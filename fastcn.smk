"""
Title: SnakeFastCN
Desc: Snakemake pipeline for copy number determination using Illumina reads
- This pipeline is based on Jeff Kidd's copy number calculation pipeline
- This version contains an update in combine_depth to check for the case of one url
"""

import os
import re
import pandas as pd

# -----------
# CONFIG FILE
# -----------

FASTCN = "/share/dennislab/programs/fastCN_new"
MRSFAST = "/share/dennislab/programs/mrsfast"

SAMPLE = config['sample']
URLS = config['urls']

UNMASKED_REF = config['unmasked_ref']
CHROM_SIZES = config['chrom_sizes']
MASKED_REF = config['masked_ref']
GC_EXCLUDE_BED = config['gc_exclude_bed']
GC_GAPS_BED = config['gc_gaps_bed']

WINDOWS = config['windows']
WINDOWSIZE = '1kb'

# -----------------
# VARIABLES PARSING
# -----------------

def parse_filename(filename):
  accession = re.split("_|\.", filename, maxsplit = 2)[0]
  extension = re.split("_|\.", filename, maxsplit = 2)[2]
  return accession, extension

urlsfile = pd.read_csv(URLS, delimiter = '\n', names = ['url'])
urls = urlsfile['url'].tolist()
filenames = [os.path.basename(x) for x in urls if "_1" in x or "_2" in x]
extensions = dict([parse_filename(x) for x in filenames])
accessions = list(extensions.keys())

# ---------------
# WORKFLOW TARGET
# ---------------

rule all:
  input:
    expand("bigBed/{smp}.depth."+WINDOWSIZE+".bed.CN.bb", smp = SAMPLE)

# -------------
# DATA DOWNLOAD
# -------------

rule download_fastq:
  output:
    expand("fastq/{filename}", filename = filenames)
  threads: 10
  shell:
    '''
    cat {URLS} | xargs -n 1 -P {threads} wget --directory-prefix=fastq
    '''

# ------------------   
# REFERENCE INDEXING
# ------------------

rule fai_unmasked:
  input:
    fasta = UNMASKED_REF
  output:
    fai = UNMASKED_REF+".fai"
  shell:
    '''
    samtools faidx {input.fasta}
    '''
   
rule index_masked:
  input:
    fasta = MASKED_REF
  output:
    index = MASKED_REF+".index"
  shell:
    '''
    {MRSFAST}/mrsfast --index {input.fasta}
    '''
  
# ------------------------
# GC CORRECTION GENERATION
# ------------------------

rule gc_correction_ref:
  input:
    fasta = UNMASKED_REF,
    exclude = GC_EXCLUDE_BED,
    gaps = GC_GAPS_BED
  output:
    bin = UNMASKED_REF+".GC_control.bin"
  shell:
    '''
    {FASTCN}/GC_control_gen {input.fasta} {input.exclude} {input.gaps} 400 {output.bin}
    '''

# --------------------
# MAPPING TO REFERENCE
# --------------------

def accession_forward(wildcards):
  return "fastq/"+wildcards.acc+"_1."+extensions.get(wildcards.acc)

def accession_reverse(wildcards):
  return "fastq/"+wildcards.acc+"_2."+extensions.get(wildcards.acc)

rule mapping:
  input:
    forward = accession_forward,
    reverse = accession_reverse,
    reference = MASKED_REF
  output:
    "mapping/{acc}.sam.gz"
  params:
    "mapping/{acc}"
  threads: 5
  shell:
    '''
    {FASTCN}/extract-from-fastq36-pair.py --in1 {input.forward} --in2 {input.reverse} | {MRSFAST}/mrsfast --search {input.reference} --seq /dev/fd/0 --disable-nohits --mem 16 --threads {threads} -e 2 --outcomp -o {params}
    '''

# -----------------------
# GC-CORRECTED READ-DEPTH
# -----------------------

rule gc_correction_sam:
  input:
    alignment = "mapping/{acc}.sam.gz",
    reference = UNMASKED_REF+".fai",
    gccontrol = UNMASKED_REF+".GC_control.bin"
  output:
    bpdepth = "binary/{acc}.bin"
  params:
    "binary/{acc}"
  shell:
    '''
    zcat {input.alignment} | {FASTCN}/SAM_GC_correction {input.reference} {input.gccontrol} /dev/fd/0 {params}
    '''

# ------------
# COMBINE RUNS
# ------------

rule combine_depth:
  input:
    expand("binary/{acc}.bin", acc = accessions)
  output:
    "binary/{smp}.bin.gz"
  params:
    accessions = " ".join(expand("binary/{acc}.bin", acc = accessions)),
    combined = "binary/{smp}.bin"
  run:
    if len(accessions) <= 1: # check to see if there is only one url
      shell("mv {params.accessions} {params.combined}")

    else:
      shell("depth_combine -H {params.accessions} > {params.combined}")

    shell("gzip {params.combined}")

# ------------------------------------
# CONVERT READ-DEPTH FROM BP TO WINDOW
# ------------------------------------

rule make_windows:
  input:
    
  output:
  
  shell:
    '''
    
    '''

rule bp_to_windows:
  input:
    bpdepth = "binary/{smp}.bin.gz",
    chromlen = UNMASKED_REF+".fai", 
    windows = UNMASKED+"."+WINDOWSIZE+".bed"
  output:
    windepth = "windows/{smp}.depth."+WINDOWSIZE+".bed"
  shell:
    '''
    perbp-to-windows.py --depth {input.bpdepth} --out {output.windepth} --chromlen {input.chromlen} --windows {input.windows}
    '''

# -------------------------
# READ-DEPTH TO COPY NUMBER
# -------------------------

rule depth_to_cn:
  input:
     depth = "windows/{smp}.depth."+WINDOWSIZE+".bed",
     auto = UNMASKED_REF+"."+WINDOWSIZE+".autoControl.bed",
     chrx = UNMASKED_REF+"."+WINDOWSIZE+".chrXnonParControl.bed"
  output:
     "windows/{smp}.depth."+WINDOWSIZE+".bed.CN.bed"
  shell:
    '''
    {FASTCN}/depth-to-cn.py --in {input.depth} --autocontrol {input.auto} --chrX {input.chrx}
    '''

# -----------------
# BIGBED CONVERSION
# -----------------

rule bed2bigBed:
  input:
    bedGraph = "windows/{smp}.depth."+WINDOWSIZE+".bed.CN.bed",
    chromsizes = CHROM_SIZES
  output:
    bed9 = temp("bigBed/{smp}.depth."+WINDOWSIZE+".bed.CN.bed9"),
    sorted = temp("bigBed/{smp}.depth."+WINDOWSIZE+".bed.CN.srt.bed9"),
    bigbed = "bigBed/{smp}.depth."+WINDOWSIZE+".bed.CN.bb"
  shell:
    '''
    python3 scripts/bedToBed9.py {input.bedGraph} {output.bed9} 
    sort -k1,1 -k2,2n {output.bed9} > {output.sorted}
    bedToBigBed -type=bed9 {output.sorted} {input.chromsizes} {output.bigbed}  
    '''

