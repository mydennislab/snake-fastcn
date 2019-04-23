"""
Title: SnakeFastCN
Desc: Snakemake pipeline for copy number determination using Illumina reads
- This pipeline is based in Jeff Kidd's copy number calculation pipeline
- This version contains an update in combine_depth to check for the case of one url
- This version skips fast q download step for files where fast q was manually downloaded
"""

import os
import re
import pandas as pd

# -----------
# CONFIG FILE
# -----------

SAMPLE = config['sample']
URLS = config['urls']
REFERENCE_PATH = config['reference_path']
CHROM_SIZES = config['chrom_sizes']

WINDOWSIZE = '1kb'

# -----------------
# VARIABLES PARSING
# -----------------

def parse_filename(filename):
  accession = re.split("_1|_2", filename, maxsplit = 2)[0]
  extension = re.split("_1|_2", filename, maxsplit = 2)[1]
  return accession, extension

urlsfile = pd.read_csv(URLS, delimiter = '\n', names = ['url'])
urls = urlsfile['url'].tolist()
filenames = urls
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

# rule download_fastq:
  #output:
   # expand("fastq/{filename}", filename = filenames)
 # threads: 10
 # shell:
  #  '''
   # cat {URLS} | xargs -n 1 -P 10 wget --directory-prefix=fastq
   # '''

# --------------------
# REFERENCE PROCESSING
# --------------------

rule gc_correction_ref:
  input:
    fasta = REFERENCE_PATH+"/ref-WMDUST/GRCh38_BSM_WMDUST_unmasked.fa",
    exclude = REFERENCE_PATH+"/ref-WMDUST/GRCh38_BSM_WMDUST.exclude.bed.sort2",
    gaps = REFERENCE_PATH+"/ref-WMDUST/GRCh38_BSM_WMDUST.gaps.bed.slop36.sorted.merged.sort2"
  output:
    bin = REFERENCE_PATH+"/ref-WMDUST/GRCh38_BSM_WMDUST.GC_control.bin"
  shell:
    '''
    GC_control_gen {input.fasta} {input.exclude} {input.gaps} 400 {output.bin}
    '''

rule index_masked:
  input:
    fasta = REFERENCE_PATH+"/ref-WMDUST-masked/GRCh38_BSM_WMDUST_masked.fa"
  output:
    index = REFERENCE_PATH+"/ref-WMDUST-masked/GRCh38_BSM_WMDUST_masked.fa.index"
  shell:
    '''
    /share/dennislab/programs/mrsfast/mrsfast --index {input.fasta}
    '''

# --------------------
# MAPPING TO REFERENCE
# --------------------

def accession_forward(wildcards):
  return "fastq/"+wildcards.acc+"_1"+extensions.get(wildcards.acc)

def accession_reverse(wildcards):
  return "fastq/"+wildcards.acc+"_2"+extensions.get(wildcards.acc)

rule mapping:
  input:
    forward = accession_forward,
    reverse = accession_reverse,
    reference = REFERENCE_PATH+"/ref-WMDUST-masked/GRCh38_BSM_WMDUST_masked.fa"
  output:
    "mapping/{acc}.sam.gz"
  params:
    "mapping/{acc}"
  threads: 5
  shell:
    '''
    extract-from-fastq36-pair.py --in1 {input.forward} --in2 {input.reverse} | /share/dennislab/programs/mrsfast/mrsfast --search {input.reference} --seq /dev/fd/0 --disable-nohits --mem 16 --threads {threads} -e 2 --outcomp -o {params}
    '''

# -----------------------
# GC-CORRECTED READ-DEPTH
# -----------------------

rule gc_correction_sam:
  input:
    alignment = "mapping/{acc}.sam.gz",
    reference = REFERENCE_PATH+"/ref-WMDUST/GRCh38_BSM_WMDUST_unmasked.fa.fai",
    gccontrol = REFERENCE_PATH+"/ref-WMDUST/GRCh38_BSM_WMDUST.GC_control.bin"
  output:
    bpdepth = "binary/{acc}.bin"
  params:
    "binary/{acc}"
  shell:
    '''
    zcat {input.alignment} | SAM_GC_correction {input.reference} {input.gccontrol} /dev/fd/0 {params}
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

rule bp_to_windows:
  input:
    bpdepth = "binary/{smp}.bin.gz",
    chromlen = REFERENCE_PATH+"/ref-WMDUST/GRCh38_BSM_WMDUST_unmasked.fa.fai",
    windows = REFERENCE_PATH+"/windows-WMDUST/GRCh38_BSM_WMDUST."+WINDOWSIZE+".bed"
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
     auto = REFERENCE_PATH+"/windows-WMDUST/GRCh38_BSM_WMDUST."+WINDOWSIZE+".autoControl.bed",
     chrx = REFERENCE_PATH+"/windows-WMDUST/GRCh38_BSM_WMDUST."+WINDOWSIZE+".chrXnonParControl.bed"
  output:
     "windows/{smp}.depth."+WINDOWSIZE+".bed.CN.bed"
  shell:
    '''
    depth-to-cn.py --in {input.depth} --autocontrol {input.auto} --chrX {input.chrx}
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