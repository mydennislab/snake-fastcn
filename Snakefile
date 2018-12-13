"""
Title: SnakeFastCN
Desc: Snakemake pipeline for copy number determination using Illumina reads
- This pipeline is based in Jeff Kidd's copy number calculation pipeline
"""

import os
import re
import pandas as pd

# -----------
# Config-file
# -----------

TSV = config['tsv']
REFERENCE_PATH = config["reference_path"]
CHROM_SIZES = config["chrom_sizes"]

# -----------------------
# Variables and functions
# -----------------------

df = pd.read_csv(TSV, delimiter = '\t')
URLS = df['url'].tolist() # urls of fastq files to download
FILENAMES = [os.path.basename(x) for x in URLS if "_1" in x or "_2" in x] # fastq filenames to download
SAMPLES = list(set(df['Sample'].tolist())) # individuals in tsv
ACCESSIONS = list(set([x.split("_")[0] for x in FILENAMES]))
EXTENSION = list(set([x.split(".",1)[1] for x in FILENAMES]))

# ---------------
# WORKFLOW SET-UP
# ---------------

TARGETS = expand("windows/{acc}.depth.3kb.bed.CN.bb", acc = ACCESSIONS)

rule all:
  input:
    TARGETS

# -------------
# Data download
# -------------

rule download_fastq:
  input:
    TSV
  output:
    expand("fastq/{filename}", filename = FILENAMES)
  shell:
    '''
    cut -f1 {input} | tail -n +2 | xargs -P4 wget -P fastq
    '''

# -----------------------
# Reference pre-treatment
# -----------------------

rule gc_correction_ref:
  input:
    fasta = REFERENCE_PATH+"/ref/GRCh38_BSM.fa",
    exclude = REFERENCE_PATH+"/ref/toexclude.bed.sorted.merge.sorted2",
    gaps = REFERENCE_PATH+"/ref/GRCh38_bsm.gaps.bed.slop36.sorted.merged.sorted2"
  output:
    bin = REFERENCE_PATH+"/ref/GRCh38_BSM.GC_control.bin"
  shell:
    '''
    GC_control_gen {input.fasta} {input.exclude} {input.gaps} 400 {output.bin}
    ''' 

rule index_masked:
  input:
    fasta = REFERENCE_PATH+"/masked/GRCh38_BSM.fa"
  output:
    index = REFERENCE_PATH+"/masked/GRCh38_BSM.fa.index"
  shell:
    '''
    /share/dennislab/programs/mrsfast/mrsfast --index {input.fasta}
    '''

# --------------------
# Mapping to reference
# --------------------

rule mapping:
  input:
    forward = expand("fastq/{{acc}}_1.{ext}", ext = EXTENSION),
    reverse = expand("fastq/{{acc}}_2.{ext}", ext = EXTENSION),
    reference = REFERENCE_PATH+"/masked/GRCh38_BSM.fa",
    gccontrol = REFERENCE_PATH+"/ref/GRCh38_BSM.GC_control.bin"
  output:
    "mapping/{acc}.sam.gz"
  params:
    "mapping/{acc}"
  threads: 5
  shell:
    '''
    extract-from-fastq36-pair.py --in1 {input.forward} --in2 {input.reverse} | /share/dennislab/programs/mrsfast/mrsfast --search {input.reference} --seq /dev/fd/0 --disable-nohits --mem 16 --threads {threads} -e 2 --outcomp -o {params}
    '''

# -----------------------------------
# GC-corrected read-depth calculation 
# -----------------------------------

rule gc_correction_sam:
  input:
    alignment = "mapping/{acc}.sam.gz",
    reference = REFERENCE_PATH+"/ref/GRCh38_BSM.fa.fai",
    gccontrol = REFERENCE_PATH+"/ref/GRCh38_BSM.GC_control.bin"
  output:
    bpdepth = "binary/{acc}.bin.gz"
  params:
    "binary/{acc}"
  shell:
    '''
    zcat {input.alignment} | SAM_GC_correction {input.reference} {input.gccontrol} /dev/fd/0 {params}
    gzip {params}.bin
    '''

rule bp_to_windows:
  input:
    bpdepth = "binary/{acc}.bin.gz",
    chromlen = REFERENCE_PATH+"/ref/GRCh38_BSM.fa.fai", 
    windows = REFERENCE_PATH+"/windows/GRCh38_bsm.3kb.bed"
  output:
    windepth = "windows/{acc}.depth.3kb.bed"
  shell:
    '''
    perbp-to-windows.py --depth {input.bpdepth} --out {output.windepth} --chromlen {input.chromlen} --windows {input.windows}
    '''

# -------------------------------------
# Read-depth to copy number calculation
# -------------------------------------

rule depth_to_cn:
  input:
     depth = "windows/{acc}.depth.3kb.bed",
     auto = REFERENCE_PATH+"/windows/GRCh38_bsm.3kb.bed.autoControl",
     chrx = REFERENCE_PATH+"/windows/GRCh38_bsm.3kb.bed.chrXnonParControl"
  output:
     "windows/{acc}.depth.3kb.bed.CN.bed"
  shell:
    '''
    depth-to-cn.py --in {input.depth} --autocontrol {input.auto} --chrX {input.chrx}
    '''

# ------------------------------------------------
# Merge samples by individual  + bigBed conversion
# ------------------------------------------------

rule bed2bigBed:
  input:
    bed = "windows/{acc}.depth.3kb.bed.CN.bed",
    chromsizes = CHROM_SIZES
  output:
    sorted = temp("windows/{acc}.depth.3kb.bed.CN.srt.bed"),
    bigbed = "windows/{acc}.depth.3kb.bed.CN.bb"
  shell:
    '''
    sort -k1,1 -k2,2n {input.bed} > {output.sorted}
    bedToBigBed {output.sorted} {input.chromsizes} {output.bigbed}  
    '''
 
