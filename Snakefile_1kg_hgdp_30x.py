"""
Title: SnakeFastCN
Desc: Snakemake pipeline for copy number determination using Illumina reads
- This pipeline is based on Jeff Kidd's copy number calculation pipeline
- This version contains updates to run specifically on the HGDP and 1KG 30x high coverage data
  located at /mnt/datasets/dennislab/1kg /mnt/datasets/dennislab/hgdp
- These files are named as {sample}.cram, and are single end, so I removed the filename parsing and combine depth steps
"""

import os
import re
import pandas as pd

# -----------
# CONFIG FILE
# -----------

SAMPLE = [line. rstrip('\n') for line in open(config['samples'])]
REFERENCE_PATH = config['reference_path']
CHROM_SIZES = config['chrom_sizes']

WINDOWSIZE = config['window_size'] # must be 1kb or 3kb


# ---------------
# WORKFLOW TARGET
# ---------------

rule all:
  input:
    expand("bigBed/{smp}.depth."+WINDOWSIZE+".bed.CN.bb", smp = SAMPLE)

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
# CONVERT TO FASTQ
# --------------------

rule cram_convert:
  input:
    "cram/{smp}.cram"
  output:
    temp("fastq/{smp}.fastq.gz")
  shell:
    """
    samtools fastq {input} | gzip -c > {output}
    """

# --------------------
# MAPPING TO REFERENCE
# --------------------

rule mapping:
  input:
    fastq = "fastq/{smp}.fastq.gz",
    reference = REFERENCE_PATH+"/ref-WMDUST-masked/GRCh38_BSM_WMDUST_masked.fa"
  output:
    "mapping/{smp}.sam.gz"
  params:
    "mapping/{smp}"
  threads: 5
  shell:
    '''
    extract-from-fastq36.py --in {input.fastq} | /share/dennislab/programs/mrsfast/mrsfast --search {input.reference} --seq /dev/fd/0 --disable-nohits --mem 16 --threads {threads} -e 2 --outcomp -o {params}
    '''

# -----------------------
# GC-CORRECTED READ-DEPTH
# -----------------------

rule gc_correction_sam:
  input:
    alignment = "mapping/{smp}.sam.gz",
    reference = REFERENCE_PATH+"/ref-WMDUST/GRCh38_BSM_WMDUST_unmasked.fa.fai",
    gccontrol = REFERENCE_PATH+"/ref-WMDUST/GRCh38_BSM_WMDUST.GC_control.bin"
  output:
    bpdepth = "binary/{smp}.bin.gz"
  params:
    "binary/{smp}"
  shell:
    '''
    zcat {input.alignment} | SAM_GC_correction {input.reference} {input.gccontrol} /dev/fd/0 {params}
    gzip {params}.bin
    '''


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
