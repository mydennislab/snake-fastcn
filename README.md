# Snakemake file to determine copy number based on read-depth 

This pipeline uses Jeff Kidd's tool FastCN tool to predict copy number in a region based on read-depth of Illumina reads.

## FastCN installation

```bash
git clone https://github.com/KiddLab/fastCN.git
cd fastCN
g++ -o GC_control_gen GC_control_gen.cc
g++ -o SAM_GC_correction SAM_GC_correction.cc
```

## Set-up environment

Before running the Snakafile you need to have in your path:
- fastCN 
- MrsFast
- bedToBigBed
- Python 2 with pandas, numpy and matplotlib libraries

The best way to do this is to create just a Conda environment:
```bash
conda create -n snakecn python=2.7 pandas numpy matplotlib ucsc-bedToBigBed
```

FastCN and MrsFast should be manually added to your path. 

Example: activating the environment and adding fastCN and MrsFast to the path:
```bash
source activate snakecn
export PATH="/share/dennislab/programs/fastCN:/share/dennislab/programs/mrsfast/:$PATH"
```

## Download reference

Using Jeff Kidd's reference:
```bash
wget http://guest:kiddlab@kiddlabshare.umms.med.umich.edu/shared-data/public-data/fastCN/GRCh38_BSM_fastCN.tgz
tar -xvzf GRCh38_BSM_fastCN.tgz
rm GRCh38_BSM_fastCN.tgz
```

You also need a file containing chrom sizes. 

> We generate a custom file containing chrom sizes for that reference.

## Running pipeline with Snakemake

This pipeline only needs a TSV file containing fastq ulrs as well as the sample/individual name (like the ones that you can download from 1K genomes).

Then just run:
```bash
snakemake -p --config tsv="filename.tsv" reference_path="path/to/referece" chrom_sizes="path/to/chromsizes" 
```

Example: (running with 10 cores maximum)
```bash
snakemake -p --config tsv=igsr_NA18507_undefined.tsv reference_path=/share/dennislab/databases/assemblies/GRCh38/GRCh38_BSM_fastCN chrom_sizes=/share/dennislab/databases/assemblies/GRCh38/GRCh38_BSM_fastCN/ref/GRCh38_BSM.chromsizes -j 15
```

Overview of rules:

![Rule graph](figs/rulegraph.png)

## Limitations and future directions

Limitations:
- It works for paired-end reads only
- It considers that all files in TSV have the same extension 

Future directions:
- Merging copy number from the same individual

