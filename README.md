# Snakemake file to determine copy number based on read-depth 

This pipeline uses Jeff Kidd's tool FastCN tool to predict copy number in a region based on read-depth of Illumina reads.

## FastCN installation

```bash
git clone https://github.com/KiddLab/fastCN.git
cd fastCN
g++ -o GC_control_gen GC_control_gen.cc
g++ -o SAM_GC_correction SAM_GC_correction.cc
gcc -std=c99 depth_combine.c -O3 -o depth_combine
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

> Python scripts must have execution permissions set (change with `chmod`)

## Download reference

Using Jeff Kidd's reference:
```bash
wget http://guest:kiddlab@kiddlabshare.umms.med.umich.edu/shared-data/public-data/fastCN/GRCh38_BSM_fastCN.tgz
tar -xvzf GRCh38_BSM_fastCN.tgz
rm GRCh38_BSM_fastCN.tgz
```

You also need a file containing chrom sizes. 

> We generated a custom file containing chrom sizes for that reference.

## Running pipeline with Snakemake

This pipeline needs a file containing complete link addresses for all fastq files related to that sample. 

Example:
```
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA18507/sequence_read/ERR002346_2.filt.fastq.gz
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA18507/sequence_read/ERR002351_2.filt.fastq.gz
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA18507/sequence_read/ERR002346.filt.fastq.gz
...
``` 

Then just run:
```bash
snakemake -p --config sample="sample_name" urls="filename.urls" reference_path="path/to/referece" chrom_sizes="path/to/chromsizes" 
```

Example: (running with 10 cores maximum)
```bash
snakemake -p --config sample=NA18507 urls=NA18507.urls reference_path=/share/dennislab/databases/assemblies/GRCh38/GRCh38_BSM_fastCN chrom_sizes=/share/dennislab/databases/assemblies/GRCh38/GRCh38_BSM_fastCN/ref/GRCh38_BSM.chromsizes -j 10
```

## Limitations

- This pipeline uses only paired-end reads with the extensions "_1" and "_2". This can be modified in the future (if we want to use single-end reads).

