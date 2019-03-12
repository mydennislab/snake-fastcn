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
export PATH=/share/dennislab/programs/new_miniconda3/bin/:$PATH
source activate snakecn
export PATH="/share/dennislab/programs/fastCN:/share/dennislab/programs/mrsfast/:$PATH"
```

> Python scripts must have execution permissions set (change with `chmod`)

## Download reference

Using Jeff Kidd's reference (3 Kb windows):
```bash
wget http://guest:kiddlab@kiddlabshare.umms.med.umich.edu/shared-data/public-data/fastCN/GRCh38_BSM_fastCN.tgz
tar -xvzf GRCh38_BSM_fastCN.tgz
rm GRCh38_BSM_fastCN.tgz
```

You also need a file containing chrom sizes. We generated a custom file containing chrom sizes for that reference.

We also created a referece with window size 1 Kb as follows:
```bash
intersectBed -v -wao  -a GRCh38_bsm.1kb.bed -b GRCh38-badRegions.bed.sorted.merge >  GRCh38_bsm.1kb.control
awk '{if($1=="chrX"){print}}' GRCh38_bsm.1kb.control > GRCh38_bsm.1kb.bed.chrXnonParControl
awk '{if($1 !~ "_" && $1 != "chrX" && $1 != "chrM" && $1 != "chrY" && $1 != "chrEBV"){print}}' GRCh38_bsm.1kb.control > GRCh38_bsm.1kb.bed.autoControl
```

> The custom reference genome for this analysis is located in: `/share/dennislab/databases/assemblies/GRCh38/GRCh38_BSM_fastCN`

## Running pipeline with Snakemake

This pipeline needs a file containing complete link addresses for all fastq files related to that sample.

Example:
```
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA18507/sequence_read/ERR002346_2.filt.fastq.gz
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA18507/sequence_read/ERR002346_1.filt.fastq.gz
```

Then just run:
```bash
snakemake -p --config sample="sample_name" urls="filename.urls" reference_path="path/to/referece" chrom_sizes="path/to/chromsizes"
```

Example: (running with 10 cores maximum)
```bash
snakemake -p --config sample=NA18507 urls=NA18507.urls reference_path=/share/dennislab/databases/assemblies/GRCh38/GRCh38_BSM_fastCN chrom_sizes=/share/dennislab/databases/assemblies/GRCh38/GRCh38_BSM_fastCN/ref/GRCh38_BSM.chromsizes -j 10
```

## Pipeline overview

![Rulegraph](figs/rulegraph.png)

## Limitations

- This pipeline uses only paired-end reads with the extensions "_1" and "_2". This can be modified in the future (if we want to use single-end reads).
- This pipeline is **not paralog specific**

## Extra scripts for genotyping

Additional scripts to genotype copy number for certain genes of interest are provided in the scripts folder. 
- `genotype_cn.py` receives a bed file with regions for copy number genotyping, sample name and its respective bed file with copy number estimates
- `genotype_cn_global.py` receives a bed file with regions for copy number genotyping, and path that contains copy number estimates for one or several individuals. The script will automatically read all files with extension "CN.bed" and output a tsv matrix with copy number variants for each individual

To run these scripts create a python3 environment with pandas installed.

Example genotyping:
```bash
python scripts/genotype_cn.py --sample NA18507 --copynumber windows/NA18507.depth.1kb.bed.CN.bed --genes data/genotypable_regions.bed 
python scripts/genotype_cn_global.py --path windows/ --genes data/genotypable_regions.bed
```

Results will be stored in `NA18507_cnv.tsv` and `samples_cnv.tsv` respectively.
