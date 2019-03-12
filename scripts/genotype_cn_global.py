import pandas as pd
import glob
import re

def parse_filename(filename):
    no_path = re.split("/", filename)[-1]
    no_ext = re.split("\.", no_path)[0]
    return no_ext

def main():
    import argparse
    p = argparse.ArgumentParser("Genotypes copy number for all samples in a folder")
    p.add_argument('--path', help = "Path to copy number files")
    p.add_argument('--genes', help = "Bed file with genes and coordinates to genotype")
    args = p.parse_args()

    # Reads all files in folder with CN.bed extension
    mylist = [f for f in glob.glob(args.path+"/*.CN.bed")]

    # Reads genotypable region in bed format
    genes = pd.read_csv(args.genes, delimiter = '\t', names = ['chrom','chromStart','chromEnd','geneName'])
    output = genes

    # For each sample, envxtracts copy number for each region and calculates its mean
    # Borders of the region are left out
    for file in mylist:
        sample = parse_filename(file)
        copynumber = pd.read_csv(file, delimiter = "\t", names = ['chrom','chromStart','chromEnd','copyNumber'])
        cnv = list()
        for index, row in genes.iterrows():
            a = copynumber[copynumber['chrom'] == row['chrom']]
            b = a[(a['chromStart'] >= row['chromStart']) & (a['chromEnd'] <= row['chromEnd'])]
            c = b['copyNumber'].mean()
            cnv.append(c)
        output[sample] = pd.Series(cnv)

    # Saves output
    output.to_csv("samples_cnv.tsv", sep = '\t', header = True, index = False)

if __name__ == '__main__':
    main()
