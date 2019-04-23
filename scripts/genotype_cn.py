import pandas as pd
import numpy as np

def main():
    import argparse
    p = argparse.ArgumentParser("Genotypes copy number in a sample by averaging copy number in a region")
    p.add_argument('--sample', help = "Sample name")
    p.add_argument('--genes', help = "Bed file with genes and coordinates to genotype")
    p.add_argument('--copynumber', help = "Bed file with sample's copy numbers per window")
    args = p.parse_args()

    # Reads input files into a pandas data frame
    genes = pd.read_csv(args.genes, delimiter = '\t', names = ['chrom','chromStart','chromEnd','geneName'])
    copynumber = pd.read_csv(args.copynumber, delimiter = "\t", names = ['chrom','chromStart','chromEnd','copyNumber'])
    output = genes
    cnv = list()

    # Extracts copy number for each region and calculates its mean
    # Borders of the region are left out
    for index, row in genes.iterrows():
        a = copynumber[copynumber['chrom'] == row['chrom']]
        b = a[(a['chromStart'] >= row['chromStart']) & (a['chromEnd'] <= row['chromEnd'])]
        c = b['copyNumber'].mean()
        cnv.append(c)
    output['raw_cnv'] = pd.Series(cnv)
    output['round_cnv'] = pd.Series(round(x) for x in cnv)

    # Saves output
    output.to_csv(args.sample+"_cnv.tsv", sep = '\t', header = False, index = False)

if __name__ == '__main__':
    main()
