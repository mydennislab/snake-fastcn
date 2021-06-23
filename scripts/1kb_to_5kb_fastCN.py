## Author: Mira Mastoras mnmastoras@ucdavis.edu
## Description: This script combines 5 consecutive bed windows and averages their copy number (column 4)
## Usage: python 1kb_to_5kb_fastCN.py amargosa.depth.1000.bed.CN.bed > amargosa.5kb.bed

import pandas as pd
import numpy as np
import sys

bedfile = pd.read_csv(sys.argv[1], sep='\t', header=None)

# function to combine 5 consecutive windows in bedfile and average their CN
def combine_windows(bed):
    nrows=len(bed.index)
    index = 0
    if nrows == 1: # if window spans the entire chromosome, just print it
        bed_print=bed.iloc[0].tolist()
        print(*bed_print, sep="\t")
    for n in range(nrows-1):
        if index == nrows: # special case where chrom has exact number of windows divisible by 5
            break
        elif index + 4 >= nrows: # case of the last window in chromsome being less than 5
            offset = nrows-index # offset is new number of windows to combine (less than 5 bc we hit the end of the chrom)
            bedsum = 0
            for i in range(offset):
                bedsum += bed.iloc[index + i,3]
            CN = bedsum / offset
            start = bed.iloc[index,1]
            end = bed.iloc[index+offset-1, 2]
            chr = bed.iloc[index,0]
            print(chr, "\t", start, "\t", end, "\t", CN)
            break
        else:
            CN = (bed.iloc[index,3] + bed.iloc[index+1,3] + bed.iloc[index+2,3]+ bed.iloc[index+3,3]+ bed.iloc[index+4,3]) / 5
            start = bed.iloc[index,1]
            end = bed.iloc[index+4, 2]
            chr = bed.iloc[index,0]
            print(chr, "\t", start, "\t", end, "\t", CN)

            index += 5

# Get chromosome names
bedfile[0]=bedfile[0].astype(str)
chroms = np.unique(bedfile[0])

# combine windows for each chromosome:

for chr in chroms:
    combine_windows(bedfile.loc[bedfile[0]==chr])
