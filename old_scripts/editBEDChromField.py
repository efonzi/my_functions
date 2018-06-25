# The WES reference BED files from Illumina have the chromosome names as 'chr#', 'chrX', 'chrY', 'chrM' and before running the pipeline we must change them to '1', '2', .. , 'Y', 'MT'

# usage: python3 editBEDChromField.py arg1 [arg2]
# arg1 is mandatory and is the name of the BED file
# arg2 is optional and is an integer defining the number of the line to skip when loading the file

import pandas as pd
import sys

args = sys.argv

# get name of BED file from command line
bed_name = args[1]

# some files have a first line that we have to skip
# if no 2nd argument is provided, skip no lines
# otherwise skip the line number provided in the 2nd argument
rows_to_skip = None
if len(args) > 2:
    rows_to_skip = [int(r) for r in args[2]]

# load file
bed = pd.read_table(bed_name, sep='\t', header=None, dtype='str', skiprows=rows_to_skip)

# change 'chr#' to '#' and 'chrM' to 'MT'
bed[0] = ['MT' if c == 'chrM' else c[3:] for c in bed[0]]

# cut extension
bed_name = bed_name[:-4]

# write to file
bed.to_csv(bed_name + '_edited.bed', sep='\t', header=False, index=False)
