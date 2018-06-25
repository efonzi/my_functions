#!/bin/bash

# to run CIRCOS package
# run from terminal like: sh script/directory (1) (2) (3)
#(1) = working directory
#(2) = name of desired ".conf" file (that contains the instructions to build the plot)
#(3) = circos package directory up to "bin" folder
#(4) = desired name for output file


# the "working/directory" must contain a ".conf" file and a karyotype file whose name must match the name of the karyotype file indicated inside the "circos.conf" file
# the output will be written in working/directory

# directories
INPATH=$1
OUTPATH=$1
CIRCOSPATH=$3

# ".conf" name
CIRCOS=$2

# output file name
OUTNAME=$4

# "circos.conf" file contains structure info to build the CIRCOS (ideograms, spacing, etc.)
# use options -outputdir -outputfile to specify outdirectory and name of output file
$CIRCOSPATH/circos -conf $INPATH/$CIRCOS.conf -outputdir $OUTPATH -outputfile $OUTNAME
