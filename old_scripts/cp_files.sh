#!/bin/bash

# BASH SCRIPT TO CREATE  A LIST OF SOFT LINKS TO FILES IN DIR <indir> TO <outdir>
#
# USAGE
# cp_files.sh /media/DATA1/DATA/SEQ/Data/20150703_Exome/20150703_Exome/Project_AML /media/DATA1/DATA/SEQ/Project/20150703_Exome/Analysis/fastq 
#

if [ "$#" -ne 2 ]; then
    echo "Illegal number of parameters."
    echo "Usage: cp_files.sh SOURCE DEST" 
else

    ls $1 > dirList

    while read dirName; do
        echo "Uncompressing .fastq.gz files in $dirName .."  
        gzip -kd  $1/$dirName/*.fastq.gz
	
	echo "Creating links to .fastq files .."
	ln -s $1/$dirName/*.fastq $2/

	echo "Done."
    done < dirList
fi
