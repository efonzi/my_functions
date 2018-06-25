#!/bin/bash


# Copyright (C) 2015  Marco Manfrini, PhD

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


### SAMTOOLS short indels pipeline ###

 
# Input files are pre-processed as GATK data-cleanup pipeline reccomandations
# usage example: ./varscan_paired.sh /media/eugenio/DATA2/NGS/WES/Projects/LAL/first_batch/Analysis

####################

# START SCRIPT

####################

# Dirs and files

if [ "$#" -ne 1 ]; then
    echo "Illegal number of parameters"
    echo "Usage: ./varscan_paired.sh BASEDIR"
    exit 1
fi

echo "basedir: $1"

inpath=$1/Data-cleanup/BQSR
outpath=$1/Results/
refpath=/media/DATA1/DATA/Genomes

reference=human_g1k_v37.fasta
exome_targeted_regions=truseq-exome-targeted-regions-manifest-v1-2_edited.bed

# Software

SAMTOOLS=/media/DATA1/DATA/Software/jre1.6.0_45/bin
ANNOVAR=/media/DATA1/DATA/Software/annovar/
VARSCAN=/media/DATA1/DATA/Software/
 

# Redirecting

exec 2>&1 > $outpath/Logs/varscan_err_log.txt

echo "------------------------------------------\n"
echo "--                                      --\n"
echo "-- Starting analysis                    --\n"
echo "--                                      --\n"
echo "------------------------------------------\n"


#Preparing files list

ls -l $inpath | awk '{print $9}' | grep -E .bam  > $outpath/samtoolsList

k=0
while read line; do

        bamFiles+=($line)
        ((k++)) 

done < $outpath/samtoolsList

for (( c=0; c<$k; c++ )); do
    
    # SAMTOOLS VARSCAN

    echo ${bamFiles[c]}

    samtools mpileup -B -q 1 -f $refpath/$reference --positions $refpath/$exome_targeted_regions $inpath/${bamFiles[c]} > $outpath/${bamFiles[c]}.mpileup

done

ls -l $outpath | awk '{print $9}' | grep -E .mpileup  > $outpath/varscanList

k=0
while read line; do

        mpileupFiles+=($line)
        ((k++)) 

done < $outpath/varscanList

#k=1 TEST PURPOSE ONLY
j=1
for (( c=0; c<$k; c++ )); do

    if (( $c % 2 ==  0 )); then

	echo "Comparing tumor vs normal samples.."

	echo $outpath/${mpileupFiles[c]}

	echo $outpath/${mpileupFiles[c+1]}

	java -jar $VARSCAN/VarScan.v2.3.9.jar somatic $outpath/${mpileupFiles[c]} $outpath/${mpileupFiles[c+1]} --output-snp $outpath/${mpileupFiles[c+1]}.snp.pairs.tsv --output-indel $outpath/${mpileupFiles[c+1]}.indels.pairs.tsv --min-avg-qual 15 --strand-filter 1 --min-var-freq 0.05 --somatic-p-value 0.05

	((j++))

    fi

done

# INTRODUCTED IN SAMTOOLS MPILEUP PARAMETERS -B AND -q 1
# CALL TO VARSCAN SOMATIC WITH --min-avg-qual 15 --strand-filter 1 --min-var-freq 0.05 --somatic-p-value 0.05. BEFORE WAS --p-value 0.01 --somatic-p-value 0.01
