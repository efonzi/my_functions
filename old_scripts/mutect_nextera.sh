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


### MUTECT pipeline ###

 
# Input files are pre-processed as GATK data-cleanup pipeline reccomandations
# command example: ./mutect.sh /media/eugenio/DATA2/NGS/WES/Projects/LAL/first_batch/Analysis


####################

# START SCRIPT

####################

# Dirs and files

if [ "$#" -ne 1 ]; then
    echo "Illegal number of parameters"
    echo "Usage: ./mutect.sh BASEDIR"
    exit 1
fi

echo "basedir: $1"

inpath=$1/Data-cleanup/BQSR
outpath=$1/Results
refpath=~/giovanni/wes_lal

reference=human_g1k_v37.fasta
exome_targeted_regions=nexterarapidcapture_expandedexome_targetedregions_fixed.bed
dbsnp=dbsnp_138.b37.vcf
cosmicdb=b37_cosmic_v54_120711.vcf

# Software

JAVA6=/home/condivisi/packages/jdk/jdk1.6.0_35/jre/bin
MUTECT=/home/condivisi/softwares/mutect

# Redirecting

exec 2>&1 > $outpath/Logs/err_log.txt

echo "------------------------------------------\n"
echo "--                                      --\n"
echo "-- Starting analysis                    --\n"
echo "--                                      --\n"
echo "------------------------------------------\n"


# Preparing files list

ls -l $inpath | awk '{print $9}' | grep -E .bam | grep -E germline > $outpath/normal_samplesList

ls -l $inpath | awk '{print $9}' | grep -E .bam | grep -E tumor > $outpath/tumor_samplesList

k=0
while read line; do
        normalbam+=($line)
        ((k++)) 
done < $outpath/normal_samplesList

k=0
while read line; do
        tumorbam+=($line)
        ((k++)) 
done < $outpath/tumor_samplesList


# Mutect

for (( c=0; c<$k; c++ )); do
    
    echo ${normalbam[c]}
    echo ${tumorbam[c]}
    #echo call_stats_$((c+1))

    $JAVA6/java -Xmx2g -jar $MUTECT/muTect-1.1.4.jar --analysis_type MuTect --reference_sequence $refpath/$reference --cosmic $refpath/$cosmicdb --dbsnp $refpath/$dbsnp --intervals $refpath/$exome_targeted_regions --input_file:normal $inpath/${normalbam[c]}  --input_file:tumor $inpath/${tumorbam[c]}  --out $outpath/call_stats_$((c+1)).out --coverage_file $outpath/coverage.wig.$((c+1)).txt 
    
    
done


