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

# CREATES DIRS STRUCTURE FOR GENOME ANAYSIS BASED ON GATK AND MUTECT PIPELINES
#
# USAGE: create_dir_tree.sh <path to dir>
#
#Example: ./create_dir_tree.sh /media/DATA1/DATA/SEQ/Projects/2015_07_23_Exome                                                                                                                                                                     
###START SCRIPT

if [ "$#" -ne 1 ]; then
    echo "Illegal number of parameters."
    echo "Usage: create_dir_tree.sh <path_to_dir>"
    exit
fi

cd "$1"

mkdir Analysis

cd Analysis
mkdir fastq
mkdir Coverage
mkdir Data-cleanup
mkdir Variant-Discovery
mkdir Results


mkdir ./Variant-Discovery/Logs
mkdir ./Results/Logs

cd Data-cleanup
mkdir Aligned
mkdir BQSR
mkdir IndelRealigned
mkdir Logs
mkdir Marked
mkdir Sorted
mkdir Trimmed


### END SCRIPT
