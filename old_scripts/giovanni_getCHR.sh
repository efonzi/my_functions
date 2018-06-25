#!/bin/bash

### given a list of gene symbols, query NCBI and obtain chromosome location

### example usage: getCHR.bash list.of.genes.txt

list=$(cat $1)

#echo $list
#> chr.csv

for gene in $list
do
#	echo $gene
	esearch -db gene -query "($gene[gene]) AND (Homo sapiens[orgn])" | efetch -db clinvar -format docsum | xtract -pattern DocumentSummary -element Chromosome -element MapLocation | head -n1 #>> chr.csv

done
