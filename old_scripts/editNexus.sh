#!/bin/bash

INPATH=$1
OUTPATH=$1
#suffix=ed.txt

ls $INPATH | grep "124" | grep -Ev "annot|BP|CC|MF" > $OUTPATH/scriptList

k=0
while read line; do
    sheet+=($line)
    ((k++))
done < $OUTPATH/scriptList

for (( c=0; c<$k; c++ )); do

    sheetName=$(echo ${sheet[c]} | cut -d '.' -f1)
    awk 'NR>7' $INPATH/${sheet[c]} |
    awk -F\t '{OFS="\t"}{for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = "."};{print}' |
    awk -F\t '{OFS="\t"}{if(match($15,"mi")) {} else $15="."; print}' |
    awk '{gsub(/ - /, "_")};{print}' |
    awk '{sub(/Region/, "chr\tstart\tend")};{print}' |
    awk '{gsub(/Gene Symbols/, "Gene Symbol")};{print}' |
    awk -F\t 'NR>1 {OFS="\t"}{sub(/:/, "\t", $1)};{print}' |
    awk -F\t 'NR>1 {OFS="\t"}{sub(/,/, "", $2)};{print}' |
    awk -F\t 'NR>1 {OFS="\t"}{sub(/,/, "", $2)};{print}' |
    awk -F\t 'NR>1 {OFS="\t"}{sub(/-/, "\t", $2)};{print}' |
    awk -F\t 'NR>1 {OFS="\t"}{sub(/,/, "", $3)};{print}' |
    awk -F\t 'NR>1 {OFS="\t"}{sub(/,/, "", $3)};{print}' |
    awk -F\t ' {if(match($15,", ")) {split($15,a,", ") ; for (i=1; i<=length(a); i++) {sub($15, a[i]) ; print}} else print}' > $OUTPATH/${sheetName}_ed.txt

done

rm $OUTPATH/scriptList
unset sheet





ls $INPATH | grep "annot" > $OUTPATH/scriptList

k=0
while read line; do
    sheet+=($line)
    ((k++))
done < $OUTPATH/scriptList

for (( c=0; c<$k; c++ )); do

    sheetName=$(echo ${sheet[c]} | cut -d '.' -f1)
    awk 'NR>1' $INPATH/${sheet[c]} |
    awk -F\t '{OFS="\t"}{for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = "."};{print}' > $OUTPATH/${sheetName}_ed.txt

done

rm $OUTPATH/scriptList
