#!/bin/bash

WD=$1

FILES=$WD/*CEL

echo $FILES

for f in $FILES
do
  echo $f
  echo cac
done


while read line
do
  echo $line
done <  $WD/mcf_cel_name_change.tsv
