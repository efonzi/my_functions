#!/bin/bash

home='/volume1'
#dest=$home/'Fontana_Marconi_et_al_2017_1'
dest=$home/'180613_chromot_snparray'
source=$home/'DATA/unibo/SNParray'
prova=$home/'prova'

for i in $( ls $dest | grep CEL ); do
    echo $source/$i
    rm $dest/$i
    ln -s $source/$i $dest/$i
done
