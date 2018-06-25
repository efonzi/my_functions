#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Illegal number of parameters"
    echo "Usage: ./data_cleanup.sh BASEDIR"
    exit 1
fi


if [ $1 = "truseq" ]; then
	p=trus
	echo $p
elif [ $1 = "nextera" ]; then
	p=next
	echo $p
else
	echo "error"
	exit
fi

echo "last line"

#####################################################################

#for i in $( ls ~/aneuploidy/circos/170621/ | grep conf ); do
#        echo $i
#done
#
#for i in $( ls ~/aneuploidy/circos/170621/ | grep conf ); do
#	mv ~/aneuploidy/circos/170621/$i ~/aneuploidy/circos/170621/${i#170621_}
#done"
