#!/bin/bash

# first argument must be the working directory
# second argument must be the desired file extension
# third argument must be the file containing the name conversion 
# (old names in the first column and new names in the second column)
# the conversion file must be located in working directory
WD=$1
EXT=$2
CONV=$3

# create directory for output
mkdir $WD/output

# loop through lines of the conversion file
while read line
do
  old=$(echo $line | awk '{print $1}') #extract old name from first column of conversion file
  new=$(echo $line | awk '{print $2}') #extract new name from second column of conversion file

  # check if a file named like "old name" is present in the working directory
  if [ -f $WD/$old$EXT ]
  # if present, change its name to the new name and move it to the "output" dir
  then 
    mv $WD/$old$EXT $WD/output/$new$EXT
  # if not present, print error message
  else 
    echo "cannot find "$old$EXT
    echo "          "
  fi

done <  $WD/$CONV


echo "###############"
echo "end of script"
echo "###############"

