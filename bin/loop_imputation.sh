#!/bin/bash

if [ $# -lt 4 ] ; then
  echo "<study> <use X?[0|1]> <startchunk> <endchunk>"
  exit 1
fi

study=$1
useX=$2
startchunk=$3
endchunk=$4

if [ $useX -eq 1 ] ; then
  chrs='1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X'
elif [ $useX -eq 0 ] ; then
  chrs='1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22'
else
  echo use X is 1 or 0
  exit 1
fi


for chr in $chrs
do
  run_imputation.sh $study $chr $startchunk $endchunk
done
