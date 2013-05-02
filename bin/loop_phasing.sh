#!/bin/bash

if [ $# -lt 4 ] ; then
  echo "Usage <dbname> <use X?[1|0]> <person_chunk_start> <person_chunk_end>"
  exit 1
fi

dbname=$1
useX=$2
chunk_start=$3
chunk_end=$4

if [ $useX -eq 1 ] ; then
  chroms='1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X'
elif [ $useX -eq 0 ] ; then
  chroms='22'
  #chroms='1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22'
else
  echo use X is 1 or 0
  exit 1
fi

maxpart=20
for chr in $chroms
do
  run_phasing.sh $dbname $chr $chunk_start $chunk_end
done
