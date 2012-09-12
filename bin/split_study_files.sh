#!/bin/bash

if [ $# -lt 2 ] ; then
  echo "<study> <chunks>"
  exit 1
fi
study=$1
chunks=$2

chrs='1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X'
for chr in $chrs
do
  for((chunk=0;chunk<$chunks;++chunk))
  do
    basefile=$study'/chr.'$chr'.chunk.'$chunk
    gunzip -c -v $basefile'.gz' | ./split_file.pl $basefile'.part' 10000 0 
    gzip -v -f $basefile'.part'*
  done
done
