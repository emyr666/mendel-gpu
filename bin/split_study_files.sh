#!/bin/bash

if [ $# -lt 3 ] ; then
  echo "<study> <chunks> <maleonly=1>"
  exit 1
fi
study=$1
chunks=$2
maleonly=$3
header=0

chrs='1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X'
for chr in $chrs
do
  if [ $chr == 'X' && $maleonly == '1' ] ; then
    echo Haploid X
    $haploid=1
  else
    $haploid=0
  fi
  for((chunk=0;chunk<$chunks;++chunk))
  do
    basefile=$study'/chr.'$chr'.chunk.'$chunk
    gunzip -c -v $basefile'.gz' | ./split_file.pl $basefile'.part' 10000 $header $haploid
    gzip -v -f $basefile'.part'*
  done
done
