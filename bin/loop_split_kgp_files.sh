#!/bin/bash

if [ $# -lt 2 ] ; then
  echo "<minimum minor allele counts> <regionsize>"
  exit 1
fi
mincounts=$1
regionsize=$2

chrs='1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X_nonPAR X_PAR1 X_PAR2'

for chr in $chrs do
  split_kgp_file.sh $chr $mincounts $regionsize
done
