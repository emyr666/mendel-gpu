#!/bin/bash

# set the following to the minimum minor allele counts observed at a site to avoid including extremely rare SNPs in the imputation/phasing
mincounts=1
# set the following to the maximum size in base pairs for each sub region
regionsize=100000

chrs='1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X_nonPAR X_PAR1 X_PAR2'

for chr in $chrs
do
  split_kgp_file.sh $chr $mincounts $regionsize
done
