#!/bin/bash

if [ $# -lt 3 ] ; then
  echo "<chr> <minimum minor allele counts> <regionsize>"
  exit 1
fi
chr=$1
mincounts=$2
regionsize=$3


cd $IMPUTATION_DATA
echo Unzipping raw data
gunzip -c $KGP_DIR'/ALL_1000G_phase1integrated_v3_chr'$chr'_impute.hap.gz' > temp.hap
gunzip -c $KGP_DIR'/ALL_1000G_phase1integrated_v3_chr'$chr'_impute.legend.gz' > temp.legend
echo Pruning rare SNPs
prune_kgp_file.pl $mincounts temp.legend temp.hap $chr'.legend' $chr'.hap'
rm -f $KGP_DIR'/'*'.chr'$chr'.part'*'.gz'
echo Splitting into regions
split_file.pl $KGP_DIR'/hap.chr'$chr'.part' $regionsize 0 < $chr'.hap'
gzip -f $KGP_DIR'/hap.chr'$chr'.part'*
split_file.pl $KGP_DIR'/legend.chr'$chr'.part' $regionsize 1 < $chr'.legend'
gzip -f $KGP_DIR'/legend.chr'$chr'.part'*
