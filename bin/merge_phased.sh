#!/bin/bash

if [ $# -lt 2 ] ; then
  echo "<study> <chr>"
  exit 1
fi

study=$1
chr=$2

chunks=10

phased_results=$IMPUTATION_DATA'/'$study'_phased'
imputed_results=$IMPUTATION_DATA'/'$study'_imputed'
snplist_temp=$phased_results'/snplist.merged.'$chr
phasing_temp=$phased_results'/phasing.merged.'$chr
# MERGE ALL THE PHASED OUTPUTS INTO INPUTFILES
rm -f $snplist_temp $phasing_temp.*
parts=30
for((part=1;part<$parts;++part))
do
  snpfile=$IMPUTATION_DATA'/'$study'/snplist.chr.'$chr'.chunk.0.part'$part
  if [ -f $snpfile ] ; then
    #echo $snpfile
    # concatenate all chr parts's IDs into this file
    cat $snpfile  >> $snplist_temp
    snps=`cat $snpfile|wc -l`
    # concatenate all chr parts's phasing data into each chunk
    for((chunk=0;chunk<$chunks;++chunk))
    do
      phasefile=$phased_results'/phasing.chr'$chr'.part'$part'.chunk'$chunk
      if [ -f $phasefile'.gz' ] ; then
        echo "Unzipping $phasefile"
        cols=`gunzip -c  $phasefile'.gz' | head -n1 | wc -w`
        persons=`echo "($cols-1)/4"|bc`
        echo Making hap posteriors for $snps SNPs and $persons persons.
        convert2log=0
        confidence=1
        gunzip -c $phasefile'.gz' | make_haplo_posteriors $persons $snps $convert2log $confidence >> $phasing_temp.$chunk
      fi
    done
  fi
done
gzip -f $phasing_temp.*
