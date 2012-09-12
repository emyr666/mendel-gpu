#!/bin/bash

if [ $# -lt 4 ] ; then
  echo "Usage <dbname> <chrom> <person_chunk_start> <person_chunk_end>"
  exit 1
fi

dbname=$1
chr=$2
chunk_start=$3
chunk_end=$4

maxpart=20
ln -fs $IMPUTATION_SOFTWARE/src/kernels.c 
for((chunk=$chunk_start;chunk<=$chunk_end;++chunk))
do
  basefile=$IMPUTATION_DATA'/'$dbname'/chr.'$chr'.chunk.'$chunk
  #gunzip -c -v $basefile'.gz' | split_file.pl $basefile'.part' 10000 0
  #gzip -v -f $basefile'.part'*
  for((part=1;part<=$maxpart;++part))
  do
    if [[ -f $IMPUTATION_DATA'/'$dbname'/chr.'$chr'.chunk.'$chunk'.part'$part'.gz' && ! -f $IMPUTATION_DATA'/'$dbname'_phased/results/phasing.chr'$chr'.part'$part'.chunk'$chunk'.gz' ]] ; then
      echo Prepping input for Chr $chr Region $part Subjectchunk $chunk
      gunzip -c $IMPUTATION_DATA'/'$dbname'/chr.'$chr'.chunk.'$chunk'.part'$part'.gz' > study.mec
      # do conversion
      prepare_denovo 2 study.mec
      compress_float 3 < out.genotype.data > genotype.stream
      cut -f2 out.genotype.info  |sed '1d' > $IMPUTATION_DATA'/'$dbname'/snplist.chr.'$chr'.chunk.'$chunk'.part'$part
      # do phasing
      people=`head -n1 out.genotype.data|cut -f1`
      snps=`head -n1 out.genotype.data|cut -f2`
      total_regions=`echo "$snps/10000"|bc`
      remainder=`echo "$snps%10000"|bc`
      if [ $remainder -gt 0 ] ; then
        let total_regions=$total_regions+1
      fi
      #echo $people $snps $total_regions $remainder
      cat settings.phasing | sed "s/people/$people/" | sed "s/snps/$snps/" | sed "s/total_regions/$total_regions/" > settings.txt
      echo Beginning phasing
      impute 1>impute.out 2>impute.debug
      # post process files
      cat POSTERIORS |gzip -c - > $IMPUTATION_DATA'/'$dbname'_phased/phasing.chr'$chr'.part'$part'.chunk'$chunk'.gz'
    fi
  done
  #rm -f $basefile'.part'*
done
