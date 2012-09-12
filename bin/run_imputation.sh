#!/bin/bash

if [ $# -lt 4 ] ; then
  echo "<study> <chr> <startchunk> <endchunk>"
  exit 1
fi

study=$1
chr=$2
startchunk=$3
endchunk=$4

chr_code=$chr
if [[ $chr == 'X_nonPAR' || $chr == 'X_PAR1' || $chr == 'X_PAR2' ]] ; then
  chr_code='X'
fi
phased_results=$IMPUTATION_DATA'/'$study'_phased'
imputed_results=$IMPUTATION_DATA'/'$study'_imputed'

ln -fs $IMPUTATION_SOFTWARE/src/kernels.c
mkdir -p $imputed_results

parts=30
# INTERSECT THE INPUTFILES WITH THE REF HAPLOTYPES
for((part=1;part<$parts;++part))
do
  if [[ -f $KGPDIR'/legend.chr'$chr'.part'$part'.gz' ]] ; then
    gunzip -c $KGPDIR'/legend.chr'$chr'.part'$part'.gz' > kgp.legend
    gunzip -c $KGPDIR'/hap.chr'$chr'.part'$part'.gz' > kgp.hap
    for((chunk=$startchunk;chunk<=$endchunk;++chunk))
    do
      merged_gz=$IMPUTATION_DATA'/'$study'_phased/phasing.merged.'$chr_code'.'$chunk'.gz'
      if [[ -f $merged_gz && ! -f $imputed_results'/quality.chr'$chr'.part'$part'.chunk'$chunk'.gz' ]] ; then
        echo Prepping inputfiles for Chr $chr_code Region $part Subjectchunk $chunk
        gunzip -c $merged_gz > study.phasing
        cp $IMPUTATION_DATA'/'$study'_phased/snplist.merged.'$chr_code study.snplist
        prep_swiftgpu_input.sh $chr_code kgp.legend kgp.hap phasing study
        cut -f2 out.genotype.info |sed '1d' > $imputed_results'/snplist.chr'$chr'.part'$part
        # do imputation
        people=`head -n1 out.genotype.data|cut -f1`
        snps=`head -n1 out.genotype.data|cut -f2`
        total_regions=`echo "$snps/10000"|bc`
        remainder=`echo "$snps%10000"|bc`
        if [ $remainder -gt 0 ] ; then
          let total_regions=$total_regions+1
        fi
        echo Running haploid imputation on $people persons $snps SNPs $total_regions regions with remainder $remainder
        cat settings.imputation | sed "s/people/$people/" | sed "s/snps/$snps/" | sed "s/total_regions/$total_regions/" > settings.txt
        impute 1>impute.out 2>impute.debug
        # post process files
        cat QUALITY |gzip -c - > $imputed_results'/quality.chr'$chr'.part'$part'.chunk'$chunk'.gz'
        cat DOSAGES |gzip -c - > $imputed_results'/dosages.chr'$chr'.part'$part'.chunk'$chunk'.gz'
        cat GENOTYPES |gzip -c - > $imputed_results'/genotypes.chr'$chr'.part'$part'.chunk'$chunk'.gz'
      else
         echo "Skipping because $imputed_results'/quality.chr'$chr'.part'$part'.chunk'$chunk'.gz' exists"
      fi
    done
  fi
done # loop through chr parts
