#!/bin/bash

# set this to location of swift-gpu/bin directory
#scripts=$HOME'/gary/code/haplotyping/google/swift-gpu/bin'

if [ $# -lt 5 ] ; then
  echo "Inputformat: <chr><legend><hapfile> [plink|mendel|vcf|mec|phasing] <basefilename>";
  exit 1
fi

chr=$1
legend=$2
hapfile=$3
format=$4
basename=$5
geno_dim=3
if [ $format = 'plink' ] ; then
  prepare_inputfiles $chr $legend $hapfile 0 $basename.fam $basename.bim $basename.bed
elif [ $format = 'mec' ] ; then
  prepare_inputfiles $chr $legend $hapfile 2 $basename.mec 
elif [ $format = 'phasing' ] ; then
  prepare_inputfiles $chr $legend $hapfile 3 $basename.snplist $basename.phasing
  geno_dim=4
elif [ $format = 'mendel' ] ; then
  convertmeta_mendel2plink.pl $basename.ped $basename.def $basename.fam $basename.bim
  prepare_inputfiles $chr $legend $hapfile 0 $basename.fam $basename.bim $basename.bed
elif [ $format = 'vcf' ] ; then
  gunzip -c $basename.vcf.gz |  gawk -f $IMPUTATION_SOFTWARE/bin/convertvcf.awk > $basename.glf
  cols=`head -n1 $basename.glf|wc -w`
  echo "Total cols in glf is $cols"
  cut -f1-5 $basename.glf > glf.info
  cut -f6-$cols $basename.glf > glf.data
  echo Running prepare_inputfiles
  prepare_inputfiles $chr $legend $hapfile 1 glf.info glf.data
  #rm -f glf.info glf.data $basename.glf
else
  exit 1
fi

# below will store $basename data in compressed format for swift-gpu
echo Compressing floats
compress_float $geno_dim < out.genotype.data > genotype.stream
