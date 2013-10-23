#!/bin/bash

faidx='simref2.fasta.fai'
mergedsam='merged.sam'
rm -f $faidx $mergedsam
subjects=100
tar -xzvf variants3.tar.gz
for((subject=0;subject<$subjects;++subject))
do
  sed "s/SUBJECT/$subject/" $faidx.template >> $faidx
  ./make_bam.sh $subject
  cat 'compressed_subject'$subject'.sam' >> $mergedsam
done
rm -f compressed*sam
samtools view -bT simref2.fasta merged.sam > merged.bam
samtools index merged.bam
