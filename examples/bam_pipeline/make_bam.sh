#!/bin/bash

if [ $# -lt 1 ] ; then
  echo "<Subject index>"
  exit 1
fi
subject=$1

# tuning parameters
topsnps=100
regionlen='11e3'
coverage=4
echo "Running bam workflow for subject $subject for top $topsnps SNPs with coverage: $coverage region length: $regionlen."
# input locations
bin_dir='../../src/bam_utility'
variants_dir='variants3'
reffasta='simref2.fasta'
snpfile='snpinfo3.txt'
# output files
vcf_dir='vcffiles'
compressed_bam_dir='../bam/bamfiles'
subject_fasta='subject'$subject'.fasta'
subject_reads='subject'$subject
sorted='subject'$subject'_sorted'
compressed='compressed_subject'$subject
subject_vcf='subject'$subject'.vcf'

gunzip -v -c ref_fasta.gz > $reffasta
toprows=`expr $topsnps + 1`
paste $snpfile $variants_dir/variant.$subject|head -n $toprows | $bin_dir/make_subject_fasta $regionlen $reffasta > $subject_fasta
samtools faidx $subject_fasta
./art_illumina --paired --in $subject_fasta --out 'subject'$subject --len 100 --fcov $coverage --mflen 200 --sdev 3.5 -sam
rm -f *.aln *.fq 
rm -f $subject_fasta*
$bin_dir/postprocess_art_sam $subject < $subject_reads.sam | samtools view -bT $reffasta - > $subject_reads.bam
samtools sort $subject_reads.bam $sorted
samtools index $sorted.bam
rm $subject_reads.?am
head -n $toprows $snpfile | sed '1d'  | sed 's/\ /\t/g' | awk '{print $1,$3,"0",$2,$4,$5;}' | sed 's/\ /\t/g' > bamsnps.bim 
tempsnplist=$snpfile
# make compressed form
$bin_dir/compress $subject $tempsnplist $sorted.bam $reffasta 2> compressed.debug > $compressed.sam
rm compressed.debug 
#samtools view -bT $reffasta $compressed.sam > $compressed_bam_dir/$compressed.bam
#rm $compressed.sam
#samtools index $compressed_bam_dir/$compressed.bam
# make VCF pass2
samtools mpileup -uf $reffasta $sorted.bam | bcftools view - | $bin_dir/filter_vcf $tempsnplist > $subject_vcf
rm $sorted.bam*
rm $reffasta

./bgzip -f $subject_vcf
./tabix -p vcf $subject_vcf.gz
mv $subject_vcf* $vcf_dir
