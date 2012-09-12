In this directory we have small examples extracted from real data of Chr 22 from the KGP.  The files are described below:

REFERENCE HAPLOTYPE FILES
*************************

kgp.haplotypes:

This contains phased haplotypes with alleles 0 and 1.  The rows are SNPs, and columns are haplotypes (2 * individuals).  There are 1092 individuals and 70 SNPs in this example.

kgp.legend:

This is meta data describing the SNPs.

FILES DERIVED FROM GWAS CHIP STUDY
****************************************

study.bim and study.def:

Meta data on the each SNP. 
For PLINK compliant format (.bim): Chr,ID,Genetic position,Physical position,Ref allele, Other allele
For MENDEL compliant format (.def): ID,Chr,Physical position

study.fam and study.ped

Meta data on each subject.
For PLINK compliant format (.fam): Ped ID, Indiv ID, Father ID, Mother ID, Sex (1=male,2=female), Affection (-9=unknown)
For MENDEL compliant format (.ped): Ped ID, Indiv ID, Father ID, Mother ID, Sex (1=male,2=female), Affection (-9=unknown)

study.bed

This is a binary file that stores genotypes at the bit level.  More information can be found at the PLINK website.

FILES DERIVED FROM GWAS LOW-PASS SEQUENCING STUDY
*************************************************

chr22.vcf.gz:

The Variant Calling Format file contains many bits of information. The relevant bits are extracted by piping the contents through the convertvcf.awk script as demonstrated in convert.sh
