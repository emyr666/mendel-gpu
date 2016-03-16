# Introduction #

In this tutorial we describe how one can quickly impute missing genotypes using the 1000 Genomes Project data as external information for informing haplotype structure.

We have included some utilities to convert files in standard formats (e.g. PLINK, MENDEL, VCF) into native formats.  You can run these utilities from a directory that contains sample input files.  These are located at $IMPUTATION\_HOME/examples/rawinput where $IMPUTATION\_HOME is the location of the installation.  In this directory you will find files extracted from real data of Chr 22 from the KGP.  Descriptions of the files and how to process them using the included scripts are described below:

# Reference haplotypes #

kgp.haplotypes:

This contains phased haplotypes with alleles 0 and 1.  The rows are SNPs, and columns are haplotypes.  There are 1092 individuals and 70 SNPs in this example.

kgp.legend:

This is meta data describing the SNPs.

# Data from GWAS chip #

study.bed

This is a binary file that stores genotypes at the bit level.  More information can be found at the PLINK website.

study.bim and study.def:

Meta data on the each SNP.
For PLINK compliant format (.bim): Chr,ID,Genetic position,Physical position,Ref allele, Other allele
For MENDEL compliant format (.def): ID,Chr,Physical position

study.fam and study.ped

Meta data on each subject.
For PLINK compliant format (.fam): Ped ID, Indiv ID, Father ID, Mother ID, Sex (1=male,2=female), Affection (-9=unknown)
For MENDEL compliant format (.ped): Ped ID, Indiv ID, Father ID, Mother ID, Sex (1=male,2=female), Affection (-9=unknown)

**Usage from MENDEL format input: _prep\_swiftgpu\_input.sh 22 kgp.legend kgp.haplotypes mendel study_**

**Usage from PLINK format input: _prep\_swiftgpu\_input.sh 22 kgp.legend kgp.haplotypes plink study_**

# Low-coverage resequencing study #

study.vcf.gz:

The Variant Calling Format file contains many bits of information. The relevant bits are extracted by piping the contents through the convertvcf.awk script which is called prep\_swiftgpu\_input.sh

