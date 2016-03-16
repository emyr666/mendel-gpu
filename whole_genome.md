# Introduction #

Be sure you have completed the steps in the [Configuration Page](configuration.md).  There are four major components to this pipeline. The first two involve prepping the reference and the study data.  The last two involve pre-phasing and imputation.  You'll want to make a workspace directory and launch all your processes from here.  Let's say it is in /scratch/test.

# Details #

## Splitting the reference panels into sub-regions ##

In order to fit datasets into memory, our analyses will be run sequentially on sub-regions.  First, we'll let Mendel-GPU chunk up the reference haplotypes, in our case 1000 Genomes haplotypes, into sub-regions.  You can configure how this script works by editing two parameters near the top of the script, or leave the defaults alone.  An excerpt is shown as
```
# set the following to the minimum minor allele counts observed at a site to avoid including extremely rare SNPs in the imputation/phasing
mincounts=1
# set the following to the maximum size in base pairs for each sub region
regionsize=100000
```

To run the script, execute:
```
[garykche@epigraph test]$ loop_split_kgp_files.sh 
Unzipping raw data for Chr 1
Pruning rare SNPs
Splitting into regions
Unzipping raw data for Chr 2
Pruning rare SNPs
Splitting into regions
...
Unzipping raw data for Chr 22
Pruning rare SNPs
Splitting into regions
Unzipping raw data for Chr X_nonPAR
Pruning rare SNPs
Splitting into regions
Unzipping raw data for Chr X_PAR1
Pruning rare SNPs
Splitting into regions
Unzipping raw data for Chr X_PAR2
Pruning rare SNPs
Splitting into regions
[garykche@epigraph test]$
```

Note that this step only has to be done once. You can subsequently repeat the following steps for each study you wish to impute.

## Fetching study data from the database ##

Now run loop\_fetch.sh to fetch all data from a desired study (e.g. labc) and _automatically\_chunk this data up across subjects (if it's a large study) and SNPs.  You can configure how big each dimension is by setting parameters at the top of this script.
```
# set the following parameter to the maximum number of subjects per chunk
chunksize=1000
# set the following parameter to the maximum size of a sub-region in base pairs
regionsize=10000
```_

The script is now run as:
```
[garykche@epigraph test]$ loop_fetch.sh 
Usage <study name (database name)> <include X? valid values are [1|0]>
[garykche@epigraph test]$ loop_fetch.sh labc 1
Fetching data on labc for Chr 1
Fetching data on labc for Chr 2
...
Fetching data on labc for Chr 22
Fetching data on labc for Chr X
[garykche@epigraph test]$
```

## Pre-phase the study data ##

We generally want to phase the typed SNPs into haplotypes before we impute into a much larger reference panel.  This saves time and is also recommended by authors of other imputation programs.  The following scripts handles intersecting the KGP and the study data, and launches the Mendel-GPU core code for inferring phase.
```
[garykche@epigraph test]$ loop_phasing.sh 
Usage <dbname> <use X?[1|0]> <person_chunk_start> <person_chunk_end>
[garykche@epigraph test]$ loop_phasing.sh labc 1 0 1
Prepping input for Chr 22 Region 1 Subjectchunk 0
Beginning phasing
Prepping input for Chr 22 Region 1 Subjectchunk 1
Beginning phasing
Prepping input for Chr X Region 1 Subjectchunk 0
Beginning phasing
Prepping input for Chr X Region 2 Subjectchunk 0
Beginning phasing
Prepping input for Chr X Region 1 Subjectchunk 1
Beginning phasing
Prepping input for Chr X Region 2 Subjectchunk 1
Beginning phasing
[garykche@epigraph test]$ 
```
Notice that this script asks for the first person chunk and the last person chunk to run.  In this manner, if you had multiple CPUs, or GPUs, each processor can be tasked with running the whole-genome-phasing procedure.

After loop\_phasing.sh is completed on all subject chunks, you will want to merge these files in preparation for the fast step of haploid imputation.  Run this as:
```
[garykche@epigraph test]$ loop_merge_phased.sh 
<study>
[garykche@epigraph test]$ loop_merge_phased.sh labc
[garykche@epigraph test]$
```

## Do the imputation using reference panels ##

This is very analogous to the previous step.  Note that if you have all males in the sample, it is recommended that for the X (haploid) chromosome, valid genotypes are coded as either 0s or 2s.  A setting of '1' for the third parameter in the script enables this.

```
[garykche@epigraph test]$ loop_imputation.sh
<study> <use X?[0|1]> <haploid X?[0|1]> <startchunk> <endchunk>
[garykche@epigraph test]$ loop_imputation.sh labc 1 0 0 1
[garykche@epigraph test]$ 
```