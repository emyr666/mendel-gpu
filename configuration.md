# Introduction #

Mendel-GPU includes many utilities that facilitate whole-genome imputation using 1000 Genomes Project (KGP) data.  However, you will need to have a mysql server installed either on the analysis machine or a remote machine.  This will contain the positions of the SNP IDs and their positions according the KGP.

# Environment variables #

You will need to set a few environment variables in your startup script under your home directory.  Below is an excerpt from my .bash\_profile startup script.

```
export IMPUTATION_SOFTWARE=/home/garykche/gary/code/haplotyping/google/swift-gpu
export IMPUTATION_DATA=/scratch/swift-gpu
export KGP_DIR=$IMPUTATION_DATA'/kgp_v3'
export DBHOST=sailfish.usc.edu
export DBUSER=garyc
export DBPW=mypassword
export MASTERDB=genome_builds
export MASTERTABLE=grch37_unique
```

The first two variables are somewhat straightforward.  The former variable defines where you installed the program, and the latter variable defines where you plan for the imputation data to reside.

The remaining variables are only necessary for whole-genome imputation using our pipeline:

The third variable defines where to find the KGP legend and haplotype files.  You can download these from the IMPUTE2 website, and save them in $KGP\_DIR.

The last few database relevant variabltes inform scripts as to where to find data from the mysql database.  In our case, the database server is at host sailfish.usc.edu, the username with read access to all tables is garyc, and his password is mypassword.  The MASTERDB and MASTERTABLE refer to the name of the database and the table that KGP SNP IDs are located in for the purpose of sorting study genotype data in position order.

You'll want to make a copy of $IMPUTATION\_SOFTWARE/settings.imputation and $IMPUTATION\_SOFTWARE/settings.phasing into the directory you plan to run Mendel-GPU in.  These provide default values that you can change for customization.

# GPU support #

If a GPU is present, you will want to make sure you install the OpenCL drivers.  Depending on the make of your GPU, you can obtain more information at the nVidia or ATI website.  Once the driver is installed, test to see that the OpenCL runtime can see the GPU on your system.  In this example, we installed the ATI OpenCL driver, and ran CLInfo, which lists the available devices, along with their specifications:

```
[garykche@epigraph x86_64]$ pwd
/home/garykche/ati-stream-sdk-v2.2-lnx64/samples/opencl/bin/x86_64
[garykche@epigraph x86_64]$ ./CLInfo 
Number of platforms:				 2
  Platform Profile:				 FULL_PROFILE
  Platform Version:				 OpenCL 1.0 CUDA 4.0.1
  Platform Name:					 NVIDIA CUDA
  Platform Vendor:				 NVIDIA Corporation
  Platform Extensions:			 cl_khr_byte_addressable_store cl_khr_icd cl_khr_gl_sharing cl_nv_compiler_options cl_nv_device_attribute_query cl_nv_pragma_unroll 
  Platform Profile:				 FULL_PROFILE
  Platform Version:				 OpenCL 1.1 AMD-APP-SDK-v2.4 (595.10)
  Platform Name:					 AMD Accelerated Parallel Processing
  Platform Vendor:				 Advanced Micro Devices, Inc.
  Platform Extensions:			 cl_khr_icd cl_amd_event_callback cl_amd_offline_devices


  Platform Name:					 NVIDIA CUDA
Number of devices:				 2
  Device Type:					 CL_DEVICE_TYPE_GPU
  Device ID:					 4318
  Max compute units:				 14
  Max work items dimensions:			 3
    Max work items[0]:				 1024
    Max work items[1]:				 1024
    Max work items[2]:				 64
  Max work group size:				 1024
  Preferred vector width char:			 1
  Preferred vector width short:			 1
  Preferred vector width int:			 1
  Preferred vector width long:			 1
  Preferred vector width float:			 1
  Preferred vector width double:		 1
  Max clock frequency:				 1147Mhz
  Address bits:					 32
  Max memory allocation:			 704495616
  Image support:				 Yes
  Max number of images read arguments:	 128
  Max number of images write arguments:	 8
  Max image 2D width:			 4096
  Max image 2D height:			 32768
  Max image 3D width:			 2048
  Max image 3D height:	 2048
  Max image 3D depth:			 2048
  Max samplers within kernel:		 16
  Max size of kernel argument:			 4352
  Alignment (bits) of base address:		 4096
  Minimum alignment (bytes) for any datatype:	 128
  Single precision floating point capability
    Denorms:					 Yes
    Quiet NaNs:					 Yes
    Round to nearest even:			 Yes
    Round to zero:				 Yes
    Round to +ve and infinity:			 Yes
    IEEE754-2008 fused multiply-add:		 Yes
  Cache type:					 Read/Write
  Cache line size:				 128
  Cache size:					 229376
  Global memory size:				 2817982464
  Constant buffer size:				 65536
  Max number of constant args:			 9
  Local memory type:				 Scratchpad
  Local memory size:				 49152
  Profiling timer resolution:			 1000
  Device endianess:				 Little
  Available:					 Yes
  Compiler available:				 Yes
  Execution capabilities:				 
    Execute OpenCL kernels:			 Yes
    Execute native function:			 No
  Queue properties:				 
    Out-of-Order:				 Yes
    Profiling :					 Yes
  Platform ID:					 0x3eb61e0
  Name:						 Tesla C2050
  Vendor:					 NVIDIA Corporation
  Driver version:				 270.27
  Profile:					 FULL_PROFILE
  Version:					 OpenCL 1.0 CUDA
  Extensions:					 cl_khr_byte_addressable_store cl_khr_icd cl_khr_gl_sharing cl_nv_compiler_options cl_nv_device_attribute_query cl_nv_pragma_unroll  cl_khr_global_int32_base_atomics cl_khr_global_int32_extended_atomics cl_khr_local_int32_base_atomics cl_khr_local_int32_extended_atomics cl_khr_fp64 
  Device Type:					 CL_DEVICE_TYPE_GPU
  Device ID:					 4318
  Max compute units:				 14
  Max work items dimensions:			 3
    Max work items[0]:				 1024
    Max work items[1]:				 1024
    Max work items[2]:				 64
  Max work group size:				 1024
  Preferred vector width char:			 1
  Preferred vector width short:			 1
  Preferred vector width int:			 1
  Preferred vector width long:			 1
  Preferred vector width float:			 1
  Preferred vector width double:		 1
  Max clock frequency:				 1147Mhz
  Address bits:					 32
  Max memory allocation:			 704495616
  Image support:				 Yes
  Max number of images read arguments:	 128
  Max number of images write arguments:	 8
  Max image 2D width:			 4096
  Max image 2D height:			 32768
  Max image 3D width:			 2048
  Max image 3D height:	 2048
  Max image 3D depth:			 2048
  Max samplers within kernel:		 16
  Max size of kernel argument:			 4352
  Alignment (bits) of base address:		 4096
  Minimum alignment (bytes) for any datatype:	 128
  Single precision floating point capability
    Denorms:					 Yes
    Quiet NaNs:					 Yes
    Round to nearest even:			 Yes
    Round to zero:				 Yes
    Round to +ve and infinity:			 Yes
    IEEE754-2008 fused multiply-add:		 Yes
  Cache type:					 Read/Write
  Cache line size:				 128
  Cache size:					 229376
  Global memory size:				 2817982464
  Constant buffer size:				 65536
  Max number of constant args:			 9
  Local memory type:				 Scratchpad
  Local memory size:				 49152
  Profiling timer resolution:			 1000
  Device endianess:				 Little
  Available:					 Yes
  Compiler available:				 Yes
  Execution capabilities:				 
    Execute OpenCL kernels:			 Yes
    Execute native function:			 No
  Queue properties:				 
    Out-of-Order:				 Yes
    Profiling :					 Yes
  Platform ID:					 0x3eb61e0
  Name:						 Tesla C2050
  Vendor:					 NVIDIA Corporation
  Driver version:				 270.27
  Profile:					 FULL_PROFILE
  Version:					 OpenCL 1.0 CUDA
  Extensions:					 cl_khr_byte_addressable_store cl_khr_icd cl_khr_gl_sharing cl_nv_compiler_options cl_nv_device_attribute_query cl_nv_pragma_unroll  cl_khr_global_int32_base_atomics cl_khr_global_int32_extended_atomics cl_khr_local_int32_base_atomics cl_khr_local_int32_extended_atomics cl_khr_fp64 


Error : atomics mismatch!
Error : Bytes mismatch!
Error : glSharing mismatch!
Error : images mismatch!
Error : printf mismatch!
Error : deviceAttributeQuery mismatch!
Failed!
  Platform Name:					 AMD Accelerated Parallel Processing
Number of devices:				 1
  Device Type:					 CL_DEVICE_TYPE_CPU
  Device ID:					 4098
  Max compute units:				 24
  Max work items dimensions:			 3
    Max work items[0]:				 1024
    Max work items[1]:				 1024
    Max work items[2]:				 1024
  Max work group size:				 1024
  Preferred vector width char:			 16
  Preferred vector width short:			 8
  Preferred vector width int:			 4
  Preferred vector width long:			 2
  Preferred vector width float:			 4
  Preferred vector width double:		 0
  Max clock frequency:				 1600Mhz
  Address bits:					 64
  Max memory allocation:			 12658999296
  Image support:				 Yes
  Max number of images read arguments:	 128
  Max number of images write arguments:	 8
  Max image 2D width:			 8192
  Max image 2D height:			 8192
  Max image 3D width:			 2048
  Max image 3D height:	 2048
  Max image 3D depth:			 2048
  Max samplers within kernel:		 16
  Max size of kernel argument:			 4096
  Alignment (bits) of base address:		 1024
  Minimum alignment (bytes) for any datatype:	 128
  Single precision floating point capability
    Denorms:					 Yes
    Quiet NaNs:					 Yes
    Round to nearest even:			 Yes
    Round to zero:				 Yes
    Round to +ve and infinity:			 Yes
    IEEE754-2008 fused multiply-add:		 No
  Cache type:					 Read/Write
  Cache line size:				 0
  Cache size:					 0
  Global memory size:				 50635997184
  Constant buffer size:				 65536
  Max number of constant args:			 8
  Local memory type:				 Global
  Local memory size:				 32768
  Profiling timer resolution:			 999848
  Device endianess:				 Little
  Available:					 Yes
  Compiler available:				 Yes
  Execution capabilities:				 
    Execute OpenCL kernels:			 Yes
    Execute native function:			 Yes
  Queue properties:				 
    Out-of-Order:				 No
    Profiling :					 Yes
  Platform ID:					 0x2b8b12e92800
  Name:						 Intel(R) Xeon(R) CPU           X5650  @ 2.67GHz
  Vendor:					 GenuineIntel
  Driver version:				 2.0
  Profile:					 FULL_PROFILE
  Version:					 OpenCL 1.1 AMD-APP-SDK-v2.4 (595.10)
  Extensions:					 cl_khr_fp64 cl_amd_fp64 cl_khr_global_int32_base_atomics cl_khr_global_int32_extended_atomics cl_khr_local_int32_base_atomics cl_khr_local_int32_extended_atomics cl_khr_int64_base_atomics cl_khr_int64_extended_atomics cl_khr_byte_addressable_store cl_khr_gl_sharing cl_ext_device_fission cl_amd_device_attribute_query cl_amd_vec3 cl_amd_media_ops cl_amd_popcnt cl_amd_printf 


Error : atomics mismatch!
Error : Bytes mismatch!
Error : glSharing mismatch!
Error : images mismatch!
Error : printf mismatch!
Error : deviceAttributeQuery mismatch!
Failed!
[garykche@epigraph x86_64]$ 
```

Note here, that it lists two platforms, and within each platform a list of devices.  The order the platforms and devices are listed will matter, when you configure the settings for Mendel-GPU.  In the excerpt above, Platform 0 is the GPU platform, and Platform 1 is the CPU platform.  In principle, Mendel-GPU will run on a multi-core CPU platform but we have not tested it extensively there, and the architectural differences with current CPUs may not lead to substantial speedups on the CPU.

Based on what you discovered from CLInfo above, you will want to set the PLATFORM\_ID and DEVICE\_ID fields of the settings.imputation and settings.phasing files to appropriate values.  In our case, we could have two working directories that apply to the two GPUs on our system, where for each of them, the settings file has PLATFORM\_ID set to 0, and DEVICE\_ID set to 0 and 1 for each respective directory.

# Compiling the program #

If GPU support is not enabled, skip this paragraph.  For GPU support, change to $IMPUTATION\_SOFTWARE/src and open the file locations.gpu.  In there, you will want to set CL\_ROOT to the directory where you installed the OpenCL SDK. Then, in makefile, comment out the line that reads "include locations.baseline" and uncomment out the line that reads "include locations.gpu".

Finally, compile the program by running "make" at $IMPUTATION\_SOFTWARE/src.  You can also compile the utilities necessary for the pipeline utilities, by running "make utility".

# Section on configuring for whole-genome pipeline #

## Initializing the KGP database ##

You can create MASTERTABLE by using the script below:

```
use MASTERDB

CREATE TABLE `MASTERTABLE` (  
`rsnumber` varchar(32) NOT NULL default '',  
`chrom` varchar(2) default NULL, 
`position` int(11) default NULL, 
PRIMARY KEY  (`rsnumber`), 
KEY `chrom` (`chrom`,`position`));
```

where you replace MASTERDB and MASTERTABLE with their respective values.

You will want to load the values of the KGP SNPs and their positions into MASTER TABLE using the LOAD DATA command.  An example is:

LOAD DATA INFILE '/home/garyc/kgp\_snps.txt' into table MASTERTABLE;

where kgp\_snps.txt is tab delimited with rows looking like:
```
rs100     1    100312
rs102     1    100402
...
```

Please refer to the MySQL documentation for more details on this syntax.

## Initializing the study database ##

Suppose we have a study called mystudy.  We will create a database called mystudy:

```
CREATE DATABASE mystudy;
USE mystudy;
```

Now we need to create two tables that Mendel-GPU will read, in order to chunk up data for whole genome imputation.  First, create a sample list:

```
CREATE TABLE `person` (
  `seq` mediumint unsigned NOT NULL default '0',
  `id` varchar(100) default NULL,
  `sex` char(1) default NULL,
  `perm` mediumint(8) unsigned default NULL,
  PRIMARY KEY  (`seq`),
  KEY `index_perm` (`perm`)
)
```

Below is an excerpt of example data:

```
mysql> select * from person limit 5;
+-----+---------+------+------+
| seq | id      | sex  | perm |
+-----+---------+------+------+
|   0 | 5112641 | F    |  837 | 
|   1 | 5117125 | F    | 4197 | 
|   2 | 5117343 | F    | 1486 | 
|   3 | 5118531 | F    |  151 | 
|   4 | 5119660 | F    | 2415 | 
+-----+---------+------+------+
5 rows in set (0.05 sec)

mysql> 
```

The perm column simply allows MENDEL-GPU to permute the order of the samples when chunking up data across subjects.  This can reduce the possibility of batch effects.  The lowest value for perm must be zero and the highest value is N-1 where N is the sample size.  It is optional and you can simply copy seq into perm if you wish.

The second table contains all the genotypes.  Create the genotype table as:

```
CREATE TABLE `genotype` (
  `rs` varchar(20) NOT NULL default '',
  `vector` mediumtext,
  PRIMARY KEY  (`rs`)
) ENGINE=MyISAM
```

A tab delimited text file that you could import into genotype would look something like

```
rs10231 03123
```

In reality, the second column as many more characters (thousands), each character representing the genotype of the ordered individuals, ordered by the sequence defined in the person table.  We define missing genotype as 0, homozygous reference (allele corresponding to reference genome allele) genotype as 1, hetero genotype as 2, and homozygous other (other allele) genotype as 3.

Here we only import one SNP, where the 0-th (first) person (ID 5112641) has missing genotype, and 4-th (fifth) person (5119660) has homozygous other genotype.

Please visit the [tutorial](whole_genome.md) for whole-genome analysis for more information on how to do large scale phasing, and imputation.