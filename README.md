SNPRelate: Parallel computing toolset for relatedness and principal component analysis of SNP data
===

Version: 1.1.1

[![Build Status](https://travis-ci.org/zhengxwen/SNPRelate.png)](https://travis-ci.org/zhengxwen/SNPRelate)


## Features

Genome-wide association studies are widely used to investigate the genetic basis of diseases and traits, but they pose many computational challenges. We developed SNPRelate (R package for multi-core symmetric multiprocessing computer architectures) to accelerate two key computations on SNP data: principal component analysis (PCA) and relatedness analysis using identity-by-descent measures. The kernels of our algorithms are written in C/C++ and highly optimized.

## Wiki
[Wiki Page](https://github.com/zhengxwen/SNPRelate/wiki)

## Bioconductor:

[Development Version](http://www.bioconductor.org/packages/devel/bioc/html/SNPRelate.html)


## News v1.1.1:

	* fix bug in snpgdsVCF2GDS when 'method="biallelic.only"'
	* add 'snpgdsVCF2GDS_R' for the R implementation
	* fix bug in snpgdsBED2GDS if 'family=TRUE'

## News
```gdsfmt_1.1.2``` should be installed immediately, if you see the error like
```
Invalid Zip Deflate Stream operation 'Seek'!
```

## News in Bioconductor (v1.0.0) compared to CRAN version:

	* fully support long vectors (>= R v3.0)
	* >5x speedup in the function 'snpgdsVCF2GDS'
	* SNP GDS format allows character-type chromosome codes
	* add a new argument 'ref.allele' in 'snpgdsVCF2GDS'
	* add new functions 'snpgdsOpen' and 'snpgdsClose'
	* add a new function 'snpgdsTranspose' to transpose the genotypic matrix
	* add a new function 'snpgdsAlleleSwitch' to switch alleles if needed
	* add a new function 'snpgdsApartSelection'
	* add a new function 'snpgdsGEN2GDS' to import Oxford GEN data
	* use NA instead of 3 as missing value in 'snpgdsGetGeno'
	* add a new argument 'snpfirstdim' in the function 'snpgdsGDS2BED'
	* add a new argument 'with.id' in the functions 'snpgdsSNPRateFreq' and 'snpgdsSampMissRate'
	* return a numeric vector instead of data.frame in 'snpgdsLDpair'
	* add estimating nine Jacquard's coefficients in 'snpgdsIBDMLE'
	* fix the memory issues reported by valgrind
	

## Installation


* Bioconductor repository:
```
source("http://bioconductor.org/biocLite.R")
library(BiocInstaller)
BiocInstaller::useDevel()

biocLite("SNPRelate")
```


* Development version from Github:
```
library("devtools")
install_github("zhengxwen/gdsfmt")
install_github("zhengxwen/SNPRelate")
```
The `install_github()` approach requires that you build from source, i.e. `make` and compilers must be installed on your system -- see the R FAQ for your operating system; you may also need to install dependencies manually.


* Install the package from the source code:
[gdsfmt](https://github.com/zhengxwen/gdsfmt/tarball/master),
[SNPRelate](https://github.com/zhengxwen/SNPRelate/tarball/master)
```
wget --no-check-certificate https://github.com/zhengxwen/gdsfmt/tarball/master -O gdsfmt_latest.tar.gz
wget --no-check-certificate https://github.com/zhengxwen/SNPRelate/tarball/master -O SNPRelate_latest.tar.gz
** Or **
curl -L https://github.com/zhengxwen/gdsfmt/tarball/master/ -o gdsfmt_latest.tar.gz
curl -L https://github.com/zhengxwen/SNPRelate/tarball/master/ -o SNPRelate_latest.tar.gz

** Install **
R CMD INSTALL gdsfmt_latest.tar.gz
R CMD INSTALL SNPRelate_latest.tar.gz
```


* Old CRAN version (v0.9.19) from r-forge repository:
```
install.packages("gdsfmt", repos="http://R-Forge.R-project.org")
install.packages("SNPRelate", repos="http://R-Forge.R-project.org")
```
