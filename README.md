SNPRelate: Parallel computing toolset for relatedness and principal component analysis of SNP data
===

Version: 0.99.1

[![Build Status](https://travis-ci.org/zhengxwen/SNPRelate.png)](https://travis-ci.org/zhengxwen/SNPRelate)


## Features

Genome-wide association studies are widely used to investigate the genetic basis of diseases and traits, but they pose many computational challenges. We developed SNPRelate (R package for multi-core symmetric multiprocessing computer architectures) to accelerate two key computations on SNP data: principal component analysis (PCA) and relatedness analysis using identity-by-descent measures. The kernels of our algorithms are written in C/C++ and highly optimized.

## Wiki
[Wiki Page](https://github.com/zhengxwen/SNPRelate/wiki)

## News in Bioconductor version:

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

* Development version from Github (v0.99.1):
```
library("devtools")
install_github("zhengxwen/gdsfmt")
install_github("zhengxwen/SNPRelate")
```
The `install_github()` approach requires that you build from source, i.e. `make` and compilers must be installed on your system -- see the R FAQ for your operating system; you may also need to install dependencies manually.

* Install the packages (gdsfmt and SNPRelate) from the source code (v0.99.0):
[gdsfmt](https://codeload.github.com/zhengxwen/gdsfmt/tar.gz/v1.0.5)
and
[SNPRelate](https://codeload.github.com/zhengxwen/SNPRelate/tar.gz/v0.99.0)
```
wget https://codeload.github.com/zhengxwen/gdsfmt/tar.gz/v1.0.5 -O gdsfmt_1.0.5.tar.gz
wget https://codeload.github.com/zhengxwen/SNPRelate/tar.gz/v0.99.0 -O SNPRelate_0.99.0.tar.gz
** Or **
curl https://codeload.github.com/zhengxwen/gdsfmt/tar.gz/v1.0.5 -o gdsfmt_1.0.5.tar.gz
curl https://codeload.github.com/zhengxwen/SNPRelate/tar.gz/v0.99.0 -o SNPRelate_0.99.0.tar.gz

** Install **
R CMD INSTALL gdsfmt_1.0.5.tar.gz
R CMD INSTALL SNPRelate_0.99.0.tar.gz
```
