SNPRelate: Parallel computing toolset for relatedness and principal component analysis of SNP data
====

![GPLv3](http://www.gnu.org/graphics/gplv3-88x31.png)
[GNU General Public License, GPLv3](http://www.gnu.org/copyleft/gpl.html)

[![Availability](http://www.bioconductor.org/shields/availability/release/SNPRelate.svg)](http://www.bioconductor.org/packages/release/bioc/html/SNPRelate.html)
[![Years-in-BioC](http://www.bioconductor.org/shields/years-in-bioc/SNPRelate.svg)](http://www.bioconductor.org/packages/release/bioc/html/SNPRelate.html)
[![Build Status](https://travis-ci.org/zhengxwen/SNPRelate.png)](https://travis-ci.org/zhengxwen/SNPRelate)
[![Build status](https://ci.appveyor.com/api/projects/status/odo1jcrxg65k748g?svg=true)](https://ci.appveyor.com/project/zhengxwen/snprelate)
[![Comparison is done across all Bioconductor packages over the last 6 months](http://www.bioconductor.org/shields/downloads/SNPRelate.svg)](http://www.bioconductor.org/packages/release/bioc/html/SNPRelate.html)
[![codecov.io](https://codecov.io/github/Bioconductor-mirror/SNPRelate/coverage.svg?branch=master)](https://codecov.io/github/Bioconductor-mirror/SNPRelate?branch=master)


## Features

Genome-wide association studies are widely used to investigate the genetic basis of diseases and traits, but they pose many computational challenges. We developed SNPRelate (R package for multi-core symmetric multiprocessing computer architectures) to accelerate two key computations on SNP data: principal component analysis (PCA) and relatedness analysis using identity-by-descent measures. The kernels of our algorithms are written in C/C++ and highly optimized.

The GDS format offers the efficient operations specifically designed for integers with two bits, since a SNP could occupy only two bits. The SNP GDS format in this package is also used by the [GWASTools](http://bioconductor.org/packages/GWASTools) package with the support of S4 classes and generic functions. The extended GDS format is implemented in the [SeqArray](https://github.com/zhengxwen/SeqArray) package to support the storage of single nucleotide variation (SNV), insertion/deletion polymorphism (indel) and structural variation calls.


## Bioconductor

Release Version: v1.12.1

[http://www.bioconductor.org/packages/release/bioc/html/SNPRelate.html](http://www.bioconductor.org/packages/release/bioc/html/SNPRelate.html)

Development Version: v1.13.3

[http://www.bioconductor.org/packages/devel/bioc/html/SNPRelate.html](http://www.bioconductor.org/packages/devel/bioc/html/SNPRelate.html)


## News

* Supports the [SeqArray](http://bioconductor.org/packages/release/bioc/html/SeqArray.html) GDS format, see [the vignette](http://www.bioconductor.org/packages/release/bioc/vignettes/SeqArray/inst/doc/R_Integration.html#integration-with-snprelate).


## Tutorials

[http://corearray.sourceforge.net/tutorials/SNPRelate](http://corearray.sourceforge.net/tutorials/SNPRelate)

[http://www.bioconductor.org/packages/devel/bioc/vignettes/SNPRelate/inst/doc/SNPRelateTutorial.html](http://www.bioconductor.org/packages/devel/bioc/vignettes/SNPRelate/inst/doc/SNPRelateTutorial.html)


## Citation

Zheng X, Levine D, Shen J, Gogarten SM, Laurie C, Weir BS (2012). A High-performance Computing Toolset for Relatedness and Principal Component Analysis of SNP Data. *Bioinformatics*. [DOI: 10.1093/bioinformatics/bts606](http://dx.doi.org/10.1093/bioinformatics/bts606).

Zheng X, Gogarten S, Lawrence M, Stilp A, Conomos M, Weir BS, Laurie C, Levine D (2017). SeqArray -- A storage-efficient high-performance data format for WGS variant calls. *Bioinformatics*. [DOI: 10.1093/bioinformatics/btx145](http://dx.doi.org/10.1093/bioinformatics/btx145).


## Installation

* Bioconductor repository:
```R
source("http://bioconductor.org/biocLite.R")
biocLite("SNPRelate")
```

* Development version from Github:
```R
library("devtools")
install_github("zhengxwen/gdsfmt")
install_github("zhengxwen/SNPRelate")
```
The `install_github()` approach requires that you build from source, i.e. `make` and compilers must be installed on your system -- see the [R FAQ](http://cran.r-project.org/faqs.html) for your operating system; you may also need to install dependencies manually.


* Install the package from the source code:
[gdsfmt](https://github.com/zhengxwen/gdsfmt/tarball/master),
[SNPRelate](https://github.com/zhengxwen/SNPRelate/tarball/master)
```sh
wget --no-check-certificate https://github.com/zhengxwen/gdsfmt/tarball/master -O gdsfmt_latest.tar.gz
wget --no-check-certificate https://github.com/zhengxwen/SNPRelate/tarball/master -O SNPRelate_latest.tar.gz
R CMD INSTALL gdsfmt_latest.tar.gz
R CMD INSTALL SNPRelate_latest.tar.gz

## Or
curl -L https://github.com/zhengxwen/gdsfmt/tarball/master/ -o gdsfmt_latest.tar.gz
curl -L https://github.com/zhengxwen/SNPRelate/tarball/master/ -o SNPRelate_latest.tar.gz
R CMD INSTALL gdsfmt_latest.tar.gz
R CMD INSTALL SNPRelate_latest.tar.gz
```



## Implementation with Intel Intrinsics

### Implementation Table:

| Function             | No SIMD | SSE2 | AVX | AVX2 | AVX-512 |
|:---------------------|:-------:|:----:|:---:|:----:|:-------:|
| snpgdsDiss [»](http://zhengxwen.github.io/SNPRelate/release/help/snpgdsDiss.html)           | X |
| snpgdsEIGMIX [»](http://zhengxwen.github.io/SNPRelate/release/help/snpgdsEIGMIX.html)        | X | X | X |
| snpgdsGRM [»](http://zhengxwen.github.io/SNPRelate/release/help/snpgdsGRM.html)           | X | X | X |
| snpgdsIBDKING [»](http://zhengxwen.github.io/SNPRelate/release/help/snpgdsIBDKING.html)       | X | X |   | X |
| snpgdsIBDMoM [»](http://zhengxwen.github.io/SNPRelate/release/help/snpgdsIBDMoM.html)        | X |
| snpgdsIBS [»](http://zhengxwen.github.io/SNPRelate/release/help/snpgdsIBS.html)           | X | X |
| snpgdsIBSNum [»](http://zhengxwen.github.io/SNPRelate/release/help/snpgdsIBSNum.html)        | X | X |
| snpgdsIndivBeta [»](http://zhengxwen.github.io/SNPRelate/release/help/snpgdsIndivBeta.html)     | X | X | P | X |
| snpgdsPCA [»](http://zhengxwen.github.io/SNPRelate/release/help/snpgdsPCA.html)           | X | X | X |
| snpgdsPCACorr [»](http://zhengxwen.github.io/SNPRelate/release/help/snpgdsPCACorr.html)       | X |
| snpgdsPCASampLoading [»](http://zhengxwen.github.io/SNPRelate/release/help/snpgdsPCASampLoading.html) | X |
| snpgdsPCASNPLoading [»](http://zhengxwen.github.io/SNPRelate/release/help/snpgdsPCASNPLoading.html) | X |
| [...](http://zhengxwen.github.io/SNPRelate/release/help/00Index.html) |

`X: fully supported;  .: partially supported;  P: POPCNT instruction.`


### Install the package from the source code with the support of Intel SIMD Intrinsics:

You have to customize the package compilation, see: [CRAN: Customizing-package-compilation](http://cran.r-project.org/doc/manuals/r-release/R-admin.html#Customizing-package-compilation)

Change `~/.R/Makevars` to, assuming GNU Compilers (gcc/g++) or Clang compiler (clang++) are installed:
```sh
## for C code
CFLAGS=-g -O3 -march=native -mtune=native
## for C++ code
CXXFLAGS=-g -O3 -march=native -mtune=native
```


## Implementation with OpenCL

In progress ...
