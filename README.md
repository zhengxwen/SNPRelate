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

The GDS format offers the efficient operations specifically designed for integers with two bits, since a SNP could occupy only two bits. The SNP GDS format in this package is also used by the [GWASTools](http://bioconductor.org/packages/GWASTools) package with the support of S4 classes and generic functions. The extended GDS format is implemented in the [SeqArray](https://github.com/zhengxwen/SeqArray) package to support the storage of single nucleotide variation (SNV), insertion/deletion polymorphism (indel) and structural variation calls. It is strongly suggested to use [SeqArray](https://github.com/zhengxwen/SeqArray) for large-scale whole-exome and whole-genome sequencing variant data instead of [SNPRelate](https://github.com/zhengxwen/SNPRelate).


## Bioconductor

Release Version: v1.30.0

[http://www.bioconductor.org/packages/SNPRelate](http://www.bioconductor.org/packages/SNPRelate)


## News

* See [package news](NEWS).


## Tutorials

[http://www.bioconductor.org/packages/release/bioc/vignettes/SNPRelate/inst/doc/SNPRelate.html](http://www.bioconductor.org/packages/release/bioc/vignettes/SNPRelate/inst/doc/SNPRelate.html)


## Citations

Zheng X, Levine D, Shen J, Gogarten SM, Laurie C, Weir BS (2012). A High-performance Computing Toolset for Relatedness and Principal Component Analysis of SNP Data. *Bioinformatics*. [DOI: 10.1093/bioinformatics/bts606](http://dx.doi.org/10.1093/bioinformatics/bts606).

Zheng X, Gogarten S, Lawrence M, Stilp A, Conomos M, Weir BS, Laurie C, Levine D (2017). SeqArray -- A storage-efficient high-performance data format for WGS variant calls. *Bioinformatics*. [DOI: 10.1093/bioinformatics/btx145](http://dx.doi.org/10.1093/bioinformatics/btx145).


## Installation

* Bioconductor repository:
```R
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("SNPRelate")
```

* Development version from Github (for developers/testers only):
```R
library("devtools")
install_github("zhengxwen/gdsfmt")
install_github("zhengxwen/SNPRelate")
```
The `install_github()` approach requires that you build from source, i.e. `make` and compilers must be installed on your system -- see the [R FAQ](https://cran.r-project.org/faqs.html) for your operating system; you may also need to install dependencies manually.



## Implementation with Intel Intrinsics

| Functions             | No SIMD | SSE2 | AVX | AVX2 | AVX-512 |
|:----------------------|:-------:|:----:|:---:|:----:|:-------:|
| snpgdsDiss [»](https://rdrr.io/bioc/SNPRelate/man/snpgdsDiss.html)                    | X |
| snpgdsEIGMIX [»](https://rdrr.io/bioc/SNPRelate/man/snpgdsEIGMIX.html)                 | X | X | X |
| snpgdsGRM [»](https://rdrr.io/bioc/SNPRelate/man/snpgdsGRM.html)                       | X | X | X | . |
| snpgdsIBDKING [»](https://rdrr.io/bioc/SNPRelate/man/snpgdsIBDKING.html)               | X | X |   | X |
| snpgdsIBDMoM [»](https://rdrr.io/bioc/SNPRelate/man/snpgdsIBDMoM.html)                 | X |
| snpgdsIBS [»](https://rdrr.io/bioc/SNPRelate/man/snpgdsIBS.html)                       | X | X |
| snpgdsIBSNum [»](https://rdrr.io/bioc/SNPRelate/man/snpgdsIBSNum.html)                 | X | X |
| snpgdsIndivBeta [»](https://rdrr.io/bioc/SNPRelate/man/snpgdsIndivBeta.html)           | X | X | P | X |
| snpgdsPCA [»](https://rdrr.io/bioc/SNPRelate/man/snpgdsPCA.html)                       | X | X | X |
| snpgdsPCACorr [»](https://rdrr.io/bioc/SNPRelate/man/snpgdsPCACorr.html)               | X |
| snpgdsPCASampLoading [»](https://rdrr.io/bioc/SNPRelate/man/snpgdsPCASampLoading.html) | X |
| snpgdsPCASNPLoading [»](https://rdrr.io/bioc/SNPRelate/man/snpgdsPCASNPLoading.html)   | X |
| [...](http://rdrr.io/bioc/SNPRelate/man) |

`X: fully supported;  .: partially supported;  P: POPCNT instruction.`


### Install the package from the source code with the support of Intel SIMD Intrinsics:

You have to customize the package compilation, see: [CRAN: Customizing-package-compilation](https://cran.r-project.org/doc/manuals/r-release/R-admin.html#Customizing-package-compilation)

Change `~/.R/Makevars` to, assuming GNU Compilers (gcc/g++) or Clang compiler (clang++) are installed:
```sh
## for C code
CFLAGS=-g -O3 -march=native -mtune=native
## for C++ code
CXXFLAGS=-g -O3 -march=native -mtune=native
```
