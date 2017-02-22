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

## Bioconductor:

Release Version: v1.8.0

[http://www.bioconductor.org/packages/release/bioc/html/SNPRelate.html](http://www.bioconductor.org/packages/release/bioc/html/SNPRelate.html)

Development Version: v1.9.5

[http://www.bioconductor.org/packages/devel/bioc/html/SNPRelate.html](http://www.bioconductor.org/packages/devel/bioc/html/SNPRelate.html)


## News

* Supports the [SeqArray](http://bioconductor.org/packages/release/bioc/html/SeqArray.html) GDS format, see [the vignette](http://www.bioconductor.org/packages/release/bioc/vignettes/SeqArray/inst/doc/R_Integration.html#integration-with-snprelate).


## Tutorials

[http://corearray.sourceforge.net/tutorials/SNPRelate](http://corearray.sourceforge.net/tutorials/SNPRelate)

[http://www.bioconductor.org/packages/devel/bioc/vignettes/SNPRelate/inst/doc/SNPRelateTutorial.html](http://www.bioconductor.org/packages/devel/bioc/vignettes/SNPRelate/inst/doc/SNPRelateTutorial.html)


## Citation

Zheng X, Levine D, Shen J, Gogarten SM, Laurie C, Weir BS (2012). A High-performance Computing Toolset for Relatedness and Principal Component Analysis of SNP Data. *Bioinformatics*. [DOI: 10.1093/bioinformatics/bts606](http://dx.doi.org/10.1093/bioinformatics/bts606).


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


* Old version (<= v0.9.19) from [R-Forge](http://R-Forge.R-project.org) repository:
```R
install.packages("gdsfmt", repos="http://R-Forge.R-project.org")
install.packages("SNPRelate", repos="http://R-Forge.R-project.org")
```



## Implementation with Intel Intrinsics

| Function             | No SIMD | SSE2 | AVX | AVX2 | AVX-512 |
|:---------------------|:-------:|:----:|:---:|:----:|:-------:|
| snpgdsDiss           | X |
| snpgdsEIGMIX         | X | X | X |
| snpgdsGRM            | X | X | X |
| snpgdsIBDKING        | X | X |   | X |
| snpgdsIBDMoM         | X |
| snpgdsIBS            | X | X |
| snpgdsIBSNum         | X | X |
| snpgdsIndivBeta      | X | X | P | X |
| snpgdsPCA            | X | X | X |
| snpgdsPCACorr        | X |
| snpgdsPCASampLoading | X |
| snpgdsPCASNPLoading  | X |

`X`: fully supported;  `.`: partially supported; `P`: POPCNT instruction.
