SNPRelate: Parallel computing toolset for relatedness and principal component analysis of SNP data
===

Version: 0.9.20

[![Build Status](https://travis-ci.org/zhengxwen/SNPRelate.png)](https://travis-ci.org/zhengxwen/SNPRelate)


## Features

Genome-wide association studies are widely used to investigate the genetic basis of diseases and traits, but they pose many computational challenges. We developed SNPRelate (R package for multi-core symmetric multiprocessing computer architectures) to accelerate two key computations on SNP data: principal component analysis (PCA) and relatedness analysis using identity-by-descent measures. The kernels of our algorithms are written in C/C++ and highly optimized.

## Installation

* Development version from Github:
```
library("devtools")
install_github("zhengxwen/gdsfmt")
install_github("zhengxwen/SNPRelate")
```
The `install_github()` approach requires that you build from source, i.e. `make` and compilers must be installed on your system -- see the R FAQ for your operating system; you may also need to install dependencies manually.
* Nearly up-to-date development binaries from `gdsfmt` r-forge repository:
```
install.packages(c("gdsfmt", "SNPRelate"),
   repos=c("http://gdsfmt.r-forge.r-project.org/repos",
          getOption("repos")[["CRAN"]]))
```

* Install the packages (gdsfmt and SNPRelate) from the source code:
[gdsfmt](https://github.com/zhengxwen/gdsfmt/archive/v1.0.5.tar.gz)
and
[SNPRelate](https://github.com/zhengxwen/SNPRelate/archive/v0.9.20.tar.gz)
```
R CMD INSTALL gdsfmt-1.0.5.tar.gz
R CMD INSTALL SNPRelate-0.9.20.tar.gz
```
