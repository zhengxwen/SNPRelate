SNPRelate: Parallel computing toolset for relatedness and principal component analysis of SNP data
====

## Features

Genome-wide association studies are widely used to investigate the genetic basis of diseases and traits, but they pose many computational challenges. We developed SNPRelate (R package for multi-core symmetric multiprocessing computer architectures) to accelerate two key computations on SNP data: principal component analysis (PCA) and relatedness analysis using identity-by-descent measures. The kernels of our algorithms are written in C/C++ and highly optimized

## Installation

* Development version from Github:
```
library("devtools")
install_github("SNPRelate", user="zhengxwen")
```
The `install_github()` approach requires that you build from source, i.e. `make` and compilers must be installed on your system -- see the R FAQ for your operating system; you may also need to install dependencies manually.
* Nearly up-to-date development binaries from `gdsfmt` r-forge repository:
```
install.packages("SNPRelate",
   repos=c("http://SNPRelate.r-forge.r-project.org/repos",
          getOption("repos")[["CRAN"]]))
```
