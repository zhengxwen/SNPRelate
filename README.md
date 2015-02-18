SNPRelate: Parallel computing toolset for relatedness and principal component analysis of SNP data
====

[![Build Status](https://travis-ci.org/zhengxwen/SNPRelate.png)](https://travis-ci.org/zhengxwen/SNPRelate)


## Features

Genome-wide association studies are widely used to investigate the genetic basis of diseases and traits, but they pose many computational challenges. We developed SNPRelate (R package for multi-core symmetric multiprocessing computer architectures) to accelerate two key computations on SNP data: principal component analysis (PCA) and relatedness analysis using identity-by-descent measures. The kernels of our algorithms are written in C/C++ and highly optimized.


## Bioconductor:

Release Version: v1.0.1

[http://www.bioconductor.org/packages/release/bioc/html/SNPRelate.html](http://www.bioconductor.org/packages/release/bioc/html/SNPRelate.html)

Development Version: v1.1.9

[http://www.bioconductor.org/packages/devel/bioc/html/SNPRelate.html](http://www.bioconductor.org/packages/devel/bioc/html/SNPRelate.html)


## Tutorials

[http://corearray.sourceforge.net/tutorials/SNPRelate](http://corearray.sourceforge.net/tutorials/SNPRelate)

[http://www.bioconductor.org/packages/devel/bioc/vignettes/SNPRelate/inst/doc/SNPRelateTutorial.html](http://www.bioconductor.org/packages/devel/bioc/vignettes/SNPRelate/inst/doc/SNPRelateTutorial.html)


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


* Old version (<=v0.9.19) from [R-Forge](http://R-Forge.R-project.org) repository:
```R
install.packages("gdsfmt", repos="http://R-Forge.R-project.org")
install.packages("SNPRelate", repos="http://R-Forge.R-project.org")
```
