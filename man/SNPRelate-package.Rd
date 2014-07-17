\name{SNPRelate-package}
\alias{SNPRelate-package}
\alias{SNPRelate}
\docType{package}
\title{Parallel Computing Toolset for Genome-Wide Association Studies}
\description{
    Genome-wide association studies are widely used to investigate
the genetic basis of diseases and traits, but they pose many computational
challenges. We developed SNPRelate (R package for multi-core symmetric
multiprocessing computer architectures) to accelerate two key computations
on SNP data: principal component analysis (PCA) and relatedness analysis
using identity-by-descent measures. The kernels of our algorithms are
written in C/C++ and highly optimized.
}
\details{
    \tabular{ll}{
        Package: \tab SNPRelate \cr
        Type: \tab Package \cr
        License: \tab GPL version 3 \cr
        Depends: \tab gdsfmt (>= 1.0.4) \cr
    }
    The genotypes stored in GDS format can be analyzed by the R functions in
    SNPRelate, which utilize the multi-core feature of machine for a single
    computer.

    Webpage: \url{http://github.com/zhengxwen/SNPRelate},
        \url{http://corearray.sourceforge.net/}

    Tutorial: \url{http://corearray.sourceforge.net/materials/SNPRelate/}
}
\references{
    Zheng X, Levine D, Shen J, Gogarten SM, Laurie C, Weir BS.
    A High-performance Computing Toolset for Relatedness and Principal
    Component Analysis of SNP Data.
    Bioinformatics (2012); doi: 10.1093/bioinformatics/bts610
}
\author{
    Xiuwen Zheng \email{zhengxwen@gmail.com}
}

\examples{
####################################################################
# Convert the PLINK BED file to the GDS file
#

# PLINK BED files
bed.fn <- system.file("extdata", "plinkhapmap.bed", package="SNPRelate")
bim.fn <- system.file("extdata", "plinkhapmap.bim", package="SNPRelate")
fam.fn <- system.file("extdata", "plinkhapmap.fam", package="SNPRelate")
# convert
snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, "HapMap.gds")


####################################################################
# Principal Component Analysis
#

# open
genofile <- snpgdsOpen("HapMap.gds")

RV <- snpgdsPCA(genofile)
plot(RV$eigenvect[,2], RV$eigenvect[,1], xlab="PC 2", ylab="PC 1",
    col=rgb(0,0,150, 50, maxColorValue=255), pch=19)

# close the file
snpgdsClose(genofile)


####################################################################
# Identity-By-Descent (IBD) Analysis
#

# open
genofile <- snpgdsOpen(snpgdsExampleFileName())

RV <- snpgdsIBDMoM(genofile)
flag <- lower.tri(RV$k0)
plot(RV$k0[flag], RV$k1[flag], xlab="k0", ylab="k1",
    col=rgb(0,0,150, 50, maxColorValue=255), pch=19)
abline(1, -1, col="red", lty=4)

# close the file
snpgdsClose(genofile)


####################################################################
# Identity-By-State (IBS) Analysis
#

# open
genofile <- snpgdsOpen(snpgdsExampleFileName())

RV <- snpgdsIBS(genofile)
m <- 1 - RV$ibs
colnames(m) <- rownames(m) <- RV$sample.id
GeneticDistance <- as.dist(m[1:45, 1:45])
HC <- hclust(GeneticDistance, "ave")
plot(HC)

# close the file
snpgdsClose(genofile)


####################################################################
# Linkage Disequilibrium (LD) Analysis
#

# open an example dataset (HapMap)
genofile <- snpgdsOpen(snpgdsExampleFileName())

snpset <- read.gdsn(index.gdsn(genofile, "snp.id"))[1:200]
L1 <- snpgdsLDMat(genofile, snp.id=snpset, method="composite", slide=-1)

# plot
image(abs(L1$LD), col=terrain.colors(64))

# close the file
snpgdsClose(genofile)
}

\keyword{GDS}
\keyword{GWAS}
