\name{SNPRelate-package}
\alias{SNPRelate-package}
\alias{SNPRelate}
\docType{package}
\title{Parallel Computing Toolset for genome-wide association studies}
\description{
	A High-performance Computing Package for Relatedness, Linkage Disequilibrium and
Principal Component Analysis in Genome-wide Association Studies
}
\details{
\tabular{ll}{
Package: \tab SNPRelate\cr
Type: \tab Package\cr
Version: \tab 0.9.1\cr
Date: \tab 2011-08-11\cr
License: \tab LGPL version 3\cr
Depends: \tab gdsfmt (>= 0.9.5)\cr
}
	The genotypes stored in GDS format can be analyzed by the R functions in SNPRelate, which
	utilize the multi-core feature of machine for a single computer.
}
\author{Xiuwen Zheng \email{zhengx@u.washington.edu}}
\examples{
####################################################################
# Convert the PLINK BED file to the GDS file
#

# PLINK BED files
path <- paste(.path.package("SNPRelate"), "/extdata", sep="")
bed.fn <- paste(path, "/plinkhapmap.bed", sep="")
bim.fn <- paste(path, "/plinkhapmap.bim", sep="")
fam.fn <- paste(path, "/plinkhapmap.fam", sep="")
# convert
snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, "HapMap.gds")


####################################################################
# Principal Component Analysis
#

# open
genofile <- openfn.gds("HapMap.gds")

RV <- snpgdsPCA(genofile, num.thread=2)
plot(RV$eigenvect[,2], RV$eigenvect[,1], xlab="PC 2", ylab="PC 1",
	col=rgb(0,0,150, 50, maxColorValue=255), pch=19)

# close
closefn.gds(genofile)


####################################################################
# Identity-By-Descent (IBD) Analysis
#

# open
genofile <- openfn.gds(snpgdsExampleFileName())

RV <- snpgdsIBDMoM(genofile, num.thread=2)
flag <- lower.tri(RV$k0)
plot(RV$k0[flag], RV$k1[flag], xlab="k0", ylab="k1",
	col=rgb(0,0,150, 50, maxColorValue=255), pch=19)
abline(1, -1, col="red", lty=4)

# close
closefn.gds(genofile)


####################################################################
# Identity-By-State (IBS) Analysis
#

# open
genofile <- openfn.gds(snpgdsExampleFileName())

RV <- snpgdsIBS(genofile, num.thread=2)
m <- 1 - RV$ibs
colnames(m) <- rownames(m) <- RV$sample.id
GeneticDistance <- as.dist(m[1:45, 1:45])
HC <- hclust(GeneticDistance, "ave")
plot(HC)

# close
closefn.gds(genofile)


####################################################################
# Linkage Disequilibrium (LD) Analysis
#

# open an example dataset (HapMap)
genofile <- openfn.gds(snpgdsExampleFileName())

snpset <- read.gdsn(index.gdsn(genofile, "snp.id"))[1:200]
L1 <- snpgdsLDMat(genofile, snp.id=snpset, method="composite")

# plot
image(abs(L1$LD), col=terrain.colors(64))

# close the genotype file
closefn.gds(genofile)
}

\keyword{gds}
\keyword{GWAS}
\keyword{genetics}
