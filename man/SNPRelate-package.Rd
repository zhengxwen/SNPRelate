\name{SNPRelate-package}
\alias{SNPRelate-package}
\alias{SNPRelate}
\docType{package}
\title{Parallel Computing Toolset for genome-wide association studies}
\description{
	A high-performance computing toolset for relatedness and principal component
analysis in Genome-wide Association Studies
}
\details{
	\tabular{ll}{
		Package: \tab SNPRelate \cr
		Type: \tab Package \cr
		Version: \tab 0.9.14 \cr
		Date: \tab 2013-09-04 \cr
		License: \tab GPL version 3 \cr
		Depends: \tab gdsfmt (>= 0.9.7) \cr
	}
	The genotypes stored in GDS format can be analyzed by the R functions in SNPRelate,
	which utilize the multi-core feature of machine for a single computer.

	Webpage: \url{http://corearray.sourceforge.net/}

	Tutorial: \url{http://corearray.sourceforge.net/tutorials/SNPRelate/}

	Forums:	\url{http://sourceforge.net/projects/corearray/forums}
}
\references{
	Zheng X, Levine D, Shen J, Gogarten SM, Laurie C, Weir BS.
	A High-performance Computing Toolset for Relatedness and Principal Component
	Analysis of SNP Data.
	Bioinformatics (2012); doi: 10.1093/bioinformatics/bts610
}
\author{Xiuwen Zheng \email{zhengx@u.washington.edu}}

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
genofile <- openfn.gds("HapMap.gds")

RV <- snpgdsPCA(genofile)
plot(RV$eigenvect[,2], RV$eigenvect[,1], xlab="PC 2", ylab="PC 1",
	col=rgb(0,0,150, 50, maxColorValue=255), pch=19)

# close
closefn.gds(genofile)


####################################################################
# Identity-By-Descent (IBD) Analysis
#

# open
genofile <- openfn.gds(snpgdsExampleFileName())

RV <- snpgdsIBDMoM(genofile)
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

RV <- snpgdsIBS(genofile)
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
L1 <- snpgdsLDMat(genofile, snp.id=snpset, method="composite", slide=-1)

# plot
image(abs(L1$LD), col=terrain.colors(64))

# close the genotype file
closefn.gds(genofile)
}

\keyword{multicore}
\keyword{gds}
\keyword{GWAS}
\keyword{genetics}
