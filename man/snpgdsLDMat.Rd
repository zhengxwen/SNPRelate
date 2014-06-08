\name{snpgdsLDMat}
\alias{snpgdsLDMat}
\title{
	Linkage Disequilibrium (LD) analysis
}
\description{
	Return a LD matrix for SNP pairs.
}
\usage{
snpgdsLDMat(gdsobj, sample.id = NULL, snp.id = NULL,
	method = c("composite", "r", "dprime", "corr"), num.thread = 1, verbose = TRUE)
}
\arguments{
	\item{gdsobj}{the \code{gdsclass} object in the \link{gdsfmt} package}
	\item{sample.id}{a vector of sample id specifying selected samples; if NULL, all samples are used}
	\item{snp.id}{a vector of snp id specifying selected SNPs; if NULL, all SNPs are used}
	\item{method}{"composite", "r", "dprime", "corr", see details}
	\item{num.thread}{the number of CPU cores used}
	\item{verbose}{if TRUE, show information}
}
\details{
	Four methods can be used to calculate linkage disequilibrium values:
	"composite" for LD composite measure, "r" for r square, "dprime" for D', and "corr" for
	correlation coefficient. The method "corr" is equivalent to "composite", when SNP genotypes
	are coded as: 0 -- BB, 1 -- AB, 2 -- AA.
}
\value{
	Return a list:
	\item{sample.id}{the sample ids used in the analysis}
	\item{snp.id}{the SNP ids used in the analysis}
	\item{LD}{a matrix of LD values}
}
\references{
	Weir B: Inferences about linkage disequilibrium. Biometrics 1979; 35: 235-254.

	Weir B: Genetic Data Analysis II. Sunderland, MA: Sinauer Associates, 1996.

	Weir BS, Cockerham CC: Complete characterization of disequilibrium at two loci;
		in Feldman MW (ed): Mathematical Evolutionary Theory.
		Princeton, NJ: Princeton University Press, 1989.
}
\author{Xiuwen Zheng \email{zhengx@u.washington.edu}}
\seealso{
	\code{\link{snpgdsLDpair}}, \code{\link{snpgdsLDpruning}}
}

\examples{
# open an example dataset (HapMap)
genofile <- openfn.gds(snpgdsExampleFileName())

snpset <- read.gdsn(index.gdsn(genofile, "snp.id"))[1:200]
L1 <- snpgdsLDMat(genofile, snp.id=snpset, method="composite")
L2 <- snpgdsLDMat(genofile, snp.id=snpset, method="r")
L3 <- snpgdsLDMat(genofile, snp.id=snpset, method="dprime")
L4 <- snpgdsLDMat(genofile, snp.id=snpset, method="corr")

# plot
image(abs(L1$LD))

# close the genotype file
closefn.gds(genofile)
}

\keyword{gds}
\keyword{GWAS}
\keyword{LD}
