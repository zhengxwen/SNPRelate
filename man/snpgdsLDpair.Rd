\name{snpgdsLDpair}
\alias{snpgdsLDpair}
\title{
	Linkage Disequilibrium (LD)
}
\description{
	Return a LD value between snp1 and snp2.
}
\usage{
snpgdsLDpair(snp1, snp2, method = c("composite", "r", "dprime", "corr"))
}
\arguments{
	\item{snp1}{a vector of SNP genotypes (0 -- BB, 1 -- AB, 2 -- AA)}
	\item{snp2}{a vector of SNP genotypes (0 -- BB, 1 -- AB, 2 -- AA)}
	\item{method}{"composite", "r", "dprime", "corr", see details}
}
\details{
	Four methods can be used to calculate linkage disequilibrium values:
"composite" for LD composite measure, "r" for r square, "dprime" for D', and "corr" for
correlation coefficient. The method "corr" is equivalent to "composite", when SNP genotypes
are coded as: 0 -- BB, 1 -- AB, 2 -- AA.
}
\value{
	A measure of linkage disequilibrium.
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
	\code{\link{snpgdsLDMat}}, \code{\link{snpgdsLDpruning}}
}

\examples{
# open an example dataset (HapMap)
genofile <- openfn.gds(snpgdsExampleFileName())

snp1 <- read.gdsn(index.gdsn(genofile, "genotype"), start=c(1,1), count=c(1,-1))
snp2 <- read.gdsn(index.gdsn(genofile, "genotype"), start=c(2,1), count=c(1,-1))

snpgdsLDpair(snp1, snp2)

# close the genotype file
closefn.gds(genofile)
}

\keyword{gds}
\keyword{GWAS}
\keyword{LD}
