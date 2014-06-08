\name{snpgdsGetGeno}
\alias{snpgdsGetGeno}
\title{
	To get a genotype matrix
}
\description{
	To get a genotype matrix from a specified GDS file
}
\usage{
snpgdsGetGeno(gdsobj, sample.id=NULL, snp.id=NULL, snpfirstdim=NULL, verbose=TRUE)
}
\arguments{
	\item{gdsobj}{a GDS file object (\code{\link[gdsfmt]{gds.class}})}
	\item{sample.id}{a vector of sample id specifying selected samples; if NULL, all samples are used}
	\item{snp.id}{a vector of snp id specifying selected SNPs; if NULL, all SNPs are used}
	\item{snpfirstdim}{if TRUE, genotypes are stored in the individual-major mode,
		(i.e, list all SNPs for the first individual, and then list all SNPs for
		the second individual, etc); if NULL, determine automatically}
	\item{verbose}{if TRUE, show information}
}
\value{
	Return an integer matrix.
}
\author{Xiuwen Zheng}
\seealso{
	\code{\link{snpgdsCreateGeno}}, \code{\link{snpgdsCreateGenoSet}},
	\code{\link{snpgdsCombineGeno}}
}

\examples{
# open an example dataset (HapMap)
genofile <- openfn.gds(snpgdsExampleFileName())

snpset <- sample(read.gdsn(index.gdsn(genofile, "snp.id")), 1000)
mat <- snpgdsGetGeno(genofile, snp.id=snpset)
mat[!is.element(mat, c(0,1,2))] <- NA
dim(mat)

# close the file
closefn.gds(genofile)
}

\keyword{gds}
\keyword{GWAS}
\keyword{LD}
