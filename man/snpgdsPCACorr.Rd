\name{snpgdsPCACorr}
\alias{snpgdsPCACorr}
\title{
	SNP correlation in principal component analysis
}
\description{
	To calculate the SNP correlations between eigenvactors and SNP genotypes
}
\usage{
snpgdsPCACorr(pcaobj, gdsobj, snp.id=NULL, eig.which=NULL, num.thread=1, verbose=TRUE)
}
\arguments{
	\item{pcaobj}{the \code{snpgdsPCAClass} object returned from the function \link{snpgdsPCA}}
	\item{gdsobj}{a GDS file object (\code{\link[gdsfmt]{gds.class}})}
	\item{snp.id}{a vector of snp id specifying selected SNPs; if NULL, all SNPs are used}
	\item{eig.which}{a vector of integers, to specify which eigenvectors to be used}
	\item{num.thread}{the number of CPU cores used}
	\item{verbose}{if TRUE, show information}
}
\value{
	Return a list:
	\item{sample.id}{the sample ids used in the analysis}
	\item{snp.id}{the SNP ids used in the analysis}
	\item{snpcorr}{a matrix of correlation coefficients, "# of eigenvectors" x "# of SNPs"}
}
\references{
	Patterson N, Price AL, Reich D (2006) Population structure and eigenanalysis. PLoS Genetics 2:e190.
}
\author{Xiuwen Zheng}
\seealso{
	\code{\link{snpgdsPCA}}, \code{\link{snpgdsPCASampLoading}}, \code{\link{snpgdsPCASNPLoading}}
}

\examples{
# open an example dataset (HapMap)
genofile <- openfn.gds(snpgdsExampleFileName())
# get chromosome index
chr <- read.gdsn(index.gdsn(genofile, "snp.chromosome"))

RV <- snpgdsPCA(genofile)
CORR <- snpgdsPCACorr(RV, genofile, eig.which=1:4)
plot(abs(CORR$snpcorr[3,]), xlab="SNP Index", ylab="PC 3", col=chr)

# close the genotype file
closefn.gds(genofile)
}

\keyword{gds}
\keyword{PCA}
\keyword{GWAS}
