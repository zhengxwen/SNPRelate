\name{snpgdsPCA}
\alias{snpgdsPCA}
\title{
	Principal Component Analysis (PCA) on SNP genotype data
}
\description{
	To calculate the eigenvectors and eigenvalues for principal component analysis in GWAS.
}
\usage{
snpgdsPCA(gdsobj, sample.id = NULL, snp.id = NULL, autosome.only = TRUE,
	remove.monosnp = TRUE, maf = NaN, missing.rate = NaN, eigen.cnt = 32,
	num.thread = 1, bayesian = FALSE, need.genmat = FALSE, genmat.only = FALSE,
	verbose = TRUE)
}
\arguments{
	\item{gdsobj}{a GDS file object (\code{\link[gdsfmt]{gds.class}})}
	\item{sample.id}{a vector of sample id specifying selected samples; if NULL, all samples are used}
	\item{snp.id}{a vector of snp id specifying selected SNPs; if NULL, all SNPs are used}
	\item{autosome.only}{if TRUE, use autosomal SNPs only}
	\item{remove.monosnp}{if TRUE, remove monomorphic SNPs}
	\item{maf}{to use the SNPs with ">= maf" only; if NaN, no MAF threshold}
	\item{missing.rate}{to use the SNPs with "<= missing.rate" only; if NaN, no missing threshold}
	\item{eigen.cnt}{output the number of eigenvectors; if eigen.cnt <= 0, then return all eigenvectors}
	\item{num.thread}{the number of CPU cores used}
	\item{bayesian}{if TRUE, use bayesian normalization}
	\item{need.genmat}{if TRUE, return the genetic covariance matrix}
	\item{genmat.only}{return the genetic covariance matrix only, do not compute the eigenvalues and eigenvectors}
	\item{verbose}{if TRUE, show information}
}
\details{
	The minor allele frequency and missing rate for each SNP passed in \code{snp.id} are
calculated over all the samples in \code{sample.id}.
}
\value{
	Return a \code{snpgdsPCAClass} object, and it is a list:
	\item{sample.id}{the sample ids used in the analysis}
	\item{snp.id}{the SNP ids used in the analysis}
	\item{eigenval}{eigenvalues}
	\item{eigenvect}{eigenvactors, "# of samples" x "eigen.cnt"}
	\item{TraceXTX}{the trace of the genetic covariance matrix}
	\item{Bayesian}{whether use bayerisan normalization}
	\item{genmat}{the genetic covariance matrix}
}
\references{
	Patterson N, Price AL, Reich D (2006) Population structure and eigenanalysis. PLoS Genetics 2:e190.

	Price AL, Patterson NJ, Plenge RM, Weinblatt ME, Shadick NA, Reich D (2006)
		Principal components analysis corrects for stratification in genome-wide association studies.
		Nat Genet. 38, 904-909.
}
\author{Xiuwen Zheng \email{zhengx@u.washington.edu}}
\seealso{
	\code{\link{snpgdsPCACorr}}, \code{\link{snpgdsPCASampLoading}}, \code{\link{snpgdsPCASNPLoading}}
}

\examples{
# open an example dataset (HapMap)
genofile <- openfn.gds(snpgdsExampleFileName())

# get population information
#   or pop_code <- scan("pop.txt", what=character()), if it is stored in a text file "pop.txt"
pop_code <- read.gdsn(index.gdsn(genofile, "sample.annot/pop.group"))

RV <- snpgdsPCA(genofile)
plot(RV$eigenvect[,2], RV$eigenvect[,1], col=as.integer(factor(pop_code)),
	xlab="PC 2", ylab="PC 1")
legend("topleft", legend=levels(factor(pop_code)), pch="o", col=1:4)

# close the genotype file
closefn.gds(genofile)
}

\keyword{multicore}
\keyword{gds}
\keyword{PCA}
\keyword{GWAS}
