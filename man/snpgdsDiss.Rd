\name{snpgdsDiss}
\alias{snpgdsDiss}
\title{
	Individual dissimilarity analysis
}
\description{
	Calculate the individual dissimilarities for each pair of individuals
}
\usage{
snpgdsDiss(gdsobj, sample.id=NULL, snp.id=NULL, autosome.only=TRUE,
	remove.monosnp=TRUE, maf=NaN, missing.rate=NaN, num.thread=1, verbose=TRUE)
}
\arguments{
	\item{gdsobj}{a GDS file object (\code{\link[gdsfmt]{gds.class}})}
	\item{sample.id}{a vector of sample id specifying selected samples; if NULL, all samples are used}
	\item{snp.id}{a vector of snp id specifying selected SNPs; if NULL, all SNPs are used}
	\item{autosome.only}{if TRUE, use autosomal SNPs only}
	\item{remove.monosnp}{if TRUE, remove monomorphic SNPs}
	\item{maf}{to use the SNPs with ">= maf" only; if NaN, no MAF threshold}
	\item{missing.rate}{to use the SNPs with "<= missing.rate" only; if NaN, no missing threshold}
	\item{num.thread}{the number of CPU cores used}
	\item{verbose}{if TRUE, show information}
}
\details{
	The minor allele frequency and missing rate for each SNP passed in \code{snp.id} are
calculated over all the samples in \code{sample.id}.

	The details will be described in future.
}
\value{
	Return a class "snpgdsDissClass":
	\item{sample.id}{the sample ids used in the analysis}
	\item{snp.id}{the SNP ids used in the analysis}
	\item{diss}{a matrix of individual dissimilarity}
}
\author{Xiuwen Zheng}
\seealso{
	\code{\link{snpgdsHCluster}}
}

\keyword{multicore}
\keyword{gds}
\keyword{GWAS}
\keyword{IBS}
