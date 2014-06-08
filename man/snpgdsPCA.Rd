\name{snpgdsPCA}
\alias{snpgdsPCA}
\title{
	Principal Component Analysis (PCA) on SNP genotype data
}
\description{
	To calculate the eigenvectors and eigenvalues for principal component analysis
in GWAS.
}
\usage{
snpgdsPCA(gdsobj, sample.id = NULL, snp.id = NULL, autosome.only = TRUE,
	remove.monosnp = TRUE, maf = NaN, missing.rate = NaN, eigen.cnt = 32,
	num.thread = 1, bayesian = FALSE, need.genmat = FALSE,
	genmat.only = FALSE, verbose = TRUE)
}
\arguments{
	\item{gdsobj}{a GDS file object (\code{\link[gdsfmt]{gds.class}})}
	\item{sample.id}{a vector of sample id specifying selected samples;
		if NULL, all samples are used}
	\item{snp.id}{a vector of snp id specifying selected SNPs; if NULL,
		all SNPs are used}
	\item{autosome.only}{if TRUE, use autosomal SNPs only}
	\item{remove.monosnp}{if TRUE, remove monomorphic SNPs}
	\item{maf}{to use the SNPs with ">= maf" only; if NaN, no MAF threshold}
	\item{missing.rate}{to use the SNPs with "<= missing.rate" only; if NaN,
		no missing threshold}
	\item{eigen.cnt}{output the number of eigenvectors; if eigen.cnt <= 0, then
		return all eigenvectors}
	\item{num.thread}{the number of CPU cores used}
	\item{bayesian}{if TRUE, use bayesian normalization}
	\item{need.genmat}{if TRUE, return the genetic covariance matrix}
	\item{genmat.only}{return the genetic covariance matrix only, do not compute
		the eigenvalues and eigenvectors}
	\item{verbose}{if TRUE, show information}
}
\details{
	The minor allele frequency and missing rate for each SNP passed in \code{snp.id}
are calculated over all the samples in \code{sample.id}.
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
	Patterson N, Price AL, Reich D (2006) Population structure and eigenanalysis.
		PLoS Genetics 2:e190.

	Price AL, Patterson NJ, Plenge RM, Weinblatt ME, Shadick NA, Reich D (2006)
		Principal components analysis corrects for stratification in genome-wide
		association studies. Nat Genet. 38, 904-909.
}
\author{Xiuwen Zheng}
\seealso{
	\code{\link{snpgdsPCACorr}}, \code{\link{snpgdsPCASampLoading}},
	\code{\link{snpgdsPCASNPLoading}}
}

\examples{
# open an example dataset (HapMap)
genofile <- openfn.gds(snpgdsExampleFileName())

# get population information
#   or pop_code <- scan("pop.txt", what=character()), if it is stored in a text file "pop.txt"
pop_code <- read.gdsn(index.gdsn(genofile, "sample.annot/pop.group"))

# get sample id
samp.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

# run PCA
RV <- snpgdsPCA(genofile)

# eigenvalues
RV$eigenval

# variance proportion
head(RV$eigenval / sum(RV$eigenval))
# [1] 0.122763042 0.058776373 0.010073833 0.009404113 0.008397200 0.007333046

# make a data.frame
tab <- data.frame(sample.id = samp.id, pop = factor(pop_code),
	EV1 = RV$eigenvect[,1],    # the first eigenvector
	EV2 = RV$eigenvect[,2],    # the second eigenvector
	stringsAsFactors = FALSE)
head(tab)
#   sample.id pop         EV1         EV2
# 1   NA19152 YRI  0.08405438  0.01243429
# 2   NA19139 YRI  0.08363225  0.01099277
# 3   NA18912 YRI  0.08104325  0.01196670
# 4   NA19160 YRI  0.08682153  0.01465503
# 5   NA07034 CEU -0.03096588 -0.07717809
# 6   NA07055 CEU -0.03212054 -0.08159248

# draw
plot(tab$EV2, tab$EV1, col=as.integer(tab$pop),
	xlab="eigenvector 2", ylab="eigenvector 1")
legend("topleft", legend=levels(tab$pop), pch="o", col=1:4)


# close the genotype file
closefn.gds(genofile)
}

\keyword{multicore}
\keyword{gds}
\keyword{PCA}
\keyword{GWAS}
