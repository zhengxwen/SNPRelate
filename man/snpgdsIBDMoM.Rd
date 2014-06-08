\name{snpgdsIBDMoM}
\alias{snpgdsIBDMoM}
\title{
	PLINK method of moment (MoM) for the Identity-By-Descent (IBD) Analysis
}
\description{
	Calculate three IBD coefficients for non-inbred individual pairs by
	PLINK method of moment (MoM).
}
\usage{
snpgdsIBDMoM(gdsobj, sample.id = NULL, snp.id = NULL, autosome.only = TRUE,
	remove.monosnp = TRUE, maf = NaN, missing.rate = NaN, kinship.constraint = FALSE,
	num.thread = 1, verbose = TRUE)
}
\arguments{
	\item{gdsobj}{the \code{\link[gdsfmt]{gdsclass}} object in the \link{gdsfmt} package}
	\item{sample.id}{a vector of sample id specifying selected samples; if NULL, all samples are used}
	\item{snp.id}{a vector of snp id specifying selected SNPs; if NULL, all SNPs are used}
	\item{autosome.only}{if TRUE, use autosomal SNPs only}
	\item{remove.monosnp}{if TRUE, remove monomorphic SNPs}
	\item{maf}{to use the SNPs with ">= maf" only; if NaN, no MAF threshold}
	\item{missing.rate}{to use the SNPs with "<= missing.rate" only; if NaN, no missing threshold}
	\item{kinship.constraint}{if TRUE, constrict IBD coefficients ($k_0,k_1,k_2$) in the geneloical region ($2 k_0 k_1 >= k_2^2$)}
	\item{num.thread}{the number of CPU cores used}
	\item{verbose}{if TRUE, show information}
}
\details{
	PLINK IBD estimator is a moment estimator, and it is computationally efficient relative to
MLE method.
	The minor allele frequency and missing rate for each SNP passed in snp.id are calculated
over all the samples in sample.id.
}
\value{
	Return a list:
	\item{sample.id}{the sample ids used in the analysis}
	\item{snp.id}{the SNP ids used in the analysis}
	\item{k0}{IBD coefficient, the probability of sharing ZERO IBD}
	\item{k1}{IBD coefficient, the probability of sharing ONE IBD}
}
\references{
	Purcell S, Neale B, Todd-Brown K, Thomas L, Ferreira MAR, Bender D, Maller J, Sklar P,
	de Bakker PIW, Daly MJ & Sham PC. 2007. PLINK: a toolset for whole-genome association and
	population-based linkage analysis. American Journal of Human Genetics, 81.
}
\author{Xiuwen Zheng \email{zhengx@u.washington.edu}}
\seealso{
	\code{\link{snpgdsIBDMLE}}, \code{\link{snpgdsIBDMLELogLik}}
}

\examples{
# open an example dataset (HapMap)
genofile <- openfn.gds(snpgdsExampleFileName())

# CEU population
CEU.id <- read.gdsn(index.gdsn(genofile, "sample.id"))[
	read.gdsn(index.gdsn(genofile, c("sample.annot", "pop.group")))=="CEU"]
pibd <- snpgdsIBDMoM(genofile, sample.id=c(CEU.id))
flag <- lower.tri(pibd$k0)
plot(NaN, xlim=c(0,1), ylim=c(0,1), xlab="k0", ylab="k1")
lines(c(0,1), c(1,0), col="red", lty=3)
points(pibd$k0[flag], pibd$k1[flag])

# YRI population
YRI.id <- read.gdsn(index.gdsn(genofile, "sample.id"))[
	read.gdsn(index.gdsn(genofile, c("sample.annot", "pop.group")))=="YRI"]
pibd <- snpgdsIBDMoM(genofile, sample.id=c(YRI.id))
flag <- lower.tri(pibd$k0)
plot(NaN, xlim=c(0,1), ylim=c(0,1), xlab="k0", ylab="k1")
lines(c(0,1), c(1,0), col="red", lty=3)
points(pibd$k0[flag], pibd$k1[flag])

# close the genotype file
closefn.gds(genofile)
}

\keyword{gds}
\keyword{GWAS}
\keyword{IBD}
