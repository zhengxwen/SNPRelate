\name{snpgdsIBDKING}
\alias{snpgdsIBDKING}
\title{
	KING method of moment for the identity-by-descent (IBD) analysis
}
\description{
	Calculate IBD coefficients by KING method of moment.
}
\usage{
snpgdsIBDKING(gdsobj, sample.id=NULL, snp.id=NULL, autosome.only=TRUE,
	remove.monosnp=TRUE, maf=NaN, missing.rate=NaN,
	method = c("robust.across.family", "robust.within.family", "homo"),
	kinship=TRUE, num.thread=1, verbose=TRUE)
}
\arguments{
	\item{gdsobj}{a GDS file object (\code{\link[gdsfmt]{gds.class}})}
	\item{sample.id}{a vector of sample id specifying selected samples; if NULL,
		all samples are used}
	\item{snp.id}{a vector of snp id specifying selected SNPs; if NULL, all SNPs are used}
	\item{autosome.only}{if TRUE, use autosomal SNPs only}
	\item{remove.monosnp}{if TRUE, remove monomorphic SNPs}
	\item{maf}{to use the SNPs with ">= maf" only; if NaN, no MAF threshold}
	\item{missing.rate}{to use the SNPs with "<= missing.rate" only; if NaN, no
		missing threshold}
	\item{method}{should be "robust.across.family" (default),
		"robust.within.family" or "homo", see details}
	\item{kinship}{if \code{TRUE}, output the estimated kinship coefficients}
	\item{num.thread}{the number of CPU cores used}
	\item{verbose}{if TRUE, show information}
}
\details{
	KING IBD estimator is a moment estimator, and it is computationally efficient
relative to MLE method. The method includes \code{"robust.across.family"} -- robust
relationship inference across family in the presence of population substructure
(according to Equation 11 in Manichaikul et al. 2010); \code{"robust.within.family"} --
robust relationship inference within family in the presence of population substructure
(according to Equation 9 in Manichaikul et al. 2010); \code{"homo"} -- relationship
inference in a homogeneous population (according to Equation 5 in Manichaikul
et al. 2010). The method of \code{"robust.across.family"} always returns lower or equal
kinship coefficients than what \code{"robust.within.family"} does. The probability of
no sharing IBD (k0) is calculated using Equation 2 in Manichaikul et al. 2010.
\code{\link{snpgdsIBSNum}} could be used to calculate the proportions of SNPs sharing
zero IBS.

	The minor allele frequency and missing rate for each SNP passed in \code{snp.id}
are calculated over all the samples in \code{sample.id}.
}
\value{
	Return a list:
	\item{sample.id}{the sample ids used in the analysis}
	\item{snp.id}{the SNP ids used in the analysis}
	\item{k0}{IBD coefficient, the probability of sharing ZERO IBD}
	\item{k1}{IBD coefficient, the probability of sharing ONE IBD}
	\item{kinship}{the estimated kinship coefficients, if the parameter \code{kinship=TRUE}}
}
\references{
	Manichaikul A, Mychaleckyj JC, Rich SS, Daly K, Sale M, Chen WM.
	Robust relationship inference in genome-wide association studies.
	Bioinformatics. 2010 Nov 15;26(22):2867-73. doi: 10.1093/bioinformatics/btq559.
	Epub 2010 Oct 5.
}
\author{Xiuwen Zheng \email{zhengx@u.washington.edu}}
\seealso{
	\code{\link{snpgdsIBDMLE}}, \code{\link{snpgdsIBDMoM}}
}

\examples{
# open an example dataset (HapMap)
genofile <- openfn.gds(snpgdsExampleFileName())

# CEU population
samp.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
CEU.id <- samp.id[read.gdsn(index.gdsn(genofile, "sample.annot/pop.group"))=="CEU"]


####  KING-homo: relationship inference in a homogeneous population

ibd.homo <- snpgdsIBDKING(genofile, sample.id=CEU.id, method="homo")
names(ibd.homo)
# "sample.id" "snp.id"    "afreq"     "k0"        "k1"        "kinship"

flag <- lower.tri(ibd.homo$k0)
plot(ibd.homo$k0[flag], ibd.homo$kinship[flag],
	xlab="Pr(IBD=0)", ylab="Estimated Kinship Coefficient (KING-homo)")


####  KING-robust: relationship inference in the presence of population stratification
####               robust relationship inference across family

ibd.robust <- snpgdsIBDKING(genofile, sample.id=CEU.id)
names(ibd.robust)
# "sample.id" "snp.id"    "afreq"     "k0"        "k1"        "kinship"

ibs <- snpgdsIBSNum(genofile, sample.id=CEU.id)
# proportion of SNPs with zero IBS
ibs0 <- ibs$ibs0 / (ibs$ibs0 + ibs$ibs1 + ibs$ibs2)

flag <- lower.tri(ibd.robust$kinship)
plot(ibs0[flag], ibd.robust$kinship[flag],
	xlab="Proportion of Zero IBS", ylab="Estimated Kinship Coefficient (KING-robust)")


# close the genotype file
closefn.gds(genofile)
}

\keyword{multicore}
\keyword{gds}
\keyword{GWAS}
\keyword{IBD}
