\name{snpgdsSNPRateFreq}
\alias{snpgdsSNPRateFreq}
\title{
	Allele Frequency, Minor Allele Frequency, Missing Rate of SNPs
}
\description{
	Calculate the allele frequency, minor allele frequency and missing rate per SNP
}
\usage{
snpgdsSNPRateFreq(gdsobj, sample.id=NULL, snp.id=NULL)
}
\arguments{
	\item{gdsobj}{a GDS file object (\code{\link[gdsfmt]{gds.class}})}
	\item{sample.id}{a vector of sample id specifying selected samples; if NULL, all samples are used}
	\item{snp.id}{a vector of snp id specifying selected SNPs; if NULL, all SNPs are used}
}
\value{
	Return a list:
	\item{AlleleFreq}{the allele frequencies}
	\item{MinorFreq}{the minor allele frequencies}
	\item{MissingRate}{the missing rates}
}
\author{Xiuwen Zheng}
\seealso{
	\code{\link{snpgdsSampMissrate}}
}

\examples{
# open an example dataset (HapMap)
genofile <- openfn.gds(snpgdsExampleFileName())

RV <- snpgdsSNPRateFreq(genofile)
hist(RV$AlleleFreq, breaks=128)
summary(RV$MissingRate)

# close the genotype file
closefn.gds(genofile)
}

\keyword{gds}
\keyword{GWAS}
