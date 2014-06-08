\name{snpgdsSelectSNP}
\alias{snpgdsSelectSNP}
\title{
	SNP selection
}
\description{
	Create a list of candidate SNPs based on specified criteria
}
\usage{
snpgdsSelectSNP(gdsobj, sample.id=NULL, snp.id=NULL, autosome.only=TRUE,
	remove.monosnp=TRUE, maf=NaN, missing.rate=NaN, verbose=TRUE)
}
\arguments{
	\item{gdsobj}{a GDS file object (\code{\link[gdsfmt]{gds.class}})}
	\item{sample.id}{a vector of sample id specifying selected samples; if NULL, all samples are used}
	\item{snp.id}{a vector of snp id specifying selected SNPs; if NULL, all SNPs are used}
	\item{autosome.only}{if TRUE, use autosomal SNPs only}
	\item{remove.monosnp}{if TRUE, remove monomorphic SNPs}
	\item{maf}{to use the SNPs with ">= maf" only; if NaN, no any MAF threshold}
	\item{missing.rate}{to use the SNPs with "<= missing.rate" only; if NaN, no any missing threshold}
	\item{verbose}{if TRUE, show information}
}
\value{
	Return a list of snp ids.
}
\author{Xiuwen Zheng}
\seealso{
	\code{\link{snpgdsSampMissrate}}, \code{\link{snpgdsSNPRateFreq}},
	\code{\link{snpgdsLDpruning}}
}

\examples{
# open an example dataset (HapMap)
genofile <- openfn.gds(snpgdsExampleFileName())

snpset <- snpgdsSelectSNP(genofile, maf=0.05, missing.rate=0.95)

# close the genotype file
closefn.gds(genofile)
}

\keyword{gds}
\keyword{GWAS}
