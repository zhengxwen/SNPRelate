\name{snpgdsSampMissrate}
\alias{snpgdsSampMissrate}
\title{
	Missing Rate of Samples
}
\description{
	Return the missing fraction for each sample
}
\usage{
snpgdsSampMissrate(gdsobj, sample.id = NULL, snp.id = NULL)
}
\arguments{
	\item{gdsobj}{the \code{\link[gdsfmt]{gds.class}} object in the \link{gdsfmt} package}
	\item{sample.id}{a vector of sample id specifying selected samples; if NULL, all samples are used}
	\item{snp.id}{a vector of snp id specifying selected SNPs; if NULL, all SNPs are used}
}
\value{
	a vector of numeric values
}
\author{Xiuwen Zheng \email{zhengx@u.washington.edu}}
\seealso{
	\code{\link{snpgdsSNPRateFreq}}
}

\examples{
# open an example dataset (HapMap)
genofile <- openfn.gds(snpgdsExampleFileName())

RV <- snpgdsSampMissrate(genofile)
summary(RV)

# close the genotype file
closefn.gds(genofile)
}

\keyword{gds}
\keyword{GWAS}
