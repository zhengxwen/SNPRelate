\name{snpgdsSampMissRate}
\alias{snpgdsSampMissRate}
\title{
    Missing Rate of Samples
}
\description{
    Return the missing fraction for each sample
}
\usage{
snpgdsSampMissRate(gdsobj, sample.id=NULL, snp.id=NULL, with.id=FALSE)
}
\arguments{
    \item{gdsobj}{an object of class \code{\link{SNPGDSFileClass}},
        a SNP GDS file}
    \item{sample.id}{a vector of sample id specifying selected samples;
        if \code{NULL}, all samples will be used}
    \item{snp.id}{a vector of snp id specifying selected SNPs;
        if \code{NULL}, all SNPs will be used}
    \item{with.id}{if \code{TRUE}, the returned value with sample id}
}
\value{
    A vector of numeric values.
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{snpgdsSNPRateFreq}}
}

\examples{
# open an example dataset (HapMap)
genofile <- snpgdsOpen(snpgdsExampleFileName())

RV <- snpgdsSampMissRate(genofile)
summary(RV)

# close the genotype file
snpgdsClose(genofile)
}

\keyword{GDS}
\keyword{GWAS}
