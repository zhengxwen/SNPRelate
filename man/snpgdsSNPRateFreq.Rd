\name{snpgdsSNPRateFreq}
\alias{snpgdsSNPRateFreq}
\title{
    Allele Frequency, Minor Allele Frequency, Missing Rate of SNPs
}
\description{
    Calculate the allele frequency, minor allele frequency and missing rate
per SNP.
}
\usage{
snpgdsSNPRateFreq(gdsobj, sample.id=NULL, snp.id=NULL, with.id=FALSE)
}
\arguments{
    \item{gdsobj}{an object of class \code{\link{SNPGDSFileClass}},
        a SNP GDS file}
    \item{sample.id}{a vector of sample id specifying selected samples;
        if \code{NULL}, all samples will be used}
    \item{snp.id}{a vector of snp id specifying selected SNPs;
        if \code{NULL}, all SNPs will be used}
    \item{with.id}{if \code{TRUE}, return sample and SNP IDs}
}
\value{
    Return a list:
    \item{AlleleFreq}{allele frequencies}
    \item{MinorFreq}{minor allele frequencies}
    \item{MissingRate}{missing rates}
    \item{sample.id}{sample id, if \code{with.id=TRUE}}
    \item{snp.id}{SNP id, if \code{with.id=TRUE}}
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{snpgdsSampMissRate}}
}

\examples{
# open an example dataset (HapMap)
genofile <- snpgdsOpen(snpgdsExampleFileName())

RV <- snpgdsSNPRateFreq(genofile)
hist(RV$AlleleFreq, breaks=128)
summary(RV$MissingRate)

# close the file
snpgdsClose(genofile)
}

\keyword{GDS}
\keyword{GWAS}
