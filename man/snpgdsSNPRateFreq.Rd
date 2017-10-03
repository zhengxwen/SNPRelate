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
snpgdsSNPRateFreq(gdsobj, sample.id=NULL, snp.id=NULL, with.id=FALSE,
    with.sample.id=FALSE, with.snp.id=FALSE)
}
\arguments{
    \item{gdsobj}{an object of class \code{\link{SNPGDSFileClass}},
        a SNP GDS file}
    \item{sample.id}{a vector of sample id specifying selected samples;
        if \code{NULL}, all samples will be used}
    \item{snp.id}{a vector of snp id specifying selected SNPs;
        if \code{NULL}, all SNPs will be used}
    \item{with.id}{if \code{TRUE}, return both sample and SNP IDs}
    \item{with.sample.id}{if \code{TRUE}, return sample IDs}
    \item{with.snp.id}{if \code{TRUE}, return SNP IDs}
}
\value{
    Return a list:
    \item{AlleleFreq}{allele frequencies}
    \item{MinorFreq}{minor allele frequencies}
    \item{MissingRate}{missing rates}
    \item{sample.id}{sample id, if \code{with.id=TRUE} or
        \code{with.sample.id=TRUE}}
    \item{snp.id}{SNP id, if \code{with.id=TRUE} or \code{with.snp.id=TRUE}}
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{snpgdsSampMissRate}}
}

\examples{
# open an example dataset (HapMap)
genofile <- snpgdsOpen(snpgdsExampleFileName())

RV <- snpgdsSNPRateFreq(genofile, with.snp.id=TRUE)
head(data.frame(RV))

hist(RV$AlleleFreq, breaks=128)
summary(RV$MissingRate)

# close the file
snpgdsClose(genofile)
}

\keyword{GDS}
\keyword{GWAS}
