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
    \item{gdsobj}{an object of class \code{\link{SNPGDSFileClass}},
        a SNP GDS file}
    \item{sample.id}{a vector of sample id specifying selected samples;
        if \code{NULL}, all samples will be used}
    \item{snp.id}{a vector of snp id specifying selected SNPs;
        if \code{NULL}, all SNPs will be used}
    \item{autosome.only}{if \code{TRUE}, use autosomal SNPs only}
    \item{remove.monosnp}{if \code{TRUE}, remove monomorphic SNPs}
    \item{maf}{to use the SNPs with ">= maf" only; if \code{NaN},
        no any MAF threshold}
    \item{missing.rate}{to use the SNPs with "<= missing.rate" only;
        if \code{NaN}, no any missing threshold}
    \item{verbose}{if \code{TRUE}, show information}
}
\value{
    Return a list of snp ids.
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{snpgdsSampMissRate}}, \code{\link{snpgdsSNPRateFreq}},
    \code{\link{snpgdsLDpruning}}
}

\examples{
# open an example dataset (HapMap)
genofile <- snpgdsOpen(snpgdsExampleFileName())

snpset <- snpgdsSelectSNP(genofile, maf=0.05, missing.rate=0.95)
length(snpset)
# 7502

# close the genotype file
snpgdsClose(genofile)
}

\keyword{GDS}
\keyword{GWAS}
