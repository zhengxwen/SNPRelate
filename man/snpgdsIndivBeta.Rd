\name{snpgdsIndivBeta}
\alias{snpgdsIndivBeta}
\title{
    Individual inbreeding and relatedness estimation (beta estimator)
}
\description{
    Calculate individual inbreeding and relatedness estimation (beta estimator)
using SNP genotype data.
}
\usage{
snpgdsIndivBeta(gdsobj, sample.id=NULL, snp.id=NULL, autosome.only=TRUE,
    remove.monosnp=TRUE, maf=NaN, missing.rate=NaN, method=c("weighted"),
    inbreeding=TRUE, num.thread=1L, with.id=TRUE, verbose=TRUE)
}
\arguments{
    \item{gdsobj}{an object of class \code{\link{SNPGDSFileClass}},
        a SNP GDS file}
    \item{sample.id}{a vector of sample id specifying selected samples;
        if NULL, all samples are used}
    \item{snp.id}{a vector of snp id specifying selected SNPs; if NULL,
        all SNPs are used}
    \item{autosome.only}{if \code{TRUE}, use autosomal SNPs only; if it is a
        numeric or character value, keep SNPs according to the specified
        chromosome}
    \item{remove.monosnp}{if TRUE, remove monomorphic SNPs}
    \item{maf}{to use the SNPs with ">= maf" only; if NaN, no MAF threshold}
    \item{missing.rate}{to use the SNPs with "<= missing.rate" only; if NaN,
        no missing threshold}
    \item{method}{"weighted" estimator}
    \item{inbreeding}{\code{TRUE}, the diagonal is a vector of inbreeding
        coefficients; otherwise, individual variance estimates}
    \item{num.thread}{the number of (CPU) cores used; if \code{NA}, detect
        the number of cores automatically}
    \item{with.id}{if \code{TRUE}, the returned value with \code{sample.id}
        and \code{sample.id}}
    \item{verbose}{if \code{TRUE}, show information}
}
\value{
    Return a list if \code{with.id = TRUE}:
        \item{sample.id}{the sample ids used in the analysis}
        \item{snp.id}{the SNP ids used in the analysis}
        \item{beta}{the genetic relationship matrix; different methods might
            have different meanings and interpretation for estimates}
    If \code{with.id = FALSE}, this function returns the genetic relationship
    matrix without sample and SNP IDs.
}
\references{
    Weir BS, Zheng X.
    SNPs and SNVs in Forensic Science.
    Forensic Science International: Genetics Supplement Series. 2015.
    doi:10.1016/j.fsigss.2015.09.106
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{snpgdsGRM}}, \code{\link{snpgdsIndInb}},
    \code{\link{snpgdsFst}}
}

\examples{
library(SNPRelate)
# open an example dataset (HapMap)
genofile <- snpgdsOpen(snpgdsExampleFileName())

beta <- snpgdsIndivBeta(genofile, with.id=FALSE)
beta[1:10, 1:10]


# close the file
snpgdsClose(genofile)
}

\keyword{GDS}
\keyword{GWAS}
