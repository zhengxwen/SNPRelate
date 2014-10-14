\name{snpgdsGRM}
\alias{snpgdsGRM}
\title{
    Genetic Relationship Matrix (GRM) for SNP genotype data
}
\description{
    Eigen-analysis on the Genetic Relationship Matrix (GRM) for SNP
genotype data.
}
\usage{
snpgdsGRM(gdsobj, sample.id=NULL, snp.id=NULL, autosome.only=TRUE,
    remove.monosnp=TRUE, maf=NaN, missing.rate=NaN, num.thread=1L,
    with.id=TRUE, verbose=TRUE)
}
\arguments{
    \item{gdsobj}{an object of class \code{\link{SNPGDSFileClass}},
        a SNP GDS file}
    \item{sample.id}{a vector of sample id specifying selected samples;
        if NULL, all samples are used}
    \item{snp.id}{a vector of snp id specifying selected SNPs; if NULL,
        all SNPs are used}
    \item{autosome.only}{if TRUE, use autosomal SNPs only}
    \item{remove.monosnp}{if TRUE, remove monomorphic SNPs}
    \item{maf}{to use the SNPs with ">= maf" only; if NaN, no MAF threshold}
    \item{missing.rate}{to use the SNPs with "<= missing.rate" only; if NaN,
        no missing threshold}
    \item{num.thread}{the number of (CPU) cores used; if \code{NA}, detect
        the number of cores automatically}
    \item{with.id}{if \code{TRUE}, the returned value with \code{sample.id} and
        \code{sample.id}}
    \item{verbose}{if \code{TRUE}, show information}
}
\value{
    Return a list if \code{with.id = TRUE}:
        \item{sample.id}{the sample ids used in the analysis}
        \item{snp.id}{the SNP ids used in the analysis}
        \item{GRM}{the genetic relationship matrix}
    If \code{with.id = FALSE}, return the genetic relationship matrix without
    sample and SNP IDs.
}
\references{
    Yang, J., Lee, S. H., Goddard, M. E. & Visscher, P. M.
    GCTA: a tool for genome-wide complex trait analysis.
    American journal of human genetics 88, 76-82 (2011).
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{snpgdsIndInb}}, \code{\link{snpgdsFst}}
}

\examples{
# open an example dataset (HapMap)
genofile <- snpgdsOpen(snpgdsExampleFileName())

# Genetic relationship matrix
grm <- snpgdsGRM(genofile)

# close the file
snpgdsClose(genofile)
}

\keyword{GDS}
\keyword{GWAS}
