\name{snpgdsSlidingWindow}
\alias{snpgdsSlidingWindow}
\title{
    Sliding window
}
\description{
    Apply a user-defined function with a sliding window.
}
\usage{
snpgdsSlidingWindow(gdsobj, sample.id=NULL, snp.id=NULL,
    FUN=NULL, winsize=100000L, shift=10000L, unit=c("basepair", "locus"),
    autosome.only=FALSE, remove.monosnp=TRUE, maf=NaN, missing.rate=NaN,
    with.id=c("snp.id", "snp.id.in.window", "none"),
    num.thread=1, verbose=TRUE, ...)
}
\arguments{
    \item{gdsobj}{an object of class \code{\link{SNPGDSFileClass}},
        a SNP GDS file}
    \item{sample.id}{a vector of sample id specifying selected samples;
        if NULL, all samples are used}
    \item{snp.id}{a vector of snp id specifying selected SNPs; if NULL,
        all SNPs are used}
    \item{FUN}{the user-defined function}
    \item{winsize}{the size of sliding window}
    \item{shift}{the shift of sliding window}
    \item{unit}{\code{"basepair"} -- \code{winsize} and \code{shift} are
        applied with SNP coordinate of basepair; \code{"locus"} --
        \code{winsize} and \code{shift} are applied with SNP index}
    \item{autosome.only}{if \code{TRUE}, use autosomal SNPs only; if it is a
        numeric or character value, keep SNPs according to the specified
        chromosome}
    \item{remove.monosnp}{if TRUE, remove monomorphic SNPs}
    \item{maf}{to use the SNPs with ">= maf" only; if NaN, no MAF threshold}
    \item{missing.rate}{to use the SNPs with "<= missing.rate" only; if NaN,
        no missing threshold}
    \item{with.id}{"snp.id", "snp.id.in.window" or "none"}
    \item{num.thread}{the number of (CPU) cores used; if \code{NA}, detect
        the number of cores automatically}
    \item{verbose}{if TRUE, show information}
    \item{...}{optional arguments to \code{FUN}}
}
\value{
    Return a list
}
\author{Xiuwen Zheng}

\examples{
# open an example dataset (HapMap)
genofile <- snpgdsOpen(snpgdsExampleFileName())

# sliding windows
rv <- snpgdsSlidingWindow(genofile, winsize=500000, shift=100000,
    FUN=function(...) NULL)

# plot
plot(rv$chr1.num, ylab="# of SNPs in the sliding window")

# close the genotype file
snpgdsClose(genofile)
}

\keyword{GDS}
\keyword{GWAS}
