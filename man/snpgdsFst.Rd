\name{snpgdsFst}
\alias{snpgdsFst}
\title{
    F-statistics (fixation indices)
}
\description{
    Calculate relatedness measures F-statistics (also known as fixation
indices) for given populations
}
\usage{
snpgdsFst(gdsobj, population, method=c("W&H02", "W&C84"), sample.id=NULL,
    snp.id=NULL, autosome.only=TRUE, remove.monosnp=TRUE, maf=NaN,
    missing.rate=NaN, with.id=FALSE, verbose=TRUE)
}
\arguments{
    \item{gdsobj}{an object of class \code{\link{SNPGDSFileClass}},
        a SNP GDS file}
    \item{population}{a factor, indicating population information for each
        individual}
    \item{method}{\code{"W&H02"} -- relative beta estimator in Weir&Hill 2002,
        \code{"W&C84"} -- Fst estimator in Weir&Cockerham 1984}
    \item{sample.id}{a vector of sample id specifying selected samples;
        if NULL, all samples are used}
    \item{snp.id}{a vector of snp id specifying selected SNPs;
        if NULL, all SNPs are used}
    \item{autosome.only}{if \code{TRUE}, use autosomal SNPs only; if it is a
        numeric or character value, keep SNPs according to the specified
        chromosome}
    \item{remove.monosnp}{if TRUE, remove monomorphic SNPs}
    \item{maf}{to use the SNPs with ">= maf" only; if NaN, no MAF threshold}
    \item{missing.rate}{to use the SNPs with "<= missing.rate" only;
        if NaN, no missing threshold}
    \item{with.id}{if \code{TRUE}, the returned value with \code{sample.id} and
        \code{sample.id}}
    \item{verbose}{if \code{TRUE}, show information}
}
\details{
    The minor allele frequency and missing rate for each SNP passed in
\code{snp.id} are calculated over all the samples in \code{sample.id}.
}
\value{
    Return a list:
    \item{sample.id}{the sample ids used in the analysis}
    \item{snp.id}{the SNP ids used in the analysis}
    \item{Fst}{Fst estimator}
    \item{Beta}{Beta matrix}
}

\references{
    Weir, B. S. & Cockerham, C. C.
    Estimating F-statistics for the analysis of population structure. (1984).

    Weir, B. S. & Hill, W. G. Estimating F-statistics.
    Annual review of genetics 36, 721-50 (2002).
}
\author{Xiuwen Zheng}

\examples{
# open an example dataset (HapMap)
genofile <- snpgdsOpen(snpgdsExampleFileName())

group <- as.factor(read.gdsn(index.gdsn(
    genofile, "sample.annot/pop.group")))

# Fst estimation
snpgdsFst(genofile, population=group, method="W&H02")

# or
snpgdsFst(genofile, population=group, method="W&C84")


# close the genotype file
snpgdsClose(genofile)
}

\keyword{GDS}
\keyword{GWAS}
