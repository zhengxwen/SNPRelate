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
snpgdsFst(gdsobj, population, method=c("W&C84", "W&H02"), sample.id=NULL,
    snp.id=NULL, autosome.only=TRUE, remove.monosnp=TRUE, maf=NaN,
    missing.rate=NaN, with.id=FALSE, verbose=TRUE)
}
\arguments{
    \item{gdsobj}{an object of class \code{\link{SNPGDSFileClass}},
        a SNP GDS file}
    \item{population}{a factor, indicating population information for each
        individual}
    \item{method}{\code{"W&C84"} -- Fst estimator in Weir & Cockerham 1984
        (by default), \code{"W&H02"} -- relative beta estimator in
        Weir & Hill 2002, see details}
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

    The \code{"W&H02"} option implements the calculation in Buckleton et. al.
2016.
}
\value{
    Return a list:
    \item{sample.id}{the sample ids used in the analysis}
    \item{snp.id}{the SNP ids used in the analysis}
    \item{Fst}{weighted Fst estimate}
    \item{MeanFst}{the average of Fst estimates across SNPs}
    \item{FstSNP}{a vector of Fst for each SNP}
    \item{Beta}{Beta matrix}
}

\references{
    Weir, BS. & Cockerham, CC.
    Estimating F-statistics for the analysis of population structure. (1984).

    Weir, BS. & Hill, WG. Estimating F-statistics.
    Annual review of genetics 36, 721-50 (2002).

    Population-specific FST values for forensic STR markers: A worldwide survey.
    Buckleton J, Curran J, Goudet J, Taylor D, Thiery A, Weir BS.
    Forensic Sci Int Genet. 2016 Jul;23:91-100.
    doi: 10.1016/j.fsigen.2016.03.004.
}
\author{Xiuwen Zheng}

\examples{
# open an example dataset (HapMap)
genofile <- snpgdsOpen(snpgdsExampleFileName())

group <- as.factor(read.gdsn(index.gdsn(
    genofile, "sample.annot/pop.group")))

# Fst estimation
v <- snpgdsFst(genofile, population=group, method="W&C84")
v$Fst
v$MeanFst
summary(v$FstSNP)

# or
v <- snpgdsFst(genofile, population=group, method="W&H02")
v$Fst
v$MeanFst
v$Beta
summary(v$FstSNP)

# close the genotype file
snpgdsClose(genofile)
}

\keyword{GDS}
\keyword{GWAS}
