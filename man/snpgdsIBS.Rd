\name{snpgdsIBS}
\alias{snpgdsIBS}
\title{
    Identity-By-State (IBS) proportion
}
\description{
    Calculate the fraction of identity by state for each pair of samples
}
\usage{
snpgdsIBS(gdsobj, sample.id = NULL, snp.id = NULL, autosome.only = TRUE,
    remove.monosnp = TRUE, maf = NaN, missing.rate = NaN,
    num.thread = 1, verbose = TRUE)
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
    \item{missing.rate}{to use the SNPs with "<= missing.rate" only;
        if NaN, no missing threshold}
    \item{num.thread}{the number of (CPU) cores used; if \code{NA}, detect
        the number of cores automatically}
    \item{verbose}{if TRUE, show information}
}
\details{
    The minor allele frequency and missing rate for each SNP passed in
\code{snp.id} are calculated over all the samples in \code{sample.id}.

    The values of the IBS matrix range from ZERO to ONE.
}
\value{
    Return a list (class "snpgdsIBSClass"):
    \item{sample.id}{the sample ids used in the analysis}
    \item{snp.id}{the SNP ids used in the analysis}
    \item{ibs}{a matrix of IBS proportion, "# of samples" x "# of samples"}
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{snpgdsIBSNum}}
}

\examples{
# open an example dataset (HapMap)
genofile <- snpgdsOpen(snpgdsExampleFileName())

# perform identity-by-state calculations
ibs <- snpgdsIBS(genofile)

# perform multidimensional scaling analysis on
# the genome-wide IBS pairwise distances:
loc <- cmdscale(1 - ibs$ibs, k = 2)
x <- loc[, 1]; y <- loc[, 2]
race <- as.factor(read.gdsn(index.gdsn(genofile, "sample.annot/pop.group")))
plot(x, y, col=race, xlab = "", ylab = "", main = "cmdscale(IBS Distance)")
legend("topleft", legend=levels(race), text.col=1:nlevels(race))

# close the file
snpgdsClose(genofile)
}

\keyword{GDS}
\keyword{GWAS}
