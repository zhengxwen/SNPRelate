\name{snpgdsIBSNum}
\alias{snpgdsIBSNum}
\title{
    Identity-By-State (IBS)
}
\description{
    Calculate the number of SNPs for identity by state for each pair of
samples.
}
\usage{
snpgdsIBSNum(gdsobj, sample.id = NULL, snp.id = NULL, autosome.only = TRUE,
    remove.monosnp = TRUE, maf = NaN, missing.rate = NaN, num.thread = 1,
    verbose = TRUE)
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
    \item{verbose}{if TRUE, show information}
}
\details{
    The minor allele frequency and missing rate for each SNP passed in
\code{snp.id} are calculated over all the samples in \code{sample.id}.
}
\value{
    Return a list (n is the number of samples):
    \item{sample.id}{the sample ids used in the analysis}
    \item{snp.id}{the SNP ids used in the analysis}
    \item{ibs0}{a n-by-n matrix, the number of SNPs sharing 0 IBS}
    \item{ibs1}{a n-by-n matrix, the number of SNPs sharing 1 IBS}
    \item{ibs2}{a n-by-n matrix, the number of SNPs sharing 2 IBS}
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{snpgdsIBS}}
}

\examples{
# open an example dataset (HapMap)
genofile <- snpgdsOpen(snpgdsExampleFileName())

RV <- snpgdsIBSNum(genofile)
pop <- read.gdsn(index.gdsn(genofile, "sample.annot/pop.group"))
L <- order(pop)
image(RV$ibs0[L, L]/length(RV$snp.id))

# close the genotype file
snpgdsClose(genofile)
}

\keyword{GDS}
\keyword{GWAS}
