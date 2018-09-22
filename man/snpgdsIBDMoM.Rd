\name{snpgdsIBDMoM}
\alias{snpgdsIBDMoM}
\title{
    PLINK method of moment (MoM) for the Identity-By-Descent (IBD) Analysis
}
\description{
    Calculate three IBD coefficients for non-inbred individual pairs by
    PLINK method of moment (MoM).
}
\usage{
snpgdsIBDMoM(gdsobj, sample.id=NULL, snp.id=NULL, autosome.only=TRUE,
    remove.monosnp=TRUE, maf=NaN, missing.rate=NaN, allele.freq=NULL,
    kinship=FALSE, kinship.constraint=FALSE, num.thread=1L, useMatrix=FALSE,
    verbose=TRUE)
}
\arguments{
    \item{gdsobj}{an object of class \code{\link{SNPGDSFileClass}},
        a SNP GDS file}
    \item{sample.id}{a vector of sample id specifying selected samples;
        if \code{NULL}, all samples are used}
    \item{snp.id}{a vector of snp id specifying selected SNPs; if \code{NULL},
        all SNPs are used}
    \item{autosome.only}{if \code{TRUE}, use autosomal SNPs only; if it is a
        numeric or character value, keep SNPs according to the specified
        chromosome}
    \item{remove.monosnp}{if \code{TRUE}, remove monomorphic SNPs}
    \item{maf}{to use the SNPs with ">= maf" only; if NaN, no MAF threshold}
    \item{missing.rate}{to use the SNPs with "<= missing.rate" only;
        if NaN, no missing threshold}
    \item{allele.freq}{to specify the allele frequencies; if NULL, determine
        the allele frequencies from \code{gdsobj} using the specified samples;
        if \code{snp.id} is specified, \code{allele.freq} should have
        the same order as \code{snp.id}}
    \item{kinship}{if \code{TRUE}, output the estimated kinship coefficients}
    \item{kinship.constraint}{if TRUE, constrict IBD coefficients
        ($k_0,k_1,k_2$) in the geneloical region ($2 k_0 k_1 >= k_2^2$)}
    \item{num.thread}{the number of (CPU) cores used; if \code{NA}, detect
        the number of cores automatically}
    \item{useMatrix}{if \code{TRUE}, use \code{Matrix::dspMatrix} to store
        the output square matrix to save memory}
    \item{verbose}{if TRUE, show information}
}
\details{
    PLINK IBD estimator is a moment estimator, and it is computationally
efficient relative to MLE method. In the PLINK method of moment, a correction
factor based on allele counts is used to adjust for sampling. However, if
allele frequencies are specified, no correction factor is conducted since the
specified allele frequencies are assumed to be known without sampling.

    The minor allele frequency and missing rate for each SNP passed in
\code{snp.id} are calculated over all the samples in \code{sample.id}.
}
\value{
    Return a list:
    \item{sample.id}{the sample ids used in the analysis}
    \item{snp.id}{the SNP ids used in the analysis}
    \item{k0}{IBD coefficient, the probability of sharing ZERO IBD}
    \item{k1}{IBD coefficient, the probability of sharing ONE IBD}
    \item{kinship}{the estimated kinship coefficients, if the parameter
        \code{kinship=TRUE}}
}
\references{
    Purcell S, Neale B, Todd-Brown K, Thomas L, Ferreira MAR, Bender D,
    Maller J, Sklar P, de Bakker PIW, Daly MJ & Sham PC. 2007.
    PLINK: a toolset for whole-genome association and population-based
    linkage analysis. American Journal of Human Genetics, 81.
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{snpgdsIBDMLE}}, \code{\link{snpgdsIBDMLELogLik}}
}

\examples{
# open an example dataset (HapMap)
genofile <- snpgdsOpen(snpgdsExampleFileName())

#########################################################
# CEU population

CEU.id <- read.gdsn(index.gdsn(genofile, "sample.id"))[
    read.gdsn(index.gdsn(genofile, "sample.annot/pop.group"))=="CEU"]
pibd <- snpgdsIBDMoM(genofile, sample.id=CEU.id)
names(pibd)

flag <- lower.tri(pibd$k0)
plot(NaN, xlim=c(0,1), ylim=c(0,1), xlab="k0", ylab="k1")
lines(c(0,1), c(1,0), col="red", lty=3)
points(pibd$k0[flag], pibd$k1[flag])

# select a set of pairs of individuals
d <- snpgdsIBDSelection(pibd, kinship.cutoff=1/8)
head(d)


#########################################################
# YRI population

YRI.id <- read.gdsn(index.gdsn(genofile, "sample.id"))[
    read.gdsn(index.gdsn(genofile, "sample.annot/pop.group"))=="YRI"]
pibd <- snpgdsIBDMoM(genofile, sample.id=YRI.id)
flag <- lower.tri(pibd$k0)
plot(NaN, xlim=c(0,1), ylim=c(0,1), xlab="k0", ylab="k1")
lines(c(0,1), c(1,0), col="red", lty=3)
points(pibd$k0[flag], pibd$k1[flag])


# specify the allele frequencies
afreq <- snpgdsSNPRateFreq(genofile, sample.id=YRI.id)$AlleleFreq
aibd <- snpgdsIBDMoM(genofile, sample.id=YRI.id, allele.freq=afreq)
flag <- lower.tri(aibd$k0)
plot(NaN, xlim=c(0,1), ylim=c(0,1), xlab="k0", ylab="k1")
lines(c(0,1), c(1,0), col="red", lty=3)
points(aibd$k0[flag], aibd$k1[flag])

# analysis on a subset
subibd <- snpgdsIBDMoM(genofile, sample.id=YRI.id[1:25], allele.freq=afreq)
summary(c(subibd$k0 - aibd$k0[1:25, 1:25]))
# ZERO
summary(c(subibd$k1 - aibd$k1[1:25, 1:25]))
# ZERO


# close the genotype file
snpgdsClose(genofile)
}

\keyword{GDS}
\keyword{GWAS}
