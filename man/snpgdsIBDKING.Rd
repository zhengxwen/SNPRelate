\name{snpgdsIBDKING}
\alias{snpgdsIBDKING}
\title{
    KING method of moment for the identity-by-descent (IBD) analysis
}
\description{
    Calculate IBD coefficients by KING method of moment.
}
\usage{
snpgdsIBDKING(gdsobj, sample.id=NULL, snp.id=NULL, autosome.only=TRUE,
    remove.monosnp=TRUE, maf=NaN, missing.rate=NaN,
    type=c("KING-robust", "KING-homo"), family.id=NULL, num.thread=1L,
    useMatrix=FALSE, verbose=TRUE)
}
\arguments{
    \item{gdsobj}{an object of class \code{\link{SNPGDSFileClass}},
        a SNP GDS file}
    \item{sample.id}{a vector of sample id specifying selected samples; if NULL,
        all samples are used}
    \item{snp.id}{a vector of snp id specifying selected SNPs;
        if NULL, all SNPs are used}
    \item{autosome.only}{if \code{TRUE}, use autosomal SNPs only; if it is a
        numeric or character value, keep SNPs according to the specified
        chromosome}
    \item{remove.monosnp}{if TRUE, remove monomorphic SNPs}
    \item{maf}{to use the SNPs with ">= maf" only; if NaN, no MAF threshold}
    \item{missing.rate}{to use the SNPs with "<= missing.rate" only; if NaN, no
        missing threshold}
    \item{type}{\code{"KING-robust"} -- relationship inference in the presence
        of population stratification; \code{"KING-homo"} -- relationship
        inference in a homogeneous population}
    \item{family.id}{if \code{NULL}, all individuals are treated as singletons;
        if family id is given, within- and between-family relationship are
        estimated differently. If \code{sample.id=NULL}, \code{family.id}
        should have the same length as "sample.id" in the GDS file, otherwise
        \code{family.id} should have the same length and order as the argument
        \code{sample.id}}
    \item{num.thread}{the number of (CPU) cores used; if \code{NA}, detect
        the number of cores automatically}
    \item{useMatrix}{if \code{TRUE}, use \code{Matrix::dspMatrix} to store
        the output square matrix to save memory}
    \item{verbose}{if TRUE, show information}
}
\details{
    KING IBD estimator is a moment estimator, and it is computationally
efficient relative to MLE method. The approaches include \code{"KING-robust"}
-- robust relationship inference within or across families in the presence of
population substructure, and \code{"KING-homo"} -- relationship inference in
a homogeneous population.

    With \code{"KING-robust"}, the function would return the proportion of SNPs
with zero IBS (\code{IBS0}) and kinship coefficient (\code{kinship}). With
\code{"KING-homo"} it would return the probability of sharing one IBD
(\code{k1}) and the probability of sharing zero IBD (\code{k0}).

    The minor allele frequency and missing rate for each SNP passed in
\code{snp.id} are calculated over all the samples in \code{sample.id}.
}
\value{
    Return a list:
    \item{sample.id}{the sample ids used in the analysis}
    \item{snp.id}{the SNP ids used in the analysis}
    \item{k0}{IBD coefficient, the probability of sharing zero IBD}
    \item{k1}{IBD coefficient, the probability of sharing one IBD}
    \item{IBS0}{proportion of SNPs with zero IBS}
    \item{kinship}{the estimated kinship coefficients, if the parameter
        \code{kinship=TRUE}}
}
\references{
    Manichaikul A, Mychaleckyj JC, Rich SS, Daly K, Sale M, Chen WM.
    Robust relationship inference in genome-wide association studies.
    Bioinformatics. 2010 Nov 15;26(22):2867-73.
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{snpgdsIBDMLE}}, \code{\link{snpgdsIBDMoM}}
}

\examples{
# open an example dataset (HapMap)
genofile <- snpgdsOpen(snpgdsExampleFileName())

# CEU population
samp.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
CEU.id <- samp.id[
    read.gdsn(index.gdsn(genofile, "sample.annot/pop.group"))=="CEU"]



####  KING-robust:
####  relationship inference in the presence of population stratification
####               robust relationship inference across family

ibd.robust <- snpgdsIBDKING(genofile, sample.id=CEU.id)
names(ibd.robust)
# [1] "sample.id" "snp.id"    "afreq"     "IBS0"      "kinship"

# select a set of pairs of individuals
dat <- snpgdsIBDSelection(ibd.robust, 1/32)
head(dat)

plot(dat$IBS0, dat$kinship, xlab="Proportion of Zero IBS",
    ylab="Estimated Kinship Coefficient (KING-robust)")


# using Matrix
ibd.robust <- snpgdsIBDKING(genofile, sample.id=CEU.id, useMatrix=TRUE)
is(ibd.robust$IBS0)  # dspMatrix
is(ibd.robust$kinship)  # dspMatrix



####  KING-robust:
####  relationship inference in the presence of population stratification
####               within- and between-family relationship inference

# incorporate with pedigree information
family.id <- read.gdsn(index.gdsn(genofile, "sample.annot/family.id"))
family.id <- family.id[match(CEU.id, samp.id)]

ibd.robust2 <- snpgdsIBDKING(genofile, sample.id=CEU.id, family.id=family.id)
names(ibd.robust2)

# select a set of pairs of individuals
dat <- snpgdsIBDSelection(ibd.robust2, 1/32)
head(dat)

plot(dat$IBS0, dat$kinship, xlab="Proportion of Zero IBS",
    ylab="Estimated Kinship Coefficient (KING-robust)")



####  KING-homo: relationship inference in a homogeneous population

ibd.homo <- snpgdsIBDKING(genofile, sample.id=CEU.id, type="KING-homo")
names(ibd.homo)
# "sample.id" "snp.id"    "afreq"     "k0"        "k1"

# select a subset of pairs of individuals
dat <- snpgdsIBDSelection(ibd.homo, 1/32)
head(dat)

plot(dat$k0, dat$kinship, xlab="Pr(IBD=0)",
    ylab="Estimated Kinship Coefficient (KING-homo)")


# using Matrix
ibd.homo <- snpgdsIBDKING(genofile, sample.id=CEU.id, type="KING-homo",
    useMatrix=TRUE)
is(ibd.homo$k0)  # dspMatrix
is(ibd.homo$k1)  # dspMatrix


# close the genotype file
snpgdsClose(genofile)
}

\keyword{GDS}
\keyword{GWAS}
