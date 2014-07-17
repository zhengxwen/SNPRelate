\name{snpgdsIBDMLELogLik}
\alias{snpgdsIBDMLELogLik}
\title{
    Log likelihood for MLE method in the Identity-By-Descent (IBD) Analysis
}
\description{
    Calculate the log likelihood values from maximum likelihood estimation.
}
\usage{
snpgdsIBDMLELogLik(gdsobj, ibdobj, k0 = NaN, k1 = NaN,
    relatedness=c("", "self", "fullsib", "offspring",
    "halfsib", "cousin", "unrelated"))
}
\arguments{
    \item{gdsobj}{an object of class \code{\link{SNPGDSFileClass}},
        a SNP GDS file}
    \item{ibdobj}{the \code{snpgdsIBDClass} object returned from
        \link{snpgdsIBDMLE}}
    \item{k0}{specified IBD coefficient}
    \item{k1}{specified IBD coefficient}
    \item{relatedness}{specify a relatedness, otherwise use the values of
        k0 and k1}
}
\details{
    If (\code{relatedness} == "") and (k0 == NaN or k1 == NaN), then return
the log likelihood values for each (k0, k1) stored in ibdobj. \\
    If (\code{relatedness} == "") and (k0 != NaN) and (k1 != NaN), then return
the log likelihood values for a specific IBD coefficient (k0, k1). \\
    If \code{relatedness} is: "self", then k0 = 0, k1 = 0; "fullsib", then
k0 = 0.25, k1 = 0.5; "offspring", then k0 = 0, k1 = 1; "halfsib", then
k0 = 0.5, k1 = 0.5; "cousin", then k0 = 0.75, k1 = 0.25; "unrelated", then
k0 = 1, k1 = 0.
}
\value{
    Return a n-by-n matrix of log likelihood values, where n is the
number of samples.
}
\references{
    Milligan BG. 2003. Maximum-likelihood estimation of relatedness.
    Genetics 163:1153-1167.

    Weir BS, Anderson AD, Hepler AB. 2006.
        Genetic relatedness analysis: modern data and new challenges.
        Nat Rev Genet. 7(10):771-80.

    Choi Y, Wijsman EM, Weir BS. 2009.
        Case-control association testing in the presence of unknown
        relationships.
        Genet Epidemiol 33(8):668-78.
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{snpgdsIBDMLE}}, \code{\link{snpgdsIBDMoM}}
}

\examples{
# open an example dataset (HapMap)
genofile <- snpgdsOpen(snpgdsExampleFileName())

YRI.id <- read.gdsn(index.gdsn(genofile, "sample.id"))[
    read.gdsn(index.gdsn(genofile, "sample.annot/pop.group"))=="YRI"]
YRI.id <- YRI.id[1:30]

# SNP pruning
set.seed(10)
snpset <- snpgdsLDpruning(genofile, sample.id=YRI.id, maf=0.05,
    missing.rate=0.05)
snpset <- sample(unlist(snpset), 250)
mibd <- snpgdsIBDMLE(genofile, sample.id=YRI.id, snp.id=snpset)
names(mibd)

# select a set of pairs of individuals
d <- snpgdsIBDSelection(mibd, kinship.cutoff=1/8)
head(d)


# log likelihood

loglik <- snpgdsIBDMLELogLik(genofile, mibd)
loglik0 <- snpgdsIBDMLELogLik(genofile, mibd, relatedness="unrelated")

# likelihood ratio test
p.value <- pchisq(loglik - loglik0, 1, lower.tail=FALSE)


flag <- lower.tri(mibd$k0)
plot(NaN, xlim=c(0,1), ylim=c(0,1), xlab="k0", ylab="k1")
lines(c(0,1), c(1,0), col="red", lty=3)
points(mibd$k0[flag], mibd$k1[flag])

# specify the allele frequencies
afreq <- snpgdsSNPRateFreq(genofile, sample.id=YRI.id,
    snp.id=snpset)$AlleleFreq
subibd <- snpgdsIBDMLE(genofile, sample.id=YRI.id[1:25], snp.id=snpset,
    allele.freq=afreq)
summary(c(subibd$k0 - mibd$k0[1:25, 1:25]))
# ZERO
summary(c(subibd$k1 - mibd$k1[1:25, 1:25]))
# ZERO


# close the genotype file
snpgdsClose(genofile)
}

\keyword{GDS}
\keyword{GWAS}
