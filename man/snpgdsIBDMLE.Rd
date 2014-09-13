\name{snpgdsIBDMLE}
\alias{snpgdsIBDMLE}
\title{
    Maximum likelihood estimation (MLE) for the Identity-By-Descent (IBD)
Analysis
}
\description{
    Calculate the three IBD coefficients (k0, k1, k2) for non-inbred
individual pairs by Maximum Likelihood Estimation.
}
\usage{
snpgdsIBDMLE(gdsobj, sample.id=NULL, snp.id=NULL, autosome.only=TRUE,
    remove.monosnp=TRUE, maf=NaN, missing.rate=NaN, kinship=FALSE,
    kinship.constraint=FALSE, allele.freq=NULL,
    method=c("EM", "downhill.simplex", "Jacquard"), max.niter=1000L,
    reltol=sqrt(.Machine$double.eps), coeff.correct=TRUE,
    out.num.iter=TRUE, num.thread=1, verbose=TRUE)
}
\arguments{
    \item{gdsobj}{an object of class \code{\link{SNPGDSFileClass}},
        a SNP GDS file}
    \item{sample.id}{a vector of sample id specifying selected samples;
        if NULL, all samples are used}
    \item{snp.id}{a vector of snp id specifying selected SNPs;
        if NULL, all SNPs are used}
    \item{autosome.only}{if TRUE, use autosomal SNPs only}
    \item{remove.monosnp}{if TRUE, remove monomorphic SNPs}
    \item{maf}{to use the SNPs with ">= maf" only; if NaN, no any MAF threshold}
    \item{missing.rate}{to use the SNPs with "<= missing.rate" only;
        if NaN, no any missing threshold}
    \item{kinship}{if \code{TRUE}, output the estimated kinship coefficients}
    \item{kinship.constraint}{if TRUE, constrict IBD coefficients
        ($k_0,k_1,k_2$) in the geneloical region ($2 k_0 k_1 >= k_2^2$)}
    \item{allele.freq}{to specify the allele frequencies; if NULL, determine
        the allele frequencies from \code{gdsobj} using the specified samples;
        if \code{snp.id} is specified, \code{allele.freq} should have
        the same order as \code{snp.id}}
    \item{method}{"EM", "downhill.simplex", "Jacquard", see details}
    \item{max.niter}{the maximum number of iterations}
    \item{reltol}{relative convergence tolerance; the algorithm stops if it
        is unable to reduce the value of log likelihood by a factor of
        $reltol * (abs(log likelihood with the initial parameters) + reltol)$
        at a step.}
    \item{coeff.correct}{\code{TRUE} by default, see details}
    \item{out.num.iter}{if TRUE, output the numbers of iterations}
    \item{num.thread}{the number of (CPU) cores used; if \code{NA}, detect
        the number of cores automatically}
    \item{verbose}{if TRUE, show information}
}
\details{
    The minor allele frequency and missing rate for each SNP passed in
\code{snp.id} are calculated over all the samples in \code{sample.id}.

    The PLINK moment estimates are used as the initial values in the algorithm
of searching maximum value of log likelihood function. Two numeric approaches
can be used: one is Expectation-Maximization (EM) algorithm, and the other is
Nelder-Mead method or downhill simplex method. Generally, EM algorithm is more
robust than downhill simplex method. \code{"Jacquard"} refers to the estimation
of nine Jacquard's coefficients.

    If \code{coeff.correct} is \code{TRUE}, the final point that is found by
searching algorithm (EM or downhill simplex) is used to compare the six points
(fullsib, offspring, halfsib, cousin, unrelated), since any numeric approach
might not reach the maximum position after a finit number of steps. If any of
these six points has a higher value of log likelihood, the final point will
be replaced by the best one.

    Although MLE estimates are more reliable than MoM, MLE is much more
computationally intensive than MoM, and might not be feasible to estimate
pairwise relatedness for a large dataset.
}
\value{
    Return a \code{snpgdsIBDClass} object, which is a list:
    \item{sample.id}{the sample ids used in the analysis}
    \item{snp.id}{the SNP ids used in the analysis}
    \item{afreq}{the allele frequencies used in the analysis}
    \item{k0}{IBD coefficient, the probability of sharing ZERO IBD, if
        \code{method="EM"} or \code{"downhill.simplex"}}
    \item{k1}{IBD coefficient, the probability of sharing ONE IBD, if
        \code{method="EM"} or \code{"downhill.simplex"}}
    \item{D1, ..., D8}{Jacquard's coefficients, if \code{method="Jacquard"}}
    \item{kinship}{the estimated kinship coefficients, if the parameter
        \code{kinship=TRUE}}
}
\references{
    Milligan BG. 2003.
    Maximum-likelihood estimation of relatedness. Genetics 163:1153-1167.

    Weir BS, Anderson AD, Hepler AB. 2006.
    Genetic relatedness analysis: modern data and new challenges.
    Nat Rev Genet. 7(10):771-80.

    Choi Y, Wijsman EM, Weir BS. 2009.
    Case-control association testing in the presence of unknown relationships.
    Genet Epidemiol 33(8):668-78.

    Jacquard, A. Structures Genetiques des Populations
    (Masson & Cie, Paris, 1970); English translation available in
    Charlesworth, D. & Chalesworth, B. Genetics of Human Populations
    (Springer, New York, 1974).
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{snpgdsIBDMLELogLik}}, \code{\link{snpgdsIBDMoM}}
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
