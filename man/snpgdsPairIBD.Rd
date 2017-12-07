\name{snpgdsPairIBD}
\alias{snpgdsPairIBD}
\title{
    Calculate Identity-By-Descent (IBD) Coefficients
}
\description{
    Calculate the three IBD coefficients (k0, k1, k2) for non-inbred
individual pairs by Maximum Likelihood Estimation (MLE) or PLINK
Method of Moment (MoM).
}
\usage{
snpgdsPairIBD(geno1, geno2, allele.freq,
    method=c("EM", "downhill.simplex", "MoM", "Jacquard"),
    kinship.constraint=FALSE, max.niter=1000L, reltol=sqrt(.Machine$double.eps),
    coeff.correct=TRUE, out.num.iter=TRUE, verbose=TRUE)
}
\arguments{
    \item{geno1}{the SNP genotypes for the first individual,
        0 -- BB, 1 -- AB, 2 -- AA, other values -- missing}
    \item{geno2}{the SNP genotypes for the second individual,
        0 -- BB, 1 -- AB, 2 -- AA, other values -- missing}
    \item{allele.freq}{the allele frequencies}
    \item{method}{"EM", "downhill.simplex", "MoM" or "Jacquard", see details}
    \item{kinship.constraint}{if TRUE, constrict IBD coefficients
        ($k_0,k_1,k_2$) in the genealogical region ($2 k_0 k_1 >= k_2^2$)}
    \item{max.niter}{the maximum number of iterations}
    \item{reltol}{relative convergence tolerance; the algorithm stops if
        it is unable to reduce the value of log likelihood by a factor of
        $reltol * (abs(log likelihood with the initial parameters) + reltol)$
        at a step.}
    \item{coeff.correct}{\code{TRUE} by default, see details}
    \item{out.num.iter}{if TRUE, output the numbers of iterations}
    \item{verbose}{if TRUE, show information}
}
\details{
    If \code{method = "MoM"}, then PLINK Method of Moment without a
allele-count-based correction factor is conducted. Otherwise, two numeric
approaches for maximum likelihood estimation can be used: one is
Expectation-Maximization (EM) algorithm, and the other is Nelder-Mead method
or downhill simplex method. Generally, EM algorithm is more robust than
downhill simplex method. \code{"Jacquard"} refers to the estimation of nine
Jacquard's coefficients.

    If \code{coeff.correct} is \code{TRUE}, the final point that is found by
searching algorithm (EM or downhill simplex) is used to compare the six points
(fullsib, offspring, halfsib, cousin, unrelated), since any numeric approach
might not reach the maximum position after a finit number of steps. If any of
these six points has a higher value of log likelihood, the final point will be
replaced by the best one.
}
\value{
    Return a \code{data.frame}:
    \item{k0}{IBD coefficient, the probability of sharing ZERO IBD}
    \item{k1}{IBD coefficient, the probability of sharing ONE IBD}
    \item{loglik}{the value of log likelihood}
    \item{niter}{the number of iterations}
}
\references{
    Milligan BG. 2003. Maximum-likelihood estimation of relatedness.
    Genetics 163:1153-1167.

    Weir BS, Anderson AD, Hepler AB. 2006.
    Genetic relatedness analysis: modern data and new challenges.
    Nat Rev Genet. 7(10):771-80.

    Choi Y, Wijsman EM, Weir BS. 2009.
    Case-control association testing in the presence of unknown relationships.
    Genet Epidemiol 33(8):668-78.

    Purcell S, Neale B, Todd-Brown K, Thomas L, Ferreira MAR, Bender D,
    Maller J, Sklar P, de Bakker PIW, Daly MJ & Sham PC. 2007.
    PLINK: a toolset for whole-genome association and population-based linkage
    analysis. American Journal of Human Genetics, 81.
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{snpgdsPairIBDMLELogLik}}, \code{\link{snpgdsIBDMLE}},
    \code{\link{snpgdsIBDMLELogLik}}, \code{\link{snpgdsIBDMoM}}
}

\examples{
# open an example dataset (HapMap)
genofile <- snpgdsOpen(snpgdsExampleFileName())

YRI.id <- read.gdsn(index.gdsn(genofile, "sample.id"))[
    read.gdsn(index.gdsn(genofile, "sample.annot/pop.group"))=="YRI"]

# SNP pruning
set.seed(10)
snpset <- snpgdsLDpruning(genofile, sample.id=YRI.id, maf=0.05,
    missing.rate=0.05)
snpset <- unname(sample(unlist(snpset), 250))

# the number of samples
n <- 25

# specify allele frequencies
RF <- snpgdsSNPRateFreq(genofile, sample.id=YRI.id, snp.id=snpset,
    with.id=TRUE)
summary(RF$AlleleFreq)

subMLE <- snpgdsIBDMLE(genofile, sample.id=YRI.id[1:n], snp.id=RF$snp.id,
    allele.freq=RF$AlleleFreq)
subMoM <- snpgdsIBDMoM(genofile, sample.id=YRI.id[1:n], snp.id=RF$snp.id,
    allele.freq=RF$AlleleFreq)
subJac <- snpgdsIBDMLE(genofile, sample.id=YRI.id[1:n], snp.id=RF$snp.id,
    allele.freq=RF$AlleleFreq, method="Jacquard")



########################

# genotype matrix
mat <- snpgdsGetGeno(genofile, sample.id=YRI.id[1:n], snp.id=snpset,
    snpfirstdim=TRUE)

rv <- NULL
for (i in 2:n)
{
    rv <- rbind(rv, snpgdsPairIBD(mat[,1], mat[,i], RF$AlleleFreq, "EM"))
    print(snpgdsPairIBDMLELogLik(mat[,1], mat[,i], RF$AlleleFreq,
        relatedness="unrelated", verbose=TRUE))
}
rv
summary(rv$k0 - subMLE$k0[1, 2:n])
summary(rv$k1 - subMLE$k1[1, 2:n])
# ZERO

rv <- NULL
for (i in 2:n)
    rv <- rbind(rv, snpgdsPairIBD(mat[,1], mat[,i], RF$AlleleFreq, "MoM"))
rv
summary(rv$k0 - subMoM$k0[1, 2:n])
summary(rv$k1 - subMoM$k1[1, 2:n])
# ZERO

rv <- NULL
for (i in 2:n)
    rv <- rbind(rv, snpgdsPairIBD(mat[,1], mat[,i], RF$AlleleFreq, "Jacquard"))
rv
summary(rv$D1 - subJac$D1[1, 2:n])
summary(rv$D2 - subJac$D2[1, 2:n])
# ZERO

# close the genotype file
snpgdsClose(genofile)
}

\keyword{GDS}
\keyword{GWAS}
