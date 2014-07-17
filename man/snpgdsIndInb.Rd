\name{snpgdsIndInb}
\alias{snpgdsIndInb}
\title{
    Individual Inbreeding Coefficients
}
\description{
    To calculate individual inbreeding coefficients using SNP genotype data
}
\usage{
snpgdsIndInb(gdsobj, sample.id=NULL, snp.id=NULL,
    autosome.only=TRUE, remove.monosnp=TRUE, maf=NaN, missing.rate=NaN,
    method=c("mom.weir", "mom.visscher", "mle"), allele.freq=NULL,
    out.num.iter=TRUE, reltol=.Machine$double.eps^0.75, verbose=TRUE)
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
    \item{maf}{to use the SNPs with ">= maf" only; if NaN, no MAF threshold}
    \item{missing.rate}{to use the SNPs with "<= missing.rate" only; if NaN,
        no missing threshold}
    \item{method}{see details}
    \item{allele.freq}{to specify the allele frequencies; if NULL, the allele
        frequencies are estimated from the given samples}
    \item{out.num.iter}{output the numbers of iterations}
    \item{reltol}{relative convergence tolerance used in MLE; the algorithm
        stops if it is unable to reduce the value of log likelihood by a
        factor of $reltol * (abs(log likelihood with the initial parameters)
        + reltol)$ at a step.}
    \item{verbose}{if TRUE, show information}
}
\details{
    The \code{method} can be:
    "mom.weir": a modified Visscher's estimator, proposed by Bruce Weir;
    "mom.visscher": Visscher's estimator described in Yang et al. (2010);
    "mle": the maximum likelihood estimation.
}
\value{
    Return estimated inbreeding coefficient.
}
\references{
    Yang J, Benyamin B, McEvoy BP, Gordon S, Henders AK, Nyholt DR, Madden PA,
    Heath AC, Martin NG, Montgomery GW, Goddard ME, Visscher PM. 2010.
    Common SNPs explain a large proportion of the heritability for human height.
    Nat Genet. 42(7):565-9. Epub 2010 Jun 20.
}
\author{Xiuwen Zheng}

\examples{
# open an example dataset (HapMap)
genofile <- snpgdsOpen(snpgdsExampleFileName())

rv <- snpgdsIndInb(genofile, method="mom.visscher")
head(rv$inbreeding)
summary(rv$inbreeding)

# close the genotype file
snpgdsClose(genofile)
}

\keyword{GDS}
\keyword{GWAS}
