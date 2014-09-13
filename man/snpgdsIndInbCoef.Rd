\name{snpgdsIndInbCoef}
\alias{snpgdsIndInbCoef}
\title{
    Individual Inbreeding Coefficient
}
\description{
    To calculate an individual inbreeding coefficient using SNP genotype data
}
\usage{
snpgdsIndInbCoef(x, p, method = c("mom.weir", "mom.visscher", "mle"),
    reltol=.Machine$double.eps^0.75)
}
\arguments{
    \item{x}{SNP genotypes}
    \item{p}{allele frequencies}
    \item{method}{see details}
    \item{reltol}{relative convergence tolerance used in MLE; the algorithm
        stops if it is unable to reduce the value of log likelihood by a
        factor of $reltol * (abs(log likelihood with the initial parameters)
        + reltol)$ at a step.}
}
\details{
    The \code{method} can be:
    \code{"mom.weir"}: a modified Visscher's estimator, proposed by Bruce Weir;
    \code{"mom.visscher"}: Visscher's estimator described in Yang et al. (2010);
    \code{"mle"}: the maximum likelihood estimation.
}
\value{
    Return estimated inbreeding coefficient.
}
\references{
    Yang J, Benyamin B, McEvoy BP, Gordon S, Henders AK, Nyholt DR, Madden PA,
    Heath AC, Martin NG, Montgomery GW, Goddard ME, Visscher PM. 2010.
    Common SNPs explain a large proportion of the heritability for human
    height. Nat Genet. 42(7):565-9. Epub 2010 Jun 20.
}
\author{Xiuwen Zheng}

\examples{
# open an example dataset (HapMap)
genofile <- snpgdsOpen(snpgdsExampleFileName())

chr1 <- read.gdsn(index.gdsn(genofile, "snp.id"))[
    read.gdsn(index.gdsn(genofile, "snp.chromosome"))==1]
chr1idx <- match(chr1, read.gdsn(index.gdsn(genofile, "snp.id")))

AF <- snpgdsSNPRateFreq(genofile)
g <- read.gdsn(index.gdsn(genofile, "genotype"), start=c(1,1), count=c(-1,1))

snpgdsIndInbCoef(g[chr1idx], AF$AlleleFreq[chr1idx], method="mom.weir")
snpgdsIndInbCoef(g[chr1idx], AF$AlleleFreq[chr1idx], method="mom.visscher")
snpgdsIndInbCoef(g[chr1idx], AF$AlleleFreq[chr1idx], method="mle")

# close the genotype file
snpgdsClose(genofile)
}

\keyword{GDS}
\keyword{GWAS}
