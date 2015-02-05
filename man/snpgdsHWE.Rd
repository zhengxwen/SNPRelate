\name{snpgdsHWE}
\alias{snpgdsHWE}
\title{
    Statistical test of Hardy-Weinberg Equilibrium
}
\description{
    Calculate the p-values for the exact SNP test of Hardy-Weinberg Equilibrium
}
\usage{
snpgdsHWE(gdsobj, sample.id=NULL, snp.id=NULL, with.id=FALSE)
}
\arguments{
    \item{gdsobj}{an object of class \code{\link{SNPGDSFileClass}},
        a SNP GDS file}
    \item{sample.id}{a vector of sample id specifying selected samples;
        if \code{NULL}, all samples will be used}
    \item{snp.id}{a vector of snp id specifying selected SNPs;
        if \code{NULL}, all SNPs will be used}
    \item{with.id}{if \code{TRUE}, the returned value with sample and SNP IDs}
}
\value{
    If \code{with.id=FALSE}, return a vector of numeric values (p-value);
otherwise, return a list with three components "pvalue", "sample.id" and
"snp.id".
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{snpgdsSNPRateFreq}}
}

\examples{
# open an example dataset (HapMap)
genofile <- snpgdsOpen(snpgdsExampleFileName())

RV <- snpgdsHWE(genofile)
summary(RV)

# close the genotype file
snpgdsClose(genofile)
}

\keyword{GDS}
\keyword{GWAS}
