\name{snpgdsHWE}
\alias{snpgdsHWE}
\title{
    Statistical test of Hardy-Weinberg Equilibrium
}
\description{
    Calculate the p-values for the exact SNP test of Hardy-Weinberg
Equilibrium.
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
\references{
    Wigginton, J. E., Cutler, D. J. & Abecasis, G. R.
    A note on exact tests of Hardy-Weinberg equilibrium.
    Am. J. Hum. Genet. 76, 887-93 (2005).
}
\author{Xiuwen Zheng, Janis E. Wigginton}
\seealso{
    \code{\link{snpgdsSNPRateFreq}}
}

\examples{
# open an example dataset (HapMap)
genofile <- snpgdsOpen(snpgdsExampleFileName())

# Japanese samples
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
pop <- read.gdsn(index.gdsn(genofile, "sample.annot/pop.group"))
(samp.sel <- sample.id[pop=="JPT"])
samp.sel <- samp.sel[nchar(samp.sel) == 7]

# chromosome 1
snp.id <- snpgdsSelectSNP(genofile, sample.id=samp.sel, autosome.only=1L)

# HWE test
p <- snpgdsHWE(genofile, sample.id=samp.sel, snp.id=snp.id)
summary(p)

# QQ plot
plot(-log10((1:length(p))/length(p)), -log10(p[order(p)]),
	xlab="-log10(expected P)", ylab="-log10(observed P)", main="QQ plot")
abline(a=0, b=1, col="blue")

# close the genotype file
snpgdsClose(genofile)
}

\keyword{GDS}
\keyword{GWAS}
