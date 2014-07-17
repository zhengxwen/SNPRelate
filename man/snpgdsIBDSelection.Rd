\name{snpgdsIBDSelection}
\alias{snpgdsIBDSelection}
\title{
    Get a table of IBD coefficients
}
\description{
    Return a data frame with IBD coefficients.
}
\usage{
snpgdsIBDSelection(ibdobj, kinship.cutoff=NaN, samp.sel=NULL)
}
\arguments{
    \item{ibdobj}{an object of \code{snpgdsIBDClass} returned by
        \code{\link{snpgdsIBDMLE}} or \code{\link{snpgdsIBDMoM}}}
    \item{kinship.cutoff}{select the individual pairs with kinship coefficients
        >= kinship.cutoff; no filter if \code{kinship.cutoff = NaN}}
    \item{samp.sel}{a logical vector or integer vector to specify selection
        of samples}
}
\value{
    Return a data.frame:
    \item{ID1}{the id of the first individual}
    \item{ID2}{the id of the second individual}
    \item{k0}{the probability of sharing ZERO alleles}
    \item{k1}{the probability of sharing ONE alleles}
    \item{kinship}{kinship coefficient}
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{snpgdsIBDMLE}}, \code{\link{snpgdsIBDMoM}},
    \code{\link{snpgdsIBDKING}}
}

\examples{
# open an example dataset (HapMap)
genofile <- snpgdsOpen(snpgdsExampleFileName())

# YRI population
YRI.id <- read.gdsn(index.gdsn(genofile, "sample.id"))[
    read.gdsn(index.gdsn(genofile, "sample.annot/pop.group"))=="YRI"]
pibd <- snpgdsIBDMoM(genofile, sample.id=YRI.id)
flag <- lower.tri(pibd$k0)
plot(NaN, xlim=c(0,1), ylim=c(0,1), xlab="k0", ylab="k1")
lines(c(0,1), c(1,0), col="red", lty=3)
points(pibd$k0[flag], pibd$k1[flag])

# close the genotype file
snpgdsClose(genofile)

# IBD coefficients
dat <- snpgdsIBDSelection(pibd, 1/32)
head(dat)
#       ID1     ID2          k0        k1    kinship
# 1 NA19152 NA19154 0.010749154 0.9892508 0.24731271
# 2 NA19152 NA19093 0.848207777 0.1517922 0.03794806
# 3 NA19139 NA19138 0.010788047 0.9770181 0.25035144
# 4 NA19139 NA19137 0.012900661 0.9870993 0.24677483
# 5 NA18912 NA18914 0.008633077 0.9913669 0.24784173
# 6 NA19160 NA19161 0.008635754 0.9847777 0.24948770
}

\keyword{GDS}
\keyword{GWAS}
