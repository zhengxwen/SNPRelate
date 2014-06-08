\name{snpgdsIBDSelection}
\alias{snpgdsIBDSelection}
\title{
	Get a table of IBD coefficients
}
\description{
	Return a data frame with IBD coefficients.
}
\usage{
snpgdsIBDSelection(ibdobj, kinship.cutoff=-1, samp.sel=NULL)
}
\arguments{
	\item{ibdobj}{an object of \code{snpgdsIBDClass} returned by \code{\link{snpgdsIBDMLE}} or \code{\link{snpgdsIBDMoM}}}
	\item{kinship.cutoff}{select the individual pairs with kinship coefficients >= kinship.cutoff}
	\item{samp.sel}{a logical vector or integer vector to specify selection of samples}
}
\value{
	Return a data.frame:
	\item{sample1}{the id of the first sample}
	\item{sample2}{the id of the second sample}
	\item{k0}{the probability of sharing ZERO alleles}
	\item{k1}{the probability of sharing ONE alleles}
	\item{kinshipcoeff}{kinship coefficient}
}
\author{Xiuwen Zheng \email{zhengx@u.washington.edu}}
\seealso{
	\code{\link{snpgdsIBDMLE}}, \code{\link{snpgdsIBDMoM}}
}

\examples{
# open an example dataset (HapMap)
genofile <- openfn.gds(snpgdsExampleFileName())

# YRI population
YRI.id <- read.gdsn(index.gdsn(genofile, "sample.id"))[
	read.gdsn(index.gdsn(genofile, c("sample.annot", "pop.group")))=="YRI"]
pibd <- snpgdsIBDMoM(genofile, sample.id=YRI.id)
flag <- lower.tri(pibd$k0)
plot(NaN, xlim=c(0,1), ylim=c(0,1), xlab="k0", ylab="k1")
lines(c(0,1), c(1,0), col="red", lty=3)
points(pibd$k0[flag], pibd$k1[flag])

# close the genotype file
closefn.gds(genofile)

# IBD coefficient list
dat <- snpgdsIBDSelection(pibd, 1/32)
#        sample1     sample2          k0          k1 kinshipcoeff
# 1      NA19152     NA19154 0.010796493 0.981324805   0.24927055
# 2      NA19139     NA19138 0.000000000 1.000000000   0.25000000
# 3      NA19139     NA19137 0.008646155 0.989673572   0.24825853
# 4      NA18912     NA18914 0.006489406 0.993510594   0.24837765
# 5      NA19160     NA19161 0.006033827 0.977803758   0.25253215
# ......
}

\keyword{gds}
\keyword{GWAS}
\keyword{IBD}
