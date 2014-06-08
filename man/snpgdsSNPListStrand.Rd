\name{snpgdsSNPListStrand}
\alias{snpgdsSNPListStrand}
\title{
	Switch allele strand.
}
\description{
	To get a logical vector, indicating whether allele references of \code{snplist2}
need to be switched, with respect to \code{snplist1}.
}
\usage{
snpgdsSNPListStrand(snplist1, snplist2)
}
\arguments{
	\item{snplist1}{the first SNP list object ``snpgdsSNPListClass''}
	\item{snplist2}{the second SNP list object ``snpgdsSNPListClass''}
}
\value{
	a logical vector, where TRUE indiciates the allele references need to be switched, and NA
indiciates that locus is not in the common snp list.
}
\author{Xiuwen Zheng \email{zhengx@u.washington.edu}}
\seealso{
	\code{\link{snpgdsSNPList}}, \code{\link{snpgdsSNPListStrand}}
}

\examples{
# open an example dataset (HapMap)
genofile <- openfn.gds(snpgdsExampleFileName())

# to get a snp list object
snplist1 <- snpgdsSNPList(genofile)
snplist2 <- snpgdsSNPList(genofile)

# a common snp list
L <- snpgdsSNPListStrand(snplist1, snplist2)

table(L, exclude=NULL)

# close the genotype file
closefn.gds(genofile)
}

\keyword{gds}
\keyword{GWAS}
