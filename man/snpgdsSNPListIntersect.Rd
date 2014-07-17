\name{snpgdsSNPListIntersect}
\alias{snpgdsSNPListIntersect}
\title{
    Get a common SNP list between two SNP list objects
}
\description{
    Get a common SNP list by comparing their rs id, chromosome indices and
positions.
}
\usage{
snpgdsSNPListIntersect(snplist1, snplist2)
}
\arguments{
    \item{snplist1}{the first SNP list object \code{\link{snpgdsSNPListClass}}}
    \item{snplist2}{the second SNP list object \code{\link{snpgdsSNPListClass}}}
}
\value{
    Return an object of \code{snpgdsSNPListClass} including the following
components:
    \item{rs.id}{SNP id}
    \item{chromosome}{SNP chromosome index}
    \item{position}{SNP physical position in basepair}
    \item{allele}{reference / non-ref alleles from the first SNP list object}
    \item{afreq}{allele frequency from the first SNP list object}
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{snpgdsSNPList}}, \code{\link{snpgdsSNPListStrand}}
}

\examples{
# open an example dataset (HapMap)
genofile <- snpgdsOpen(snpgdsExampleFileName())

# to get a snp list object
snplist1 <- snpgdsSNPList(genofile)
snplist2 <- snpgdsSNPList(genofile)

# a common snp list
snplist <- snpgdsSNPListIntersect(snplist1, snplist2)

summary(snplist$afreq)

# close the file
snpgdsClose(genofile)
}

\keyword{GDS}
\keyword{GWAS}
