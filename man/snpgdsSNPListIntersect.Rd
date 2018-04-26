\name{snpgdsSNPListIntersect}
\alias{snpgdsSNPListIntersect}
\title{
    Get a common SNP list between/among SNP list objects
}
\description{
    Get a common SNP list by comparing their snp id, chromosome, positions
and allele frequency if needed.
}
\usage{
snpgdsSNPListIntersect(snplist1, snplist2, ..., method=c("position", "exact"),
    na.rm=TRUE, same.strand=FALSE, verbose=TRUE)
}
\arguments{
    \item{snplist1}{the SNP list object \code{\link{snpgdsSNPListClass}}}
    \item{snplist2}{the SNP list object \code{\link{snpgdsSNPListClass}}}
    \item{...}{the other SNP list objects}
    \item{method}{\code{"exact"}: matching by all snp.id, chromosomes, positions
        and alleles; \code{"position"}: matching by chromosomes and positions}
    \item{na.rm}{if TRUE, remove mismatched alleles}
    \item{same.strand}{if TRUE, assuming the alleles on the same strand}
    \item{verbose}{if TRUE, show information}
}
\value{
    Return a list of \code{snpgdsSNPListClass} including the following
components:
    \item{idx1}{the indices of common SNPs in the first GDS file}
    \item{idx2}{the indices of common SNPs in the second GDS file}
    \item{idx...}{}
    \item{idx\code{n}}{the indices of common SNPs in the n-th GDS file}
    \item{flag2}{an integer vector, flip flag for each common SNP for the
        second GDS file (assuming a value \code{v}):
        \code{bitwAnd(v, 1)}: 0 -- no flip of allele names,
            1 -- flip of allele names;
        \code{bitwAnd(v, 2)}: 0 -- on the same strand,
            2 -- on the different strands, comparing with the first GDS file;
        \code{bitwAnd(v, 4)}: 0 -- no strand ambiguity,
            4 -- ambiguous allele names, determined by allele frequencies;
        \code{NA} -- mismatched allele names
        (there is no NA if \code{na.rm=TRUE})}
    \item{flag...}{}
    \item{flag\code{n}}{flip flag for each common SNP for the n-th GDS file}
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{snpgdsSNPList}}
}

\examples{
# open an example dataset (HapMap)
genofile <- snpgdsOpen(snpgdsExampleFileName())

# to get a snp list object
snplist1 <- snpgdsSNPList(genofile)
snplist2 <- snpgdsSNPList(genofile)

# a common snp list, exactly matching
v <- snpgdsSNPListIntersect(snplist1, snplist2)
names(v)
# "idx1" "idx2"

# a common snp list, matching by position
v <- snpgdsSNPListIntersect(snplist1, snplist2, method="pos")
names(v)
# "idx1"  "idx2"  "flag2"

table(v$flag2, exclude=NULL)


# close the file
snpgdsClose(genofile)
}

\keyword{GDS}
\keyword{GWAS}
