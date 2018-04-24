\name{snpgdsSNPList}
\alias{snpgdsSNPList}
\title{
    Create a SNP list object
}
\description{
    A list object of SNP information including rs, chr, pos, allele and allele
frequency.
}
\usage{
snpgdsSNPList(gdsobj, sample.id=NULL)
}
\arguments{
    \item{gdsobj}{an object of class \code{\link{SNPGDSFileClass}},
        a SNP GDS file}
    \item{sample.id}{a vector of sample id specifying selected samples;
        if NULL, all samples are used}
}
\value{
    Return an object of \code{snpgdsSNPListClass} including the following
components:
    \item{snp.id}{SNP id}
    \item{chromosome}{SNP chromosome index}
    \item{position}{SNP physical position in basepair}
    \item{allele}{reference / non-ref alleles}
    \item{afreq}{allele frequency}
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{snpgdsSNPListIntersect}}
}

\examples{
# open an example dataset (HapMap)
genofile <- snpgdsOpen(snpgdsExampleFileName())

# to get a snp list object
snplist <- snpgdsSNPList(genofile)
head(snplist)

# close the file
snpgdsClose(genofile)
}

\keyword{GDS}
\keyword{GWAS}
