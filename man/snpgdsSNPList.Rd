\name{snpgdsSNPList}
\alias{snpgdsSNPList}
\title{
	Create a SNP list object
}
\description{
	A list object of SNP information including rs, chr, pos, allele and allele frequency.
}
\usage{
snpgdsSNPList(gdsobj, sample.id=NULL)
}
\arguments{
	\item{gdsobj}{the \code{\link[gdsfmt]{gds.class}} object in the \link{gdsfmt} package}
	\item{sample.id}{a vector of sample id specifying selected samples; if NULL, all samples are used}
}
\value{
	Return an object of ``snpgdsSNPListClass'' including the following components:
	\item{rs.id}{SNP id}
	\item{chromosome}{SNP chromosome index}
	\item{position}{SNP physical position in basepair}
	\item{allele}{reference / non-ref alleles}
	\item{afreq}{allele frequency}
}
\author{Xiuwen Zheng \email{zhengx@u.washington.edu}}
\seealso{
	\code{\link{snpgdsSNPListIntersect}}, \code{\link{snpgdsSNPListStrand}}
}

\examples{
# open an example dataset (HapMap)
genofile <- openfn.gds(snpgdsExampleFileName())

# to get a snp list object
snplist <- snpgdsSNPList(genofile)

# close the genotype file
closefn.gds(genofile)
}

\keyword{gds}
\keyword{GWAS}
