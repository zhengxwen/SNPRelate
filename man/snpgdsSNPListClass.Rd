\name{snpgdsSNPListClass}
\alias{snpgdsSNPListClass}
\title{the class of a SNP list}
\description{
    the class of a SNP list, and its instance is returned from
\code{\link{snpgdsSNPList}}.
}

\value{
    Return an object of ``snpgdsSNPListClass'' including the following
components:
    \item{rs.id}{SNP id}
    \item{chromosome}{SNP chromosome index}
    \item{position}{SNP physical position in basepair}
    \item{allele}{reference / non-ref alleles}
    \item{afreq}{allele frequency}
}

\author{Xiuwen Zheng}
\seealso{\code{\link{snpgdsSNPList}}, \code{\link{snpgdsSNPListIntersect}}}

\keyword{GDS}
\keyword{GWAS}
