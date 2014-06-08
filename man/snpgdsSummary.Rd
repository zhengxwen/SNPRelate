\name{snpgdsSummary}
\alias{snpgdsSummary}
\title{
	Summary of GDS genotype file
}
\description{
	Print the information stored in the gds object
}
\usage{
snpgdsSummary(gds, show=TRUE)
}
\arguments{
	\item{gds}{the \code{\link[gdsfmt]{gds.class}} object in the \link{gdsfmt} package, or a file name}
	\item{show}{if TRUE, show information}
}
\value{
	Return a list:
	\item{sample.id}{the IDs of valid samples}
	\item{snp.id}{the IDs of valid SNPs}
}
\author{Xiuwen Zheng \email{zhengx@u.washington.edu}}

\examples{
snpgdsSummary(snpgdsExampleFileName())
}

\keyword{gds}
\keyword{GWAS}
