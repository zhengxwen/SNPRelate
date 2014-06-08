\name{snpgdsOption}
\alias{snpgdsOption}
\title{
	Option Settings
}
\description{
	Return an option list used by the SNPRelate package or a GDS file
}
\usage{
snpgdsOption(gdsobj=NULL, autosome.start=1, autosome.end=22)
}
\arguments{
	\item{gdsobj}{the \code{\link[gdsfmt]{gdsclass}} object in the \link{gdsfmt} package}
	\item{autosome.start}{the starting index of autosome}
	\item{autosome.end}{the ending index of autosome}
}
\author{Xiuwen Zheng \email{zhengx@u.washington.edu}}

\examples{
# open an example dataset (HapMap)
genofile <- openfn.gds(snpgdsExampleFileName())

snpgdsOption(genofile)

# close the genotype file
closefn.gds(genofile)
}

\keyword{gds}
\keyword{GWAS}
