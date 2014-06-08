\name{snpgdsOption}
\alias{snpgdsOption}
\title{
	Option settings: chromosome coding, etc
}
\description{
	Return an option list used by the SNPRelate package or a GDS file
}
\usage{
snpgdsOption(gdsobj=NULL, autosome.start=1, autosome.end=22, ...)
}
\arguments{
	\item{gdsobj}{a GDS file object (\code{\link[gdsfmt]{gds.class}})}
	\item{autosome.start}{the starting index of autosome}
	\item{autosome.end}{the ending index of autosome}
	\item{...}{optional arguments for new chromosome coding}
}
\author{Xiuwen Zheng}

\examples{

# define the new chromosomes 'Z' and 'W'
snpgdsOption(Z=27, W=28)


# open an example dataset (HapMap)
genofile <- openfn.gds(snpgdsExampleFileName())

snpgdsOption(genofile)

# close the genotype file
closefn.gds(genofile)
}

\keyword{gds}
\keyword{GWAS}
