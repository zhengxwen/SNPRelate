\name{snpgdsOption}
\alias{snpgdsOption}
\title{
    Option settings: chromosome coding, etc
}
\description{
    Return an option list used by the SNPRelate package or a GDS file
}
\usage{
snpgdsOption(gdsobj=NULL, autosome.start=1L, autosome.end=22L, ...)
}
\arguments{
    \item{gdsobj}{an object of class \code{\link{SNPGDSFileClass}},
        a SNP GDS file}
    \item{autosome.start}{the starting index of autosome}
    \item{autosome.end}{the ending index of autosome}
    \item{...}{optional arguments for new chromosome coding}
}
\author{Xiuwen Zheng}

\examples{

# define the new chromosomes 'Z' and 'W'
snpgdsOption(Z=27L, W=28L)


# open an example dataset (HapMap)
genofile <- snpgdsOpen(snpgdsExampleFileName())

snpgdsOption(genofile)

# close the genotype file
snpgdsClose(genofile)
}

\keyword{GDS}
\keyword{GWAS}
