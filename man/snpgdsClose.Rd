\name{snpgdsClose}
\alias{snpgdsClose}
\title{Close the SNP GDS File}
\description{
    Close the SNP GDS file
}
\usage{
snpgdsClose(gdsobj)
}
\arguments{
    \item{gdsobj}{an object of class \code{\link{SNPGDSFileClass}},
        a SNP GDS file}
}
\value{
    None.
}
\details{
    It is suggested to call \code{snpgdsClose} instead of
\code{\link[gdsfmt]{closefn.gds}}.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{snpgdsOpen}}
}

\examples{
# open an example dataset (HapMap)
genofile <- snpgdsOpen(snpgdsExampleFileName())

genofile

# close the file
snpgdsClose(genofile)
}

\keyword{GDS}
\keyword{GWAS}
