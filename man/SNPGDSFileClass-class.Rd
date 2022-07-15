\name{SNPGDSFileClass}
\docType{class}
\alias{SNPGDSFileClass-class}
\alias{SNPGDSFileClass}

\title{SNPGDSFileClass}

\description{
    A \code{SNPGDSFileClass} object provides access to a GDS file containing
genome-wide SNP data. It extends the class \code{\link{gds.class}} in
the gdsfmt package.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{snpgdsOpen}}, \code{\link{snpgdsClose}}
}

\examples{
# open an example dataset (HapMap)
genofile <- snpgdsOpen(snpgdsExampleFileName())
genofile

class(genofile)
# "SNPGDSFileClass" "gds.class"

# close the file
snpgdsClose(genofile)
}

\keyword{GDS}
\keyword{GWAS}
