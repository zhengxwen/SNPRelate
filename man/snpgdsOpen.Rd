\name{snpgdsOpen}
\alias{snpgdsOpen}
\title{Open a SNP GDS File}
\description{
    Open a SNP GDS file
}
\usage{
snpgdsOpen(filename, readonly=TRUE, allow.duplicate=FALSE, allow.fork=TRUE)
}
\arguments{
    \item{filename}{the file name}
    \item{readonly}{whether read-only or not}
    \item{allow.duplicate}{if \code{TRUE}, it is allowed to open a GDS file
        with read-only mode when it has been opened in the same R session,
        see \code{\link[gdsfmt]{openfn.gds}}}
    \item{allow.fork}{\code{TRUE} for parallel environment using forking,
        see \code{\link[gdsfmt]{openfn.gds}}}
}
\value{
    Return an object of class \code{\link{SNPGDSFileClass}}.
}
\details{
    It is strongly suggested to call \code{snpgdsOpen} instead of
\code{\link[gdsfmt]{openfn.gds}}, since \code{snpgdsOpen} will perform
internal checking for data integrality.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{snpgdsClose}}
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
