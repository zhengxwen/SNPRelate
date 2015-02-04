\name{snpgdsTranspose}
\alias{snpgdsTranspose}
\title{
    Transpose genotypic matrix
}
\description{
    Transpose the genotypic matrix if needed.
}
\usage{
snpgdsTranspose(gds.fn, snpfirstdim=FALSE, compress=NULL, optimize=TRUE,
    verbose=TRUE)
}
\arguments{
    \item{gds.fn}{the file name of SNP GDS format}
    \item{snpfirstdim}{if \code{TRUE}, genotypes are stored in snp-by-sample;
        if \code{FALSE}, sample-by-snp mode; if \code{NA}, force to transpose
        the SNP matrix}
    \item{compress}{the compression mode for SNP genotypes, optional values
        are defined in the function of \code{add.gdsn}; if \code{NULL}, to use
        the compression mode}
    \item{optimize}{if \code{TRUE}, call \code{\link[gdsfmt]{cleanup.gds}}
        after transposing}
    \item{verbose}{if TRUE, show information}
}
\value{
    None.
}

\author{Xiuwen Zheng}
\examples{
# the file name of SNP GDS
(fn <- snpgdsExampleFileName())

# copy the file
file.copy(fn, "test.gds", overwrite=TRUE)

# summary
snpgdsSummary("test.gds")

# transpose the SNP matrix
snpgdsTranspose("test.gds", snpfirstdim=TRUE)

# summary
snpgdsSummary("test.gds")


# delete the temporary file
unlink("test.gds", force=TRUE)
}

\keyword{GDS}
\keyword{GWAS}
