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
    \item{gds}{a GDS file name, or an object of class
        \code{\link{SNPGDSFileClass}}}
    \item{show}{if TRUE, show information}
}
\value{
    Return a list:
    \item{sample.id}{the IDs of valid samples}
    \item{snp.id}{the IDs of valid SNPs}
}
\author{Xiuwen Zheng}

\examples{
snpgdsSummary(snpgdsExampleFileName())
}

\keyword{GDS}
\keyword{GWAS}
