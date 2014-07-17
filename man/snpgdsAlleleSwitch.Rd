\name{snpgdsAlleleSwitch}
\alias{snpgdsAlleleSwitch}
\title{
    Allele-switching
}
\description{
    Switch alleles according to the reference if needed.
}
\usage{
snpgdsAlleleSwitch(gdsobj, A.allele, verbose=TRUE)
}
\arguments{
    \item{gdsobj}{an object of class \code{\link{SNPGDSFileClass}},
        a SNP GDS file}
    \item{A.allele}{characters, referring to A allele}
    \item{verbose}{if TRUE, show information}
}
\value{
    A logical vector with \code{TRUE} indicating allele-switching and
\code{NA} when it is unable to determine. \code{NA} occurs when
\code{A.allele = NA} or \code{A.allele} is not in the list of alleles.
}

\author{Xiuwen Zheng}
\examples{
# the file name of SNP GDS
(fn <- snpgdsExampleFileName())

# copy the file
file.copy(fn, "test.gds", overwrite=TRUE)

# open the SNP GDS file
genofile <- snpgdsOpen("test.gds", readonly=FALSE)

# allelic information
allele <- read.gdsn(index.gdsn(genofile, "snp.allele"))
allele.list <- strsplit(allele, "/")

A.allele <- sapply(allele.list, function(x) { x[1] })
B.allele <- sapply(allele.list, function(x) { x[2] })

set.seed(1000)
flag <- rep(FALSE, length(A.allele))
flag[sample.int(length(A.allele), 50, replace=TRUE)] <- TRUE

A.allele[flag] <- B.allele[flag]
A.allele[sample.int(length(A.allele), 10, replace=TRUE)] <- NA
table(A.allele, exclude=NULL)


# allele switching
z <- snpgdsAlleleSwitch(genofile, A.allele)

table(z, exclude=NULL)


# close the file
snpgdsClose(genofile)


# delete the temporary file
unlink("test.gds", force=TRUE)
}

\keyword{GDS}
\keyword{GWAS}
