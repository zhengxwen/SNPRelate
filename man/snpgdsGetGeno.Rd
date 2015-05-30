\name{snpgdsGetGeno}
\alias{snpgdsGetGeno}
\title{
    To get a genotype matrix
}
\description{
    To get a genotype matrix from a specified GDS file
}
\usage{
snpgdsGetGeno(gdsobj, sample.id=NULL, snp.id=NULL, snpfirstdim=NULL,
    with.id=FALSE, verbose=TRUE)
}
\arguments{
    \item{gdsobj}{an object of class \code{\link{SNPGDSFileClass}},
        a SNP GDS file; or characters to specify the file name of SNP GDS}
    \item{sample.id}{a vector of sample id specifying selected samples;
        if NULL, all samples are used}
    \item{snp.id}{a vector of snp id specifying selected SNPs;
        if NULL, all SNPs are used}
    \item{snpfirstdim}{if \code{TRUE}, genotypes are stored in the
        individual-major mode, (i.e, list all SNPs for the first individual,
        and then list all SNPs for the second individual, etc); if \code{NULL},
        determine automatically}
    \item{with.id}{if \code{TRUE}, return \code{sample.id} and \code{snp.id}}
    \item{verbose}{if TRUE, show information}
}
\value{
    The function returns an integer matrix with values 0, 1, 2 or NA
representing the number of reference allele when \code{with.id=FALSE};
or \code{list(genotype, sample.id, snp.id)} when \code{with.id=TRUE}. The
orders of sample and SNP IDs in the genotype matrix are actually consistent
with \code{sample.id} and \code{snp.id} in the GDS file, which may not be
as the same as the arguments \code{sampel.id} and \code{snp.id} specified
by users.
}
\author{Xiuwen Zheng}

\examples{
# open an example dataset (HapMap)
genofile <- snpgdsOpen(snpgdsExampleFileName())

set.seed(1000)
snpset <- sample(read.gdsn(index.gdsn(genofile, "snp.id")), 1000)

mat1 <- snpgdsGetGeno(genofile, snp.id=snpset, snpfirstdim=TRUE)
dim(mat1)
# 1000  279
table(c(mat1), exclude=NULL)

mat2 <- snpgdsGetGeno(genofile, snp.id=snpset, snpfirstdim=FALSE)
dim(mat2)
# 279 1000
table(c(mat2), exclude=NULL)

# close the file
snpgdsClose(genofile)
}

\keyword{GDS}
\keyword{GWAS}
