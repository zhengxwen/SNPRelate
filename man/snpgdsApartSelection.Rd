\name{snpgdsApartSelection}
\alias{snpgdsApartSelection}
\title{Select SNPs with a basepair distance}
\description{
    Randomly selects SNPs for which each pair is at least as far apart as the
specified basepair distance.
}
\usage{
snpgdsApartSelection(chromosome, position, min.dist=100000,
    max.n.snp.perchr=-1, verbose=TRUE) 
}
\arguments{
    \item{chromosome}{chromosome codes}
    \item{position}{SNP positions in base pair}
    \item{min.dist}{A numeric value to specify minimum distance required
        (in basepairs)}
    \item{max.n.snp.perchr}{A numeric value specifying the maximum number
        of SNPs to return per chromosome, "-1" means no number limit}
    \item{verbose}{if TRUE, show information}
}
\value{
    A logical vector indicating which SNPs were selected.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{snpgdsLDpruning}}
}

\examples{
# open an example dataset (HapMap)
genofile <- snpgdsOpen(snpgdsExampleFileName())
genofile

chr <- read.gdsn(index.gdsn(genofile, "snp.chromosome"))
pos <- read.gdsn(index.gdsn(genofile, "snp.position"))

set.seed(1000)
flag <- snpgdsApartSelection(chr, pos, min.dist=250000, verbose=TRUE)
table(flag)

# close the genotype file
snpgdsClose(genofile)
}

\keyword{GDS}
\keyword{GWAS}
