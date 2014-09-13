\name{snpgdsSlidingWindow}
\alias{snpgdsSlidingWindow}
\title{
    Sliding window
}
\description{
    Apply a user-defined function with a sliding window.
}
\usage{
snpgdsSlidingWindow(snp.id, position, chromosome=NULL, win.size, shift, FUN,
    param=c("id", "index+id", "index+id+pos"), verbose=TRUE, ...)
}
\arguments{
    \item{snp.id}{a vector, SNP IDs}
    \item{position}{SNP positions}
    \item{chromosome}{chromosome coding}
    \item{win.size}{the size of sliding window}
    \item{shift}{the shift of sliding window}
    \item{FUN}{the user-defined function}
    \item{param}{"id" -- calling \code{FUN(snp.id, ...)};
        "index+id" -- calling \code{FUN(index, snp.id, ...)};
        "index+id+pos" -- calling \code{FUN(index, snp.id, pos, ...)}}
    \item{verbose}{if TRUE, show information}
    \item{...}{optional arguments to \code{FUN}}
}
\value{
    Return a list
}
\author{Xiuwen Zheng}

\examples{
# open an example dataset (HapMap)
genofile <- snpgdsOpen(snpgdsExampleFileName())

set.seed(1000)
snpset <- snpgdsLDpruning(genofile)
names(snpset)
#  [1] "chr1"  "chr2"  "chr3"  "chr4"  "chr5"  "chr6"  "chr7"  "chr8"  "chr9"
# [10] "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18"
# ......
head(snpset$chr1)
# [1] 1 2 3 4 5 6

# get SNP ids
snp.id <- unlist(snpset)

# close the genotype file
snpgdsClose(genofile)
}

\keyword{GDS}
\keyword{GWAS}
