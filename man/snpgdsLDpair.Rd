\name{snpgdsLDpair}
\alias{snpgdsLDpair}
\title{
    Linkage Disequilibrium (LD)
}
\description{
    Return a LD value between snp1 and snp2.
}
\usage{
snpgdsLDpair(snp1, snp2, method = c("composite", "r", "dprime", "corr"))
}
\arguments{
    \item{snp1}{a vector of SNP genotypes (0 -- BB, 1 -- AB, 2 -- AA)}
    \item{snp2}{a vector of SNP genotypes (0 -- BB, 1 -- AB, 2 -- AA)}
    \item{method}{"composite", "r", "dprime", "corr", see details}
}
\details{
    Four methods can be used to calculate linkage disequilibrium values:
"composite" for LD composite measure, "r" for R coefficient (by EM algorithm
assuming HWE, it could be negative), "dprime" for D', and "corr" for
correlation coefficient. The method "corr" is equivalent to "composite", when
SNP genotypes are coded as: 0 -- BB, 1 -- AB, 2 -- AA.
}
\value{
    Return a numeric vector:
    \item{ld}{a measure of linkage disequilibrium}

    if \code{method = "r" or "dprime"},
    \item{pA_A}{haplotype frequency of AA, the first locus is A and the
        second locus is A}
    \item{pA_B}{haplotype frequency of AB, the first locus is A and the
        second locus is B}
    \item{pB_A}{haplotype frequency of BA, the first locus is B and the
        second locus is A}
    \item{pB_B}{haplotype frequency of BB, the first locus is B and the
        second locus is B}
}
\references{
    Weir B: Inferences about linkage disequilibrium.
    Biometrics 1979; 35: 235-254.

    Weir B: Genetic Data Analysis II. Sunderland, MA: Sinauer Associates, 1996.

    Weir BS, Cockerham CC:
        Complete characterization of disequilibrium at two loci;
        in Feldman MW (ed): Mathematical Evolutionary Theory.
        Princeton, NJ: Princeton University Press, 1989.
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{snpgdsLDMat}}, \code{\link{snpgdsLDpruning}}
}

\examples{
# open an example dataset (HapMap)
genofile <- snpgdsOpen(snpgdsExampleFileName())

snp1 <- read.gdsn(index.gdsn(genofile, "genotype"), start=c(1,1), count=c(1,-1))
snp2 <- read.gdsn(index.gdsn(genofile, "genotype"), start=c(2,1), count=c(1,-1))

snpgdsLDpair(snp1, snp2, method = "composite")
snpgdsLDpair(snp1, snp2, method = "r")
snpgdsLDpair(snp1, snp2, method = "dprime")
snpgdsLDpair(snp1, snp2, method = "corr")

# close the genotype file
snpgdsClose(genofile)
}

\keyword{GDS}
\keyword{GWAS}
