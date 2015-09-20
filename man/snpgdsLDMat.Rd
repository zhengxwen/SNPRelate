\name{snpgdsLDMat}
\alias{snpgdsLDMat}
\title{
    Linkage Disequilibrium (LD) analysis
}
\description{
    Return a LD matrix for SNP pairs.
}
\usage{
snpgdsLDMat(gdsobj, sample.id=NULL, snp.id=NULL, slide=250L,
    method=c("composite", "r", "dprime", "corr", "cov"), mat.trim=FALSE,
    num.thread=1L, with.id=TRUE, verbose=TRUE)
}
\arguments{
    \item{gdsobj}{an object of class \code{\link{SNPGDSFileClass}},
        a SNP GDS file}
    \item{sample.id}{a vector of sample id specifying selected samples;
        if NULL, all samples are used}
    \item{snp.id}{a vector of snp id specifying selected SNPs;
        if NULL, all SNPs are used}
    \item{slide}{# of SNPs, the size of sliding window, see details}
    \item{method}{"composite", "r", "dprime", "corr", "cov", see details}
    \item{mat.trim}{if \code{TRUE}, trim the matrix when \code{slide > 0}:
        the function returns a \code{"num_slide x (n_snp - slide)"} matrix}
    \item{num.thread}{the number of (CPU) cores used; if \code{NA}, detect
        the number of cores automatically}
    \item{with.id}{if \code{TRUE}, the returned value with \code{sample.id} and
        \code{sample.id}}
    \item{verbose}{if TRUE, show information}
}
\details{
    Four methods can be used to calculate linkage disequilibrium values:
"composite" for LD composite measure, "r" for R coefficient (by EM algorithm
assuming HWE, it could be negative), "dprime" for D', and "corr" for
correlation coefficient. The method "corr" is equivalent to "composite",
when SNP genotypes are coded as: 0 -- BB, 1 -- AB, 2 -- AA.

    If \code{slide <= 0}, the function returns a n-by-n LD matrix where the
value of i row and j column is LD of i and j SNPs. If \code{slide > 0}, it
returns a m-by-n LD matrix where n is the number of SNPs, m is the size of
sliding window, and the value of i row and j column is LD of j and j+i SNPs.
}
\value{
    Return a list:
    \item{sample.id}{the sample ids used in the analysis}
    \item{snp.id}{the SNP ids used in the analysis}
    \item{LD}{a matrix of LD values}
    \item{slide}{the size of sliding window}
}
\references{
    Weir B: Inferences about linkage disequilibrium. Biometrics 1979;
    35: 235-254.

    Weir B: Genetic Data Analysis II. Sunderland, MA: Sinauer Associates, 1996.

    Weir BS, Cockerham CC:
    Complete characterization of disequilibrium at two loci;
    in Feldman MW (ed): Mathematical Evolutionary Theory.
    Princeton, NJ: Princeton University Press, 1989.
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{snpgdsLDpair}}, \code{\link{snpgdsLDpruning}}
}

\examples{
# open an example dataset (HapMap)
genofile <- snpgdsOpen(snpgdsExampleFileName())

# chromosome 15
snpset <- read.gdsn(index.gdsn(genofile, "snp.id"))[
    read.gdsn(index.gdsn(genofile, "snp.chromosome")) == 15]
length(snpset)

# LD matrix without sliding window
ld.noslide <- snpgdsLDMat(genofile, snp.id=snpset, slide=-1, method="composite")
# plot
image(t(ld.noslide$LD^2), col=terrain.colors(16))

# LD matrix with a sliding window
ld.slide <- snpgdsLDMat(genofile, snp.id=snpset, method="composite")
# plot
image(t(ld.slide$LD^2), col=terrain.colors(16))


# close the genotype file
snpgdsClose(genofile)
}

\keyword{GDS}
\keyword{GWAS}
