\name{snpgdsGRM}
\alias{snpgdsGRM}
\title{
    Genetic Relationship Matrix (GRM) for SNP genotype data
}
\description{
    Calculate Genetic Relationship Matrix (GRM) using SNP genotype data.
}
\usage{
snpgdsGRM(gdsobj, sample.id=NULL, snp.id=NULL, autosome.only=TRUE,
    remove.monosnp=TRUE, maf=NaN, missing.rate=NaN,
    method=c("Visscher", "Eigenstrat", "EIGMIX", "W&Z15"), num.thread=1L,
    with.id=TRUE, verbose=TRUE)
}
\arguments{
    \item{gdsobj}{an object of class \code{\link{SNPGDSFileClass}},
        a SNP GDS file}
    \item{sample.id}{a vector of sample id specifying selected samples;
        if NULL, all samples are used}
    \item{snp.id}{a vector of snp id specifying selected SNPs; if NULL,
        all SNPs are used}
    \item{autosome.only}{if \code{TRUE}, use autosomal SNPs only; if it is a
        numeric or character value, keep SNPs according to the specified
        chromosome}
    \item{remove.monosnp}{if TRUE, remove monomorphic SNPs}
    \item{maf}{to use the SNPs with ">= maf" only; if NaN, no MAF threshold}
    \item{missing.rate}{to use the SNPs with "<= missing.rate" only; if NaN,
        no missing threshold}
    \item{method}{"Visscher" -- genetic relationship matrix defined in CGTA;
        "Eigenstrat" -- genetic covariance matrix in EIGENSTRAT;
        "EIGMIX" -- coancestry matrix defined in Zheng & Weir (2015),
        "W&Z15" -- individual beta estimator}
    \item{num.thread}{the number of (CPU) cores used; if \code{NA}, detect
        the number of cores automatically}
    \item{with.id}{if \code{TRUE}, the returned value with \code{sample.id}
        and \code{sample.id}}
    \item{verbose}{if \code{TRUE}, show information}
}
\value{
    Return a list if \code{with.id = TRUE}:
        \item{sample.id}{the sample ids used in the analysis}
        \item{snp.id}{the SNP ids used in the analysis}
        \item{grm}{the genetic relationship matrix}
    If \code{with.id = FALSE}, this function returns the genetic relationship
    matrix (GRM) without sample and SNP IDs.
}
\references{
    Yang, J., Lee, S. H., Goddard, M. E. & Visscher, P. M.
    GCTA: a tool for genome-wide complex trait analysis.
    American journal of human genetics 88, 76-82 (2011).

    Patterson, N., Price, A. L. & Reich, D.
    Population structure and eigenanalysis. PLoS Genet. 2, e190 (2006).

    Zheng X, Weir BS.
    Eigenanalysis on SNP Data with an Interpretation of Identity by Descent.
    2015. Theoretical Population Biology.

    Weir BS, Zheng X.
    SNPs and SNVs in Forensic Science. 2015.
    Forensic Science International: Genetics Supplement Series.
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{snpgdsIndInb}}, \code{\link{snpgdsFst}}
}

\examples{
# open an example dataset (HapMap)
genofile <- snpgdsOpen(snpgdsExampleFileName())


rv <- snpgdsGRM(genofile, method="Visscher")
eig <- eigen(rv$grm)  # Eigendecomposition

pop <- factor(read.gdsn(index.gdsn(genofile, "sample.annot/pop.group")))
plot(eig$vectors[,1], eig$vectors[,2], col=pop)
legend("topleft", legend=levels(pop), pch=19, col=1:4)


# close the file
snpgdsClose(genofile)
}

\keyword{GDS}
\keyword{GWAS}
