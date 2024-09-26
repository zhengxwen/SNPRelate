\name{snpgdsGRM}
\alias{snpgdsGRM}
\title{
    Genetic Relationship Matrix (GRM) for SNP genotype data
}
\description{
    Calculate Genetic Relationship Matrix (GRM) using SNP genotype data.
}
\usage{
snpgdsGRM(gdsobj, sample.id=NULL, snp.id=NULL,
    autosome.only=TRUE, remove.monosnp=TRUE, maf=NaN, missing.rate=NaN,
    method=c("GCTA", "Eigenstrat", "EIGMIX", "Weighted", "Corr", "IndivBeta"),
    num.thread=1L, useMatrix=FALSE, out.fn=NULL, out.prec=c("double", "single"),
    out.compress="LZMA_RA", with.id=TRUE, verbose=TRUE)
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
    \item{method}{"GCTA" -- genetic relationship matrix defined in CGTA;
        "Eigenstrat" -- genetic covariance matrix in EIGENSTRAT;
        "EIGMIX" -- two times coancestry matrix defined in Zheng&Weir (2016),
        "Weighted" -- weighted GCTA, as the same as "EIGMIX",
        "Corr" -- Scaled GCTA GRM (dividing each i,j element by the
            product of the square root of the i,i and j,j elements),
        "IndivBeta" -- two times individual beta estimate relative to the
            minimum of beta; see details}
    \item{num.thread}{the number of (CPU) cores used; if \code{NA}, detect
        the number of cores automatically}
    \item{useMatrix}{if \code{TRUE}, use \code{Matrix::dspMatrix} to store
        the output square matrix to save memory}
    \item{out.fn}{NULL for no GDS output, or a file name}
    \item{out.prec}{double or single precision for storage}
    \item{out.compress}{the compression method for storing the GRM matrix in
        the GDS file}
    \item{with.id}{if \code{TRUE}, the returned value with \code{sample.id}
        and \code{sample.id}}
    \item{verbose}{if \code{TRUE}, show information}
}
\details{
    "GCTA": the genetic relationship matrix in GCTA is defined as
$G_ij = avg_l [(g_il - 2*p_l*(g_jl - 2*p_l) / 2*p_l*(1 - p_l)]$ for
individuals i,j and locus l;

    "Eigenstrat": the genetic covariance matrix in EIGENSTRAT
$G_ij = avg_l [(g_il - 2*p_l)*(g_jl - 2*p_l) / 2*p_l*(1 - p_l)]$ for
individuals i,j and locus l; the missing genotype is imputed by the dosage
mean of that locus.

    "EIGMIX" / "Weighted": it is the same as `2 * snpgdsEIGMIX(, ibdmat=TRUE,
diagadj=FALSE)$ibd`:
$G_ij = [sum_l (g_il - 2*p_l)*(g_jl - 2*p_l)] / [sum_l 2*p_l*(1 - p_l)]$
for individuals i,j and locus l;

    "IndivBeta": `beta = snpgdsIndivBeta(, inbreeding=TRUE)`
(Weir&Goudet, 2017), and beta-based GRM is
$grm_ij = 2 * (beta_ij - beta_min) / (1 - beta_min)$ for $i!=j$,
$grm_ij = 1 + (beta_i - beta_min) / (1 - beta_min)$ for $i=j$.
It is relative to the minimum value of beta estimates.
}
\value{
    Return a list if \code{with.id = TRUE}:
        \item{sample.id}{the sample ids used in the analysis}
        \item{snp.id}{the SNP ids used in the analysis}
        \item{method}{characters, the method used}
        \item{grm}{the genetic relationship matrix; different methods might
            have different meanings and interpretation for estimates}

    If \code{with.id = FALSE}, this function returns the genetic relationship
    matrix (GRM) without sample and SNP IDs.
}
\references{
    Patterson, N., Price, A. L. & Reich, D.
    Population structure and eigenanalysis. PLoS Genet. 2, e190 (2006).

    Yang, J., Lee, S. H., Goddard, M. E. & Visscher, P. M.
    GCTA: a tool for genome-wide complex trait analysis.
    American journal of human genetics 88, 76-82 (2011).

    Zheng X, Weir BS.
    Eigenanalysis on SNP Data with an Interpretation of Identity by Descent.
    Theoretical Population Biology. 2016 Feb;107:65-76.
    doi: 10.1016/j.tpb.2015.09.004

    Weir BS, Zheng X.
    SNPs and SNVs in Forensic Science.
    Forensic Science International: Genetics Supplement Series. 2015.
    doi:10.1016/j.fsigss.2015.09.106

    Weir BS, Goudet J. A Unified Characterization of Population Structure
    and Relatedness. Genetics. 2017 Aug;206(4):2085-2103.
    doi: 10.1534/genetics.116.198424.
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{snpgdsPCA}}, \code{\link{snpgdsEIGMIX}},
    \code{\link{snpgdsIndivBeta}},
    \code{\link{snpgdsIndInb}}, \code{\link{snpgdsFst}},
    \code{\link{snpgdsMergeGRM}}
}

\examples{
# open an example dataset (HapMap)
genofile <- snpgdsOpen(snpgdsExampleFileName())

rv <- snpgdsGRM(genofile, method="GCTA")
eig <- eigen(rv$grm)  # Eigen-decomposition

# output to a GDS file
snpgdsGRM(genofile, method="GCTA", out.fn="test.gds")


pop <- factor(read.gdsn(index.gdsn(genofile, "sample.annot/pop.group")))
plot(eig$vectors[,1], eig$vectors[,2], col=pop)
legend("topleft", legend=levels(pop), pch=19, col=1:4)


# close the file
snpgdsClose(genofile)

# delete the temporary file
unlink("test.gds", force=TRUE)
}

\keyword{GDS}
\keyword{GWAS}
