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
    method=c("GCTA", "Eigenstrat", "EIGMIX", "Weighted", "IndivBeta"),
    num.thread=1L, with.id=TRUE, verbose=TRUE)
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
        "EIGMIX" -- two times coancestry matrix defined in Zheng & Weir (2015),
        "Weighted" -- weighted GCTA, as the same as "EIGMIX",
        "IndivBeta" -- two times individual beta estimate; see details}
    \item{num.thread}{the number of (CPU) cores used; if \code{NA}, detect
        the number of cores automatically}
    \item{with.id}{if \code{TRUE}, the returned value with \code{sample.id}
        and \code{sample.id}}
    \item{verbose}{if \code{TRUE}, show information}
}
\details{
    "GCTA": the genetic relationship matrix in GCTA is defined as
G_{ij} = avg_l [(g_{il} - 2*p_l)*(g_{jl} - 2*p_l) / 2*p_l*(1 - p_l)] for
individuals i,j and locus l;

    "Eigenstrat": the genetic covariance matrix in EIGENSTRAT
G_{ij} = avg_l [(g_{il} - 2*p_l)*(g_{jl} - 2*p_l) / 2*p_l*(1 - p_l)] for
individuals i,j and locus l; it is identical as "GCTA" if there is no missing
genotype.

    "EIGMIX" / "Weighted": it is the same as `2 * snpgdsEIGMIX(, ibdmat=TRUE,
diagadj=FALSE)$ibd`:
G_{ij} = [sum_l (g_{il} - 2*p_l)*(g_{jl} - 2*p_l)] / [sum_l 2*p_l*(1 - p_l)]
for individuals i,j and locus l;

    "IndivBeta": it is the same as `2 * snpgdsIndivBeta(, with.id=FALSE)`.
}
\value{
    Return a list if \code{with.id = TRUE}:
        \item{sample.id}{the sample ids used in the analysis}
        \item{snp.id}{the SNP ids used in the analysis}
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
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{snpgdsPCA}}, \code{\link{snpgdsEIGMIX}},
    \code{\link{snpgdsIndivBeta}},
    \code{\link{snpgdsIndInb}}, \code{\link{snpgdsFst}}
}

\examples{
# open an example dataset (HapMap)
genofile <- snpgdsOpen(snpgdsExampleFileName())


rv <- snpgdsGRM(genofile, method="GCTA")
eig <- eigen(rv$grm)  # Eigen-decomposition

pop <- factor(read.gdsn(index.gdsn(genofile, "sample.annot/pop.group")))
plot(eig$vectors[,1], eig$vectors[,2], col=pop)
legend("topleft", legend=levels(pop), pch=19, col=1:4)


# close the file
snpgdsClose(genofile)
}

\keyword{GDS}
\keyword{GWAS}
