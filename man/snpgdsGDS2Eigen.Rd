\name{snpgdsGDS2Eigen}
\alias{snpgdsGDS2Eigen}
\title{Conversion from GDS to Eigen (EIGENSTRAT)}
\description{
    Convert a GDS file to an EIGENSTRAT file.
}
\usage{
snpgdsGDS2Eigen(gdsobj, eigen.fn, sample.id=NULL, snp.id=NULL, verbose=TRUE)
}
\arguments{
    \item{gdsobj}{an object of class \code{\link{SNPGDSFileClass}},
        a SNP GDS file}
    \item{eigen.fn}{the file name of EIGENSTRAT}
    \item{sample.id}{a vector of sample id specifying selected samples;
        if NULL, all samples are used}
    \item{snp.id}{a vector of snp id specifying selected SNPs; if NULL,
        all SNPs are used}
    \item{verbose}{if TRUE, show information}
}
\value{
    None.
}
\details{
    GDS -- Genomic Data Structures, the extended file name used for storing
genetic data, and the file format used in the \link{gdsfmt} package.

    Eigen -- the text format used in EIGENSTRAT.
}
\references{
    Patterson N, Price AL, Reich D (2006) Population structure and
        eigenanalysis. PLoS Genetics 2:e190.

    Price AL, Patterson NJ, Plenge RM, Weinblatt ME, Shadick NA, Reich D (2006)
        Principal components analysis corrects for stratification in
        genome-wide association studies.
        Nat Genet. 38, 904-909.

    \url{http://corearray.sourceforge.net/}
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{snpgdsGDS2PED}}
}

\examples{
# open an example dataset (HapMap)
genofile <- snpgdsOpen(snpgdsExampleFileName())

snpset <- snpgdsSelectSNP(genofile, missing.rate=0.95)
snpgdsGDS2Eigen(genofile, eigen.fn="tmpeigen", snp.id=snpset)

# close the genotype file
snpgdsClose(genofile)
}

\keyword{GDS}
\keyword{GWAS}
