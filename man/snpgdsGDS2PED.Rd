\name{snpgdsGDS2PED}
\alias{snpgdsGDS2PED}
\title{Conversion from GDS to PED}
\description{
    Convert a GDS file to a PLINK ped file.
}
\usage{
snpgdsGDS2PED(gdsobj, ped.fn, sample.id=NULL, snp.id=NULL, use.snp.rsid=TRUE,
    format=c("A/G/C/T", "A/B", "1/2"), verbose=TRUE)
}
\arguments{
    \item{gdsobj}{a GDS file object (\code{\link[gdsfmt]{gds.class}})}
    \item{ped.fn}{the file name of output}
    \item{sample.id}{a vector of sample id specifying selected samples;
        if NULL, all samples are used}
    \item{snp.id}{a vector of snp id specifying selected SNPs; if NULL, all
        SNPs are used}
    \item{use.snp.rsid}{if \code{TRUE}, use "snp.rs.id" instead of "snp.id" if
        available}
    \item{format}{specify the coding: "A/G/C/T" -- allelic codes stored in
        "snp.allele" of the GDS file; "A/B" -- A and B codes; "1/2" -- 1 and
        2 codes}
    \item{verbose}{if TRUE, show information}
}
\value{
    None.
}
\details{
    GDS -- Genomic Data Structures, the extended file name used for storing
        genetic data, and the file format used in the \link{gdsfmt} package.

    PED -- the PLINK text ped format.
}
\references{
    Purcell S, Neale B, Todd-Brown K, Thomas L, Ferreira MAR, Bender D,
    Maller J, Sklar P, de Bakker PIW, Daly MJ & Sham PC. 2007.
    PLINK: a toolset for whole-genome association and
    population-based linkage analysis. American Journal of Human Genetics, 81.

    \url{http://corearray.sourceforge.net/}
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{snpgdsGDS2BED}}
}

\examples{
# open an example dataset (HapMap)
genofile <- snpgdsOpen(snpgdsExampleFileName())

snpset <- snpgdsSelectSNP(genofile, missing.rate=0.95)
snpgdsGDS2PED(genofile, ped.fn="tmp", snp.id=snpset)

# close the genotype file
snpgdsClose(genofile)
}

\keyword{GDS}
\keyword{GWAS}
