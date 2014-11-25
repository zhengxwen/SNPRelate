\name{snpgdsGDS2BED}
\alias{snpgdsGDS2BED}
\title{Conversion from GDS to PLINK BED}
\description{
    Convert a GDS file to a PLINK binary ped file.
}
\usage{
snpgdsGDS2BED(gdsobj, bed.fn, sample.id=NULL, snp.id=NULL, snpfirstdim=NULL,
    verbose=TRUE)
}
\arguments{
    \item{gdsobj}{an object of class \code{\link{SNPGDSFileClass}},
        a SNP GDS file; or characters, the file name of GDS}
    \item{bed.fn}{the file name of output}
    \item{sample.id}{a vector of sample id specifying selected samples;
        if NULL, all samples are used}
    \item{snp.id}{a vector of snp id specifying selected SNPs;
        if NULL, all SNPs are used}
    \item{snpfirstdim}{if TRUE, genotypes are stored in the individual-major
        mode, (i.e, list all SNPs for the first individual, and then list all
        SNPs for the second individual, etc); if NULL, determine
        automatically}
    \item{verbose}{if TRUE, show information}
}
\value{
    None.
}
\details{
    GDS -- Genomic Data Structures, the extended file name used for storing
genetic data, and the file format used in the \link{gdsfmt} package.

    BED -- the PLINK binary ped format.
}
\references{
    Purcell S, Neale B, Todd-Brown K, Thomas L, Ferreira MAR, Bender D,
    Maller J, Sklar P, de Bakker PIW, Daly MJ & Sham PC. 2007.
    PLINK: a toolset for whole-genome association and population-based linkage
    analysis. American Journal of Human Genetics, 81.

    \url{http://corearray.sourceforge.net/}
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{snpgdsBED2GDS}}, \code{\link{snpgdsGDS2PED}}
}

\examples{
# open an example dataset (HapMap)
genofile <- snpgdsOpen(snpgdsExampleFileName())

snpset <- snpgdsSelectSNP(genofile, missing.rate=0.95)
snpgdsGDS2BED(genofile, bed.fn="test", snp.id=snpset)

# close the genotype file
snpgdsClose(genofile)
}

\keyword{GDS}
\keyword{GWAS}
