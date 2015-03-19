\name{snpgdsPED2GDS}
\alias{snpgdsPED2GDS}
\title{Conversion from PLINK PED to GDS}
\description{
    Convert a PLINK PED text file to a GDS file.
}
\usage{
snpgdsPED2GDS(ped.fn, map.fn, out.gdsfn, family=TRUE, snpfirstdim=FALSE,
    compress.annotation="ZIP_RA.max", compress.geno="", verbose=TRUE)
}
\arguments{
    \item{ped.fn}{the file name of PED file, genotype information}
    \item{map.fn}{the file name of MAP file}
    \item{out.gdsfn}{the output GDS file}
    \item{family}{if \code{TRUE}, to include family information in the
        sample annotation}
    \item{snpfirstdim}{if TRUE, genotypes are stored in the individual-major
        mode, (i.e, list all SNPs for the first individual, and then list all
        SNPs for the second individual, etc)}
    \item{compress.annotation}{the compression method for the GDS variables,
        except "genotype"; optional values are defined in the function
        \code{add.gdsn}}
    \item{compress.geno}{the compression method for "genotype"; optional
        values are defined in the function \code{add.gdsn}}
    \item{verbose}{if TRUE, show information}
}
\value{
    None.
}
\details{
    GDS -- Genomic Data Structures, the extended file name used for storing
        genetic data, and the file format is used in the \link{gdsfmt} package.

    PED -- PLINK PED format.
}
\references{
    Purcell S, Neale B, Todd-Brown K, Thomas L, Ferreira MAR, Bender D,
    Maller J, Sklar P, de Bakker PIW, Daly MJ & Sham PC. 2007.
    PLINK: a toolset for whole-genome association and population-based linkage
    analysis. American Journal of Human Genetics, 81.
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{snpgdsGDS2PED}}, \code{\link{snpgdsBED2GDS}},
    \code{\link{snpgdsGDS2BED}}
}

\examples{
# open
genofile <- snpgdsOpen(snpgdsExampleFileName())

snpgdsGDS2PED(genofile, "tmp")

# close
snpgdsClose(genofile)


# PED ==> GDS
snpgdsPED2GDS("tmp.ped", "tmp.map", "test.gds")


# delete the temporary file
unlink(c("tmp.ped", "tmp.map", "test.gds"), force=TRUE)
}

\keyword{GDS}
\keyword{GWAS}
