\name{snpgdsBED2GDS}
\alias{snpgdsBED2GDS}
\title{Conversion from PLINK BED to GDS}
\description{
    Convert a PLINK binary ped file to a GDS file.
}
\usage{
snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, out.gdsfn, family=FALSE,
    snpfirstdim=NA, compress.annotation="LZMA_RA", compress.geno="",
    option=NULL, cvt.chr=c("int", "char"), cvt.snpid=c("auto", "int"),
    verbose=TRUE)
}
\arguments{
    \item{bed.fn}{the file name of binary file, genotype information}
    \item{fam.fn}{the file name of first six columns of \code{".ped"}; if it is
        missing, ".fam" is added to \code{bed.fn}}
    \item{bim.fn}{the file name of extended MAP file: two extra columns =
        allele names; if it is missing, ".bim" is added to \code{bim.fn}}
    \item{out.gdsfn}{the output file name of GDS file}
    \item{family}{if \code{TRUE}, to include family information in the
        sample annotation}
    \item{snpfirstdim}{if TRUE, genotypes are stored in the individual-major
        mode, (i.e, list all SNPs for the first individual, and then list all
        SNPs for the second individual, etc); \code{NA}, the dimension is
        determined by the BED file}
    \item{compress.annotation}{the compression method for the GDS variables,
        except "genotype"; optional values are defined in the function
        \code{add.gdsn}}
    \item{compress.geno}{the compression method for "genotype"; optional
        values are defined in the function \code{add.gdsn}}
    \item{option}{\code{NULL} or an object from \code{\link{snpgdsOption}},
        see details}
    \item{cvt.chr}{\code{"int"} -- chromosome code in the GDS file is integer;
        \code{"char"} -- chromosome code in the GDS file is character}
    \item{cvt.snpid}{\code{"int"} -- to create an integer \code{snp.id}
        starting from 1; \code{"auto"} -- if SNP IDs in the PLINK file are not
        unique, to create an an integer \code{snp.id}, otherwise to use SNP
        IDs for \code{snp.id}}
    \item{verbose}{if TRUE, show information}
}
\value{
    Return the file name of GDS format with an absolute path.
}
\details{
    GDS -- Genomic Data Structures, the extended file name used for storing
        genetic data, and the file format is used in the \link{gdsfmt} package.

    BED -- the PLINK binary ped format.

    The user could use \code{option} to specify the range of code for autosomes.
For humans there are 22 autosomes (from 1 to 22), but dogs have 38 autosomes.
Note that the default settings are used for humans. The user could call
\code{option = snpgdsOption(autosome.end=38)} for importing the BED file of dog.
It also allow define new chromosome coding, e.g.,
\code{option = snpgdsOption(Z=27)}.
}
\references{
    Purcell S, Neale B, Todd-Brown K, Thomas L, Ferreira MAR, Bender D,
    Maller J, Sklar P, de Bakker PIW, Daly MJ & Sham PC. 2007.
    PLINK: a toolset for whole-genome association and population-based linkage
    analysis. American Journal of Human Genetics, 81.
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{snpgdsOption}}, \code{\link{snpgdsPED2GDS}},
    \code{\link{snpgdsGDS2PED}}
}

\examples{
# PLINK BED files
bed.fn <- system.file("extdata", "plinkhapmap.bed.gz", package="SNPRelate")
fam.fn <- system.file("extdata", "plinkhapmap.fam.gz", package="SNPRelate")
bim.fn <- system.file("extdata", "plinkhapmap.bim.gz", package="SNPRelate")

# convert
snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, "HapMap.gds")

# open
genofile <- snpgdsOpen("HapMap.gds")
genofile

# close
snpgdsClose(genofile)


# delete the temporary file
unlink("HapMap.gds", force=TRUE)
}

\keyword{GDS}
\keyword{GWAS}
