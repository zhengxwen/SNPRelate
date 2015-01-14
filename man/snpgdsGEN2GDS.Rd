\name{snpgdsGEN2GDS}
\alias{snpgdsGEN2GDS}
\title{Conversion from Oxford GEN format to GDS}
\description{
    Convert an Oxford GEN file (text format) to a GDS file.
}
\usage{
snpgdsGEN2GDS(gen.fn, sample.fn, out.fn, chr.code=NULL,
    call.threshold=0.9, version=c(">=2.0", "<=1.1.5"),
    snpfirstdim=FALSE, compress.annotation="ZIP_RA.max", compress.geno="",
    verbose=TRUE)
}
\arguments{
    \item{gen.fn}{the file name of Oxford GEN text file(s), it could be a
        vector indicate merging all files}
    \item{sample.fn}{the file name of sample annotation}
    \item{out.fn}{the output GDS file}
    \item{chr.code}{a vector of chromosome code according to \code{gen.fn},
        indicating chromosomes. It could be either numeric or character-type}
    \item{call.threshold}{the threshold to determine missing genotypes}
    \item{version}{either \code{">=2.0"} or \code{"<=1.1.5"}, see details}
    \item{snpfirstdim}{if TRUE, genotypes are stored in the individual-major
        mode, (i.e, list all SNPs for the first individual, and then list all
        SNPs for the second individual, etc)}
    \item{compress.annotation}{the compression mode of the nodes stored,
        except "genotype"; optional value are defined in the function of
        \code{add.gdsn}}
    \item{compress.geno}{the data compression mode for the variable "genotype",
        optional values are defined in the function of \code{add.gdsn}}
    \item{verbose}{if TRUE, show information}
}
\value{
    Return the file name of GDS format with an absolute path.
}
\details{
    GDS -- Genomic Data Structures, the extended file name used for storing
        genetic data, and the file format is used in the \link{gdsfmt} package.

    NOTE : the sample file format (\code{sample.fn}) has changed with the
release of SNPTEST v2. Specifically, the way in which covariates and phenotypes
are coded on the second line of the header file has changed. \code{version} has
to be specified, and the function uses ">=2.0" by default.
}
\references{
    \url{http://www.stats.ox.ac.uk/~marchini/software/gwas/file_format.html}
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{snpgdsBED2GDS}}, \code{\link{snpgdsVCF2GDS}}
}

\examples{
cat("running snpgdsGEN2GDS ...\n")
\dontrun{
snpgdsGEN2GDS("test.gen", "test.sample", "output.gds", chr.code=1)
}
}

\keyword{GDS}
\keyword{GWAS}
