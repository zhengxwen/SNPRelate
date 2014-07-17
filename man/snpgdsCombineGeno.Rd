\name{snpgdsCombineGeno}
\alias{snpgdsCombineGeno}
\title{
    Merge SNP datasets
}
\description{
    To merge GDS files of SNP genotypes into a single GDS file
}
\usage{
snpgdsCombineGeno(gds.fn, out.fn, sample.id=NULL, snpobj=NULL,
    name.prefix=NULL, snpfirstdim=TRUE, compress.annotation="ZIP.MAX",
    compress.geno="", other.vars=NULL, verbose=TRUE)
}
\arguments{
    \item{gds.fn}{a list of SNP GDS files to be merged}
    \item{out.fn}{the name of output GDS file}
    \item{sample.id}{NULL, or a list. If it is a list, specify sample ids
        for each SNP GDS file}
    \item{snpobj}{specify a \code{\link{snpgdsSNPListClass}} object, used for
        strand switch; if NULL, the strand information of the first SNP GDS
        file is used}
    \item{name.prefix}{NULL, a character vector (added to sample ids for each
        GDS file)}
    \item{snpfirstdim}{if TRUE, genotypes are stored in the individual-major
        mode, (i.e, list all SNPs for the first individual, and then list all
        SNPs for the second individual, etc)}
    \item{compress.annotation}{the compression method for the variables except
        \code{genotype}}
    \item{compress.geno}{the compression method for the variable
        \code{genotype}}
    \item{other.vars}{a list object storing other variables}
    \item{verbose}{if TRUE, show information}
}
\details{
    The typical variables specified in \code{other.vars} are ``sample.annot''
and ``snp.annot'', which are data.frame objects.
}
\value{
    None.
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{snpgdsCreateGeno}}, \code{\link{snpgdsCreateGenoSet}}
}

\examples{
# get the file name of a gds file
fn <- snpgdsExampleFileName()

# combine
snpgdsCombineGeno(c(fn, fn), "test.gds")

snpgdsSummary("test.gds")
}

\keyword{GDS}
\keyword{GWAS}
