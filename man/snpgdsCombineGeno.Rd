\name{snpgdsCombineGeno}
\alias{snpgdsCombineGeno}
\title{
    Merge SNP datasets
}
\description{
    To merge GDS files of SNP genotypes into a single GDS file
}
\usage{
snpgdsCombineGeno(gds.fn, out.fn, method=c("position", "exact"),
    compress.annotation="ZIP_RA.MAX", compress.geno="",
    same.strand=FALSE, snpfirstdim=FALSE, verbose=TRUE)
}
\arguments{
    \item{gds.fn}{a character vector of GDS file names to be merged}
    \item{out.fn}{the name of output GDS file}
    \item{method}{\code{"exact"}: matching by all snp.id, chromosomes, positions
        and alleles; \code{"position"}: matching by chromosomes and positions}
    \item{compress.annotation}{the compression method for the variables except
        \code{genotype}}
    \item{compress.geno}{the compression method for the variable
        \code{genotype}}
    \item{same.strand}{if TRUE, assuming the alleles on the same strand}
    \item{snpfirstdim}{if TRUE, genotypes are stored in the individual-major
        mode, (i.e, list all SNPs for the first individual, and then list all
        SNPs for the second individual, etc)}
    \item{verbose}{if TRUE, show information}
}
\details{
    This function calls \code{\link{snpgdsSNPListIntersect}} internally to
determine the common SNPs. Allele definitions are taken from the first GDS file.
}
\value{
    None.
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{snpgdsCreateGeno}}, \code{\link{snpgdsCreateGenoSet}},
    \code{\link{snpgdsSNPList}}, \code{\link{snpgdsSNPListIntersect}}
}

\examples{
# get the file name of a gds file
fn <- snpgdsExampleFileName()

f <- snpgdsOpen(fn)
samp_id <- read.gdsn(index.gdsn(f, "sample.id"))
snp_id <- read.gdsn(index.gdsn(f, "snp.id"))
geno <- read.gdsn(index.gdsn(f, "genotype"), start=c(1,1), count=c(-1, 3000))
snpgdsClose(f)


# split the GDS file with different samples
snpgdsCreateGenoSet(fn, "t1.gds", sample.id=samp_id[1:10],
    snp.id=snp_id[1:3000])
snpgdsCreateGenoSet(fn, "t2.gds", sample.id=samp_id[11:30],
    snp.id=snp_id[1:3000])

# combine with different samples
snpgdsCombineGeno(c("t1.gds", "t2.gds"), "test.gds", same.strand=TRUE)
f <- snpgdsOpen("test.gds")
g <- read.gdsn(index.gdsn(f, "genotype"))
snpgdsClose(f)

identical(geno[1:30, ], g)  # TRUE


# split the GDS file with different SNPs
snpgdsCreateGenoSet(fn, "t1.gds", snp.id=snp_id[1:100])
snpgdsCreateGenoSet(fn, "t2.gds", snp.id=snp_id[101:300])

# combine with different SNPs
snpgdsCombineGeno(c("t1.gds", "t2.gds"), "test.gds")
f <- snpgdsOpen("test.gds")
g <- read.gdsn(index.gdsn(f, "genotype"))
snpgdsClose(f)

identical(geno[, 1:300], g)  # TRUE


# delete the temporary files
unlink(c("t1.gds", "t2.gds", "t3.gds", "t4.gds", "test.gds"), force=TRUE)
}

\keyword{GDS}
\keyword{GWAS}
