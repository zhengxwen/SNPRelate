\name{snpgdsCreateGenoSet}
\alias{snpgdsCreateGenoSet}
\title{
    Create a SNP genotype dataset from a GDS file
}
\description{
    To create a GDS file of genotypes from a specified GDS file.
}
\usage{
snpgdsCreateGenoSet(src.fn, dest.fn, sample.id=NULL, snp.id=NULL,
    snpfirstdim=NULL, compress.annotation="ZIP.max", compress.geno="",
    verbose=TRUE)
}
\arguments{
    \item{src.fn}{the file name of a specified GDS file}
    \item{dest.fn}{the file name of output GDS file}
    \item{sample.id}{a vector of sample id specifying selected samples;
        if NULL, all samples are used}
    \item{snp.id}{a vector of snp id specifying selected SNPs; if NULL,
        all SNPs are used}
    \item{snpfirstdim}{if TRUE, genotypes are stored in the individual-major
        mode, (i.e, list all SNPs for the first individual, and then list all
        SNPs for the second individual, etc)}
    \item{compress.annotation}{the compression method for the variables except
        \code{genotype}}
    \item{compress.geno}{the compression method for the variable
        \code{genotype}}
    \item{verbose}{if TRUE, show information}
}
\value{
    None.
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{snpgdsCreateGeno}}, \code{\link{snpgdsCombineGeno}}
}

\examples{
# open an example dataset (HapMap)
(genofile <- snpgdsOpen(snpgdsExampleFileName()))
# +     [  ] *
# |--+ sample.id    { FStr8 279 ZIP(23.10%) }
# |--+ snp.id   { Int32 9088 ZIP(34.76%) }
# |--+ snp.rs.id    { FStr8 9088 ZIP(42.66%) }
# |--+ snp.position { Int32 9088 ZIP(94.73%) }
# |--+ snp.chromosome   { UInt8 9088 ZIP(0.94%) } *
# |--+ snp.allele   { FStr8 9088 ZIP(14.45%) }
# |--+ genotype { Bit2 9088x279 } *
# |--+ sample.annot [ data.frame ] *
# |  |--+ sample.id { FStr8 279 ZIP(23.10%) }
# |  |--+ family.id { FStr8 279 ZIP(28.37%) }
# |  |--+ geneva.id { Int32 279 ZIP(80.29%) }
# |  |--+ father.id { FStr8 279 ZIP(12.98%) }
# |  |--+ mother.id { FStr8 279 ZIP(12.86%) }
# |  |--+ plate.id  { FStr8 279 ZIP(1.29%) }
# |  |--+ sex   { FStr8 279 ZIP(28.32%) }
# |  |--+ pop.group { FStr8 279 ZIP(7.89%) }

set.seed(1000)
snpset <- unlist(snpgdsLDpruning(genofile))
length(snpset)
# 6547

# close the file
snpgdsClose(genofile)

snpgdsCreateGenoSet(snpgdsExampleFileName(), "test.gds", snp.id=snpset)

####################################################
# check

(gfile <- snpgdsOpen("test.gds"))
# +     [  ] *
# |--+ sample.id    { VStr8 279 ZIP(29.89%) }
# |--+ snp.id   { Int32 6547 ZIP(34.89%) }
# |--+ snp.rs.id    { VStr8 6547 ZIP(40.52%) }
# |--+ snp.position { Int32 6547 ZIP(94.85%) }
# |--+ snp.chromosome   { Int32 6547 ZIP(0.41%) }
# |--+ snp.allele   { VStr8 6547 ZIP(11.51%) }
# |--+ genotype { Bit2 6547x279 } *

# close the file
snpgdsClose(gfile)


unlink("test.gds", force=TRUE)
}

\keyword{GDS}
\keyword{GWAS}
