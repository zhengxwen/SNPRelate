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
    snpfirstdim=NULL, compress.annotation="ZIP_RA.max", compress.geno="",
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
# +    [  ] *
# |--+ sample.id   { VStr8 279 ZIP(29.9%), 679B }
# |--+ snp.id   { Int32 9088 ZIP(34.8%), 12.3K }
# |--+ snp.rs.id   { VStr8 9088 ZIP(40.1%), 36.2K }
# |--+ snp.position   { Int32 9088 ZIP(94.7%), 33.6K }
# |--+ snp.chromosome   { UInt8 9088 ZIP(0.94%), 85B } *
# |--+ snp.allele   { VStr8 9088 ZIP(11.3%), 4.0K }
# |--+ genotype   { Bit2 279x9088, 619.0K } *
# \--+ sample.annot   [ data.frame ] *
#    |--+ family.id   { VStr8 279 ZIP(34.4%), 514B }
#    |--+ father.id   { VStr8 279 ZIP(31.5%), 220B }
#    |--+ mother.id   { VStr8 279 ZIP(30.9%), 214B }
#    |--+ sex   { VStr8 279 ZIP(17.0%), 95B }
#    \--+ pop.group   { VStr8 279 ZIP(6.18%), 69B }

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
# +    [  ] *
# |--+ sample.id   { Str8 279 ZIP_ra(31.2%), 715B }
# |--+ snp.id   { Int32 6547 ZIP_ra(34.9%), 8.9K }
# |--+ snp.rs.id   { Str8 6547 ZIP_ra(41.5%), 27.1K }
# |--+ snp.position   { Int32 6547 ZIP_ra(94.9%), 24.3K }
# |--+ snp.chromosome   { Int32 6547 ZIP_ra(0.45%), 124B }
# |--+ snp.allele   { Str8 6547 ZIP_ra(11.5%), 3.0K }
# \--+ genotype   { Bit2 279x6547, 446.0K } *

# close the file
snpgdsClose(gfile)


unlink("test.gds", force=TRUE)
}

\keyword{GDS}
\keyword{GWAS}
