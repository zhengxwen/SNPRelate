\name{snpgdsCreateGeno}
\alias{snpgdsCreateGeno}
\title{
    Create a SNP genotype dataset from a matrix
}
\description{
    To create a GDS file of genotypes from a matrix.
}
\usage{
snpgdsCreateGeno(gds.fn, genmat, sample.id=NULL, snp.id=NULL, snp.rs.id=NULL,
    snp.chromosome=NULL, snp.position=NULL, snp.allele=NULL, snpfirstdim=TRUE,
    compress.annotation="ZIP.max", compress.geno="", other.vars=NULL)
}
\arguments{
    \item{gds.fn}{the file name of gds}
    \item{genmat}{a matrix of genotypes}
    \item{sample.id}{the sample ids, which should be unique}
    \item{snp.id}{the SNP ids, which should be unique}
    \item{snp.rs.id}{the rs ids for SNPs, which can be not unique}
    \item{snp.chromosome}{the chromosome indices}
    \item{snp.position}{the SNP positions in basepair}
    \item{snp.allele}{the reference/non-reference alleles}
    \item{snpfirstdim}{if TRUE, genotypes are stored in the individual-major
        mode, (i.e, list all SNPs for the first individual, and then list
        all SNPs for the second individual, etc)}
    \item{compress.annotation}{the compression method for the variables except
        \code{genotype}}
    \item{compress.geno}{the compression method for the variable
        \code{genotype}}
    \item{other.vars}{a list object storing other variables}
}
\details{
    There are possible values stored in the variable \code{genmat}: 0, 1, 2
and other values. ``0'' indicates two B alleles, ``1'' indicates one A allele
and one B allele, ``2'' indicates two A alleles, and other values indicate a
missing genotype.

    If \code{snpfirstdim} is \code{TRUE}, then \code{genmat} should be ``# of
SNPs X # of samples''; if \code{snpfirstdim} is \code{FALSE}, then
\code{genmat} should be ``# of samples X # of SNPs''.

    The typical variables specified in \code{other.vars} are ``sample.annot''
and ``snp.annot'', which are data.frame objects.
}
\value{
    None.
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{snpgdsCreateGenoSet}}, \code{\link{snpgdsCombineGeno}}
}
\examples{
# load data
data(hapmap_geno)

# create a gds file
with(hapmap_geno, snpgdsCreateGeno("test.gds", genmat=genotype,
    sample.id=sample.id, snp.id=snp.id, snp.chromosome=snp.chromosome,
    snp.position=snp.position, snp.allele=snp.allele, snpfirstdim=TRUE))

# open the gds file
genofile <- snpgdsOpen("test.gds")

RV <- snpgdsPCA(genofile)
plot(RV$eigenvect[,2], RV$eigenvect[,1], xlab="PC 2", ylab="PC 1")

# close the file
snpgdsClose(genofile)
}

\keyword{GDS}
\keyword{GWAS}
