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
	snpfirstdim=NULL, compress.annotation="ZIP.max", compress.geno="", verbose=TRUE)
}
\arguments{
	\item{src.fn}{the file name of a specified GDS file}
	\item{dest.fn}{the file name of output GDS file}
	\item{sample.id}{a vector of sample id specifying selected samples; if NULL, all samples are used}
	\item{snp.id}{a vector of snp id specifying selected SNPs; if NULL, all SNPs are used}
	\item{snpfirstdim}{if TRUE, genotypes are stored in the individual-major mode,
		(i.e, list all SNPs for the first individual, and then list all SNPs for
		the second individual, etc)}
	\item{compress.annotation}{the compression method for the variables except \code{genotype}}
	\item{compress.geno}{the compression method for the variable \code{genotype}}
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
(genofile <- openfn.gds(snpgdsExampleFileName()))
# + 	[   ]
# |--+ sample.id	[ dFStr8 279 ZIP(23.10%) ]
# |--+ snp.id	[ dInt32 45440 ZIP(34.60%) ]
# |--+ snp.rs.id	[ dFStr8 45440 ZIP(41.05%) ]
# |--+ snp.position	[ dInt32 45440 ZIP(94.25%) ]
# |--+ snp.chromosome	[ dInt32 45440 ZIP(0.14%) ]
# |--+ snp.allele	[ dFStr8 45440 ZIP(13.36%) ]
# |--+ genotype	[ dBit2 45440x279 ]
# |--+ sample.annot	[   ] *
# |  |--+ sex	[ dFStr8 279 ZIP(28.32%) ] *
# |  |--+ pop.group	[ dFStr8 279 ZIP(7.89%) ] *

snpset <- unlist(snpgdsLDpruning(genofile))
length(snpset)
# 43518

# close the genotype file
closefn.gds(genofile)

snpgdsCreateGenoSet(snpgdsExampleFileName(), "test.gds", snp.id=snpset)

####################################################
# check

(gfile <- openfn.gds("test.gds"))
# + 	[   ]
# |--+ sample.id	[ dFStr8 279 ZIP(23.10%) ]
# |--+ snp.id	[ dInt32 43518 ZIP(34.61%) ]
# |--+ snp.rs.id	[ dFStr8 43518 ZIP(41.09%) ]
# |--+ snp.position	[ dInt32 43518 ZIP(94.16%) ]
# |--+ snp.chromosome	[ dInt32 43518 ZIP(0.14%) ]
# |--+ snp.allele	[ dFStr8 43518 ZIP(13.37%) ]
# |--+ genotype	[ dBit2 43518x279 ] *

closefn.gds(gfile)
}

\keyword{gds}
\keyword{GWAS}
\keyword{LD}
