\name{snpgdsVCF2GDS}
\alias{snpgdsVCF2GDS}
\title{Conversion from VCF to GDS}
\description{
	Extract SNP data from a VCF sequence file
}
\usage{
snpgdsVCF2GDS(vcf.fn, outfn.gds, nblock=1024, compress.annotation="ZIP.max",
	verbose=TRUE)
}
\arguments{
	\item{vcf.fn}{the file name of VCF format}
	\item{outfn.gds}{the output gds file}
	\item{nblock}{the buffer lines}
	\item{compress.annotation}{the compression flag of the nodes stored, except
		"genotype"; the string value is defined in the function of \code{add.gdsn}}
	\item{verbose}{if TRUE, show information}
}
\value{
	None.
}
\details{
	GDS -- Genomic Data Structures, the extended file name used for storing genetic data,
		and the file format used in the \link{gdsfmt} package.

	VCF -- The Variant Call Format (VCF), which is a generic format for storing DNA
		polymorphism data such as SNPs, insertions, deletions and structural variants,
		together with rich annotations.
}
\references{
	The variant call format and VCFtools.
	Danecek P, Auton A, Abecasis G, Albers CA, Banks E, DePristo MA, Handsaker RE,
	Lunter G, Marth GT, Sherry ST, McVean G, Durbin R; 1000 Genomes Project Analysis Group.
	Bioinformatics. 2011 Aug 1;27(15):2156-8. Epub 2011 Jun 7.

	\url{http://corearray.sourceforge.net/}
}
\author{Xiuwen Zheng \email{zhengx@u.washington.edu}}
\seealso{
	\code{\link{snpgdsBED2GDS}}, \code{\link{snpgdsGDS2PED}}
}

\examples{
# The VCF file
vcf.fn <- system.file("extdata", "sequence.vcf", package="SNPRelate")
readLines(vcf.fn)

snpgdsVCF2GDS(vcf.fn, "test.gds", verbose=TRUE)
snpgdsSummary("test.gds")


# open "test.gds"
(genofile <- openfn.gds("test.gds"))

read.gdsn(index.gdsn(genofile, "sample.id"))
read.gdsn(index.gdsn(genofile, "snp.rs.id"))
read.gdsn(index.gdsn(genofile, "genotype"))

# close the genotype file
closefn.gds(genofile)
}

\keyword{gds}
\keyword{GWAS}
