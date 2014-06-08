\name{snpgdsBED2GDS}
\alias{snpgdsBED2GDS}
\title{Conversion from BED to GDS}
\description{
	Convert a PLINK binary ped file to a GDS file.
}
\usage{
snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, out.gdsfn, compress.annotation="ZIP.max", verbose=TRUE)
}
\arguments{
	\item{bed.fn}{the file name of binary file, genotype information}
	\item{fam.fn}{the file name of first six columns of *.ped}
	\item{bim.fn}{the file name of extended MAP file: two extra cols = allele names}
	\item{out.gdsfn}{the output gds file}
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

	BED -- the PLINK binary ped format.
}
\references{
	Purcell S, Neale B, Todd-Brown K, Thomas L, Ferreira MAR, Bender D, Maller J, Sklar P,
	de Bakker PIW, Daly MJ & Sham PC. 2007. PLINK: a toolset for whole-genome association and
	population-based linkage analysis. American Journal of Human Genetics, 81.

	\url{http://corearray.sourceforge.net/}
}
\author{Xiuwen Zheng \email{zhengx@u.washington.edu}}
\seealso{
	\code{\link{snpgdsBED2GDS}}, \code{\link{snpgdsGDS2PED}}
}

\examples{
# PLINK BED files
bed.fn <- system.file("extdata", "plinkhapmap.bed", package="SNPRelate")
bim.fn <- system.file("extdata", "plinkhapmap.bim", package="SNPRelate")
fam.fn <- system.file("extdata", "plinkhapmap.fam", package="SNPRelate")

# convert
snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, "HapMap.gds")

# open
genofile <- openfn.gds("HapMap.gds")
genofile

# close
closefn.gds(genofile)
}

\keyword{gds}
\keyword{GWAS}
