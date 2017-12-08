\name{snpgdsMergeGRM}
\alias{snpgdsMergeGRM}
\title{
    Merge Multiple Genetic Relationship Matrices (GRM)
}
\description{
    Combine multiple genetic relationship matrices with weighted averaging.
}
\usage{
snpgdsMergeGRM(filelist, out.fn=NULL, out.prec=c("double", "single"),
    out.compress="LZMA_RA", weight=NULL, verbose=TRUE)
}
\arguments{
    \item{filelist}{a character vector, list of GDS file names}
    \item{out.fn}{NULL, return a GRM object; or characters, the output GDS
        file name}
    \item{out.prec}{double or single precision for storage}
    \item{out.compress}{the compression method used in the GDS file}
    \item{weight}{NULL, weights proportional to the numbers of SNPs;
        a numeric vector, or a logical vector (FALSE for excluding some GRMs
        with a negative weight, weights proportional to the numbers of SNPs)}
    \item{verbose}{if \code{TRUE}, show information}
}
\details{
    The final GRM is the weighted averaged matrix combining multiple GRMs.
The merged GRM may not be identical to the GRM calculated using full SNPs, due
to missing genotypes or the internal weighting strategy of the specified GRM
calculation.
}
\value{
    None or a GRM object if \code{out.fn=NULL}.
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{snpgdsGRM}}
}

\examples{
# open an example dataset (HapMap)
genofile <- snpgdsOpen(snpgdsExampleFileName())

snpid <- read.gdsn(index.gdsn(genofile, "snp.id"))
snpid <- snpid[snpgdsSNPRateFreq(genofile)$MissingRate == 0]

grm <- snpgdsGRM(genofile, snp.id=snpid, method="GCTA")


# save two GRMs
set1 <- grm$snp.id[1:(length(grm$snp.id)/2)]
set2 <- setdiff(grm$snp.id, set1)
snpgdsGRM(genofile, method="GCTA", snp.id=set1, out.fn="tmp1.gds")
snpgdsGRM(genofile, method="GCTA", snp.id=set2, out.fn="tmp2.gds")

# merge GRMs and export to a new GDS file
snpgdsMergeGRM(c("tmp1.gds", "tmp2.gds"), "tmp.gds")

# return the GRM
grm2 <- snpgdsMergeGRM(c("tmp1.gds", "tmp2.gds"))


# check
f <- openfn.gds("tmp.gds")
m <- read.gdsn(index.gdsn(f, "grm"))
closefn.gds(f)

summary(c(m - grm$grm))  # ~zero
summary(c(m - grm2$grm))  # zero


# close the file
snpgdsClose(genofile)

# delete the temporary file
unlink(c("tmp1.gds", "tmp2.gds", "tmp.gds"), force=TRUE)
}

\keyword{GDS}
\keyword{GWAS}
