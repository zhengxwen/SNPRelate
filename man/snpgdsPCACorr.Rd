\name{snpgdsPCACorr}
\alias{snpgdsPCACorr}
\title{
    PC-correlated SNPs in principal component analysis
}
\description{
    To calculate the SNP correlations between eigenvactors and SNP genotypes
}
\usage{
snpgdsPCACorr(pcaobj, gdsobj, snp.id=NULL, eig.which=NULL, num.thread=1L,
    with.id=TRUE, outgds=NULL, verbose=TRUE)
}
\arguments{
    \item{pcaobj}{a \code{snpgdsPCAClass} object returned from the function
        \link{snpgdsPCA}, a \code{snpgdsEigMixClass} from
        \link{snpgdsEIGMIX}, or an eigenvector matrix with row names
        (sample id)}
    \item{gdsobj}{an object of class \code{\link{SNPGDSFileClass}},
        a SNP GDS file}
    \item{snp.id}{a vector of snp id specifying selected SNPs;
        if NULL, all SNPs are used}
    \item{eig.which}{a vector of integers, to specify which eigenvectors
        to be used}
    \item{num.thread}{the number of (CPU) cores used; if \code{NA}, detect
        the number of cores automatically}
    \item{with.id}{if \code{TRUE}, the returned value with \code{sample.id} and
        \code{sample.id}}
    \item{outgds}{\code{NULL} or a character of file name for exporting
        correlations to a GDS file, see details}
    \item{verbose}{if TRUE, show information}
}
\details{
    If an output file name is specified via \code{outgds}, "sample.id",
"snp.id" and "correlation" will be stored in the GDS file. The GDS node
"correlation" is a matrix of correlation coefficients, and it is stored with
the format of packed real number ("packedreal16" preserving 4 digits, 0.0001
is the smallest number greater zero, see \link{add.gdsn}).
}
\value{
    Return a list if \code{outgds=NULL},
    \item{sample.id}{the sample ids used in the analysis}
    \item{snp.id}{the SNP ids used in the analysis}
    \item{snpcorr}{a matrix of correlation coefficients,
        "# of eigenvectors" x "# of SNPs"}
}
\references{
    Patterson N, Price AL, Reich D (2006)
    Population structure and eigenanalysis. PLoS Genetics 2:e190.
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{snpgdsPCA}}, \code{\link{snpgdsPCASampLoading}},
    \code{\link{snpgdsPCASNPLoading}}
}

\examples{
# open an example dataset (HapMap)
genofile <- snpgdsOpen(snpgdsExampleFileName())
# get chromosome index
chr <- read.gdsn(index.gdsn(genofile, "snp.chromosome"))

pca <- snpgdsPCA(genofile)
cr <- snpgdsPCACorr(pca, genofile, eig.which=1:4)
plot(abs(cr$snpcorr[3,]), xlab="SNP Index", ylab="PC 3", col=chr)


# output to a gds file if limited memory
snpgdsPCACorr(pca, genofile, eig.which=1:4, outgds="test.gds")

(f <- openfn.gds("test.gds"))
m <- read.gdsn(index.gdsn(f, "correlation"))
closefn.gds(f)

# check
summary(c(m - cr$snpcorr))  # should < 1e-4


# close the file
snpgdsClose(genofile)

# delete the temporary file
unlink("test.gds", force=TRUE)
}

\keyword{PCA}
\keyword{GDS}
\keyword{GWAS}
