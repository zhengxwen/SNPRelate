\name{snpgdsPCASNPLoading}
\alias{snpgdsPCASNPLoading}
\title{
    SNP loadings in principal component analysis
}
\description{
    To calculate the SNP loadings in Principal Component Analysis
}
\usage{
snpgdsPCASNPLoading(pcaobj, gdsobj, num.thread=1L, verbose=TRUE)
}
\arguments{
    \item{pcaobj}{the snpgdsPCAClass object returned from the function
        \link{snpgdsPCA}}
    \item{gdsobj}{an object of class \code{\link{SNPGDSFileClass}},
        a SNP GDS file}
    \item{num.thread}{the number of (CPU) cores used; if \code{NA}, detect
        the number of cores automatically}
    \item{verbose}{if TRUE, show information}
}
\details{
    Calculate the SNP loadings (or SNP eigenvectors) from the principal
component analysis conducted in \code{snpgdsPCA}.
}
\value{
    Return a \code{snpgdsPCASNPLoading} object, which is a list:
    \item{sample.id}{the sample ids used in the analysis}
    \item{snp.id}{the SNP ids used in the analysis}
    \item{eigenval}{eigenvalues}
    \item{snploading}{the SNP loadings, or SNP eigenvectors}
    \item{TraceXTX}{the trace of the genetic covariance matrix}
    \item{Bayesian}{whether use bayerisan normalization}
    \item{avefreq}{the allele frequency used in \code{snpgdsPCA}}
    \item{scale}{internal parameter}
}

\references{
    Patterson N, Price AL, Reich D (2006)
    Population structure and eigenanalysis. PLoS Genetics 2:e190.

    Price AL, Patterson NJ, Plenge RM, Weinblatt ME, Shadick NA, Reich D (2006)
    Principal components analysis corrects for stratification in genome-wide
    association studies. Nat Genet. 38, 904-909.

    Zhu, X., Li, S., Cooper, R. S., and Elston, R. C. (2008).
    A unified association analysis approach for family and unrelated samples
    correcting for stratification. Am J Hum Genet, 82(2), 352-365.
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{snpgdsPCA}}, \code{\link{snpgdsPCASampLoading}},
    \code{\link{snpgdsPCACorr}}
}

\examples{
# open an example dataset (HapMap)
genofile <- snpgdsOpen(snpgdsExampleFileName())

PCARV <- snpgdsPCA(genofile, eigen.cnt=8)
SnpLoad <- snpgdsPCASNPLoading(PCARV, genofile)

names(SnpLoad)
# [1] "sample.id"  "snp.id"     "eigenval"   "snploading" "TraceXTX"
# [6] "Bayesian"   "avefreq"    "scale"
dim(SnpLoad$snploading)
# [1]     8 8722

plot(SnpLoad$snploading[1,], type="h", ylab="PC 1")

# close the genotype file
snpgdsClose(genofile)
}

\keyword{PCA}
\keyword{GDS}
\keyword{GWAS}
