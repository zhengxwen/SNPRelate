\name{snpgdsPCASampLoading}
\alias{snpgdsPCASampLoading}
\title{
    Project individuals onto existing principal component axes
}
\description{
    To calculate the sample eigenvectors using the specified SNP loadings
}
\usage{
snpgdsPCASampLoading(loadobj, gdsobj, sample.id=NULL, num.thread=1L,
    verbose=TRUE)
}
\arguments{
    \item{loadobj}{a \code{snpgdsPCASNPLoadingClass} or
        \code{snpgdsEigMixSNPLoadingClass} object returned from
        \link{snpgdsPCASNPLoading}}
    \item{gdsobj}{an object of class \code{\link{SNPGDSFileClass}},
        a SNP GDS file}
    \item{sample.id}{a vector of sample id specifying selected samples;
        if NULL, all samples are used}
    \item{num.thread}{the number of CPU cores used}
    \item{verbose}{if TRUE, show information}
}
\details{
    The \code{sample.id} are usually different from the samples used in the
calculation of SNP loadings.
}
\value{
    Returns a \code{snpgdsPCAClass} object, and it is a list:
    \item{sample.id}{the sample ids used in the analysis}
    \item{snp.id}{the SNP ids used in the analysis}
    \item{eigenval}{eigenvalues}
    \item{eigenvect}{eigenvactors, ``# of samples'' x ``eigen.cnt''}
    \item{TraceXTX}{the trace of the genetic covariance matrix}
    \item{Bayesian}{whether use bayerisan normalization}

    Or returns a \code{snpgdsEigMixClass} object, and it is a list:
    \item{sample.id}{the sample ids used in the analysis}
    \item{snp.id}{the SNP ids used in the analysis}
    \item{eigenval}{eigenvalues}
    \item{eigenvect}{eigenvactors, ``# of samples'' x ``eigen.cnt''}
    \item{afreq}{allele frequencies}
}
\references{
    Patterson N, Price AL, Reich D (2006)
    Population structure and eigenanalysis. PLoS Genetics 2:e190.

    Zhu, X., Li, S., Cooper, R. S., and Elston, R. C. (2008).
    A unified association analysis approach for family and unrelated samples
    correcting for stratification. Am J Hum Genet, 82(2), 352-365.
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{snpgdsPCA}}, \code{\link{snpgdsPCACorr}},
    \code{\link{snpgdsPCASNPLoading}}
}

\examples{
# open an example dataset (HapMap)
genofile <- snpgdsOpen(snpgdsExampleFileName())

sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

# first PCA
pca <- snpgdsPCA(genofile, eigen.cnt=8)
snp_load <- snpgdsPCASNPLoading(pca, genofile)

# calculate sample eigenvectors from SNP loadings
samp_load <- snpgdsPCASampLoading(snp_load, genofile, sample.id=sample.id[1:100])

diff <- pca$eigenvect[1:100,] - samp_load$eigenvect
summary(c(diff))
# ~ ZERO


# combine eigenvectors
allpca <- list(
    sample.id = c(pca$sample.id, samp_load$sample.id),
    snp.id = pca$snp.id,
    eigenval = c(pca$eigenval, samp_load$eigenval),
    eigenvect = rbind(pca$eigenvect, samp_load$eigenvect),
    varprop = c(pca$varprop, samp_load$varprop),
    TraceXTX = pca$TraceXTX
)
class(allpca) <- "snpgdsPCAClass"
allpca


# close the genotype file
snpgdsClose(genofile)
}

\keyword{PCA}
\keyword{GDS}
\keyword{GWAS}
