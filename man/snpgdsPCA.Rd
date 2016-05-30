\name{snpgdsPCA}
\alias{snpgdsPCA}
\title{
    Principal Component Analysis (PCA) on SNP genotype data
}
\description{
    To calculate the eigenvectors and eigenvalues for principal component
analysis in GWAS.
}
\usage{
snpgdsPCA(gdsobj, sample.id=NULL, snp.id=NULL,
    autosome.only=TRUE, remove.monosnp=TRUE, maf=NaN, missing.rate=NaN,
    algorithm=c("exact", "randomized"),
    eigen.cnt=ifelse(identical(algorithm, "randomized"), 16L, 32L),
    num.thread=1L, bayesian=FALSE, need.genmat=FALSE,
    genmat.only=FALSE, eigen.method=c("DSPEVX", "DSPEV"),
    aux.dim=eigen.cnt*2L, iter.num=10L, verbose=TRUE)
}
\arguments{
    \item{gdsobj}{an object of class \code{\link{SNPGDSFileClass}},
        a SNP GDS file}
    \item{sample.id}{a vector of sample id specifying selected samples;
        if NULL, all samples are used}
    \item{snp.id}{a vector of snp id specifying selected SNPs; if NULL,
        all SNPs are used}
    \item{autosome.only}{if \code{TRUE}, use autosomal SNPs only; if it is a
        numeric or character value, keep SNPs according to the specified
        chromosome}
    \item{remove.monosnp}{if TRUE, remove monomorphic SNPs}
    \item{maf}{to use the SNPs with ">= maf" only; if NaN, no MAF threshold}
    \item{missing.rate}{to use the SNPs with "<= missing.rate" only; if NaN,
        no missing threshold}
    \item{eigen.cnt}{output the number of eigenvectors; if eigen.cnt <= 0, then
        return all eigenvectors}
    \item{algorithm}{"exact", traditional exact calculation; "randomized",
        fast PCA with randomized algorithm introduced in Galinsky et al. 2016}
    \item{num.thread}{the number of (CPU) cores used; if \code{NA}, detect
        the number of cores automatically}
    \item{bayesian}{if TRUE, use bayesian normalization}
    \item{need.genmat}{if TRUE, return the genetic covariance matrix}
    \item{genmat.only}{return the genetic covariance matrix only, do not
        compute the eigenvalues and eigenvectors}
    \item{eigen.method}{"DSPEVX" -- compute the top \code{eigen.cnt}
        eigenvalues and eigenvectors using LAPACK::DSPEVX; "DSPEV" -- to be
        compatible with SNPRelate_1.1.6 or earlier, using LAPACK::DSPEV;
        "DSPEVX" is significantly faster than "DSPEV" if only top principal
        components are of interest}
    \item{aux.dim}{auxiliary dimension used in fast randomized algorithm}
    \item{iter.num}{iteration number used in fast randomized algorithm}
    \item{verbose}{if TRUE, show information}
}
\details{
    The minor allele frequency and missing rate for each SNP passed in
\code{snp.id} are calculated over all the samples in \code{sample.id}.
}
\value{
    Return a \code{snpgdsPCAClass} object, and it is a list:
    \item{sample.id}{the sample ids used in the analysis}
    \item{snp.id}{the SNP ids used in the analysis}
    \item{eigenval}{eigenvalues}
    \item{eigenvect}{eigenvactors, "# of samples" x "eigen.cnt"}
    \item{varprop}{variance proportion for each principal component}
    \item{TraceXTX}{the trace of the genetic covariance matrix}
    \item{Bayesian}{whether use bayerisan normalization}
    \item{genmat}{the genetic covariance matrix}
}
\references{
    Patterson N, Price AL, Reich D. Population structure and eigenanalysis.
PLoS Genet. 2006 Dec;2(12):e190.

	Galinsky KJ, Bhatia G, Loh PR, Georgiev S, Mukherjee S, Patterson NJ,
Price AL. Fast Principal-Component Analysis Reveals Convergent Evolution of
ADH1B in Europe and East Asia. Am J Hum Genet. 2016 Mar 3;98(3):456-72.
doi: 10.1016/j.ajhg.2015.12.022. Epub 2016 Feb 25.
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{snpgdsPCACorr}}, \code{\link{snpgdsPCASampLoading}},
    \code{\link{snpgdsPCASNPLoading}}
}

\examples{
# open an example dataset (HapMap)
genofile <- snpgdsOpen(snpgdsExampleFileName())

# run PCA
RV <- snpgdsPCA(genofile)

# eigenvalues
head(RV$eigenval)

# variance proportion (%)
head(round(RV$varprop*100, 2))
# [1] 12.23  5.84  1.01  0.95  0.84  0.74


####  there is no population information  ####

# make a data.frame
tab <- data.frame(sample.id = RV$sample.id,
    EV1 = RV$eigenvect[,1],    # the first eigenvector
    EV2 = RV$eigenvect[,2],    # the second eigenvector
    stringsAsFactors = FALSE)
head(tab)
#   sample.id         EV1         EV2
# 1   NA19152 -0.08411287 -0.01226860
# 2   NA19139 -0.08360644 -0.01085849
# 3   NA18912 -0.08110808 -0.01184524
# 4   NA19160 -0.08680864 -0.01447106
# 5   NA07034  0.03109761  0.07709255
# 6   NA07055  0.03228450  0.08155730

# draw
plot(tab$EV2, tab$EV1, xlab="eigenvector 2", ylab="eigenvector 1")



####  there are population information  ####

# get population information
#   or pop_code <- scan("pop.txt", what=character())
#   if it is stored in a text file "pop.txt"
pop_code <- read.gdsn(index.gdsn(genofile, "sample.annot/pop.group"))

# get sample id
samp.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

# assume the order of sample IDs is as the same as population codes
cbind(samp.id, pop_code)
#        samp.id       pop_code
#   [1,] "NA19152"     "YRI"   
#   [2,] "NA19139"     "YRI"   
#   [3,] "NA18912"     "YRI"   
#   [4,] "NA19160"     "YRI"   
#   [5,] "NA07034"     "CEU"   
#   ...

# make a data.frame
tab <- data.frame(sample.id = RV$sample.id,
    pop = factor(pop_code)[match(RV$sample.id, samp.id)],
    EV1 = RV$eigenvect[,1],    # the first eigenvector
    EV2 = RV$eigenvect[,2],    # the second eigenvector
    stringsAsFactors = FALSE)
head(tab)
#   sample.id pop         EV1         EV2
# 1   NA19152 YRI -0.08411287 -0.01226860
# 2   NA19139 YRI -0.08360644 -0.01085849
# 3   NA18912 YRI -0.08110808 -0.01184524
# 4   NA19160 YRI -0.08680864 -0.01447106
# 5   NA07034 CEU  0.03109761  0.07709255
# 6   NA07055 CEU  0.03228450  0.08155730

# draw
plot(tab$EV2, tab$EV1, col=as.integer(tab$pop),
    xlab="eigenvector 2", ylab="eigenvector 1")
legend("bottomright", legend=levels(tab$pop), pch="o", col=1:4)


# close the file
snpgdsClose(genofile)
}

\keyword{PCA}
\keyword{GDS}
\keyword{GWAS}
