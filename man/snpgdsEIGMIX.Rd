\name{snpgdsEIGMIX}
\alias{snpgdsEIGMIX}
\alias{plot.snpgdsEigMixClass}
\title{
    Eigen-analysis on SNP genotype data
}
\description{
    Eigen-analysis on IBD matrix based SNP genotypes.
}
\usage{
snpgdsEIGMIX(gdsobj, sample.id=NULL, snp.id=NULL, autosome.only=TRUE,
    remove.monosnp=TRUE, maf=NaN, missing.rate=NaN, num.thread=1L,
    eigen.cnt=32L, diagadj=TRUE, ibdmat=FALSE, verbose=TRUE)
\method{plot}{snpgdsEigMixClass}(x, eig=c(1L,2L), ...)
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
    \item{remove.monosnp}{if \code{TRUE}, remove monomorphic SNPs}
    \item{maf}{to use the SNPs with ">= maf" only; if NaN, no MAF threshold}
    \item{missing.rate}{to use the SNPs with "<= missing.rate" only; if NaN,
        no missing threshold}
    \item{num.thread}{the number of (CPU) cores used; if \code{NA}, detect
        the number of cores automatically}
    \item{eigen.cnt}{output the number of eigenvectors; if \code{eigen.cnt < 0},
        returns all eigenvectors; if \code{eigen.cnt==0}, no eigen calculation}
    \item{diagadj}{\code{TRUE} for diagonal adjustment by default}
    \item{ibdmat}{if \code{TRUE}, returns the IBD matrix}
    \item{verbose}{if \code{TRUE}, show information}
    \item{x}{a \code{snpgdsEigMixClass} object}
    \item{eig}{indices of eigenvectors, like \code{1:2} or \code{1:4}}
    \item{...}{the arguments passed to or from other methods, like
        \code{pch}, \code{col}}
}
\value{
    Return a \code{snpgdsEigMixClass} object, and it is a list:
    \item{sample.id}{the sample ids used in the analysis}
    \item{snp.id}{the SNP ids used in the analysis}
    \item{eigenval}{eigenvalues}
    \item{eigenvect}{eigenvactors, "# of samples" x "eigen.cnt"}
    \item{afreq}{allele frequencies}
    \item{ibd}{the IBD matrix when \code{ibdmat=TRUE}}
    \item{diagadj}{the argument \code{diagadj}}
}
\references{
    Zheng X, Weir BS.
    Eigenanalysis on SNP Data with an Interpretation of Identity by Descent.
    Theoretical Population Biology. 2016 Feb;107:65-76.
    doi: 10.1016/j.tpb.2015.09.004
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{snpgdsAdmixProp}}, \code{\link{snpgdsAdmixPlot}},
    \code{\link{snpgdsPCA}}, \code{\link{snpgdsPCASNPLoading}}, 
    \code{\link{snpgdsPCASampLoading}} 
}

\examples{
# open an example dataset (HapMap)
genofile <- snpgdsOpen(snpgdsExampleFileName())

# get population information
#   or pop_code <- scan("pop.txt", what=character())
#   if it is stored in a text file "pop.txt"
pop_code <- read.gdsn(index.gdsn(genofile, "sample.annot/pop.group"))

# get sample id
samp.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

# run eigen-analysis
RV <- snpgdsEIGMIX(genofile)
RV

# eigenvalues
RV$eigenval

# make a data.frame
tab <- data.frame(sample.id = samp.id, pop = factor(pop_code),
    EV1 = RV$eigenvect[,1],    # the first eigenvector
    EV2 = RV$eigenvect[,2],    # the second eigenvector
    stringsAsFactors = FALSE)
head(tab)

# draw
plot(tab$EV2, tab$EV1, col=as.integer(tab$pop),
    xlab="eigenvector 2", ylab="eigenvector 1")
legend("topleft", legend=levels(tab$pop), pch="o", col=1:4)


# define groups
groups <- list(CEU = samp.id[pop_code == "CEU"],
    YRI = samp.id[pop_code == "YRI"],
    CHB = samp.id[is.element(pop_code, c("HCB", "JPT"))])

prop <- snpgdsAdmixProp(RV, groups=groups)

# draw
plot(prop[, "YRI"], prop[, "CEU"], col=as.integer(tab$pop),
    xlab = "Admixture Proportion from YRI",
    ylab = "Admixture Proportion from CEU")
abline(v=0, col="gray25", lty=2)
abline(h=0, col="gray25", lty=2)
abline(a=1, b=-1, col="gray25", lty=2)
legend("topright", legend=levels(tab$pop), pch="o", col=1:4)


# close the genotype file
snpgdsClose(genofile)
}

\keyword{GDS}
\keyword{GWAS}
