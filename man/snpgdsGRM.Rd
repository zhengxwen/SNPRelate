\name{snpgdsGRM}
\alias{snpgdsGRM}
\title{
    Genetic Relationship Matrix (GRM) for SNP genotype data
}
\description{
    Eigen-analysis on the Genetic Relationship Matrix (GRM) for SNP
genotype data.
}
\usage{
snpgdsGRM(gdsobj, sample.id=NULL, snp.id=NULL,
    autosome.only=TRUE, remove.monosnp=TRUE, maf=NaN, missing.rate=NaN,
    num.thread=1L, eigen.cnt=32, need.GRM=FALSE, GRM.only=FALSE,
    verbose=TRUE)
}
\arguments{
    \item{gdsobj}{an object of class \code{\link{SNPGDSFileClass}},
        a SNP GDS file}
    \item{sample.id}{a vector of sample id specifying selected samples;
        if NULL, all samples are used}
    \item{snp.id}{a vector of snp id specifying selected SNPs; if NULL,
        all SNPs are used}
    \item{autosome.only}{if TRUE, use autosomal SNPs only}
    \item{remove.monosnp}{if TRUE, remove monomorphic SNPs}
    \item{maf}{to use the SNPs with ">= maf" only; if NaN, no MAF threshold}
    \item{missing.rate}{to use the SNPs with "<= missing.rate" only; if NaN,
        no missing threshold}
    \item{num.thread}{the number of (CPU) cores used; if \code{NA}, detect
        the number of cores automatically}
    \item{eigen.cnt}{output the number of eigenvectors; if eigen.cnt <= 0, then
        return all eigenvectors}
    \item{need.GRM}{if TRUE, return the GRM matrix}
    \item{GRM.only}{return the IBD matrix only, do not compute the
        eigenvalues and eigenvectors}
    \item{verbose}{if TRUE, show information}
}
\value{
    Return a \code{snpgdsGRMClass} object, and it is a list:
    \item{sample.id}{the sample ids used in the analysis}
    \item{snp.id}{the SNP ids used in the analysis}
    \item{eigenval}{eigenvalues}
    \item{eigenvect}{eigenvactors, "# of samples" x "eigen.cnt"}
    \item{ibdmat}{the IBD matrix}
}
\references{
    Yang, J., Lee, S. H., Goddard, M. E. & Visscher, P. M.
    GCTA: a tool for genome-wide complex trait analysis.
    American journal of human genetics 88, 76-82 (2011).
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{snpgdsPCA}}, \code{\link{snpgdsEIGMIX}}
}

\examples{
# open an example dataset (HapMap)
genofile <- snpgdsOpen(snpgdsExampleFileName())

# run PCA
RV <- snpgdsGRM(genofile)

####  there is no population information  ####

# make a data.frame
tab <- data.frame(sample.id = RV$sample.id,
    EV1 = RV$eigenvect[,1],    # the first eigenvector
    EV2 = RV$eigenvect[,2],    # the second eigenvector
    stringsAsFactors = FALSE)
head(tab)
#   sample.id         EV1         EV2
#   ...               ...         ...

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
#   ...               ...         ...

# draw
plot(tab$EV2, tab$EV1, col=as.integer(tab$pop),
    xlab="eigenvector 2", ylab="eigenvector 1")
legend("topleft", legend=levels(tab$pop), pch="o", col=1:4)


# close the file
snpgdsClose(genofile)
}

\keyword{GDS}
\keyword{GWAS}
