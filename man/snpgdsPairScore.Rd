\name{snpgdsPairScore}
\alias{snpgdsPairScore}
\title{
    Genotype Score for Pairs of Individuals
}
\description{
    Calculate the genotype score for pairs of individuals based on
identity-by-state (IBS) measure
}
\usage{
snpgdsPairScore(gdsobj, sample1.id, sample2.id, snp.id=NULL,
    method=c("IBS", "GVH", "HVG"), type=c("avg+sd", "matrix"), with.id=TRUE,
    verbose=TRUE)
}
\arguments{
    \item{gdsobj}{an object of class \code{\link{SNPGDSFileClass}},
        a SNP GDS file}
    \item{sample1.id}{a vector of sample id specifying selected samples;
        if NULL, all samples are used}
    \item{sample2.id}{a vector of sample id specifying selected samples;
        if NULL, all samples are used}
    \item{snp.id}{a vector of snp id specifying selected SNPs; if NULL,
        all SNPs are used}
    \item{method}{\code{"IBS"} -- identity-by-state score; \code{"GVH"} -- ;
        \code{"HVG"} -- }
    \item{type}{\code{"avg+sd"} -- ; \code{"matrix"} -- }
    \item{with.id}{if \code{TRUE}, returns "sample.id" and "snp.id"; see Value}
    \item{verbose}{if TRUE, show information}
}
\details{
    XX.
}
\value{
    Return a list:
    \item{sample.id}{the sample ids used in the analysis,
        if \code{with.id=TRUE}}
    \item{snp.id}{the SNP ids used in the analysis, if \code{with.id=TRUE}}
    \item{score}{a matrix of genotype score over SNPs:
        if \code{type="avg+sd"}, a 2-by-\code{# of SNPs} matrix with the first
        row for average scores and the second row for standard deviation;
        if \code{type="matrix"}, a \code{# of pairs}-by-\code{# of SNPs}
        matrix with rows for pairs of individuals}
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{snpgdsIBS}}
}

\examples{
# open an example dataset (HapMap)
genofile <- snpgdsOpen(snpgdsExampleFileName())

# autosomal SNPs
selsnp <- snpgdsSelectSNP(genofile, autosome.only=TRUE, remove.monosnp=FALSE)

# sample ID
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
father.id <- read.gdsn(index.gdsn(genofile, "sample.annot/father.id"))

offspring.id <- sample.id[father.id != ""]
father.id <- father.id[father.id != ""]


# calculate average genotype scores
z1 <- snpgdsPairScore(genofile, offspring.id, father.id, snp.id=selsnp,
    method="IBS", type="avg+sd")
names(z1)
mean(z1$score["Avg",])
mean(z1$score["SD",])

plot(z1$score["Avg",], pch=20, cex=0.75, xlab="SNP Index", ylab="IBS score")


# calculate a matrix of genotype scores over samples and SNPs
z2 <- snpgdsPairScore(genofile, offspring.id, father.id, snp.id=selsnp,
    method="IBS", type="matrix")
dim(z2$score)


# close the file
snpgdsClose(genofile)
}

\keyword{GDS}
\keyword{GWAS}
