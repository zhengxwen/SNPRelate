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
    method=c("IBS", "GVH", "HVG"), type=c("per.pair", "per.snp", "matrix"),
    with.id=TRUE, verbose=TRUE)
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
    \item{method}{\code{"IBS"} -- identity-by-state score, \code{"GVH"}
        or \code{"HVG"}, see Details}
    \item{type}{\code{"per.pair"}, \code{"per.snp"} or \code{"matrix"},
        see Value}
    \item{with.id}{if \code{TRUE}, returns "sample.id" and "snp.id"; see Value}
    \item{verbose}{if TRUE, show information}
}
\details{
    \tabular{ccccccc}{
        Patient (sample1.id) \tab Coded Genotype \tab
        Donor (sample2.id)   \tab Coded Genotype \tab
        IBS \tab GVH \tab HVG \cr

        AA \tab 0 \tab AA \tab 0 \tab 2 \tab 0 \tab 0 \cr
        AA \tab 0 \tab AB \tab 1 \tab 1 \tab 0 \tab 1 \cr
        AA \tab 0 \tab BB \tab 2 \tab 0 \tab 2 \tab 2 \cr
        AB \tab 1 \tab AA \tab 0 \tab 1 \tab 1 \tab 0 \cr
        AB \tab 1 \tab AB \tab 1 \tab 2 \tab 0 \tab 0 \cr
        AB \tab 1 \tab BB \tab 2 \tab 1 \tab 1 \tab 0 \cr
        BB \tab 2 \tab AA \tab 0 \tab 0 \tab 2 \tab 2 \cr
        BB \tab 2 \tab AB \tab 1 \tab 1 \tab 0 \tab 1 \cr
        BB \tab 2 \tab BB \tab 2 \tab 2 \tab 0 \tab 0 \cr
    }
}
\value{
    Return a list:
    \item{sample.id}{the sample ids used in the analysis,
        if \code{with.id=TRUE}}
    \item{snp.id}{the SNP ids used in the analysis, if \code{with.id=TRUE}}
    \item{score}{a matrix of genotype score:
        if \code{type="per.pair"}, a \code{# of pairs}-by-3 matrix with the
        first column for average scores, the second column for standard
        deviation and the third column for the valid number of SNPs;
        if \code{type="per.snp"}, a 3-by-\code{# of SNPs} matrix with the first
        row for average scores, the second row for standard deviation and
        the third row for the valid number of individual pairs;
        if \code{type="matrix"}, a \code{# of pairs}-by-\code{# of SNPs}
        matrix with rows for pairs of individuals}
}
\references{
    Warren, E. H., Zhang, X. C., Li, S., Fan, W., Storer, B. E., Chien, J. W.,
Boeckh, M. J., et al. (2012). Effect of MHC and non-MHC donor/recipient
genetic disparity on the outcome of allogeneic HCT. Blood, 120(14), 2796-806.
doi:10.1182/blood-2012-04-347286
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
    method="IBS", type="per.pair")
names(z1)
z1$score


# calculate average genotype scores
z2 <- snpgdsPairScore(genofile, offspring.id, father.id, snp.id=selsnp,
    method="IBS", type="per.snp")
names(z2)
mean(z2$score["Avg",])
mean(z2$score["SD",])

plot(z2$score["Avg",], pch=20, cex=0.75, xlab="SNP Index", ylab="IBS score")


# calculate a matrix of genotype scores over samples and SNPs
z3 <- snpgdsPairScore(genofile, offspring.id, father.id, snp.id=selsnp,
    method="IBS", type="matrix")
dim(z3$score)


# close the file
snpgdsClose(genofile)
}

\keyword{GDS}
\keyword{GWAS}
