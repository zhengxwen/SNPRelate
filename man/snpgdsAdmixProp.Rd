\name{snpgdsAdmixProp}
\alias{snpgdsAdmixProp}
\title{
    Estimate ancestral proportions from the eigen-analysis
}
\description{
    Estimate ancestral (admixture) proportions based on the eigen-analysis.
}
\usage{
snpgdsAdmixProp(eigobj, groups, bound=FALSE)
}
\arguments{
    \item{eigobj}{an object of \code{snpgdsEigMixClass} from
        \code{\link{snpgdsEIGMIX}}, or an object of \code{snpgdsPCAClass}
        from \code{\link{snpgdsPCA}}}
    \item{groups}{a list of sample IDs, such like \code{groups = list(
        CEU = c("NA0101", "NA1022", ...), YRI = c("NAxxxx", ...),
        Asia = c("NA1234", ...))}}
    \item{bound}{if \code{TRUE}, the estimates are bounded in [0, 1], and
        the sum of proportions is one; \code{bound=FALSE} for unbiased
        estimates}
}
\details{
    The minor allele frequency and missing rate for each SNP passed in
\code{snp.id} are calculated over all the samples in \code{sample.id}.
}
\value{
    Return a matrix of ancestral proportions with rows for study individuals
(\code{rownames()} is sample ID).
}

\references{
    Zheng X, Weir BS.
    Eigenanalysis on SNP Data with an Interpretation of Identity by Descent.
    Theoretical Population Biology. 2015 Oct 23. pii: S0040-5809(15)00089-1.
    doi: 10.1016/j.tpb.2015.09.004. [Epub ahead of print]
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{snpgdsEIGMIX}}, \code{\link{snpgdsPCA}},
    \code{\link{snpgdsAdmixPlot}}
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
legend("bottomleft", legend=levels(tab$pop), pch="o", col=1:4)


# define groups
groups <- list(CEU = samp.id[pop_code == "CEU"],
    YRI = samp.id[pop_code == "YRI"],
    CHB = samp.id[is.element(pop_code, c("HCB", "JPT"))])

prop <- snpgdsAdmixProp(RV, groups=groups)
head(prop)

# draw
plot(prop[, "YRI"], prop[, "CEU"], col=as.integer(tab$pop),
    xlab = "Admixture Proportion from YRI",
    ylab = "Admixture Proportion from CEU")
abline(v=0, col="gray25", lty=2)
abline(h=0, col="gray25", lty=2)
abline(a=1, b=-1, col="gray25", lty=2)
legend("topright", legend=levels(tab$pop), pch="o", col=1:4)


# draw
snpgdsAdmixPlot(prop, group=pop_code)


# close the genotype file
snpgdsClose(genofile)
}

\keyword{GDS}
\keyword{GWAS}
