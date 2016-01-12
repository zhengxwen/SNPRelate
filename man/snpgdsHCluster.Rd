\name{snpgdsHCluster}
\alias{snpgdsHCluster}
\title{
    Hierarchical cluster analysis
}
\description{
    Perform hierarchical cluster analysis on the dissimilarity matrix.
}
\usage{
snpgdsHCluster(dist, sample.id=NULL, need.mat=TRUE, hang=0.25)
}
\arguments{
    \item{dist}{an object of "snpgdsDissClass" from \code{\link{snpgdsDiss}},
        an object of "snpgdsIBSClass" from \code{\link{snpgdsIBS}}, or
        a square matrix for dissimilarity}
    \item{sample.id}{to specify sample id, only work if dist is a matrix}
    \item{need.mat}{if TRUE, store the dissimilarity matrix in the result}
    \item{hang}{The fraction of the plot height by which labels should hang
        below the rest of the plot. A negative value will cause the labels to
        hang down from 0.}
}
\details{
    Call the function \code{\link{hclust}} to perform hierarchical cluster
analysis, using \code{method="average"}.
}
\value{
    Return a list (class "snpgdsHCClass"):
    \item{sample.id}{the sample ids used in the analysis}
    \item{hclust}{an object returned from \code{\link{hclust}}}
    \item{dendrogram}{}
    \item{dist}{the dissimilarity matrix, if \code{need.mat = TRUE}}
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{snpgdsIBS}}, \code{\link{snpgdsDiss}},
    \code{\link{snpgdsCutTree}}
}

\examples{
# open an example dataset (HapMap)
genofile <- snpgdsOpen(snpgdsExampleFileName())

pop.group <- read.gdsn(index.gdsn(genofile, "sample.annot/pop.group"))
pop.group <- as.factor(pop.group)
pop.level <- levels(pop.group)

diss <- snpgdsDiss(genofile)
hc <- snpgdsHCluster(diss)
rv <- snpgdsCutTree(hc)
rv

# call 'plot' to draw a dendrogram
plot(rv$dendrogram, leaflab="none", main="HapMap Phase II")


# the distribution of Z scores
snpgdsDrawTree(rv, type="z-score", main="HapMap Phase II")

# draw dendrogram
snpgdsDrawTree(rv, main="HapMap Phase II",
    edgePar=list(col=rgb(0.5,0.5,0.5, 0.75), t.col="black"))


# close the file
snpgdsClose(genofile)
}

\keyword{GDS}
\keyword{GWAS}
