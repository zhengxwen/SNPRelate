\name{snpgdsCutTree}
\alias{snpgdsCutTree}
\title{
    Determine clusters of individuals
}
\description{
    To determine sub groups of individuals using a specified dendrogram from
hierarchical cluster analysis
}
\usage{
snpgdsCutTree(hc, z.threshold=15, outlier.n=5, n.perm = 5000, samp.group=NULL,
    col.outlier="red", col.list=NULL, pch.outlier=4, pch.list=NULL,
    label.H=FALSE, label.Z=TRUE, verbose=TRUE)
}
\arguments{
    \item{hc}{an object of \code{\link{snpgdsHCluster}}}
    \item{z.threshold}{the threshold of Z score to determine whether split the
        node or not}
    \item{outlier.n}{the cluster with size less than or equal to
        \code{outlier.n} is considered as outliers}
    \item{n.perm}{the times for permutation}
    \item{samp.group}{if \code{NULL}, determine groups by Z score; if a vector
        of factor, assign each individual in dendrogram with respect to
        \code{samp.group}}
    \item{col.outlier}{the color of outlier}
    \item{col.list}{the list of colors for different clusters}
    \item{pch.outlier}{plotting 'character' for outliers}
    \item{pch.list}{plotting 'character' for different clusters}
    \item{label.H}{if TRUE, plotting heights in a dendrogram}
    \item{label.Z}{if TRUE, plotting Z scores in a dendrogram}
    \item{verbose}{if TRUE, show information}
}
\details{
    The details will be described in future.
}
\value{
    Return a list:
    \item{sample.id}{the sample ids used in the analysis}
    \item{z.threshold}{the threshold of Z score to determine whether split
        the node or not}
    \item{outlier.n}{the cluster with size less than or equal to
        \code{outlier.n} is considered as outliers}
    \item{samp.order}{the order of samples in the dendrogram}
    \item{samp.group}{a vector of factor, indicating the group of each
        individual}
    \item{dmat}{a matrix of pairwise group dissimilarity}
    \item{dendrogram}{the dendrogram of individuals}
    \item{merge}{a data.frame of \code{(z, n1, n2)} describing each
        combination: \code{z}, the Z score; \code{n1}, the size of the first
        cluster; \code{n2}, the size of the second cluster}
    \item{clust.count}{the counts for clusters}
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{snpgdsHCluster}}, \code{\link{snpgdsDrawTree}},
    \code{\link{snpgdsIBS}}, \code{\link{snpgdsDiss}}
}

\examples{
# open an example dataset (HapMap)
genofile <- snpgdsOpen(snpgdsExampleFileName())

pop.group <- as.factor(read.gdsn(index.gdsn(
    genofile, "sample.annot/pop.group")))
pop.level <- levels(pop.group)

diss <- snpgdsDiss(genofile)
hc <- snpgdsHCluster(diss)

# close the genotype file
snpgdsClose(genofile)



###################################################################
# cluster individuals
#

set.seed(100)
rv <- snpgdsCutTree(hc, label.H=TRUE, label.Z=TRUE)

# the distribution of Z scores
snpgdsDrawTree(rv, type="z-score", main="HapMap Phase II")

# draw dendrogram
snpgdsDrawTree(rv, main="HapMap Phase II",
    edgePar=list(col=rgb(0.5,0.5,0.5, 0.75), t.col="black"))


###################################################################
# or cluster individuals by ethnic information
#

rv2 <- snpgdsCutTree(hc, samp.group=pop.group)

# cluster individuals by Z score, specifying 'clust.count'
snpgdsDrawTree(rv2, rv$clust.count, main="HapMap Phase II",
    edgePar = list(col=rgb(0.5,0.5,0.5, 0.75), t.col="black"),
    labels = c("YRI", "CHB/JPT", "CEU"), y.label=0.1)
legend("bottomleft", legend=levels(pop.group), col=1:nlevels(pop.group),
    pch=19, ncol=4, bg="white")



###################################################################
# zoom in ...
#

snpgdsDrawTree(rv2, rv$clust.count, dend.idx = c(1),
    main="HapMap Phase II -- YRI",
    edgePar=list(col=rgb(0.5,0.5,0.5, 0.75), t.col="black"),
    y.label.kinship=TRUE)

snpgdsDrawTree(rv2, rv$clust.count, dend.idx = c(2,2),
    main="HapMap Phase II -- CEU",
    edgePar=list(col=rgb(0.5,0.5,0.5, 0.75), t.col="black"),
    y.label.kinship=TRUE)

snpgdsDrawTree(rv2, rv$clust.count, dend.idx = c(2,1),
    main="HapMap Phase II -- CHB/JPT",
    edgePar=list(col=rgb(0.5,0.5,0.5, 0.75), t.col="black"),
    y.label.kinship=TRUE)
}

\keyword{GDS}
\keyword{GWAS}
