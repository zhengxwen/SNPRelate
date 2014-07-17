\name{snpgdsDrawTree}
\alias{snpgdsDrawTree}
\title{
    Draw a dendrogram
}
\description{
    To draw a dendrogram or the distribution of Z scores
}
\usage{
snpgdsDrawTree(obj, clust.count=NULL, dend.idx=NULL,
    type=c("dendrogram", "z-score"), yaxis.height=TRUE, yaxis.kinship=TRUE,
    y.kinship.baseline=NaN, y.label.kinship=FALSE, outlier.n=NULL,
    shadow.col=c(rgb(0.5, 0.5, 0.5, 0.25), rgb(0.5, 0.5, 0.5, 0.05)),
    outlier.col=rgb(1, 0.50, 0.50, 0.5), leaflab="none",
    labels=NULL, y.label=0.2, ...)
}
\arguments{
    \item{obj}{an object returned by \code{\link{snpgdsCutTree}}}
    \item{clust.count}{the counts for clusters, drawing shadows}
    \item{dend.idx}{the index of sub tree, plot obj$dendrogram[[dend.idx]],
        or NULL for the whole tree}
    \item{type}{"dendrogram", draw a dendrogram; or "z-score", draw the
        distribution of Z score}
    \item{yaxis.height}{if TRUE, draw the left Y axis: height of tree}
    \item{yaxis.kinship}{if TRUE, draw the right Y axis: kinship coefficient}
    \item{y.kinship.baseline}{the baseline value of kinship; if NaN, it is the
        height of the first split from top in a dendrogram; only works when
        \code{yaxis.kinship = TRUE}}
    \item{y.label.kinship}{if TRUE, show 'PO/FS' etc on the right axis}
    \item{outlier.n}{the cluster with size less than or equal to
        \code{outlier.n} is considered as outliers; if \code{NULL}, let
        \code{outlier.n = obj$outlier.n}}
    \item{shadow.col}{two colors for shadow}
    \item{outlier.col}{the colors for outliers}
    \item{leaflab}{a string specifying how leaves are labeled. The default
        \code{"perpendicular"} write text vertically (by default).
        \code{"textlike"} writes text horizontally (in a rectangle), and
        \code{"none"} suppresses leaf labels.}
    \item{labels}{the legend for different regions}
    \item{y.label}{y positions of labels}
    \item{...}{Arguments to be passed to the method \code{"plot(, ...)"},
        such as graphical parameters.}
}
\details{
    The details will be described in future.
}
\value{
    None.
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{snpgdsCutTree}}
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


# split
set.seed(100)
rv <- snpgdsCutTree(hc, label.H=TRUE, label.Z=TRUE)

# draw dendrogram
snpgdsDrawTree(rv, main="HapMap Phase II",
    edgePar=list(col=rgb(0.5,0.5,0.5, 0.75), t.col="black"))
}

\keyword{GDS}
\keyword{GWAS}
