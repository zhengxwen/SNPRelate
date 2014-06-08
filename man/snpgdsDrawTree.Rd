\name{snpgdsDrawTree}
\alias{snpgdsDrawTree}
\title{
	Draw a dendrogram
}
\description{
	To draw a dendrogram or the distribution of Z scores
}
\usage{
snpgdsDrawTree(obj, clust.count=NULL, type=c("dendrogram", "z-score"),
	outlier.n=NULL, shadow.col = c(rgb(0.5, 0.5, 0.5, 0.25), rgb(0.5, 0.5, 0.5, 0.05)),
	outlier.col = rgb(1, 0.50, 0.50, 0.5), leaflab="none", labels=NULL, y.label=0.2,
	...)
}
\arguments{
	\item{obj}{an object returned by \code{\link{snpgdsCutTree}}}
	\item{clust.count}{the counts for clusters, drawing shadows}
	\item{type}{"dendrogram", draw a dendrogram; or "z-score", draw the distribution of Z score}
	\item{outlier.n}{the cluster with size less than or equal to \code{outlier.n} is
		considered as outliers; if \code{NULL}, let \code{outlier.n = obj$outlier.n}}
	\item{shadow.col}{two colors for shadow}
	\item{outlier.col}{the colors for outliers}
	\item{leaflab}{a string specifying how leaves are labeled. The default
		\code{"perpendicular"} write text vertically (by default). \code{"textlike"}
		writes text horizontally (in a rectangle), and \code{"none"} suppresses leaf labels.}
	\item{labels}{the legend for different regions}
	\item{y.label}{y positions of labels}
	\item{...}{Arguments to be passed to the method \code{"plot(...)"}, such as graphical
		parameters.}
}
\details{
	The details will be described in future.
}
\value{
	None.
}
\author{Xiuwen Zheng \email{zhengx@u.washington.edu}}
\seealso{
	\code{\link{snpgdsCutTree}}
}

\keyword{multicore}
\keyword{gds}
\keyword{GWAS}
\keyword{IBS}
