\name{snpgdsCutTree}
\alias{snpgdsCutTree}
\title{
	Determine sub groups of individuals
}
\description{
	To determine sub groups of individuals using a specified tree from hierarchical
cluster analysis
}
\usage{
snpgdsCutTree(hc, z.threshold=10, method=c("permutation", "exact"),
	n.perm = 1000, samp.group=NULL, verbose=TRUE)
}
\arguments{
	\item{hc}{an object of "snpgdsHCClass"}
	\item{z.threshold}{the threshold of Z score to determine whether split the node or not}
	\item{method}{"permutation", "exact"}
	\item{n.perm}{the times for permutation}
	\item{samp.group}{if \code{NULL}, determine groups by Z score; if a vector of factor,
		assign each individual in dendrogram with respect to \code{samp.group}}
	\item{verbose}{if TRUE, show information}
}
\details{
	The details will be described in future.
}
\value{
	Return a list (class "snpgdsIBSClass"):
	\item{sample.id}{the sample ids used in the analysis}
	\item{samp.group}{a vector of factor, indicating the group of each individual}
	\item{group.dist}{a matrix of pairwise group dissimilarity}
	\item{dendrogram}{the dendrogram of individuals of sub groups}
}
\author{Xiuwen Zheng \email{zhengx@u.washington.edu}}
\seealso{
	\code{\link{snpgdsHCluster}},
	\code{\link{snpgdsIBS}}, \code{\link{snpgdsDiss}}
}

\examples{
# open an example dataset (HapMap)
genofile <- openfn.gds(snpgdsExampleFileName())

pop.group <- as.factor(read.gdsn(index.gdsn(genofile, c("sample.annot", "pop.group"))))
pop.level <- levels(pop.group)

diss <- snpgdsDiss(genofile, num.thread=2)
hc <- snpgdsHCluster(diss)
rv <- snpgdsCutTree(hc)
plot(rv$dendrogram, leaflab="none", main="HapMap Phase II")

rv2 <- snpgdsCutTree(hc, samp.group=as.factor(pop.group))
plot(rv2$dendrogram, leaflab="none", main="HapMap Phase II")

# close the genotype file
closefn.gds(genofile)
}

\keyword{multicore}
\keyword{gds}
\keyword{GWAS}
\keyword{IBS}
