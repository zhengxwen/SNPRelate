\name{snpgdsDiss}
\alias{snpgdsDiss}
\title{
    Individual dissimilarity analysis
}
\description{
    Calculate the individual dissimilarities for each pair of individuals.
}
\usage{
snpgdsDiss(gdsobj, sample.id=NULL, snp.id=NULL, autosome.only=TRUE,
    remove.monosnp=TRUE, maf=NaN, missing.rate=NaN, num.thread=1, verbose=TRUE)
}
\arguments{
    \item{gdsobj}{an object of class \code{\link{SNPGDSFileClass}},
        a SNP GDS file}
    \item{sample.id}{a vector of sample id specifying selected samples;
        if NULL, all samples are used}
    \item{snp.id}{a vector of snp id specifying selected SNPs;
        if NULL, all SNPs are used}
    \item{autosome.only}{if \code{TRUE}, use autosomal SNPs only; if it is a
        numeric or character value, keep SNPs according to the specified
        chromosome}
    \item{remove.monosnp}{if TRUE, remove monomorphic SNPs}
    \item{maf}{to use the SNPs with ">= maf" only; if NaN, no MAF threshold}
    \item{missing.rate}{to use the SNPs with "<= missing.rate" only;
        if NaN, no missing threshold}
    \item{num.thread}{the number of (CPU) cores used; if \code{NA}, detect
        the number of cores automatically}
    \item{verbose}{if TRUE, show information}
}
\details{
    The minor allele frequency and missing rate for each SNP passed in
\code{snp.id} are calculated over all the samples in \code{sample.id}.

    \code{snpgdsDiss()} returns \code{1 - beta_ij} which is formally described
in Weir&Goudet (2017).
}
\value{
    Return a class "snpgdsDissClass":
    \item{sample.id}{the sample ids used in the analysis}
    \item{snp.id}{the SNP ids used in the analysis}
    \item{diss}{a matrix of individual dissimilarity}
}

\references{
    Zheng, Xiuwen. 2013. Statistical Prediction of HLA Alleles and
Relatedness Analysis in Genome-Wide Association Studies. PhD dissertation,
the department of Biostatistics, University of Washington.

    Weir BS, Zheng X. SNPs and SNVs in Forensic Science. 2015.
Forensic Science International: Genetics Supplement Series.

    Weir BS, Goudet J. A Unified Characterization of Population Structure and
Relatedness. Genetics. 2017 Aug;206(4):2085-2103.
doi: 10.1534/genetics.116.198424.
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{snpgdsHCluster}}
}

\examples{
# open an example dataset (HapMap)
genofile <- snpgdsOpen(snpgdsExampleFileName())

pop.group <- as.factor(read.gdsn(index.gdsn(
    genofile, "sample.annot/pop.group")))
pop.level <- levels(pop.group)

diss <- snpgdsDiss(genofile)
diss
hc <- snpgdsHCluster(diss)
names(hc)
plot(hc$dendrogram)

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
