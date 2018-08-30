\name{snpgdsLDpruning}
\alias{snpgdsLDpruning}
\title{
    Linkage Disequilibrium (LD) based SNP pruning
}
\description{
    Recursively removes SNPs within a sliding window
}
\usage{
snpgdsLDpruning(gdsobj, sample.id=NULL, snp.id=NULL, autosome.only=TRUE,
    remove.monosnp=TRUE, maf=NaN, missing.rate=NaN,
    method=c("composite", "r", "dprime", "corr"), slide.max.bp=500000L,
    slide.max.n=NA, ld.threshold=0.2, start.pos=c("random", "first", "last"),
    num.thread=1L, verbose=TRUE)
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
    \item{method}{"composite", "r", "dprime", "corr", see details}
    \item{slide.max.bp}{the maximum basepairs in the sliding window}
    \item{slide.max.n}{the maximum number of SNPs in the sliding window}
    \item{ld.threshold}{the LD threshold}
    \item{start.pos}{"random": a random starting position; "first": start from
        the first position; "last": start from the last position}
    \item{num.thread}{the number of (CPU) cores used; if \code{NA}, detect
        the number of cores automatically}
    \item{verbose}{if TRUE, show information}
}
\details{
    The minor allele frequency and missing rate for each SNP passed in
\code{snp.id} are calculated over all the samples in \code{sample.id}.

    Four methods can be used to calculate linkage disequilibrium values:
"composite" for LD composite measure, "r" for R coefficient (by EM algorithm
assuming HWE, it could be negative), "dprime" for D', and "corr" for
correlation coefficient. The method "corr" is equivalent to "composite",
when SNP genotypes are coded as: 0 -- BB, 1 -- AB, 2 -- AA. The argument
\code{ld.threshold} is the absolute value of measurement.

    It is useful to generate a pruned subset of SNPs that are in approximate
linkage equilibrium with each other. The function \code{snpgdsLDpruning}
recursively removes SNPs within a sliding window based on the pairwise
genotypic correlation. SNP pruning is conducted chromosome by chromosome,
since SNPs in a chromosome can be considered to be independent with the other
chromosomes.

    The pruning algorithm on a chromosome is described as follows (n is
the total number of SNPs on that chromosome):

    1) Randomly select a starting position i (\code{start.pos="random"}), 
i=1 if \code{start.pos="first"}, or i=last if \code{start.pos="last"}; and let
the current SNP set S=\{ i \};

    2) For each right position j from i+1 to n: if any LD between j and k is
greater than \code{ld.threshold}, where k belongs to S, and both of j and k
are in the sliding window, then skip j; otherwise, let S be S + \{ j \};

    3) For each left position j from i-1 to 1: if any LD between j and k is
greater than \code{ld.threshold}, where k belongs to S, and both of j and k
are in the sliding window, then skip j; otherwise, let S be S + \{ j \};

    4) Output S, the final selection of SNPs.
}
\value{
    Return a list of SNP IDs stratified by chromosomes.
}
\references{
    Weir B: Inferences about linkage disequilibrium.
    Biometrics 1979; 35: 235-254.

    Weir B: Genetic Data Analysis II. Sunderland, MA: Sinauer Associates, 1996.

    Weir BS, Cockerham CC:
    Complete characterization of disequilibrium at two loci;
    in Feldman MW (ed): Mathematical Evolutionary Theory.
    Princeton, NJ: Princeton University Press, 1989.
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{snpgdsLDMat}}, \code{\link{snpgdsLDpair}}
}

\examples{
# open an example dataset (HapMap)
genofile <- snpgdsOpen(snpgdsExampleFileName())

set.seed(1000)
snpset <- snpgdsLDpruning(genofile)
names(snpset)
#  [1] "chr1"  "chr2"  "chr3"  "chr4"  "chr5"  "chr6"  "chr7"  "chr8"  "chr9"
# [10] "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18"
# ......
head(snpset$chr1)
# [1] 1 2 3 4 5 6

# get SNP ids
snp.id <- unlist(unname(snpset))

# close the genotype file
snpgdsClose(genofile)
}

\keyword{GDS}
\keyword{GWAS}
