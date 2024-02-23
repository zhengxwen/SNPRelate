\name{snpgdsAdmixPlot}
\alias{snpgdsAdmixPlot}
\alias{snpgdsAdmixTable}
\title{
    Plot Ancestry Proportions
}
\description{
    Plot the admixture proportions according to their ancestries.
}
\usage{
snpgdsAdmixPlot(propmat, group=NULL, col=NULL, multiplot=TRUE, showgrp=TRUE,
    shownum=TRUE, ylim=TRUE, na.rm=TRUE)
snpgdsAdmixTable(propmat, group, sort=FALSE)
}
\arguments{
    \item{propmat}{a sample-by-ancestry matrix of proportion estimates,
        returned from \code{\link{snpgdsAdmixProp}()}}
    \item{group}{a character vector of a factor according to the rows
        in \code{propmat}}
    \item{col}{specify colors; if \code{group} is not specified, it is a color
        for each sample; otherwise specify colors for the groups}
    \item{multiplot}{single plot or multiple plots}
    \item{showgrp}{show group names in the plot; applicable when \code{group}
        is used}
    \item{shownum}{\code{TRUE}: show the number of each group on the X-axis
        in the figure; applicable when \code{group} is used}
    \item{ylim}{\code{TRUE}: y-axis is limited to [0, 1];
        \code{FALSE}: \code{ylim <- range(propmat)};
        a 2-length numeric vector: \code{ylim} used in \code{plot()}}
    \item{na.rm}{\code{TRUE}: remove the sample(s) according to the
        missing value(s) in \code{group}}
    \item{sort}{\code{TRUE}: rearranges the rows of proportion matrices
        into descending order}
}
\details{
    The minor allele frequency and missing rate for each SNP passed in
\code{snp.id} are calculated over all the samples in \code{sample.id}.
}
\value{
    \code{snpgdsAdmixPlot()}: none.

    \code{snpgdsAdmixTable()}: a list of \code{data.frame} consisting of
        \code{group, num, mean, sd, min, max}
}

\references{
    Zheng X, Weir BS.
    Eigenanalysis on SNP Data with an Interpretation of Identity by Descent.
    Theoretical Population Biology. 2015 Oct 23. pii: S0040-5809(15)00089-1.
    doi: 10.1016/j.tpb.2015.09.004.
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{snpgdsEIGMIX}}, \code{\link{snpgdsAdmixProp}}
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

# define groups
groups <- list(CEU = samp.id[pop_code == "CEU"],
    YRI = samp.id[pop_code == "YRI"],
    CHB = samp.id[is.element(pop_code, c("HCB", "JPT"))])

prop <- snpgdsAdmixProp(RV, groups=groups, bound=TRUE)

# draw
snpgdsAdmixPlot(prop, group=pop_code)

# use user-defined colors for the groups
snpgdsAdmixPlot(prop, group=pop_code, multiplot=FALSE, col=c(3,2,4))

snpgdsAdmixTable(prop, group=pop_code)

# close the genotype file
snpgdsClose(genofile)
}

\keyword{GDS}
\keyword{GWAS}
