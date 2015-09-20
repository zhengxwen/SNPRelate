#############################################################
#
# DESCRIPTION: LD matrix
#

library(RUnit)
library(SNPRelate)


#############################################################
# test function
#

test.LD_Matrix <- function()
{
    # the file name of SNP GDS
    fn <- snpgdsExampleFileName()

    f <- snpgdsOpen(fn)
    snpset <- read.gdsn(index.gdsn(f, "snp.id"), start=1, count=1000)

    # SNP matrix
    geno <- snpgdsGetGeno(f, snp.id=snpset, verbose=FALSE)

    c1 <- snpgdsLDMat(f, snp.id=snpset, method="cov", slide=-1, with.id=FALSE)
    c2 <- cov(geno, use="pairwise.complete.obs")
    checkEquals(c1, c2, "covariance matrix")

    c1 <- snpgdsLDMat(f, snp.id=snpset, method="corr", slide=-1, with.id=FALSE)
    c2 <- suppressWarnings(cor(geno, use="pairwise.complete.obs"))
    checkEquals(c1, c2, "covariance matrix")

    # close the GDS file
    snpgdsClose(f)
}
