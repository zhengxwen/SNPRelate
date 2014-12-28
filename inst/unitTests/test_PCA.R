#############################################################
#
# DESCRIPTION: test PCA
#

library(RUnit)
library(SNPRelate)


# open a GDS file
# genofile <- snpgdsOpen(snpgdsExampleFileName())
# samp.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
# pca <- snpgdsPCA(genofile, sample.id=samp.id[1:90], need.genmat=TRUE)
# save(pca, file="Validate.PCA.RData", compress="xz")
# snpgdsClose(genofile)




#############################################################
# test function
#

.test.PCA <- function()
{
    valid.dta <- get(load(system.file(
        "unitTests", "valid", "Validate.PCA.RData", package="SNPRelate")))

    # open a GDS file
    genofile <- snpgdsOpen(snpgdsExampleFileName(), allow.duplicate=TRUE)

    samp.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

    # run on one core
    pca.1 <- snpgdsPCA(genofile, sample.id=samp.id[1:90],
        num.thread=1, need.genmat=TRUE, verbose=FALSE)
    checkEquals(pca.1$genmat, valid.dta$genmat, "PCA (one core)")

    # run on 4 cores
    pca.4 <- snpgdsPCA(genofile, sample.id=samp.id[1:90],
        num.thread=4, need.genmat=TRUE, verbose=FALSE)
    checkEquals(pca.4$genmat, valid.dta$genmat, "PCA (four cores)")

    # run on 16 cores
    pca.16 <- snpgdsPCA(genofile, sample.id=samp.id[1:90],
        num.thread=16, need.genmat=TRUE, verbose=FALSE)
    checkEquals(pca.16$genmat, valid.dta$genmat, "PCA (16 cores)")

    # close the file
    snpgdsClose(genofile)
}
