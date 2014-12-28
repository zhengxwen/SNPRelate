#############################################################
#
# DESCRIPTION: test IBS
#

library(RUnit)
library(SNPRelate)


# open a GDS file
# genofile <- snpgdsOpen(snpgdsExampleFileName())
# samp.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
# ibs <- snpgdsIBS(genofile, sample.id=samp.id[1:90])
# save(ibs, file="Validate.IBS.RData", compress="xz")
# snpgdsClose(genofile)




#############################################################
# test function
#

.test.IBS <- function()
{
    valid.dta <- get(load(system.file(
        "unitTests", "valid", "Validate.IBS.RData", package="SNPRelate")))

    # open a GDS file
    genofile <- snpgdsOpen(snpgdsExampleFileName(), allow.duplicate=TRUE)
    samp.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

    # run on one core
    ibs.1 <- snpgdsIBS(genofile, sample.id=samp.id[1:90],
        num.thread=1, verbose=FALSE)
    checkEquals(ibs.1, valid.dta, "IBS (one core)")

    # run on 4 cores
    ibs.4 <- snpgdsIBS(genofile, sample.id=samp.id[1:90],
        num.thread=4, verbose=FALSE)
    checkEquals(ibs.4, valid.dta, "IBS (four cores)")

    # run on 16 cores
    ibs.16 <- snpgdsIBS(genofile, sample.id=samp.id[1:90],
        num.thread=16, verbose=FALSE)
    checkEquals(ibs.16, valid.dta, "IBS (16 cores)")

    # close the file
    snpgdsClose(genofile)
}
