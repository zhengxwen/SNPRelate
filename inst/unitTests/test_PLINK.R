#############################################################
#
# DESCRIPTION: test PLINK Method of Moment
#

library(RUnit)
library(SNPRelate)


# open a GDS file
# genofile <- snpgdsOpen(snpgdsExampleFileName())
# samp.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
# ibd <- snpgdsIBDMoM(genofile, sample.id=samp.id[1:90])
# save(ibd, file="Validate.MoM.RData", compress="xz")
# snpgdsClose(genofile)




#############################################################
# test function
#

test.PLINK.MoM <- function()
{
    valid.dta <- get(load(system.file(
        "unitTests", "valid", "Validate.MoM.RData", package="SNPRelate")))

    # open a GDS file
    genofile <- snpgdsOpen(snpgdsExampleFileName(), allow.duplicate=TRUE)
    samp.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

    # run on one core
    ibd.1 <- snpgdsIBDMoM(genofile, sample.id=samp.id[1:90],
        num.thread=1, verbose=FALSE)
    checkEquals(ibd.1, valid.dta, "PLINK MoM (one core)")

    # run on one core
    ibd.2 <- snpgdsIBDMoM(genofile, sample.id=samp.id[1:90],
        num.thread=2, verbose=FALSE)
    checkEquals(ibd.2, valid.dta, "PLINK MoM (two cores)")

    # close the file
    snpgdsClose(genofile)
}
