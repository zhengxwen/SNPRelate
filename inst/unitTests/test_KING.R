#############################################################
#
# DESCRIPTION: test the KING IBD/IBS method
#

library(RUnit)
library(SNPRelate)


# open a GDS file
# genofile <- snpgdsOpen(snpgdsExampleFileName())
# samp.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
# v1 <- snpgdsIBDKING(genofile, sample.id=samp.id[1:60], type="KING-robust")
# v2 <- snpgdsIBDKING(genofile, sample.id=samp.id[1:60], type="KING-homo")
# .king <- list(v1, v2)
# save(.king, file="Validate.KING.RData", compress="xz")
# snpgdsClose(genofile)




#############################################################
# test function
#

test.KING <- function()
{
    valid.dta <- get(load(system.file(
        "unitTests", "valid", "Validate.KING.RData", package="SNPRelate")))

    # open a GDS file
    genofile <- snpgdsOpen(snpgdsExampleFileName(), allow.duplicate=TRUE)
    samp.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

    # run on one core
    king.1 <- snpgdsIBDKING(genofile, sample.id=samp.id[1:60],
        type="KING-robust", num.thread=1L, verbose=FALSE)
    checkEquals(king.1, valid.dta[[1L]], "KING robust MoM (one core)")

    king.2 <- snpgdsIBDKING(genofile, sample.id=samp.id[1:60],
        type="KING-homo", num.thread=1L, verbose=FALSE)
    checkEquals(king.2, valid.dta[[2L]], "KING homo MoM (one core)")

    # run on two cores
    king.1 <- snpgdsIBDKING(genofile, sample.id=samp.id[1:60],
        type="KING-robust", num.thread=2L, verbose=FALSE)
    checkEquals(king.1, valid.dta[[1L]], "KING robust MoM (two cores)")

    king.2 <- snpgdsIBDKING(genofile, sample.id=samp.id[1:60],
        type="KING-homo", num.thread=2L, verbose=FALSE)
    checkEquals(king.2, valid.dta[[2L]], "KING homo MoM (two cores)")

    # close the file
    snpgdsClose(genofile)
}
