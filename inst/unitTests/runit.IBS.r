#############################################################
#
# DESCRIPTION: test IBS
#

library(RUnit)
library(SNPRelate)


# open a GDS file
# genofile <- openfn.gds(snpgdsExampleFileName())

# ibs <- snpgdsIBS(genofile)
# save(ibs, file="Validate.IBS.RData")

# closefn.gds(genofile)






#############################################################
# test function
#

test.IBS <- function()
{
	valid.dta <- get(load(
		system.file("unitTests", "valid", "Validate.IBS.RData", package="SNPRelate")))

	# open a GDS file
	genofile <- openfn.gds(snpgdsExampleFileName())

	# run on one core
	ibs.1 <- snpgdsIBS(genofile, num.thread=1)
	checkEquals(ibs.1, valid.dta, "IBS (one core)")

	# run on four core
	ibs.4 <- snpgdsIBS(genofile, num.thread=4)
	checkEquals(ibs.4, valid.dta, "IBS (four cores)")

	# run on 16 core
	ibs.16 <- snpgdsIBS(genofile, num.thread=16)
	checkEquals(ibs.16, valid.dta, "IBS (16 cores)")

	# close the file
	closefn.gds(genofile)
}
