#############################################################
#
# DESCRIPTION: test PLINK Method of Moment
#

library(RUnit)
library(SNPRelate)


# open a GDS file
# genofile <- openfn.gds(snpgdsExampleFileName())

# ibd <- snpgdsIBDMoM(genofile)
# save(ibd, file="Validate.MoM.RData")

# closefn.gds(genofile)






#############################################################
# test function
#

test.PLINK.MoM <- function()
{
	valid.dta <- get(load(
		system.file("unitTests", "valid", "Validate.MoM.RData", package="SNPRelate")))

	# open a GDS file
	genofile <- openfn.gds(snpgdsExampleFileName())

	# run on one core
	ibd.1 <- snpgdsIBDMoM(genofile, num.thread=1)
	checkEquals(ibd.1, valid.dta, "PLINK MoM (one core)")

	# run on four core
	ibd.4 <- snpgdsIBDMoM(genofile, num.thread=4)
	checkEquals(ibd.4, valid.dta, "PLINK MoM (four cores)")

	# run on 16 core
	ibd.16 <- snpgdsIBDMoM(genofile, num.thread=16)
	checkEquals(ibd.16, valid.dta, "PLINK MoM (16 cores)")

	# close the file
	closefn.gds(genofile)
}
