#############################################################
#
# DESCRIPTION: test PCA
#

library(RUnit)
library(SNPRelate)


# open a GDS file
# genofile <- openfn.gds(snpgdsExampleFileName())

# pca <- snpgdsPCA(genofile, need.genmat=TRUE)
# save(pca, file="Validate.PCA.RData")

# closefn.gds(genofile)






#############################################################
# test function
#

test.PCA <- function()
{
	valid.dta <- get(load(
		system.file("unitTests", "valid", "Validate.PCA.RData", package="SNPRelate")))

	# open a GDS file
	genofile <- openfn.gds(snpgdsExampleFileName())

	# run on one core
	pca.1 <- snpgdsPCA(genofile, num.thread=1, need.genmat=TRUE)
	checkEquals(pca.1, valid.dta, "PCA (one core)")

	# run on four core
	pca.4 <- snpgdsPCA(genofile, num.thread=4, need.genmat=TRUE)
	checkEquals(pca.4, valid.dta, "PCA (four cores)")

	# run on 16 core
	pca.16 <- snpgdsPCA(genofile, num.thread=16, need.genmat=TRUE)
	checkEquals(pca.16, valid.dta, "PCA (16 cores)")

	# close the file
	closefn.gds(genofile)
}
