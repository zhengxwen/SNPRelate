#############################################################
#
# DESCRIPTION: test PCA
#

library(RUnit)
library(SNPRelate)


CreatePCA <- function()
{
	# open a GDS file
	f <- snpgdsOpen(snpgdsExampleFileName())

	samp.id <- read.gdsn(index.gdsn(f, "sample.id"))
	pca <- snpgdsPCA(f, sample.id=samp.id[1:90], need.genmat=TRUE)
	pca$corr <- round(snpgdsPCACorr(pca, f, eig.which=1:2)$snpcorr, 3)

	save(pca, file="Validate.PCA.RData", compress="xz")

	snpgdsClose(f)
	invisible()
}



#############################################################
# test function
#

test.PCA <- function()
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

	corr.1 <- round(snpgdsPCACorr(pca.1, genofile, eig.which=1:2,
		num.thread=1)$snpcorr, 3)
	checkEquals(corr.1, valid.dta$corr, "PCA correlation (one core)")


	# run on one core
	pca.2 <- snpgdsPCA(genofile, sample.id=samp.id[1:90],
		num.thread=2, need.genmat=TRUE, verbose=FALSE)
	checkEquals(pca.2$genmat, valid.dta$genmat, "PCA (two cores)")

	corr.2 <- round(snpgdsPCACorr(pca.2, genofile, eig.which=1:2,
		num.thread=2)$snpcorr, 3)
	checkEquals(corr.2, valid.dta$corr, "PCA correlation (two cores)")


	# close the file
	snpgdsClose(genofile)
}
