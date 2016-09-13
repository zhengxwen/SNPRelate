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
	pca <- snpgdsPCA(f, sample.id=samp.id[1:90], need.genmat=TRUE, eigen.cnt=8L)
	corr <- round(snpgdsPCACorr(pca, f, eig.which=1:2)$snpcorr, 3)

	SnpLoad <- snpgdsPCASNPLoading(pca, f)
	snploading <- round(SnpLoad$snploading, 3)

	SL <- snpgdsPCASampLoading(SnpLoad, f, sample.id=samp.id[1:100])

	.rv <- list(genmat = pca$genmat, corr = corr,
		snploading = snploading,
		samploading = round(SL$eigenvect, 4))
	save(.rv, file="Validate.PCA.RData", compress="xz")

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
	on.exit({ snpgdsClose(genofile) })

	samp.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

	# run on one core
	pca <- snpgdsPCA(genofile, sample.id=samp.id[1:90],
		num.thread=1, need.genmat=TRUE, eigen.cnt=8L, verbose=FALSE)
	checkEquals(pca$genmat, valid.dta$genmat, "PCA (one core)")

	corr <- round(snpgdsPCACorr(pca, genofile, eig.which=1:2,
		num.thread=1)$snpcorr, 3)
	checkEquals(corr, valid.dta$corr, "PCA correlation (one core)")

	SnpLoad <- snpgdsPCASNPLoading(pca, genofile)
	snploading <- round(SnpLoad$snploading, 3)
	checkEquals(snploading, valid.dta$snploading, "PCA SNP loading (one core)")

	SL <- snpgdsPCASampLoading(SnpLoad, genofile, sample.id=samp.id[1:100])
	checkEquals(round(SL$eigenvect, 4), valid.dta$samploading,
		"PCA sample loading (one core)")


	# run on one core
	pca <- snpgdsPCA(genofile, sample.id=samp.id[1:90],
		num.thread=2, need.genmat=TRUE, eigen.cnt=8L, verbose=FALSE)
	checkEquals(pca$genmat, valid.dta$genmat, "PCA (two cores)")

	corr <- round(snpgdsPCACorr(pca, genofile, eig.which=1:2,
		num.thread=2)$snpcorr, 3)
	checkEquals(corr, valid.dta$corr, "PCA correlation (two cores)")

	SnpLoad <- snpgdsPCASNPLoading(pca, genofile)
	snploading <- round(SnpLoad$snploading, 3)
	checkEquals(snploading, valid.dta$snploading, "PCA SNP loading (two cores)")

	SL <- snpgdsPCASampLoading(SnpLoad, genofile, sample.id=samp.id[1:100])
	checkEquals(round(SL$eigenvect, 4), valid.dta$samploading,
		"PCA sample loading (two cores)")
}
