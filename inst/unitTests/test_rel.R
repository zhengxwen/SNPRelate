#############################################################
#
# DESCRIPTION: test units for relatedness
#

library(RUnit)
library(SNPRelate)


CreateIBS <- function()
{
	genofile <- snpgdsOpen(snpgdsExampleFileName())
	samp.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
	ibs <- snpgdsIBS(genofile, sample.id=samp.id[1:90])
	save(ibs, file="Validate.IBS.RData", compress="xz")
	snpgdsClose(genofile)
}


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


CreatePLINK <- function()
{
	genofile <- snpgdsOpen(snpgdsExampleFileName())
	samp.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
	ibd <- snpgdsIBDMoM(genofile, sample.id=samp.id[1:90])
	save(ibd, file="Validate.MoM.RData", compress="xz")
	snpgdsClose(genofile)
}


CreateKING <- function()
{
	genofile <- snpgdsOpen(snpgdsExampleFileName())
	samp.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
	v1 <- snpgdsIBDKING(genofile, sample.id=samp.id[1:60], type="KING-robust")
	v2 <- snpgdsIBDKING(genofile, sample.id=samp.id[1:60], type="KING-homo")
	.king <- list(v1, v2)
	save(.king, file="Validate.KING.RData", compress="xz")
	snpgdsClose(genofile)
}


CreateIndivBeta <- function()
{
	genofile <- snpgdsOpen(snpgdsExampleFileName())
	samp.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
	.beta <- snpgdsIndivBeta(genofile, sample.id=samp.id[1:90])
	save(.beta, file="Validate.Beta.RData", compress="xz")
	snpgdsClose(genofile)
}


CreateEIGMIX <- function()
{
	genofile <- snpgdsOpen(snpgdsExampleFileName())
	samp.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
	.eigmix <- snpgdsEIGMIX(genofile, sample.id=samp.id[1:90], eigen.cnt=0,
		ibdmat=TRUE)$ibd
	save(.eigmix, file="Validate.EIGMIX.RData", compress="xz")
	snpgdsClose(genofile)
}




##############################################################################
# test functions
#

test.IBS <- function()
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

	ibs.1 <- snpgdsIBS(genofile, sample.id=samp.id[1:90],
		num.thread=1, useMatrix=TRUE, verbose=FALSE)
	z <- as.matrix(ibs.1$ibs)
	dimnames(z) <- NULL
	checkEquals(z, valid.dta$ibs, "IBS (one core, Matrix)")

	# run on two cores
	ibs.2 <- snpgdsIBS(genofile, sample.id=samp.id[1:90],
		num.thread=2, verbose=FALSE)
	checkEquals(ibs.2, valid.dta, "IBS (two cores)")

	# close the file
	snpgdsClose(genofile)
}



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

	corrX <- snpgdsPCACorr(pca, genofile, eig.which=1:2, num.thread=1)$snpcorr
	corr <- round(corrX, 3)
	checkEquals(corr, valid.dta$corr, "PCA correlation (one core)")

	snpgdsPCACorr(pca, genofile, eig.which=1:2, num.thread=1, outgds="test.gds")
	f <- openfn.gds("test.gds")
	corr <- read.gdsn(index.gdsn(f, "correlation"))
	closefn.gds(f)
	checkEquals(corr, round(corrX, 4), "PCA correlation (one core)")

	SnpLoad <- snpgdsPCASNPLoading(pca, genofile)
	snploading <- round(SnpLoad$snploading, 3)
	checkEquals(snploading, valid.dta$snploading, "PCA SNP loading (one core)")

	SL <- snpgdsPCASampLoading(SnpLoad, genofile, sample.id=samp.id[1:100])
	checkEquals(round(SL$eigenvect, 4), valid.dta$samploading,
		"PCA sample loading (one core)")


	# run on two cores
	pca <- snpgdsPCA(genofile, sample.id=samp.id[1:90],
		num.thread=2, need.genmat=TRUE, eigen.cnt=8L, verbose=FALSE)
	checkEquals(pca$genmat, valid.dta$genmat, "PCA (two cores)")

	corrX <- snpgdsPCACorr(pca, genofile, eig.which=1:2, num.thread=2)$snpcorr
	corr <- round(corrX, 3)
	checkEquals(corr, valid.dta$corr, "PCA correlation (two cores)")

	snpgdsPCACorr(pca, genofile, eig.which=1:2, num.thread=2, outgds="test.gds")
	f <- openfn.gds("test.gds")
	corr <- read.gdsn(index.gdsn(f, "correlation"))
	closefn.gds(f)
	checkEquals(corr, round(corrX, 4), "PCA correlation (one core)")

	SnpLoad <- snpgdsPCASNPLoading(pca, genofile)
	snploading <- round(SnpLoad$snploading, 3)
	checkEquals(snploading, valid.dta$snploading, "PCA SNP loading (two cores)")

	SL <- snpgdsPCASampLoading(SnpLoad, genofile, sample.id=samp.id[1:100])
	checkEquals(round(SL$eigenvect, 4), valid.dta$samploading,
		"PCA sample loading (two cores)")


	# delete the temporary file
	unlink("test.gds", force=TRUE)
}



test.PLINK.MoM <- function()
{
	valid.dta <- get(load(system.file(
		"unitTests", "valid", "Validate.MoM.RData", package="SNPRelate")))

	# open a GDS file
	genofile <- snpgdsOpen(snpgdsExampleFileName(), allow.duplicate=TRUE)
	samp.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

	# run on one core
	ibd.1 <- snpgdsIBDMoM(genofile, sample.id=samp.id[1:90],
		num.thread=1L, verbose=FALSE)
	checkEquals(ibd.1, valid.dta, "PLINK MoM (one core)")

	ibd <- snpgdsIBDMoM(genofile, sample.id=samp.id[1:90],
		num.thread=1L, useMatrix=TRUE, verbose=FALSE)
	checkEquals(as.numeric(ibd.1$k0), as.numeric(ibd$k0), "PLINK MoM (one core)")
	checkEquals(as.numeric(ibd.1$k1), as.numeric(ibd$k1), "PLINK MoM (one core)")

	# run on two cores
	ibd.2 <- snpgdsIBDMoM(genofile, sample.id=samp.id[1:90],
		num.thread=2L, verbose=FALSE)
	checkEquals(ibd.2, valid.dta, "PLINK MoM (two cores)")

	ibd <- snpgdsIBDMoM(genofile, sample.id=samp.id[1:90],
		num.thread=1L, useMatrix=TRUE, verbose=FALSE)
	checkEquals(as.numeric(ibd.2$k0), as.numeric(ibd$k0), "PLINK MoM (two cores)")
	checkEquals(as.numeric(ibd.2$k1), as.numeric(ibd$k1), "PLINK MoM (two cores)")

	# close the file
	snpgdsClose(genofile)
}



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

	king.1 <- snpgdsIBDKING(genofile, sample.id=samp.id[1:60],
		type="KING-robust", num.thread=1L, useMatrix=TRUE, verbose=FALSE)
	king.1$IBS0 <- as.matrix(king.1$IBS0)
	dimnames(king.1$IBS0) <- NULL
	king.1$kinship <- as.matrix(king.1$kinship)
	dimnames(king.1$kinship) <- NULL
	checkEquals(king.1, valid.dta[[1L]], "KING robust MoM (one core, Matrix)")

	king.2 <- snpgdsIBDKING(genofile, sample.id=samp.id[1:60],
		type="KING-homo", num.thread=1L, verbose=FALSE)
	checkEquals(king.2, valid.dta[[2L]], "KING homo MoM (one core)")

	king.2 <- snpgdsIBDKING(genofile, sample.id=samp.id[1:60],
		type="KING-homo", num.thread=1L, verbose=FALSE)
	king.2$k0 <- as.matrix(king.2$k0)
	dimnames(king.2$k0) <- NULL
	king.2$k1 <- as.matrix(king.2$k1)
	dimnames(king.2$k1) <- NULL
	checkEquals(king.2, valid.dta[[2L]], "KING homo MoM (one core, Matrix)")

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



test.IndivBeta <- function()
{
	valid.dta <- get(load(system.file(
		"unitTests", "valid", "Validate.Beta.RData", package="SNPRelate")))

	# open a GDS file
	genofile <- snpgdsOpen(snpgdsExampleFileName(), allow.duplicate=TRUE)
	samp.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

	# run on one core
	beta.1 <- snpgdsIndivBeta(genofile, sample.id=samp.id[1:90],
		num.thread=1, verbose=FALSE)
	checkEquals(beta.1$beta, valid.dta$beta, "Individual Beta (one core)")

	beta.1 <- snpgdsIndivBeta(genofile, sample.id=samp.id[1:90],
		num.thread=1, useMatrix=TRUE, verbose=FALSE)
	beta.1$beta <- as.matrix(beta.1$beta)
	dimnames(beta.1$beta) <- NULL
	checkEquals(beta.1$beta, valid.dta$beta, "Individual Beta (one core, Matrix)")

	# run on two cores
	beta.2 <- snpgdsIndivBeta(genofile, sample.id=samp.id[1:90],
		num.thread=2, verbose=FALSE)
	checkEquals(beta.2$beta, valid.dta$beta, "Individual Beta (two cores)")

	# close the file
	snpgdsClose(genofile)
}



test.EIGMIX <- function()
{
	valid.dta <- get(load(system.file(
		"unitTests", "valid", "Validate.EIGMIX.RData", package="SNPRelate")))

	# open a GDS file
	genofile <- snpgdsOpen(snpgdsExampleFileName(), allow.duplicate=TRUE)
	samp.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

	# run on one core
	eigmix.1 <- snpgdsEIGMIX(genofile, sample.id=samp.id[1:90], ibdmat=TRUE,
		num.thread=1, verbose=FALSE)
	checkEquals(eigmix.1$ibd, valid.dta, "EIGMIX (one core)")

	# run on two cores
	eigmix.2 <- snpgdsEIGMIX(genofile, sample.id=samp.id[1:90], ibdmat=TRUE,
		num.thread=2, verbose=FALSE)
	checkEquals(eigmix.2$ibd, valid.dta, "EIGMIX (two cores)")

	# close the file
	snpgdsClose(genofile)
}
