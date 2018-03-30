#############################################################
#
# DESCRIPTION: test the calculation GRM matrix
#

library(RUnit)
library(SNPRelate)


#############################################################
# test function
#

test.merge.GCTA.grm <- function()
{
	# open an example dataset (HapMap)
	f <- snpgdsOpen(snpgdsExampleFileName())

	# there is no missing genotype
	snpid <- read.gdsn(index.gdsn(f, "snp.id"))
	snpid <- snpid[snpgdsSNPRateFreq(f)$MissingRate == 0]
	# split the SNP set
	snp1 <- snpid[1:1000]
	snp2 <- snpid[1001:3000]
	snp3 <- setdiff(snpid, c(snp1, snp2))

	# run
	snpgdsGRM(f, snp.id=snp1, method="GCTA", out.fn="tmp1.gds")
	snpgdsGRM(f, snp.id=snp2, method="GCTA", out.fn="tmp2.gds")
	snpgdsGRM(f, snp.id=snp3, method="GCTA", out.fn="tmp3.gds")
	# merge GRMs and export to a new GDS file
	snpgdsMergeGRM(c("tmp1.gds", "tmp2.gds", "tmp3.gds"), "tmp.gds")

	# run using all SNPs
	grm <- snpgdsGRM(f, method="GCTA", snp.id=snpid)
	# close the file
	snpgdsClose(f)

	# check
	f <- openfn.gds("tmp.gds")
	m <- read.gdsn(index.gdsn(f, "grm"))
	closefn.gds(f)

    # check
    checkEquals(m, grm$grm, "check the merged GCTA GRM")

	# delete the temporary file
	unlink(c("tmp1.gds", "tmp2.gds", "tmp3.gds", "tmp.gds"), force=TRUE)
}


test.merge.beta.grm <- function()
{
	# open an example dataset (HapMap)
	f <- snpgdsOpen(snpgdsExampleFileName())

	# there is no missing genotype
	snpid <- read.gdsn(index.gdsn(f, "snp.id"))
	snpid <- snpid[snpgdsSNPRateFreq(f)$MissingRate == 0]
	# split the SNP set
	snp1 <- snpid[1:1000]
	snp2 <- snpid[1001:3000]
	snp3 <- setdiff(snpid, c(snp1, snp2))

	# run
	snpgdsGRM(f, snp.id=snp1, method="IndivBeta", out.fn="tmp1.gds")
	snpgdsGRM(f, snp.id=snp2, method="IndivBeta", out.fn="tmp2.gds")
	snpgdsGRM(f, snp.id=snp3, method="IndivBeta", out.fn="tmp3.gds")
	# merge GRMs and export to a new GDS file
	snpgdsMergeGRM(c("tmp1.gds", "tmp2.gds", "tmp3.gds"), "tmp.gds")

	# run using all SNPs
	grm <- snpgdsGRM(f, method="IndivBeta", snp.id=snpid)
	# close the file
	snpgdsClose(f)

	# check
	f <- openfn.gds("tmp.gds")
	m <- read.gdsn(index.gdsn(f, "grm"))
	closefn.gds(f)

    # check
    checkEquals(m, grm$grm, "check the merged beta-based GRM")

	# delete the temporary file
	unlink(c("tmp1.gds", "tmp2.gds", "tmp3.gds", "tmp.gds"), force=TRUE)
}
