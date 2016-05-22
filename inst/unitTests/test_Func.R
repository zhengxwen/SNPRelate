#############################################################
#
# DESCRIPTION: test IBS
#

library(RUnit)
library(SNPRelate)


#############################################################
# test function
#

test.AlleleFreq <- function()
{
	# open the SNP GDS file
	genofile <- snpgdsOpen(snpgdsExampleFileName())
	on.exit({ snpgdsClose(genofile) })

	# get genotype
	geno <- snpgdsGetGeno(genofile, snpfirstdim=FALSE, verbose=FALSE)

	af <- colMeans(geno, na.rm=TRUE) * 0.5
	maf <- pmin(af, 1 - af)
	mr <- colMeans(is.na(geno))
	x <- snpgdsSNPRateFreq(genofile)

    checkEquals(af, x$AlleleFreq, "allele frequency")
    checkEquals(mr, x$MissingRate, "missing rate")
    checkEquals(maf, x$MinorFreq, "minor allele frequency")
}


test.Allele_Switching <- function()
{
    # the file name of SNP GDS
    (fn <- snpgdsExampleFileName())

    # copy the file
    file.copy(fn, "test.gds", overwrite=TRUE)

    # open the SNP GDS file
    genofile <- snpgdsOpen("test.gds",readonly=FALSE, allow.duplicate=TRUE)

    # get genotype
    g1 <- snpgdsGetGeno(genofile, verbose=FALSE)

    # allelic information
    allele <- read.gdsn(index.gdsn(genofile, "snp.allele"))
    allele.list <- strsplit(allele, "/")

    A1 <- A.allele <- sapply(allele.list, function(x) { x[1] })
    A2 <- B.allele <- sapply(allele.list, function(x) { x[2] })

    set.seed(1000)
    flag <- rep(FALSE, length(A.allele))
    flag[sample.int(length(A.allele), 250, replace=TRUE)] <- TRUE

    A.allele[flag] <- B.allele[flag]
    A.allele[sample.int(length(A.allele), 10, replace=TRUE)] <- NA

    # allele switching
    flag <- snpgdsAlleleSwitch(genofile, A.allele, verbose=FALSE)

    # close the file
    snpgdsClose(genofile)

    # get genotype and alleles

    f <- snpgdsOpen("test.gds", allow.duplicate=TRUE)
    g2 <- snpgdsGetGeno(f, verbose=FALSE)
    alt <- read.gdsn(index.gdsn(f, "snp.allele"))
    snpgdsClose(f)


    ######    ######

    flag[is.na(flag)] <- FALSE

    g1[,flag] <- 2L - g1[,flag]
    checkEquals(g1, g2, "allele switching")

    alt1 <- paste(A1, A2, sep="/")
    alt1[flag] <- paste(A2, A1, sep="/")[flag]
    checkEquals(alt, alt1, "allele switching")


    # delete the temporary file
    unlink("test.gds", force=TRUE)
}
