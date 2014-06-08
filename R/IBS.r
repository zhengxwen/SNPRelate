#######################################################################
#
# Package name: SNPRelate
#
# Description:
#     A High-performance Computing Toolset for Relatedness and
# Principal Component Analysis of SNP Data
#
# Author: Xiuwen Zheng
# Email: zhengx@u.washington.edu
#


#######################################################################
# Identity-By-State (IBS) analysis
#######################################################################

#######################################################################
# To calculate the identity-by-state (IBS) matrix for SNP genotypes
#
# INPUT:
#   gdsobj -- an object of gds file
#   sample.id -- a vector of sample.id; if NULL, to use all samples
#   snp.id -- a vector of snp.id; if NULL, to use all SNPs
#   autosome.only -- whether only use autosomal SNPs
#   remove.monosnp -- whether remove monomorphic snps or not
#   maf -- the threshold of minor allele frequencies, keeping ">= maf"
#   missing.rate -- the threshold of missing rates, keeping "<= missing.rate"
#   verbose -- show information
#

snpgdsIBS <- function(gdsobj, sample.id=NULL, snp.id=NULL,
	autosome.only=TRUE, remove.monosnp=TRUE, maf=NaN, missing.rate=NaN,
	num.thread=1, verbose=TRUE)
{
	# check
	stopifnot(inherits(gdsobj, "gds.class"))
	stopifnot(is.numeric(num.thread) & (num.thread>0))
	stopifnot(is.logical(verbose))

	# samples
	sample.ids <- read.gdsn(index.gdsn(gdsobj, "sample.id"))
	if (!is.null(sample.id))
	{
		n.tmp <- length(sample.id)
		sample.id <- sample.ids %in% sample.id
		n.samp <- sum(sample.id);
		if (n.samp != n.tmp)
			stop("Some of sample.id do not exist!")
		if (n.samp <= 0)
			stop("No sample in the working dataset.")
		sample.ids <- sample.ids[sample.id]
	}

	if (verbose)
		cat("Identity-By-State (IBS) analysis on SNP genotypes:\n");

	# SNPs
	snp.ids <- read.gdsn(index.gdsn(gdsobj, "snp.id"))
	if (!is.null(snp.id))
	{
		n.tmp <- length(snp.id)
		snp.id <- snp.ids %in% snp.id
		n.snp <- sum(snp.id)
		if (n.snp != n.tmp)
			stop("Some of snp.id do not exist!")
		if (n.snp <= 0)
			stop("No SNP in the working dataset.")
		if (autosome.only)
		{
			chr <- read.gdsn(index.gdsn(gdsobj, "snp.chromosome"))
			opt <- snpgdsOption(gdsobj)
			snp.id <- snp.id & (chr %in% c(opt$autosome.start : opt$autosome.end))
			if (verbose)
			{
				tmp <- n.snp - sum(snp.id)
				if (tmp > 0) cat("Removing", tmp, "non-autosomal SNPs\n")
			}
		}
		snp.ids <- snp.ids[snp.id]
	} else {
		if (autosome.only)
		{
			chr <- read.gdsn(index.gdsn(gdsobj, "snp.chromosome"))
			opt <- snpgdsOption(gdsobj)
			snp.id <- chr %in% c(opt$autosome.start : opt$autosome.end)
			snp.ids <- snp.ids[snp.id]
			if (verbose)
				cat("Removing", length(chr) - length(snp.ids), "non-autosomal SNPs\n")
		}
	}

	# call C codes
	# set genotype working space
	node <- .C("gnrSetGenoSpace", as.integer(index.gdsn(gdsobj, "genotype")),
		as.logical(sample.id), as.logical(!is.null(sample.id)),
		as.logical(snp.id), as.logical(!is.null(snp.id)),
		n.snp=integer(1), n.samp=integer(1),
		err=integer(1), NAOK=TRUE, PACKAGE="SNPRelate")
	if (node$err != 0) stop(snpgdsErrMsg())

	# call allele freq. and missing rates
	if (remove.monosnp || is.finite(maf) || is.finite(missing.rate))
	{
		if (!is.finite(maf)) maf <- -1;
		if (!is.finite(missing.rate)) missing.rate <- 2;
		# call
		rv <- .C("gnrSelSNP_Base", as.logical(remove.monosnp),
			as.double(maf), as.double(missing.rate),
			out.num=integer(1), out.snpflag = logical(node$n.snp),
			err=integer(1), NAOK=TRUE, PACKAGE="SNPRelate")
		if (rv$err != 0) stop(snpgdsErrMsg())
		snp.ids <- snp.ids[rv$out.snpflag]
		# show
		if (verbose)
			cat("Removing", rv$out.num, "SNPs (monomorphic, < MAF, or > missing rate)\n")
	}

	# get the dimension of SNP genotypes
	node <- .C("gnrGetGenoDim", n.snp=integer(1), n.samp=integer(1),
		NAOK=TRUE, PACKAGE="SNPRelate")

	if (verbose)
	{
		cat("Working space:", node$n.samp, "samples,", node$n.snp, "SNPs\n");
		if (num.thread <= 1)
			cat("\tUsing", num.thread, "CPU core.\n")
		else
			cat("\tUsing", num.thread, "CPU cores.\n")
	}

	# call the C function
	rv <- list(sample.id = sample.ids, snp.id = snp.ids,
		ibs = .Call("gnrIBSAve", as.logical(verbose), TRUE,
			as.integer(num.thread), PACKAGE="SNPRelate"))
	class(rv) <- "snpgdsIBSClass"

	return(rv)
}


#######################################################################
# To calculate the identity-by-state (IBS) matrix for SNP genotypes
#
# INPUT:
#   gdsobj -- an object of gds file
#   sample.id -- a vector of sample.id; if NULL, to use all samples
#   snp.id -- a vector of snp.id; if NULL, to use all SNPs
#   autosome.only -- whether only use autosomal SNPs
#   remove.monosnp -- whether remove monomorphic snps or not
#   maf -- the threshold of minor allele frequencies, keeping ">= maf"
#   missing.rate -- the threshold of missing rates, keeping "<= missing.rate"
#   verbose -- show information
#

snpgdsIBSNum <- function(gdsobj, sample.id=NULL, snp.id=NULL,
	autosome.only=TRUE, remove.monosnp=TRUE, maf=NaN, missing.rate=NaN,
	num.thread=1, verbose=TRUE)
{
	# check
	stopifnot(inherits(gdsobj, "gds.class"))
	stopifnot(is.numeric(num.thread) & (num.thread>0))
	stopifnot(is.logical(verbose))

	# samples
	sample.ids <- read.gdsn(index.gdsn(gdsobj, "sample.id"))
	if (!is.null(sample.id))
	{
		n.tmp <- length(sample.id)
		sample.id <- sample.ids %in% sample.id
		n.samp <- sum(sample.id);
		if (n.samp != n.tmp)
			stop("Some of sample.id do not exist!")
		if (n.samp <= 0)
			stop("No sample in the working dataset.")
		sample.ids <- sample.ids[sample.id]
	}

	if (verbose)
		cat("Identity-By-State (IBS) analysis on SNP genotypes:\n");

	# SNPs
	snp.ids <- read.gdsn(index.gdsn(gdsobj, "snp.id"))
	if (!is.null(snp.id))
	{
		n.tmp <- length(snp.id)
		snp.id <- snp.ids %in% snp.id
		n.snp <- sum(snp.id)
		if (n.snp != n.tmp)
			stop("Some of snp.id do not exist!")
		if (n.snp <= 0)
			stop("No SNP in the working dataset.")
		if (autosome.only)
		{
			chr <- read.gdsn(index.gdsn(gdsobj, "snp.chromosome"))
			opt <- snpgdsOption(gdsobj)
			snp.id <- snp.id & (chr %in% c(opt$autosome.start : opt$autosome.end))
			if (verbose)
			{
				tmp <- n.snp - sum(snp.id)
				if (tmp > 0) cat("Removing", tmp, "non-autosomal SNPs\n")
			}
		}
		snp.ids <- snp.ids[snp.id]
	} else {
		if (autosome.only)
		{
			chr <- read.gdsn(index.gdsn(gdsobj, "snp.chromosome"))
			opt <- snpgdsOption(gdsobj)
			snp.id <- chr %in% c(opt$autosome.start : opt$autosome.end)
			snp.ids <- snp.ids[snp.id]
			if (verbose)
				cat("Removing", length(chr) - length(snp.ids), "non-autosomal SNPs\n")
		}
	}

	# call C codes
	# set genotype working space
	node <- .C("gnrSetGenoSpace", as.integer(index.gdsn(gdsobj, "genotype")),
		as.logical(sample.id), as.logical(!is.null(sample.id)),
		as.logical(snp.id), as.logical(!is.null(snp.id)),
		n.snp=integer(1), n.samp=integer(1),
		err=integer(1), NAOK=TRUE, PACKAGE="SNPRelate")
	if (node$err != 0) stop(snpgdsErrMsg())

	# call allele freq. and missing rates
	if (remove.monosnp || is.finite(maf) || is.finite(missing.rate))
	{
		if (!is.finite(maf)) maf <- -1;
		if (!is.finite(missing.rate)) missing.rate <- 2;
		# call
		rv <- .C("gnrSelSNP_Base", as.logical(remove.monosnp),
			as.double(maf), as.double(missing.rate),
			out.num=integer(1), out.snpflag = logical(node$n.snp),
			err=integer(1), NAOK=TRUE, PACKAGE="SNPRelate")
		if (rv$err != 0) stop(snpgdsErrMsg())
		snp.ids <- snp.ids[rv$out.snpflag]
		# show
		if (verbose)
			cat("Removing", rv$out.num, "SNPs (monomorphic, < MAF, or > missing rate)\n")
	}

	# get the dimension of SNP genotypes
	node <- .C("gnrGetGenoDim", n.snp=integer(1), n.samp=integer(1),
		NAOK=TRUE, PACKAGE="SNPRelate")

	if (verbose)
	{
		cat("Working space:", node$n.samp, "samples,", node$n.snp, "SNPs\n");
		if (num.thread <= 1)
			cat("\tUsing", num.thread, "CPU core.\n")
		else
			cat("\tUsing", num.thread, "CPU cores.\n")
	}

	# call the C function
	rv <- .Call("gnrIBSNum", as.logical(verbose), TRUE,
		as.integer(num.thread), PACKAGE="SNPRelate")
	names(rv) <- c("ibs0", "ibs1", "ibs2")

	# return
	rv <- list(sample.id = sample.ids, snp.id = snp.ids,
		ibs0 = rv$ibs0, ibs1 = rv$ibs1, ibs2 = rv$ibs2)
	return(rv)
}
