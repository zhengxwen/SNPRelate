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
# Linkage Disequilibrium (LD) analysis
#######################################################################

#######################################################################
# To calculate LD for a pair of SNPs
#
# INPUT:
#   snp1 -- an array of snp genotypes at the first locus
#   snp2 -- an array of snp genotypes at the second locus
#   method -- "composite"  Composite LD coefficients (by default)
#             "r"          R coefficient (by EM algorithm assuming HWE)
#             "dprime"     D' coefficient
#             "corr"       Correlation coefficient (BB, AB, AA are codes as 0, 1, 2)
#

snpgdsLDpair <- function(snp1, snp2, method=c("composite", "r", "dprime", "corr"))
{
	# check
	stopifnot(is.integer(snp1))
	stopifnot(is.integer(snp2))
	stopifnot(length(snp1) == length(snp2))
	stopifnot(is.character(method))

	method <- match(method[1], c("composite", "r", "dprime", "corr"))
	if (is.na(method))
		stop("method should be one of \"composite\", \"r\", \"dprime\" and \"corr\"")

	# call
	rv <- .C("gnrLDpair", as.integer(snp1), as.integer(snp2), length(snp1), method,
		out=double(1), pA_A=double(1), pA_B=double(1), pB_A=double(1), pB_B=double(1),
		err=integer(1), NAOK=TRUE, PACKAGE="SNPRelate")
	if (rv$err != 0) stop(snpgdsErrMsg())

	# output
	if (method %in% c(2, 3))
		return(data.frame(ld=rv$out, pA_A=rv$pA_A, pA_B=rv$pA_B, pB_A=rv$pB_A, pB_B=rv$pB_B))
	else
		return(data.frame(ld=rv$out))
}


#######################################################################
# To calculate LD for each pair of SNPs in the region
#
# INPUT:
#   gdsobj -- an object of gds file
#   sample.id -- a vector of sample.id; if NULL, to use all samples
#   snp.id -- a vector of snp.id; if NULL, to use all SNPs
#   slide -- the size of sliding windows,
#            if slide <= 0, no sliding window, output a n-by-n LD matrix
#            otherwise, output a slide-by-n LD matrix
#   method -- "composite"  Composite LD coefficients (by default)
#             "r"          R coefficient (by EM algorithm assuming HWE)
#             "dprime"     D' coefficient
#             "corr"       Correlation coefficient (BB, AB, AA are codes as 0, 1, 2)
#   num.thread -- the number of threads to be used
#   verbose -- show information
#

snpgdsLDMat <- function(gdsobj, sample.id=NULL, snp.id=NULL,
	slide=250, method=c("composite", "r", "dprime", "corr"),
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
		snp.ids <- snp.ids[snp.id]
	}

	# method
	method <- match(method[1], c("composite", "r", "dprime", "corr"))
	if (is.na(method))
		stop("method should be one of \"composite\", \"r\", \"dprime\" and \"corr\"")

	# call C codes
	# set genotype working space
	node <- .C("gnrSetGenoSpace", as.integer(index.gdsn(gdsobj, "genotype")),
		as.logical(sample.id), as.logical(!is.null(sample.id)),
		as.logical(snp.id), as.logical(!is.null(snp.id)),
		n.snp=integer(1), n.samp=integer(1),
		err=integer(1), NAOK=TRUE, PACKAGE="SNPRelate")
	if (node$err != 0) stop(snpgdsErrMsg())

	slide <- as.integer(slide)
	if (is.na(slide)) slide <- as.integer(-1)
	if (slide > node$n.snp) slide <- node$n.snp

	if (verbose)
	{
		cat("Linkage Disequilibrium (LD) analysis on SNP genotypes:\n");
		cat("Working space:", node$n.samp, "samples,", node$n.snp, "SNPs\n");
		if (num.thread <= 1)
			cat("\tUsing", num.thread, "CPU core.\n")
		else
			cat("\tUsing", num.thread, "CPU cores.\n")
		if (slide > 0)
			cat("\tSliding window size:", slide, "\n")
	}

	# call parallel IBD
	rv <- .C("gnrLDMat", method, verbose, TRUE, as.integer(num.thread), slide,
		LD = { if (slide <= 0)
					matrix(NaN, nrow=node$n.snp, ncol=node$n.snp)
				else
					matrix(NaN, nrow=slide, ncol=node$n.snp) },
		err = integer(1), NAOK=TRUE, PACKAGE="SNPRelate")
	if (rv$err != 0) stop(snpgdsErrMsg())

	# return
	rv <- list(sample.id=sample.ids, snp.id=snp.ids, LD=rv$LD, slide=slide)
	return(rv)
}



#######################################################################
# To prune SNPs based on LD
#
# INPUT:
#   gdsobj -- an object of gds file
#   sample.id -- a vector of sample.id; if NULL, to use all samples
#   snp.id -- a vector of snp.id; if NULL, to use all SNPs
#   autosome.only -- whether only use autosomal SNPs
#   remove.monosnp -- whether remove monomorphic snps or not
#   maf -- the threshold of minor allele frequencies, keeping ">= maf"
#   missing.rate -- the threshold of missing rates, keeping "<= missing.rate"
#   method -- "composite"  Composite LD coefficients (by default)
#             "r"          R coefficient (by EM algorithm assuming HWE)
#             "dprime"     D' coefficient
#             "corr"       Correlation coefficient (BB, AB, AA are codes as 0, 1, 2)
#   num.thread -- the number of threads to be used
#   verbose -- show information
#

snpgdsLDpruning <- function(gdsobj, sample.id=NULL, snp.id=NULL,
	autosome.only=TRUE, remove.monosnp=TRUE, maf=NaN, missing.rate=NaN,
	method=c("composite", "r", "dprime", "corr"),
	slide.max.bp=500000, slide.max.n=NA, ld.threshold=0.2,
	num.thread=1, verbose=TRUE)
{
	# check
	stopifnot(inherits(gdsobj, "gds.class"))
	stopifnot(is.na(slide.max.bp) | is.numeric(slide.max.bp))
	stopifnot(is.na(slide.max.n) | is.numeric(slide.max.n))
	stopifnot(is.numeric(ld.threshold) & is.finite(ld.threshold))
	stopifnot(is.numeric(num.thread) & (num.thread>0))
	if (num.thread > 1)
		warning("The current version of 'snpgdsLDpruning' does not support multi processes.")
	stopifnot(is.logical(verbose))

	if (verbose)
	{
		cat("SNP pruning based on LD:\n")
		bp <- slide.max.bp; mn <- slide.max.n
		if (!is.finite(bp)) bp <- Inf
		if (!is.finite(mn)) mn <- Inf
		cat(sprintf("\tSliding window: %g basepairs, %g SNPs\n", bp, mn))
		cat(sprintf("\t|LD| threshold: %g\n", ld.threshold))
	}

	if (!is.finite(slide.max.bp))
		slide.max.bp <- .Machine$double.xmax
	if (!is.finite(slide.max.n))
		slide.max.n <- .Machine$integer.max

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

	# SNPs
	total.snp.ids <- snp.ids <- read.gdsn(index.gdsn(gdsobj, "snp.id"))
	chr <- read.gdsn(index.gdsn(gdsobj, "snp.chromosome"))
	position <- as.integer(read.gdsn(index.gdsn(gdsobj, "snp.position")))
	if (!all(length(snp.ids)==length(chr), length(snp.ids)==length(position)))
		stop("snp.id, snp.chromosome and snp.position should have the same length!")
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
			opt <- snpgdsOption(gdsobj)
			snp.id <- chr %in% c(opt$autosome.start : opt$autosome.end)
			snp.ids <- snp.ids[snp.id]
			if (verbose)
				cat("Removing", length(chr) - length(snp.ids), "non-autosomal SNPs\n")
		}
	}

	# method
	method <- match(method[1], c("composite", "r", "dprime", "corr"))
	if (is.na(method))
		stop("method should be one of \"composite\", \"r\", \"dprime\" and \"corr\"")

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
		cat("Working space:", node$n.samp, "samples,", node$n.snp, "SNPs\n");

	# for-loop each chromosome
	ntotal <- 0; res <- list()
	snp.flag <- total.snp.ids %in% snp.ids

	for (ch in setdiff(unique(chr), c(0, NA)))
	{
		flag <- snp.flag & (chr == ch)
		n.tmp <- sum(flag)
		if (n.tmp > 0)
		{
			# set genotype working space
			node <- .C("gnrSetGenoSpace", as.integer(index.gdsn(gdsobj, "genotype")),
				as.logical(sample.id), as.logical(!is.null(sample.id)),
				flag, TRUE, n.snp=integer(1), n.samp=integer(1),
				err=integer(1), NAOK=TRUE, PACKAGE="SNPRelate")
			if (node$err != 0) stop(snpgdsErrMsg())

			# call LD prune for this chromosome
			startidx <- sample(1:n.tmp, 1)
			rv <- .C("gnrLDpruning", as.integer(startidx-1), position[flag],
				as.integer(slide.max.bp), as.integer(slide.max.n), as.double(ld.threshold),
				method, out_snp = logical(node$n.snp),
				err = integer(1), NAOK=TRUE, PACKAGE="SNPRelate")
			if (rv$err != 0) stop(snpgdsErrMsg())

			# output
			L <- rep(FALSE, length(total.snp.ids))
			L[flag] <- rv$out_snp
			res[[paste("chr", ch, sep="")]] <- total.snp.ids[L]
			ntotal <- ntotal + sum(rv$out_snp)

			# information
			if (verbose)
			{
				ntmp <- sum(rv$out_snp); ntot <- sum(chr == ch)
				cat(sprintf("Chromosome %d: %0.2f%%, %d/%d\n", ch, 100*ntmp/ntot, ntmp, ntot))
			}
		}
	}

	if (verbose)
		cat(sprintf("%d SNPs are selected in total.\n", ntotal))

	# return
	return(res)
}
