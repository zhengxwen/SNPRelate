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
# Principal Component Analysis
#######################################################################

#######################################################################
# To conduct Principal Component Analysis
#
# INPUT:
#   gdsobj -- an object of gds file
#   sample.id -- a vector of sample.id; if NULL, to use all samples
#   snp.id -- a vector of snp.id; if NULL, to use all SNPs
#   autosome.only -- whether only use autosomal SNPs
#   remove.monosnp -- whether remove monomorphic snps or not
#   maf -- the threshold of minor allele frequencies, keeping ">= maf"
#   missing.rate -- the threshold of missing rates, keeping "<= missing.rate"
#   eigen.cnt -- the number of eigenvectors output
#   num.thread -- the number of threads
#   bayesian -- if TRUE, to use Bayesian adjustment
#   need.genmat -- if TRUE, return genetic covariance matrix
#   verbose -- show information, if TRUE
#

snpgdsPCA <- function(gdsobj, sample.id=NULL, snp.id=NULL,
	autosome.only=TRUE, remove.monosnp=TRUE, maf=NaN, missing.rate=NaN,
	eigen.cnt=32, num.thread=1, bayesian=FALSE, need.genmat=FALSE, genmat.only=FALSE,
	verbose=TRUE)
{
	# check
	stopifnot(inherits(gdsobj, "gds.class"))
	stopifnot(is.numeric(num.thread) & (num.thread>0))
	stopifnot(is.numeric(eigen.cnt))
	stopifnot(is.logical(bayesian))
	stopifnot(is.logical(need.genmat))
	stopifnot(is.logical(genmat.only))
	stopifnot(is.logical(verbose))
	if (genmat.only) need.genmat <- TRUE

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
		cat("Principal Component Analysis (PCA) on SNP genotypes:\n");

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
				if (tmp > 0) cat("Removing", tmp, "non-autosomal SNPs.\n")
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
				cat("Removing", length(chr) - length(snp.ids), "non-autosomal SNPs.\n")
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

	if (eigen.cnt <= 0) eigen.cnt <- node$n.samp

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

	# call parallel PCA
	rv <- .C("gnrPCA", as.integer(eigen.cnt), as.integer(num.thread),
		as.logical(bayesian), as.logical(need.genmat), as.logical(genmat.only),
		as.logical(verbose), TRUE, eigenval = double(node$n.samp),
		eigenvect = matrix(NaN, nrow=node$n.samp, ncol=eigen.cnt),
		TraceXTX = double(1),
		genmat = switch(as.integer(need.genmat)+1, double(0),
			matrix(NaN, nrow=node$n.samp, ncol=node$n.samp)),
		err = integer(1), NAOK=TRUE, PACKAGE="SNPRelate")
	if (rv$err != 0) stop(snpgdsErrMsg())

	# return
	rv <- list(sample.id = sample.ids, snp.id = snp.ids,
		eigenval = rv$eigenval, eigenvect = rv$eigenvect, TraceXTX = rv$TraceXTX,
		Bayesian = bayesian, genmat = rv$genmat)
	class(rv) <- "snpgdsPCAClass"
	return(rv)
}



#######################################################################
# To calculate SNP correlations from principal component analysis
#
# INPUT:
#   pcaobj -- a "snpgdsPCAClass" object from the function "snpgds.pca"
#   gdsobj -- an object of gds file
#   snp.id -- a vector of snp.id; if NULL, to use all SNPs
#   eig.which -- specify which eigenvectors to be used
#   num.thread -- the number of threads
#   verbose -- show information
#

snpgdsPCACorr <- function(pcaobj, gdsobj, snp.id=NULL, eig.which=NULL,
	num.thread=1, verbose=TRUE)
{
	# check
	stopifnot(inherits(pcaobj, "snpgdsPCAClass"))
	stopifnot(inherits(gdsobj, "gds.class"))
	stopifnot(is.numeric(num.thread) & (num.thread>0))
	stopifnot(is.logical(verbose))

	# samples
	sample.ids <- read.gdsn(index.gdsn(gdsobj, "sample.id"))
	sample.id <- sample.ids %in% pcaobj$sample.id
	sample.ids <- pcaobj$sample.id

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

	if (is.null(eig.which))
		eig.which <- 1:dim(pcaobj$eigenvect)[2]

	# call C codes
	# set genotype working space
	node <- .C("gnrSetGenoSpace", as.integer(index.gdsn(gdsobj, "genotype")),
		as.integer(sample.id), as.integer(!is.null(sample.id)),
		as.integer(snp.id), as.integer(!is.null(snp.id)),
		n.snp=integer(1), n.samp=integer(1),
		err=integer(1), NAOK=TRUE, PACKAGE="SNPRelate")
	if (node$err != 0) stop(snpgdsErrMsg())

	if (verbose)
	{
		cat("SNP correlations:\n")
		cat("Working space:", node$n.samp, "samples,", node$n.snp, "SNPs\n");
		if (num.thread <= 1)
			cat("\tUsing", num.thread, "CPU core.\n")
		else
			cat("\tUsing", num.thread, "CPU cores.\n")
		cat("\tUsing the top", dim(pcaobj$eigenvect)[2], "eigenvectors.\n")
	}

	# call parallel PCA
	dm <- as.integer(c(dim(pcaobj$eigenvect)[1], length(eig.which)))
	rv <- .C("gnrPCACorr", dm, pcaobj$eigenvect[, eig.which],
		as.integer(num.thread), verbose, TRUE,
		snpcorr = matrix(NaN, nrow=length(eig.which), ncol=node$n.snp),
		err=integer(1), NAOK=TRUE, PACKAGE="SNPRelate")
	if (rv$err != 0) stop(snpgdsErrMsg())

	# return
	rv <- list(sample.id=sample.ids, snp.id=snp.ids, snpcorr=rv$snpcorr)
	return(rv)
}



#######################################################################
# To calculate SNP loadings from principal component analysis
#
# INPUT:
#   pcaobj -- a "snpgdsPCAClass" object from the function "snpgds.pca"
#   gdsobj -- an object of gds file
#   num.thread -- the number of threads
#   verbose -- show information
#

snpgdsPCASNPLoading <- function(pcaobj, gdsobj, num.thread=1, verbose=TRUE)
{
	# check
	stopifnot(inherits(pcaobj, "snpgdsPCAClass"))
	stopifnot(inherits(gdsobj, "gds.class"))
	stopifnot(is.numeric(num.thread) & (num.thread>0))
	stopifnot(is.logical(verbose))

	# samples
	sample.ids <- read.gdsn(index.gdsn(gdsobj, "sample.id"))
	sample.id <- sample.ids %in% pcaobj$sample.id
	sample.ids <- pcaobj$sample.id

	# SNPs
	snp.ids <- read.gdsn(index.gdsn(gdsobj, "snp.id"))
	snp.id <- snp.ids %in% pcaobj$snp.id
	snp.ids <- pcaobj$snp.id

	# call C codes
	# set genotype working space
	node <- .C("gnrSetGenoSpace", as.integer(index.gdsn(gdsobj, "genotype")),
		as.integer(sample.id), as.integer(!is.null(sample.id)),
		as.integer(snp.id), as.integer(!is.null(snp.id)),
		n.snp=integer(1), n.samp=integer(1),
		err=integer(1), NAOK=TRUE, PACKAGE="SNPRelate")
	if (node$err != 0) stop(snpgdsErrMsg())

	if (verbose)
	{
		cat("SNP loadings:\n")
		cat("Working space:", node$n.samp, "samples,", node$n.snp, "SNPs\n");
		if (num.thread <= 1)
			cat("\tUsing", num.thread, "CPU core.\n")
		else
			cat("\tUsing", num.thread, "CPU cores.\n")
		cat("\tUsing the top", dim(pcaobj$eigenvect)[2], "eigenvectors.\n")
	}

	# call parallel PCA
	rv <- .C("gnrPCASNPLoading", pcaobj$eigenval, dim(pcaobj$eigenvect),
		pcaobj$eigenvect, pcaobj$TraceXTX, as.integer(num.thread),
		pcaobj$Bayesian, verbose, TRUE,
		snploading = matrix(NaN, nrow=dim(pcaobj$eigenvect)[2], ncol=node$n.snp),
		afreq=double(node$n.snp), scale=double(node$n.snp),
		err=integer(1), NAOK=TRUE, PACKAGE="SNPRelate")
	if (rv$err != 0) stop(snpgdsErrMsg())

	# return
	rv <- list(sample.id=sample.ids, snp.id=snp.ids, eigenval=pcaobj$eigenval,
		snploading=rv$snploading, TraceXTX=pcaobj$TraceXTX,
		Bayesian=pcaobj$Bayesian, avefreq=rv$afreq, scale=rv$scale)
	class(rv) <- "snpgdsPCASNPLoadingClass"
	return(rv)
}



#######################################################################
# To calculate sample loadings from SNP loadings in principal component analysis
#
# INPUT:
#   loadobj -- a "snpgdsPCASNPLoading" object from the function "snpgds.pca"
#   gdsobj -- an object of gds file
#   sample.id -- a vector of sample.id; if NULL, to use all samples
#   num.thread -- the number of threads
#   verbose -- show information
#

snpgdsPCASampLoading <- function(loadobj, gdsobj, sample.id=NULL,
	num.thread=1, verbose=TRUE)
{
	# check
	stopifnot(inherits(loadobj, "snpgdsPCASNPLoadingClass"))
	stopifnot(inherits(gdsobj, "gds.class"))
	stopifnot(is.numeric(num.thread) & (num.thread>0))
	stopifnot(is.logical(verbose))

	# samples
	sample.ids <- read.gdsn(index.gdsn(gdsobj, "sample.id"))
	if (!is.null(sample.id))
	{
		sample.id <- sample.ids %in% sample.id
		n.samp <- sum(sample.id);
		if (n.samp <= 0)
			stop("No sample in the working dataset.")
		sample.ids <- sample.ids[sample.id]
	}

	# SNPs
	snp.ids <- read.gdsn(index.gdsn(gdsobj, "snp.id"))
	snp.id <- snp.ids %in% loadobj$snp.id
	snp.ids <- loadobj$snp.id

	# call C codes
	# set genotype working space
	node <- .C("gnrSetGenoSpace", as.integer(index.gdsn(gdsobj, "genotype")),
		as.integer(sample.id), as.integer(!is.null(sample.id)),
		as.integer(snp.id), as.integer(!is.null(snp.id)),
		n.snp=integer(1), n.samp=integer(1),
		err=integer(1), NAOK=TRUE, PACKAGE="SNPRelate")
	if (node$err != 0) stop(snpgdsErrMsg())

	eigcnt <- dim(loadobj$snploading)[1]
	if (verbose)
	{
		cat("Sample loadings:\n")
		cat("Working space:", node$n.samp, "samples,", node$n.snp, "SNPs\n");
		if (num.thread <= 1)
			cat("\tUsing", num.thread, "CPU core.\n")
		else
			cat("\tUsing", num.thread, "CPU cores.\n")
		cat("\tUsing the top", eigcnt, "eigenvectors.\n")
	}

	# call parallel PCA
	rv <- .C("gnrPCASampLoading", length(loadobj$sample.id), loadobj$eigenval,
		eigcnt, as.double(loadobj$snploading), loadobj$TraceXTX,
		loadobj$avefreq, loadobj$scale, as.integer(num.thread), verbose, TRUE,
		eigenvect = matrix(NaN, nrow=node$n.samp, ncol=eigcnt),
		err=integer(1), NAOK=TRUE, PACKAGE="SNPRelate")
	if (rv$err != 0) stop(snpgdsErrMsg())

	# return
	rv <- list(sample.id = sample.ids, snp.id = loadobj$snp.ids,
		eigenval = loadobj$eigenval, eigenvect = rv$eigenvect,
		TraceXTX = loadobj$TraceXTX,
		Bayesian = loadobj$Bayesian, genmat = NULL)
	class(rv) <- "snpgdsPCAClass"
	return(rv)
}
