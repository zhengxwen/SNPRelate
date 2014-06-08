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
# Default human chromosome coding:
#          autosome                        -> 1 .. 22
#     X    X chromosome                    -> 23
#     XY   Pseudo-autosomal region of X    -> 24
#     Y    Y chromosome                    -> 25
#     MT   Mitochondrial                   -> 26
#######################################################################





#######################################################################
# Summary Descriptive Statistics
#######################################################################

#######################################################################
# To calculate the missing rate and allele frequency for each SNP
#
# INPUT:
#   gdsobj -- an object of gds file
#   sample.id -- a vector of sample.id; if NULL, to use all samples
#   snp.id -- a vector of snp.id; if NULL, to use all SNPs
#

snpgdsSNPRateFreq <- function(gdsobj, sample.id=NULL, snp.id=NULL)
{
	# check
	stopifnot(inherits(gdsobj, "gds.class"))
	# samples
	if (!is.null(sample.id))
	{
		n.tmp <- length(sample.id)
		sample.id <- read.gdsn(index.gdsn(gdsobj, "sample.id")) %in% sample.id
		n.samp <- sum(sample.id);
		if (n.samp != n.tmp)
			stop("Some of sample.id do not exist!")
		if (n.samp <= 0)
			stop("No sample in the working dataset.")
	}
	# SNPs
	if (!is.null(snp.id))
	{
		n.tmp <- length(snp.id)
		snp.id <- read.gdsn(index.gdsn(gdsobj, "snp.id")) %in% snp.id
		n.snp <- sum(snp.id)
		if (n.snp != n.tmp)
			stop("Some of snp.id do not exist!")
		if (n.snp <= 0)
			stop("No SNP in the working dataset.")
	}

	# call C codes
	# set genotype working space
	rv <- .C("gnrSetGenoSpace", as.integer(index.gdsn(gdsobj, "genotype")),
		as.logical(sample.id), as.logical(!is.null(sample.id)),
		as.logical(snp.id), as.logical(!is.null(snp.id)),
		n.snp=integer(1), n.samp=integer(1),
		err=integer(1), NAOK=TRUE, PACKAGE="SNPRelate")
	if (rv$err != 0) stop(snpgdsErrMsg())

	# call allele freq. and missing rates
	rv <- .C("gnrSNPFreq", AF=double(rv$n.snp), MF=double(rv$n.snp),
		MR=double(rv$n.snp), err=integer(1), NAOK=TRUE, PACKAGE="SNPRelate")
	if (rv$err != 0) stop(snpgdsErrMsg())
	list(AlleleFreq=rv$AF, MinorFreq=rv$MF, MissingRate=rv$MR)
}


#######################################################################
# To calculate the missing rate for each sample
#
# INPUT:
#   gdsobj -- an object of gds file
#   sample.id -- a vector of sample.id; if NULL, to use all samples
#   snp.id -- a vector of snp.id; if NULL, to use all SNPs
#

snpgdsSampMissrate <- function(gdsobj, sample.id=NULL, snp.id=NULL)
{
	# check
	stopifnot(inherits(gdsobj, "gds.class"))
	# samples
	if (!is.null(sample.id))
	{
		n.tmp <- length(sample.id)
		sample.id <- read.gdsn(index.gdsn(gdsobj, "sample.id")) %in% sample.id
		n.samp <- sum(sample.id);
		if (n.samp != n.tmp)
			stop("Some of sample.id do not exist!")
		if (n.samp <= 0)
			stop("No sample in the working dataset.")
	}
	# SNPs
	if (!is.null(snp.id))
	{
		n.tmp <- length(snp.id)
		snp.id <- read.gdsn(index.gdsn(gdsobj, "snp.id")) %in% snp.id
		n.snp <- sum(snp.id)
		if (n.snp != n.tmp)
			stop("Some of snp.id do not exist!")
		if (n.snp <= 0)
			stop("No SNP in the working dataset.")
	}

	# call C codes
	# set genotype working space
	rv <- .C("gnrSetGenoSpace", as.integer(index.gdsn(gdsobj, "genotype")),
		as.logical(sample.id), as.logical(!is.null(sample.id)),
		as.logical(snp.id), as.logical(!is.null(snp.id)),
		n.snp=integer(1), n.samp=integer(1),
		err=integer(1), NAOK=TRUE, PACKAGE="SNPRelate")
	if (rv$err != 0) stop(snpgdsErrMsg())

	# call allele freq. and missing rates
	rv <- .C("gnrSampFreq", MR=double(rv$n.samp), err=integer(1), NAOK=TRUE, PACKAGE="SNPRelate")
	rv$MR
}


#######################################################################
# Return a list of candidate SNPs
#
# INPUT:
#   gdsobj -- an object of gds file
#   sample.id -- a vector of sample.id; if NULL, to use all samples
#   snp.id -- a vector of snp.id; if NULL, to use all SNPs
#   autosome.only -- whether only use autosomal SNPs
#   remove.monosnp -- whether remove monomorphic snps or not
#   maf -- the threshold of minor allele frequencies, keeping ">= maf"
#   missing.rate -- the threshold of missing rates, keeping "<= missing.rate"
#

snpgdsSelectSNP <- function(gdsobj, sample.id=NULL, snp.id=NULL,
	autosome.only=TRUE, remove.monosnp=TRUE, maf=NaN, missing.rate=NaN, verbose=TRUE)
{
	# check
	stopifnot(inherits(gdsobj, "gds.class"))

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

	# output
	return(snp.ids)
}




#######################################################################
# Individual Inbreeding Coefficients
#######################################################################

#######################################################################
# To calculate individual inbreeding coefficient
#
# INPUT:
#   x -- an array of snp genotypes
#   p -- allele frequencies
#   method -- "mom", "mle"
#

snpgdsIndInbCoef <- function(x, p, method=c("mom.weir", "mom.visscher", "mle"),
	reltol=.Machine$double.eps^0.75)
{
	# check
	stopifnot(length(x) == length(p))
	stopifnot(method %in% c("mom.weir", "mom.visscher", "mle"))
	method <- method[1]
	x[!(x %in% c(0,1,2))] <- NA

	if (method == "mom.weir")
	{
		num <- x*x - (1+2*p)*x + 2*p*p
		den <- 2*p*(1-p)
		flag <- is.finite(num) & is.finite(den)
		return(sum(num[flag]) / sum(den[flag]))
	} else if (method == "mom.visscher")
	{
		d <- (x*x - (1+2*p)*x + 2*p*p) / (2*p*(1-p))
		return(mean(d[is.finite(d)]))
	} else if (method == "mle")
	{
		rv <- .C("gnrIndInbCoef", length(x), as.integer(x), as.double(p), as.double(reltol),
			out=double(1), err = integer(1), NAOK=TRUE, PACKAGE="SNPRelate")
		if (rv$err != 0) stop(snpgdsErrMsg())
		return(rv$out)
	}
}


#######################################################################
# To calculate individual inbreeding coefficients
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

snpgdsIndInb <- function(gdsobj, sample.id=NULL, snp.id=NULL,
	autosome.only=TRUE, remove.monosnp=TRUE, maf=NaN, missing.rate=NaN,
	method=c("mom.weir", "mom.visscher", "mle"),
	allele.freq=NULL, out.num.iter=TRUE, reltol=.Machine$double.eps^0.75, verbose=TRUE)
{
	# check
	stopifnot(inherits(gdsobj, "gds.class"))
	stopifnot(method[1] %in% c("mom.weir", "mom.visscher", "mle"))
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
		cat("Estimate individual inbreeding coefficients:\n");

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
		if (!is.null(allele.freq))
			allele.freq <- allele.freq[rv$out.snpflag]
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
	}

	# call allele freq.
	if (is.null(allele.freq))
	{
		rv <- .C("gnrSNPFreq", AF=double(node$n.snp), MF=double(node$n.snp),
			MR=double(node$n.snp), err=integer(1), NAOK=TRUE, PACKAGE="SNPRelate")
		if (rv$err != 0) stop(snpgdsErrMsg())
		allele.freq <- rv$AF
	}

	# call individual inbreeding coefficients
	r <- .C("gnrIndInb", allele.freq,
		as.integer(match(method[1], c("mom.weir", "mom.visscher", "mle"))),
		as.double(reltol), coeff=double(node$n.samp), iternum=integer(node$n.samp),
		err=integer(1), NAOK=TRUE, PACKAGE="SNPRelate")
	if (r$err != 0) stop(snpgdsErrMsg())

	# output
	rv <- list(sample.id = sample.ids, snp.id = snp.ids,
		inbreeding = r$coeff)
	if (out.num.iter & method[1]=="mle")
		rv$out.num.iter <- r$iternum
	return(rv)
}




#######################################################################
# Genetic Distance analysis
#######################################################################

#######################################################################
# To calculate the genetic dissimilarity matrix for SNP genotypes
#
# INPUT:
#   gdsobj -- an object of gds file
#   sample.id -- a vector of sample.id; if NULL, to use all samples
#   snp.id -- a vector of snp.id; if NULL, to use all SNPs
#   autosome.only -- whether only use autosomal SNPs
#   remove.monosnp -- whether remove monomorphic snps or not
#   maf -- the threshold of minor allele frequencies, keeping ">= maf"
#   missing.rate -- the threshold of missing rates, keeping "<= missing.rate"
#   num.thread -- the number of threads to be used
#   verbose -- show information
#

snpgdsDiss <- function(gdsobj, sample.id=NULL, snp.id=NULL, autosome.only=TRUE,
	remove.monosnp=TRUE, maf=NaN, missing.rate=NaN, num.thread=1, verbose=TRUE)
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
		cat("Individual dissimilarity analysis on SNP genotypes:\n");

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
		cat("\tUse", num.thread, "CPU cores.\n")
	}

	# call C function
	rv <- .C("gnrDiss", as.logical(verbose), TRUE, as.integer(num.thread),
		diss = matrix(NaN, ncol=node$n.samp, nrow=node$n.samp),
		err = integer(1), NAOK=TRUE, PACKAGE="SNPRelate")
	if (rv$err != 0) stop(snpgdsErrMsg())

	# return
	rv <- list(sample.id = sample.ids, snp.id = snp.ids, diss = rv$diss)
	class(rv) <- "snpgdsDissClass"
	return(rv)
}



#######################################################################
# To perform hierarchical cluster analysis
#
# INPUT:
#   dist -- the dissimilarity matrix
#   sample.id -- only work if dist is a matrix, to specify sample id
#   need.mat -- if TRUE, store the dissimilarity matrix in the result
#   hang -- see plot.hclust
#

snpgdsHCluster <- function(dist, sample.id=NULL, need.mat=TRUE, hang=0.25)
{
	# check
	stopifnot(is.matrix(dist) | inherits(dist, "snpgdsDissClass") |
		inherits(dist, "snpgdsIBSClass"))

	if (inherits(dist, "snpgdsDissClass"))
	{
		sample.id <- dist$sample.id
		dist <- dist$diss
	}
	if (inherits(dist, "snpgdsIBSClass"))
	{
		sample.id <- dist$sample.id
		dist <- 1 - dist$ibs
	}

	if (is.null(sample.id))
	{
		stopifnot(nrow(dist) == ncol(dist))
		if (is.null(colnames(dist)) | is.null(rownames(dist)))
			stop("Please specify 'sample.id'.")
		sample.id <- colnames(dist)
	} else {
		stopifnot(nrow(dist) == length(sample.id))
		stopifnot(ncol(dist) == length(sample.id))
		colnames(dist) <- sample.id
		rownames(dist) <- sample.id
	}

	# run
	hc <- hclust(as.dist(dist), method="average")
	obj <- as.dendrogram(hc, hang=hang)

	rv <- list(sample.id = sample.id, hclust = hc, dendrogram = obj)
	if (need.mat) rv$dist <- dist
	class(rv) <- "snpgdsHCClass"
	rv
}



#######################################################################
# To determine sub groups of individuals
#
# INPUT:
#   hc -- a "snpgdsHCClass" object
#   verbose -- show information
#

snpgdsCutTree <- function(hc, z.threshold=15, outlier.n=5, n.perm = 5000,
	samp.group=NULL, col.outlier="red", col.list=NULL, pch.outlier=4, pch.list=NULL,
	label.H=FALSE, label.Z=TRUE, verbose=TRUE)
{
	# check
	stopifnot(inherits(hc, "snpgdsHCClass"))
	stopifnot(is.finite(z.threshold))
	stopifnot(is.numeric(n.perm))
	stopifnot(is.logical(label.H))
	stopifnot(is.logical(label.Z))
	stopifnot(is.logical(verbose))
	stopifnot(n.perm >= 50)

	if (is.null(hc$dist))
		stop("`hc' should have a matrix of dissimilarity.")

	auto.cluster <- is.null(samp.group)
	if (verbose)
	{
		if (auto.cluster)
		{
			cat(sprintf("Determine groups by permutation (Z threshold: %g, outlier threshold: %d):\n",
				z.threshold, outlier.n))
		}
	}


	# result
	ans <- list(sample.id = hc$sample.id)
	ans$z.threshold <- z.threshold
	ans$outlier.n <- outlier.n
	ans$samp.order <- hc$hclust$order

	if (!auto.cluster)
	{
		stopifnot(is.factor(samp.group))
		stopifnot(length(samp.group) == length(hc$sample.id))
		merge <- NULL

	} else {
		# determine the number of groups

		# check
		stopifnot(!is.null(hc$dist))
		stopifnot(!is.null(hc$hclust$merge))

		n <- dim(hc$dist)[1]
		rv <- .C("gnrDistPerm", n, as.double(hc$dist),
			as.integer(hc$hclust$merge), as.integer(n.perm), as.double(z.threshold),
			OutZ = double(n-1), OutN1 = integer(n-1), OutN2 = integer(n-1),
			OutGrp = integer(n), err=integer(1), NAOK=TRUE, PACKAGE="SNPRelate")
		if (rv$err != 0) stop(snpgdsErrMsg())

		merge <- data.frame(z=rv$OutZ, n1=rv$OutN1, n2=rv$OutN2)

		if (is.finite(outlier.n))
		{
			tab <- table(rv$OutGrp)
			x <- as.integer(names(tab)[tab <= outlier.n])
			flag <- rv$OutGrp %in% x

			samp.group <- sprintf("G%03d", rv$OutGrp)
			samp.group[flag] <- sprintf("Outlier%03d", rv$OutGrp[flag])
			samp.group <- as.factor(samp.group)

			n.g <- length(tab) - length(x)
			n.o <- length(x)
			nm <- character()

			if (n.g > 0) nm <- c(nm, sprintf("G%03d", 1:n.g))
			if (n.o > 0) nm <- c(nm, sprintf("Outlier%03d", 1:n.o))
			levels(samp.group) <- nm

		} else {
			# not detect outliers

			samp.group <- as.factor(sprintf("G%03d", rv$OutGrp))
			n <- levels(samp.group)
			n <- sprintf("G%03d", 1:length(n))
			levels(samp.group) <- n		
		}
	}


	# create a new dendrogram
	add.point <- function(n, sample.id, subgroup)
	{
		if(is.leaf(n))
		{
			idx <- match(attr(n, "label"), sample.id)
			if (as.character(subgroup[idx]) == "outlier")
			{
				attr(n, "nodePar") <- list(pch=pch.outlier, col=col.outlier)
			} else {
				idx <- as.integer(subgroup[idx])
				attr(n, "nodePar") <- list(pch=pch.list[idx], col=col.list[idx])
			}
		} else 
		{
			if (!is.null(merge))
			{
				if (label.H | label.Z)
				{
					x1 <- attr(n[[1]], "members"); x2 <- attr(n[[2]], "members")
					w1 <- (x1 == merge$n1) & (x2 == merge$n2)
					w2 <- (x2 == merge$n1) & (x1 == merge$n2)
					w <- which(w1 | w2)
					if (length(w) > 0)
					{
						if (merge$z[w[1]] >= z.threshold)
						{
							if (label.H)
							{
								if (label.Z)
								{
									attr(n, "edgetext") <- sprintf("H: %0.2f, Z: %0.1f",
										attr(n, "height"), merge$z[w[1]])
								} else {
									attr(n, "edgetext") <- sprintf("H: %0.2f",
										attr(n, "height"))
								}
							} else {
								if (label.Z)
								{
									attr(n, "edgetext") <- sprintf("Z: %0.1f", merge$z[w[1]])
								}
							}
						}
					}
				}
			}
		}
		n
	}

	ans$samp.group <- samp.group

	n <- nlevels(samp.group); grps <- levels(samp.group)
	dmat <- matrix(0.0, nrow=n, ncol=n)
	for (i in 1:n)
	{
		m <- hc$dist[grps[i] == samp.group, grps[i] == samp.group]
		ms <- sum(grps[i] == samp.group)
		mflag <- matrix(TRUE, nrow=ms, ncol=ms)
		diag(mflag) <- FALSE
		dmat[i, i] <- mean(m[mflag], na.rm=TRUE)
		
		if (i < n)
		{
			for (j in (i+1):n)
			{
				dmat[i, j] <- dmat[j, i] <-
					mean(hc$dist[grps[i] == samp.group, grps[j] == samp.group], na.rm=TRUE)
			}
		}
	}
	colnames(dmat) <- rownames(dmat) <- grps
	ans$dmat <- dmat

	if (is.null(pch.list))
	{
		pch.list <- rep(20, n)
	} else {
		pch.list <- rep(pch.list, (n %/% length(pch.list))+1)
	}
	if (is.null(col.list))
	{
		col.list <- 1:n
	} else {
		col.list <- rep(col.list, (n %/% length(col.list))+1)
	}

	ans$dendrogram <- dendrapply(hc$dendrogram, add.point, sample.id=hc$sample.id,
		subgroup=samp.group)
	ans$merge <- merge

	if (!is.null(merge))
	{
		cluster <- samp.group[hc$hclust$order]
		ans$clust.count <- table(cluster)[unique(cluster)]
	} else {
		ans$clust.count <- NULL
	}

	if (verbose)
		cat(sprintf("Create %d groups.\n", n))

	ans
}


#######################################################################
# Draw a dendrogram
#
# INPUT:
#   hc -- a "snpgdsHCClass" object
#   verbose -- show information
#

snpgdsDrawTree <- function(obj, clust.count=NULL, dend.idx=NULL,
	type=c("dendrogram", "z-score"), yaxis.height=TRUE, yaxis.kinship=TRUE,
	y.kinship.baseline=NaN, y.label.kinship=FALSE,
	outlier.n=NULL, shadow.col = c(rgb(0.5, 0.5, 0.5, 0.25), rgb(0.5, 0.5, 0.5, 0.05)),
	outlier.col=rgb(1, 0.50, 0.50, 0.5), leaflab="none", labels=NULL, y.label=0.2,
	...)
{
	# check
	stopifnot(is.null(dend.idx) | is.numeric(dend.idx))
	type <- match.arg(type)
	stopifnot(is.numeric(y.kinship.baseline))

	if (type == "dendrogram")
	{
		stopifnot(!is.null(obj$dendrogram))
		stopifnot(is.null(outlier.n) | is.numeric(outlier.n))

		# initialize ...
		if (is.null(clust.count))
			clust.count <- obj$clust.count
		if (is.null(outlier.n))
			outlier.n <- obj$outlier.n

		# draw
		if (!is.null(dend.idx))
		{
			den <- obj$dendrogram[[dend.idx]]
			x.offset <- 0
			for (i in 1:length(dend.idx))
			{
				if (dend.idx[i] == 2)
				{
					IX <- dend.idx[1:i]
					IX[i] <- 1
					x.offset <- x.offset + attr(obj$dendrogram[[IX]], "member")
				}
			}
		} else {
			den <- obj$dendrogram
			x.offset <- 0
		}

		par(mar=c(4, 4, 4, 4))
		oldpar <- par(mgp=c(5, 1, 0))  # to avoid ylab
		plot(den, leaflab=leaflab, axes=F, ...)
		par(oldpar)

		# draw left y-axis
		if (yaxis.height)
		{
			axis(side=2, line=0)
			tmp <- list(...)
			if (!is.null(tmp$ylab))
				ylab <- tmp$ylab
			else
				ylab <- "individual dissimilarity"
			mtext(ylab, side=2, line=2.5)
		}

		# draw right y-axis
		if (yaxis.kinship)
		{
			if (is.finite(y.kinship.baseline))
			{
				y.kinship.baseline <- y.kinship.baseline[1]
			} else {
				y.kinship.baseline <- attr(den, "height")
			}

			ym <- pretty(c(0, 1))
			axis(side=4, (1 - ym) * y.kinship.baseline, ym, line=0)
			mtext("kinship coefficient", 4, line=2.5)
		}

		# draw others
		if (!is.null(clust.count))
		{
			m <- c(0, cumsum(clust.count))
			jj <- 1; k <- 1
			for (i in 1:length(clust.count))
			{
				if (clust.count[i] > outlier.n)
				{
					rect(m[i] + 0.5 - x.offset, par("usr")[3L],
						m[i+1] + 0.5 - x.offset, par("usr")[4L],
						col = shadow.col[jj], border = NA)
					jj <- 3 - jj
					if (!is.null(labels[k]))
						text((m[i]+m[i+1])/2 - x.offset, y.label, labels[k])
					k <- k + 1
				} else {
					rect(m[i] + 0.5 - x.offset, par("usr")[3L],
						m[i+1] + 0.5 - x.offset, par("usr")[4L],
						col = outlier.col, border = NA)
				}
			}
		}

		# draw kinship label
		if (yaxis.kinship & y.label.kinship)
		{
			# identical twins
			h1 <- (1 - 0.5)*y.kinship.baseline
			abline(h=h1, lty=2, col="gray")

			# parent–child / full-siblings
			h2 <- (1 - 0.25)*y.kinship.baseline
			abline(h=h2, lty=2, col="gray")

			# parent–child / full-siblings
			h3 <- (1 - 1/8)*y.kinship.baseline
			abline(h=h3, lty=2, col="gray")

			# first cousins
			h4 <- (1 - 1/16)*y.kinship.baseline
			abline(h=h4, lty=2, col="gray")

			axis(side=4, c(h1, h2, h3, h4), c("twins", "PC/FS", "DFC/HS", "FC"),
				tick=FALSE, line=-0.75, las=2, cex.axis=0.75, col.axis="gray25")
		}

	} else if (type == "z-score")
	{
		# the distribution of Z scores
		if (is.null(obj$merge))
			stop("There is no Z score in this object.")

		y <- obj$merge[,1]
		y <- y[order(y, decreasing=TRUE)]

		plot(y, xlab="the order of Z score", ylab="Z score", type="b", pch="+", log="x", ...)
		abline(h=15, col="gray", lty=2)
	}

	invisible()
}







#######################################################################
# SNP functions
#######################################################################

#######################################################################
# To get a list of SNP information including rs, chr, pos, allele
#   and allele frequency
#
# INPUT:
#   gdsobj -- an object of gds file
#   sample.id -- a set of sample used in calculating the allele frequencies
#

snpgdsSNPList <- function(gdsobj, sample.id=NULL)
{
	# check
	stopifnot(inherits(gdsobj, "gds.class"))

	# rs id
	if (is.null(index.gdsn(gdsobj, "snp.rs.id", silent=TRUE)))
		rs.id <- read.gdsn(index.gdsn(gdsobj, "snp.id"))
	else
		rs.id <- read.gdsn(index.gdsn(gdsobj, "snp.rs.id"))
	# chromosome
	chromosome <- read.gdsn(index.gdsn(gdsobj, "snp.chromosome"))
	# position
	position <- read.gdsn(index.gdsn(gdsobj, "snp.position"))
	# allele
	allele <- read.gdsn(index.gdsn(gdsobj, "snp.allele"))
	# allele freq.
	afreq <- snpgdsSNPRateFreq(gdsobj, sample.id=sample.id)$AlleleFreq

	rv <- list(rs.id = rs.id, chromosome = chromosome,
		position = position, allele = allele, afreq = afreq)
	class(rv) <- "snpgdsSNPListClass"
	return(rv)
}


#######################################################################
# To get a common list of SNPs from SNP objects, and return
# snp alleles from the first snp object
#
# INPUT:
#   snplist1 -- the first object of snpgdsSNPListClass
#   snplist2 -- the second object of snpgdsSNPListClass
#

snpgdsSNPListIntersect <- function(snplist1, snplist2)
{
	# check
	stopifnot(inherits(snplist1, "snpgdsSNPListClass"))
	stopifnot(inherits(snplist2, "snpgdsSNPListClass"))

	s1 <- paste(snplist1$rs.id, snplist1$chromosome, snplist1$position, sep="-")
	s2 <- paste(snplist2$rs.id, snplist2$chromosome, snplist2$position, sep="-")
	s <- intersect(s1, s2)
	flag <- s1 %in% s

	rv <- list(rs.id = snplist1$rs.id[flag],
		chromosome = snplist1$chromosome[flag], position = snplist1$position[flag],
		allele = snplist1$allele[flag], afreq = snplist1$afreq[flag])
	class(rv) <- "snpgdsSNPListClass"
	return(rv)
}


#######################################################################
# To get a vector of logical variables, indicating whether genotypes
# need to be converted in snplist2.
#
# INPUT:
#   snplist1 -- the first object of snpgdsSNPListClass
#   snplist2 -- the second object of snpgdsSNPListClass
#

snpgdsSNPListStrand <- function(snplist1, snplist2, same.strand=FALSE)
{
	# check
	stopifnot(inherits(snplist1, "snpgdsSNPListClass"))
	stopifnot(inherits(snplist2, "snpgdsSNPListClass"))
	stopifnot(is.logical(same.strand))

	s1 <- paste(snplist1$rs.id, snplist1$chromosome, snplist1$position, sep="-")
	s2 <- paste(snplist2$rs.id, snplist2$chromosome, snplist2$position, sep="-")
	s <- intersect(s1, s2)
	I1 <- match(s, s1); I2 <- match(s, s2)

	# call
	rv <- .C("gnrAlleleStrand", snplist1$allele, snplist1$afreq, I1,
		snplist2$allele, snplist2$afreq, I2,
		same.strand, length(s), out=logical(length(s)),
		out.n.ambiguity=integer(1), out.n.mismatching=integer(1),
		err=integer(1), NAOK=TRUE, PACKAGE="SNPRelate")
	if (rv$err != 0) stop(snpgdsErrMsg())

	# result
	res <- logical(length(s2)); res[I2] <- rv$out
	res[setdiff(1:length(s2), I2)] <- NA
	return(res)
}




#######################################################################
# GDS file management
#######################################################################

#######################################################################
# To summarize the gds file
#
# INPUT:
#   gds -- an object of gds file, or a file names
#   show -- print information on screen
#

snpgdsSummary <- function(gds, show=TRUE)
{
	# check
	stopifnot(inherits(gds, "gds.class") | is.character(gds))

	# open ...
	if (is.character(gds))
	{
		gds.tmp <- openfn.gds(gds)
	} else {
		gds.tmp <- gds
	}

	#
	# checking ...
	#

	warn.flag <- FALSE

	# check sample id
	dm <- objdesp.gdsn(index.gdsn(gds.tmp, "sample.id"))$dim
	if (length(dm) != 1)
	{
		print(gds.tmp)
		if (is.character(gds)) closefn.gds(gds.tmp)
		stop("Invalid dimension of `sample.id'.")
	}
	samp.id <- read.gdsn(index.gdsn(gds.tmp, "sample.id"))
	if (length(samp.id) != length(unique(samp.id)))
	{
		warning("sample.id is not unique!")
		samp.id <- unique(samp.id)
		warn.flag <- TRUE
	}
	n.samp <- dm[1]

	# check snp id
	dm <- objdesp.gdsn(index.gdsn(gds.tmp, "snp.id"))$dim
	if (length(dm) != 1)
	{
		print(gds.tmp)
		if (is.character(gds)) closefn.gds(gds.tmp)
		stop("Invalid dimension of `snp.id'.")
	}
	n.snp <- dm[1]
	snp.id <- read.gdsn(index.gdsn(gds.tmp, "snp.id"))
	if (length(snp.id) != length(unique(snp.id)))
	{
		warning("snp.id is not unique!")
		warn.flag <- TRUE
		snp.flag <- rep(FALSE, n.snp)
		snp.flag[match(unique(snp.id), snp.id)] <- TRUE
	} else {
		snp.flag <- rep(TRUE, n.snp)
	}

	# check snp position
	dm <- objdesp.gdsn(index.gdsn(gds.tmp, "snp.position"))$dim
	if ((length(dm) != 1) | (dm[1] != n.snp))
	{
		print(gds.tmp)
		if (is.character(gds)) closefn.gds(gds.tmp)
		stop("Invalid dimension of `snp.position'.")
	}
	snp.pos <- read.gdsn(index.gdsn(gds.tmp, "snp.position"))
	snp.pos[!is.finite(snp.pos)] <- -1
	if (any(snp.pos <= 0))
	{
		message("Some values of snp.position are invalid (should be > 0)!")
		warn.flag <- TRUE
		snp.flag <- snp.flag & (snp.pos > 0)
	}

	# check snp chromosome
	dm <- objdesp.gdsn(index.gdsn(gds.tmp, "snp.chromosome"))$dim
	if ((length(dm) != 1) | (dm[1] != n.snp))
	{
		print(gds.tmp)
		if (is.character(gds)) closefn.gds(gds.tmp)
		stop("Invalid dimension of `snp.chromosome'.")
	}
	snp.chr <- read.gdsn(index.gdsn(gds.tmp, "snp.chromosome"))
	snp.chr[!is.finite(snp.chr)] <- -1
	flag <- (snp.chr >= 1)
	if (any(!flag))
	{
		message("Some values of snp.chromosome are invalid (should be finite and >= 1)!")
		warn.flag <- TRUE
		snp.flag <- snp.flag & flag
	}

	# check snp allele
	if (!is.null(index.gdsn(gds.tmp, "snp.allele", silent=TRUE)))
	{
		dm <- objdesp.gdsn(index.gdsn(gds.tmp, "snp.allele"))$dim
		if ((length(dm) != 1) | (dm[1] != n.snp))
		{
			print(gds.tmp)
			if (is.character(gds)) closefn.gds(gds.tmp)
			stop("Invalid dimension of `snp.allele'.")
		}
		snp.allele <- read.gdsn(index.gdsn(gds.tmp, "snp.allele"))
		snp.allele[is.na(snp.allele)] <- "?/?"
		flag <- sapply(strsplit(snp.allele, "/"),
			function(x)
			{
				if (length(x) == 2)
				{
					all(x %in% c("A", "G", "C", "T"))
				} else {
					FALSE
				}
			}
		)
		if (any(!flag))
		{
			s <- as.character((snp.allele[!flag])[1])
			message(sprintf("Some of snp.allele are not standard! E.g., %s", s))
			warn.flag <- TRUE
			snp.flag <- snp.flag & flag
		}
	}

	# check genotype
	dm <- objdesp.gdsn(index.gdsn(gds.tmp, "genotype"))$dim
	if (length(dm) != 2)
	{
		print(gds.tmp)
		if (is.character(gds)) closefn.gds(gds.tmp)
		stop("Invalid dimension of `genotype'.")
	}
	lv <- get.attr.gdsn(index.gdsn(gds.tmp, "genotype"))
	if ("sample.order" %in% names(lv))
	{
		if (dm[1]!=n.samp | dm[2]!=n.snp)
		{
			print(gds.tmp)
			if (is.character(gds)) closefn.gds(gds.tmp)
			stop("Invalid dimension of `genotype'.")
		}
	} else {
		if (dm[2]!=n.samp | dm[1]!=n.snp)
		{
			print(gds.tmp)
			if (is.character(gds)) closefn.gds(gds.tmp)
			stop("Invalid dimension of `genotype'.")
		}
	}

	# print
	if (show)
	{
		if (inherits(gds, "gds.class"))
			cat("The file name:", gds$filename, "\n")
		else
			cat("The file name:", gds, "\n")
		cat("The total number of samples:", n.samp, "\n")
		cat("The total number of SNPs:", n.snp, "\n")
		if ("sample.order" %in% names(lv))
		{
			cat("SNP genotypes are stored in SNP-major mode.\n")
		} else {
			cat("SNP genotypes are stored in individual-major mode.\n")
		}
		if (warn.flag)
		{
			cat("The number of valid samples:", length(samp.id), "\n")
			cat("The number of valid SNPs:", sum(snp.flag), "\n")
		}
	}

#	if (warn.flag)
#	{
#		warning("Call `snpgdsCreateGenoSet' to create a valid set of genotypes, using the returned sample.id and snp.id.")
#	}

	warn.flag <- FALSE
	snp.chr <- snp.chr[snp.flag]
	snp.pos <- snp.pos[snp.flag]
	chrtab <- setdiff(unique(snp.chr), 0)
	for (chr in chrtab)
	{
		pos <- snp.pos[snp.chr == chr]
		if (length(pos) > 0)
		{
			if (any(order(pos) != seq_len(length(pos))))
			{
				message(sprintf("The SNP positions are not in ascending order on chromosome %d.", chr))
				warn.flag <- TRUE
				break
			}
		}
	}

	# check -- sample annotation
	if (!is.null(index.gdsn(gds.tmp, "sample.annot", silent=TRUE)))
	{
		# sample.id
		if (!is.null(index.gdsn(gds.tmp, "sample.annot/sample.id", silent=TRUE)))
		{
			s <- read.gdsn(index.gdsn(gds.tmp, "sample.annot/sample.id"))
			if (length(s) != length(samp.id))
				warning("Invalid length of `sample.annot/sample.id'.")
			if (any(s != samp.id))
				warning("Invalid `sample.annot/sample.id'.")
		}

		# others
		lst <- ls.gdsn(index.gdsn(gds.tmp, "sample.annot"))
		for (n in lst)
		{
			dm <- objdesp.gdsn(index.gdsn(gds.tmp, index=c("sample.annot", n)))$dim
			if (!is.null(dm))
			{
				if (dm[1] != length(samp.id))
					warning(sprintf("Invalid `sample.annot/%s'.", n))
			}
		}
	}

	# check -- snp annotation
	if (!is.null(index.gdsn(gds.tmp, "snp.annot", silent=TRUE)))
	{
		if (!is.null(index.gdsn(gds.tmp, "snp.annot/snp.id", silent=TRUE)))
		{
			s <- read.gdsn(index.gdsn(gds.tmp, "snp.annot/snp.id"))
			if (length(s) != length(snp.id))
				warning("Invalid length of `snp.annot/snp.id'.")
			if (any(s != snp.id))
				warning("Invalid `snp.annot/snp.id'.")
		}

		# others
		lst <- ls.gdsn(index.gdsn(gds.tmp, "snp.annot"))
		for (n in lst)
		{
			dm <- objdesp.gdsn(index.gdsn(gds.tmp, index=c("snp.annot", n)))$dim
			if (!is.null(dm))
			{
				if (dm[1] != length(snp.id))
					warning(sprintf("Invalid `snp.annot/%s'.", n))
			}
		}
	}

	# close ...
	if (is.character(gds)) closefn.gds(gds.tmp)

	invisible(list(sample.id = samp.id, snp.id = snp.id[snp.flag]))
}


#######################################################################
# To get a subset of genotypes from a gds file
#
# INPUT:
#   gdsobj -- an object of gds file
#   sample.id -- a vector of sample.id; if NULL, to use all samples
#   snp.id -- a vector of snp.id; if NULL, to use all SNPs
#   snpfirstdim -- if TRUE, indicate store in individual-major order; NULL, by default
#   verbose -- show information
#

snpgdsGetGeno <- function(gdsobj, sample.id=NULL, snp.id=NULL,
	snpfirstdim=NULL, verbose=TRUE)
{
	# check
	stopifnot(inherits(gdsobj, "gds.class"))
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

	# call C codes
	# set genotype working space
	node <- .C("gnrSetGenoSpace", as.integer(index.gdsn(gdsobj, "genotype")),
		as.logical(sample.id), as.logical(!is.null(sample.id)),
		as.logical(snp.id), as.logical(!is.null(snp.id)),
		n.snp=integer(1), n.samp=integer(1),
		err=integer(1), NAOK=TRUE, PACKAGE="SNPRelate")
	if (node$err != 0) stop(snpgdsErrMsg())

	# get the dimension of SNP genotypes
	node <- .C("gnrGetGenoDim", n.snp=integer(1), n.samp=integer(1),
		NAOK=TRUE, PACKAGE="SNPRelate")

	# snp order
	if (is.null(snpfirstdim))
	{
		snpfirstdim <- TRUE
		rd <- names(get.attr.gdsn(index.gdsn(gdsobj, "genotype")))
		if ("snp.order" %in% rd) snpfirstdim <- TRUE
		if ("sample.order" %in% rd) snpfirstdim <- FALSE
	}
	if (snpfirstdim)
	{
		n1 <- node$n.snp; n2 <- node$n.samp
	} else {
		n1 <- node$n.samp; n2 <- node$n.snp
	}

	if (verbose)
	{
		cat("genotype matrix:", node$n.samp, "samples,", node$n.snp, "SNPs\n");
		cat("Whether the SNP is the first dimension: ", snpfirstdim, "\n", sep="");
	}

	# get the dimension of SNP genotypes
	rv <- .C("gnrCopyGenoMem",
		geno = matrix(as.integer(0), nrow=n1, ncol=n2),
		snpfirstdim, err=integer(1), NAOK=TRUE, PACKAGE="SNPRelate")
	if (rv$err != 0) stop(snpgdsErrMsg())

	rv$geno
}


#######################################################################
# To create a gds file for SNP genotypes
#
# INPUT:
#   gds.fn -- the file name of SNP genotypes
#   genmat -- a genotype matrix
#   sample.id -- a vector of sample id (should be unique)
#   snp.id -- a vector of snp id (should be unique)
#   snp.rs.id -- rs.id for snps
#   snp.chromosome -- the chromosome indeces (1..22)
#   snp.position -- the positions in basepair
#   snp.allele -- the reference/non-reference alleles
#   snpfirstdim -- if TRUE, indicate store in individual-major order;
#   compress.annotation -- the compression method for sample and snp annotations
#   compress.geno -- the compression method for genotypes
#   other.vars -- a list of variables
#

snpgdsCreateGeno <- function(gds.fn, genmat,
	sample.id=NULL, snp.id=NULL, snp.rs.id=NULL, snp.chromosome=NULL, snp.position=NULL,
	snp.allele=NULL, snpfirstdim=TRUE, compress.annotation="ZIP.max", compress.geno="",
	other.vars=NULL)
{
	# check
	stopifnot(is.matrix(genmat))
	stopifnot(is.numeric(genmat))
	if (snpfirstdim)
	{
		n.snp <- dim(genmat)[1]; n.samp <- dim(genmat)[2]
	} else {
		n.snp <- dim(genmat)[2]; n.samp <- dim(genmat)[1]
	}

	if (!is.null(sample.id))
	{
		stopifnot(n.samp == length(sample.id))
		if (anyDuplicated(sample.id) > 0) stop("sample.id is not unique!")
	} else
		sample.id <- 1:n.samp
	if (!is.null(snp.id))
	{
		stopifnot(n.snp == length(snp.id))
		if (anyDuplicated(snp.id) > 0) stop("snp.id is not unique!")
	} else
		snp.id <- 1:n.snp
	if (!is.null(snp.rs.id))
		stopifnot(n.snp == length(snp.rs.id))
	if (!is.null(snp.chromosome))
		stopifnot(n.snp == length(snp.chromosome))
	if (!is.null(snp.position))
		stopifnot(n.snp == length(snp.position))
	stopifnot(is.null(other.vars) | is.list(other.vars))

	# create a gds file
	gfile <- createfn.gds(gds.fn)
	add.gdsn(gfile, "sample.id", sample.id, compress=compress.annotation, closezip=TRUE)

	add.gdsn(gfile, "snp.id", snp.id, compress=compress.annotation, closezip=TRUE)
	if (!is.null(snp.rs.id))
		add.gdsn(gfile, "snp.rs.id", snp.rs.id, compress=compress.annotation, closezip=TRUE)

	if (is.null(snp.position))
		snp.position <- as.integer(1:n.snp)
	add.gdsn(gfile, "snp.position", snp.position, compress=compress.annotation, closezip=TRUE)

	if (is.null(snp.chromosome))
		snp.chromosome <- as.integer(rep(1, n.snp))
	add.gdsn(gfile, "snp.chromosome", snp.chromosome, compress=compress.annotation, closezip=TRUE)

	if (!is.null(snp.allele))
		add.gdsn(gfile, "snp.allele", snp.allele, compress=compress.annotation, closezip=TRUE)

	# add genotype
	genmat[is.na(genmat)] <- 3
	genmat[!(genmat %in% c(0,1,2))] <- 3
	node.geno <- add.gdsn(gfile, "genotype", genmat, storage="bit2",
		compress=compress.geno)
	if (snpfirstdim)
		put.attr.gdsn(node.geno, "snp.order")
	else
		put.attr.gdsn(node.geno, "sample.order")

	# other variables
	if (!is.null(other.vars))
	{
		for (i in 1:length(other.vars))
		{
			nm <- names(other.vars)[i]
			add.gdsn(gfile, nm, val=other.vars[[i]], compress=compress.annotation)
		}
	}

	# close the gds file
	closefn.gds(gfile)

	# return
	return(invisible(NULL))
}


#######################################################################
# To create a gds file from a specified gds file
#
# INPUT:
#   src.fn -- the file name of source gds file
#   dest.fn -- the file name of destination gds file
#   sample.id -- a vector of sample id (should be unique)
#   snp.id -- a vector of snp id (should be unique)
#   snpfirstdim -- if TRUE, indicate store in individual-major order; NULL, by default
#   compress.annotation -- the compression method for sample and snp annotations
#   compress.geno -- the compression method for genotypes
#   verbose -- show information
#

snpgdsCreateGenoSet <- function(src.fn, dest.fn, sample.id=NULL, snp.id=NULL,
	snpfirstdim=NULL, compress.annotation="ZIP.max", compress.geno="", verbose=TRUE)
{
	# check
	stopifnot(is.character(src.fn))
	stopifnot(is.character(dest.fn))
	stopifnot(is.logical(snpfirstdim) | is.null(snpfirstdim))

	if (verbose)
		cat("Create a GDS genotype file:\n");

	# open and create gds files
	srcobj <- openfn.gds(src.fn)
	destobj <- createfn.gds(dest.fn)

	# samples
	sample.ids <- read.gdsn(index.gdsn(srcobj, "sample.id"))
	if (!is.null(sample.id))
	{
		n.tmp <- length(sample.id)
		sample.id <- sample.ids %in% sample.id
		n.samp <- sum(sample.id);
		if (n.samp != n.tmp)
		{
			closefn.gds(srcobj); closefn.gds(destobj)
			stop("Some of sample.id do not exist!")
		}
		if (n.samp <= 0)
		{
			closefn.gds(srcobj); closefn.gds(destobj)
			stop("No sample in the working dataset.")
		}
		sample.ids <- sample.ids[sample.id]
	}
	# SNPs
	snp.ids <- read.gdsn(index.gdsn(srcobj, "snp.id"))
	if (!is.null(snp.id))
	{
		n.tmp <- length(snp.id)
		snp.id <- snp.ids %in% snp.id
		n.snp <- sum(snp.id)
		if (n.snp != n.tmp)
		{
			closefn.gds(srcobj); closefn.gds(destobj)
			stop("Some of snp.id do not exist!")
		}
		if (n.snp <= 0)
		{
			closefn.gds(srcobj); closefn.gds(destobj)
			stop("No SNP in the working dataset.")
		}
		snp.ids <- snp.ids[snp.id]
	}

	# call C codes
	# set genotype working space
	node <- .C("gnrSetGenoSpace", as.integer(index.gdsn(srcobj, "genotype")),
		as.logical(sample.id), as.logical(!is.null(sample.id)),
		as.logical(snp.id), as.logical(!is.null(snp.id)),
		n.snp=integer(1), n.samp=integer(1),
		err=integer(1), NAOK=TRUE, PACKAGE="SNPRelate")
	if (node$err != 0) stop(snpgdsErrMsg())

	if (verbose)
	{
		cat(sprintf("The new dataset consists of %d samples and %d SNPs\n",
			length(sample.ids), length(snp.ids)))
	}
	# write to the destination file
	# sample.id
	add.gdsn(destobj, "sample.id", sample.ids, compress=compress.annotation, closezip=TRUE)
	if (verbose) cat("\twrite sample.id\n");
	# snp.id
	add.gdsn(destobj, "snp.id", snp.ids, compress=compress.annotation, closezip=TRUE)
	if (verbose) cat("\twrite snp.id\n");
	# snp.rs.id
	if (!is.null(index.gdsn(srcobj, "snp.rs.id", silent=TRUE)))
	{
		rs.id <- read.gdsn(index.gdsn(srcobj, "snp.rs.id"))
		if (!is.null(snp.id)) rs.id <- rs.id[snp.id]
		add.gdsn(destobj, "snp.rs.id", rs.id, compress=compress.annotation, closezip=TRUE)
		if (verbose) cat("\twrite snp.rs.id\n");
	}
	# snp.position
	pos <- read.gdsn(index.gdsn(srcobj, "snp.position"))
	if (!is.null(snp.id)) pos <- pos[snp.id]
	add.gdsn(destobj, "snp.position", pos, compress=compress.annotation, closezip=TRUE)
	if (verbose) cat("\twrite snp.position\n");
	# snp.chromosome
	chr <- read.gdsn(index.gdsn(srcobj, "snp.chromosome"))
	if (!is.null(snp.id)) chr <- chr[snp.id]
	add.gdsn(destobj, "snp.chromosome", chr, compress=compress.annotation, closezip=TRUE)
	if (verbose) cat("\twrite snp.chromosome\n");
	# snp.allele
	if (!is.null(index.gdsn(srcobj, "snp.allele", silent=TRUE)))
	{
		allele <- read.gdsn(index.gdsn(srcobj, "snp.allele"))
		if (!is.null(snp.id)) allele <- allele[snp.id]
		add.gdsn(destobj, "snp.allele", allele, compress=compress.annotation, closezip=TRUE)
		if (verbose) cat("\twrite snp.allele\n");
	}

	# snp order
	if (is.null(snpfirstdim))
	{
		snpfirstdim <- TRUE
		rd <- names(get.attr.gdsn(index.gdsn(srcobj, "genotype")))
		if ("snp.order" %in% rd) snpfirstdim <- TRUE
		if ("sample.order" %in% rd) snpfirstdim <- FALSE
	}
	if (verbose)
	{
		if (snpfirstdim)
		{
			cat("SNP genotypes are stored in individual-major mode.\n")
		} else {
			cat("SNP genotypes are stored in SNP-major mode.\n")
		}
	}

	# write genotypes to the destination file
	# set genotype working space
	node <- .C("gnrSetGenoSpace", as.integer(index.gdsn(srcobj, "genotype")),
		as.logical(sample.id), as.logical(!is.null(sample.id)),
		as.logical(snp.id), as.logical(!is.null(snp.id)),
		n.snp=integer(1), n.samp=integer(1),
		err=integer(1), NAOK=TRUE, PACKAGE="SNPRelate")
	if (node$err != 0) stop(snpgdsErrMsg())

	if (snpfirstdim)
	{
		gGeno <- add.gdsn(destobj, "genotype", storage="bit2",
			valdim=c(node$n.snp, node$n.samp), compress="")
		put.attr.gdsn(gGeno, "snp.order")
	} else {
		gGeno <- add.gdsn(destobj, "genotype", storage="bit2",
			valdim=c(node$n.samp, node$n.snp), compress="")
		put.attr.gdsn(gGeno, "sample.order")
	}

	rv <- .C("gnrCopyGeno", as.integer(gGeno),
		snpfirstdim, err=integer(1), NAOK=TRUE, PACKAGE="SNPRelate")
	if (rv$err != 0) stop(snpgdsErrMsg())


	# close the gds files
	closefn.gds(srcobj)
	closefn.gds(destobj)

	# return
	return(invisible(NULL))
}


#######################################################################
# To merge gds files of SNP genotypes into a single gds file
#
# INPUT:
#   gds.fn -- the files name of SNP genotypes
#   out.fn -- the file name of output
#   sample.id -- a list of sample id (should be unique)
#   snpobj --
#   name.prefix --
#   snpfirstdim -- if TRUE, indicate store in individual-major order;
#   compress.annotation -- the compression method for annotations
#   compress.geno -- the compression method for genotypes
#   other.vars --
#

snpgdsCombineGeno <- function(gds.fn, out.fn,
	sample.id=NULL, snpobj=NULL, name.prefix=NULL,
	snpfirstdim=TRUE, compress.annotation="ZIP.MAX", compress.geno="",
	other.vars=NULL, verbose=TRUE)
{
	# check
	stopifnot(is.character(gds.fn))
	if (!is.null(sample.id))
	{
		stopifnot(is.list(sample.id))
		stopifnot(length(gds.fn) == length(sample.id))
	}
	if (!is.null(name.prefix))
	{
		stopifnot(is.character(name.prefix))
		stopifnot(length(gds.fn) == length(name.prefix))
	}
	stopifnot(is.null(snpobj) | inherits(snpobj, "snpgdsSNPListClass"))
	stopifnot(is.logical(snpfirstdim))

	# samples
	total.sampid <- NULL
	for (i in 1:length(gds.fn))
	{
		gdsobj <- openfn.gds(gds.fn[i])
		sample.ids <- read.gdsn(index.gdsn(gdsobj, "sample.id"))
		sampid <- sample.id[[i]]
		if (!is.null(sampid))
		{
			n.tmp <- length(sampid)
			sampid <- sample.ids %in% sampid
			n.samp <- sum(sampid)
			if (n.samp != n.tmp)
			{
				closefn.gds(gdsobj)
				stop("Some of sample.id do not exist!")
			}
			if (n.samp <= 0)
			{
				closefn.gds(gdsobj)
				stop("No sample in the working dataset.")
			}
			sample.ids <- sample.ids[sampid]
		}
		closefn.gds(gdsobj)
		# sample id
		if (is.null(name.prefix))
			total.sampid <- c(total.sampid, sample.ids)
		else
			total.sampid <- c(total.sampid, paste(name.prefix[i], sample.ids, sep="."))
	}

	# check whether total.sampid is unique
	if (length(unique(total.sampid)) != length(total.sampid))
	{
		if (!is.null(name.prefix))
			stop("Unable to make sample id unique, please specify correct `name.prefix'.")

		snpgdsCombineGeno(gds.fn=gds.fn, out.fn=out.fn,
			sample.id=sample.id, snpobj=snpobj,
			name.prefix = sprintf("p%02d", 1:length(gds.fn)),
			snpfirstdim=snpfirstdim,
			compress.annotation=compress.annotation, compress.geno=compress.geno,
			other.vars=other.vars, verbose=verbose)

		if (verbose)
			warning("Has specified the value of `name.prefix' to make sample id unique.")
		return(invisible(NULL))
	}

	# SNPs
	for (i in 1:length(gds.fn))
	{
		gdsobj <- openfn.gds(gds.fn[i])
		s <- snpgdsSNPList(gdsobj)
		if (is.null(snpobj))
		{
			# get unique snp id
			tmps <- paste(s$rs.id, s$chromosome, s$position, sep="-")
			i <- match(unique(tmps), tmps)
			snpobj <- s
			snpobj$rs.id <- snpobj$rs.id[i]
			snpobj$chromosome <- snpobj$chromosome[i]
			snpobj$position <- snpobj$position[i]
			snpobj$allele <- snpobj$allele[i]
			snpobj$afreq <- snpobj$afreq[i]
		} else
			snpobj <- snpgdsSNPListIntersect(snpobj, s)
		closefn.gds(gdsobj)
		# check
		if (length(snpobj$rs.id) <= 0)
			stop("There is no common SNP.")
	}

	# create a gds file
	gfile <- createfn.gds(out.fn)
	if (verbose)
	{
		cat("Create", out.fn, "with", length(total.sampid), "samples and",
			length(snpobj$rs.id), "SNPs\n")
	}

	add.gdsn(gfile, "sample.id", total.sampid, compress=compress.annotation, closezip=TRUE)
	if (length(snpobj$rs.id) == length(unique(snpobj$rs.id)))
	{
		add.gdsn(gfile, "snp.id", snpobj$rs.id,
			compress=compress.annotation, closezip=TRUE)
	} else {
		add.gdsn(gfile, "snp.id", as.integer(1:length(snpobj$rs.id)),
			compress=compress.annotation, closezip=TRUE)
		add.gdsn(gfile, "snp.rs.id", snpobj$rs.id, compress=compress.annotation,
			closezip=TRUE)
	}
	add.gdsn(gfile, "snp.position", snpobj$position, compress=compress.annotation,
		closezip=TRUE)
	add.gdsn(gfile, "snp.chromosome", snpobj$chromosome, compress=compress.annotation,
		closezip=TRUE)
	add.gdsn(gfile, "snp.allele", snpobj$allele, compress=compress.annotation,
		closezip=TRUE)

	# add genotype
	if (snpfirstdim)
	{
		node.geno <- add.gdsn(gfile, "genotype", valdim=c(length(snpobj$rs.id), 0),
			storage="bit2", compress=compress.geno)
		put.attr.gdsn(node.geno, "snp.order")
	} else {
		node.geno <- add.gdsn(gfile, "genotype", valdim=c(length(total.sampid), 0),
			storage="bit2", compress=compress.geno)
		put.attr.gdsn(node.geno, "sample.order")
	}

	for (i in 1:length(gds.fn))
	{
		# open the file
		gdsobj <- openfn.gds(gds.fn[i])

		# samples
		sample.ids <- read.gdsn(index.gdsn(gdsobj, "sample.id"))
		sampid <- sample.id[[i]]
		if (!is.null(sampid))
			sampid <- sample.ids %in% sampid
		# SNPs
		L <- snpgdsSNPListStrand(snpobj, snpgdsSNPList(gdsobj))

		if (verbose)
		{
			cat("\tOpen the gds file ", gds.fn[i], ".\n", sep="")
			cat("\t\t", sum(L, na.rm=TRUE), " strands of SNP loci need to be switched.\n", sep="")
		}

		# set genotype working space
		rv <- .C("gnrSetGenoSpace", as.integer(index.gdsn(gdsobj, "genotype")),
			as.logical(sampid), as.integer(!is.null(sampid)),
			as.logical(!is.na(L)), as.integer(TRUE),
			n.snp=integer(1), n.samp=integer(1),
			err=integer(1), NAOK=TRUE, PACKAGE="SNPRelate")

		# check
		if (rv$err != 0) stop(snpgdsErrMsg())
		if (rv$n.snp != length(snpobj$position))
			stop("Invalid snp alleles.")

		# write genotypes
		rv <- .C("gnrAppendGenoSpaceStrand", as.integer(node.geno), snpfirstdim,
			L[!is.na(L)], err=integer(1), NAOK=TRUE, PACKAGE="SNPRelate")
		if (rv$err != 0) stop(snpgdsErrMsg())

		# close the file
		closefn.gds(gdsobj)
	}

	# other variables
	if (!is.null(other.vars))
	{
		for (i in 1:length(other.vars))
		{
			nm <- names(other.vars)[i]
			add.gdsn(gfile, nm, val=other.vars[[i]], compress=compress.annotation)
		}
	}

	# close the gds file
	closefn.gds(gfile)
	# return
	return(invisible(NULL))
}







#######################################################################
# Plot functions
#######################################################################









#######################################################################
# To get the error message
#

snpgdsErrMsg <- function()
{
	rv <- .C("gnrErrMsg", msg=character(1), NAOK=TRUE, PACKAGE="SNPRelate")
	rv$msg
}


#######################################################################
# To get the file name of an example
#
snpgdsExampleFileName <- function()
{
	system.file("extdata", "hapmap_geno.gds", package="SNPRelate")
}







#######################################################################
# Conversion
#######################################################################
#     X    X chromosome                    -> 23
#     XY   Pseudo-autosomal region of X    -> 24
#     Y    Y chromosome                    -> 25
#     MT   Mitochondrial                   -> 26

#######################################################################
# Convert a GDS file to a PLINK ped file
#
# INPUT:
#   gdsobj -- an object of gds file
#   ped.fn -- the file name of ped format with the extended name
#   sample.id -- a vector of sample.id; if NULL, to use all samples
#   snp.id -- a vector of snp.id; if NULL, to use all SNPs
#   use.snp.rsid -- if TRUE, use "snp.rs.id" instead of "snp.id" if available
#   format -- "A/G/C/T" -- allelic codes, "A/B" -- A or B, "1/2" -- 1 or 2
#   verbose -- show information
#

snpgdsGDS2PED <- function(gdsobj, ped.fn, sample.id=NULL, snp.id=NULL,
	use.snp.rsid=TRUE, format=c("A/G/C/T", "A/B", "1/2"), verbose=TRUE)
{
	# check
	stopifnot(inherits(gdsobj, "gds.class"))
	stopifnot(is.character(ped.fn))
	format <- match.arg(format)

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
		cat("Converting from GDS to PLINK PED:\n")

	# SNPs
	total.snp.ids <- read.gdsn(index.gdsn(gdsobj, "snp.id"))
	snp.ids <- total.snp.ids
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

	# format code
	snp.idx <- match(snp.ids, total.snp.ids)
	if (format == "A/G/C/T")
	{
		n <- index.gdsn(gdsobj, "snp.allele", silent=TRUE)
		if (is.null(n))
			stop("There is no 'snp.allele' variable in the GDS file.")
		al <- read.gdsn(n)
		if (length(al) != length(total.snp.ids))
			stop("Invalid 'snp.allele' in the GDS file.")
		al <- al[snp.idx]
		fmt.code <- 1L
	} else if (format == "A/B")
	{
		al <- character(0)
		fmt.code <- 2L
	} else if (format == "1/2")
	{
		al <- character(0)
		fmt.code <- 3L
	} else
		stop("Invalid 'format'.")

	# output a MAP file
	tmp.snp.id <- snp.ids
	if (use.snp.rsid)
	{
		if (!is.null(index.gdsn(gdsobj, "snp.rs.id", silent=TRUE)))
		{
			tmp.snp.id <- read.gdsn(index.gdsn(gdsobj, "snp.rs.id"))[snp.idx]
		}
	}
	xchr <- as.character(read.gdsn(index.gdsn(gdsobj, "snp.chromosome")))[snp.idx]
	xchr[xchr=="23"] <- "X"; xchr[xchr=="25"] <- "Y"
	xchr[xchr=="24"] <- "XY"; xchr[xchr=="26"] <- "MT"
	D <- data.frame(chr = xchr, rs = tmp.snp.id,
		gen = rep(0, length(snp.idx)),
		base = read.gdsn(index.gdsn(gdsobj, "snp.position"))[snp.idx],
		stringsAsFactors = FALSE)
	write.table(D, file=paste(ped.fn, ".map", sep=""), sep="\t",
		quote=FALSE, row.names=FALSE, col.names=FALSE)
	if (verbose)
		cat("\tOutput a MAP file DONE.\n");

	# output a PED file
	if (verbose)
		cat("\tOutput a PED file ...\n");

	# set genotype working space
	node <- .C("gnrSetGenoSpace", as.integer(index.gdsn(gdsobj, "genotype")),
		as.logical(sample.id), as.logical(!is.null(sample.id)),
		as.logical(snp.id), as.logical(!is.null(snp.id)),
		n.snp=integer(1), n.samp=integer(1),
		err=integer(1), NAOK=TRUE, PACKAGE="SNPRelate")
	if (node$err != 0) stop(snpgdsErrMsg())

	# run the C code
	rv <- .C("gnrConvGDS2PED", paste(ped.fn, ".ped", sep=""),
		as.character(sample.ids), al, fmt.code, as.logical(verbose),
		err=integer(1), NAOK=TRUE, PACKAGE="SNPRelate")
	if (rv$err != 0) stop(snpgdsErrMsg())

	return(invisible(NULL))
}


#######################################################################
# Convert a GDS file to a PLINK binary ped file
#
# INPUT:
#   gdsobj -- an object of gds file
#   bed.fn -- the file name of binary ped format with the extended name
#   sample.id -- a vector of sample.id; if NULL, to use all samples
#   snp.id -- a vector of snp.id; if NULL, to use all SNPs
#   verbose -- show information
#

snpgdsGDS2BED <- function(gdsobj, bed.fn, sample.id=NULL, snp.id=NULL, verbose=TRUE)
{
	# check
	stopifnot(inherits(gdsobj, "gds.class"))
	stopifnot(is.character(bed.fn) & (length(bed.fn)==1))
	stopifnot(is.logical(verbose) & (length(verbose)==1))

	# samples
	total.samp.ids <- read.gdsn(index.gdsn(gdsobj, "sample.id"))
	sample.ids <- total.samp.ids
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
		cat("Converting from GDS to PLINK binary PED:\n")

	# SNPs
	total.snp.ids <- read.gdsn(index.gdsn(gdsobj, "snp.id"))
	snp.ids <- total.snp.ids
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

	# output a bim file
	snp.idx <- match(snp.ids, total.snp.ids)
	xchr <- as.character(read.gdsn(index.gdsn(gdsobj, "snp.chromosome")))[snp.idx]

	opt <- snpgdsOption(gdsobj)
	for (i in 1:length(opt$chromosome.code))
	{
		xchr[ xchr == opt$chromosome.code[[i]] ] <- names(opt$chromosome.code)[i]
	}
	xchr[is.na(xchr)] <- "0"

	if ((opt$autosome.start==1) & (opt$autosome.end==22))
	{
		# PLINK: Chromosome codes
		#   The autosomes should be coded 1 through 22.
		#   The following other codes can be used to specify other chromosome types:
		# X    X chromosome                    -> 23
		# Y    Y chromosome                    -> 24
		# XY   Pseudo-autosomal region of X    -> 25
		# MT   Mitochondrial                   -> 26

		xchr[xchr=="X"] <- "23"; xchr[xchr=="Y"] <- "24"; xchr[xchr=="XY"] <- "25"
		xchr[xchr=="M"] <- "26"; xchr[xchr=="MT"] <- "26"
	}

	if (!is.null(index.gdsn(gdsobj, "snp.allele", silent=TRUE)))
	{
		allele <- read.gdsn(index.gdsn(gdsobj, "snp.allele"))
		s <- unlist(strsplit(allele, "/"))
		ref <- s[seq(1, length(s), 2)]; ref <- ref[snp.idx]
		nonref <- s[seq(2, length(s), 2)]; nonref <- nonref[snp.idx]
	} else {
		warning("There is no allele information in the GDS file. ``A/B'' is used for the last two columns.")
		ref <- rep("A", length(snp.idx))
		nonref <- rep("B", length(snp.idx))
	}

	D <- data.frame(chr = xchr, rs = snp.ids,
		gen = rep(0, length(snp.idx)),
		base = read.gdsn(index.gdsn(gdsobj, "snp.position"))[snp.idx],
		A1 = ref, A2 = nonref,
		stringsAsFactors = FALSE)
	write.table(D, file=paste(bed.fn, ".bim", sep=""), sep="\t",
		quote=FALSE, row.names=FALSE, col.names=FALSE)
	if (verbose)
		cat("\tOutput a bim file DONE.\n");

	# output a fam file
	samp.idx <- match(sample.ids, total.samp.ids)
	D <- data.frame(fam = rep(0, length(samp.idx)), ind = sample.ids,
		fat = rep(0, length(samp.idx)), mot = rep(0, length(samp.idx)),
		sex = rep(0, length(samp.idx)), pheno = rep(-9, length(samp.idx)),
		stringsAsFactors = FALSE)
	write.table(D, file=paste(bed.fn, ".fam", sep=""), sep="\t",
		quote=FALSE, row.names=FALSE, col.names=FALSE)

	# output a BED file
	if (verbose)
		cat("\tOutput a BED file ...\n");

	# set genotype working space
	node <- .C("gnrSetGenoSpace", as.integer(index.gdsn(gdsobj, "genotype")),
		as.logical(sample.id), as.logical(!is.null(sample.id)),
		as.logical(snp.id), as.logical(!is.null(snp.id)),
		n.snp=integer(1), n.samp=integer(1),
		err=integer(1), NAOK=TRUE, PACKAGE="SNPRelate")
	if (node$err != 0) stop(snpgdsErrMsg())

	# Whether SNP major order or not
	snpfirstdim <- TRUE
	rd <- names(get.attr.gdsn(index.gdsn(gdsobj, "genotype")))
	if ("snp.order" %in% rd) snpfirstdim <- TRUE
	if ("sample.order" %in% rd) snpfirstdim <- FALSE

	rv <- .C("gnrConvGDS2BED", paste(bed.fn, ".bed", sep=""),
		snpfirstdim, verbose, err=integer(1),
		NAOK=TRUE, PACKAGE="SNPRelate")
	if (rv$err != 0) stop(snpgdsErrMsg())

	return(invisible(NULL))
}


#######################################################################
# Convert a GDS file to a PLINK ped file
#
# INPUT:
#   bed.fn -- binary file, genotype information
#   fam.fn -- first six columns of mydata.ped
#   bim.fn -- extended MAP file: two extra cols = allele names
#   out.gdsn -- the output gds file
#   compress.annotation -- "ZIP.max"
#   option -- allows new chromosome coding
#   verbose -- show information
#

snpgdsBED2GDS <- function(bed.fn, fam.fn, bim.fn, out.gdsfn, family=FALSE,
	compress.annotation="ZIP.max", option=NULL, verbose=TRUE)
{
	# check
	stopifnot(is.character(bed.fn) & (length(bed.fn)==1))
	stopifnot(is.character(fam.fn) & (length(fam.fn)==1))
	stopifnot(is.character(bim.fn) & (length(bim.fn)==1))
	stopifnot(is.logical(verbose))

	bed.fn <- normalizePath(bed.fn, mustWork=FALSE)
	fam.fn <- normalizePath(fam.fn, mustWork=FALSE)
	bim.fn <- normalizePath(bim.fn, mustWork=FALSE)

	if (verbose)
		cat("Start snpgdsBED2GDS ...\n")

	# detec bed.fn
	bed <- .C("gnrBEDFlag", bed.fn, snporder=integer(1), err=integer(1),
		NAOK=TRUE, PACKAGE="SNPRelate")
	if (bed$err != 0) stop(snpgdsErrMsg())
	if (verbose)
	{
		cat("\topen", bed.fn)
		if (bed$snporder == 0)
			cat(" in the individual-major mode\n")
		else
			cat(" in the SNP-major mode\n")
	}

	# read fam.fn
	famD <- read.table(fam.fn, header=FALSE, stringsAsFactors=FALSE)
	names(famD) <- c("FamilyID", "InvID", "PatID", "MatID", "Sex", "Pheno")
	if (length(unique(famD$InvID)) == dim(famD)[1])
	{
		sample.id <- famD$InvID
	} else {
		sample.id <- paste(famD$FamilyID, famD$InvID, sep="-")
		if (length(unique(sample.id)) != dim(famD)[1])
			stop("IDs in PLINK bed are not unique!")
	}
	if (verbose)
		cat("\topen", fam.fn, "DONE.\n")

	# read bim.fn
	bimD <- read.table(bim.fn, header=FALSE, stringsAsFactors=FALSE)
	names(bimD) <- c("chr", "snp.id", "map", "pos", "allele1", "allele2")

	# chromosome
	if (is.null(option)) option <- snpgdsOption()
	chrcode <- option$chromosome.code
	chr <- bimD$chr
	for (i in names(chrcode))
		chr[bimD$chr == i] <- chrcode[[i]]
	chr <- as.integer(chr)
	chr[is.na(chr)] <- 0

	# snp.id
	if (length(unique(bimD$snp.id)) == dim(bimD)[1])
	{
		snp.id <- bimD$snp.id; snp.rs.id <- NULL
	} else {
		snp.id <- 1:dim(bimD)[1]; snp.rs.id <- bimD$snp.id
	}

	if (verbose)
		cat("\topen", bim.fn, "DONE.\n")

	# create GDS file
	gfile <- createfn.gds(out.gdsfn)

	# add "sample.id"
	add.gdsn(gfile, "sample.id", sample.id, compress=compress.annotation, closezip=TRUE)
	# add "snp.id"
	add.gdsn(gfile, "snp.id", snp.id, compress=compress.annotation, closezip=TRUE)
	# add "snp.rs.id"
	if (!is.null(snp.rs.id))
		add.gdsn(gfile, "snp.rs.id", snp.rs.id, compress=compress.annotation, closezip=TRUE)
	# add "snp.position"
	add.gdsn(gfile, "snp.position", bimD$pos, compress=compress.annotation, closezip=TRUE)
	# add "snp.chromosome"
	v.chr <- add.gdsn(gfile, "snp.chromosome", chr, storage="uint8",
		compress=compress.annotation, closezip=TRUE)
	# add "snp.allele"
	add.gdsn(gfile, "snp.allele", paste(bimD$allele1, bimD$allele2, sep="/"),
		compress=compress.annotation, closezip=TRUE)

	# snp.chromosome
	put.attr.gdsn(v.chr, "autosome.start", option$autosome.start)
	put.attr.gdsn(v.chr, "autosome.end", option$autosome.end)
	for (i in 1:length(option$chromosome.code))
	{
		put.attr.gdsn(v.chr, names(option$chromosome.code)[i],
			option$chromosome.code[[i]])
	}

	# sync file
	sync.gds(gfile)

	nSamp <- dim(famD)[1]; nSNP <- dim(bimD)[1]
	if (verbose)
	{
		cat(date(), "\tstore sample id, snp id, position, and chromosome.\n")
		cat(sprintf("\tstart writing: %d samples, %d SNPs ...\n", nSamp, nSNP))
	}

	# add "gonetype", 2 bits to store one genotype
	if (bed$snporder == 0)
	{
		gGeno <- add.gdsn(gfile, "genotype", storage="bit2", valdim=c(nSNP, nSamp))
		put.attr.gdsn(gGeno, "snp.order")
	} else {
		gGeno <- add.gdsn(gfile, "genotype", storage="bit2", valdim=c(nSamp, nSNP))
		put.attr.gdsn(gGeno, "sample.order")
	}

	rv <- .C("gnrConvBED2GDS", bed.fn, as.integer(gGeno), verbose, err=integer(1),
		NAOK=TRUE, PACKAGE="SNPRelate")
	if (rv$err != 0) stop(snpgdsErrMsg())

	# sync file
	sync.gds(gfile)

	# add "sample.annot"
	sex <- rep("", length(sample.id))
	sex[famD$Sex==1] <- "M"; sex[famD$Sex==2] <- "F"

	if (family)
	{
		samp.annot <- data.frame(father=famD$PatID, mother=famD$MatID,
			sex=sex, phenotype=famD$Pheno, stringsAsFactors=FALSE)
	} else {
		samp.annot <- data.frame(sex=sex, phenotype=famD$Pheno, stringsAsFactors=FALSE)
	}
	add.gdsn(gfile, "sample.annot", samp.annot, compress=compress.annotation, closezip=TRUE)

	# close files
	closefn.gds(gfile)

	if (verbose) cat(date(), "\tDone.\n")

	return(invisible(NULL))
}


#######################################################################
# To convert a genotype GDS file to Eigenstrat format
#
# INPUT:
#   gdsobj -- an object of gds file
#   eigen.fn -- the file name of Eigenstrat format with the extended name
#   sample.id -- a vector of sample.id; if NULL, to use all samples
#   snp.id -- a vector of snp.id; if NULL, to use all SNPs
#   verbose -- show progress
#
# OUTPUT:
#    *.eigenstratgeno -- "genotype file in EIGENSTRAT format"
#    *.snp -- "snp file"
#    *.ind -- "individual file"
#

snpgdsGDS2Eigen <- function(gdsobj, eigen.fn, sample.id=NULL, snp.id=NULL, verbose=TRUE)
{
	# check
	stopifnot(inherits(gdsobj, "gds.class"))
	stopifnot(is.character(eigen.fn))

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
	} else
		sample.id <- rep(TRUE, length(sample.ids))

	if (verbose)
		cat("Converting from GDS to EIGENSOFT:\n")

	# SNPs
	total.snp.ids <- read.gdsn(index.gdsn(gdsobj, "snp.id"))
	snp.ids <- total.snp.ids
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
	} else
		snp.id <- rep(TRUE, length(snp.ids))

	# making the "*.snp" file ...
	tmpD <- data.frame(
		snpid = read.gdsn(index.gdsn(gdsobj, "snp.id"))[snp.id],
		chrom = read.gdsn(index.gdsn(gdsobj, "snp.chromosome"))[snp.id],
		map = rep(0.0, sum(snp.id)),
		pos = read.gdsn(index.gdsn(gdsobj, "snp.position"))[snp.id],
		stringsAsFactors = FALSE
	)
	write.table(tmpD, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE,
		file = paste(eigen.fn, ".snp", sep=""))
	if (verbose)
		cat("\tsave to *.snp:", dim(tmpD)[1], "snps\n")

	# making the "*.ind" file ...
	sex <- try(read.gdsn(index.gdsn(gdsobj, "sample.annot/sex")), TRUE)
	if (class(sex) == "try-error")
	{
		sex <- rep("U", sum(sample.id))
	} else {
		sex <- as.character(sex)[sample.id]
		if (!all(sex %in% c("F", "M"), na.rm=TRUE))
			stop("The gender variable in GDS file should be either \"M\" or \"F\".")
		sex[is.na(sex)] <- "U"
	}
	tmpD <- data.frame(
		sampid = read.gdsn(index.gdsn(gdsobj, "sample.id"))[sample.id],
		gender = sex, label = rep("control", sum(sample.id)),
		stringsAsFactors = FALSE
	)
	write.table(tmpD, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE,
		file = paste(eigen.fn, ".ind", sep=""))
	if (verbose)
		cat("\tsave to *.ind:", dim(tmpD)[1], "samples\n")

	# making the "*.eigenstratgeno" file ...

	# set genotype working space
	node <- .C("gnrSetGenoSpace", as.integer(index.gdsn(gdsobj, "genotype")),
		as.logical(sample.id), as.logical(!is.null(sample.id)),
		as.logical(snp.id), as.logical(!is.null(snp.id)),
		n.snp=integer(1), n.samp=integer(1),
		err=integer(1), NAOK=TRUE, PACKAGE="SNPRelate")
	if (node$err != 0) stop(snpgdsErrMsg())

	rv <- .C("gnrConvGDS2EIGEN", paste(eigen.fn, ".eigenstratgeno", sep=""),
		as.logical(verbose), err=integer(1), NAOK=TRUE, PACKAGE="SNPRelate")
	if (rv$err != 0) stop(snpgdsErrMsg())

	if (verbose) cat("Done.\n")

	return(invisible(NULL))
}


#######################################################################
# Convert a VCF (sequence) file to a GDS file (extract SNP data)
#
# INPUT:
#   vcf.fn -- the file name of VCF format
#   outfn.gds -- the output gds file
#   nblock -- the number of lines in buffer
#   method -- biallelic SNPs, or copy number of variants
#   compress.annotation -- the compression method for sample and snp annotations
#   verbose -- show information
#

snpgdsVCF2GDS <- function(vcf.fn, outfn.gds, nblock=1024,
	method = c("biallelic.only", "copy.num.of.ref"),
	compress.annotation="ZIP.max", snpfirstdim=FALSE, option = NULL,
	verbose=TRUE)
{
	# check
	stopifnot(is.character(vcf.fn))
	stopifnot(is.character(outfn.gds))
	stopifnot(is.logical(snpfirstdim) & (length(snpfirstdim)==1))

	method <- match.arg(method)
	if (is.null(option)) option <- snpgdsOption()



	######################################################################
	# Scan VCF file -- get sample id

	scan.vcf.sampid <- function(fn)
	{
		# open the vcf file
		opfile <- file(fn, open="r")

		# read header
		fmtstr <- substring(readLines(opfile, n=1), 3)
		samp.id <- NULL
		while (length(s <- readLines(opfile, n=1)) > 0)
		{
			if (substr(s, 1, 6) == "#CHROM")
			{
				samp.id <- scan(text=s, what=character(0), sep="\t", quiet=TRUE)[-c(1:9)]
				break
			}
		}
		if (is.null(samp.id))
		{
			close(opfile)
			stop("Error VCF format: invalid sample id!")
		}

		# close the file
		close(opfile)

		return(samp.id)
	}


	######################################################################
	# Scan VCF file -- get marker information

	scan.vcf.marker <- function(fn, method)
	{
		if (verbose)
			cat(sprintf("\tfile: %s\n", fn))

		# total number of rows and columns
		Cnt <- count.fields(fn, sep="\t")
		# check
		if (any(Cnt != Cnt[1]))
			stop(sprintf("The file (%s) has different numbers of columns.", fn))

		line.cnt <- length(Cnt)
		col.cnt <- max(Cnt)
		if (verbose)
			cat(sprintf("\tcontent: %d rows x %d columns\n", line.cnt, col.cnt))

		# open the vcf file
		opfile <- file(fn, open="r")

		# read header
		fmtstr <- substring(readLines(opfile, n=1), 3)
		while (length(s <- readLines(opfile, n=1)) > 0)
		{
			if (substr(s, 1, 6) == "#CHROM")
				break
		}

		# init ...
		chr <- character(line.cnt); position <- integer(line.cnt)
		snpidx <- integer(line.cnt); snp.rs <- character(line.cnt)
		snp.allele <- character(line.cnt)
		snp.cnt <- 0; var.cnt <- 0

		if (method == "biallelic.only")
		{
			while (length(s <- readLines(opfile, n=nblock)) > 0)
			{
				for (i in 1:length(s))
				{
					var.cnt <- var.cnt + 1
					ss <- scan(text=s[i], what=character(0), sep="\t", quiet=TRUE, n=5)
					if (all(ss[c(4,5)] %in% c("A", "G", "C", "T", "a", "g", "c", "t")))
					{
						snp.cnt <- snp.cnt + 1
						chr[snp.cnt] <- ss[1]
						position[snp.cnt] <- as.integer(ss[2])
						snpidx[snp.cnt] <- var.cnt
						snp.rs[snp.cnt] <- ss[3]
						snp.allele[snp.cnt] <- paste(ss[4], ss[5], sep="/")
					}
				}
			}
		} else {
			while (length(s <- readLines(opfile, n=nblock)) > 0)
			{
				for (i in 1:length(s))
				{
					var.cnt <- var.cnt + 1
					ss <- scan(text=s[i], what=character(0), sep="\t", quiet=TRUE, n=5)
					snp.cnt <- snp.cnt + 1
					chr[snp.cnt] <- ss[1]
					position[snp.cnt] <- as.integer(ss[2])
					snpidx[snp.cnt] <- var.cnt
					snp.rs[snp.cnt] <- ss[3]
					snp.allele[snp.cnt] <- paste(ss[4], ss[5], sep="/")
				}
			}
		}

		# close the file
		close(opfile)


		# chromosomes
		chr <- chr[1:snp.cnt]
		flag <- match(chr, names(option$chromosome.code))
		chr[!is.na(flag)] <- unlist(option$chromosome.code)[ flag[!is.na(flag)] ]
		chr <- suppressWarnings(as.integer(chr))
		chr[is.na(chr)] <- -1

		snp.allele <- gsub(".", "/", snp.allele[1:snp.cnt], fixed=TRUE)
		list(chr = chr, position = position[1:snp.cnt],
			snpidx = snpidx[1:snp.cnt], snp.rs = snp.rs[1:snp.cnt],
			snp.allele = snp.allele
		)
	}


	######################################################################
	# Scan VCF file -- get marker information

	scan.vcf.geno <- function(fn, gGeno, method, start)
	{
		# matching codes
		geno.str <- c("0|0", "0|1", "1|0", "1|1", "0/0", "0/1", "1/0", "1/1",
			"0", "1",
			"0|0|0", "0|0|1", "0|1|0", "0|1|1", "1|0|0", "1|0|1", "1|1|0", "1|1|1",
			"0/0/0", "0/0/1", "0/1/0", "0/1/1", "1/0/0", "1/0/1", "1/1/0", "1/1/1")
		geno.code <- as.integer(c(2, 1, 1, 0, 2, 1, 1, 0,
			1, 0,
			2, 1, 1, 1, 1, 1, 1, 0, 
			2, 1, 1, 1, 1, 1, 1, 0))

		# open the vcf file
		opfile <- file(fn, open="r")
		# read header
		fmtstr <- substring(readLines(opfile, n=1), 3)
		while (length(s <- readLines(opfile, n=1)) > 0)
		{
			if (substr(s, 1, 6) == "#CHROM")
				break
		}

		# scan
		snp.cnt <- start

		if (method == "biallelic.only")
		{
			while (length(s <- readLines(opfile, n=nblock)) > 0)
			{
				gx <- NULL
				for (i in 1:length(s))
				{
					ss <- scan(text=s[i], what=character(0), sep="\t", quiet=TRUE, n=5)
					if (all(ss[c(4,5)] %in% c("A", "G", "C", "T", "a", "g", "c", "t")))
					{
						ss <- scan(text=s[i], what=character(0), sep="\t", quiet=TRUE)[-c(1:9)]
						ss <- sapply(strsplit(ss, ":"), FUN = function(x) x[1])
						x <- match(ss, geno.str)
						x <- geno.code[x]
						x[is.na(x)] <- as.integer(3)
						gx <- cbind(gx, x)
					}
				}
				if (!is.null(gx))
				{
					if (snpfirstdim)
						write.gdsn(gGeno, t(gx), start=c(snp.cnt,1), count=c(ncol(gx),-1))
					else {
						print(snp.cnt)
						write.gdsn(gGeno, gx, start=c(1,snp.cnt), count=c(-1,ncol(gx)))
					}
					snp.cnt <- snp.cnt + ncol(gx)
				}
			}
		} else {
			while (length(s <- readLines(opfile, n=nblock)) > 0)
			{
				gx <- NULL
				for (i in 1:length(s))
				{
					ss <- scan(text=s[i], what=character(0), sep="\t", quiet=TRUE)[-c(1:9)]
					x <- sapply(strsplit(ss, ":"), FUN = function(x) {
						a <- unlist(strsplit(x[1], ""))
						if (any(a == "."))
							NA
						else
							sum(a == "0")
					})
					x[x > 2] <- 2
					x[is.na(x)] <- as.integer(3)
					gx <- cbind(gx, x)
				}
				if (!is.null(gx))
				{
					if (snpfirstdim)
						write.gdsn(gGeno, t(gx), start=c(snp.cnt,1), count=c(ncol(gx),-1))
					else
						write.gdsn(gGeno, gx, start=c(1,snp.cnt), count=c(-1,ncol(gx)))
					snp.cnt <- snp.cnt + ncol(gx)
				}
			}
		}

		# close the file
		close(opfile)
		
		snp.cnt - start
	}



	######################################################################
	######################################################################
	# Starting ...
	######################################################################
	######################################################################

	if (verbose)
	{
		cat("Start snpgdsVCF2GDS ...\n")
		if (method == "biallelic.only")
			cat("\tExtracting bi-allelic and polymorhpic SNPs.\n")
		else
			cat("\tStoring dosage of the reference allele for all variant sites, including bi-allelic SNPs, multi-allelic SNPs, indels and structural variants.\n")
		cat("\tScanning ...\n")
	}


	####################################
	# sample.id

	sample.id <- NULL
	for (fn in vcf.fn)
	{
		s <- scan.vcf.sampid(fn)
		if (!is.null(sample.id))
		{
			if (length(sample.id) != length(s))
				stop("All VCF files should have the same sample id.")
			if (any(sample.id != s))
				stop("All VCF files should have the same sample id.")
		} else
			sample.id <- s
	}


	####################################
	# genetic markers

	all.chr <- integer()
	all.position <- integer()
	all.snpidx <- integer()
	all.snp.rs <- character()
	all.snp.allele <- character()

	for (fn in vcf.fn)
	{
		v <- scan.vcf.marker(fn, method)
		
		all.chr <- c(all.chr, v$chr)
		all.position <- c(all.position, v$position)
		all.snpidx <- c(all.snpidx, length(all.snpidx) + v$snpidx)
		all.snp.rs <- c(all.snp.rs, v$snp.rs)
		all.snp.allele <- c(all.snp.allele, v$snp.allele)
	}


	####################################
	# genetic variants

	nSamp <- length(sample.id)
	nSNP <- length(all.chr)
	if (verbose)
	{
		cat(date(), "\tstore sample id, snp id, position, and chromosome.\n")
		cat(sprintf("\tstart writing: %d samples, %d SNPs ...\n", nSamp, nSNP))
	}


	######################################################################
	# create GDS file
	#
	gfile <- createfn.gds(outfn.gds)

	# add "sample.id"
	add.gdsn(gfile, "sample.id", sample.id, compress=compress.annotation, closezip=TRUE)
	# add "snp.id"
	add.gdsn(gfile, "snp.id", as.integer(all.snpidx), compress=compress.annotation, closezip=TRUE)
	# add "snp.rs.id"
	add.gdsn(gfile, "snp.rs.id", all.snp.rs, compress=compress.annotation, closezip=TRUE)
	# add "snp.position"
	add.gdsn(gfile, "snp.position", all.position, compress=compress.annotation, closezip=TRUE)
	# add "snp.chromosome"
	v.chr <- add.gdsn(gfile, "snp.chromosome", all.chr, storage="int32", compress=compress.annotation, closezip=TRUE)
	# add "snp.allele"
	add.gdsn(gfile, "snp.allele", all.snp.allele, compress=compress.annotation, closezip=TRUE)

	# snp.chromosome
	put.attr.gdsn(v.chr, "autosome.start", option$autosome.start)
	put.attr.gdsn(v.chr, "autosome.end", option$autosome.end)
	for (i in 1:length(option$chromosome.code))
	{
		put.attr.gdsn(v.chr, names(option$chromosome.code)[i],
			option$chromosome.code[[i]])
	}

	# sync file
	sync.gds(gfile)

	# add "gonetype", 2 bits to store one genotype
	if (snpfirstdim)
	{
		gGeno <- add.gdsn(gfile, "genotype", storage="bit2", valdim=c(nSNP, nSamp))
		put.attr.gdsn(gGeno, "snp.order")
	} else {
		gGeno <- add.gdsn(gfile, "genotype", storage="bit2", valdim=c(nSamp, nSNP))
		put.attr.gdsn(gGeno, "sample.order")
	}
	# sync file
	sync.gds(gfile)


	####################################
	# genetic genotypes

	snp.start <- 1
	for (fn in vcf.fn)
	{
		if (verbose)
			cat(sprintf("\tfile: %s\n", fn))
		s <- scan.vcf.geno(fn, gGeno, method, start=snp.start)
		snp.start <- snp.start + s		
		sync.gds(gfile)  # sync file
	}


	# close files
	closefn.gds(gfile)

	if (verbose) cat(date(), "\tDone.\n")

	return(invisible(NULL))
}




#######################################################################
# SNPRelate Option
#

snpgdsOption <- function(gdsobj=NULL, autosome.start=1, autosome.end=22, ...)
{
	ans <- list(
		autosome.start = as.integer(autosome.start),  # the starting index of autosome
		autosome.end   = as.integer(autosome.end)     # the ending idex of autosome
	)

	if (!is.null(gdsobj))
	{
		# ignore the arguments "..."

		# chromosome
		n <- index.gdsn(gdsobj, "snp.chromosome")
		lst <- get.attr.gdsn(n)

		if (!is.null(lst$autosome.start))
		{
			ans$autosome.start <- lst$autosome.start
			lst <- lst[-match("autosome.start", names(lst))]
		}
		if (!is.null(lst$autosome.end))
		{
			ans$autosome.end <- lst$autosome.end
			lst <- lst[-match("autosome.end", names(lst))]
		}

		ns <- names(lst)
		if (!("X" %in% ns))       # X chromosome
			lst$X <- as.integer(ans$autosome.end + 1)
		if (!("XY" %in% ns))      # Pseudo-autosomal region of X
			lst$XY <- as.integer(ans$autosome.end + 2)
		if (!("Y" %in% ns))       # Y chromosome
			lst$Y <- as.integer(ans$autosome.end + 3)
		if (!("M" %in% ns))       # Mitochondrial
			lst$M <- as.integer(ans$autosome.end + 4)
		if (!("MT" %in% ns))      # Mitochondrial
			lst$MT = as.integer(ans$autosome.end + 4)

		ans$chromosome.code <- lst

		# SNP genotype
		ans$file$filename <- gdsobj$filename

		n <- index.gdsn(gdsobj, "genotype")
		lst <- get.attr.gdsn(n)
		if ("sample.order" %in% names(lst))
			ans$file$geno.dim <- "sample-by-snp"
		else
			ans$file$geno.dim <- "snp-by-sample"

	} else {
		# incorporate the arguments "..."

		lst <- list(...)
		lst <- lst[names(lst) != ""]

		if (length(lst) <= 0)
		{
			ans$chromosome.code = list(
				X  = as.integer(autosome.end + 1),    # X chromosome
				XY = as.integer(autosome.end + 2),    # Pseudo-autosomal region of X
				Y  = as.integer(autosome.end + 3),    # Y chromosome
				M  = as.integer(autosome.end + 4),    # Mitochondrial
				MT = as.integer(autosome.end + 4)     # Mitochondrial
			)
		} else {
			ans$chromosome.code <- lst
		}
	}

	ans
}






#######################################################################
# Internal R library functions
#######################################################################

.onAttach <- function(lib, pkg)
{
	# get the filename of the dynamic-link library
	lib.fn <- as.character(getLoadedDLLs()$gdsfmt[[2]])
	if (length(lib.fn) != 1)
		stop("The package 'gdsfmt' was not installed correctly.")

	# init SNPRelate
	rv <- .C("gnrInit", lib.fn, err=character(1), sse=integer(1),
		PACKAGE="SNPRelate")
	if (rv$err != "") stop(rv$err)

	# information
	packageStartupMessage("SNPRelate: 0.9.19")
	if (rv$sse != 0)
		packageStartupMessage("Supported by Streaming SIMD Extensions 2 (SSE2).")

	TRUE
}

.Last.lib <- function(libpath)
{
	# finalize SNPRelate
	rv <- .C("gnrDone", PACKAGE="SNPRelate")
	TRUE
}
