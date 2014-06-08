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
# Identity-by-Descent (IBD) analysis
#######################################################################

#######################################################################
# To calculate the identity-by-descent (IBD) matrix (PLINK Moment) for SNP genotypes
#
# INPUT:
#   gdsobj -- an object of gds file
#   sample.id -- a vector of sample.id; if NULL, to use all samples
#   snp.id -- a vector of snp.id; if NULL, to use all SNPs
#   autosome.only -- whether only use autosomal SNPs
#   remove.monosnp -- whether remove monomorphic snps or not
#   maf -- the threshold of minor allele frequencies, keeping ">= maf"
#   missing.rate -- the threshold of missing rates, keeping "<= missing.rate"#
#   allele.freq -- specify the allele frequencies
#   kinship -- if TRUE, output estimated kinship coefficients
#   kinship.constraint -- constrict IBD coeff in the geneloical region
#   num.thread -- the number of threads to be used
#   verbose -- show information
#

snpgdsIBDMoM <- function(gdsobj, sample.id=NULL, snp.id=NULL, autosome.only=TRUE,
	remove.monosnp=TRUE, maf=NaN, missing.rate=NaN, allele.freq=NULL,
	kinship=FALSE, kinship.constraint=FALSE, num.thread=1, verbose=TRUE)
{
	# check
	stopifnot(inherits(gdsobj, "gds.class"))
	stopifnot(is.numeric(num.thread) & (num.thread>0))
	stopifnot(is.logical(verbose))
	stopifnot(is.null(allele.freq) | is.numeric(allele.freq))

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
		cat("Identity-By-Descent analysis (PLINK method of moment) on SNP genotypes:\n");

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
		if (!is.null(allele.freq))
		{
			if (n.snp != length(allele.freq))
				stop("`length(allele.freq)' should equal to the number of SNPs.")
		}
		if (autosome.only)
		{
			chr <- read.gdsn(index.gdsn(gdsobj, "snp.chromosome"))
			opt <- snpgdsOption(gdsobj)
			chridx <- chr[snp.id]
			snp.id <- snp.id & (chr %in% c(opt$autosome.start : opt$autosome.end))
			if (!is.null(allele.freq))
				allele.freq <- allele.freq[chridx %in% 1:22]
			if (verbose)
			{
				tmp <- n.snp - sum(snp.id)
				if (tmp > 0) cat("Removing", tmp, "non-autosomal SNPs\n")
			}
		}
		snp.ids <- snp.ids[snp.id]
	} else {
		if (!is.null(allele.freq))
		{
			if (length(snp.ids) != length(allele.freq))
				stop("`length(allele.freq)' should equal to the number of SNPs.")
		}
		if (autosome.only)
		{
			chr <- read.gdsn(index.gdsn(gdsobj, "snp.chromosome"))
			opt <- snpgdsOption(gdsobj)
			snp.id <- chr %in% c(opt$autosome.start : opt$autosome.end)
			snp.ids <- snp.ids[snp.id]
			if (!is.null(allele.freq))
				allele.freq <- allele.freq[snp.id]
			if (verbose)
				cat("Removing", length(chr) - length(snp.ids), "non-autosomal SNPs\n")
		}
	}

	# check
	if (!is.null(allele.freq))
	{
		cat(sprintf("Specifying allele frequencies, mean: %0.3f, sd: %0.3f\n",
			mean(allele.freq, na.rm=TRUE), sd(allele.freq, na.rm=TRUE)))
		cat("*** A correction factor based on allele counts is not used, since the allele frequencies are specified.\n")
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
		if (is.null(allele.freq))
		{
			rv <- .C("gnrSelSNP_Base", as.logical(remove.monosnp),
				as.double(maf), as.double(missing.rate),
				out.num=integer(1), out.snpflag = logical(node$n.snp),
				err=integer(1), NAOK=TRUE, PACKAGE="SNPRelate")
		} else {
			rv <- .C("gnrSelSNP_Base_Ex", as.double(allele.freq), as.logical(remove.monosnp),
				as.double(maf), as.double(missing.rate),
				out.num=integer(1), out.snpflag = logical(node$n.snp),
				err=integer(1), NAOK=TRUE, PACKAGE="SNPRelate")
		}
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
		if (num.thread <= 1)
			cat("\tUsing", num.thread, "CPU core.\n")
		else
			cat("\tUsing", num.thread, "CPU cores.\n")
	}

	# call the C function
	rv <- .Call("gnrIBD_PLINK", verbose, TRUE, as.integer(num.thread),
		as.double(allele.freq), !is.null(allele.freq),
		as.logical(kinship.constraint), PACKAGE="SNPRelate")
	names(rv) <- c("k0", "k1", "afreq")

	# return
	rv <- list(sample.id = sample.ids, snp.id = snp.ids,
		afreq = rv$afreq, k0 = rv$k0, k1 = rv$k1)
	if (kinship)
		rv$kinship <- 0.5*(1 - rv$k0 - rv$k1) + 0.25*rv$k1
	rv$afreq[rv$afreq < 0] <- NaN
	class(rv) <- "snpgdsIBDClass"
	return(rv)
}



#######################################################################
# To calculate the identity-by-descent (IBD) matrix (MLE) for SNP genotypes
#
# INPUT:
#   gdsobj -- an object of gds file
#   sample.id -- a vector of sample.id; if NULL, to use all samples
#   snp.id -- a vector of snp.id; if NULL, to use all SNPs
#   autosome.only -- whether only use autosomal SNPs
#   remove.monosnp -- whether remove monomorphic snps or not
#   maf -- the threshold of minor allele frequencies, keeping ">= maf"
#   missing.rate -- the threshold of missing rates, keeping "<= missing.rate"
#   kinship -- if TRUE, output estimated kinship coefficients
#   kinship.constraint --
#   out.num.iter = FALSE
#   num.thread -- the number of threads to be used
#   verbose -- show information
#

snpgdsIBDMLE <- function(gdsobj, sample.id=NULL, snp.id=NULL, autosome.only=TRUE,
	remove.monosnp=TRUE, maf=NaN, missing.rate=NaN, kinship=FALSE,
	kinship.constraint=FALSE, allele.freq=NULL, method=c("EM", "downhill.simplex"),
	max.niter=1000, reltol=sqrt(.Machine$double.eps), coeff.correct=TRUE,
	out.num.iter=TRUE, num.thread=1, verbose=TRUE)
{
	# check
	stopifnot(inherits(gdsobj, "gds.class"))
	stopifnot(is.numeric(num.thread) & (num.thread>0))
	stopifnot(is.logical(out.num.iter))
	stopifnot(is.logical(verbose))
	stopifnot(method %in% c("EM", "downhill.simplex"))
	stopifnot(is.null(allele.freq) | is.numeric(allele.freq))

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
		cat("Identity-By-Descent analysis (MLE) on SNP genotypes:\n");

	# SNPs
	snp.ids <- read.gdsn(index.gdsn(gdsobj, "snp.id"))
	if (!is.null(snp.id))
	{
		n.tmp <- length(snp.id)
		snp.id <- snp.ids %in% snp.id
		n.snp <- sum(snp.id)
		if (n.snp != n.tmp)
			stop("Some of `snp.id' do not exist!")
		if (n.snp <= 0)
			stop("No SNP in the working dataset.")
		if (!is.null(allele.freq))
		{
			if (n.snp != length(allele.freq))
				stop("`length(allele.freq)' should equal to the number of SNPs.")
		}
		if (autosome.only)
		{
			chr <- read.gdsn(index.gdsn(gdsobj, "snp.chromosome"))
			opt <- snpgdsOption(gdsobj)
			chridx <- chr[snp.id]
			snp.id <- snp.id & (chr %in% c(opt$autosome.start : opt$autosome.end))
			if (!is.null(allele.freq))
				allele.freq <- allele.freq[chridx %in% 1:22]
			if (verbose)
			{
				tmp <- n.snp - sum(snp.id)
				if (tmp > 0) cat("Removing", tmp, "non-autosomal SNPs\n")
			}
		}
		snp.ids <- snp.ids[snp.id]
	} else {
		if (!is.null(allele.freq))
		{
			if (length(snp.ids) != length(allele.freq))
				stop("`length(allele.freq)' should equal to the number of SNPs.")
		}
		if (autosome.only)
		{
			chr <- read.gdsn(index.gdsn(gdsobj, "snp.chromosome"))
			opt <- snpgdsOption(gdsobj)
			snp.id <- chr %in% c(opt$autosome.start : opt$autosome.end)
			snp.ids <- snp.ids[snp.id]
			if (!is.null(allele.freq))
				allele.freq <- allele.freq[snp.id]
			if (verbose)
				cat("Removing", length(chr) - length(snp.ids), "non-autosomal SNPs\n")
		}
	}

	# check
	if (!is.null(allele.freq))
	{
		cat(sprintf("Specifying allele frequencies, mean: %0.3f, sd: %0.3f\n",
			mean(allele.freq, na.rm=TRUE), sd(allele.freq, na.rm=TRUE)))
	}

	# method
	if (method[1] == "EM")
		method <- 0
	else if (method[1] == "downhill.simplex")
		method <- 1
	else
		stop("Invalid MLE method!")

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
		if (is.null(allele.freq))
		{
			rv <- .C("gnrSelSNP_Base", as.logical(remove.monosnp),
				as.double(maf), as.double(missing.rate),
				out.num=integer(1), out.snpflag = logical(node$n.snp),
				err=integer(1), NAOK=TRUE, PACKAGE="SNPRelate")
		} else {
			rv <- .C("gnrSelSNP_Base_Ex", as.double(allele.freq), as.logical(remove.monosnp),
				as.double(maf), as.double(missing.rate),
				out.num=integer(1), out.snpflag = logical(node$n.snp),
				err=integer(1), NAOK=TRUE, PACKAGE="SNPRelate")
		}
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
		if (num.thread <= 1)
			cat("\tUsing", num.thread, "CPU core.\n")
		else
			cat("\tUsing", num.thread, "CPU cores.\n")
	}

	# get internal info
	sz <- .C("gnrIBD_SizeInt", size=integer(1), nsnp4 = integer(1),
		NAOK=TRUE, PACKAGE="SNPRelate")

	# call parallel IBD
	rv <- .C("gnrIBD_MLE", as.double(allele.freq), !is.null(allele.freq),
		as.logical(kinship.constraint),
		as.integer(max.niter), as.double(reltol), as.logical(coeff.correct),
		as.integer(method), verbose, TRUE,
		as.integer(num.thread), out.num.iter,
		integer(sz$size), double(sz$nsnp4),
		k0 = matrix(NaN, ncol=node$n.samp, nrow=node$n.samp),
		k1 = matrix(NaN, ncol=node$n.samp, nrow=node$n.samp),
		afreq = double(node$n.snp),
		niter = switch(out.num.iter+1, integer(0),
			matrix(as.integer(NA), ncol=node$n.samp, nrow=node$n.samp)),
		err = integer(1), NAOK=TRUE, PACKAGE="SNPRelate")
	if (rv$err != 0) stop(snpgdsErrMsg())

	# return
	rv <- list(sample.id=sample.ids, snp.id=snp.ids, afreq=rv$afreq,
		k0=rv$k0, k1=rv$k1, niter=rv$niter)
	if (kinship)
		rv$kinship <- 0.5*(1 - rv$k0 - rv$k1) + 0.25*rv$k1
	rv$afreq[rv$afreq < 0] <- NaN
	class(rv) <- "snpgdsIBDClass"
	return(rv)
}



#######################################################################
# To calculate the identity-by-descent (IBD) matrix (MLE) for SNP genotypes
#
# INPUT:
#   gdsobj -- an object of gds file
#   ibdobj -- an object of snpgdsIBDClass
#

snpgdsIBDMLELogLik <- function(gdsobj, ibdobj, k0=NaN, k1=NaN,
	relatedness=c("", "self", "fullsib", "offspring", "halfsib", "cousin", "unrelated"))
{
	# check
	stopifnot(inherits(gdsobj, "gds.class"))
	stopifnot(inherits(ibdobj, "snpgdsIBDClass"))
	stopifnot(is.numeric(k0) & length(k0)==1)
	stopifnot(is.numeric(k1) & length(k1)==1)
	stopifnot(relatedness %in% c("", "self", "fullsib", "offspring", "halfsib",
		"cousin", "unrelated"))

	# samples
	sample.ids <- read.gdsn(index.gdsn(gdsobj, "sample.id"))
	sample.id <- sample.ids %in% ibdobj$sample.id
	# SNPs
	snp.ids <- read.gdsn(index.gdsn(gdsobj, "snp.id"))
	snp.id <- snp.ids %in% ibdobj$snp.id

	# call C codes
	# set genotype working space
	node <- .C("gnrSetGenoSpace", as.integer(index.gdsn(gdsobj, "genotype")),
		as.integer(sample.id), as.integer(!is.null(sample.id)),
		as.integer(snp.id), as.integer(!is.null(snp.id)),
		n.snp=integer(1), n.samp=integer(1),
		err=integer(1), NAOK=TRUE, PACKAGE="SNPRelate")
	if (node$err != 0) stop(snpgdsErrMsg())

	# relatedness
	relatedness <- relatedness[1]
	if (relatedness == "self")
	{
		k0 <- 0; k1 <- 0
	} else if (relatedness == "fullsib")
	{
		k0 <- 0.25; k1 <- 0.5
	} else if (relatedness == "offspring")
	{
		k0 <- 0; k1 <- 1
	} else if (relatedness == "halfsib")
	{
		k0 <- 0.5; k1 <- 0.5
	} else if (relatedness == "cousin")
	{
		k0 <- 0.75; k1 <- 0.25
	} else if (relatedness == "unrelated")
	{
		k0 <- 1; k1 <- 0
	}

	# call log likelihood
	sz <- .C("gnrIBD_SizeInt", size=integer(1), nsnp4 = integer(1),
		NAOK=TRUE, PACKAGE="SNPRelate")
	if (is.finite(k0) & is.finite(k1))
	{
		rv <- .C("gnrIBD_LogLik_k01", ibdobj$afreq, as.double(k0), as.double(k1),
			integer(sz$size), double(sz$nsnp4),
			loglik = matrix(NaN, ncol=node$n.samp, nrow=node$n.samp),
			err=integer(1), NAOK=TRUE, PACKAGE="SNPRelate")
	} else {
		rv <- .C("gnrIBD_LogLik", ibdobj$afreq, ibdobj$k0, ibdobj$k1,
			integer(sz$size), double(sz$nsnp4),
			loglik = matrix(NaN, ncol=node$n.samp, nrow=node$n.samp),
			err=integer(1), NAOK=TRUE, PACKAGE="SNPRelate")
	}
	if (rv$err != 0) stop(snpgdsErrMsg())

	# return
	return(rv$loglik)
}



#######################################################################
# To calculate the identity-by-descent (IBD) for a pair of SNP genotypes
#   using MLE
#
# INPUT:
#   geno1 -- SNP profile of the first individual
#   geno2 -- SNP profile of the second individual
#   allele.freq --  allele frequencies
#   method -- the searching algorithm
#   kinship.constraint -- restrict the IBD coefficients
#   max.niter -- the maximum number of iterations
#   reltol -- the relative tolerance
#   coeff.correct -- internal use
#

snpgdsPairIBD <- function(geno1, geno2, allele.freq,
	method=c("EM", "downhill.simplex", "MoM"), kinship.constraint=FALSE, max.niter=1000,
	reltol=sqrt(.Machine$double.eps), coeff.correct=TRUE,
	out.num.iter=TRUE, verbose=TRUE)
{
	# check
	stopifnot(is.vector(geno1) & is.numeric(geno1))
	stopifnot(is.vector(geno2) & is.numeric(geno2))
	stopifnot(is.vector(allele.freq) & is.numeric(allele.freq))
	stopifnot(length(geno1) == length(geno2))
	stopifnot(length(geno1) == length(allele.freq))
	stopifnot(is.logical(kinship.constraint))
	stopifnot(is.logical(coeff.correct))
	method <- match.arg(method)

	# method
	if (method == "EM")
		method <- 0
	else if (method == "downhill.simplex")
		method <- 1
	else if (method == "MoM")
		method <- -1
	else
		stop("Invalid MLE method!")

	allele.freq[!is.finite(allele.freq)] <- -1
	flag <- (0 <= allele.freq) & (allele.freq <= 1)
	if (sum(flag) < length(geno1))
	{
		if (verbose)
		{
			cat(sprintf(
				"IBD MLE for %d SNPs in total, after removing loci with invalid allele frequencies.\n",
				sum(flag)))
		}
		geno1 <- geno1[flag]; geno2 <- geno2[flag]
		allele.freq <- allele.freq[flag]
	}

	# call C codes
	rv <- .C("gnrPairIBD", length(geno1), as.integer(geno1), as.integer(geno2),
		as.double(allele.freq), as.logical(kinship.constraint), as.integer(max.niter),
		as.double(reltol), as.logical(coeff.correct), as.integer(method),
		out.k0 = double(1), out.k1 = double(1), out.loglik = double(1),
		out.niter = integer(1), double(3*length(geno1) + 3*4),
		err = integer(1), NAOK=TRUE, PACKAGE="SNPRelate")
	if (rv$err != 0) stop(snpgdsErrMsg())

	# return
	ans <- data.frame(k0=rv$out.k0, k1=rv$out.k1, loglik=rv$out.loglik)
	if (out.num.iter) ans$niter <- rv$out.niter
	ans
}



#######################################################################
# To calculate the identity-by-descent (IBD) matrix (MLE) for SNP genotypes
#
# INPUT:
#   geno1 -- SNP profile of the first individual
#   geno2 -- SNP profile of the second individual
#   allele.freq --  allele frequencies
#   k0 -- the first IBD coefficient
#   k1 -- the second IBD coefficient
#

snpgdsPairIBDMLELogLik <- function(geno1, geno2, allele.freq, k0=NaN, k1=NaN,
	relatedness=c("", "self", "fullsib", "offspring", "halfsib", "cousin", "unrelated"),
	verbose=TRUE)
{
	# check
	stopifnot(is.vector(geno1) & is.numeric(geno1))
	stopifnot(is.vector(geno2) & is.numeric(geno2))
	stopifnot(is.vector(allele.freq) & is.numeric(allele.freq))
	stopifnot(length(geno1) == length(geno2))
	stopifnot(length(geno1) == length(allele.freq))
	stopifnot(is.numeric(k0))
	stopifnot(is.numeric(k1))
	stopifnot(is.character(relatedness))

	allele.freq[!is.finite(allele.freq)] <- -1
	flag <- (0 <= allele.freq) & (allele.freq <= 1)
	if (sum(flag) < length(geno1))
	{
		if (verbose)
		{
			cat(sprintf(
				"IBD MLE for %d SNPs in total, after removing loci with invalid allele frequencies.\n",
				sum(flag)))
		}
		geno1 <- geno1[flag]; geno2 <- geno2[flag]
		allele.freq <- allele.freq[flag]
	}

	# relatedness
	relatedness <- relatedness[1]
	if (relatedness == "self")
	{
		k0 <- 0; k1 <- 0
	} else if (relatedness == "fullsib")
	{
		k0 <- 0.25; k1 <- 0.5
	} else if (relatedness == "offspring")
	{
		k0 <- 0; k1 <- 1
	} else if (relatedness == "halfsib")
	{
		k0 <- 0.5; k1 <- 0.5
	} else if (relatedness == "cousin")
	{
		k0 <- 0.75; k1 <- 0.25
	} else if (relatedness == "unrelated")
	{
		k0 <- 1; k1 <- 0
	}

	# call C codes
	rv <- .C("gnrPairIBDLogLik", length(geno1), as.integer(geno1), as.integer(geno2),
		as.double(allele.freq), as.double(k0), as.double(k1), out.loglik = double(1),
		double(3*length(geno1)), err = integer(1), NAOK=TRUE, PACKAGE="SNPRelate")
	if (rv$err != 0) stop(snpgdsErrMsg())

	# return
	rv$out.loglik
}



#######################################################################
# Identity-by-Descent (IBD) analysis using KING robust estimat
#######################################################################

#######################################################################
# To calculate the identity-by-descent (IBD) matrix (KING) for SNP genotypes
#
# INPUT:
#   gdsobj -- an object of gds file
#   sample.id -- a vector of sample.id; if NULL, to use all samples
#   snp.id -- a vector of snp.id; if NULL, to use all SNPs
#   autosome.only -- whether only use autosomal SNPs
#   remove.monosnp -- whether remove monomorphic snps or not
#   maf -- the threshold of minor allele frequencies, keeping ">= maf"
#   missing.rate -- the threshold of missing rates, keeping "<= missing.rate"
#   type -- "KING-home" or "KING-robust"
#   family.id -- NULL or a list of family id
#   num.thread -- the number of threads to be used
#   verbose -- show information
#

snpgdsIBDKING <- function(gdsobj, sample.id=NULL, snp.id=NULL, autosome.only=TRUE,
	remove.monosnp=TRUE, maf=NaN, missing.rate=NaN,
	type=c("KING-robust", "KING-homo"), family.id=NULL,
	num.thread=1, verbose=TRUE)
{
	# check
	stopifnot(inherits(gdsobj, "gds.class"))
	stopifnot(is.logical(autosome.only))
	stopifnot(is.logical(remove.monosnp))
	type <- match.arg(type)
	stopifnot(is.null(family.id) | is.vector(family.id))
	stopifnot(is.numeric(num.thread) & (num.thread>0))
	stopifnot(is.logical(verbose))
	allele.freq <- NULL

	# samples
	sample.ids <- read.gdsn(index.gdsn(gdsobj, "sample.id"))
	if (!is.null(sample.id))
	{
		tmp <- sample.id
		n.tmp <- length(sample.id)
		sample.id <- sample.ids %in% sample.id
		n.samp <- sum(sample.id);
		if (n.samp != n.tmp)
			stop("Some of sample.id do not exist!")
		if (n.samp <= 0)
			stop("No sample in the working dataset.")
		sample.ids <- sample.ids[sample.id]

		# family id
		if (is.vector(family.id))
		{
			if (length(tmp) != length(family.id))
				stop("'family.id' should have the same length and order as 'sample.id'.")
			family.id <- family.id[match(sample.ids, tmp)]
		}
	} else {
		if (is.vector(family.id))
		{
			if (length(sample.ids) != length(family.id))
				stop("'family.id' should have the same length as 'sample.id' in the GDS file")
		}
	}

	if (verbose)
		cat("Identity-By-Descent analysis (KING method of moment) on SNP genotypes:\n");

	# family id
	if (is.vector(family.id))
	{
		if (is.character(family.id))
			family.id[family.id == ""] <- NA
		family.id <- as.factor(family.id)
		if (verbose & (type=="KING-robust"))
		{
			cat("# of families: ", nlevels(family.id),
				", and within- and between-family relationship are estimated differently.\n",
				sep="")
		}
	} else {
		if (verbose & (type=="KING-robust"))
		{
			cat("No family is specified, and all individuals are treated as singletons.\n")
		}
		family.id <- rep(NA, length(sample.ids))
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
		if (!is.null(allele.freq))
		{
			if (n.snp != length(allele.freq))
				stop("`length(allele.freq)' should equal to the number of SNPs.")
		}
		if (autosome.only)
		{
			chr <- read.gdsn(index.gdsn(gdsobj, "snp.chromosome"))
			opt <- snpgdsOption(gdsobj)
			chridx <- chr[snp.id]
			snp.id <- snp.id & (chr %in% c(opt$autosome.start : opt$autosome.end))
			if (!is.null(allele.freq))
				allele.freq <- allele.freq[chridx %in% 1:22]
			if (verbose)
			{
				tmp <- n.snp - sum(snp.id)
				if (tmp > 0) cat("Removing", tmp, "non-autosomal SNPs\n")
			}
		}
		snp.ids <- snp.ids[snp.id]
	} else {
		if (!is.null(allele.freq))
		{
			if (length(snp.ids) != length(allele.freq))
				stop("`length(allele.freq)' should equal to the number of SNPs.")
		}
		if (autosome.only)
		{
			chr <- read.gdsn(index.gdsn(gdsobj, "snp.chromosome"))
			opt <- snpgdsOption(gdsobj)
			snp.id <- chr %in% c(opt$autosome.start : opt$autosome.end)
			snp.ids <- snp.ids[snp.id]
			if (!is.null(allele.freq))
				allele.freq <- allele.freq[snp.id]
			if (verbose)
				cat("Removing", length(chr) - length(snp.ids), "non-autosomal SNPs\n")
		}
	}

	# check
	if (!is.null(allele.freq))
	{
		cat(sprintf("Specifying allele frequencies, mean: %0.3f, sd: %0.3f\n",
			mean(allele.freq, na.rm=TRUE), sd(allele.freq, na.rm=TRUE)))
	}

	# call C codes
	# set genotype working space
	node <- .C("gnrSetGenoSpace", as.integer(index.gdsn(gdsobj, "genotype")),
		as.logical(sample.id), as.logical(!is.null(sample.id)),
		as.logical(snp.id), as.logical(!is.null(snp.id)),
		n.snp=integer(1), n.samp=integer(1),
		err=integer(1), NAOK=TRUE, PACKAGE="SNPRelate")
	if (node$err != 0) stop(snpgdsErrMsg())

	# call allele frequencies and missing rates
	if (remove.monosnp || is.finite(maf) || is.finite(missing.rate))
	{
		if (!is.finite(maf)) maf <- -1;
		if (!is.finite(missing.rate)) missing.rate <- 2;

		# call
		if (is.null(allele.freq))
		{
			rv <- .C("gnrSelSNP_Base", as.logical(remove.monosnp),
				as.double(maf), as.double(missing.rate),
				out.num=integer(1), out.snpflag = logical(node$n.snp),
				err=integer(1), NAOK=TRUE, PACKAGE="SNPRelate")
		} else {
			rv <- .C("gnrSelSNP_Base_Ex", as.double(allele.freq), as.logical(remove.monosnp),
				as.double(maf), as.double(missing.rate),
				out.num=integer(1), out.snpflag = logical(node$n.snp),
				err=integer(1), NAOK=TRUE, PACKAGE="SNPRelate")
		}
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
		if (num.thread <= 1)
			cat("\tUsing", num.thread, "CPU core.\n")
		else
			cat("\tUsing", num.thread, "CPU cores.\n")
	}

	if (type == "KING-homo")
	{
		if (verbose)
			cat("Relationship inference in a homogeneous population.\n")

		# call the C function
		rv <- .Call("gnrIBD_KING_Homo", as.logical(verbose), TRUE,
			as.integer(num.thread), PACKAGE="SNPRelate")

		rv <- list(sample.id=sample.ids, snp.id=snp.ids, afreq=NULL,
			k0=rv[[1]], k1=rv[[2]])

	} else if (type == "KING-robust")
	{
		if (verbose)
			cat("Relationship inference in the presence of population stratification.\n")

		# call the C function
		rv <- .Call("gnrIBD_KING_Robust", as.logical(verbose), TRUE,
			as.integer(num.thread), as.integer(family.id), PACKAGE="SNPRelate")

		rv <- list(sample.id=sample.ids, snp.id=snp.ids, afreq=NULL,
			IBS0=rv[[1]], kinship=rv[[2]])
	} else
		stop("Invalid 'type'.")

	# return
	rv$afreq[rv$afreq < 0] <- NaN
	class(rv) <- "snpgdsIBDClass"
	return(rv)
}




#######################################################################
#
#######################################################################

#######################################################################
# Return a data.frame of pairs of individuals with IBD coefficients
#
# INPUT:
#   ibdobj -- an object of snpgdsIBDClass
#   kinship.cutoff -- a threshold of kinship coefficient
#   samp.sel -- a logical vector, or integer vector
#
# OUTPUT:
#   a data.frame of "ID1", "ID2", "k0", "k1", "kinship"
#

snpgdsIBDSelection <- function(ibdobj, kinship.cutoff=NaN, samp.sel=NULL)
{
	# check
	stopifnot(inherits(ibdobj, "snpgdsIBDClass"))
	stopifnot(is.numeric(kinship.cutoff))
	stopifnot(is.null(samp.sel) | is.logical(samp.sel) | is.numeric(samp.sel))
	if (is.logical(samp.sel))
		stopifnot(length(samp.sel) == length(ibdobj$sample.id))

	# the variables in the output
	ns <- setdiff(names(ibdobj), c("sample.id", "snp.id", "afreq"))

	# subset
	if (!is.null(samp.sel))
	{
		ibdobj$sample.id <- ibdobj$sample.id[samp.sel]
		for (i in ns)
			ibdobj[[i]] <- ibdobj[[i]][samp.sel, samp.sel]
	}

	if (is.null(ibdobj$kinship))
	{
		ibdobj$kinship <- (1 - ibdobj$k0 - ibdobj$k1)*0.5 + ibdobj$k1*0.25
		ns <- c(ns, "kinship")
	}
	if (is.finite(kinship.cutoff))
	{
		flag <- lower.tri(ibdobj$kinship) & (ibdobj$kinship >= kinship.cutoff)
		flag[is.na(flag)] <- FALSE
	} else
		flag <- lower.tri(ibdobj$kinship)

	# output
	n <- length(ibdobj$sample.id)
	ans <- data.frame(
		ID1 = .Call("gnrIBDSelSampID_List1", ibdobj$sample.id, flag,
			PACKAGE="SNPRelate"),
		ID2 = .Call("gnrIBDSelSampID_List2", ibdobj$sample.id, flag,
			PACKAGE="SNPRelate"),
		stringsAsFactors=FALSE)
	for (i in ns)
		ans[[i]] <- ibdobj[[i]][flag]

	ans
}
