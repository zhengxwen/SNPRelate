#######################################################################
#
# Package name: SNPRelate
#
# Description:
#     A High-performance Computing Toolset for Relatedness and
# Principal Component Analysis of SNP Data
#
# Copyright (C) 2011 - 2015        Xiuwen Zheng
# License: GPL-3
# Email: zhengxwen@gmail.com
#


#######################################################################
# Principal Component Analysis
#######################################################################

#######################################################################
# To conduct Principal Component Analysis
#

snpgdsPCA <- function(gdsobj, sample.id=NULL, snp.id=NULL,
    autosome.only=TRUE, remove.monosnp=TRUE, maf=NaN, missing.rate=NaN,
    eigen.cnt=32L, num.thread=1L, bayesian=FALSE, need.genmat=FALSE,
    genmat.only=FALSE, eigen.method=c("DSPEVX", "DSPEV"), verbose=TRUE)
{
    # check
    ws <- .InitFile2(
        cmd="Principal Component Analysis (PCA) on SNP genotypes:",
        gdsobj=gdsobj, sample.id=sample.id, snp.id=snp.id,
        autosome.only=autosome.only, remove.monosnp=remove.monosnp,
        maf=maf, missing.rate=missing.rate, num.thread=num.thread,
        verbose=verbose)

    stopifnot(is.numeric(eigen.cnt))
    stopifnot(is.logical(bayesian))
    stopifnot(is.logical(need.genmat))
    stopifnot(is.logical(genmat.only))
    if (genmat.only) need.genmat <- TRUE
    if (eigen.cnt <= 0L) eigen.cnt <- ws$n.samp

	eigen.method <- match.arg(eigen.method)

    # call parallel PCA
    rv <- .Call(gnrPCA, eigen.cnt, ws$num.thread,
        bayesian, need.genmat, genmat.only, eigen.method, verbose)

    # return
    rv <- list(sample.id = ws$sample.id, snp.id = ws$snp.id,
        eigenval = rv[[3]], eigenvect = rv[[4]],
        varprop = rv[[3]] / rv[[5]],
        TraceXTX = rv[[1]], Bayesian = bayesian, genmat = rv[[2]])
    class(rv) <- "snpgdsPCAClass"
    return(rv)
}



#######################################################################
# To calculate SNP correlations from principal component analysis
#

snpgdsPCACorr <- function(pcaobj, gdsobj, snp.id=NULL, eig.which=NULL,
    num.thread=1L, verbose=TRUE)
{
    # check
    stopifnot(inherits(pcaobj, "snpgdsPCAClass"))
    ws <- .InitFile(gdsobj, sample.id=pcaobj$sample.id, snp.id=snp.id)

    snp.id <- read.gdsn(index.gdsn(gdsobj, "snp.id"))
    if (!is.null(ws$snp.flag))
        snp.id <- snp.id[ws$snp.flag]

    stopifnot(is.numeric(num.thread) & (num.thread>0))
    stopifnot(is.logical(verbose))

    if (is.null(eig.which))
    {
        eig.which <- 1:dim(pcaobj$eigenvect)[2]
    } else {
        stopifnot(is.vector(eig.which))
        stopifnot(is.numeric(eig.which))
        eig.which <- as.integer(eig.which)
    }

    if (verbose)
    {
        cat("SNP correlations:\n")
        cat("Working space:", ws$n.samp, "samples,", ws$n.snp, "SNPs\n");
        if (num.thread <= 1)
            cat("\tUsing", num.thread, "(CPU) core.\n")
        else
            cat("\tUsing", num.thread, "(CPU) cores.\n")
        cat("\tUsing the top", dim(pcaobj$eigenvect)[2], "eigenvectors.\n")
    }

    # call C function
    rv <- .Call(gnrPCACorr, length(eig.which), pcaobj$eigenvect[,eig.which],
        as.integer(num.thread), verbose)

    # return
    list(sample.id=pcaobj$sample.id, snp.id=snp.id, snpcorr=rv)
}



#######################################################################
# To calculate SNP loadings from principal component analysis
#

snpgdsPCASNPLoading <- function(pcaobj, gdsobj, num.thread=1L, verbose=TRUE)
{
    # check
    stopifnot(inherits(pcaobj, "snpgdsPCAClass"))
    ws <- .InitFile(gdsobj, sample.id=pcaobj$sample.id, snp.id=pcaobj$snp.id)
    stopifnot(is.numeric(num.thread) & (num.thread>0))
    stopifnot(is.logical(verbose))

    if (verbose)
    {
        cat("SNP loading:\n")
        cat("Working space:", ws$n.samp, "samples,", ws$n.snp, "SNPs\n");
        if (num.thread <= 1)
            cat("\tUsing", num.thread, "(CPU) core.\n")
        else
            cat("\tUsing", num.thread, "(CPU) cores.\n")
        cat("\tUsing the top", dim(pcaobj$eigenvect)[2], "eigenvectors.\n")
    }

    # call parallel PCA
    rv <- .Call(gnrPCASNPLoading, pcaobj$eigenval, dim(pcaobj$eigenvect),
        pcaobj$eigenvect, pcaobj$TraceXTX, as.integer(num.thread),
        pcaobj$Bayesian, verbose)

    # return
    rv <- list(sample.id=pcaobj$sample.id, snp.id=pcaobj$snp.id,
        eigenval=pcaobj$eigenval, snploading=rv[[1]],
        TraceXTX=pcaobj$TraceXTX, Bayesian=pcaobj$Bayesian,
        avefreq=rv[[2]], scale=rv[[3]])
    class(rv) <- "snpgdsPCASNPLoadingClass"
    return(rv)
}



#######################################################################
# To calculate sample loadings from SNP loadings in PCA
#

snpgdsPCASampLoading <- function(loadobj, gdsobj, sample.id=NULL,
    num.thread=1L, verbose=TRUE)
{
    # check
    stopifnot(inherits(loadobj, "snpgdsPCASNPLoadingClass"))

    ws <- .InitFile(gdsobj, sample.id=sample.id, snp.id=loadobj$snp.id)

    sample.id <- read.gdsn(index.gdsn(gdsobj, "sample.id"))
    if (!is.null(ws$samp.flag))
        sample.id <- sample.id[ws$samp.flag]

    stopifnot(is.numeric(num.thread) & (num.thread>0))
    stopifnot(is.logical(verbose))

    eigcnt <- dim(loadobj$snploading)[1]
    if (verbose)
    {
        cat("Sample loading:\n")
        cat("Working space:", ws$n.samp, "samples,", ws$n.snp, "SNPs\n");
        if (num.thread <= 1)
            cat("\tUsing", num.thread, "(CPU) core.\n")
        else
            cat("\tUsing", num.thread, "(CPU) cores.\n")
        cat("\tUsing the top", eigcnt, "eigenvectors.\n")
    }

    # call C function
    rv <- .Call(gnrPCASampLoading, length(loadobj$sample.id),
        loadobj$eigenval, eigcnt, loadobj$snploading, loadobj$TraceXTX,
        loadobj$avefreq, loadobj$scale, as.integer(num.thread), verbose)

    # return
    rv <- list(sample.id = sample.id, snp.id = loadobj$snp.id,
        eigenval = loadobj$eigenval, eigenvect = rv,
        TraceXTX = loadobj$TraceXTX,
        Bayesian = loadobj$Bayesian, genmat = NULL)
    class(rv) <- "snpgdsPCAClass"
    return(rv)
}



#######################################################################
# Eigen-Analysis on SNP genotypes
#

snpgdsEIGMIX <- function(gdsobj, sample.id=NULL, snp.id=NULL,
    autosome.only=TRUE, remove.monosnp=TRUE, maf=NaN, missing.rate=NaN,
    num.thread=1L, eigen.cnt=32, need.ibdmat=FALSE, ibdmat.only=FALSE,
    method=c("Coancestry", "GRM"), verbose=TRUE)
{
    # check and initialize ...
    ws <- .InitFile2(cmd="Eigen-analysis on SNP genotypes:",
        gdsobj=gdsobj, sample.id=sample.id, snp.id=snp.id,
        autosome.only=autosome.only, remove.monosnp=remove.monosnp,
        maf=maf, missing.rate=missing.rate, num.thread=num.thread,
        verbose=verbose)

    stopifnot(is.numeric(eigen.cnt) & is.vector(eigen.cnt))
    stopifnot(length(eigen.cnt) == 1)
    eigen.cnt <- as.integer(eigen.cnt)
    if (eigen.cnt < 1)
        eigen.cnt <- as.integer(ws$n.samp)

    stopifnot(is.logical(need.ibdmat) & is.vector(need.ibdmat))
    stopifnot(length(need.ibdmat) == 1)

    stopifnot(is.logical(ibdmat.only) & is.vector(ibdmat.only))
    stopifnot(length(ibdmat.only) == 1)

    method <- match.arg(method)
    diag.adjustment <- (method == "Coancestry")

    in.method <- c("RatioOfAverages", "AverageOfRatios")
    in.method <- as.character(in.method)[1]
    in.method <- match(in.method, c("RatioOfAverages", "AverageOfRatios"))


    #########################################################
    # call eigen-analysis

    rv <- .Call(gnrEIGMIX, eigen.cnt, ws$num.thread, need.ibdmat,
        ibdmat.only, in.method, diag.adjustment, verbose)

    # return
    rv <- list(sample.id = ws$sample.id, snp.id = ws$snp.id,
        eigenval = rv$eigenval, eigenvect = rv$eigenvect,
        ibdmat = rv$ibdmat)
    class(rv) <- "snpgdsEigMixClass"
    return(rv)
}



#######################################################################
# Admixture proportion from eigen-analysis
#

snpgdsAdmixProp <- function(eigobj, groups, bound=FALSE)
{
    # check
    stopifnot( inherits(eigobj, "snpgdsEigMixClass") |
        inherits(eigobj, "snpgdsPCAClass") )
    # 'sample.id' and 'eigenvect' should exist
    stopifnot(!is.null(eigobj$sample.id))
    stopifnot(is.matrix(eigobj$eigenvect))

    stopifnot(is.list(groups))
    stopifnot(length(groups) > 1)
    if (length(groups) > (ncol(eigobj$eigenvect)+1))
    {
        stop("`eigobj' should have more eigenvectors than ",
            "what is specified in `groups'.")
    }

    grlist <- NULL
    for (i in 1:length(groups))
    {
        if (!is.vector(groups[[i]]) & !is.factor(groups[[i]]))
        {
            stop(
                "`groups' should be a list of sample IDs ",
                "with respect to multiple groups."
            )
        }
        if (any(!(groups[[i]] %in% eigobj$sample.id)))
        {
            stop(sprintf(
                "`groups[[%d]]' includes sample(s) not existing ",
                "in `eigobj$sample.id'.", i))
        }

        if (any(groups[[i]] %in% grlist))
            warning("There are some overlapping between group sample IDs.")
        grlist <- c(grlist, groups[[i]])
    }

    stopifnot(is.logical(bound) & is.vector(bound))
    stopifnot(length(bound) == 1)

    # calculate ...

    E <- eigobj$eigenvect[, 1:(length(groups)-1)]
    if (is.vector(E)) E <- matrix(E, ncol=1)
    mat <- NULL
    for (i in 1:length(groups))
    {
        k <- match(groups[[i]], eigobj$sample.id)
        Ek <- E[k, ]
        if (is.vector(Ek))
            mat <- rbind(mat, mean(Ek))
        else
            mat <- rbind(mat, colMeans(Ek))
    }

    # check
    if (any(is.na(mat)))
        stop("The eigenvectors should not have missing value!")

    T.P <- mat[length(groups), ]
    T.R <- solve(mat[-length(groups), ] -
        matrix(T.P, nrow=length(T.P), ncol=length(T.P), byrow=TRUE))

    new.p <- (E - matrix(T.P, nrow=dim(E)[1], ncol=length(T.P),
        byrow=TRUE)) %*% T.R
    new.p <- cbind(new.p, 1 - rowSums(new.p))
    colnames(new.p) <- names(groups)
    rownames(new.p) <- eigobj$sample.id

    # whether bounded
    if (bound)
    {
        new.p[new.p < 0] <- 0
        r <- 1.0 / rowSums(new.p)
        new.p <- new.p * r
    }

    new.p
}
