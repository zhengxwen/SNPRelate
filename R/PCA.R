#######################################################################
#
# Package name: SNPRelate
#
# Description:
#     A High-performance Computing Toolset for Relatedness and
# Principal Component Analysis of SNP Data
#
# Copyright (C) 2011 - 2017        Xiuwen Zheng
# License: GPL-3
# Email: zhengxwen@gmail.com
#


#######################################################################
# Principal Component Analysis
#######################################################################

#######################################################################
# Conduct Principal Component Analysis (PCA)
#

snpgdsPCA <- function(gdsobj, sample.id=NULL, snp.id=NULL,
    autosome.only=TRUE, remove.monosnp=TRUE, maf=NaN, missing.rate=NaN,
    algorithm=c("exact", "randomized"),
    eigen.cnt=ifelse(identical(algorithm, "randomized"), 16L, 32L),
    num.thread=1L, bayesian=FALSE, need.genmat=FALSE,
    genmat.only=FALSE, eigen.method=c("DSPEVX", "DSPEV"),
    aux.dim=eigen.cnt*2L, iter.num=10L, verbose=TRUE)
{
    # check
    ws <- .InitFile2(
        cmd="Principal Component Analysis (PCA) on genotypes:",
        gdsobj=gdsobj, sample.id=sample.id, snp.id=snp.id,
        autosome.only=autosome.only, remove.monosnp=remove.monosnp,
        maf=maf, missing.rate=missing.rate, num.thread=num.thread,
        verbose=verbose)

    stopifnot(is.numeric(eigen.cnt))
    stopifnot(is.logical(bayesian))
    stopifnot(is.logical(need.genmat))
    stopifnot(is.logical(genmat.only))
    algorithm <- match.arg(algorithm)
    eigen.method <- match.arg(eigen.method)

    if (genmat.only) need.genmat <- TRUE
    if (eigen.cnt <= 0L) eigen.cnt <- ws$n.samp

    # call parallel PCA
    param <- list(bayesian=bayesian, need.genmat=need.genmat,
        genmat.only=genmat.only, eigen.method=eigen.method,
        aux.dim=aux.dim, iter.num=iter.num)
    if (algorithm == "randomized")
        param$aux.mat <- rnorm(aux.dim * ws$n.samp)
    rv <- .Call(gnrPCA, eigen.cnt, algorithm, ws$num.thread, param, verbose)

    # return
    if (algorithm == "exact")
    {
        rv <- list(sample.id = ws$sample.id, snp.id = ws$snp.id,
            eigenval = rv[[3L]], eigenvect = rv[[4L]],
            varprop = rv[[3L]] / rv[[5L]],
            TraceXTX = rv[[1L]], Bayesian = bayesian, genmat = rv[[2L]])
        class(rv) <- "snpgdsPCAClass"
    } else if (algorithm == "randomized")
    {
        rv <- list(sample.id = ws$sample.id, snp.id = ws$snp.id,
            eigenval = rv[[1L]], eigenvect = t(rv[[2]][seq_len(eigen.cnt), ]),
            varprop = NaN,
            TraceXTX = NaN, Bayesian = FALSE)
        class(rv) <- "snpgdsPCAClass"
    }
    return(rv)
}



#######################################################################
# To calculate SNP correlations from principal component analysis
#

snpgdsPCACorr <- function(pcaobj, gdsobj, snp.id=NULL, eig.which=NULL,
    num.thread=1L, with.id=TRUE, outgds=NULL, verbose=TRUE)
{
    # check
    stopifnot(inherits(pcaobj, "snpgdsPCAClass") |
        inherits(pcaobj, "snpgdsEigMixClass") | is.matrix(pcaobj))
    if (is.matrix(pcaobj))
    {
        sampid <- rownames(pcaobj)
        if (is.null(sampid))
            stop("'rownames(pcaobj)' should be sample id.")
        eigenvect <- pcaobj
    } else {
        sampid <- pcaobj$sample.id
        eigenvect <- pcaobj$eigenvect
    }

    stopifnot(is.null(outgds) | is.character(outgds))
    if (is.character(outgds))
    {
        stopifnot(length(outgds) == 1L)
        with.id <- TRUE
    }

    ws <- .InitFile(gdsobj, sample.id=sampid, snp.id=snp.id, with.id=with.id)

    stopifnot(is.numeric(num.thread), length(num.thread)==1L, num.thread>0L)
    if (length(sampid) != nrow(eigenvect))
    {
        stop("Internal error: the number of samples should be ",
            "equal to the number of rows in 'eigenvect'.")
    }
    stopifnot(is.logical(verbose), length(verbose)==1L)

    if (is.null(eig.which))
    {
        eig.which <- seq_len(dim(eigenvect)[2L])
    } else {
        stopifnot(is.vector(eig.which))
        stopifnot(is.numeric(eig.which))
        eig.which <- as.integer(eig.which)
    }

    if (verbose)
    {
        cat("SNP Correlation:\n")
        cat("Working space:", ws$n.samp, "samples,", ws$n.snp, "SNPs\n");
        cat("    using ", num.thread, " (CPU) core", .plural(num.thread), "\n",
            sep="")
        cat("    using the top", length(eig.which), "eigenvectors\n")
    }

    gds <- NULL
    if (is.character(outgds))
    {
        f <- createfn.gds(outgds)
        on.exit({ closefn.gds(f) })
        add.gdsn(f, "sample.id", sampid, compress="LZMA_RA", closezip=TRUE)
        add.gdsn(f, "snp.id", ws$snp.id, compress="LZMA_RA", closezip=TRUE)
        sync.gds(f)
        gds <- add.gdsn(f, "correlation", storage="packedreal16",
            valdim=c(length(eig.which), 0L), compress="LZMA_RA")
        if (verbose)
            cat("Creating '", outgds, "' ...\n", sep="")
    }

    # call C function
    rv <- .Call(gnrPCACorr, length(eig.which), eigenvect[,eig.which],
        num.thread, gds, verbose)

    # return
    if (is.null(outgds))
    {
        if (isTRUE(with.id))
            rv <- list(sample.id=sampid, snp.id=ws$snp.id, snpcorr=rv)
	    rv
    } else
        invisible()
}



#######################################################################
# To calculate SNP loadings from principal component analysis
#

snpgdsPCASNPLoading <- function(pcaobj, gdsobj, num.thread=1L, verbose=TRUE)
{
    # check
    stopifnot(inherits(pcaobj, "snpgdsPCAClass") |
        inherits(pcaobj, "snpgdsEigMixClass"))
    ws <- .InitFile(gdsobj, sample.id=pcaobj$sample.id, snp.id=pcaobj$snp.id)
    stopifnot(is.numeric(num.thread), num.thread > 0L)
    stopifnot(is.logical(verbose), length(verbose)==1L)
    stopifnot(!is.null(pcaobj$eigenval), !is.null(pcaobj$eigenvect))

    if (verbose)
    {
        cat("SNP loading:\n")
        cat("Working space:", ws$n.samp, "samples,", ws$n.snp, "SNPs\n");
        cat("    using ", num.thread, " (CPU) core", .plural(num.thread), "\n",
            sep="")
        cat("    using the top", dim(pcaobj$eigenvect)[2L], "eigenvectors\n")
    }

    # computing SNP loading
    if (inherits(pcaobj, "snpgdsPCAClass"))
    {
        # call C function
        rv <- .Call(gnrPCASNPLoading, pcaobj$eigenval, pcaobj$eigenvect,
            pcaobj$TraceXTX, num.thread, pcaobj$Bayesian, verbose)
        rv <- list(sample.id=pcaobj$sample.id, snp.id=pcaobj$snp.id,
            eigenval=pcaobj$eigenval, snploading=rv[[1L]],
            TraceXTX=pcaobj$TraceXTX, Bayesian=pcaobj$Bayesian,
            avgfreq=rv[[2L]], scale=rv[[3L]])
        class(rv) <- "snpgdsPCASNPLoadingClass"

    } else {
        if (isTRUE(pcaobj$diagadj))
        {
            stop(
    "Please run `snpgdsEIGMIX(, diagadj=FALSE)` for projecting new samples.")
        }
        # call C function
        mat <- .Call(gnrEigMixSNPLoading, pcaobj$eigenval, pcaobj$eigenvect,
            pcaobj$afreq, num.thread, verbose)
        rv <- list(sample.id=pcaobj$sample.id, snp.id=pcaobj$snp.id,
            eigenval=pcaobj$eigenval, snploading=mat,
            afreq=pcaobj$afreq)
        class(rv) <- "snpgdsEigMixSNPLoadingClass"
    }
    return(rv)
}



#######################################################################
# To calculate sample loadings from SNP loadings in PCA
#

snpgdsPCASampLoading <- function(loadobj, gdsobj, sample.id=NULL,
    num.thread=1L, verbose=TRUE)
{
    # check
    stopifnot(inherits(loadobj, "snpgdsPCASNPLoadingClass") |
        inherits(loadobj, "snpgdsEigMixSNPLoadingClass"))
    ws <- .InitFile(gdsobj, sample.id=sample.id, snp.id=loadobj$snp.id)

    sample.id <- read.gdsn(index.gdsn(gdsobj, "sample.id"))
    if (!is.null(ws$samp.flag))
        sample.id <- sample.id[ws$samp.flag]
    stopifnot(is.numeric(num.thread), length(num.thread)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L)

    eigcnt <- nrow(loadobj$snploading)
    if (verbose)
    {
        cat("Sample loading:\n")
        cat("Working space:", ws$n.samp, "samples,", ws$n.snp, "SNPs\n")
        cat("    using ", num.thread, " (CPU) core", .plural(num.thread), "\n",
            sep="")
        cat("    using the top", eigcnt, "eigenvectors\n")
    }

    if (inherits(loadobj, "snpgdsPCASNPLoadingClass"))
    {
        # prepare post-eigenvectors
        ss <- (length(loadobj$sample.id) - 1) / loadobj$TraceXTX
        sqrt_eigval <- sqrt(ss / loadobj$eigenval[1:eigcnt])
        sload <- loadobj$snploading * sqrt_eigval

        # call C function
        mm <- .Call(gnrPCASampLoading, eigcnt, sload, loadobj$avgfreq,
            loadobj$scale, num.thread, verbose)

        # return
        rv <- list(sample.id = sample.id, snp.id = loadobj$snp.id,
            eigenval = loadobj$eigenval, eigenvect = mm,
            TraceXTX = loadobj$TraceXTX,
            Bayesian = loadobj$Bayesian, genmat = NULL)
        class(rv) <- "snpgdsPCAClass"

    } else {
        # prepare post-eigenvectors
        sqrt_eigval <- sqrt(1 / loadobj$eigenval[1:eigcnt])
        sload <- loadobj$snploading * sqrt_eigval

        # call C function
        mm <- .Call(gnrEigMixSampLoading, sload, loadobj$afreq, num.thread,
            verbose)

        # return
        rv <- list(sample.id = sample.id, snp.id = loadobj$snp.id,
            eigenval = loadobj$eigenval, eigenvect = mm,
            afreq = loadobj$afreq)
        class(rv) <- "snpgdsEigMixClass"
    }
    return(rv)
}



#######################################################################
# Eigen-Analysis on genotypes
#

snpgdsEIGMIX <- function(gdsobj, sample.id=NULL, snp.id=NULL,
    autosome.only=TRUE, remove.monosnp=TRUE, maf=NaN, missing.rate=NaN,
    num.thread=1L, eigen.cnt=32L, diagadj=TRUE, ibdmat=FALSE, verbose=TRUE)
{
    # check and initialize ...
    ws <- .InitFile2(cmd="Eigen-analysis on genotypes:",
        gdsobj=gdsobj, sample.id=sample.id, snp.id=snp.id,
        autosome.only=autosome.only, remove.monosnp=remove.monosnp,
        maf=maf, missing.rate=missing.rate, num.thread=num.thread,
        verbose=verbose)

    stopifnot(is.numeric(eigen.cnt), length(eigen.cnt)==1L)
    if (eigen.cnt < 0L)
        eigen.cnt <- ws$n.samp
    stopifnot(is.logical(diagadj), length(diagadj)==1L)
    stopifnot(is.logical(ibdmat), length(ibdmat)==1L)

    # call eigen-analysis
    param <- list(diagadj=diagadj, ibdmat=ibdmat)
    rv <- .Call(gnrEigMix, eigen.cnt, ws$num.thread, param, verbose)

    # return
    rv <- list(sample.id = ws$sample.id, snp.id = ws$snp.id,
        eigenval = rv[[1L]], eigenvect = rv[[2L]], afreq = rv[[3L]],
        ibd = rv[[4L]], diagadj=diagadj)
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



snpgdsAdmixPlot <- function(propmat, group=NULL, col=NULL, multiplot=TRUE,
    showgrp=TRUE, shownum=TRUE, ylim=TRUE, na.rm=TRUE)
{
    # check
    stopifnot(is.numeric(propmat), is.matrix(propmat))
    stopifnot(is.null(group) | is.vector(group) | is.factor(group))
    if (!is.null(group))
        stopifnot(nrow(propmat) == length(group))
    stopifnot(is.null(col) | is.vector(col))
    stopifnot(is.logical(multiplot), length(multiplot)==1L)
    stopifnot(is.logical(showgrp), length(showgrp)==1L)
    stopifnot(is.logical(shownum), length(shownum)==1L)
    stopifnot(is.logical(ylim) | is.numeric(ylim))
    if (is.numeric(ylim))
        stopifnot(length(ylim) == 2L)
    stopifnot(is.logical(na.rm), length(na.rm)==1L)

    if (is.logical(ylim))
    {
        if (isTRUE(ylim))
            ylim <- c(0, 1)
        else
            ylim <- range(propmat, na.rm=TRUE)
    }

    if (!is.null(group))
    {
        if (anyNA(group) & !isTRUE(na.rm))
        {
            group <- as.character(group)
            group[is.na(group)] <- "<NA>"
        }
        grp_name <- levels(factor(group))
        idx <- list()
        for (n in grp_name)
        {
            i <- which(group == n)
            i <- i[order(propmat[i, 1L], decreasing=TRUE)]
            idx <- c(idx, list(i))
        }
        propmat <- propmat[unlist(idx), ]
        grp_len <- lengths(idx, use.names=FALSE)
        xl <- c(0, cumsum(grp_len))
        x <- xl[-1L] - 0.5*grp_len
    }

    if (multiplot)
    {
        opar <- par(mfrow=c(ncol(propmat), 1L), mar=c(1.25, 5, 1.75, 2),
            oma=c(0, 0, 4, 0))
        on.exit(par(opar))
        grp <- colnames(propmat)
        if (is.null(grp))
            grp <- paste("group", seq_len(ncol(propmat)))
        ylab <- paste("Prop. of", grp)

        for (i in seq_len(ncol(propmat)))
        {
            barplot(unname(propmat[, i]), space=0, border=NA, ylab=ylab[i],
                ylim=ylim, col=col)
            lines(c(1, nrow(propmat)), c(0, 0))
            lines(c(1, nrow(propmat)), c(1, 1))
            if (!is.null(group))
            {
                abline(v=xl, col="blue")
                if (showgrp)
                    text(x, 0.5, labels=grp_name, srt=30)
                if ((i == 1L) & shownum)
                {
                    axis(1, c(0, x, nrow(propmat)),
                        c("", as.character(lengths(idx)), ""), cex.axis=0.75)
                }
            }
        }
    } else {
        if (is.null(col)) col <- rainbow(ncol(propmat))
        barplot(t(unname(propmat)), col=col,
            xlab="Individual #", ylab="Ancestry", space=0, border=NA)
        if (!is.null(group))
        {
            abline(v=xl, col="black")
            if (showgrp)
                text(x, 0.5, labels=grp_name, srt=30)
            if (shownum)
            {
                axis(1, c(0, x, nrow(propmat)),
                    c("", as.character(lengths(idx)), ""), cex.axis=0.75)
            }
        }
    }

    invisible()
}



snpgdsAdmixTable <- function(propmat, group, sort=FALSE)
{
    # check
    stopifnot(is.numeric(propmat), is.matrix(propmat))
    stopifnot(is.vector(group) | is.factor(group), nrow(propmat)==length(group))
    stopifnot(is.logical(sort), length(sort)==1L)

    ans <- vector("list", ncol(propmat))
    names(ans) <- colnames(propmat)
    for (i in seq_along(ans))
    {
        rv <- NULL
        for (grp in levels(factor(group)))
        {
            x <- group == grp
            x[is.na(x)] <- FALSE
            if (any(x))
            {
                y <- propmat[x, i]
                v <- data.frame(group=grp, num=sum(x),
                    mean = mean(y, na.rm=TRUE),  sd  = sd(y, na.rm=TRUE),
                    min  = min(y, na.rm=TRUE),   max = max(y, na.rm=TRUE),
                    stringsAsFactors=FALSE)
                rv <- rbind(rv, v)
            }
        }
        if (sort)
            rv <- rv[order(rv$mean, decreasing=TRUE), ]
        ans[[i]] <- rv
    }
    ans
}



#######################################################################
# plot PCA results
#

plot.snpgdsPCAClass <- function(x, eig=c(1L, 2L), ...)
{
    stopifnot(inherits(x, "snpgdsPCAClass"))
    stopifnot(is.numeric(eig), length(eig) >= 2L)

    if (length(eig) == 2L)
    {
        v <- x$varprop[eig] * 100
        plot(x$eigenvect[,eig[1L]], x$eigenvect[,eig[2L]],
            xlab=sprintf("Eigenvector %d (%.1f%%)", eig[1L], v[eig[1L]]),
            ylab=sprintf("Eigenvector %d (%.1f%%)", eig[2L], v[eig[2L]]),
            ...)
    } else {
        pairs(x$eigenvect[, eig],
            labels=sprintf("Eig %d\n(%.1f%%)", eig, x$varprop[eig]*100),
            gap=0.2, ...)
    }

    invisible()
}


plot.snpgdsEigMixClass <- function(x, eig=c(1L, 2L), ...)
{
    stopifnot(inherits(x, "snpgdsEigMixClass"))
    stopifnot(is.numeric(eig), length(eig) >= 2L)

    if (length(eig) == 2L)
    {
        plot(x$eigenvect[,eig[1L]], x$eigenvect[,eig[2L]],
            xlab=sprintf("Eigenvector %d", eig[1L]),
            ylab=sprintf("Eigenvector %d", eig[2L]),
            ...)
    } else {
        pairs(x$eigenvect[, eig], labels=sprintf("Eig %d", eig),
            gap=0.2, ...)
    }

    invisible()
}
