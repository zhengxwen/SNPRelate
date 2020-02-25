#######################################################################
#
# Package name: SNPRelate
#
# Description:
#     A High-performance Computing Toolset for Relatedness and
# Principal Component Analysis of SNP Data
#
# Copyright (C) 2011 - 2020        Xiuwen Zheng
# License: GPL-3
# Email: zhengxwen@gmail.com
#


#######################################################################
# Linkage Disequilibrium (LD) analysis
#######################################################################

#######################################################################
# To calculate LD for a pair of SNPs
#

snpgdsLDpair <- function(snp1, snp2,
    method=c("composite", "r", "dprime", "corr"))
{
    # check
    stopifnot(is.numeric(snp1), is.vector(snp1))
    stopifnot(is.numeric(snp2), is.vector(snp2))
    stopifnot(length(snp1)==length(snp2))

    method <- match.arg(method)
    method <- match(method, c("composite", "r", "dprime", "corr"))

    # call
    rv <- .Call(gnrLDpair, as.integer(snp1), as.integer(snp2), method)

    # output
    if (method %in% c(2L, 3L))
    {
        names(rv) <- c("ld", "pA_A", "pA_B", "pB_A", "pB_B")
    } else {
        rv <- rv[1L]
        names(rv) <- "ld"
    }
    rv
}



#######################################################################
# To calculate LD for each pair of SNPs in the region
#

snpgdsLDMat <- function(gdsobj, sample.id=NULL, snp.id=NULL, slide=250L,
    method=c("composite", "r", "dprime", "corr", "cov"), mat.trim=FALSE,
    num.thread=1L, with.id=TRUE, verbose=TRUE)
{
    # check
    ws <- .InitFile(gdsobj, sample.id=sample.id, snp.id=snp.id, with.id=with.id)

    stopifnot(is.numeric(slide), length(slide)==1L)
    stopifnot(is.numeric(num.thread), num.thread > 0L)
    stopifnot(is.logical(mat.trim), length(mat.trim)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L)

    # method
    method <- match.arg(method)
    method <- match(method, c("composite", "r", "dprime", "corr", "cov"))

    slide <- as.integer(slide)
    if (is.na(slide)) slide <- -1L
    if (slide > ws$n.snp) slide <- ws$n.snp

    if (verbose)
    {
        cat("Linkage Disequilibrium (LD) estimation on genotypes:\n")
        .cat("    # of samples: ", .pretty(ws$n.samp))
        .cat("    # of SNPs: ", .pretty(ws$n.snp))
        .cat("    using ", num.thread, " thread", .plural(num.thread))
        if (slide > 0L)
            .cat("    sliding window size: ", slide)
        .cat("    method: ",
            c("composite", "R", "D'", "correlation", "covariance")[method])
    }

    # call C function
    m <- .Call(gnrLDMat, method, slide, mat.trim, num.thread, verbose)

    # output
    if (with.id)
        m <- list(sample.id=ws$sample.id, snp.id=ws$snp.id, LD=m, slide=slide)
    m
}



#######################################################################
# To prune SNPs based on LD
#

snpgdsLDpruning <- function(gdsobj, sample.id=NULL, snp.id=NULL,
    autosome.only=TRUE, remove.monosnp=TRUE, maf=NaN, missing.rate=NaN,
    method=c("composite", "r", "dprime", "corr"),
    slide.max.bp=500000L, slide.max.n=NA, ld.threshold=0.2,
    start.pos=c("random", "first", "last"), num.thread=1L, verbose=TRUE)
{
    # check
    ws <- .InitFile2(
        cmd=paste(ifelse(inherits(gdsobj, "SeqVarGDSClass"), "SNV", "SNP"),
            "pruning based on LD:"),
        gdsobj=gdsobj, sample.id=sample.id, snp.id=snp.id,
        autosome.only=autosome.only, remove.monosnp=remove.monosnp,
        maf=maf, missing.rate=missing.rate, num.thread=num.thread,
        verbose=verbose)

    stopifnot(is.na(slide.max.bp) | is.numeric(slide.max.bp))
    stopifnot(is.na(slide.max.n) | is.numeric(slide.max.n))
    stopifnot(is.numeric(ld.threshold), is.finite(ld.threshold),
        length(ld.threshold)==1L)
    stopifnot(is.numeric(num.thread), num.thread > 0L)
    if (num.thread > 1L)
    {
        warning("The current version of 'snpgdsLDpruning()' ",
            "does not support multi-threading.", immediate.=TRUE)
    }
    stopifnot(is.logical(verbose), length(verbose)==1L)

    start.pos <- match.arg(start.pos)

    if (verbose)
    {
        bp <- slide.max.bp; mn <- slide.max.n
        if (!is.finite(bp)) bp <- Inf
        if (!is.finite(mn)) mn <- Inf
        cat("    sliding window:", prettyNum(bp, ",", scientific=FALSE),
            "basepairs,", prettyNum(mn, ",", scientific=FALSE), "SNPs\n")
        cat(sprintf("    |LD| threshold: %g\n", ld.threshold))
    }

    if (!is.finite(slide.max.bp))
        slide.max.bp <- .Machine$double.xmax
    if (!is.finite(slide.max.n))
        slide.max.n <- .Machine$integer.max

    # method
    method <- match(method[1L], c("composite", "r", "dprime", "corr"))
    if (is.na(method))
    {
        stop("method should be one of ",
            "\"composite\", \"r\", \"dprime\" and \"corr\"")
    }
    if (verbose)
        .cat("    method: ", c("composite", "R", "D'", "correlation")[method])

    if (!inherits(gdsobj, "SeqVarGDSClass"))
    {
        total.snp.ids <- read.gdsn(index.gdsn(gdsobj, "snp.id"))
        total.samp.ids <- read.gdsn(index.gdsn(gdsobj, "sample.id"))
        chr <- read.gdsn(index.gdsn(gdsobj, "snp.chromosome"))
        position <- read.gdsn(index.gdsn(gdsobj, "snp.position"))
    } else {
        total.snp.ids <- read.gdsn(index.gdsn(gdsobj, "variant.id"))
        total.samp.ids <- read.gdsn(index.gdsn(gdsobj, "sample.id"))
        chr <- read.gdsn(index.gdsn(gdsobj, "chromosome"))
        position <- read.gdsn(index.gdsn(gdsobj, "position"))
    }

    # for-loop each chromosome
    ntotal <- 0L
    res <- list()
    snp.flag <- total.snp.ids %in% ws$snp.id
    samp.flag <- total.samp.ids %in% ws$sample.id

    if (is.numeric(chr))
        chrset <- setdiff(unique(chr), c(0L, NA))
    else if (is.character(chr))
        chrset <- setdiff(unique(chr), c("", NA))
    else
        stop("Unknown format of 'snp.chromosome'!")

    for (ch in chrset)
    {
        flag <- snp.flag & (chr == ch)
        flag[is.na(flag)] <- FALSE
        n.tmp <- sum(flag)
        if (n.tmp > 0)
        {
            # set genotype working space
            if (!inherits(gdsobj, "SeqVarGDSClass"))
            {
                .Call(gnrSetGenoSpace, index.gdsn(gdsobj, "genotype"),
                    samp.flag, flag)
            } else {
                .Call(gnrSetSeqSpace, gdsobj, samp.flag, flag)
            }

            # call LD prune for this chromosome
            startidx <- switch(start.pos,
                random = sample.int(n.tmp, 1L),
                first  = 1L,
                last   = n.tmp)
            rv <- .Call(gnrLDpruning, startidx, position[flag],
                slide.max.bp, slide.max.n, ld.threshold, method)

            # output
            L <- rep(FALSE, length(total.snp.ids))
            L[flag] <- rv
            res[[paste("chr", ch, sep="")]] <- total.snp.ids[L]
            ntotal <- ntotal + sum(rv)

            # information
            if (verbose)
            {
                ntmp <- sum(rv); ntot <- sum(chr == ch)
                cat(sprintf("Chromosome %s: %0.2f%%, %s/%s\n",
                    as.character(ch), 100*ntmp/ntot,
                    prettyNum(ntmp, ",", scientific=FALSE),
                    prettyNum(ntot, ",", scientific=FALSE)))
            }
        }
    }

    if (verbose)
    {
        cat(prettyNum(ntotal, ",", scientific=FALSE),
            "markers are selected in total.\n")
    }

    # return
    return(res)
}



#######################################################################
# Randomly selects SNPs for which each pair is at least as far apart
#    as the specified basepair distance
#

snpgdsApartSelection <- function(chromosome, position, min.dist=100000,
    max.n.snp.perchr=-1, verbose=TRUE) 
{
    # check
    stopifnot(is.vector(chromosome))
    stopifnot(is.vector(position) & is.numeric(position))
    if (length(chromosome) != length(position))
        stop("The lengths of 'chomosome' and 'position' do not match.")

    stopifnot(is.numeric(min.dist) & is.vector(min.dist))
    stopifnot(length(min.dist)==1L)

    stopifnot(is.numeric(max.n.snp.perchr) & is.vector(max.n.snp.perchr))
    stopifnot(length(max.n.snp.perchr)==1L)

    stopifnot(is.logical(verbose) & is.vector(verbose))
    stopifnot(length(verbose)==1L)

    # chromosome codes
    chrlist <- unique(chromosome)

    # initialize the result
    rv <- rep(FALSE, length(chromosome))

    # for-loop
    for (chr in chrlist)
    {
        b <- (chromosome==chr)
        b[is.na(b)] <- FALSE
        pos <- position[b]
        sel <- seq_len(length(pos))
        flag <- rep(FALSE, length(pos))

        iter.num <- 0
        repeat {
            if ((length(sel)==0) | (iter.num==max.n.snp.perchr))
                break
            iter.num <- iter.num + 1
            p.i <- sel[sample.int(length(sel), 1)]
            flag[p.i] <- TRUE
            sel <- sel[abs(pos[sel] - pos[p.i]) >= min.dist]
            if (verbose & (iter.num %% 5000 == 0))
            {
                cat(date(), "\tChromosome ", chr, ", # of SNPs: ",
                    iter.num, "\n", sep="")
            }
        }
        if (verbose)
        {
            cat(date(), "\tChromosome ", chr, ", # of SNPs: ",
                    sum(flag), "\n", sep="")
        }
        rv[b] <- flag
    }
    if (verbose) 
        cat("Total # of SNPs selected:", sum(rv), "\n", sep="")

    return(rv)
}
