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
# Linkage Disequilibrium (LD) analysis
#######################################################################

#######################################################################
# To calculate LD for a pair of SNPs
#

snpgdsLDpair <- function(snp1, snp2,
    method=c("composite", "r", "dprime", "corr"))
{
    # check
    stopifnot(is.numeric(snp1) & is.vector(snp1))
    stopifnot(is.numeric(snp2) & is.vector(snp2))
    stopifnot(length(snp1) == length(snp2))

    method <- match.arg(method)
    method <- match(method, c("composite", "r", "dprime", "corr"))

    # call
    rv <- .Call(gnrLDpair, as.integer(snp1), as.integer(snp2), method)

    # output
    if (method %in% c(2, 3))
    {
        names(rv) <- c("ld", "pA_A", "pA_B", "pB_A", "pB_B")
    } else {
        rv <- rv[1]
        names(rv) <- "ld"
    }
    rv
}



#######################################################################
# To calculate LD for each pair of SNPs in the region
#

snpgdsLDMat <- function(gdsobj, sample.id=NULL, snp.id=NULL,
    slide=250, method=c("composite", "r", "dprime", "corr"),
    num.thread=1, verbose=TRUE)
{
    # check
    ws <- .InitFile(gdsobj, sample.id=sample.id, snp.id=snp.id)

    stopifnot(is.numeric(slide))
    stopifnot(is.numeric(num.thread) & (num.thread>0))
    stopifnot(is.logical(verbose))

    # method
    method <- match.arg(method)
    method <- match(method[1], c("composite", "r", "dprime", "corr"))
    if (is.na(method))
    {
        stop("method should be one of ",
            "\"composite\", \"r\", \"dprime\" and \"corr\"")
    }

    slide <- as.integer(slide)
    if (is.na(slide)) slide <- as.integer(-1)
    if (slide > ws$n.snp) slide <- ws$n.snp

    if (verbose)
    {
        cat("Linkage Disequilibrium (LD) analysis on SNP genotypes:\n");
        cat("Working space:", ws$n.samp, "samples,", ws$n.snp, "SNPs\n");
        if (num.thread <= 1)
            cat("\tUsing", num.thread, "(CPU) core.\n")
        else
            cat("\tUsing", num.thread, "(CPU) cores.\n")
        if (slide > 0)
            cat("\tSliding window size:", slide, "\n")
    }

    # call C function
    rv <- list(sample.id = NULL, snp.id = NULL,
        LD = .Call(gnrLDMat, method, slide, as.integer(num.thread), verbose),
        slide = slide)

    # output

    rv$sample.id <- read.gdsn(index.gdsn(gdsobj, "sample.id"))
    if (!is.null(ws$samp.flag))
        rv$sample.id <- rv$sample.id[ws$samp.flag]

    rv$snp.id <- read.gdsn(index.gdsn(gdsobj, "snp.id"))
    if (!is.null(ws$snp.flag))
        rv$snp.id <- rv$snp.id[ws$snp.flag]

    return(rv)
}



#######################################################################
# To prune SNPs based on LD
#

snpgdsLDpruning <- function(gdsobj, sample.id=NULL, snp.id=NULL,
    autosome.only=TRUE, remove.monosnp=TRUE, maf=NaN, missing.rate=NaN,
    method=c("composite", "r", "dprime", "corr"),
    slide.max.bp=500000, slide.max.n=NA, ld.threshold=0.2,
    num.thread=1, verbose=TRUE)
{
    # check
    ws <- .InitFile2(
        cmd="SNP pruning based on LD:",
        gdsobj=gdsobj, sample.id=sample.id, snp.id=snp.id,
        autosome.only=autosome.only, remove.monosnp=remove.monosnp,
        maf=maf, missing.rate=missing.rate, num.thread=num.thread,
        verbose=verbose)

    stopifnot(is.na(slide.max.bp) | is.numeric(slide.max.bp))
    stopifnot(is.na(slide.max.n) | is.numeric(slide.max.n))
    stopifnot(is.numeric(ld.threshold) & is.finite(ld.threshold))
    stopifnot(is.numeric(num.thread) & (num.thread>0))
    if (num.thread > 1)
    {
        warning("The current version of 'snpgdsLDpruning' ",
            "does not support multi processes.")
    }
    stopifnot(is.logical(verbose))

    if (verbose)
    {
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

    # method
    method <- match(method[1], c("composite", "r", "dprime", "corr"))
    if (is.na(method))
    {
        stop("method should be one of ",
            "\"composite\", \"r\", \"dprime\" and \"corr\"")
    }

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
    ntotal <- 0; res <- list()
    snp.flag <- total.snp.ids %in% ws$snp.id
    samp.flag <- total.samp.ids %in% ws$sample.id

    if (is.numeric(chr))
        chrset <- setdiff(unique(chr), c(0, NA))
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
            startidx <- sample(1:n.tmp, 1)

            rv <- .Call(gnrLDpruning, as.integer(startidx-1), position[flag],
                as.integer(slide.max.bp), as.integer(slide.max.n),
                as.double(ld.threshold), method)

            # output
            L <- rep(FALSE, length(total.snp.ids))
            L[flag] <- rv
            res[[paste("chr", ch, sep="")]] <- total.snp.ids[L]
            ntotal <- ntotal + sum(rv)

            # information
            if (verbose)
            {
                ntmp <- sum(rv); ntot <- sum(chr == ch)
                cat(sprintf("Chromosome %s: %0.2f%%, %d/%d\n",
                    as.character(ch), 100*ntmp/ntot, ntmp, ntot))
            }
        }
    }

    if (verbose)
        cat(sprintf("%d SNPs are selected in total.\n", ntotal))

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
    stopifnot(length(min.dist) == 1)

    stopifnot(is.numeric(max.n.snp.perchr) & is.vector(max.n.snp.perchr))
    stopifnot(length(max.n.snp.perchr) == 1)

    stopifnot(is.logical(verbose) & is.vector(verbose))
    stopifnot(length(verbose) == 1)

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
