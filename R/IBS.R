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
# Identity-By-State (IBS) analysis
#######################################################################

#######################################################################
# To calculate the identity-by-state (IBS) matrix for SNP genotypes
#

snpgdsIBS <- function(gdsobj, sample.id=NULL, snp.id=NULL,
    autosome.only=TRUE, remove.monosnp=TRUE, maf=NaN, missing.rate=NaN,
    num.thread=1, verbose=TRUE)
{
    # check
    ws <- .InitFile2(
        cmd="Identity-By-State (IBS) analysis on SNP genotypes:",
        gdsobj=gdsobj, sample.id=sample.id, snp.id=snp.id,
        autosome.only=autosome.only, remove.monosnp=remove.monosnp,
        maf=maf, missing.rate=missing.rate, num.thread=num.thread,
        verbose=verbose)

    # call the C function
    ibs <- .Call(gnrIBSAve, ws$num.thread, verbose)

    rv <- list(sample.id = ws$sample.id, snp.id = ws$snp.id, ibs=ibs)
    class(rv) <- "snpgdsIBSClass"
    return(rv)
}



#######################################################################
# To calculate the identity-by-state (IBS) matrix for SNP genotypes
#

snpgdsIBSNum <- function(gdsobj, sample.id=NULL, snp.id=NULL,
    autosome.only=TRUE, remove.monosnp=TRUE, maf=NaN, missing.rate=NaN,
    num.thread=1, verbose=TRUE)
{
    # check
    ws <- .InitFile2(
        cmd="Identity-By-State (IBS) analysis on SNP genotypes:",
        gdsobj=gdsobj, sample.id=sample.id, snp.id=snp.id,
        autosome.only=autosome.only, remove.monosnp=remove.monosnp,
        maf=maf, missing.rate=missing.rate, num.thread=num.thread,
        verbose=verbose)

    # call the C function
    rv <- .Call(gnrIBSNum, ws$num.thread, verbose)

    # return
    list(sample.id = ws$sample.id, snp.id = ws$snp.id,
        ibs0 = rv[[1]], ibs1 = rv[[2]], ibs2 = rv[[3]])
}



#######################################################################
# To calculate the genotype score over SNPs given by pairs of individuals
#

snpgdsPairScore <- function(gdsobj, sample1.id, sample2.id, snp.id=NULL,
    method=c("IBS", "GVH", "HVG"), type=c("per.pair", "per.snp", "matrix"),
    with.id=TRUE, verbose=TRUE)
{
    # check
    if (anyDuplicated(sample1.id) != 0)
        stop("'sample1.id' has duplicated element(s).")
    if (anyDuplicated(sample2.id) != 0)
        stop("'sample2.id' has duplicated element(s).")
    stopifnot(length(sample1.id) == length(sample2.id))

    sample.id <- union(sample1.id, sample2.id)
    ws <- .InitFile(gdsobj, sample.id, snp.id, with.id=TRUE)

    method <- match.arg(method)
    type <- match.arg(type)
    stopifnot(is.logical(with.id))
    stopifnot(is.logical(verbose))

    if (verbose)
    {
        cat("Working space: ", ws$n.samp, " sample", .plural(ws$n.samp),
            ", ", ws$n.snp, " SNP", .plural(ws$n.snp), "\n", sep="")
        cat("Method: ", method, "\n", sep="")
    }

    # output
    if (with.id)
        ans <- list(sample.id = ws$sample.id, snp.id = ws$snp.id)
    else
        ans <- list()

    # call the C function
    ans$score <- .Call(gnrPairScore, match(sample1.id, ws$sample.id)-1L,
        match(sample2.id, ws$sample.id)-1L, method, type, verbose)
    if (type == "per.pair")
        colnames(ans$score) <- c("Avg", "SD", "Num")
    else if (type == "per.snp")
        rownames(ans$score) <- c("Avg", "SD", "Num")

    ans
}
