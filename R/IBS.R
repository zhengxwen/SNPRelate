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
