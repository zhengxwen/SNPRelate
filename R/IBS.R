#######################################################################
#
# Package name: SNPRelate
#
# Description:
#     A High-performance Computing Toolset for Relatedness and
# Principal Component Analysis of SNP Data
#
# Copyright (C) 2011 - 2016        Xiuwen Zheng
# License: GPL-3
# Email: zhengxwen@gmail.com
#


#######################################################################
# Identity-By-State (IBS) analysis
#######################################################################

#######################################################################
# Calculate the identity-by-state (IBS) matrix
#

snpgdsIBS <- function(gdsobj, sample.id=NULL, snp.id=NULL,
    autosome.only=TRUE, remove.monosnp=TRUE, maf=NaN, missing.rate=NaN,
    num.thread=1L, verbose=TRUE)
{
    # check
    ws <- .InitFile2(
        cmd="Identity-By-State (IBS) analysis on genotypes:",
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
# Calculate the identity-by-state (IBS) matrix
#

snpgdsIBSNum <- function(gdsobj, sample.id=NULL, snp.id=NULL,
    autosome.only=TRUE, remove.monosnp=TRUE, maf=NaN, missing.rate=NaN,
    num.thread=1L, verbose=TRUE)
{
    # check
    ws <- .InitFile2(
        cmd="Identity-By-State (IBS) analysis on genotypes:",
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
    method=c("IBS", "GVH", "HVG"),
    type=c("per.pair", "per.snp", "matrix", "gds.file"),
    dosage=TRUE, with.id=TRUE, output=NULL, verbose=TRUE)
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
    stopifnot(is.logical(with.id), length(with.id)==1L)
    stopifnot(is.logical(dosage), length(dosage)==1L)
    stopifnot(is.logical(verbose))

    if (type == "gds.file")
    {
        stopifnot(is.character(output) & is.vector(output))
        stopifnot(length(output) == 1L)
    } else {
        if (!is.null(output))
            stop("'output' should be NULL, if 'type' is not \"gds.file\".")
    }

    if (verbose)
    {
        cat("Working space: ", ws$n.samp, " sample", .plural(ws$n.samp),
            ", ", ws$n.snp, " SNP", .plural(ws$n.snp), "\n", sep="")
        cat("Method: ", method, "\n", sep="")
    }

    if (type == "gds.file")
    {
        ZIP <- "ZIP_RA.max"

        # create GDS file
        output <- createfn.gds(output)
        # close the file at the end
        on.exit({ closefn.gds(output) })

        # add "sample.id"
        add.gdsn(output, "sample.id", paste(sample1.id, sample2.id, sep="-"),
            compress=ZIP, closezip=TRUE)

        # SNPs
        flag <- read.gdsn(index.gdsn(gdsobj, "snp.id")) %in% ws$snp.id

        # add "snp.id"
        add.gdsn(output, "snp.id", ws$snp.id, compress=ZIP, closezip=TRUE)
        # add "snp.position"
        add.gdsn(output, "snp.position",
            read.gdsn(index.gdsn(gdsobj, "snp.position"))[flag],
            compress=ZIP, closezip=TRUE)
        # add "snp.chromosome"
        add.gdsn(output, "snp.chromosome",
            read.gdsn(index.gdsn(gdsobj, "snp.chromosome"))[flag],
            compress=ZIP, closezip=TRUE)

        # add "genotype"
        gGeno <- add.gdsn(output, "genotype", storage="bit2",
            valdim = c(length(sample1.id), 0))
        put.attr.gdsn(gGeno, "sample.order")

        # sync file
        sync.gds(output)

        if (verbose)
            cat("Output: ", output$filename, "\n", sep="")
    } else
        gGeno <- NULL

    # output
    if (with.id)
        ans <- list(sample.id = ws$sample.id, snp.id = ws$snp.id)
    else
        ans <- list()

    # call the C function
    ans$score <- .Call(gnrPairScore, match(sample1.id, ws$sample.id)-1L,
        match(sample2.id, ws$sample.id)-1L, method, type, dosage, gGeno,
        verbose)
    if (type == "per.pair")
    {
        x <- as.data.frame(ans$score)
        colnames(x) <- c("Avg", "SD", "Num")
        storage.mode(x$Num) <- "integer"
        x <- cbind(x, data.frame(
            Sample1=sample1.id, Sample2=sample2.id, stringsAsFactors=FALSE))
        ans$score <- x
    } else if (type == "per.snp")
        rownames(ans$score) <- c("Avg", "SD", "Num")

    if (type == "gds.file")
        invisible(ans)
    else
        ans
}
