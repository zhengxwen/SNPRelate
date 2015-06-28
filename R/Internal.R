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
# Internal Functions
#######################################################################

.plural <- function(num)
{
    if (num > 1) "s" else ""
}


#######################################################################
# Check whether SNP GDS file or not

.CheckFile <- function(gdsobj)
{
    if (!inherits(gdsobj, "gds.class"))
        stop("`gdsobj' should be a GDS file.")

    if (!inherits(gdsobj, "SNPGDSFileClass"))
    {
        if (!inherits(gdsobj, "SeqVarGDSClass"))
        {
            message("Hint: ",
                "it is suggested to call `snpgdsOpen' to open a SNP GDS file ",
                "instead of `openfn.gds'.")
        }
    }
}


#######################################################################
# Initialize the working space with sample and SNP IDs

.InitFile <- function(gdsobj, sample.id, snp.id, with.id=FALSE)
{
    # check
    .CheckFile(gdsobj)

    # samples
    if (!is.null(sample.id))
    {
        n.tmp <- length(sample.id)
        samp.tmp <- read.gdsn(index.gdsn(gdsobj, "sample.id"))
        sample.id <- samp.tmp %in% sample.id
        n.samp <- sum(sample.id)
        if (n.samp != n.tmp)
            stop("Some of sample.id do not exist!")
        if (n.samp <= 0)
            stop("No sample in the working dataset.")
        if (with.id)
            samp.tmp <- samp.tmp[sample.id]
    } else {
        if (with.id)
            samp.tmp <- read.gdsn(index.gdsn(gdsobj, "sample.id"))
    }

    # SNPs
    snp.var <- "snp.id"
    if (inherits(gdsobj, "SeqVarGDSClass"))
        snp.var <- "variant.id"
    if (!is.null(snp.id))
    {
        n.tmp <- length(snp.id)
        snp.tmp <- read.gdsn(index.gdsn(gdsobj, snp.var))
        snp.id <- snp.tmp %in% snp.id
        n.snp <- sum(snp.id)
        if (n.snp != n.tmp)
            stop("Some of snp.id do not exist!")
        if (n.snp <= 0)
            stop("No SNP in the working dataset.")
        if (with.id)
            snp.tmp <- snp.tmp[snp.id]
    } else {
        if (with.id)
            snp.tmp <- read.gdsn(index.gdsn(gdsobj, snp.var))
    }

    # set genotype working space
    if (!inherits(gdsobj, "SeqVarGDSClass"))
    {
        v <- .Call(gnrSetGenoSpace, index.gdsn(gdsobj, "genotype"),
            sample.id, snp.id)
    } else {
        v <- .Call(gnrSetSeqSpace, gdsobj, sample.id, snp.id)
    }

    # output
    if (with.id)
    {
        list(n.snp=v[1], n.samp=v[2], sample.id=samp.tmp, snp.id=snp.tmp)
    } else {
        list(n.snp=v[1], n.samp=v[2],
            samp.flag=sample.id, snp.flag=snp.id)
    }
}



#######################################################################
# Initialize the working space with more parameters

.InitFile2 <- function(cmd=NULL, gdsobj, sample.id, snp.id,
    autosome.only=TRUE, remove.monosnp=TRUE, maf=NaN, missing.rate=NaN,
    allele.freq=NULL, num.thread=1L, verbose=TRUE, verbose.work=TRUE,
    verbose.numthread=TRUE)
{
    # check
    stopifnot(is.null(cmd) | is.character(cmd))
    .CheckFile(gdsobj)
    stopifnot(is.null(sample.id) | is.vector(sample.id) | is.factor(sample.id))
    stopifnot(is.null(snp.id) | is.vector(snp.id) | is.factor(snp.id))

    stopifnot(is.logical(autosome.only) | is.numeric(autosome.only) |
        is.character(autosome.only))
    stopifnot(is.vector(autosome.only))
    stopifnot(length(autosome.only) == 1)

    stopifnot(is.logical(remove.monosnp))
    stopifnot(is.vector(remove.monosnp))
    stopifnot(length(remove.monosnp) == 1)

    stopifnot(is.numeric(maf))
    stopifnot(is.vector(maf))
    stopifnot(length(maf) == 1)

    stopifnot(is.numeric(missing.rate))
    stopifnot(is.vector(missing.rate))
    stopifnot(length(missing.rate) == 1)

    stopifnot(is.numeric(num.thread))
    stopifnot(is.vector(num.thread))
    stopifnot(length(num.thread) == 1)
    num.thread <- as.integer(num.thread)
    if (is.na(num.thread))
    {
        if (!requireNamespace("parallel", quietly=TRUE))
            stop("The 'parallel' package fails to load.")
        num.thread <- parallel::detectCores()
        if (is.na(num.thread))
            stop("parallel::detectCores fails to detect the number of cores.")
    }
    if (num.thread < 1)
        stop("`num.thread' should be a positive value or NA.")

    if (!is.null(allele.freq))
    {
        if (!is.numeric(allele.freq) || !is.vector(allele.freq))
            stop("'allele.freq' should be a numeric vector or NULL.")
    }

    stopifnot(is.logical(verbose))
    stopifnot(is.vector(verbose))
    stopifnot(length(verbose) == 1)


    # samples
    sample.ids <- read.gdsn(index.gdsn(gdsobj, "sample.id"))
    if (!is.null(sample.id))
    {
        n.tmp <- length(sample.id)
        sample.id <- sample.ids %in% sample.id
        n.samp <- sum(sample.id)
        if (n.samp != n.tmp)
            stop("Some of sample.id do not exist!")
        if (n.samp <= 0)
            stop("No sample in the working dataset.")
        sample.ids <- sample.ids[sample.id]
    }

    if (verbose & !is.null(cmd))
    {
        cat(cmd)
        cat("\n")
    }

    # SNPs
    snp.var <- "snp.id"
    chr.var <- "snp.chromosome"
    if (inherits(gdsobj, "SeqVarGDSClass"))
    {
        snp.var <- "variant.id"
        chr.var <- "chromosome"
    }

    snp.ids <- read.gdsn(index.gdsn(gdsobj, snp.var))
    if (!is.null(snp.id))
    {
        # specify 'snp.id'

        n.tmp <- length(snp.id)
        if (!is.null(allele.freq))
        {
            if (length(allele.freq) != n.tmp)
                stop("'length(allele.freq)' should be 'length(snp.id)'.")
            tmp.id <- snp.id
        }
        snp.id <- snp.ids %in% snp.id
        n.snp <- sum(snp.id)
        if (n.snp != n.tmp)
            stop("Some of snp.id do not exist!")
        if (n.snp <= 0)
            stop("No SNP in the working dataset.")

        if (!identical(autosome.only, FALSE))
        {
            nt <- index.gdsn(gdsobj, chr.var)
            if (identical(autosome.only, TRUE))
            {
                dt <- objdesp.gdsn(nt)
                if (dt$type == "String")
                {
                    snp.id <- snp.id & .Call(gnrChromParseNumeric, nt)
                } else {
                    opt <- snpgdsOption(gdsobj)
                    snp.id <- snp.id & .Call(gnrChromRangeNumeric, nt,
                        as.integer(opt$autosome.start),
                        as.integer(opt$autosome.end))
                }
            } else {
                snp.id <- snp.id & (read.gdsn(nt) == autosome.only)
                snp.id[is.na(snp.id)] <- FALSE
            }

            if (!is.null(allele.freq))
                allele.freq <- allele.freq[match(snp.ids[snp.id], tmp.id)]
            if (verbose)
            {
                if (identical(autosome.only, TRUE))
                {
                    m <- dt$dim - sum(snp.id)
                    cat("Excluding ", m, " SNP", .plural(m),
                        " on non-autosomes\n", sep="")
                } else {
                    m <- sum(snp.id)
                    cat("Keeping ", m, " SNP", .plural(m),
                        " according to chromosome ", autosome.only, "\n", sep="")
                }
            }
        } else {
            if (!is.null(allele.freq))
                allele.freq <- allele.freq[match(snp.ids[snp.id], tmp.id)]
        }
        snp.ids <- snp.ids[snp.id]

    } else {
        # not specify 'snp.id'

        if (!is.null(allele.freq))
        {
            if (length(allele.freq) != length(snp.ids))
                stop("'length(allele.freq)' should be the number of SNPs.")
        }

        if (!identical(autosome.only, FALSE))
        {
            nt <- index.gdsn(gdsobj, chr.var)
            if (identical(autosome.only, TRUE))
            {
                dt <- objdesp.gdsn(nt)
                if (dt$type == "String")
                {
                    snp.id <- .Call(gnrChromParseNumeric, nt)
                } else {
                    opt <- snpgdsOption(gdsobj)
                    snp.id <- .Call(gnrChromRangeNumeric, nt,
                        as.integer(opt$autosome.start),
                        as.integer(opt$autosome.end))
                }
            } else {
                snp.id <- (read.gdsn(nt) == autosome.only)
                snp.id[is.na(snp.id)] <- FALSE
            }

            if (!is.null(allele.freq))
                allele.freq <- allele.freq[snp.id]
            if (verbose)
            {
                if (identical(autosome.only, TRUE))
                {
                    m <- dt$dim - sum(snp.id)
                    cat("Excluding ", m, " SNP", .plural(m),
                        " on non-autosomes\n", sep="")
                } else {
                    m <- sum(snp.id)
                    cat("Keeping ", m, " SNP", .plural(m),
                        " according to chromosome ", autosome.only, "\n", sep="")
                }
            }
            snp.ids <- snp.ids[snp.id]
        }
    }

    # set genotype working space
    if (!inherits(gdsobj, "SeqVarGDSClass"))
    {
        .Call(gnrSetGenoSpace, index.gdsn(gdsobj, "genotype"),
            sample.id, snp.id)
    } else {
        .Call(gnrSetSeqSpace, gdsobj, sample.id, snp.id)
    }

    # call allele freq. and missing rates
    if (remove.monosnp || is.finite(maf) || is.finite(missing.rate))
    {
        t.maf <- maf; t.miss <- missing.rate

        if (!is.finite(maf)) maf <- -1;
        if (!is.finite(missing.rate)) missing.rate <- 2;

        maf <- as.double(maf)
        missing.rate <- as.double(missing.rate)

        # call
        if (is.null(allele.freq))
        {
            rv <- .Call(gnrSelSNP_Base, remove.monosnp, maf, missing.rate)
        } else {
            rv <- .Call(gnrSelSNP_Base_Ex, as.double(allele.freq),
                remove.monosnp, maf, missing.rate)
        }
        snp.ids <- snp.ids[rv[[2]]]
        if (!is.null(allele.freq))
            allele.freq <- allele.freq[rv[[2]]]

        # show
        if (verbose)
        {
            cat("Excluding ", rv[[1]], " SNP", .plural(rv[[1]]),
                " (monomorphic: ", remove.monosnp, ", < MAF: ", t.maf,
                ", or > missing rate: ", t.miss, ")\n", sep="")
        }
    }

    # get the dimension of SNP genotypes
    # dm[1] -- # of SNPs, dm[2] -- # of samples
    dm <- .Call(gnrGetGenoDim)

    if (verbose && verbose.work)
    {
        cat("Working space: ", dm[2], " sample", .plural(dm[2]), ", ",
            dm[1], " SNP", .plural(dm[1]), "\n", sep="")
        if (verbose.numthread)
        {
            cat("\tUsing ", num.thread, " (CPU) core",
                .plural(num.thread), "\n", sep="")
        }
    }

    # output
    list(sample.id = sample.ids, snp.id = snp.ids,
        n.snp = dm[1], n.samp = dm[2], allele.freq = allele.freq,
        num.thread = num.thread, verbose = verbose)
}



#######################################################################
# Open and close a connection,
# Please always call '.CloseConnection' after '.OpenConnText'

.LastStr <- function(s, len)
{
    substring(s, nchar(s)-len+1L, nchar(s))
}

.OpenConnBin <- function(filename)
{
    stopifnot(is.character(filename))
    con2 <- NULL

    if ((substr(filename, 1, 6) == "ftp://") |
        (substr(filename, 1, 7) == "http://"))
    {
        if (.LastStr(filename, 3) == ".gz")
        {
            con <- gzcon(url(filename, "rb"))
        } else
            con <- url(filename, "rb")
    } else {
        if (.LastStr(filename, 3) == ".gz")
        {
            con <- gzfile(filename, "rb")
        } else if (.LastStr(filename, 3) == ".xz")
        {
            con <- xzfile(filename, "rb")
        } else 
            con <- file(filename, "rb")
    }

    # open(con)
    list(filename=filename, con=con, con2=con2)
}

.OpenConnText <- function(filename, require.txtmode=FALSE)
{
    stopifnot(is.character(filename))
    con2 <- NULL

    if ((substr(filename, 1, 6) == "ftp://") |
        (substr(filename, 1, 7) == "http://"))
    {
        if (.LastStr(filename, 3) == ".gz")
        {
            con <- gzcon(url(filename, "rb"))
            if (require.txtmode)
            {
                fn <- tempfile(fileext=".tmpfile")
                write(readLines(con), file=fn, sep="\n")
                close(con)
                con <- fn
                con2 <- NULL
            }
        } else
            con <- url(filename, "rt")
    } else
        con <- file(filename, "rt")

    # open(con)
    list(filename=filename, con=con, con2=con2)
}

.CloseConnection <- function(conn)
{
    if (is.character(conn$con))
    {
        if (.LastStr(conn$con, 8) == ".tmpfile")
            unlink(conn$con, force=TRUE)
    } else if (inherits(conn$con, "connection"))
    {
        close(conn$con)
    }

    if (inherits(conn$con2, "connection"))
    {
        close(conn$con2)
    }
}



#######################################################################
# Internal R library functions
#######################################################################

.onAttach <- function(lib, pkg)
{
    # information
    s <- switch(.Call(gnrSSEFlag),
        " -- supported by Streaming SIMD Extensions 2 (SSE2)",
        " -- supported by Advanced Vector Extensions (AVX)"
    )
    packageStartupMessage("SNPRelate", s)
    TRUE
}
