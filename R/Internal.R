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
# Internal Functions
#######################################################################

.plural <- function(num)
{
    if (num > 1) "s" else ""
}

.pretty <- function(x)
{
    prettyNum(x, big.mark=",", scientific=FALSE)
}



#######################################################################
# Check whether SNP GDS file or not

.CheckFile <- function(gdsobj)
{
    if (is.character(gdsobj))
    {
        stop(paste0("`gdsobj' is a character variable. ",
            "Please use snpgdsOpen() to open the file '", gdsobj, "'."))
    }
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
        if (n.samp <= 0L)
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
        if (n.snp <= 0L)
        {
            if (inherits(gdsobj, "SeqVarGDSClass"))
                stop("No SNV in the working dataset.")
            else
                stop("No SNP in the working dataset.")
        }
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
        list(n.snp=v[1L], n.samp=v[2L], sample.id=samp.tmp, snp.id=snp.tmp)
    } else {
        list(n.snp=v[1L], n.samp=v[2L],
            samp.flag=sample.id, snp.flag=snp.id)
    }
}



#######################################################################
# Initialize the working space with more parameters

.seldim <- function(gdsobj)
{
    eval(parse(text="
        SeqArray::seqSummary(gdsobj, 'genotype', check='none',
            verbose=FALSE)$seldim[-1L]
    "))
}

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
    stopifnot(length(autosome.only) == 1L)

    stopifnot(is.logical(remove.monosnp), length(remove.monosnp)==1L)
    stopifnot(is.numeric(maf), length(maf)==1L)
    stopifnot(is.numeric(missing.rate), length(missing.rate)==1L)

    stopifnot(is.numeric(num.thread), length(num.thread)==1L)
    num.thread <- as.integer(num.thread)
    if (is.na(num.thread))
    {
        if (!requireNamespace("parallel", quietly=TRUE))
            stop("The 'parallel' package fails to load.")
        num.thread <- parallel::detectCores()
        if (is.na(num.thread))
            stop("parallel::detectCores() fails to detect the number of cores.")
    }
    if (num.thread < 1L)
        stop("`num.thread' should be a positive value or NA.")

    if (!is.null(allele.freq))
    {
        if (!is.numeric(allele.freq) || !is.vector(allele.freq))
            stop("'allele.freq' should be a numeric vector or NULL.")
    }

    stopifnot(is.logical(verbose), length(verbose)==1L)
    if (verbose & !is.null(cmd)) cat(cmd, "\n", sep="")


    if (inherits(gdsobj, "SeqVarGDSClass"))
    {
        # whether the function 'SeqArray::seqSetFilterCond()' exists or not
        fnok <- tryCatch({
            eval(parse(text="SeqArray::seqSetFilterCond")); TRUE
        }, error=function(e) FALSE)
        if (fnok)
        {
            SeqArray::seqResetFilter(gdsobj, verbose=FALSE)
            SeqArray::seqSetFilter(gdsobj, sample.id=sample.id,
                variant.id=snp.id, verbose=FALSE)
            n <- .seldim(gdsobj)[2L]
            if (isTRUE(autosome.only))
            {
                eval(parse(text="
                    SeqArray::seqSetFilterChrom(gdsobj, is.num=TRUE,
                        intersect=TRUE, verbose=FALSE)
                "))
                m <- .seldim(gdsobj)[2L]
                if (verbose & m<n)
                {
                    cat("Excluding ", .pretty(n-m), " SNV", .plural(n-m),
                        " on non-autosomes\n", sep="")
                }
                n <- m
            }

            # call allele freq. and missing rates
            if (remove.monosnp | is.finite(maf) | is.finite(missing.rate))
            {
                mac <- if (remove.monosnp) 1L else NA_integer_
                if (verbose)
                    cat("Calculating allele counts/frequencies ...\n")
                eval(parse(text="
                    SeqArray::seqSetFilterCond(gdsobj, maf=maf, mac=mac,
                        missing.rate=missing.rate, parallel=num.thread,
                        .progress=verbose, verbose=FALSE)
                "))
                m <- .seldim(gdsobj)[2L]
                if (verbose & m<n)
                {
                    cat("Excluding ", .pretty(n-m), " SNV", .plural(n-m),
                        " (monomorphic: ", remove.monosnp, ", MAF: ", maf,
                        ", missing rate: ", missing.rate, ")\n", sep="")
                }
                n <- m
            }

            # get the dimension of genotype
            dm <- .seldim(gdsobj)
            if (verbose & verbose.work)
            {
                cat("Working space: ",
                    .pretty(dm[1L]), " sample", .plural(dm[1L]), ", ",
                    .pretty(dm[2L]), " SNV", .plural(dm[2L]), "\n", sep="")
                if (verbose.numthread)
                {
                    cat("    using ", num.thread, " (CPU) core",
                        .plural(num.thread), "\n", sep="")
                }
            }

			if (!is.null(allele.freq))
			{
				if (length(allele.freq) != dm[3L])
				{
					stop("the length of allele.freq doest not ",
						"match with the selected variants.")
				}
			}

            sel <- SeqArray::seqGetFilter(gdsobj)
            .Call(gnrSetSeqSpace, gdsobj, sel$sample.sel, sel$variant.sel)

            sampid <- SeqArray::seqGetData(gdsobj, "sample.id")
            snpid <- SeqArray::seqGetData(gdsobj, "variant.id")

            return(list(sample.id = sampid, snp.id = snpid,
                n.snp = dm[2L], n.samp = dm[1L], allele.freq = allele.freq,
                num.thread = num.thread, verbose = verbose))
        }
    }

    # samples
    sample.ids <- read.gdsn(index.gdsn(gdsobj, "sample.id"))
    if (!is.null(sample.id))
    {
        n.tmp <- length(sample.id)
        sample.id <- sample.ids %in% sample.id
        n.samp <- sum(sample.id)
        if (n.samp != n.tmp)
            stop("Some of sample.id do not exist!")
        if (n.samp <= 0L)
            stop("No sample in the working dataset.")
        sample.ids <- sample.ids[sample.id]
    }

    # SNPs
    snp.var <- "snp.id"
    chr.var <- "snp.chromosome"
    SSS <- "SNP"
    if (inherits(gdsobj, "SeqVarGDSClass"))
    {
        snp.var <- "variant.id"
        chr.var <- "chromosome"
        SSS <- "SNV"
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
        if (n.snp <= 0L)
            stop("No ", SSS, " in the working dataset.")

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
                    cat("Excluding ", .pretty(m), " ", SSS, .plural(m),
                        " (non-autosomes or non-selection)\n", sep="")
                } else {
                    m <- sum(snp.id)
                    cat("Keeping ", .pretty(m), " ", SSS, .plural(m),
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
                stop("'length(allele.freq)' should be the number of ", SSS, ".")
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
                    cat("Excluding ", .pretty(m), " ", SSS, .plural(m),
                        " on non-autosomes\n", sep="")
                } else {
                    m <- sum(snp.id)
                    cat("Keeping ", .pretty(m), " ", SSS, .plural(m),
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
    if (remove.monosnp | is.finite(maf) | is.finite(missing.rate))
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
        snp.ids <- snp.ids[rv[[2L]]]
        if (!is.null(allele.freq))
            allele.freq <- allele.freq[rv[[2L]]]

        # show
        if (verbose)
        {
            cat("Excluding ", .pretty(rv[[1L]]), " ", SSS, .plural(rv[[1L]]),
                " (monomorphic: ", remove.monosnp, ", MAF: ", t.maf,
                ", missing rate: ", t.miss, ")\n", sep="")
        }
    }

    # get the dimension of SNP genotypes
    # dm[1] -- # of SNPs, dm[2] -- # of samples
    dm <- .Call(gnrGetGenoDim)

    if (verbose & verbose.work)
    {
        cat("Working space: ", .pretty(dm[2L]), " sample", .plural(dm[2L]),
            ", ", .pretty(dm[1L]), " ", SSS, .plural(dm[1L]), "\n", sep="")
        if (verbose.numthread)
        {
            cat("    using ", num.thread, " (CPU) core", .plural(num.thread),
                "\n", sep="")
        }
    }

    # output
    list(sample.id = sample.ids, snp.id = snp.ids,
        n.snp = dm[1L], n.samp = dm[2L], allele.freq = allele.freq,
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
