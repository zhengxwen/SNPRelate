#######################################################################
#
# Package name: SNPRelate
#
# Description:
#     A High-performance Computing Toolset for Relatedness and
# Principal Component Analysis of SNP Data
#
# Copyright (C) 2011 - 2018        Xiuwen Zheng
# License: GPL-3
#


#######################################################################
# Default human chromosome coding (from PLINK format):
#          autosome                        -> 1 .. 22
#     X    X chromosome                    -> 23
#     XY   Pseudo-autosomal region of X    -> 24
#     Y    Y chromosome                    -> 25
#     MT   Mitochondrial                   -> 26
#######################################################################


#######################################################################
# File Functions
#######################################################################

#######################################################################
# Open a SNP GDS file
#

snpgdsOpen <- function(filename, readonly=TRUE, allow.duplicate=FALSE,
    allow.fork=TRUE)
{
    ## open the GDS file
    ans <- openfn.gds(filename, readonly,
        allow.duplicate=allow.duplicate, allow.fork=allow.fork)
    on.exit({ closefn.gds(ans) })

    ##  error information
    err <- "\nInvalid SNP GDS file: "

    ##  FileFormat
    at <- get.attr.gdsn(ans$root)
    if ("FileFormat" %in% names(at))
    {
        # it does not throw any warning or error if FileFormat does not exist,
        # but it is encouraged to add this attribute
        if (identical(at$FileFormat, "SEQ_ARRAY"))
        {
            stop(err,
            "please open the file using 'seqOpen()' in the SeqArray package ",
            "since 'FileFormat=SEQ_ARRAY', or run 'SeqArray::seqSNP2GDS()' ",
            "to convert the file to GWAS SNP GDS format.")
        }
        if (!identical(at$FileFormat, "SNP_ARRAY") &&
            !identical(at$FileFormat, "IMPUTED_DOSAGE"))
        {
            stop(err, "'FileFormat' should be 'SNP_ARRAY' or 'IMPUTED_DOSAGE'.")
        }
    }

    ##  check the validation
    n <- index.gdsn(ans, "sample.id", silent=TRUE)
    if (!is.null(n))
    {
        n.samp <- objdesp.gdsn(n)$dim
        if (length(n.samp) != 1)
            stop(err, "invalid dimension of 'sample.id'.")
    } else
        stop("Invalid SNP GDS file: no 'sample.id' variable.")

    n <- index.gdsn(ans, "snp.id", silent=TRUE)
    if (!is.null(n))
    {
        n.snp <- objdesp.gdsn(n)$dim
        if (length(n.snp) != 1)
            stop(err, "invalid dimension of 'snp.id'.")
    } else
        stop(err, "no 'snp.id' variable.")

    n <- index.gdsn(ans, "snp.rs.id", silent=TRUE)
    if (!is.null(n))
    {
        if (!identical(n.snp, objdesp.gdsn(n)$dim))
        {
            stop(err,
                "inconsistent dimension between 'snp.id' and 'snp.rs.id'.")
        }
    }

    n <- index.gdsn(ans, "snp.position", silent=TRUE)
    if (!is.null(n))
    {
        if (!identical(n.snp, objdesp.gdsn(n)$dim))
        {
            stop(err,
                "inconsistent dimension between 'snp.id' and 'snp.position'.")
        }
    } else
        stop(err, "no 'snp.position' variable.")

    n <- index.gdsn(ans, "snp.chromosome", silent=TRUE)
    if (!is.null(n))
    {
        if (!identical(n.snp, objdesp.gdsn(n)$dim))
        {
            stop(err,
            "inconsistent dimension between 'snp.id' and 'snp.chromosome'.")
        }
    } else
        stop(err, "no 'snp.chromosome' variable.")

    n <- index.gdsn(ans, "snp.allele", silent=TRUE)
    if (!is.null(n))
    {
        if (!identical(n.snp, objdesp.gdsn(n)$dim))
        {
            stop(err,
                "inconsistent dimension between 'snp.id' and 'snp.allele'.")
        }
    }

    n <- index.gdsn(ans, "genotype", silent=TRUE)
    if (!is.null(n))
    {
        dm <- objdesp.gdsn(n)$dim
        if (length(dm) != 2L)
            stop(err, "'genotype' should be a matrix.")
        snpfirstdim <- TRUE
        rd <- names(get.attr.gdsn(n))
        if ("snp.order" %in% rd) snpfirstdim <- TRUE
        if ("sample.order" %in% rd) snpfirstdim <- FALSE
        if (snpfirstdim)
        {
            if ((dm[1L]!=n.snp) || (dm[2L]!=n.samp))
                stop(err, "invalid dimension of 'genotype'.")
        } else {
            if ((dm[1L]!=n.samp) || (dm[2L]!=n.snp))
                stop(err, "invalid dimension of 'genotype'.")
        }
    } else
        stop(err, "no 'genotype' variable!")

    ## output
    on.exit()
    class(ans) <- c("SNPGDSFileClass", class(ans))
    ans
}



#######################################################################
# Close the SNP GDS file
#

snpgdsClose <- function(gdsobj)
{
    if (inherits(gdsobj, "SNPGDSFileClass"))
    {
        closefn.gds(gdsobj)
    } else if (inherits(gdsobj, "gds.class"))
    {
        stop("'snpgdsClose' is only used to ",
            "close the file opened by 'snpgdsOpen'.")
    } else {
        stop("It should be a 'SNPGDSFileClass' object.")
    }
}




#######################################################################
# Summary Descriptive Statistics
#######################################################################

#######################################################################
# To calculate the missing rate and allele frequency for each SNP
#

snpgdsSNPRateFreq <- function(gdsobj, sample.id=NULL, snp.id=NULL,
    with.id=FALSE, with.sample.id=FALSE, with.snp.id=FALSE)
{
    # check
    ws <- .InitFile(gdsobj, sample.id=sample.id, snp.id=snp.id)
    stopifnot(is.logical(with.id), length(with.id)==1L)
    stopifnot(is.logical(with.sample.id), length(with.sample.id)==1L)
    stopifnot(is.logical(with.snp.id), length(with.snp.id)==1L)
    if (isTRUE(with.id))
    {
        with.snp.id <- TRUE
        with.sample.id <- TRUE
    }

    rv <- list()
    if (isTRUE(with.sample.id))
    {
        rv$sample.id <- read.gdsn(index.gdsn(gdsobj, "sample.id"))
        if (!is.null(ws$samp.flag))
            rv$sample.id <- rv$sample.id[ws$samp.flag]
    }
    if (isTRUE(with.snp.id))
    {
        rv$snp.id <- read.gdsn(index.gdsn(gdsobj, "snp.id"))
        if (!is.null(ws$snp.flag))
            rv$snp.id <- rv$snp.id[ws$snp.flag]
    }

    # call allele freq. and missing rates
    v <- .Call(gnrSNPRateFreq)
    names(v) <- c("AlleleFreq", "MinorFreq", "MissingRate")
    rv <- c(rv, v)

    rv
}



#######################################################################
# To calculate the missing rate for each sample
#

snpgdsSampMissRate <- function(gdsobj, sample.id=NULL, snp.id=NULL,
    with.id=FALSE)
{
    # check
    ws <- .InitFile(gdsobj, sample.id=sample.id, snp.id=snp.id)

    # call sample missing rates
    rv <- .Call(gnrSampFreq)

    if (with.id)
    {
        samp.id <- read.gdsn(index.gdsn(gdsobj, "sample.id"))
        if (!is.null(ws$samp.flag))
            samp.id <- samp.id[ws$samp.flag]
        names(rv) <- samp.id
    }

    rv
}


#######################################################################
# To calculate the p-value for the test of Hardy-Weinberg Equilibrium
#

snpgdsHWE <- function(gdsobj, sample.id=NULL, snp.id=NULL,
    with.id=FALSE)
{
    # check
    ws <- .InitFile(gdsobj, sample.id=sample.id, snp.id=snp.id)
    stopifnot(is.logical(with.id))

    # call C function
    rv <- .Call(gnrHWE)

    if (with.id)
    {
        rv <- list(pvalue = rv)

        rv$sample.id <- read.gdsn(index.gdsn(gdsobj, "sample.id"))
        if (!is.null(ws$samp.flag))
            rv$sample.id <- rv$sample.id[ws$samp.flag]

        rv$snp.id <- read.gdsn(index.gdsn(gdsobj, "snp.id"))
        if (!is.null(ws$snp.flag))
            rv$snp.id <- rv$snp.id[ws$snp.flag]
    }

    rv
}


#######################################################################
# Return a list of candidate SNPs
#

snpgdsSelectSNP <- function(gdsobj, sample.id=NULL, snp.id=NULL,
    autosome.only=TRUE, remove.monosnp=TRUE, maf=NaN, missing.rate=NaN,
    verbose=TRUE)
{
    # check and initialize ...
    ws <- .InitFile2(cmd=NULL,
        gdsobj=gdsobj, sample.id=sample.id, snp.id=snp.id,
        autosome.only=autosome.only, remove.monosnp=remove.monosnp,
        maf=maf, missing.rate=missing.rate, num.thread=1L,
        verbose=verbose, verbose.work=FALSE)

    # output
    ws$snp.id
}




#######################################################################
# Individual Inbreeding Coefficients
#######################################################################

#######################################################################
# To calculate individual inbreeding coefficient
#

snpgdsIndInbCoef <- function(x, p, method=c("mom.weir", "mom.visscher", "mle"),
    reltol=.Machine$double.eps^0.75)
{
    # check
    stopifnot(is.vector(x) & is.numeric(x))
    stopifnot(is.vector(p) & is.numeric(p))
    stopifnot(length(x) == length(p))

    method <- match.arg(method)
    x[!(x %in% c(0,1,2))] <- NA

    if (method == "mom.weir")
    {
        num <- x*x - (1+2*p)*x + 2*p*p
        den <- 2*p*(1-p)
        flag <- is.finite(num) & is.finite(den)
        rv <- sum(num[flag]) / sum(den[flag])
    } else if (method == "mom.visscher")
    {
        d <- (x*x - (1+2*p)*x + 2*p*p) / (2*p*(1-p))
        rv <- mean(d[is.finite(d)])
    } else if (method == "mle")
    {
        rv <- .Call(gnrIndInbCoef, x, p, reltol)
    } else {
        rv <- invisible()
    }

    return(rv)
}



#######################################################################
# To calculate individual inbreeding coefficients
#

snpgdsIndInb <- function(gdsobj, sample.id=NULL, snp.id=NULL,
    autosome.only=TRUE, remove.monosnp=TRUE, maf=NaN, missing.rate=NaN,
    method=c("mom.weir", "mom.visscher", "mle", "gcta1", "gcta2", "gcta3"),
    allele.freq=NULL, out.num.iter=TRUE, reltol=.Machine$double.eps^0.75,
    verbose=TRUE)
{
    # check
    ws <- .InitFile2(
        cmd="Estimating individual inbreeding coefficients:",
        gdsobj=gdsobj, sample.id=sample.id, snp.id=snp.id,
        autosome.only=autosome.only, remove.monosnp=remove.monosnp,
        maf=maf, missing.rate=missing.rate, allele.freq=allele.freq,
        verbose=verbose)
    method <- match.arg(method)
    stopifnot(is.logical(out.num.iter), length(out.num.iter)==1L)
    stopifnot(is.numeric(reltol), length(reltol)==1L)

    # call allele frequency
    if (is.null(allele.freq) & method=="mle")
        allele.freq <- .Call(gnrSNPFreq)

    # call individual inbreeding coefficients
    r <- .Call(gnrIndInb, allele.freq, method, reltol, out.num.iter, verbose)

    # output
    rv <- list(sample.id=ws$sample.id, snp.id=ws$snp.id, inbreeding=r[[1L]])
    if (!is.null(r[[2L]]))
        rv$out.num.iter <- r[[2L]]
    return(rv)
}



#######################################################################
# To perform hierarchical cluster analysis
#

snpgdsHCluster <- function(dist, sample.id=NULL, need.mat=TRUE, hang=0.25)
{
    # check
    stopifnot(is.matrix(dist) | inherits(dist, "snpgdsDissClass") |
        inherits(dist, "snpgdsIBSClass"))

    if (inherits(dist, "snpgdsDissClass"))
    {
        sample.id <- dist$sample.id
        dist <- dist$diss
    }
    if (inherits(dist, "snpgdsIBSClass"))
    {
        sample.id <- dist$sample.id
        dist <- 1 - dist$ibs
    }

    if (is.null(sample.id))
    {
        stopifnot(nrow(dist) == ncol(dist))
        if (is.null(colnames(dist)) | is.null(rownames(dist)))
            stop("Please specify 'sample.id'.")
        sample.id <- colnames(dist)
    } else {
        stopifnot(nrow(dist) == length(sample.id))
        stopifnot(ncol(dist) == length(sample.id))
        colnames(dist) <- sample.id
        rownames(dist) <- sample.id
    }

    # run
    hc <- hclust(as.dist(dist), method="average")
    obj <- as.dendrogram(hc, hang=hang)

    rv <- list(sample.id = sample.id, hclust = hc, dendrogram = obj)
    if (need.mat) rv$dist <- dist
    class(rv) <- "snpgdsHCClass"
    rv
}



#######################################################################
# To determine sub groups of individuals
#

snpgdsCutTree <- function(hc, z.threshold=15, outlier.n=5, n.perm=5000,
    samp.group=NULL, col.outlier="red", col.list=NULL, pch.outlier=4,
    pch.list=NULL, label.H=FALSE, label.Z=TRUE, verbose=TRUE)
{
    # check
    stopifnot(inherits(hc, "snpgdsHCClass"))
    stopifnot(is.finite(z.threshold))
    stopifnot(is.numeric(n.perm))
    stopifnot(is.logical(label.H))
    stopifnot(is.logical(label.Z))
    stopifnot(is.logical(verbose))
    stopifnot(n.perm >= 50)

    if (is.null(hc$dist))
        stop("`hc' should have a matrix of dissimilarity.")

    auto.cluster <- is.null(samp.group)
    if (verbose)
    {
        if (auto.cluster)
        {
            cat("Determine groups by permutation",
                sprintf("(Z threshold: %g, outlier threshold: %d):\n",
                z.threshold, outlier.n))
        }
    }


    # result
    ans <- list(sample.id = hc$sample.id)
    ans$z.threshold <- z.threshold
    ans$outlier.n <- outlier.n
    ans$samp.order <- hc$hclust$order

    if (!auto.cluster)
    {
        stopifnot(is.factor(samp.group))
        stopifnot(length(samp.group) == length(hc$sample.id))
        merge <- NULL

    } else {
        # determine the number of groups

        # check
        stopifnot(!is.null(hc$dist))
        stopifnot(!is.null(hc$hclust$merge))

        n <- dim(hc$dist)[1]
        rv <- .C(gnrDistPerm, n, as.double(hc$dist),
            as.integer(hc$hclust$merge), as.integer(n.perm),
            as.double(z.threshold),
            OutZ = double(n-1), OutN1 = integer(n-1), OutN2 = integer(n-1),
            OutGrp = integer(n), err=integer(1), NAOK=TRUE)
        if (rv$err != 0) stop(snpgdsErrMsg())

        merge <- data.frame(z=rv$OutZ, n1=rv$OutN1, n2=rv$OutN2)

        if (is.finite(outlier.n))
        {
            tab <- table(rv$OutGrp)
            x <- as.integer(names(tab)[tab <= outlier.n])
            flag <- rv$OutGrp %in% x

            samp.group <- sprintf("G%03d", rv$OutGrp)
            samp.group[flag] <- sprintf("Outlier%03d", rv$OutGrp[flag])
            samp.group <- as.factor(samp.group)

            n.g <- length(tab) - length(x)
            n.o <- length(x)
            nm <- character()

            if (n.g > 0) nm <- c(nm, sprintf("G%03d", 1:n.g))
            if (n.o > 0) nm <- c(nm, sprintf("Outlier%03d", 1:n.o))
            levels(samp.group) <- nm

        } else {
            # not detect outliers

            samp.group <- as.factor(sprintf("G%03d", rv$OutGrp))
            n <- levels(samp.group)
            n <- sprintf("G%03d", 1:length(n))
            levels(samp.group) <- n     
        }
    }


    # create a new dendrogram
    add.point <- function(n, sample.id, subgroup)
    {
        if(is.leaf(n))
        {
            idx <- match(attr(n, "label"), sample.id)
            if (as.character(subgroup[idx]) == "outlier")
            {
                attr(n, "nodePar") <- list(pch=pch.outlier, col=col.outlier)
            } else {
                idx <- as.integer(subgroup[idx])
                attr(n, "nodePar") <- list(pch=pch.list[idx],
                    col=col.list[idx])
            }
        } else 
        {
            if (!is.null(merge))
            {
                if (label.H | label.Z)
                {
                    x1 <- attr(n[[1]], "members")
                    x2 <- attr(n[[2]], "members")
                    w1 <- (x1 == merge$n1) & (x2 == merge$n2)
                    w2 <- (x2 == merge$n1) & (x1 == merge$n2)
                    w <- which(w1 | w2)
                    if (length(w) > 0)
                    {
                        if (merge$z[w[1]] >= z.threshold)
                        {
                            if (label.H)
                            {
                                if (label.Z)
                                {
                                    attr(n, "edgetext") <-
                                        sprintf("H: %0.2f, Z: %0.1f",
                                        attr(n, "height"), merge$z[w[1]])
                                } else {
                                    attr(n, "edgetext") <- sprintf("H: %0.2f",
                                        attr(n, "height"))
                                }
                            } else {
                                if (label.Z)
                                {
                                    attr(n, "edgetext") <-
                                        sprintf("Z: %0.1f", merge$z[w[1]])
                                }
                            }
                        }
                    }
                }
            }
        }
        n
    }

    ans$samp.group <- samp.group

    n <- nlevels(samp.group); grps <- levels(samp.group)
    dmat <- matrix(0.0, nrow=n, ncol=n)
    for (i in 1:n)
    {
        m <- hc$dist[grps[i] == samp.group, grps[i] == samp.group]
        ms <- sum(grps[i] == samp.group)
        mflag <- matrix(TRUE, nrow=ms, ncol=ms)
        diag(mflag) <- FALSE
        dmat[i, i] <- mean(m[mflag], na.rm=TRUE)
        
        if (i < n)
        {
            for (j in (i+1):n)
            {
                dmat[i, j] <- dmat[j, i] <-
                    mean(hc$dist[grps[i] == samp.group, grps[j] == samp.group],
                    na.rm=TRUE)
            }
        }
    }
    colnames(dmat) <- rownames(dmat) <- grps
    ans$dmat <- dmat

    if (is.null(pch.list))
    {
        pch.list <- rep(20, n)
    } else {
        pch.list <- rep(pch.list, (n %/% length(pch.list))+1)
    }
    if (is.null(col.list))
    {
        col.list <- 1:n
    } else {
        col.list <- rep(col.list, (n %/% length(col.list))+1)
    }

    ans$dendrogram <- dendrapply(hc$dendrogram, add.point,
        sample.id=hc$sample.id, subgroup=samp.group)
    ans$merge <- merge

    if (!is.null(merge))
    {
        cluster <- samp.group[hc$hclust$order]
        ans$clust.count <- table(cluster)[unique(cluster)]
    } else {
        ans$clust.count <- NULL
    }

    if (verbose)
        cat(sprintf("Create %d groups.\n", n))

    ans
}




#######################################################################
# SNP functions
#######################################################################

#######################################################################
# To get a list of SNP information including snp id, chr, pos, allele
#   and allele frequency
#

snpgdsSNPList <- function(gdsobj, sample.id=NULL)
{
    # check
    .CheckFile(gdsobj)

    # snp id
    snp.id <- read.gdsn(index.gdsn(gdsobj, "snp.id"))
    # chromosome
    chromosome <- read.gdsn(index.gdsn(gdsobj, "snp.chromosome"))
    # position
    position <- read.gdsn(index.gdsn(gdsobj, "snp.position"))
    # allele
    allele <- read.gdsn(index.gdsn(gdsobj, "snp.allele"))
    # allele frequency
    afreq <- snpgdsSNPRateFreq(gdsobj, sample.id=sample.id)$AlleleFreq

    # output
    rv <- data.frame(snp.id=snp.id, chromosome=chromosome, position=position,
        allele=allele, afreq=afreq, stringsAsFactors=FALSE)
    class(rv) <- c("snpgdsSNPListClass", "data.frame")
    rv
}



#######################################################################
# To get a common list of SNPs from SNP objects, and return
#     snp alleles from the first snp object
#

snpgdsSNPListIntersect <- function(snplist1, snplist2, ...,
    method=c("position", "exact"), na.rm=TRUE, same.strand=FALSE, verbose=TRUE)
{
    # check
    snplst <- list(s1=snplist1, s2=snplist2, ...)
    for (i in seq_along(snplst))
    {
        if (!inherits(snplst[[i]], "snpgdsSNPListClass"))
            stop("All objects should be 'snpgdsSNPListClass'.")
    }
    method <- match.arg(method)
    stopifnot(is.logical(na.rm), length(na.rm)==1L)
    stopifnot(is.logical(same.strand), length(same.strand)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L)

    lst <- vector("list", length(snplst))
    if (method == "position")
    {
        for (i in seq_along(snplst))
        {
            lst[[i]] <- with(snplst[[i]], paste0(chromosome, ":", position))
        }
    } else if (method == "exact")
    {
        for (i in seq_along(snplst))
        {
            lst[[i]] <- with(snplst[[i]], paste(snp.id, chromosome,
                position, allele, sep=":"))
        }
    }

    # get the common
    snp <- Reduce(intersect, lst)

    # results
    rv <- vector("list", 0L)
    for (i in seq_along(lst))
        rv[[paste0("idx", i)]] <- match(snp, lst[[i]])

    # match by position
    if (method == "position")
    {
        ii <- rv[["idx1"]]
        al1 <- snplst[[1L]]$allele[ii]
        af1 <- snplst[[1L]]$afreq[ii]
        for (i in 2:length(lst))
        {
            ii <- rv[[paste0("idx", i)]]
            rv[[paste0("flag", i)]] <- .Call(gnrAlleleStrand,
                al1, snplst[[i]]$allele[ii],
                af1, snplst[[i]]$afreq[ii], same.strand)
        }
        # remove mismatched alleles
        if (isTRUE(na.rm))
        {
            x <- rep(TRUE, length(al1))
            for (i in 2:length(lst))
                x <- x & !is.na(rv[[paste0("flag", i)]])
            if (any(!x))
            {
                for (i in seq_along(lst))
                    rv[[paste0("idx", i)]] <- rv[[paste0("idx", i)]][x]
                for (i in 2:length(lst))
                    rv[[paste0("flag", i)]] <- rv[[paste0("flag", i)]][x]
            }
        }
    }

    rv
}



#######################################################################
# GDS file management
#######################################################################

#######################################################################
# To summarize the gds file
#

snpgdsSummary <- function(gds, show=TRUE)
{
    # check
    if (is.character(gds))
    {
        gds <- snpgdsOpen(gds, allow.duplicate=TRUE)
        on.exit({ snpgdsClose(gds) })
    } else {
        .CheckFile(gds)
    }

    #
    # checking ...
    #

    warn.flag <- FALSE

    # check sample id
    dm <- objdesp.gdsn(index.gdsn(gds, "sample.id"))$dim
    if (length(dm) != 1)
    {
        print(gds)
        stop("Invalid dimension of 'sample.id'.")
    }
    samp.id <- read.gdsn(index.gdsn(gds, "sample.id"))
    if (any(duplicated(samp.id)))
    {
        warning("'sample.id' is not unique!")
        samp.id <- unique(samp.id)
        warn.flag <- TRUE
    }
    n.samp <- dm[1]

    # check snp id
    dm <- objdesp.gdsn(index.gdsn(gds, "snp.id"))$dim
    if (length(dm) != 1)
    {
        print(gds)
        stop("Invalid dimension of 'snp.id'.")
    }
    n.snp <- dm[1]
    snp.id <- read.gdsn(index.gdsn(gds, "snp.id"))
    snp.flag <- !duplicated(snp.id)
    if (!all(snp.flag))
    {
        warning("'snp.id' is not unique!")
        warn.flag <- TRUE
    }

    # check snp position
    dm <- objdesp.gdsn(index.gdsn(gds, "snp.position"))$dim
    if ((length(dm) != 1) | (dm[1] != n.snp))
    {
        print(gds)
        stop("Invalid dimension of 'snp.position'.")
    }
    snp.pos <- read.gdsn(index.gdsn(gds, "snp.position"))
    snp.pos[!is.finite(snp.pos)] <- -1
    if (any(snp.pos <= 0))
    {
        message("Some values of 'snp.position' are invalid (should be > 0)!")
        warn.flag <- TRUE
        snp.flag <- snp.flag & (snp.pos > 0)
    }

    # check snp chromosome
    dm <- objdesp.gdsn(index.gdsn(gds, "snp.chromosome"))$dim
    if ((length(dm) != 1) | (dm[1] != n.snp))
    {
        print(gds)
        stop("Invalid dimension of 'snp.chromosome'.")
    }
    snp.chr <- read.gdsn(index.gdsn(gds, "snp.chromosome"))
    if (is.numeric(snp.chr))
    {
        snp.chr[!is.finite(snp.chr)] <- -1
        flag <- (snp.chr >= 1)
        if (any(!flag))
        {
            message("Some values of 'snp.chromosome' are invalid ",
                "(should be finite and >= 0)!")
            message("Hint: specifying 'autosome.only=FALSE' in ",
                "the analysis could avoid detecting chromosome coding.")
            warn.flag <- TRUE
            snp.flag <- snp.flag & flag
        }
    }

    # check snp allele
    if (!is.null(index.gdsn(gds, "snp.allele", silent=TRUE)))
    {
        dm <- objdesp.gdsn(index.gdsn(gds, "snp.allele"))$dim
        if ((length(dm) != 1) | (dm[1] != n.snp))
        {
            print(gds)
            stop("Invalid dimension of 'snp.allele'.")
        }
        snp.allele <- read.gdsn(index.gdsn(gds, "snp.allele"))
        snp.allele[is.na(snp.allele)] <- "?/?"
        flag <- sapply(strsplit(snp.allele, "/"),
            function(x)
            {
                if (length(x) == 2)
                {
                    all(x %in% c("A", "G", "C", "T"))
                } else {
                    FALSE
                }
            }
        )
        if (any(!flag))
        {
            s <- as.character((snp.allele[!flag])[1])
            message(sprintf(
                "Some of 'snp.allele' are not standard (e.g., %s).", s))
            warn.flag <- TRUE
            snp.flag <- snp.flag & flag
        }
    }

    # check genotype
    dm <- objdesp.gdsn(index.gdsn(gds, "genotype"))$dim
    if (length(dm) != 2)
    {
        print(gds)
        stop("Invalid dimension of 'genotype'.")
    }
    lv <- get.attr.gdsn(index.gdsn(gds, "genotype"))
    if ("sample.order" %in% names(lv))
    {
        if (dm[1]!=n.samp | dm[2]!=n.snp)
        {
            print(gds)
            stop("Invalid dimension of 'genotype'.")
        }
    } else {
        if (dm[2]!=n.samp | dm[1]!=n.snp)
        {
            print(gds)
            stop("Invalid dimension of 'genotype'.")
        }
    }

    # print
    if (show)
    {
        if (inherits(gds, "gds.class"))
            cat("The file name:", gds$filename, "\n")
        else
            cat("The file name:", gds, "\n")
        cat("The total number of samples:", n.samp, "\n")
        cat("The total number of SNPs:", n.snp, "\n")
        if ("sample.order" %in% names(lv))
        {
            cat("SNP genotypes are stored in SNP-major mode (Sample X SNP).\n")
        } else {
            cat("SNP genotypes are stored in individual-major mode",
                "(SNP X Sample).\n")
        }
        if (warn.flag)
        {
            cat("The number of valid samples:", length(samp.id), "\n")
            cat("The number of biallelic unique SNPs:", sum(snp.flag), "\n")
        }
    }

#   if (warn.flag)
#   {
#       warning("Call `snpgdsCreateGenoSet' to create a valid set ",
#           "of genotypes, using the returned sample.id and snp.id.")
#   }

    warn.flag <- FALSE
    snp.chr <- snp.chr[snp.flag]
    snp.pos <- snp.pos[snp.flag]
    if (is.numeric(snp.chr))
        chrtab <- setdiff(unique(snp.chr), 0)
    else
        chrtab <- unique(snp.chr)

    for (chr in chrtab)
    {
        pos <- snp.pos[snp.chr == chr]
        if (length(pos) > 0)
        {
            if (any(order(pos) != seq_len(length(pos))))
            {
                message(sprintf(
            "The SNP positions are not in ascending order on chromosome %s.",
                    chr))
                warn.flag <- TRUE
                break
            }
        }
    }

    # check -- sample annotation
    if (!is.null(index.gdsn(gds, "sample.annot", silent=TRUE)))
    {
        # sample.id
        if (!is.null(index.gdsn(gds, "sample.annot/sample.id", silent=TRUE)))
        {
            s <- read.gdsn(index.gdsn(gds, "sample.annot/sample.id"))
            if (length(s) != length(samp.id))
                warning("Invalid length of 'sample.annot/sample.id'.")
            if (any(s != samp.id))
                warning("Invalid 'sample.annot/sample.id'.")
        }

        # others
        lst <- ls.gdsn(index.gdsn(gds, "sample.annot"))
        for (n in lst)
        {
            dm <- objdesp.gdsn(index.gdsn(gds, index=c("sample.annot", n)))$dim
            if (!is.null(dm))
            {
                if (dm[1] != length(samp.id))
                    warning(sprintf("Invalid 'sample.annot/%s'.", n))
            }
        }
    }

    # check -- snp annotation
    if (!is.null(index.gdsn(gds, "snp.annot", silent=TRUE)))
    {
        if (!is.null(index.gdsn(gds, "snp.annot/snp.id", silent=TRUE)))
        {
            s <- read.gdsn(index.gdsn(gds, "snp.annot/snp.id"))
            if (length(s) != length(snp.id))
                warning("Invalid length of 'snp.annot/snp.id'.")
            if (any(s != snp.id))
                warning("Invalid 'snp.annot/snp.id'.")
        }

        # others
        lst <- ls.gdsn(index.gdsn(gds, "snp.annot"))
        for (n in lst)
        {
            dm <- objdesp.gdsn(index.gdsn(gds, index=c("snp.annot", n)))$dim
            if (!is.null(dm))
            {
                if (dm[1] != length(snp.id))
                    warning(sprintf("Invalid 'snp.annot/%s'.", n))
            }
        }
    }

    invisible(
        list(sample.id = samp.id, snp.id = snp.id[snp.flag])
    )
}



#######################################################################
# To get a subset of genotypes from a gds file
#

snpgdsGetGeno <- function(gdsobj, sample.id=NULL, snp.id=NULL,
    snpfirstdim=NA, .snpread=NA, with.id=FALSE, verbose=TRUE)
{
    if (is.character(gdsobj))
    {
        gdsobj <- snpgdsOpen(gdsobj, readonly=TRUE, allow.duplicate=TRUE)
        on.exit({ snpgdsClose(gdsobj) })
    }

    # check
    ws <- .InitFile(gdsobj, sample.id=sample.id, snp.id=snp.id,
        with.id=with.id)

    # get genotypes
    ans <- .Call(gnrCopyGenoMem, snpfirstdim, .snpread, verbose)

    if (with.id)
        ans <- list(genotype=ans, sample.id=ws$sample.id, snp.id=ws$snp.id)
    ans
}



#######################################################################
# To create a gds file for SNP genotypes
#

snpgdsCreateGeno <- function(gds.fn, genmat, sample.id=NULL, snp.id=NULL,
    snp.rs.id=NULL, snp.chromosome=NULL, snp.position=NULL, snp.allele=NULL,
    snpfirstdim=TRUE, compress.annotation="ZIP_RA.max", compress.geno="",
    other.vars=NULL)
{
    # check
    stopifnot(is.matrix(genmat))
    stopifnot(is.numeric(genmat))
    if (snpfirstdim)
    {
        n.snp <- dim(genmat)[1]; n.samp <- dim(genmat)[2]
    } else {
        n.snp <- dim(genmat)[2]; n.samp <- dim(genmat)[1]
    }

    if (!is.null(sample.id))
    {
        if (n.samp != length(sample.id))
        {
            stop("'snpfirstdim=", snpfirstdim,
                "', but length(sample.id) is not the number of samples.")
        }
        if (anyDuplicated(sample.id) > 0) stop("sample.id is not unique!")
    } else
        sample.id <- 1:n.samp
    if (!is.null(snp.id))
    {
        if (n.snp != length(snp.id))
        {
            stop("'snpfirstdim=", snpfirstdim,
                "', but length(snp.id) is not the number of SNPs.")
        }
        if (anyDuplicated(snp.id) > 0) stop("snp.id is not unique!")
    } else
        snp.id <- 1:n.snp

    if (!is.null(snp.rs.id))
        stopifnot(n.snp == length(snp.rs.id))
    if (!is.null(snp.chromosome))
        stopifnot(n.snp == length(snp.chromosome))
    if (!is.null(snp.position))
        stopifnot(n.snp == length(snp.position))
    stopifnot(is.null(other.vars) | is.list(other.vars))

    # create a gds file
    gfile <- createfn.gds(gds.fn)

    # close the gds file
    on.exit({ closefn.gds(gfile) })

    # add a flag
    put.attr.gdsn(gfile$root, "FileFormat", "SNP_ARRAY")

    add.gdsn(gfile, "sample.id", sample.id, compress=compress.annotation,
        closezip=TRUE)

    add.gdsn(gfile, "snp.id", snp.id, compress=compress.annotation,
        closezip=TRUE)

    if (!is.null(snp.rs.id))
    {
        add.gdsn(gfile, "snp.rs.id", snp.rs.id, compress=compress.annotation,
            closezip=TRUE)
    }

    if (is.null(snp.position))
        snp.position <- as.integer(1:n.snp)
    add.gdsn(gfile, "snp.position", as.integer(snp.position),
        compress=compress.annotation, closezip=TRUE)

    if (is.null(snp.chromosome))
        snp.chromosome <- as.integer(rep(1, n.snp))
    add.gdsn(gfile, "snp.chromosome", snp.chromosome,
        compress=compress.annotation, closezip=TRUE)

    if (!is.null(snp.allele))
    {
        add.gdsn(gfile, "snp.allele", snp.allele, compress=compress.annotation,
            closezip=TRUE)
    }

    # add genotype
    genmat[is.na(genmat)] <- 3L
    genmat[!(genmat %in% c(0L,1L,2L))] <- 3L
    node.geno <- add.gdsn(gfile, "genotype", genmat, storage="bit2",
        compress=compress.geno)
    if (snpfirstdim)
        put.attr.gdsn(node.geno, "snp.order")
    else
        put.attr.gdsn(node.geno, "sample.order")

    # other variables
    if (!is.null(other.vars))
    {
        for (i in 1:length(other.vars))
        {
            nm <- names(other.vars)[i]
            add.gdsn(gfile, nm, val=other.vars[[i]],
                compress=compress.annotation)
        }
    }

    # return
    invisible()
}



#######################################################################
# To create a gds file from a specified gds file
#

snpgdsCreateGenoSet <- function(src.fn, dest.fn, sample.id=NULL, snp.id=NULL,
    snpfirstdim=NULL, compress.annotation="ZIP_RA.max", compress.geno="",
    verbose=TRUE)
{
    # check
    stopifnot(is.character(src.fn))
    stopifnot(is.character(dest.fn))
    stopifnot(is.logical(snpfirstdim) | is.null(snpfirstdim))

    if (verbose)
        cat("Create a GDS genotype file:\n");

    # open the existing gds file
    srcobj <- snpgdsOpen(src.fn, allow.duplicate=TRUE)
    on.exit({ snpgdsClose(srcobj) })

    # create a new GDS file
    destobj <- createfn.gds(dest.fn)
    on.exit({ closefn.gds(destobj) }, add=TRUE)

    # initialize the source file
    ws <- .InitFile(srcobj, sample.id=sample.id, snp.id=snp.id)

    # Sample IDs
    ws.samp.id <- read.gdsn(index.gdsn(srcobj, "sample.id"))
    if (!is.null(ws$samp.flag))
        ws.samp.id <- ws.samp.id[ws$samp.flag]
    # SNP IDs
    ws.snp.id <- read.gdsn(index.gdsn(srcobj, "snp.id"))
    if (!is.null(ws$snp.flag))
        ws.snp.id <- ws.snp.id[ws$snp.flag]

    if (verbose)
    {
        cat(sprintf("The new dataset consists of %d samples and %d SNPs\n",
            ws$n.samp, ws$n.snp))
    }


    ####  write to the destination file  ####

    # add a flag
    put.attr.gdsn(destobj$root, "FileFormat", "SNP_ARRAY")

    # sample.id
    add.gdsn(destobj, "sample.id", ws.samp.id, compress=compress.annotation,
        closezip=TRUE)
    if (verbose) cat("    write sample.id\n");

    # snp.id
    add.gdsn(destobj, "snp.id", ws.snp.id, compress=compress.annotation,
        closezip=TRUE)
    if (verbose) cat("    write snp.id\n");

    # snp.rs.id
    if (!is.null(index.gdsn(srcobj, "snp.rs.id", silent=TRUE)))
    {
        rs.id <- read.gdsn(index.gdsn(srcobj, "snp.rs.id"))
        if (!is.null(ws$snp.flag))
            rs.id <- rs.id[ws$snp.flag]
        add.gdsn(destobj, "snp.rs.id", rs.id, compress=compress.annotation,
            closezip=TRUE)
        if (verbose)
            cat("    write snp.rs.id\n");
    }

    # snp.position
    pos <- read.gdsn(index.gdsn(srcobj, "snp.position"))
    if (!is.null(ws$snp.flag))
        pos <- pos[ws$snp.flag]
    add.gdsn(destobj, "snp.position", pos, compress=compress.annotation,
        closezip=TRUE)
    if (verbose)
        cat("    write snp.position\n");

    # snp.chromosome
    chr <- read.gdsn(index.gdsn(srcobj, "snp.chromosome"))
    if (!is.null(ws$snp.flag))
        chr <- chr[ws$snp.flag]
    add.gdsn(destobj, "snp.chromosome", chr, compress=compress.annotation,
        closezip=TRUE)
    if (verbose)
        cat("    write snp.chromosome\n");

    # snp.allele
    if (!is.null(index.gdsn(srcobj, "snp.allele", silent=TRUE)))
    {
        allele <- read.gdsn(index.gdsn(srcobj, "snp.allele"))
        if (!is.null(ws$snp.flag))
            allele <- allele[ws$snp.flag]
        add.gdsn(destobj, "snp.allele", allele, compress=compress.annotation,
            closezip=TRUE)
        if (verbose)
            cat("    write snp.allele\n");
    }

    # snp order
    if (is.null(snpfirstdim))
    {
        snpfirstdim <- TRUE
        rd <- names(get.attr.gdsn(index.gdsn(srcobj, "genotype")))
        if ("snp.order" %in% rd) snpfirstdim <- TRUE
        if ("sample.order" %in% rd) snpfirstdim <- FALSE
    }

    if (verbose)
    {
        if (snpfirstdim)
        {
            cat("SNP genotypes are stored in individual-major mode",
                "(SNP X Sample).\n")
        } else {
            cat("SNP genotypes are stored in SNP-major mode (Sample X SNP).\n")
        }
    }

    if (snpfirstdim)
    {
        gGeno <- add.gdsn(destobj, "genotype", storage="bit2",
            valdim=c(ws$n.snp, ws$n.samp), compress="")
        put.attr.gdsn(gGeno, "snp.order")
    } else {
        gGeno <- add.gdsn(destobj, "genotype", storage="bit2",
            valdim=c(ws$n.samp, ws$n.snp), compress="")
        put.attr.gdsn(gGeno, "sample.order")
    }

    # call C function
    .Call(gnrCopyGeno, gGeno, snpfirstdim)

    # return
    invisible()
}



#######################################################################
# To merge gds files of SNP genotypes into a single gds file
#

snpgdsCombineGeno <- function(gds.fn, out.fn, method=c("position", "exact"),
    compress.annotation="ZIP_RA.MAX", compress.geno="ZIP_RA",
    same.strand=FALSE, snpfirstdim=FALSE, verbose=TRUE)
{
    # check
    stopifnot(is.character(gds.fn), length(gds.fn)>=2L)
    stopifnot(is.character(out.fn), length(out.fn)==1L)
    method <- match.arg(method)
    stopifnot(is.character(compress.annotation))
    stopifnot(is.character(compress.geno))
    stopifnot(is.logical(same.strand), length(same.strand)==1L)
    stopifnot(is.logical(snpfirstdim), length(snpfirstdim)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L)

    # open a list of GDS files
    gs <- vector("list", length(gds.fn))
    slst <- vector("list", length(gds.fn))
    snp1d <- logical(length(gds.fn))
    on.exit({
        for (i in seq_along(gs))
            if (!is.null(gs[[i]])) snpgdsClose(gs[[i]])
    })

    if (verbose)
        cat("Merge SNP GDS files:\n")

    for (i in seq_along(gs))
    {
        if (verbose) cat("    open '", gds.fn[i], "' ...\n", sep="")
        f <- gs[[i]] <- snpgdsOpen(gds.fn[i])
        slst[[i]] <- snpgdsSNPList(f)
        snp1d[i] <- TRUE
        rd <- names(get.attr.gdsn(index.gdsn(f, "genotype")))
        if ("snp.order" %in% rd) snp1d[i] <- TRUE
        if ("sample.order" %in% rd) snp1d[i] <- FALSE
        if (verbose)
        {
            cat("       ", prod(objdesp.gdsn(index.gdsn(f, "sample.id"))$dim),
                "samples,", prod(objdesp.gdsn(index.gdsn(f, "snp.id"))$dim),
                "SNPs\n")
        }
    }

    # samples
    samp_id <- read.gdsn(index.gdsn(gs[[1L]], "sample.id"))
    for (i in 2L:length(gs))
    {
        s <- read.gdsn(index.gdsn(gs[[i]], "sample.id"))
        samp_id <- intersect(samp_id, s)
    }

    if (length(samp_id) > 0L)
    {
        if (verbose)
            cat("Concatenating SNPs ...\n")
        # get the total number of SNPs
        n_snp <- 0L
        for (i in seq_along(gs))
            n_snp <- n_snp + prod(objdesp.gdsn(index.gdsn(gs[[i]], "snp.id"))$dim)

        # create a gds file
        gfile <- createfn.gds(out.fn)
        on.exit({ if (!is.null(gfile)) closefn.gds(gfile) }, add=TRUE)
        if (verbose)
        {
            cat("    create '", out.fn, "': ", length(samp_id), " samples, ",
                n_snp, " SNPs\n", sep="")
        }

        # genotypes in the first GDS file
        gn1 <- index.gdsn(gs[[1L]], "genotype")
        filefmt <- ifelse(objdesp.gdsn(gn1)$type=="Real", "IMPUTED_DOSAGE",
            "SNP_ARRAY")
        if (verbose)
            cat("    FileFormat = ", filefmt, "\n", sep="")

        put.attr.gdsn(gfile$root, "FileFormat", filefmt)
        add.gdsn(gfile, "sample.id", samp_id, compress=compress.annotation,
            closezip=TRUE)
        add.gdsn(gfile, "snp.id", seq_len(n_snp), compress=compress.annotation,
            closezip=TRUE)

        n <- index.gdsn(gs[[1L]], "snp.rs.id", silent=TRUE)
        if (!is.null(n))
        {
            ng <- add.gdsn(gfile, "snp.rs.id", storage="string",
                compress=compress.annotation)
            for (i in seq_along(gs))
            {
                n <- index.gdsn(gs[[i]], "snp.rs.id", silent=TRUE)
                if (is.null(n))
                {
                    n <- prod(objdesp.gdsn(index.gdsn(gs[[i]], "snp.id"))$dim)
                    append.gdsn(ng, rep("", n))
                } else {
                    append.gdsn(ng, n)
                }
            }
            readmode.gdsn(ng)
        }

        ng <- add.gdsn(gfile, "snp.position", storage="int32",
            compress=compress.annotation)
        for (i in seq_along(gs))
            append.gdsn(ng, index.gdsn(gs[[i]], "snp.position"))
        readmode.gdsn(ng)

        ng <- add.gdsn(gfile, "snp.chromosome", storage=index.gdsn(gs[[1L]],
            "snp.chromosome"), compress=compress.annotation)
        for (i in seq_along(gs))
            append.gdsn(ng, index.gdsn(gs[[i]], "snp.chromosome"))
        readmode.gdsn(ng)

        ng <- add.gdsn(gfile, "snp.allele", storage="string",
            compress=compress.annotation)
        for (i in seq_along(gs))
            append.gdsn(ng, index.gdsn(gs[[i]], "snp.allele"))
        readmode.gdsn(ng)

        sync.gds(gfile)

        # add a genotype node
        ng <- add.gdsn(gfile, "genotype", valdim=c(length(samp_id), 0L),
            storage=gn1, compress=compress.geno)
        put.attr.gdsn(ng, "sample.order")

        if (verbose)
            cat("    writing genotypes ...\n")
        for (i in seq_along(gs))
        {
            f <- gs[[i]]
            ii_samp <- match(samp_id, read.gdsn(index.gdsn(f, "sample.id")))
            idx <- 1L
            N <- prod(objdesp.gdsn(index.gdsn(f, "snp.id"))$dim)
            n <- index.gdsn(gs[[i]], "genotype")
            while (idx <= N)
            {
                k <- idx + 1024L
                if (k > N) k <- N + 1L
                L <- k - idx
                ii <- seq.int(idx, length.out=L)
                idx <- idx + L
                if (snp1d[i])
                {
                    g <- readex.gdsn(n, list(ii, ii_samp), simplify="none")
                    g <- t(g)
                } else {
                    g <- readex.gdsn(n, list(ii_samp, ii), simplify="none")
                }
                append.gdsn(ng, g)
            }
        }

    } else {
        if (verbose)
            cat("Concatenating samples (mapping to the first GDS file) ...\n")
        opt <- slst
        opt$method <- method
        opt$same.strand <- same.strand
        opt$verbose <- verbose
        snp <- do.call(snpgdsSNPListIntersect, opt)
        if (length(snp$idx1) <= 0L)
            stop("There is no common SNP.")
        n_snp <- length(snp$idx1)

        # get all samples
        samp_id <- NULL
        for (i in seq_along(gs))
        {
            samp_id <- c(samp_id, read.gdsn(index.gdsn(gs[[i]], "sample.id")))
            if (verbose)
            {
                if (i == 1L)
                {
                    cat(sprintf("    reference: %d SNPs (%.1f%%)\n", n_snp,
                        100.0*n_snp/prod(objdesp.gdsn(index.gdsn(gs[[1L]], "snp.id"))$dim)))
                } else {
                    x <- snp[[paste0("flag",i)]]
                    cat("    file ", i, ": ", sum(bitwAnd(x, 1L), na.rm=TRUE),
                        " allele flips, ", sum(bitwAnd(x, 4L), na.rm=TRUE),
                        " ambiguous locus/loci\n", sep="")
                    y <- sapply(0:7, function(v) sum(x==v, na.rm=TRUE))
                    y <- c(y, sum(is.na(x)))
                    a <- c("no flip", "flip", "no flip on different strand",
                        "flip on different strand",
                        "no flip / ambiguous strand", "flip / ambiguous strand",
                        "no flip / different and ambiguous strand",
                        "flip / different and ambiguous strand", "mismatch")
                    for (i in seq_along(y))
                    {
                        if (y[i] > 0)
                            cat("        [", a[i], "]: ", y[i], "\n", sep="")
                    }
                }
            }
        }

        # create a gds file
        gfile <- createfn.gds(out.fn)
        on.exit({ if (!is.null(gfile)) closefn.gds(gfile) }, add=TRUE)
        if (verbose)
        {
            cat("    create '", out.fn, "': ", length(samp_id), " samples, ",
                n_snp, " SNPs\n", sep="")
        }

        # genotypes in the first GDS file
        gn1 <- index.gdsn(gs[[1L]], "genotype")
        filefmt <- ifelse(objdesp.gdsn(gn1)$type=="Real", "IMPUTED_DOSAGE",
            "SNP_ARRAY")
        if (verbose)
            cat("    FileFormat = ", filefmt, "\n", sep="")

        put.attr.gdsn(gfile$root, "FileFormat", filefmt)
        add.gdsn(gfile, "sample.id", samp_id, compress=compress.annotation,
            closezip=TRUE)
        add.gdsn(gfile, "snp.id", slst[[1L]]$snp.id[snp$idx1],
            compress=compress.annotation, closezip=TRUE)

        n <- index.gdsn(gs[[1L]], "snp.rs.id", silent=TRUE)
        if (!is.null(n))
        {
            add.gdsn(gfile, "snp.rs.id", read.gdsn(n)[snp$idx1],
                compress=compress.annotation, closezip=TRUE)
        }

        add.gdsn(gfile, "snp.position", slst[[1L]]$position[snp$idx1],
            compress=compress.annotation, closezip=TRUE)
        add.gdsn(gfile, "snp.chromosome", slst[[1L]]$chromosome[snp$idx1],
            compress=compress.annotation, closezip=TRUE)
        add.gdsn(gfile, "snp.allele", slst[[1L]]$allele[snp$idx1],
            compress=compress.annotation, closezip=TRUE)
        sync.gds(gfile)

        ng <- add.gdsn(gfile, "genotype", valdim=c(length(samp_id), 0L),
            storage=gn1, compress=compress.geno)
        put.attr.gdsn(ng, "sample.order")
        if (verbose)
            cat("    writing genotypes ...\n")

        idx <- 1L
        N <- length(snp$idx1)
        while (idx <= N)
        {
            k <- idx + 1024L
            if (k > N) k <- N + 1L
            L <- k - idx
            ii <- seq.int(idx, length.out=L)
            idx <- idx + L
            geno <- NULL
            for (i in 1:length(gds.fn))
            {
                n <- index.gdsn(gs[[i]], "genotype")
                if (snp1d[i])
                {
                    g <- readex.gdsn(n, list(snp[[paste0("idx",i)]][ii], NULL),
                        simplify="none", .value=3L, .substitute=NA_integer_)
                    g <- t(g)
                } else {
                    g <- readex.gdsn(n, list(NULL, snp[[paste0("idx",i)]][ii]),
                        simplify="none", .value=3L, .substitute=NA_integer_)
                }
                if (i > 1L)
                {
                    x <- ifelse(bitwAnd(snp[[paste0("flag",i)]][ii], 1), 2L, 0L)
                    g <- abs(t(x - t(g)))
                }
                geno <- rbind(geno, g)
            }
            geno[is.na(geno)] <- 3L
            append.gdsn(ng, geno)
        }
    }

    if (snpfirstdim)
    {
        if (verbose)
            cat("    transposing the genotype matrix ...\n")
        readmode.gdsn(index.gdsn(gfile, "genotype"))
        newnode <- add.gdsn(gfile, "!genotype", storage="bit2",
            valdim=c(n_snp, 0L), compress=compress.geno)
        put.attr.gdsn(newnode, "snp.order")
        # write data
        apply.gdsn(ng, margin=1L, FUN=c, as.is="gdsnode", target.node=newnode,
            .useraw=TRUE)
        readmode.gdsn(newnode)
        # move and replace
        moveto.gdsn(newnode, ng, relpos="replace+rename")
    }

    # close file
    closefn.gds(gfile)
    gfile <- NULL
    cleanup.gds(out.fn, verbose=verbose)
    if (verbose) cat("Done.\n")

    # return
    invisible()
}



#######################################################################
# To transpose the genotypic matrix if needed
#

snpgdsTranspose <- function(gds.fn, snpfirstdim=FALSE, compress=NULL,
    optimize=TRUE, verbose=TRUE)
{
    stopifnot(is.character(gds.fn))
    stopifnot(is.na(snpfirstdim) | is.logical(snpfirstdim))
    stopifnot(is.null(compress) | is.character(compress))
    stopifnot(is.logical(optimize))
    stopifnot(is.logical(verbose))

    # open the GDS file
    gds <- snpgdsOpen(gds.fn, readonly=FALSE)
    on.exit({ snpgdsClose(gds) })

    # check dimension
    ns <- names(get.attr.gdsn(index.gdsn(gds, "genotype")))
    if ("snp.order" %in% ns) snpflag <- TRUE
    if ("sample.order" %in% ns) snpflag <- FALSE
    node <- index.gdsn(gds, "genotype")

    if (is.na(snpfirstdim))
        snpfirstdim <- !snpflag

    # dimension
    desp <- objdesp.gdsn(node)
    dm <- desp$dim
    if (verbose)
    {
        if (snpflag)
            cat(sprintf("SNP genotypes: %d samples, %d SNPs\n", dm[2], dm[1]))
        else
            cat(sprintf("SNP genotypes: %d samples, %d SNPs\n", dm[1], dm[2]))
    }

    if (snpflag != snpfirstdim)
    {
        if (verbose)
            cat("Genotype matrix is being transposed ...\n")

        # check
        dm <- rev(dm)
        if (length(dm) != 2)
            stop("Invalid 'genotype' in the GDS file!")
        dm[length(dm)] <- 0

        # compress
        if (is.null(compress))
            compress <- desp$compress

        newnode <- add.gdsn(gds, "!genotype", val=NULL,
            storage=desp$storage, valdim=dm, compress=compress)

        # write data
        apply.gdsn(node, margin=1, FUN=c, as.is="gdsnode",
            target.node=newnode, .useraw=TRUE)
        readmode.gdsn(newnode)

        if (snpfirstdim)
            put.attr.gdsn(newnode, "snp.order")
        else
            put.attr.gdsn(newnode, "sample.order")

        # move
        moveto.gdsn(newnode, node, relpos="replace+rename")

        on.exit()
        snpgdsClose(gds)

        if (optimize)
            cleanup.gds(gds.fn, verbose=verbose)

    } else {
        if (verbose)
            cat("No transposing action taken.\n")

        if (!is.null(compress))
        {
            compression.gdsn(index.gdsn(gds, "genotype"),
                compress=compress)
        }
        on.exit()
        snpgdsClose(gds)

        if (optimize)
            cleanup.gds(gds.fn, verbose=verbose)
    }

    invisible()
}



#######################################################################
# To switch SNP coding based on allelic information
#

snpgdsAlleleSwitch <- function(gdsobj, A.allele, verbose=TRUE)
{
    # check
    .CheckFile(gdsobj)

    stopifnot(is.character(A.allele))
    stopifnot(is.vector(A.allele))

    allele.node <- index.gdsn(gdsobj, "snp.allele", silent=TRUE)
    if (is.null(allele.node))
    {
        stop("There is no allelic information in the GDS file ",
            "(no 'snp.allele' variable).")
    }
    dm <- objdesp.gdsn(allele.node)$dim
    if (length(dm) != 1)
        stop("Invalid 'snp.allele' in the GDS file.")
    if (dm != length(A.allele))
    {
        stop("The length of 'A.allele' should correspond to ",
            "'snp.allele' in the GDS file.")
    }

    stopifnot(is.logical(verbose))

    # check the GDS file
    if (gdsobj$readonly)
        stop("The GDS file should allow writing.")
    geno.node <- index.gdsn(gdsobj, "genotype")
    if (objdesp.gdsn(geno.node)$compress != "")
    {
        stop("Please call 'compression.gdsn(..., compress=\"\")' to ",
            "decompress the data first.")
    }

    # new alleles
    newnode <- add.gdsn(gdsobj, "!snp.allele", storage="string", valdim=c(0),
        compress="ZIP_RA.max")

    snpfirstdim <- TRUE
    rd <- names(get.attr.gdsn(geno.node))
    if ("snp.order" %in% rd) snpfirstdim <- TRUE
    if ("sample.order" %in% rd) snpfirstdim <- FALSE

    # call C function
    flag <- .Call(gnrStrandSwitch, geno.node, allele.node, newnode,
        snpfirstdim, A.allele)

    # new alleles
    readmode.gdsn(newnode)
    moveto.gdsn(newnode, allele.node, relpos="replace+rename")

    # synchronize the GDS file
    sync.gds(gdsobj)

    if (verbose)
    {
        nTrue <- sum(flag, na.rm=TRUE)
        nNA <- sum(is.na(flag))
        cat(sprintf("Strand-switching at %d SNP locus/loci.\n", nTrue))
        cat(sprintf("Unable to determine switching at %d SNP locus/loci.\n",
            nNA))
    }

    invisible(flag)
}




#######################################################################
# Plot functions
#######################################################################

#######################################################################
# Draw a dendrogram
#

snpgdsDrawTree <- function(obj, clust.count=NULL, dend.idx=NULL,
    type=c("dendrogram", "z-score"), yaxis.height=TRUE, yaxis.kinship=TRUE,
    y.kinship.baseline=NaN, y.label.kinship=FALSE, outlier.n=NULL,
    shadow.col = c(rgb(0.5, 0.5, 0.5, 0.25), rgb(0.5, 0.5, 0.5, 0.05)),
    outlier.col=rgb(1, 0.50, 0.50, 0.5), leaflab="none",
    labels=NULL, y.label=0.2, ...)
{
    # check
    stopifnot(is.null(dend.idx) | is.numeric(dend.idx))
    type <- match.arg(type)
    stopifnot(is.numeric(y.kinship.baseline))

    if (type == "dendrogram")
    {
        stopifnot(!is.null(obj$dendrogram))
        stopifnot(is.null(outlier.n) | is.numeric(outlier.n))

        # initialize ...
        if (is.null(clust.count))
            clust.count <- obj$clust.count
        if (is.null(outlier.n))
            outlier.n <- obj$outlier.n

        # draw
        if (!is.null(dend.idx))
        {
            den <- obj$dendrogram[[dend.idx]]
            x.offset <- 0
            for (i in 1:length(dend.idx))
            {
                if (dend.idx[i] == 2)
                {
                    IX <- dend.idx[1:i]
                    IX[i] <- 1
                    x.offset <- x.offset + attr(obj$dendrogram[[IX]], "member")
                }
            }
        } else {
            den <- obj$dendrogram
            x.offset <- 0
        }

        par(mar=c(4, 4, 4, 4))
        oldpar <- par(mgp=c(5, 1, 0))  # to avoid ylab
        plot(den, leaflab=leaflab, axes=FALSE, ...)
        par(oldpar)

        # draw left y-axis
        if (yaxis.height)
        {
            axis(side=2, line=0)
            tmp <- list(...)
            if (!is.null(tmp$ylab))
                ylab <- tmp$ylab
            else
                ylab <- "individual dissimilarity"
            mtext(ylab, side=2, line=2.5)
        }

        # draw right y-axis
        if (yaxis.kinship)
        {
            if (is.finite(y.kinship.baseline))
            {
                y.kinship.baseline <- y.kinship.baseline[1]
            } else {
                y.kinship.baseline <- attr(den, "height")
            }

            ym <- pretty(c(0, 1))
            axis(side=4, (1 - ym) * y.kinship.baseline, ym, line=0)
            mtext("coancestry coefficient", 4, line=2.5)
        }

        # draw others
        if (!is.null(clust.count))
        {
            m <- c(0, cumsum(clust.count))
            jj <- 1; k <- 1
            for (i in 1:length(clust.count))
            {
                if (clust.count[i] > outlier.n)
                {
                    rect(m[i] + 0.5 - x.offset, par("usr")[3L],
                        m[i+1] + 0.5 - x.offset, par("usr")[4L],
                        col = shadow.col[jj], border = NA)
                    jj <- 3 - jj
                    if (!is.null(labels[k]))
                        text((m[i]+m[i+1])/2 - x.offset, y.label, labels[k])
                    k <- k + 1
                } else {
                    rect(m[i] + 0.5 - x.offset, par("usr")[3L],
                        m[i+1] + 0.5 - x.offset, par("usr")[4L],
                        col = outlier.col, border = NA)
                }
            }
        }

        # draw kinship label
        if (yaxis.kinship & y.label.kinship)
        {
            # identical twins
            h1 <- (1 - 0.5)*y.kinship.baseline
            abline(h=h1, lty=2, col="gray")

            # parent-child / full-siblings
            h2 <- (1 - 0.25)*y.kinship.baseline
            abline(h=h2, lty=2, col="gray")

            # parent-child / full-siblings
            h3 <- (1 - 1/8)*y.kinship.baseline
            abline(h=h3, lty=2, col="gray")

            # first cousins
            h4 <- (1 - 1/16)*y.kinship.baseline
            abline(h=h4, lty=2, col="gray")

            axis(side=4, c(h1, h2, h3, h4),
                c("twins", "PC/FS", "DFC/HS", "FC"),
                tick=FALSE, line=-0.75, las=2, cex.axis=0.75,
                col.axis="gray25")
        }

    } else if (type == "z-score")
    {
        # the distribution of Z scores
        if (is.null(obj$merge))
            stop("There is no Z score in this object.")

        y <- obj$merge[,1]
        y <- y[order(y, decreasing=TRUE)]

        plot(y, xlab="the order of Z score", ylab="Z score", type="b",
            pch="+", log="x", ...)
        abline(h=15, col="gray", lty=2)
    }

    invisible()
}



#######################################################################
# SNPRelate Option
#

snpgdsOption <- function(gdsobj=NULL, autosome.start=1L, autosome.end=22L, ...)
{
    ans <- list(
        # the starting index of autosome
        autosome.start = as.integer(autosome.start),
        # the ending idex of autosome
        autosome.end   = as.integer(autosome.end)
    )

    if (!is.null(gdsobj))
    {
        # ignore the arguments "..."

        # chromosome
        n <- index.gdsn(gdsobj, "snp.chromosome")

        if (objdesp.gdsn(n)$type == "String")
        {
            rv <- .Call(gnrChromParse, n)
            names(rv) <- c("autosome.start", "autosome.end", "coding")
            return(rv)
        }

        lst <- get.attr.gdsn(n)
        if (!is.null(lst$autosome.start))
        {
            ans$autosome.start <- lst$autosome.start
            lst <- lst[-match("autosome.start", names(lst))]
        }
        if (!is.null(lst$autosome.end))
        {
            ans$autosome.end <- lst$autosome.end
            lst <- lst[-match("autosome.end", names(lst))]
        }

        ns <- names(lst)
        if (!("X" %in% ns))       # X chromosome
            lst$X <- as.integer(ans$autosome.end + 1)
        if (!("XY" %in% ns))      # Pseudo-autosomal region of X
            lst$XY <- as.integer(ans$autosome.end + 2)
        if (!("Y" %in% ns))       # Y chromosome
            lst$Y <- as.integer(ans$autosome.end + 3)
        if (!("M" %in% ns))       # Mitochondrial
            lst$M <- as.integer(ans$autosome.end + 4)
        if (!("MT" %in% ns))      # Mitochondrial
            lst$MT = as.integer(ans$autosome.end + 4)

        ans$chromosome.code <- lst

        # SNP genotype
        ans$file$filename <- gdsobj$filename

        n <- index.gdsn(gdsobj, "genotype")
        lst <- get.attr.gdsn(n)
        if ("sample.order" %in% names(lst))
            ans$file$geno.dim <- "sample-by-snp"
        else
            ans$file$geno.dim <- "snp-by-sample"

    } else {
        # incorporate the arguments "..."

        lst <- list(...)
        lst <- lst[names(lst) != ""]

        if (length(lst) <= 0)
        {
            ans$chromosome.code = list(
                X  = as.integer(autosome.end + 1),    # X chromosome
                XY = as.integer(autosome.end + 2),    # Pseudo-autosomal X
                Y  = as.integer(autosome.end + 3),    # Y chromosome
                M  = as.integer(autosome.end + 4),    # Mitochondrial
                MT = as.integer(autosome.end + 4)     # Mitochondrial
            )
        } else {
            ans$chromosome.code <- lst
        }
    }

    ans
}



#############################################################
# Sliding Windows Analysis
#

snpgdsSlidingWindow <- function(gdsobj, sample.id=NULL, snp.id=NULL,
    FUN=NULL, winsize=100000L, shift=10000L, unit=c("basepair", "locus"),
    winstart=NULL, autosome.only=FALSE, remove.monosnp=TRUE, maf=NaN,
    missing.rate=NaN, as.is=c("list", "numeric", "array"),
    with.id=c("snp.id", "snp.id.in.window", "none"), num.thread=1L,
    verbose=TRUE, ...)
{
    # check
    ws <- .InitFile2(
        cmd="Sliding Window Analysis:",
        gdsobj=gdsobj, sample.id=sample.id, snp.id=snp.id,
        autosome.only=autosome.only, remove.monosnp=remove.monosnp,
        maf=maf, missing.rate=missing.rate, num.thread=num.thread,
        verbose=verbose)

    stopifnot(is.numeric(winsize))
    stopifnot(is.numeric(shift))
    unit <- match.arg(unit)
    as.is <- match.arg(as.is)
    with.id <- match.arg(with.id)
    winsize <- as.integer(winsize)[1]
    shift <- as.integer(shift)[1]
    stopifnot(is.finite(winsize) & is.finite(shift))

    stopifnot(is.null(winstart) | (is.numeric(winstart) & is.vector(winstart)))

    if (is.function(FUN))
    {
        FUN <- match.fun(FUN)
        if (identical(FUN, snpgdsFst))
            stop('Please use `FUN="snpgdsFst"` instead.')
        if (identical(FUN, snpgdsSNPRateFreq))
            stop('Please use `FUN="snpgdsSNPRateFreq"` instead.')
        FunIdx <- 0L
    } else if (is.character(FUN))
    {
        stopifnot(is.vector(FUN))
        stopifnot(length(FUN) == 1)
        FunList <- c("snpgdsFst", "snpgdsSNPRateFreq")
        FunIdx <- match(FUN, FunList)
        if (is.na(FunIdx))
            stop("'FUN' should be one of ", paste(FunList, collapse=","), ".")
    } else {
        stop("'FUN' should be a function, or a character.")
    }

    if (verbose)
    {
        cat("    window size: ", winsize, ", shift: ", shift, sep="")
        cat(if (unit == "basepair") " (basepair)\n" else " (locus index)\n")
    }

    # check function
    if (FunIdx > 0)
    {
        pm <- list(...)
        param <- switch(EXPR = FUN,
            snpgdsFst = {
                .paramFst(sample.id, pm$population, pm$method, ws)
            },
            snpgdsSNPRateFreq = {
                if (length(pm) > 0)
                    stop("Unused additional parameters '...'.")
            }
        )
    }


    ########    ########

    # return value
    ans <- list(sample.id = ws$sample.id)
    if (with.id %in% c("snp.id", "snp.id.in.window"))
        ans$snp.id <- ws$snp.id

    if (inherits(gdsobj, "SeqVarGDSClass"))
    {
        total.snp.ids <- read.gdsn(index.gdsn(gdsobj, "variant.id"))
        snp.flag <- total.snp.ids %in% ws$snp.id

        chr <- read.gdsn(index.gdsn(gdsobj, "chromosome"))
        snp.flag[is.na(chr)] <- FALSE

        position <- read.gdsn(index.gdsn(gdsobj, "position"))
        snp.flag[!is.finite(position)] <- FALSE
        snp.flag[position <= 0] <- FALSE
    } else {
        total.snp.ids <- read.gdsn(index.gdsn(gdsobj, "snp.id"))
        snp.flag <- total.snp.ids %in% ws$snp.id

        chr <- read.gdsn(index.gdsn(gdsobj, "snp.chromosome"))
        snp.flag[is.na(chr)] <- FALSE

        position <- read.gdsn(index.gdsn(gdsobj, "snp.position"))
        snp.flag[!is.finite(position)] <- FALSE
        snp.flag[position <= 0] <- FALSE
    }

    if (is.numeric(chr))
        chrset <- setdiff(unique(chr[snp.flag]), c(0, NA))
    else if (is.character(chr))
        chrset <- setdiff(unique(chr[snp.flag]), c("", NA))
    else
        stop("Unknown format of 'snp.chromosome'!")

    if (!is.null(winstart))
    {
        winstart <- as.integer(winstart)
        stopifnot(all(is.finite(winstart)))
        if (length(winstart) != 1)
        {
            if (length(winstart) != length(chrset))
            {
                stop("'winstart' should be specified according to ",
                    "the chromosome set (", paste(chrset, collapse =","), ")")
            }
        }
    }

    if (verbose)
        cat("Chromosome Set: ", paste(chrset, collapse =","), "\n", sep="")

    # for-loop each chromosome
    for (ch in chrset)
    {
        # specific mask for this chromosome
        chflag <- snp.flag & (chr==ch)
        sid <- total.snp.ids[chflag]
        chpos <- position[chflag]

        winst <- NULL
        if (!is.null(winstart))
        {
            winst <- winstart[1]
            if (length(winstart) > 1)
                winstart <- winstart[-1]
        }

        if (verbose)
        {
            cat(date(), ", Chromosome ", ch,
                " (", length(chpos), " SNPs)", sep="")
        }

        if (is.function(FUN))
        {
            ####  the user-defined function

            # calculate how many blocks according to sliding windows
            if (unit == "basepair")
            {
                rg <- range(chpos)
                if (is.null(winst)) winst <- rg[1]
                n <- .Call(gnrSlidingNumWin, winst, rg[2], winsize, shift)
                if (verbose) cat(",", n, "windows\n")

                if (as.is == "list")
                    rvlist <- vector("list", n)
                else
                    rvlist <- double(n)
                nlist <- integer(n)
                poslist <- double(n)
                if (with.id == "snp.id.in.window")
                    sidlist <- vector("list", n)
            
                x <- rg[1]; i <- 1L
                while (i <= n)
                {
                    k <- (x <= chpos) & (chpos < x+winsize)
                    ssid <- sid[k]
                    ppos <- chpos[k]
                    v <- FUN(ans$sample.id, ssid, ppos, ...)
                    if (as.is == "list")
                        rvlist[[i]] <- v
                    else
                        rvlist[i] <- as.double(v)[1]
                    nlist[i] <- length(ppos)
                    poslist[i] <- mean(ppos)
                    if (with.id == "snp.id.in.window")
                        sidlist[[i]] <- ssid
                    x <- x + shift
                    i <- i + 1L
                }
            } else {
                if (is.null(winst)) winst <- 1L
                n <- .Call(gnrSlidingNumWin, winst, length(sid), winsize, shift)
                if (verbose) cat(",", n, "windows\n")

                if (as.is == "list")
                    rvlist <- vector("list", n)
                else
                    rvlist <- double(n)
                nlist <- integer(n)
                poslist <- double(n)
                if (with.id == "snp.id.in.window")
                    sidlist <- vector("list", n)

                x <- 1L; i <- 1L
                while (i <= n)
                {
                    k <- seq.int(x, x+winsize-1L)
                    ssid <- sid[k]
                    ppos <- chpos[k]
                    v <- FUN(ans$sample.id, ssid, ppos, ...)
                    if (as.is == "list")
                        rvlist[[i]] <- v
                    else
                        rvlist[i] <- as.double(v)[1]
                    nlist[i] <- length(ppos)
                    poslist[i] <- mean(ppos)
                    if (with.id == "snp.id.in.window")
                        sidlist[[i]] <- ssid
                    x <- x + shift
                    i <- i + 1L
                }
            }

            ans[[paste("chr", ch, sep="")]] <- rvlist
            ans[[paste("chr", ch, ".num", sep="")]] <- nlist
            ans[[paste("chr", ch, ".pos", sep="")]] <- poslist
            ans[[paste("chr", ch, ".posrange", sep="")]] <- range(chpos)
            if (with.id == "snp.id.in.window")
                ans[[paste("chr", ch, ".snpid", sep="")]] <- sidlist

        } else {
            ####  specific functions

            v <- .Call(gnrSlidingWindow, FunIdx, winsize, shift, unit, winst,
                as.is, chflag, chpos, param, verbose)

            ans[[paste("chr", ch, ".val", sep="")]] <- v[[1]]
            ans[[paste("chr", ch, ".num", sep="")]] <- v[[2]]
            ans[[paste("chr", ch, ".pos", sep="")]] <- v[[3]]
            ans[[paste("chr", ch, ".posrange", sep="")]] <- v[[4]]
        }
    }

    if (verbose) cat(date(), "\tDone.\n")

    # output
    ans
}





#######################################################################
# To get the file name of an example
#

snpgdsExampleFileName <- function()
{
    system.file("extdata", "hapmap_geno.gds", package="SNPRelate")
}



#######################################################################
# To get the error message
#

snpgdsErrMsg <- function()
{
    .Call(gnrErrMsg)
}
