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
# Identity-by-Descent (IBD) analysis
#######################################################################

#######################################################################
# Calculate the IBD matrix (PLINK method of moment)
#

snpgdsIBDMoM <- function(gdsobj, sample.id=NULL, snp.id=NULL,
    autosome.only=TRUE, remove.monosnp=TRUE, maf=NaN, missing.rate=NaN,
    allele.freq=NULL, kinship=FALSE, kinship.constraint=FALSE, num.thread=1,
    verbose=TRUE)
{
    # check
    ws <- .InitFile2(
        cmd="IBD analysis (PLINK method of moment) on genotypes:",
        gdsobj=gdsobj, sample.id=sample.id, snp.id=snp.id,
        autosome.only=autosome.only, remove.monosnp=remove.monosnp,
        maf=maf, missing.rate=missing.rate, allele.freq=allele.freq,
        num.thread=num.thread,
        verbose=verbose)

    stopifnot(is.logical(kinship))
    stopifnot(is.logical(kinship.constraint))

    # verbose
    if (verbose & !is.null(ws$allele.freq))
    {
        cat(sprintf("Specifying allele frequencies, mean: %0.3f, sd: %0.3f\n",
            mean(ws$allele.freq, na.rm=TRUE),
            sd(ws$allele.freq, na.rm=TRUE)))
        cat("*** A correction factor based on allele count is not used,",
            "since the allele frequencies are specified.\n")
    }

    # call C function
    rv <- .Call(gnrIBD_PLINK, ws$num.thread,
        as.double(ws$allele.freq), !is.null(ws$allele.freq),
        kinship.constraint, verbose)
    names(rv) <- c("k0", "k1", "afreq")

    # return
    rv <- list(sample.id = ws$sample.id, snp.id = ws$snp.id,
        afreq = rv$afreq, k0 = rv$k0, k1 = rv$k1)
    if (kinship)
        rv$kinship <- 0.5*(1 - rv$k0 - rv$k1) + 0.25*rv$k1
    rv$afreq[rv$afreq < 0] <- NaN
    class(rv) <- "snpgdsIBDClass"
    return(rv)
}



#######################################################################
# Calculate the identity-by-descent (IBD) matrix (MLE)
#

snpgdsIBDMLE <- function(gdsobj, sample.id=NULL, snp.id=NULL,
    autosome.only=TRUE, remove.monosnp=TRUE, maf=NaN, missing.rate=NaN,
    kinship=FALSE, kinship.constraint=FALSE, allele.freq=NULL,
    method=c("EM", "downhill.simplex", "Jacquard"), max.niter=1000L,
    reltol=sqrt(.Machine$double.eps), coeff.correct=TRUE, out.num.iter=TRUE,
    num.thread=1, verbose=TRUE)
{
    # check
    ws <- .InitFile2(
        cmd="Identity-By-Descent analysis (MLE) on genotypes:",
        gdsobj=gdsobj, sample.id=sample.id, snp.id=snp.id,
        autosome.only=autosome.only, remove.monosnp=remove.monosnp,
        maf=maf, missing.rate=missing.rate, allele.freq=allele.freq,
        num.thread=num.thread, verbose=verbose)

    method <- match.arg(method)
    if (method == "EM")
        method <- 0L
    else if (method == "downhill.simplex")
        method <- 1L
    else if (method == "Jacquard")
        method <- 2L
    else
        stop("Invalid MLE method!")

    stopifnot(is.logical(kinship))
    stopifnot(is.logical(kinship.constraint))
    stopifnot(is.numeric(max.niter))
    stopifnot(is.numeric(reltol))
    stopifnot(is.logical(coeff.correct))
    stopifnot(is.logical(out.num.iter))

    # check
    if (verbose & !is.null(ws$allele.freq))
    {
        cat(sprintf("Specifying allele frequencies, mean: %0.3f, sd: %0.3f\n",
            mean(ws$allele.freq, na.rm=TRUE),
            sd(ws$allele.freq, na.rm=TRUE)))
    }

    if (method != 2L)
    {
        # call C function
        rv <- .Call(gnrIBD_MLE, ws$allele.freq,
            as.logical(kinship.constraint), as.integer(max.niter),
            as.double(reltol), as.logical(coeff.correct), method,
            out.num.iter, ws$num.thread, verbose)

        # return
        rv <- list(sample.id=ws$sample.id, snp.id=ws$snp.id, afreq=rv[[3]],
            k0=rv[[1]], k1=rv[[2]], niter=rv[[4]])
        if (kinship)
            rv$kinship <- 0.5*(1 - rv$k0 - rv$k1) + 0.25*rv$k1
        rv$afreq[rv$afreq < 0] <- NaN
        class(rv) <- "snpgdsIBDClass"

    } else {
        # call C function
        rv <- .Call(gnrIBD_MLE_Jacquard, ws$allele.freq,
            as.integer(max.niter), as.double(reltol),
            as.logical(coeff.correct), method, out.num.iter, ws$num.thread,
            verbose)

        # return
        rv <- list(sample.id=ws$sample.id, snp.id=ws$snp.id, afreq=rv[[9]],
            D1=rv[[1]], D2=rv[[2]], D3=rv[[3]], D4=rv[[4]],
            D5=rv[[5]], D6=rv[[6]], D7=rv[[7]], D8=rv[[8]],
            niter=rv[[10]])
        if (kinship)
            rv$kinship <- rv$D1 + 0.5*(rv$D3 + rv$D5 + rv$D7) + 0.25*rv$D8
        rv$afreq[rv$afreq < 0] <- NaN
        class(rv) <- "snpgdsIBDClass"
    }

    rv
}



#######################################################################
# Calculate the identity-by-descent (IBD) matrix (MLE)
#

snpgdsIBDMLELogLik <- function(gdsobj, ibdobj, k0=NaN, k1=NaN,
    relatedness=c("", "self", "fullsib", "offspring", "halfsib", "cousin",
    "unrelated"))
{
    # check
    stopifnot(inherits(ibdobj, "snpgdsIBDClass"))
    .InitFile(gdsobj, ibdobj$sample.id, ibdobj$snp.id)

    stopifnot(is.numeric(k0) & is.vector(k0))
    stopifnot(length(k0) == 1)
    stopifnot(is.numeric(k1) & is.vector(k1))
    stopifnot(length(k1) == 1)

    relatedness <- match.arg(relatedness)
    if (relatedness == "self")
    {
        k0 <- 0; k1 <- 0
    } else if (relatedness == "fullsib")
    {
        k0 <- 0.25; k1 <- 0.5
    } else if (relatedness == "offspring")
    {
        k0 <- 0; k1 <- 1
    } else if (relatedness == "halfsib")
    {
        k0 <- 0.5; k1 <- 0.5
    } else if (relatedness == "cousin")
    {
        k0 <- 0.75; k1 <- 0.25
    } else if (relatedness == "unrelated")
    {
        k0 <- 1; k1 <- 0
    }

    # call C function
    if (is.finite(k0) & is.finite(k1))
    {
        .Call(gnrIBD_LogLik_k01, ibdobj$afreq, as.double(k0), as.double(k1))
    } else {
        .Call(gnrIBD_LogLik, ibdobj$afreq, ibdobj$k0, ibdobj$k1)
    }
}



#######################################################################
# To calculate the identity-by-descent (IBD) for a pair of SNP
#   genotypes using MLE
#

snpgdsPairIBD <- function(geno1, geno2, allele.freq,
    method=c("EM", "downhill.simplex", "MoM"), kinship.constraint=FALSE,
    max.niter=1000, reltol=sqrt(.Machine$double.eps), coeff.correct=TRUE,
    out.num.iter=TRUE, verbose=TRUE)
{
    # check
    stopifnot(is.vector(geno1) & is.numeric(geno1))
    stopifnot(is.vector(geno2) & is.numeric(geno2))
    stopifnot(is.vector(allele.freq) & is.numeric(allele.freq))
    stopifnot(length(geno1) == length(geno2))
    stopifnot(length(geno1) == length(allele.freq))
    stopifnot(is.logical(kinship.constraint))
    stopifnot(is.logical(coeff.correct))
    method <- match.arg(method)

    # method
    if (method == "EM")
        method <- 0
    else if (method == "downhill.simplex")
        method <- 1
    else if (method == "MoM")
        method <- -1
    else
        stop("Invalid MLE method!")

    allele.freq[!is.finite(allele.freq)] <- -1
    flag <- (0 <= allele.freq) & (allele.freq <= 1)
    if (sum(flag) < length(geno1))
    {
        if (verbose)
        {
            cat("IBD MLE for", sum(flag), 
                "SNPs in total, after removing loci with",
                "invalid allele frequencies.\n",
            )
        }
        geno1 <- geno1[flag]; geno2 <- geno2[flag]
        allele.freq <- allele.freq[flag]
    }

    # call C code
    rv <- .Call(gnrPairIBD, as.integer(geno1), as.integer(geno2),
        as.double(allele.freq), kinship.constraint, as.integer(max.niter),
        as.double(reltol), coeff.correct, as.integer(method))

    # return
    ans <- data.frame(k0=rv[[1]], k1=rv[[2]], loglik=rv[[3]])
    if (out.num.iter) ans$niter <- rv[[4]]
    ans
}



#######################################################################
# Calculate the identity-by-descent (IBD) matrix (MLE)
#

snpgdsPairIBDMLELogLik <- function(geno1, geno2, allele.freq, k0=NaN, k1=NaN,
    relatedness=c("", "self", "fullsib", "offspring", "halfsib", "cousin",
    "unrelated"), verbose=TRUE)
{
    # check
    stopifnot(is.vector(geno1) & is.numeric(geno1))
    stopifnot(is.vector(geno2) & is.numeric(geno2))
    stopifnot(is.vector(allele.freq) & is.numeric(allele.freq))
    stopifnot(length(geno1) == length(geno2))
    stopifnot(length(geno1) == length(allele.freq))
    stopifnot(is.numeric(k0))
    stopifnot(is.numeric(k1))
    stopifnot(is.character(relatedness))

    allele.freq[!is.finite(allele.freq)] <- -1
    flag <- (0 <= allele.freq) & (allele.freq <= 1)
    if (sum(flag) < length(geno1))
    {
        if (verbose)
        {
            cat("IBD MLE for", sum(flag),
                "SNPs in total, after removing loci with",
                "invalid allele frequencies.\n",
            )
        }
        geno1 <- geno1[flag]; geno2 <- geno2[flag]
        allele.freq <- allele.freq[flag]
    }

    # relatedness
    relatedness <- relatedness[1]
    if (relatedness == "self")
    {
        k0 <- 0; k1 <- 0
    } else if (relatedness == "fullsib")
    {
        k0 <- 0.25; k1 <- 0.5
    } else if (relatedness == "offspring")
    {
        k0 <- 0; k1 <- 1
    } else if (relatedness == "halfsib")
    {
        k0 <- 0.5; k1 <- 0.5
    } else if (relatedness == "cousin")
    {
        k0 <- 0.75; k1 <- 0.25
    } else if (relatedness == "unrelated")
    {
        k0 <- 1; k1 <- 0
    }

    # call C code
    .Call(gnrPairIBDLogLik, as.integer(geno1), as.integer(geno2),
        as.double(allele.freq), as.double(k0), as.double(k1))
}



#######################################################################
# Identity-by-Descent (IBD) analysis using KING robust estimat
#######################################################################

#######################################################################
# Calculate the identity-by-descent (IBD) matrix (KING)
#

snpgdsIBDKING <- function(gdsobj, sample.id=NULL, snp.id=NULL,
    autosome.only=TRUE, remove.monosnp=TRUE, maf=NaN, missing.rate=NaN,
    type=c("KING-robust", "KING-homo"), family.id=NULL,
    num.thread=1, verbose=TRUE)
{
    # check
    ws <- .InitFile2(
        cmd="IBD analysis (KING method of moment) on genotypes:",
        gdsobj=gdsobj, sample.id=sample.id, snp.id=snp.id,
        autosome.only=autosome.only, remove.monosnp=remove.monosnp,
        maf=maf, missing.rate=missing.rate, num.thread=num.thread,
        verbose=verbose)

    type <- match.arg(type)

    # family
    stopifnot(is.null(family.id) | is.vector(family.id))
    if (!is.null(family.id))
    {
        if (ws$n.samp != length(family.id))
            stop("'length(family.id)' should be the number of samples.")
    }
    if (!is.null(sample.id))
        family.id <- family.id[match(sample.id, ws$sample.id)]

    # family id
    if (is.vector(family.id))
    {
        if (is.character(family.id))
            family.id[family.id == ""] <- NA
        family.id <- as.factor(family.id)
        if (verbose & (type=="KING-robust"))
        {
            cat("# of families: ", nlevels(family.id),
                ", and within- and between-family relationship ",
                "are estimated differently.\n",
                sep="")
        }
    } else {
        if (verbose & (type=="KING-robust"))
        {
            cat("No family is specified, and all individuals",
                "are treated as singletons.\n")
        }
        family.id <- rep(NA, ws$n.samp)
    }

    if (type == "KING-homo")
    {
        if (verbose)
            cat("Relationship inference in a homogeneous population.\n")

        # call the C function
        rv <- .Call(gnrIBD_KING_Homo, ws$num.thread, verbose)

        rv <- list(sample.id=ws$sample.id, snp.id=ws$snp.id, afreq=NULL,
            k0=rv[[1]], k1=rv[[2]])

    } else if (type == "KING-robust")
    {
        if (verbose)
        {
            cat("Relationship inference in the presence of",
                "population stratification.\n")
        }

        # call the C function
        rv <- .Call(gnrIBD_KING_Robust, as.integer(family.id),
            ws$num.thread, verbose)

        rv <- list(sample.id=ws$sample.id, snp.id=ws$snp.id, afreq=NULL,
            IBS0=rv[[1]], kinship=rv[[2]])
    } else
        stop("Invalid 'type'.")

    # return
    rv$afreq[rv$afreq < 0] <- NaN
    class(rv) <- "snpgdsIBDClass"
    return(rv)
}




#######################################################################
# Genetic dissimilarity analysis
#######################################################################

#######################################################################
# Calculate the genetic dissimilarity matrix
#

snpgdsDiss <- function(gdsobj, sample.id=NULL, snp.id=NULL, autosome.only=TRUE,
    remove.monosnp=TRUE, maf=NaN, missing.rate=NaN, num.thread=1, verbose=TRUE)
{
    # check
    ws <- .InitFile2(
        cmd="Individual dissimilarity analysis on genotypes:",
        gdsobj=gdsobj, sample.id=sample.id, snp.id=snp.id,
        autosome.only=autosome.only, remove.monosnp=remove.monosnp,
        maf=maf, missing.rate=missing.rate, num.thread=num.thread,
        verbose=verbose)

    # call C function
    d <- .Call(gnrDiss, ws$num.thread, verbose)

    # return
    ans <- list(sample.id=ws$sample.id, snp.id=ws$snp.id, diss=d)
    class(ans) <- "snpgdsDissClass"
    return(ans)
}






#######################################################################
#
#######################################################################

#######################################################################
# Return a data.frame of pairs of individuals with IBD coefficients
#

snpgdsIBDSelection <- function(ibdobj, kinship.cutoff=NaN, samp.sel=NULL)
{
    # check
    stopifnot(inherits(ibdobj, "snpgdsIBDClass"))
    stopifnot(is.numeric(kinship.cutoff))
    stopifnot(is.null(samp.sel) | is.logical(samp.sel) | is.numeric(samp.sel))
    if (is.logical(samp.sel))
        stopifnot(length(samp.sel) == length(ibdobj$sample.id))

    # the variables in the output
    ns <- setdiff(names(ibdobj), c("sample.id", "snp.id", "afreq"))

    # subset
    if (!is.null(samp.sel))
    {
        ibdobj$sample.id <- ibdobj$sample.id[samp.sel]
        for (i in ns)
            ibdobj[[i]] <- ibdobj[[i]][samp.sel, samp.sel]
    }

    if (is.null(ibdobj$kinship))
    {
        if (!is.null(ibdobj$k0) && !is.null(ibdobj$k1))
        {
            ibdobj$kinship <- (1 - ibdobj$k0 - ibdobj$k1)*0.5 + ibdobj$k1*0.25
            ns <- c(ns, "kinship")
        } else if (!is.null(ibdobj$D1))
        {
            ibdobj$kinship <- ibdobj$D1 +
                0.5*(ibdobj$D3 + ibdobj$D5 + ibdobj$D7) + 0.25*ibdobj$D8
            ns <- c(ns, "kinship")
        } else {
            if (is.finite(kinship.cutoff))
                stop("There is no kinship coefficient.")
        }
    }

    if (is.finite(kinship.cutoff))
    {
        flag <- lower.tri(ibdobj$kinship) & (ibdobj$kinship >= kinship.cutoff)
        flag[is.na(flag)] <- FALSE
    } else
        flag <- lower.tri(ibdobj$kinship)

    # output
    n <- length(ibdobj$sample.id)
    ans <- data.frame(
        ID1 = .Call(gnrIBDSelSampID_List1, ibdobj$sample.id, flag),
        ID2 = .Call(gnrIBDSelSampID_List2, ibdobj$sample.id, flag),
        stringsAsFactors=FALSE)
    for (i in ns)
        ans[[i]] <- ibdobj[[i]][flag]

    ans
}



#######################################################################
# Genetic Relatedness
#######################################################################

#######################################################################
# Genetic relationship matrix (GRM)
#

snpgdsGRM <- function(gdsobj, sample.id=NULL, snp.id=NULL,
    autosome.only=TRUE, remove.monosnp=TRUE, maf=NaN, missing.rate=NaN,
    method=c("GCTA", "Eigenstrat", "EIGMIX", "IndivBeta"), num.thread=1L,
    with.id=TRUE, verbose=TRUE)
{
    # check and initialize ...
    method <- match.arg(method)
    ws <- .InitFile2(
        cmd=paste("Genetic Relationship Matrix (GRM, ", method, "):", sep=""),
        gdsobj=gdsobj, sample.id=sample.id, snp.id=snp.id,
        autosome.only=autosome.only, remove.monosnp=remove.monosnp,
        maf=maf, missing.rate=missing.rate, num.thread=num.thread,
        verbose=verbose)
    stopifnot(is.logical(with.id))

    # call GRM C function
    rv <- .Call(gnrGRM, ws$num.thread, method, verbose)

    # return
    if (with.id)
        rv <- list(sample.id=ws$sample.id, snp.id=ws$snp.id, grm=rv)

    return(rv)
}



#######################################################################
# F_st estimation
#

.paramFst <- function(sample.id, population, method=c("W&C84", "W&H02"), ws)
{
    method <- match.arg(method)
    stopifnot(is.factor(population))

    if (is.null(sample.id))
    {
        if (length(population) != ws$n.samp)
        {
            stop("The length of 'population' should be the number of samples ",
                "in the GDS file.")
        }
    } else {
        if (length(population) != length(sample.id))
        {
            stop("The length of 'population' should be the same as ",
                "the length of 'sample.id'.")
        }
        population <- population[match(ws$sample.id, sample.id)]
    }
    if (any(is.na(population)))
        stop("'population' should not have missing values!")
    if (nlevels(population) <= 1)
        stop("There should be at least two populations!")
    if (any(table(population) < 1))
        stop("Each population should have at least one individual.")

    if (ws$verbose)
    {
        if (method == "W&C84")
            cat("Method: Weir & Cockerham, 1984\n")
        else
            cat("Method: Weir & Hill, 2002\n")
        x <- table(population)
        cat("# of Populations: ", nlevels(population), "\n    ",
            paste(sprintf("%s (%d)", names(x), x), collapse=", "),
            "\n", sep="")
    }

    list(population=population, npop=nlevels(population), method=method)
}

snpgdsFst <- function(gdsobj, population, method=c("W&C84", "W&H02"),
    sample.id=NULL, snp.id=NULL, autosome.only=TRUE, remove.monosnp=TRUE,
    maf=NaN, missing.rate=NaN, with.id=FALSE, verbose=TRUE)
{
    # check
    ws <- .InitFile2(
        cmd="Fst estimation on genotypes:",
        gdsobj=gdsobj, sample.id=sample.id, snp.id=snp.id,
        autosome.only=autosome.only, remove.monosnp=remove.monosnp,
        maf=maf, missing.rate=missing.rate, num.thread=1L,
        verbose=verbose, verbose.numthread=FALSE)

    # check
    v <- .paramFst(sample.id, population, method, ws)

    # call C function
    d <- .Call(gnrFst, v$population, v$npop, v$method)

    # return
    if (with.id)
        rv <- list(sample.id=ws$sample.id, snp.id=ws$snp.id)
    else
        rv <- list()
    rv$Fst <- d[[1L]]
    if (method == "W&C84")
    {
        rv$MeanFst <- d[[2L]]
        rv$FstSNP <- d[[3L]]
    } else {
        rv$Beta <- d[[2L]]
        colnames(rv$Beta) <- rownames(rv$Beta) <- levels(population)
    }

    rv
}



#######################################################################
# Individual inbreeding and relatedness (beta)
#

snpgdsIndivBeta <- function(gdsobj, sample.id=NULL, snp.id=NULL,
    autosome.only=TRUE, remove.monosnp=TRUE, maf=NaN, missing.rate=NaN,
    method=c("weighted"), num.thread=1L, with.id=TRUE, verbose=TRUE)
{
    # check and initialize ...
    method <- match.arg(method)
    ws <- .InitFile2(cmd="Individual Inbreeding and Relatedness:",
        gdsobj=gdsobj, sample.id=sample.id, snp.id=snp.id,
        autosome.only=autosome.only, remove.monosnp=remove.monosnp,
        maf=maf, missing.rate=missing.rate, num.thread=num.thread,
        verbose=verbose)
    stopifnot(is.logical(with.id))

    # call GRM C function
    rv <- .Call(gnrIBD_Beta, ws$num.thread, verbose)

    # return
    if (with.id)
        rv <- list(sample.id=ws$sample.id, snp.id=ws$snp.id, beta=rv)
    return(rv)
}
