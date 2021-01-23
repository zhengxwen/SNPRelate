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
#


#######################################################################
# Conversion for Human Genome:
#     X    X chromosome                    -> 23
#     XY   Pseudo-autosomal region of X    -> 24
#     Y    Y chromosome                    -> 25
#     MT   Mitochondrial                   -> 26


#######################################################################
# Convert a GDS file to PLINK PED
#

snpgdsGDS2PED <- function(gdsobj, ped.fn, sample.id=NULL, snp.id=NULL,
    use.snp.rsid=TRUE, format=c("A/G/C/T", "A/B", "1/2"), verbose=TRUE)
{
    # check
    stopifnot(inherits(gdsobj, "gds.class"))
    stopifnot(is.character(ped.fn))
    format <- match.arg(format)

    # samples
    sample.ids <- read.gdsn(index.gdsn(gdsobj, "sample.id"))
    if (!is.null(sample.id))
    {
        n.tmp <- length(sample.id)
        sample.id <- sample.ids %in% sample.id
        n.samp <- sum(sample.id);
        if (n.samp != n.tmp)
            stop("Some of sample.id do not exist!")
        if (n.samp <= 0)
            stop("No sample in the working dataset.")
        sample.ids <- sample.ids[sample.id]
    }

    if (verbose)
        cat("Converting from GDS to PLINK PED:\n")

    # SNPs
    total.snp.ids <- read.gdsn(index.gdsn(gdsobj, "snp.id"))
    snp.ids <- total.snp.ids
    if (!is.null(snp.id))
    {
        n.tmp <- length(snp.id)
        snp.id <- snp.ids %in% snp.id
        n.snp <- sum(snp.id)
        if (n.snp != n.tmp)
            stop("Some of snp.id do not exist!")
        if (n.snp <= 0)
            stop("No SNP in the working dataset.")
        snp.ids <- snp.ids[snp.id]
    }

    # format code
    snp.idx <- match(snp.ids, total.snp.ids)
    if (format == "A/G/C/T")
    {
        n <- index.gdsn(gdsobj, "snp.allele", silent=TRUE)
        if (is.null(n))
            stop("There is no 'snp.allele' variable in the GDS file.")
        al <- read.gdsn(n)
        if (length(al) != length(total.snp.ids))
            stop("Invalid 'snp.allele' in the GDS file.")
        al <- al[snp.idx]
        fmt.code <- 1L
    } else if (format == "A/B")
    {
        al <- character(0)
        fmt.code <- 2L
    } else if (format == "1/2")
    {
        al <- character(0)
        fmt.code <- 3L
    } else
        stop("Invalid 'format'.")

    # output a MAP file
    tmp.snp.id <- snp.ids
    if (use.snp.rsid)
    {
        if (!is.null(index.gdsn(gdsobj, "snp.rs.id", silent=TRUE)))
        {
            tmp.snp.id <- read.gdsn(index.gdsn(gdsobj, "snp.rs.id"))[snp.idx]
        }
    }
    xchr <- as.character(read.gdsn(index.gdsn(gdsobj,
        "snp.chromosome")))[snp.idx]
    xchr[xchr=="23"] <- "X"; xchr[xchr=="25"] <- "Y"
    xchr[xchr=="24"] <- "XY"; xchr[xchr=="26"] <- "MT"
    D <- data.frame(chr = xchr, rs = tmp.snp.id,
        gen = rep(0, length(snp.idx)),
        base = read.gdsn(index.gdsn(gdsobj, "snp.position"))[snp.idx],
        stringsAsFactors = FALSE)
    write.table(D, file=paste(ped.fn, ".map", sep=""), sep="\t",
        quote=FALSE, row.names=FALSE, col.names=FALSE)
    if (verbose)
        cat("\tOutput a MAP file DONE.\n");

    # output a PED file
    if (verbose)
        cat("\tOutput a PED file ...\n");

    # set genotype working space
    .Call(gnrSetGenoSpace, index.gdsn(gdsobj, "genotype"), sample.id, snp.id)

    # run the C code
    .Call(gnrConvGDS2PED, paste(ped.fn, ".ped", sep=""),
        as.character(sample.ids), al, fmt.code, verbose)

    # return
    invisible()
}



#######################################################################
# Convert a PLINK PED file to GDS
#

snpgdsPED2GDS <- function(ped.fn, map.fn, out.gdsfn, family=TRUE,
    snpfirstdim=FALSE, compress.annotation="ZIP_RA.max", compress.geno="",
    verbose=TRUE)
{
    # check
    stopifnot(is.character(ped.fn))
    stopifnot(is.vector(ped.fn) & (length(ped.fn)==1L))

    stopifnot(is.character(map.fn))
    stopifnot(is.vector(map.fn) & (length(map.fn)==1L))

    stopifnot(is.character(out.gdsfn))
    stopifnot(is.vector(out.gdsfn) & (length(out.gdsfn)==1L))

    stopifnot(is.character(compress.annotation))
    stopifnot(is.vector(compress.annotation) & (length(compress.annotation)==1L))

    stopifnot(is.character(compress.geno))
    stopifnot(is.vector(compress.geno) & (length(compress.geno)==1L))

    stopifnot(is.logical(family))
    stopifnot(is.logical(snpfirstdim))
    stopifnot(is.logical(verbose))

    if (verbose) cat("PLINK PED/MAP to GDS Format:\n")

    ##  read MAP file
    f <- .OpenConnText(map.fn, TRUE)
    map <- read.table(f$con, header=FALSE, stringsAsFactors=FALSE)
    .CloseConnection(f)

    nsnp <- dim(map)[1]
    if (is.numeric(map$V1))
    {
        ii <- order(map$V1, map$V4)
    } else {
        chrcode <- suppressWarnings(as.integer(map$V1))
        ii <- order(chrcode, map$V1, map$V4)
    }    
    if (is.unsorted(ii))
        cat("Hint: the SNPs are sorted and merged into the GDS file.\n")
    map <- map[ii, ]
    map$V1[is.na(map$V1)] <- 0

    if (verbose)
    {
        cat(sprintf("Import %d variant%s from '%s'\n", nsnp,
            .plural(nsnp), map.fn))
        cat("Chromosome:")
        print(table(map$V1))
    }


    # create GDS file
    gfile <- createfn.gds(out.gdsfn)
    # close the file at the end
    on.exit({ closefn.gds(gfile) })

    # add file flag
    put.attr.gdsn(gfile$root, "FileFormat", "SNP_ARRAY")
    put.attr.gdsn(gfile$root, "FileVersion", "v1.0")

    # add "sample.id"
    add.gdsn(gfile, "sample.id", valdim=0, storage="string",
        compress=compress.annotation)
    # add "snp.id"
    add.gdsn(gfile, "snp.id", seq_len(nsnp), compress=compress.annotation,
        closezip=TRUE)
    # add "snp.rs.id"
    add.gdsn(gfile, "snp.rs.id", map$V2, compress=compress.annotation,
        closezip=TRUE)
    # add "snp.position"
    add.gdsn(gfile, "snp.position", map$V4, compress=compress.annotation,
        closezip=TRUE)

    # add "snp.chromosome"
    if (is.numeric(map$V1))
    {
        var_chr <- add.gdsn(gfile, "snp.chromosome", map$V1, storage="integer",
            compress=compress.annotation, closezip=TRUE)
        option <- snpgdsOption()
        put.attr.gdsn(var_chr, "autosome.start", option$autosome.start)
        put.attr.gdsn(var_chr, "autosome.end", option$autosome.end)
        for (i in 1:length(option$chromosome.code))
        {
            put.attr.gdsn(var_chr, names(option$chromosome.code)[i],
                option$chromosome.code[[i]])
        }
    } else {
        var_chr <- add.gdsn(gfile, "snp.chromosome", map$V1, storage="string",
            compress=compress.annotation, closezip=TRUE)
    }

    # add "snp.allele"
    add.gdsn(gfile, "snp.allele", valdim=0, storage="string",
        compress=compress.annotation)

    # add "genotype"
    comp.geno <- compress.geno
    if (!snpfirstdim) comp.geno <- ""
    gGeno <- add.gdsn(gfile, "genotype", storage="bit2", valdim=c(nsnp, 0),
        compress=comp.geno)
    put.attr.gdsn(gGeno, "snp.order")

    # add family information
    if (family)
    {
        v <- addfolder.gdsn(gfile, "sample.annot")
        put.attr.gdsn(v, "R.class", "data.frame")
        add.gdsn(v, "family", valdim=0, storage="string",
            compress=compress.annotation)
        add.gdsn(v, "father", valdim=0, storage="string",
            compress=compress.annotation)
        add.gdsn(v, "mother", valdim=0, storage="string",
            compress=compress.annotation)
        add.gdsn(v, "sex", valdim=0, storage="string",
            compress=compress.annotation)
        add.gdsn(v, "phenotype", valdim=0, storage="string",
            compress=compress.annotation)
    }

    # sync file
    sync.gds(gfile)


    # read PED file
    ped1 <- .OpenConnText(ped.fn)
    ped2 <- .OpenConnText(ped.fn)
    on.exit({ .CloseConnection(ped1); .CloseConnection(ped2) }, add=TRUE)
    if (verbose)
    {
        cat("Reading '", ped.fn, "'\n", sep="")
        cat("Output: '", out.gdsfn, "'\n", sep="")
    }

    # call C function
    .Call(gnrParsePED, ped.fn, gfile$root, ii - 1L,
        readLines, ped1$con, ped2$con, new.env(), verbose)
    nsamp <- objdesp.gdsn(gGeno)$dim[2]

    if (verbose)
        cat(sprintf("Import %d sample%s\n", nsamp, .plural(nsamp)))

    on.exit({ closefn.gds(gfile) })
    .CloseConnection(ped1)
    .CloseConnection(ped2)

    if (!snpfirstdim)
    {
        if (verbose) cat("Transpose the genotypic matrix ...\n")
        tm <- add.gdsn(gfile, "~genotype", storage="bit2", valdim=c(nsamp, 0),
            compress=compress.geno)
        put.attr.gdsn(tm, "sample.order")
        apply.gdsn(gGeno, margin=1, FUN=c, as.is="gdsnode", target.node=tm)
        readmode.gdsn(tm)
        moveto.gdsn(tm, gGeno, relpos="replace+rename")
    }

    on.exit()
    closefn.gds(gfile)

    if (verbose) cat("Done.\n")

    # optimize access efficiency
    if (verbose)
        cat("Optimize the access efficiency ...\n")
    cleanup.gds(out.gdsfn, verbose=verbose)

    # return
    invisible()
}



#######################################################################
# Convert a GDS file to PLINK Binary PED (BED) file
#

snpgdsGDS2BED <- function(gdsobj, bed.fn, sample.id=NULL, snp.id=NULL,
    snpfirstdim=NULL, verbose=TRUE)
{
    if (is.character(gdsobj))
    {
        gdsobj <- snpgdsOpen(gdsobj)
        on.exit({ snpgdsClose(gdsobj) })
    }

    # check
    ws <- .InitFile(gdsobj, sample.id=sample.id, snp.id=snp.id)

    stopifnot(is.character(bed.fn) & is.vector(bed.fn))
    stopifnot(length(bed.fn)==1L)
    stopifnot(is.logical(verbose) & is.vector(verbose))
    stopifnot(length(verbose)==1L)

    # snp order
    if (is.null(snpfirstdim))
    {
        snpfirstdim <- TRUE
        rd <- names(get.attr.gdsn(index.gdsn(gdsobj, "genotype")))
        if ("snp.order" %in% rd) snpfirstdim <- TRUE
        if ("sample.order" %in% rd) snpfirstdim <- FALSE
    } else {
        stopifnot(is.logical(snpfirstdim))
    }

    if (verbose)
    {
        cat("Converting from GDS to PLINK binary PED:\n")
        cat("Working space:", ws$n.samp, "samples,", ws$n.snp, "SNPs\n");
    }


    # Sample and SNP IDs
    total.samp.ids <- read.gdsn(index.gdsn(gdsobj, "sample.id"))
    total.snp.ids <- read.gdsn(index.gdsn(gdsobj, "snp.id"))
    if (is.null(ws$samp.flag))
        ws$samp.flag <- rep(TRUE, ws$n.samp)
    if (is.null(ws$snp.flag))
        ws$snp.flag <- rep(TRUE, ws$n.snp)

    # output a bim file
    xchr <- as.character(read.gdsn(index.gdsn(gdsobj,
        "snp.chromosome")))[ws$snp.flag]

    opt <- snpgdsOption(gdsobj)
    for (i in 1:length(opt$chromosome.code))
    {
        xchr[ xchr == opt$chromosome.code[[i]] ] <-
            names(opt$chromosome.code)[i]
    }
    xchr[is.na(xchr)] <- "0"

    if ((opt$autosome.start==1) & (opt$autosome.end==22))
    {
        # PLINK: Chromosome codes
        #   The autosomes should be coded 1 through 22.
        #   The following other codes can be used to specify
        #   other chromosome types:
        # X    X chromosome                    -> 23
        # Y    Y chromosome                    -> 24
        # XY   Pseudo-autosomal region of X    -> 25
        # MT   Mitochondrial                   -> 26

        xchr[xchr=="X"]  <- "23"; xchr[xchr=="Y"] <- "24"
        xchr[xchr=="XY"] <- "25"; xchr[xchr=="M"] <- "26"
        xchr[xchr=="MT"] <- "26"
    }

    if (!is.null(index.gdsn(gdsobj, "snp.allele", silent=TRUE)))
    {
        allele <- read.gdsn(index.gdsn(gdsobj, "snp.allele"))
        s <- unlist(strsplit(allele, "/"))
        ref <- s[seq(1, length(s), 2)]
        ref <- ref[ws$snp.flag]
        nonref <- s[seq(2, length(s), 2)]
        nonref <- nonref[ws$snp.flag]
    } else {
        warning("There is no allele information in the GDS file.",
            " ``A/B'' is used for the last two columns.")
        ref <- rep("A", ws$n.snp)
        nonref <- rep("B", ws$n.snp)
    }

    D <- data.frame(chr = xchr, rs = total.snp.ids[ws$snp.flag],
        gen.base = rep(0, ws$n.snp),
        base = read.gdsn(index.gdsn(gdsobj, "snp.position"))[ws$snp.flag],
        A1 = ref, A2 = nonref,
        stringsAsFactors = FALSE)
    write.table(D, file=paste(bed.fn, ".bim", sep=""), sep="\t",
        quote=FALSE, row.names=FALSE, col.names=FALSE)
    if (verbose)
        cat("Output a BIM file.\n");

    # output a fam file
    D <- data.frame(fam = rep(0, ws$n.samp),
        ind = total.samp.ids[ws$samp.flag],
        fat = rep(0, ws$n.samp), mot = rep(0, ws$n.samp),
        sex = rep(0, ws$n.samp), pheno = rep(-9, ws$n.samp),
        stringsAsFactors = FALSE)
    write.table(D, file=paste(bed.fn, ".fam", sep=""), sep="\t",
        quote=FALSE, row.names=FALSE, col.names=FALSE)

    # output a BED file
    if (verbose)
        cat("Output a BED file ...\n");

    # call C function
    .Call(gnrConvGDS2BED, path.expand(paste(bed.fn, ".bed", sep="")),
        snpfirstdim, verbose)

    if (verbose) cat("Done.\n")
    invisible()
}



#######################################################################
# Convert a PLINK BED file to GDS
#

snpgdsBED2GDS <- function(bed.fn, fam.fn, bim.fn, out.gdsfn, family=FALSE,
    snpfirstdim=NA, compress.annotation="LZMA_RA", compress.geno="",
    option=NULL, cvt.chr=c("int", "char"), cvt.snpid=c("auto", "int"),
    verbose=TRUE)
{
    # check
    stopifnot(is.character(bed.fn), length(bed.fn)==1L)
    if (missing(fam.fn) && missing(bim.fn))
    {
        fn <- gsub("\\.bed$", "", bed.fn, ignore.case=TRUE)
        bed.fn <- paste(fn, ".bed", sep="")
        fam.fn <- paste(fn, ".fam", sep="")
        bim.fn <- paste(fn, ".bim", sep="")
    }
    stopifnot(is.character(fam.fn), length(fam.fn)==1L)
    stopifnot(is.character(bim.fn), length(bim.fn)==1L)
    stopifnot(is.character(out.gdsfn), length(out.gdsfn)==1L)

    stopifnot(is.character(compress.annotation), length(compress.annotation)==1L)
    stopifnot(is.character(compress.geno), length(compress.geno)==1L)

    cvt.chr <- match.arg(cvt.chr)
    cvt.snpid <- match.arg(cvt.snpid)

    stopifnot(is.logical(family), length(family)==1L)
    stopifnot(is.na(snpfirstdim) | is.logical(snpfirstdim))
    stopifnot(is.logical(verbose), length(verbose)==1L)

    if (verbose)
        cat("Start file conversion from PLINK BED to SNP GDS ...\n")

    ##  open and detect bed.fn  ##

    bedfile <- .OpenConnBin(bed.fn)
    on.exit({ .CloseConnection(bedfile) })
    bed.flag <- .Call(gnrConvBEDFlag, bedfile$con, readBin, new.env())
    if (verbose)
    {
        cat("    BED file: ", shQuote(bed.fn), "\n", sep="")
        s <- .pretty_size(file.size(bed.fn))
        if (bed.flag == 0L)
            cat("        individual-major mode (SNP X Sample), ", s, "\n", sep="")
        else
            cat("        SNP-major mode (Sample X SNP), ", s, "\n", sep="")
    }

    ##  read fam.fn  ##

    f <- .OpenConnText(fam.fn, TRUE)
    famD <- read.table(f$con, header=FALSE, comment.char="",
        stringsAsFactors=FALSE)
    .CloseConnection(f)

    names(famD) <- c("FamilyID", "InvID", "PatID", "MatID", "Sex", "Pheno")
    if (anyDuplicated(famD$InvID) == 0L)
    {
        sample.id <- famD$InvID
    } else {
        sample.id <- paste(famD$FamilyID, famD$InvID, sep="-")
        if (length(unique(sample.id)) != dim(famD)[1])
            stop("IDs in PLINK BED are not unique!")
    }
    if (verbose)
        cat("    FAM file: ", shQuote(fam.fn), "\n", sep="")

    ##  read bim.fn  ##
    f <- .OpenConnText(bim.fn, TRUE)
    bimD <- read.table(f$con, header=FALSE, comment.char="",
        stringsAsFactors=FALSE)
    .CloseConnection(f)
    names(bimD) <- c("chr", "snp.id", "map", "pos", "allele1", "allele2")

    # chromosome
    if (cvt.chr == "int")
    {
        if (is.null(option)) option <- snpgdsOption()
        chrcode <- option$chromosome.code
        chr <- bimD$chr
        for (i in names(chrcode))
            chr[bimD$chr == i] <- chrcode[[i]]
        chr <- as.integer(chr)
        if (any(is.na(chr)))
        {
            warning(
            "Please use cvt.chr=\"char\" for non-numeric chromosome codes, otherwise non-numeric codes are replaced by zero.",
            immediate.=TRUE)
        }
        chr[is.na(chr)] <- 0L
    } else {
        if (!is.null(option))
            stop("'option' should be NULL when 'cvt.chr=\"char\"'.")
        chr <- as.character(bimD$chr)
    }

    # snp.id
    if (cvt.snpid == "auto")
    {
        if (anyDuplicated(bimD$snp.id) == 0L)
        {
            snp.id <- bimD$snp.id
            snp.rs.id <- NULL
        } else {
            snp.id <- seq_len(dim(bimD)[1L])
            snp.rs.id <- bimD$snp.id
        }
    } else {
        snp.id <- seq_len(dim(bimD)[1L])
        snp.rs.id <- bimD$snp.id
    }

    if (verbose)
        cat("    BIM file: ", shQuote(bim.fn), "\n", sep="")

    # create GDS file
    gfile <- createfn.gds(out.gdsfn)
    # close the file at the end
    on.exit({ closefn.gds(gfile) }, add=TRUE)

    # add file flag
    put.attr.gdsn(gfile$root, "FileFormat", "SNP_ARRAY")
    put.attr.gdsn(gfile$root, "FileVersion", "v1.0")

    # add "sample.id"
    add.gdsn(gfile, "sample.id", sample.id, compress=compress.annotation,
        closezip=TRUE)
    # add "snp.id"
    add.gdsn(gfile, "snp.id", snp.id, compress=compress.annotation,
        closezip=TRUE)
    # add "snp.rs.id"
    if (!is.null(snp.rs.id))
    {
        add.gdsn(gfile, "snp.rs.id", snp.rs.id, compress=compress.annotation,
            closezip=TRUE)
    }
    # add "snp.position"
    add.gdsn(gfile, "snp.position", bimD$pos, compress=compress.annotation,
        closezip=TRUE)

    # add "snp.chromosome"
    if (cvt.chr == "int")
    {
        v.chr <- add.gdsn(gfile, "snp.chromosome", chr, storage="uint8",
            compress=compress.annotation, closezip=TRUE)

        put.attr.gdsn(v.chr, "autosome.start", option$autosome.start)
        put.attr.gdsn(v.chr, "autosome.end", option$autosome.end)
        for (i in 1:length(option$chromosome.code))
        {
            put.attr.gdsn(v.chr, names(option$chromosome.code)[i],
                option$chromosome.code[[i]])
        }
    } else {
        add.gdsn(gfile, "snp.chromosome", chr, compress=compress.annotation,
            closezip=TRUE)
    }

    # add "snp.allele"
    add.gdsn(gfile, "snp.allele", paste(bimD$allele1, bimD$allele2, sep="/"),
        compress=compress.annotation, closezip=TRUE)

    # sync file
    sync.gds(gfile)

    nSamp <- dim(famD)[1L]; nSNP <- dim(bimD)[1L]
    if (verbose)
    {
        cat(date(), "    (store sample id, snp id, position, and chromosome)\n")
        cat(sprintf("    start writing: %d samples, %d SNPs ...\n", nSamp, nSNP))
    }

    # add "gonetype", 2 bits to store one genotype
    comp.geno <- compress.geno
    transposeflag <- FALSE
    if (bed.flag == 0L)
    {
        if (identical(snpfirstdim, FALSE))
        {
            comp.geno <- ""
            transposeflag <- TRUE
        }
        gGeno <- add.gdsn(gfile, "genotype", storage="bit2",
            valdim=c(nSNP, 0L), compress=comp.geno)
        put.attr.gdsn(gGeno, "snp.order")
        n <- nSamp
    } else {
        if (identical(snpfirstdim, TRUE))
        {
            comp.geno <- ""
            transposeflag <- TRUE
        }
        gGeno <- add.gdsn(gfile, "genotype", storage="bit2",
            valdim=c(nSamp, 0L), compress=comp.geno)
        put.attr.gdsn(gGeno, "sample.order")
        n <- nSNP
    }

    # convert
    .Call(gnrConvBED2GDS, gGeno, n, bedfile$con, readBin, new.env(), verbose)

    # close the BED file
    on.exit({ closefn.gds(gfile) })
    .CloseConnection(bedfile)

    # sync file
    sync.gds(gfile)

    # add "sample.annot"
    sex <- rep("", length(sample.id))
    sex[famD$Sex==1L] <- "M"; sex[famD$Sex==2L] <- "F"

    if (family)
    {
        samp.annot <- data.frame(family=famD$FamilyID,
            father=famD$PatID, mother=famD$MatID,
            sex=sex, phenotype=famD$Pheno, stringsAsFactors=FALSE)
    } else {
        samp.annot <- data.frame(sex=sex, phenotype=famD$Pheno,
            stringsAsFactors=FALSE)
    }
    add.gdsn(gfile, "sample.annot", samp.annot, compress=compress.annotation,
        closezip=TRUE)

    if (transposeflag)
    {
        if (verbose)
            cat("Transpose the genotypic matrix ...\n")
        if (bed.flag == 0L)
        {
            tm <- add.gdsn(gfile, "~genotype", storage="bit2",
                valdim=c(nSamp, 0L), compress=compress.geno)
            put.attr.gdsn(tm, "sample.order")
        } else {
            tm <- add.gdsn(gfile, "~genotype", storage="bit2",
                valdim=c(nSNP, 0L), compress=compress.geno)
            put.attr.gdsn(tm, "snp.order")
        }
        apply.gdsn(gGeno, margin=1L, FUN=c, as.is="gdsnode", target.node=tm)
        readmode.gdsn(tm)
        moveto.gdsn(tm, gGeno, relpos="replace+rename")
    }

    on.exit()
    closefn.gds(gfile)

    if (verbose)
        cat(date(), "\tDone.\n")

    # optimize access efficiency
    if (verbose)
        cat("Optimize the access efficiency ...\n")
    cleanup.gds(out.gdsfn, verbose=verbose)

    # output
    invisible(normalizePath(out.gdsfn))
}



#######################################################################
# Convert a GDS file to Eigenstrat format
#

snpgdsGDS2Eigen <- function(gdsobj, eigen.fn, sample.id=NULL, snp.id=NULL,
    verbose=TRUE)
{
    # check
    stopifnot(inherits(gdsobj, "gds.class"))
    stopifnot(is.character(eigen.fn))

    # samples
    sample.ids <- read.gdsn(index.gdsn(gdsobj, "sample.id"))
    if (!is.null(sample.id))
    {
        n.tmp <- length(sample.id)
        sample.id <- sample.ids %in% sample.id
        n.samp <- sum(sample.id);
        if (n.samp != n.tmp)
            stop("Some of sample.id do not exist!")
        if (n.samp <= 0)
            stop("No sample in the working dataset.")
        sample.ids <- sample.ids[sample.id]
    } else
        sample.id <- rep(TRUE, length(sample.ids))

    if (verbose)
        cat("Converting from GDS to EIGENSOFT:\n")

    # SNPs
    total.snp.ids <- read.gdsn(index.gdsn(gdsobj, "snp.id"))
    snp.ids <- total.snp.ids
    if (!is.null(snp.id))
    {
        n.tmp <- length(snp.id)
        snp.id <- snp.ids %in% snp.id
        n.snp <- sum(snp.id)
        if (n.snp != n.tmp)
            stop("Some of snp.id do not exist!")
        if (n.snp <= 0)
            stop("No SNP in the working dataset.")
        snp.ids <- snp.ids[snp.id]
    } else
        snp.id <- rep(TRUE, length(snp.ids))

    # making the "*.snp" file ...
    tmpD <- data.frame(
        snpid = read.gdsn(index.gdsn(gdsobj, "snp.id"))[snp.id],
        chrom = read.gdsn(index.gdsn(gdsobj, "snp.chromosome"))[snp.id],
        map = rep(0.0, sum(snp.id)),
        pos = read.gdsn(index.gdsn(gdsobj, "snp.position"))[snp.id],
        stringsAsFactors = FALSE
    )
    write.table(tmpD, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE,
        file = paste(eigen.fn, ".snp", sep=""))
    if (verbose)
        cat("\tsave to *.snp:", dim(tmpD)[1], "snps\n")

    # making the "*.ind" file ...
    sex <- try(read.gdsn(index.gdsn(gdsobj, "sample.annot/sex")), TRUE)
    if (class(sex) == "try-error")
    {
        sex <- rep("U", sum(sample.id))
    } else {
        sex <- as.character(sex)[sample.id]
        if (!all(sex %in% c("F", "M"), na.rm=TRUE))
        {
            stop("The gender variable in GDS file should be ",
                "either \"M\" or \"F\".")
        }
        sex[is.na(sex)] <- "U"
    }
    tmpD <- data.frame(
        sampid = read.gdsn(index.gdsn(gdsobj, "sample.id"))[sample.id],
        gender = sex, label = rep("control", sum(sample.id)),
        stringsAsFactors = FALSE
    )
    write.table(tmpD, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE,
        file = paste(eigen.fn, ".ind", sep=""))
    if (verbose)
        cat("\tsave to *.ind:", dim(tmpD)[1], "samples\n")

    # making the "*.eigenstratgeno" file ...

    # set genotype working space
    .Call(gnrSetGenoSpace, index.gdsn(gdsobj, "genotype"), sample.id, snp.id)

    # call C function
    .Call(gnrConvGDS2EIGEN, paste(eigen.fn, ".eigenstratgeno", sep=""),
        verbose)

    if (verbose) cat("Done.\n")

    invisible()
}



#######################################################################
# Convert an Oxford GEN file to a GDS file
# http://www.stats.ox.ac.uk/%7Emarchini/software/gwas/file_format.html
# http://www.well.ox.ac.uk/~gav/bgen_format/bgen_format.html
#

snpgdsGEN2GDS <- function(gen.fn, sample.fn, out.fn, chr.code=NULL,
    call.threshold=0.9, version=c(">=2.0", "<=1.1.5"),
    snpfirstdim=FALSE, compress.annotation="ZIP_RA.max",
    compress.geno="", verbose=TRUE)
{
    #######################################################################
    # check

    stopifnot(is.character(gen.fn) & is.vector(gen.fn))
    stopifnot(!any(is.na(gen.fn)))

    stopifnot(is.character(sample.fn) & is.vector(sample.fn))
    stopifnot(length(sample.fn)==1L)

    stopifnot(is.character(out.fn) & is.vector(out.fn))
    stopifnot(length(out.fn)==1L)

    if (!is.null(chr.code))
    {
        stopifnot(is.vector(chr.code))
        stopifnot(is.numeric(chr.code) | is.character(chr.code))
        stopifnot(length(gen.fn) == length(chr.code))
    } else {
        stop("Please specify the argument 'chr.code', e.g., 'chr.code=1'.")
    }

    stopifnot(is.numeric(call.threshold) & is.vector(call.threshold))
    stopifnot(length(call.threshold)==1L)
    stopifnot(is.finite(call.threshold))

    version <- match.arg(version)

    stopifnot(is.character(compress.annotation))
    stopifnot(is.logical(snpfirstdim))
    stopifnot(is.character(compress.geno))
    stopifnot(is.logical(verbose))

    if (verbose)
        cat("Oxford GEN/BGEN format ---> GDS SNP format:\n")


    #######################################################################
    # running

    if (is.numeric(chr.code))
    {
        chr.code <- as.integer(chr.code)
        chr.code[is.na(chr.code)] <- 0L
    } else if (is.character(chr.code))
    {
        chr.code[is.na(chr.code)] <- ""
    }

    # read sample id
    if (version == ">=2.0")
    {
        tmp <- read.table(sample.fn, header=TRUE, nrows=1,
            stringsAsFactors=FALSE)
        samp.tab <- read.table(sample.fn, skip=2, stringsAsFactors=FALSE)
        names(samp.tab) <- names(tmp)
    } else if (version == "<=1.1.5")
    {
        samp.tab <- read.table(sample.fn, header=TRUE, stringsAsFactors=FALSE)
    }

    if (verbose)
    {
        cat(sprintf("The number of samples: %d, with SNPTEST version (%s).\n",
            dim(samp.tab)[1], version))
    }


    #######################################################################

    # create a GDS file
    gfile <- createfn.gds(out.fn)

    # close the file at the end
    on.exit(closefn.gds(gfile))

    # add file flag
    put.attr.gdsn(gfile$root, "FileFormat", "SNP_ARRAY")
    put.attr.gdsn(gfile$root, "FileVersion", "v1.0")

    # add sample id
    add.gdsn(gfile, "sample.id", samp.tab[,1], compress=compress.annotation,
        closezip=TRUE)

    # add snp.id
    add.gdsn(gfile, "snp.id", storage="string", valdim=c(0),
        compress=compress.annotation)

    # add snp.rs.id
    add.gdsn(gfile, "snp.rs.id", storage="string", valdim=c(0),
        compress=compress.annotation)

    # add position
    add.gdsn(gfile, "snp.position", storage="int32", valdim=c(0),
        compress=compress.annotation)

    # add chromosome
    add.gdsn(gfile, "snp.chromosome",
        storage = { if (is.numeric(chr.code)) "int" else "string" },
        valdim=c(0), compress=compress.annotation)

    # add allele
    add.gdsn(gfile, "snp.allele", storage="string", valdim=c(0),
        compress=compress.annotation)

    # add SNP genotypes
    cmp <- if (snpfirstdim) "" else compress.geno
    nodegeno <- add.gdsn(gfile, "genotype", storage="bit2",
        valdim=c(dim(samp.tab)[1], 0), compress=cmp)
    put.attr.gdsn(nodegeno, "sample.order")

    # add a folder for sample annotation
    add.gdsn(gfile, "sample.annot", val=samp.tab, compress=compress.annotation)


    ##################################################
    # for-loop each file
    for (i in 1:length(gen.fn))
    {
        opfile <- file(gen.fn[i], open="rt")
        on.exit({ closefn.gds(gfile); close(opfile) })

        if (verbose)
            cat("Parsing \"", gen.fn[i], "\" ...\n", sep="")

        # call C function
        n <- .Call(gnrParseGEN, gen.fn[i], gfile$root, chr.code[i],
            as.double(call.threshold),
            readLines, opfile, 1024L,  # "readLines(opfile, 1024L)"
            new.env(), verbose)

        if (verbose)
        {
            if (n > 1)
            {
                cat("\tImport ", n, " variants on chromosome ",
                    chr.code[i], ".\n", sep="")
            } else {
                cat("\tImport ", n, " variant on chromosome ",
                    chr.code[i], ".\n", sep="")
            }
        }

        on.exit()
        close(opfile)
    }

    closefn.gds(gfile)

    if (snpfirstdim)
    {
        snpgdsTranspose(out.fn, snpfirstdim=TRUE, compress=compress.geno,
            optimize=FALSE, verbose=verbose)
    }


    ##################################################
    # optimize access efficiency

    if (verbose)
        cat("Optimize the access efficiency ...\n")
    cleanup.gds(out.fn, verbose=verbose)

    # output
    invisible(normalizePath(out.fn))
}



#######################################################################
# Convert a VCF (sequence) file to GDS (extracting SNP data)
#

snpgdsVCF2GDS <- function(vcf.fn, out.fn,
    method = c("biallelic.only", "copy.num.of.ref"), snpfirstdim=FALSE,
    compress.annotation="LZMA_RA", compress.geno="",
    ref.allele=NULL, ignore.chr.prefix="chr", verbose=TRUE)
{
    # check
    stopifnot(is.character(vcf.fn) & is.vector(vcf.fn))
    stopifnot(length(vcf.fn) > 0)

    stopifnot(is.character(out.fn) & is.vector(out.fn))
    stopifnot(length(out.fn) == 1L)

    method <- match.arg(method)
    metidx <- match(method, c("biallelic.only", "copy.num.of.ref"))

    stopifnot(is.character(compress.annotation))
    stopifnot(is.logical(snpfirstdim))
    stopifnot(is.character(compress.geno))
    stopifnot(is.null(ref.allele) || is.character(ref.allele))
    if (is.character(ref.allele))
        stopifnot(is.vector(ref.allele))
    stopifnot(is.character(ignore.chr.prefix))
    stopifnot(is.logical(verbose))

    if (verbose)
    {
        cat("Start file conversion from VCF to SNP GDS ...\n")
        if (metidx == 1L)
            cat("Method: extracting biallelic SNPs\n")
        else
            cat("Method: dosage (0,1,2) of reference allele for all variant sites\n")
    }


    #######################################################################

    # get sample id from a VCF file
    VCF_SampID <- function(vcf.fn)
    {
        # open the vcf file
        opfile <- .OpenConnText(vcf.fn)
        on.exit({ .CloseConnection(opfile) })

        # read header
        samp.id <- NULL
        while (length(s <- readLines(opfile$con, n=1)) > 0)
        {
            if (substr(s, 1, 6) == "#CHROM")
            {
                samp.id <- scan(text=s, what=character(0), sep="\t",
                    quiet=TRUE)[-c(1:9)]
                break
            }
        }
        if (is.null(samp.id))
            stop("Error VCF format: invalid sample id!")

        samp.id
    }

    # read sample id
    samp.id <- NULL
    for (i in 1:length(vcf.fn))
    {
        if (is.null(samp.id))
        {
            samp.id <- VCF_SampID(vcf.fn[i])
            if (length(samp.id) <= 0)
                stop("No sample in the VCF file!")
        } else {
            tmp <- VCF_SampID(vcf.fn[i])
            if (length(samp.id) != length(tmp))
                stop(sprintf("'%s' has different sample id.", vcf.fn[i]))
            if (!identical(samp.id, tmp))
                stop(sprintf("'%s' has different sample id.", vcf.fn[i]))
        }
    }


    #######################################################################

    if (verbose)
        cat(sprintf("Number of samples: %d\n", length(samp.id)))

    # create a GDS file
    gfile <- createfn.gds(out.fn)

    # close the file at the end
    on.exit(closefn.gds(gfile))

    # add file flag
    put.attr.gdsn(gfile$root, "FileFormat", "SNP_ARRAY")
    put.attr.gdsn(gfile$root, "FileVersion", "v1.0")

    # add sample id
    add.gdsn(gfile, "sample.id", samp.id, compress=compress.annotation,
        closezip=TRUE)

    # add snp.id
    add.gdsn(gfile, "snp.id", storage="int32", valdim=c(0),
        compress=compress.annotation)

    # add snp.rs.id
    add.gdsn(gfile, "snp.rs.id", storage="string", valdim=c(0),
        compress=compress.annotation)

    # add position
    add.gdsn(gfile, "snp.position", storage="int32", valdim=c(0),
        compress=compress.annotation)

    # add chromosome
    add.gdsn(gfile, "snp.chromosome", storage="string", valdim=c(0),
        compress=compress.annotation)

    # add allele
    add.gdsn(gfile, "snp.allele", storage="string", valdim=c(0),
        compress=compress.annotation)

    # add SNP genotypes
    cmp <- if (snpfirstdim) "" else compress.geno
    nodegeno <- add.gdsn(gfile, "genotype", storage="bit2",
        valdim=c(length(samp.id), 0), compress=cmp)
    put.attr.gdsn(nodegeno, "sample.order")

    # add a folder for SNP annotation
    varAnnot <- add.gdsn(gfile, "snp.annot", storage="folder")
    add.gdsn(varAnnot, "qual", storage="float", valdim=c(0),
        compress=compress.annotation)
    add.gdsn(varAnnot, "filter", storage="string", valdim=c(0),
        compress=compress.annotation)


    ##################################################
    # initialize the internal data
    .Call(gnrParseVCF4Init)

    ##################################################
    # for-loop each file
    for (i in 1:length(vcf.fn))
    {
        opfile <- .OpenConnText(vcf.fn[i])
        on.exit({ closefn.gds(gfile); .CloseConnection(opfile) })

        if (verbose)
            cat("Parsing \"", vcf.fn[i], "\" ...\n", sep="")

        # call C function
        n <- .Call(gnrParseVCF4, vcf.fn[i], gfile$root, metidx,
            readLines, opfile$con, 1024L,  # "readLines(opfile$con, 1024L)"
            ref.allele, ignore.chr.prefix, new.env(), verbose)

        if (verbose)
        {
            if (n > 1)
                cat(sprintf("\timport %d variants.\n", n))
            else
                cat(sprintf("\timport %d variant.\n", n))
            print(nodegeno)
        }

        on.exit({ closefn.gds(gfile) })
        .CloseConnection(opfile)
    }

    on.exit()
    closefn.gds(gfile)

    if (snpfirstdim)
    {
        snpgdsTranspose(out.fn, snpfirstdim=TRUE, compress=compress.geno,
            optimize=FALSE, verbose=verbose)
    }


    ##################################################
    # optimize access efficiency

    if (verbose)
        cat("Optimize the access efficiency ...\n")
    cleanup.gds(out.fn, verbose=verbose)

    # output
    invisible(normalizePath(out.fn))
}



#######################################################################
# Convert a VCF (sequence) file to a GDS file (extract SNP data)
#
# INPUT:
#   vcf.fn -- the file name of VCF format
#   outfn.gds -- the output gds file
#   nblock -- the number of lines in buffer
#   method -- biallelic SNPs, or copy number of variants
#   compress.annotation -- the compression method for sample and snp annotations
#   verbose -- show information
#

snpgdsVCF2GDS_R <- function(vcf.fn, out.fn, nblock=1024,
    method = c("biallelic.only", "copy.num.of.ref"),
    compress.annotation="LZMA_RA", snpfirstdim=FALSE, option=NULL,
    verbose=TRUE)
{
    # check
    stopifnot(is.character(vcf.fn))
    stopifnot(is.character(out.fn))
    stopifnot(is.logical(snpfirstdim) & (length(snpfirstdim)==1L))

    method <- match.arg(method)
    if (!is.null(option) & !is.list(option))
        stop("'option' should be NULL or an object returned by 'snpgdsOption'.")


    ######################################################################
    # Scan VCF file -- get sample id

    scan.vcf.sampid <- function(fn)
    {
        # open the vcf file
        opfile <- file(fn, open="r")

        # read header
        fmtstr <- substring(readLines(opfile, n=1), 3)
        samp.id <- NULL
        while (length(s <- readLines(opfile, n=1)) > 0)
        {
            if (substr(s, 1, 6) == "#CHROM")
            {
                samp.id <- scan(text=s, what=character(0), sep="\t", quiet=TRUE)[-c(1:9)]
                break
            }
        }
        if (is.null(samp.id))
        {
            close(opfile)
            stop("Error VCF format: invalid sample id!")
        }

        # close the file
        close(opfile)

        return(samp.id)
    }


    ######################################################################
    # Scan VCF file -- get marker information

    scan.vcf.marker <- function(fn, method)
    {
        if (verbose)
            cat(sprintf("\tfile: %s\n", fn))

        # total number of rows and columns
        Cnt <- count.fields(fn, sep="\t")
        # check
        if (any(Cnt != Cnt[1]))
            stop(sprintf("The file (%s) has different numbers of columns.", fn))

        line.cnt <- length(Cnt)
        col.cnt <- max(Cnt)
        if (verbose)
            cat(sprintf("\tcontent: %d rows x %d columns\n", line.cnt, col.cnt))

        # open the vcf file
        opfile <- file(fn, open="r")

        # read header
        fmtstr <- substring(readLines(opfile, n=1), 3)
        while (length(s <- readLines(opfile, n=1)) > 0)
        {
            if (substr(s, 1, 6) == "#CHROM")
                break
        }

        # init ...
        chr <- character(line.cnt); position <- integer(line.cnt)
        snpidx <- integer(line.cnt); snp.rs <- character(line.cnt)
        snp.allele <- character(line.cnt)
        snp.cnt <- 0; var.cnt <- 0

        if (method == "biallelic.only")
        {
            while (length(s <- readLines(opfile, n=nblock)) > 0)
            {
                for (i in 1:length(s))
                {
                    var.cnt <- var.cnt + 1
                    ss <- scan(text=s[i], what=character(0), sep="\t", quiet=TRUE, n=5)
                    if (all(ss[c(4,5)] %in% c("A", "G", "C", "T", "a", "g", "c", "t")))
                    {
                        snp.cnt <- snp.cnt + 1
                        chr[snp.cnt] <- ss[1]
                        position[snp.cnt] <- as.integer(ss[2])
                        snpidx[snp.cnt] <- var.cnt
                        snp.rs[snp.cnt] <- ss[3]
                        snp.allele[snp.cnt] <- paste(ss[4], ss[5], sep="/")
                    }
                }
            }
        } else {
            while (length(s <- readLines(opfile, n=nblock)) > 0)
            {
                for (i in 1:length(s))
                {
                    var.cnt <- var.cnt + 1
                    ss <- scan(text=s[i], what=character(0), sep="\t", quiet=TRUE, n=5)
                    snp.cnt <- snp.cnt + 1
                    chr[snp.cnt] <- ss[1]
                    position[snp.cnt] <- as.integer(ss[2])
                    snpidx[snp.cnt] <- var.cnt
                    snp.rs[snp.cnt] <- ss[3]
                    snp.allele[snp.cnt] <- paste(ss[4], ss[5], sep="/")
                }
            }
        }

        # close the file
        close(opfile)


        # chromosomes
        chr <- chr[1:snp.cnt]
        if (!is.null(option))
        {
            flag <- match(chr, names(option$chromosome.code))
            chr[!is.na(flag)] <- unlist(option$chromosome.code)[ flag[!is.na(flag)] ]
            chr <- suppressWarnings(as.integer(chr))
            chr[is.na(chr)] <- -1
        }

        snp.allele <- gsub(".", "/", snp.allele[1:snp.cnt], fixed=TRUE)
        list(chr = chr, position = position[1:snp.cnt],
            snpidx = snpidx[1:snp.cnt], snp.rs = snp.rs[1:snp.cnt],
            snp.allele = snp.allele
        )
    }


    ######################################################################
    # Scan VCF file -- get marker information

    scan.vcf.geno <- function(fn, gGeno, method, start)
    {
        # matching codes
        geno.str <- c("0|0", "0|1", "1|0", "1|1", "0/0", "0/1", "1/0", "1/1",
            "0", "1",
            "0|0|0", "0|0|1", "0|1|0", "0|1|1", "1|0|0", "1|0|1", "1|1|0", "1|1|1",
            "0/0/0", "0/0/1", "0/1/0", "0/1/1", "1/0/0", "1/0/1", "1/1/0", "1/1/1")
        geno.code <- as.integer(c(2, 1, 1, 0, 2, 1, 1, 0,
            1, 0,
            2, 1, 1, 1, 1, 1, 1, 0, 
            2, 1, 1, 1, 1, 1, 1, 0))

        # open the vcf file
        opfile <- file(fn, open="r")
        # read header
        fmtstr <- substring(readLines(opfile, n=1), 3)
        while (length(s <- readLines(opfile, n=1)) > 0)
        {
            if (substr(s, 1, 6) == "#CHROM")
                break
        }

        # scan
        snp.cnt <- start

        if (method == "biallelic.only")
        {
            while (length(s <- readLines(opfile, n=nblock)) > 0)
            {
                gx <- NULL
                for (i in 1:length(s))
                {
                    ss <- scan(text=s[i], what=character(0), sep="\t", quiet=TRUE, n=5)
                    if (all(ss[c(4,5)] %in% c("A", "G", "C", "T", "a", "g", "c", "t")))
                    {
                        ss <- scan(text=s[i], what=character(0), sep="\t", quiet=TRUE)[-c(1:9)]
                        ss <- sapply(strsplit(ss, ":"), FUN = function(x) x[1])
                        x <- match(ss, geno.str)
                        x <- geno.code[x]
                        x[is.na(x)] <- as.integer(3)
                        gx <- cbind(gx, x)
                    }
                }
                if (!is.null(gx))
                {
                    if (snpfirstdim)
                        write.gdsn(gGeno, t(gx), start=c(snp.cnt,1), count=c(ncol(gx),-1))
                    else {
                        print(snp.cnt)
                        write.gdsn(gGeno, gx, start=c(1,snp.cnt), count=c(-1,ncol(gx)))
                    }
                    snp.cnt <- snp.cnt + ncol(gx)
                }
            }
        } else {
            while (length(s <- readLines(opfile, n=nblock)) > 0)
            {
                gx <- NULL
                for (i in 1:length(s))
                {
                    ss <- scan(text=s[i], what=character(0), sep="\t", quiet=TRUE)[-c(1:9)]
                    x <- sapply(strsplit(ss, ":"), FUN = function(x) {
                        a <- unlist(strsplit(x[1], ""))
                        if (any(a == "."))
                            NA
                        else
                            sum(a == "0")
                    })
                    x[x > 2] <- 2
                    x[is.na(x)] <- as.integer(3)
                    gx <- cbind(gx, x)
                }
                if (!is.null(gx))
                {
                    if (snpfirstdim)
                        write.gdsn(gGeno, t(gx), start=c(snp.cnt,1), count=c(ncol(gx),-1))
                    else
                        write.gdsn(gGeno, gx, start=c(1,snp.cnt), count=c(-1,ncol(gx)))
                    snp.cnt <- snp.cnt + ncol(gx)
                }
            }
        }

        # close the file
        close(opfile)

        snp.cnt - start
    }



    ######################################################################
    # Starting ...
    ######################################################################

    if (verbose)
    {
        cat("Start snpgdsVCF2GDS ...\n")
        if (method == "biallelic.only")
            cat("\tExtracting bi-allelic and polymorhpic SNPs.\n")
        else
            cat("\tStoring dosage of the reference allele for all variant sites, including bi-allelic SNPs, multi-allelic SNPs, indels and structural variants.\n")
        cat("\tScanning ...\n")
    }


    ####################################
    # sample.id

    sample.id <- NULL
    for (fn in vcf.fn)
    {
        s <- scan.vcf.sampid(fn)
        if (!is.null(sample.id))
        {
            if (length(sample.id) != length(s))
                stop("All VCF files should have the same sample id.")
            if (any(sample.id != s))
                stop("All VCF files should have the same sample id.")
        } else
            sample.id <- s
    }


    ####################################
    # genetic markers

    all.chr <- NULL
    all.position <- integer()
    all.snpidx <- integer()
    all.snp.rs <- character()
    all.snp.allele <- character()

    for (fn in vcf.fn)
    {
        v <- scan.vcf.marker(fn, method)

        all.chr <- c(all.chr, v$chr)
        all.position <- c(all.position, v$position)
        all.snpidx <- c(all.snpidx, length(all.snpidx) + v$snpidx)
        all.snp.rs <- c(all.snp.rs, v$snp.rs)
        all.snp.allele <- c(all.snp.allele, v$snp.allele)
    }


    ####################################
    # genetic variants

    nSamp <- length(sample.id)
    nSNP <- length(all.chr)
    if (verbose)
    {
        cat(date(), "\tstore sample id, snp id, position, and chromosome.\n")
        cat(sprintf("\tstart writing: %d samples, %d SNPs ...\n", nSamp, nSNP))
    }


    ######################################################################
    # create GDS file
    #
    gfile <- createfn.gds(out.fn)
    on.exit({ closefn.gds(gfile) })

    # add file flag
    put.attr.gdsn(gfile$root, "FileFormat", "SNP_ARRAY")
    put.attr.gdsn(gfile$root, "FileVersion", "v1.0")

    # add "sample.id"
    add.gdsn(gfile, "sample.id", sample.id, compress=compress.annotation, closezip=TRUE)
    # add "snp.id"
    add.gdsn(gfile, "snp.id", as.integer(all.snpidx), compress=compress.annotation, closezip=TRUE)
    # add "snp.rs.id"
    add.gdsn(gfile, "snp.rs.id", all.snp.rs, compress=compress.annotation, closezip=TRUE)
    # add "snp.position"
    add.gdsn(gfile, "snp.position", all.position, compress=compress.annotation, closezip=TRUE)
    # add "snp.chromosome"
    v.chr <- add.gdsn(gfile, "snp.chromosome", all.chr, compress=compress.annotation, closezip=TRUE)
    # add "snp.allele"
    add.gdsn(gfile, "snp.allele", all.snp.allele, compress=compress.annotation, closezip=TRUE)

    # snp.chromosome
    if (!is.null(option))
    {
        put.attr.gdsn(v.chr, "autosome.start", option$autosome.start)
        put.attr.gdsn(v.chr, "autosome.end", option$autosome.end)
        for (i in 1:length(option$chromosome.code))
        {
            put.attr.gdsn(v.chr, names(option$chromosome.code)[i],
                option$chromosome.code[[i]])
        }
    }

    # sync file
    sync.gds(gfile)

    # add "gonetype", 2 bits to store one genotype
    if (snpfirstdim)
    {
        gGeno <- add.gdsn(gfile, "genotype", storage="bit2", valdim=c(nSNP, nSamp))
        put.attr.gdsn(gGeno, "snp.order")
    } else {
        gGeno <- add.gdsn(gfile, "genotype", storage="bit2", valdim=c(nSamp, nSNP))
        put.attr.gdsn(gGeno, "sample.order")
    }
    # sync file
    sync.gds(gfile)


    ####################################
    # genetic genotypes

    snp.start <- 1
    for (fn in vcf.fn)
    {
        if (verbose)
            cat(sprintf("\tfile: %s\n", fn))
        s <- scan.vcf.geno(fn, gGeno, method, start=snp.start)
        snp.start <- snp.start + s      
        sync.gds(gfile)  # sync file
    }

    if (verbose) cat(date(), "\tDone.\n")

    return(invisible(NULL))
}
