#######################################################################
#
# Package name: SNPRelate
#
# Description:
#     A High-performance Computing Toolset for Relatedness and
# Principal Component Analysis of SNP Data
#
# Copyright (C) 2015 - 2017        Xiuwen Zheng
# License: GPL-3
# Email: zhengxwen@gmail.com
#



#######################################################################
# Conversion: bgl.gprobs (BEAGLE) ==> SNP Prob GDS
#

.snpgdsImport <- function(ped.fn, map.fn, out.gdsfn, family=TRUE,
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

    # add a flag
    put.attr.gdsn(gfile$root, "FileFormat", "SNP_ARRAY")

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
