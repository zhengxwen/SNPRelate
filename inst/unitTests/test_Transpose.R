#############################################################
#
# DESCRIPTION: test transposing genotypic matrix
#

library(RUnit)
library(SNPRelate)


#############################################################
# test function
#

test.Transpose_SNP_Matrix <- function()
{
    # the file name of SNP GDS
    fn <- snpgdsExampleFileName()

    # copy the file
    file.copy(fn, "test.gds", overwrite=TRUE)

    # SNP matrix
    g1 <- snpgdsGetGeno("test.gds", verbose=FALSE)

    # transpose the SNP matrix
    snpgdsTranspose("test.gds", snpfirstdim=NA, verbose=FALSE)

    # SNP matrix
    g2 <- snpgdsGetGeno("test.gds", verbose=FALSE)

    # check
    checkEquals(g1, t(g2), "transpose the SNP matrix")

    # delete the temporary file
    unlink("test.gds", force=TRUE)
}
