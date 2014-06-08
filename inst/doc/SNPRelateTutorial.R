### R code from vignette source 'SNPRelateTutorial.Rnw'

###################################################
### code chunk number 1: SNPRelateTutorial.Rnw:68-71
###################################################
# load the R packages: gdsfmt and SNPRelate
library(gdsfmt)
library(SNPRelate)


###################################################
### code chunk number 2: SNPRelateTutorial.Rnw:75-76
###################################################
snpgdsSummary(snpgdsExampleFileName())


###################################################
### code chunk number 3: SNPRelateTutorial.Rnw:79-81
###################################################
# open a GDS file
(genofile <- openfn.gds(snpgdsExampleFileName()))


###################################################
### code chunk number 4: SNPRelateTutorial.Rnw:94-96
###################################################
# get the attributes of chromosome coding
get.attr.gdsn(index.gdsn(genofile, "snp.chromosome"))


###################################################
### code chunk number 5: SNPRelateTutorial.Rnw:101-106
###################################################
# Take out genotype data for the first 3 samples and the first 5 SNPs
(g <- read.gdsn(index.gdsn(genofile, "genotype"), start=c(1,1), count=c(5,3)))

# Get the attribute of genotype
get.attr.gdsn(index.gdsn(genofile, "genotype"))


###################################################
### code chunk number 6: SNPRelateTutorial.Rnw:110-114
###################################################
# Take out snp.id
head(read.gdsn(index.gdsn(genofile, "snp.id")))
# Take out snp.rs.id
head(read.gdsn(index.gdsn(genofile, "snp.rs.id")))


###################################################
### code chunk number 7: SNPRelateTutorial.Rnw:123-129
###################################################
# Read population information
pop <- read.gdsn(index.gdsn(genofile, path="sample.annot/pop.group"))
table(pop)

# close the GDS file
closefn.gds(genofile)


###################################################
### code chunk number 8: SNPRelateTutorial.Rnw:144-159
###################################################
# load data
data(hapmap_geno)

# create a gds file
snpgdsCreateGeno("test.gds", genmat = hapmap_geno$genotype,
    sample.id = hapmap_geno$sample.id, snp.id = hapmap_geno$snp.id,
    snp.chromosome = hapmap_geno$snp.chromosome,
    snp.position = hapmap_geno$snp.position,
    snp.allele = hapmap_geno$snp.allele, snpfirstdim=TRUE)

# open the gds file
(genofile <- openfn.gds("test.gds"))

# close the genotype file
closefn.gds(genofile)


###################################################
### code chunk number 9: SNPRelateTutorial.Rnw:249-253
###################################################
# the PLINK BED file, using the example in the SNPRelate package
bed.fn <- system.file("extdata", "plinkhapmap.bed", package="SNPRelate")
bim.fn <- system.file("extdata", "plinkhapmap.bim", package="SNPRelate")
fam.fn <- system.file("extdata", "plinkhapmap.fam", package="SNPRelate")


###################################################
### code chunk number 10: SNPRelateTutorial.Rnw:256-259 (eval = FALSE)
###################################################
## bed.fn <- "C:/your_folder/your_plink_file.bed"
## bim.fn <- "C:/your_folder/your_plink_file.bim"
## fam.fn <- "C:/your_folder/your_plink_file.fam"


###################################################
### code chunk number 11: SNPRelateTutorial.Rnw:261-266
###################################################
# convert
snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, "test.gds")

# summary
snpgdsSummary("test.gds")


###################################################
### code chunk number 12: SNPRelateTutorial.Rnw:275-277
###################################################
# the VCF file, using the example in the SNPRelate package
vcf.fn <- system.file("extdata", "sequence.vcf", package="SNPRelate")


###################################################
### code chunk number 13: SNPRelateTutorial.Rnw:280-281 (eval = FALSE)
###################################################
## vcf.fn <- "C:/your_folder/your_vcf_file.vcf"


###################################################
### code chunk number 14: SNPRelateTutorial.Rnw:283-288
###################################################
# reformat
snpgdsVCF2GDS(vcf.fn, "test.gds", method="biallelic.only")

# summary
snpgdsSummary("test.gds")


###################################################
### code chunk number 15: SNPRelateTutorial.Rnw:298-300
###################################################
# open the GDS file
genofile <- openfn.gds(snpgdsExampleFileName())


###################################################
### code chunk number 16: SNPRelateTutorial.Rnw:303-309
###################################################
# get population information
#   or pop_code <- scan("pop.txt", what=character()), if it is stored in a text file "pop.txt"
pop_code <- read.gdsn(index.gdsn(genofile, path="sample.annot/pop.group"))

# display the first six values
head(pop_code)


###################################################
### code chunk number 17: SNPRelateTutorial.Rnw:318-327
###################################################
set.seed(1000)

# try different LD thresholds for sensitivity analysis
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2)
names(snpset)
head(snpset$chr1)  # snp.id

# get all selected snp id
snpset.id <- unlist(snpset)


###################################################
### code chunk number 18: SNPRelateTutorial.Rnw:336-342
###################################################
# get sample id
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

# get population information
#   or pop_code <- scan("pop.txt", what=character()), if it is stored in a text file "pop.txt"
pop_code <- read.gdsn(index.gdsn(genofile, "sample.annot/pop.group"))


###################################################
### code chunk number 19: SNPRelateTutorial.Rnw:345-360
###################################################
# run PCA
pca <- snpgdsPCA(genofile)

# make a data.frame
tab <- data.frame(sample.id = pca$sample.id,
    pop = factor(pop_code)[match(pca$sample.id, sample.id)],
    EV1 = pca$eigenvect[,1],    # the first eigenvector
    EV2 = pca$eigenvect[,2],    # the second eigenvector
    stringsAsFactors = FALSE)
head(tab)

# draw
plot(tab$EV2, tab$EV1, col=as.integer(tab$pop),
    xlab="eigenvector 2", ylab="eigenvector 1")
legend("topleft", legend=levels(tab$pop), pch="o", col=1:nlevels(tab$pop))


###################################################
### code chunk number 20: SNPRelateTutorial.Rnw:364-366
###################################################
pc.percent <- 100 * pca$eigenval[1:16]/sum(pca$eigenval)
pc.percent


###################################################
### code chunk number 21: SNPRelateTutorial.Rnw:370-372
###################################################
lbls <- paste("PC", 1:4, "\n", format(pc.percent[1:4], digits=2), "%", sep="")
pairs(pca$eigenvect[,1:4], col=tab$pop, labels=lbls)


###################################################
### code chunk number 22: SNPRelateTutorial.Rnw:376-386
###################################################
# get chromosome index
chr <- read.gdsn(index.gdsn(genofile, "snp.chromosome"))
CORR <- snpgdsPCACorr(pca, genofile, eig.which=1:4)

par( mfrow=c(3,1))
for (i in 1:3)
{
    plot(abs(CORR$snpcorr[i,]), ylim=c(0,1), xlab="SNP Index",
        ylab=paste("PC", i), col=chr, pch="+")
}


###################################################
### code chunk number 23: SNPRelateTutorial.Rnw:396-399
###################################################
# YRI samples
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
YRI.id <- sample.id[pop_code == "YRI"]


###################################################
### code chunk number 24: SNPRelateTutorial.Rnw:404-411
###################################################
# estimate IBD coefficients
ibd <- snpgdsIBDMoM(genofile, sample.id=YRI.id, snp.id=snpset.id,
    maf=0.05, missing.rate=0.05)

# make a data.frame
ibd.coeff <- snpgdsIBDSelection(ibd)
head(ibd.coeff)


###################################################
### code chunk number 25: SNPRelateTutorial.Rnw:414-417
###################################################
plot(ibd.coeff$k0, ibd.coeff$k1, xlim=c(0,1), ylim=c(0,1),
    xlab="k0", ylab="k1", main="YRI samples (MoM)")
lines(c(0,1), c(1,0), col="red", lty=2)


###################################################
### code chunk number 26: SNPRelateTutorial.Rnw:422-427 (eval = FALSE)
###################################################
## # estimate IBD coefficients
## set.seed(1000)
## snp.id <- sample(snpset.id, 5000)  # random 5000 SNPs
## ibd <- snpgdsIBDMLE(genofile, sample.id=YRI.id, snp.id=snp.id,
##     maf=0.05, missing.rate=0.05)


###################################################
### code chunk number 27: SNPRelateTutorial.Rnw:430-432 (eval = FALSE)
###################################################
## # make a data.frame
## ibd.coeff <- snpgdsIBDSelection(ibd)


###################################################
### code chunk number 28: SNPRelateTutorial.Rnw:435-438 (eval = FALSE)
###################################################
## plot(ibd.coeff$k0, ibd.coeff$k1, xlim=c(0,1), ylim=c(0,1),
##     xlab="k0", ylab="k1", main="YRI samples (MLE)")
## lines(c(0,1), c(1,0), col="red", lty=2)


###################################################
### code chunk number 29: SNPRelateTutorial.Rnw:446-457
###################################################
# incorporate with pedigree information
family.id <- read.gdsn(index.gdsn(genofile, "sample.annot/family.id"))
family.id <- family.id[match(YRI.id, sample.id)]
table(family.id)

ibd.robust <- snpgdsIBDKING(genofile, sample.id=YRI.id, family.id=family.id)
names(ibd.robust)

# pairs of individuals
dat <- snpgdsIBDSelection(ibd.robust)
head(dat)


###################################################
### code chunk number 30: SNPRelateTutorial.Rnw:460-462
###################################################
plot(dat$IBS0, dat$kinship, xlab="Proportion of Zero IBS",
    ylab="Estimated Kinship Coefficient (KING-robust)")


###################################################
### code chunk number 31: SNPRelateTutorial.Rnw:472-473
###################################################
ibs <- snpgdsIBS(genofile, num.thread=2)


###################################################
### code chunk number 32: SNPRelateTutorial.Rnw:476-480
###################################################
library(lattice)

L <- order(pop_code)
levelplot(ibs$ibs[L, L], col.regions = terrain.colors)


###################################################
### code chunk number 33: SNPRelateTutorial.Rnw:484-487
###################################################
loc <- cmdscale(1 - ibs$ibs, k = 2)
x <- loc[, 1]; y <- loc[, 2]
race <- as.factor(pop_code)


###################################################
### code chunk number 34: SNPRelateTutorial.Rnw:490-493
###################################################
plot(x, y, col=race, xlab = "", ylab = "",
	main = "Multidimensional Scaling Analysis (IBS Distance)")
legend("topleft", legend=levels(race), text.col=1:nlevels(race))


###################################################
### code chunk number 35: SNPRelateTutorial.Rnw:498-505
###################################################
set.seed(100)
ibs.hc <- snpgdsHCluster(snpgdsIBS(genofile, num.thread=2))

# to determine groups of individuals automatically
rv <- snpgdsCutTree(ibs.hc)
plot(rv$dendrogram, leaflab="none", main="HapMap Phase II")
table(rv$samp.group)


###################################################
### code chunk number 36: SNPRelateTutorial.Rnw:509-511
###################################################
# to determine groups of individuals by population information
rv2 <- snpgdsCutTree(ibs.hc, samp.group=as.factor(pop_code))


###################################################
### code chunk number 37: SNPRelateTutorial.Rnw:514-516
###################################################
plot(rv2$dendrogram, leaflab="none", main="HapMap Phase II")
legend("topright", legend=levels(race), col=1:nlevels(race), pch=19, ncol=4)


###################################################
### code chunk number 38: SNPRelateTutorial.Rnw:519-521
###################################################
# close the GDS file
closefn.gds(genofile)


