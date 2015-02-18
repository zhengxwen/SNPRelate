\name{snpgdsVCF2GDS}
\alias{snpgdsVCF2GDS}
\title{Reformat VCF file(s)}
\description{
    Reformat Variant Call Format (VCF) file(s)
}
\usage{
snpgdsVCF2GDS(vcf.fn, out.fn, method=c("biallelic.only", "copy.num.of.ref"),
    snpfirstdim=FALSE, compress.annotation="ZIP_RA.max", compress.geno="",
    ref.allele=NULL, ignore.chr.prefix="chr", verbose=TRUE)
}
\arguments{
    \item{vcf.fn}{the file name of VCF format, \code{vcf.fn} can be a vector,
        see details}
    \item{out.fn}{the file name of output GDS}
    \item{method}{either "biallelic.only" by default or "copy.num.of.ref",
        see details}
    \item{snpfirstdim}{if TRUE, genotypes are stored in the individual-major
        mode, (i.e, list all SNPs for the first individual, and then list all
        SNPs for the second individual, etc)}
    \item{compress.annotation}{the data compression mode, except the GDS
        variable "genotype"; optional values are defined in the function of
        \code{add.gdsn}}
    \item{compress.geno}{the data compression mode for the variable "genotype",
        optional values are defined in the function of \code{add.gdsn}}
    \item{ref.allele}{\code{NULL} or a character vector indicating reference
        allele (like \code{"A", "G", "T", NA, ...}) for each site where
        \code{NA} to use the original reference allele in the VCF file(s).
        The length of character vector should be the total number of variants
        in the VCF file(s).}
    \item{ignore.chr.prefix}{a vector of character, indicating the prefix of
        chromosome which should be ignored, like "chr"; it is not
        case-sensitive}
    \item{verbose}{if \code{TRUE}, show information}
}
\value{
    Return the file name of GDS format with an absolute path.
}
\details{
    GDS -- Genomic Data Structures used for storing genetic array-oriented
        data, and the file format used in the \link{gdsfmt} package.

    VCF -- The Variant Call Format (VCF), which is a generic format for
        storing DNA polymorphism data such as SNPs, insertions, deletions
        and structural variants, together with rich annotations.

    If there are more than one file names in \code{vcf.fn},
\code{snpgdsVCF2GDS} will merge all dataset together if they all contain
the same samples. It is useful to combine genetic/genomic data together if
VCF data are divided by chromosomes.

    \code{method = "biallelic.only"}: to exact bi-allelic and polymorhpic
SNP data (excluding monomorphic variants);
    \code{method = "copy.num.of.ref"}: to extract and store dosage (0, 1, 2)
of the reference allele for all variant sites, including bi-allelic SNPs,
multi-allelic SNPs, indels and structural variants.

    Haploid and triploid calls are allowed in the transfer, the variable
\code{snp.id} stores the original the row index of variants, and the
variable \code{snp.rs.id} stores the rs id.

	When \code{snp.chromosome} in the GDS file is character, SNPRelate treats
a chromosome as autosome only if it can be converted to a numeric value (
like "1", "22"). It uses "X" and "Y" for non-autosomes instead of numeric
codes. However, some software format chromosomes in VCF files with a prefix
"chr". Users should remove that prefix when importing VCF files by setting
\code{ignore.chr.prefix = "chr"}.
}
\references{
    The variant call format and VCFtools.
    Danecek P, Auton A, Abecasis G, Albers CA, Banks E, DePristo MA,
    Handsaker RE, Lunter G, Marth GT, Sherry ST, McVean G, Durbin R;
    1000 Genomes Project Analysis Group.
    Bioinformatics. 2011 Aug 1;27(15):2156-8. Epub 2011 Jun 7.

    \url{http://corearray.sourceforge.net/}
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{snpgdsBED2GDS}}
}

\examples{
# the VCF file
vcf.fn <- system.file("extdata", "sequence.vcf", package="SNPRelate")
cat(readLines(vcf.fn), sep="\n")

snpgdsVCF2GDS(vcf.fn, "test1.gds", method="biallelic.only")
snpgdsSummary("test1.gds")

snpgdsVCF2GDS(vcf.fn, "test2.gds", method="biallelic.only", snpfirstdim=TRUE)
snpgdsSummary("test2.gds")

snpgdsVCF2GDS(vcf.fn, "test3.gds", method="copy.num.of.ref", snpfirstdim=TRUE)
snpgdsSummary("test3.gds")

snpgdsVCF2GDS(vcf.fn, "test4.gds", method="copy.num.of.ref")
snpgdsSummary("test4.gds")

snpgdsVCF2GDS(vcf.fn, "test5.gds", method="copy.num.of.ref",
    ref.allele=c("A", "T", "T", "T", "A"))
snpgdsSummary("test5.gds")



# open "test1.gds"
(genofile <- snpgdsOpen("test1.gds"))

read.gdsn(index.gdsn(genofile, "sample.id"))
read.gdsn(index.gdsn(genofile, "snp.rs.id"))
read.gdsn(index.gdsn(genofile, "genotype"))

# close the file
snpgdsClose(genofile)


# open "test2.gds"
(genofile <- snpgdsOpen("test2.gds"))

read.gdsn(index.gdsn(genofile, "sample.id"))
read.gdsn(index.gdsn(genofile, "snp.rs.id"))
read.gdsn(index.gdsn(genofile, "genotype"))

# close the file
snpgdsClose(genofile)


# open "test3.gds"
(genofile <- snpgdsOpen("test3.gds"))

read.gdsn(index.gdsn(genofile, "sample.id"))
read.gdsn(index.gdsn(genofile, "snp.rs.id"))
read.gdsn(index.gdsn(genofile, "genotype"))

# close the file
snpgdsClose(genofile)


# open "test4.gds"
(genofile <- snpgdsOpen("test4.gds"))

read.gdsn(index.gdsn(genofile, "sample.id"))
read.gdsn(index.gdsn(genofile, "snp.rs.id"))
read.gdsn(index.gdsn(genofile, "snp.allele"))
read.gdsn(index.gdsn(genofile, "genotype"))

# close the file
snpgdsClose(genofile)


# open "test5.gds"
(genofile <- snpgdsOpen("test5.gds"))

read.gdsn(index.gdsn(genofile, "sample.id"))
read.gdsn(index.gdsn(genofile, "snp.rs.id"))
read.gdsn(index.gdsn(genofile, "snp.allele"))
read.gdsn(index.gdsn(genofile, "genotype"))

# close the file
snpgdsClose(genofile)


# delete the temporary files
unlink(paste("test", 1:5, ".gds", sep=""), force=TRUE)
}

\keyword{GDS}
\keyword{GWAS}
