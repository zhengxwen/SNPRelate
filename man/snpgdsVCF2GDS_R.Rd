\name{snpgdsVCF2GDS_R}
\alias{snpgdsVCF2GDS_R}
\title{Reformat a VCF file (R implementation)}
\description{
    Reformat a Variant Call Format (VCF) file
}
\usage{
snpgdsVCF2GDS_R(vcf.fn, out.fn, nblock=1024,
    method = c("biallelic.only", "copy.num.of.ref"),
    compress.annotation="ZIP.max", snpfirstdim=FALSE, option = NULL,
    verbose=TRUE)
}
\arguments{
    \item{vcf.fn}{the file name of VCF format, \code{vcf.fn} can be a vector,
        see details}
    \item{out.fn}{the output gds file}
    \item{nblock}{the buffer lines}
    \item{method}{either "biallelic.only" by default or "copy.num.of.ref",
        see details}
    \item{compress.annotation}{the compression flag of the nodes stored, except
        "genotype"; the string value is defined in the function of
        \code{\link[gdsfmt]{add.gdsn}}}
    \item{snpfirstdim}{if TRUE, genotypes are stored in the individual-major
        mode, (i.e, list all SNPs for the first individual, and then list all
        SNPs for the second individual, etc)}
    \item{option}{\code{NULL} or an object from \code{\link{snpgdsOption}},
        see details}
    \item{verbose}{if TRUE, show information}
}
\value{
    None.
}
\details{
    GDS -- Genomic Data Structures used for storing genetic array-oriented
        data, and the file format used in the \link{gdsfmt} package.

    VCF -- The Variant Call Format (VCF), which is a generic format for storing
        DNA polymorphism data such as SNPs, insertions, deletions and
        structural variants, together with rich annotations.

    If there are more than one file name in \code{vcf.fn}, \code{snpgdsVCF2GDS}
will merge all dataset together once they all contain the same samples. It is
useful to combine genetic data if VCF data are divided by chromosomes.

    \code{method = "biallelic.only"}: to exact bi-allelic and polymorhpic SNP
data; \code{method = "copy.num.of.ref"}: to extract and store dosage (0, 1, 2)
of the reference allele for all variant sites, including bi-allelic SNPs,
multi-allelic SNPs, indels and structural variants.

    Haploid and triploid calls are allowed in the transfer, the variable
\code{snp.id} stores the original the row index of variants, and the variable
\code{snp.rs.id} stores the rs id.

    The user could use \code{option} to specify the range of code for autosomes.
For humans there are 22 autosomes (from 1 to 22), but dogs have 38 autosomes.
Note that the default settings are used for humans. The user could call
\code{option = snpgdsOption(autosome.end=38)} for importing the VCF file of dog.
It also allow define new chromosome coding, e.g.,
\code{option = snpgdsOption(Z=27)}.
}
\references{
    The variant call format and VCFtools.
    Danecek P, Auton A, Abecasis G, Albers CA, Banks E, DePristo MA,
    Handsaker RE, Lunter G, Marth GT, Sherry ST, McVean G, Durbin R;
    1000 Genomes Project Analysis Group.
    Bioinformatics. 2011 Aug 1;27(15):2156-8. Epub 2011 Jun 7.
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{snpgdsVCF2GDS_R}}, \code{\link{snpgdsOption}},
    \code{\link{snpgdsBED2GDS}}
}

\examples{
# The VCF file
vcf.fn <- system.file("extdata", "sequence.vcf", package="SNPRelate")
cat(readLines(vcf.fn), sep="\n")

snpgdsVCF2GDS_R(vcf.fn, "test1.gds", method="biallelic.only")
snpgdsSummary("test1.gds")

snpgdsVCF2GDS_R(vcf.fn, "test2.gds", method="biallelic.only")
snpgdsSummary("test2.gds")

snpgdsVCF2GDS_R(vcf.fn, "test3.gds", method="copy.num.of.ref")
snpgdsSummary("test3.gds")

snpgdsVCF2GDS_R(vcf.fn, "test4.gds", method="copy.num.of.ref")
snpgdsSummary("test4.gds")
}

\keyword{gds}
\keyword{GWAS}