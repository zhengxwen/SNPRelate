\name{snpgdsPairIBDMLELogLik}
\alias{snpgdsPairIBDMLELogLik}
\title{
	Log likelihood for MLE method in the Identity-By-Descent (IBD) Analysis
}
\description{
	Calculate the log likelihood values from maximum likelihood estimation.
}
\usage{
snpgdsPairIBDMLELogLik(geno1, geno2, allele.freq, k0=NaN, k1=NaN,
	relatedness=c("", "self", "fullsib", "offspring", "halfsib",
	"cousin", "unrelated"), verbose=TRUE)
}
\arguments{
	\item{geno1}{the SNP genotypes for the first individual, 0 -- BB, 1 -- AB, 2 -- AA,
		other values -- missing}
	\item{geno2}{the SNP genotypes for the second individual, 0 -- BB, 1 -- AB, 2 -- AA,
		other values -- missing}
	\item{allele.freq}{the allele frequencies}
	\item{k0}{specified IBD coefficient}
	\item{k1}{specified IBD coefficient}
	\item{relatedness}{specify a relatedness, otherwise use the values of k0 and k1}
	\item{verbose}{if TRUE, show information}
}
\details{
	If (\code{relatedness} == "") and (k0 == NaN or k1 == NaN), then return the log likelihood values
for each (k0, k1) stored in ibdobj. \\
	If (\code{relatedness} == "") and (k0 != NaN) and (k1 != NaN), then return the log likelihood
values for a specific IBD coefficient (k0, k1). \\
	If \code{relatedness} is: "self", then k0 = 0, k1 = 0; "fullsib", then k0 = 0.25, k1 = 0.5;
"offspring", then k0 = 0, k1 = 1; "halfsib", then k0 = 0.5, k1 = 0.5;
"cousin", then k0 = 0.75, k1 = 0.25; "unrelated", then k0 = 1, k1 = 0.
}
\value{
	The value of log likelihood.
}
\references{
	Milligan BG. 2003. Maximum-likelihood estimation of relatedness. Genetics 163:1153-1167.

	Weir BS, Anderson AD, Hepler AB. 2006.
		Genetic relatedness analysis: modern data and new challenges. Nat Rev Genet. 7(10):771-80.

	Choi Y, Wijsman EM, Weir BS. 2009.
		Case-control association testing in the presence of unknown relationships.
		Genet Epidemiol 33(8):668-78.
}
\author{Xiuwen Zheng}
\seealso{
	\code{\link{snpgdsPairIBD}}, \code{\link{snpgdsIBDMLE}},
	\code{\link{snpgdsIBDMLELogLik}}, \code{\link{snpgdsIBDMoM}}
}

\keyword{gds}
\keyword{GWAS}
\keyword{IBD}
