#############################################################
#
# DESCRIPTION: run all examples in SNPRelate
#

library(SNPRelate)


#############################################################
#
# test functions
#

function.list <- c(
	"SNPGDSFileClass-class",
	"SNPRelate-package",
	"hapmap_geno",
	"snpgdsAdmixProp",
	"snpgdsAlleleSwitch",
	"snpgdsApartSelection",
	"snpgdsBED2GDS",
	"snpgdsClose",
	"snpgdsCombineGeno",
	"snpgdsCreateGeno",
	"snpgdsCreateGenoSet",
	"snpgdsCutTree",
	"snpgdsDiss",
	"snpgdsDrawTree",
	"snpgdsEIGMIX",
	"snpgdsErrMsg",
	"snpgdsExampleFileName",
	"snpgdsFst",
	"snpgdsGDS2BED",
	"snpgdsGDS2Eigen",
	"snpgdsGDS2PED",
	"snpgdsGEN2GDS",
	"snpgdsGRM",
	"snpgdsGetGeno",
	"snpgdsHCluster",
	"snpgdsHWE",
	"snpgdsIBDKING",
	"snpgdsIBDMLE",
	"snpgdsIBDMLELogLik",
	"snpgdsIBDMoM",
	"snpgdsIBDSelection",
	"snpgdsIBS",
	"snpgdsIBSNum",
	"snpgdsIndInb",
	"snpgdsIndInbCoef",
	"snpgdsLDMat",
	"snpgdsLDpair",
	"snpgdsLDpruning",
	"snpgdsOpen",
	"snpgdsOption",
	"snpgdsPCA",
	"snpgdsPCACorr",
	"snpgdsPCASNPLoading",
	"snpgdsPCASampLoading",
	"snpgdsPED2GDS",
	"snpgdsPairIBD",
	"snpgdsPairIBDMLELogLik",
	"snpgdsPairScore",
	"snpgdsSNPList",
	"snpgdsSNPListClass",
	"snpgdsSNPListIntersect",
	"snpgdsSNPListStrand",
	"snpgdsSNPRateFreq",
	"snpgdsSampMissRate",
	"snpgdsSelectSNP",
	"snpgdsSlidingWindow",
	"snpgdsSummary",
	"snpgdsTranspose",
	"snpgdsVCF2GDS",
	"snpgdsVCF2GDS_R"
)


test.examples <- function()
{
	sapply(function.list, FUN = function(func.name)
		{
			args <- list(
				topic=func.name,
				package="SNPRelate",
				echo=FALSE, verbose=FALSE, ask=FALSE
			)
			suppressWarnings(do.call(example, args))
			NULL
		})
	invisible()
}
