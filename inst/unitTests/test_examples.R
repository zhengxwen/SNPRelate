#############################################################
#
# DESCRIPTION: run all examples in SNPRelate
#

library(SNPRelate)


test.examples <- function()
{
	function.list <- readRDS(
		system.file("Meta", "Rd.rds", package="SNPRelate"))$Name

	sapply(function.list, FUN = function(func.name)
		{
			args <- list(
				topic=func.name,
				package="SNPRelate",
				echo=FALSE, verbose=FALSE, ask=FALSE
			)
			cat("FUNCTION: ", func.name, "\n", sep="")
			suppressWarnings(do.call(example, args))
			NULL
		})

	invisible()
}
