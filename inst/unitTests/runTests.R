# ===========================================================
#     _/_/_/   _/_/_/  _/_/_/_/    _/_/_/_/  _/_/_/   _/_/_/
#      _/    _/       _/             _/    _/    _/   _/   _/
#     _/    _/       _/_/_/_/       _/    _/    _/   _/_/_/
#    _/    _/       _/             _/    _/    _/   _/
# _/_/_/   _/_/_/  _/_/_/_/_/     _/     _/_/_/   _/_/
# ===========================================================
#
# runTests.r: the R interface of CoreArray library
#
# Copyright (C) 2012	Xiuwen Zheng


# load R packages
library(RUnit)
library(SNPRelate)

# define a test suite
myTestSuite <- defineTestSuite("SNPRelate examples",
	testFileRegexp = "^runit.+\\.[rR]$",
	system.file("unitTests", package = "SNPRelate"))

# run the test suite
testResult <- runTestSuite(myTestSuite)

# print detailed text protocol to standard out
printTextProtocol(testResult)

# quit
q("no")
