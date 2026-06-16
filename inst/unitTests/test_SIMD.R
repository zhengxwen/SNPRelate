#############################################################
#
# DESCRIPTION: test units for the runtime-dispatched SIMD kernels
#
# Exercises every instruction-set path available on the running CPU
# (scalar / SSE2 / AVX2 / AVX-512 on x86, NEON on AArch64) and asserts each
# produces results bit-identical to the scalar reference. On CI this is how the
# x86 SIMD kernels get executed and validated on real silicon.
#

library(RUnit)
library(SNPRelate)


# ISAs that could exist; gnrSetISA() returns the one actually activated, so an
# unsupported request returns a different name and is skipped.
.simd_candidates <- c("scalar", "sse2", "avx2", "avx512", "neon")

.set_isa <- function(isa) .Call("gnrSetISA", isa, PACKAGE="SNPRelate")


test.SIMD.dispatch <- function()
{
	# discover available ISA paths on this machine
	avail <- character(0L)
	for (isa in .simd_candidates)
		if (identical(.set_isa(isa), isa)) avail <- c(avail, isa)
	checkTrue("scalar" %in% avail, "scalar kernel must always be available")

	# open a GDS file
	genofile <- snpgdsOpen(snpgdsExampleFileName(), allow.duplicate=TRUE)
	on.exit({ .set_isa("auto"); snpgdsClose(genofile) })
	samp.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
	sid <- samp.id[1:90]

	# scalar reference (the correctness oracle)
	.set_isa("scalar")
	ref.ibs <- snpgdsIBS(genofile, sample.id=sid, missing.rate=NaN,
		num.thread=1L, verbose=FALSE)
	ref.king <- snpgdsIBDKING(genofile, sample.id=sid, type="KING-robust",
		missing.rate=NaN, num.thread=1L, verbose=FALSE)

	# every other available ISA must match the scalar reference exactly
	for (isa in setdiff(avail, "scalar"))
	{
		checkEquals(.set_isa(isa), isa, paste("activate ISA:", isa))

		v.ibs <- snpgdsIBS(genofile, sample.id=sid, missing.rate=NaN,
			num.thread=1L, verbose=FALSE)
		checkEquals(v.ibs, ref.ibs, paste("IBS dispatch:", isa))

		v.king <- snpgdsIBDKING(genofile, sample.id=sid, type="KING-robust",
			missing.rate=NaN, num.thread=1L, verbose=FALSE)
		checkEquals(v.king, ref.king, paste("KING-robust dispatch:", isa))
	}
}
