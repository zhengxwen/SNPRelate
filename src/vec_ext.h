// ===========================================================
//
// vec_ext.h: runtime-dispatched vectorized kernels (SNPRelate v2)
//
// Copyright (C) 2026    Xiuwen Zheng
//
// This file is part of SNPRelate.
//
// SNPRelate is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License Version 3 as
// published by the Free Software Foundation.
//
// SNPRelate is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with SNPRelate.
// If not, see <http://www.gnu.org/licenses/>.

// ---------------------------------------------------------------------------
// This is the modernized, C++17 replacement for the compile-time SIMD layer
// (dVect.h).  Each hot kernel is published as a function pointer that is bound
// once, at first use, to the best implementation for the running CPU.  The
// architecture-specific kernels live in their own translation units
// (vec_ext_neon.cpp, vec_ext_avx2.cpp, ...) so each can be compiled with its
// own target options without affecting the rest of the package.
//
// Design notes (vs. the 7-year-old SAIGEgds dispatcher this is modeled on):
//   * selection happens through a function-local `static const` table, so it
//     is thread-safe and lazy (C++11 "magic statics") -- no global mutable
//     state and no "remember to call vec_init_function()" footgun;
//   * one struct groups the related kernels and records the active ISA name,
//     which is exported to R for testing/observability (see gnrSIMDInfo);
//   * the override env-var SNPRELATE_FORCE_ISA lets the unit tests pin a path
//     (e.g. "def" vs "neon") and assert bit-identical results.
// ---------------------------------------------------------------------------

#ifndef _HEADER_SNP_VEC_EXT_
#define _HEADER_SNP_VEC_EXT_

#include <cstdint>
#include <cstddef>

// CoreArray's architecture / SIMD macros (COREARRAY_SIMD_NEON, COREARRAY_SIMD_
// SSE2/AVX2/..., COREARRAY_REGISTER_BIT64, the global COREARRAY_NO_SIMD switch).
// Reusing these keeps detection consistent with the rest of the package and
// makes the dispatched kernels honor CoreArray's SIMD policy automatically.
#include <CoreDEF.h>


namespace SNPvec
{
	// =====================================================================
	// IBS pairwise kernel
	//
	// Count IBS0/IBS1/IBS2 between two samples using the 1-bit "two bit-plane"
	// packed genotype representation (see GWAS::PackSNPGeno1b).  For sample i,
	// bit-plane 1 is gi[0 .. npack-1] and bit-plane 2 is gi[npack .. 2*npack-1]
	// (likewise for sample j).  The three counts are *added* to ibs0/ibs1/ibs2.
	//
	// The bit logic is identical across all implementations, so every ISA path
	// must produce exactly the same integer counts as the scalar reference.

	typedef void (*ibs_count_t)(const uint8_t *gi, const uint8_t *gj,
		size_t npack, uint32_t &ibs0, uint32_t &ibs1, uint32_t &ibs2);

	// scalar reference / portable fallback (always available)
	void ibs_count_def(const uint8_t *gi, const uint8_t *gj, size_t npack,
		uint32_t &ibs0, uint32_t &ibs1, uint32_t &ibs2);

#ifdef COREARRAY_SIMD_NEON
	// AArch64 Advanced SIMD (NEON) -- baseline on ARMv8, no runtime probe
	void ibs_count_neon(const uint8_t *gi, const uint8_t *gj, size_t npack,
		uint32_t &ibs0, uint32_t &ibs1, uint32_t &ibs2);
#endif

// x86: runtime dispatch must compile every ISA variant regardless of the
// compile-time -m flags, so it cannot key off the compile-time COREARRAY_SIMD_*
// macros. We use the raw architecture predefined macro (the same primitive
// CoreDEF.h itself keys on) and per-function target attributes (GCC/Clang) so
// AVX2/AVX-512 code is emitted even in a baseline build and selected at runtime.
#if (defined(__x86_64__) || defined(__i386__)) && !defined(COREARRAY_NO_SIMD)
#   define SNP_VEC_X86  1
#   if defined(__GNUC__) || defined(__clang__)
#       define SNP_VEC_X86_TARGET  1   // __attribute__((target(...))) available
#   endif
#endif

#ifdef SNP_VEC_X86
	void ibs_count_sse2(const uint8_t *gi, const uint8_t *gj, size_t npack,
		uint32_t &ibs0, uint32_t &ibs1, uint32_t &ibs2);
#   ifdef SNP_VEC_X86_TARGET
	void ibs_count_avx2(const uint8_t *gi, const uint8_t *gj, size_t npack,
		uint32_t &ibs0, uint32_t &ibs1, uint32_t &ibs2);
	void ibs_count_avx512(const uint8_t *gi, const uint8_t *gj, size_t npack,
		uint32_t &ibs0, uint32_t &ibs1, uint32_t &ibs2);
#   endif

	// Scalar remainder shared by the x86 SIMD kernels (the < SIMD-width tail).
	// npack is a multiple of 16 in practice, so this usually does nothing; it is
	// pure scalar (no SIMD), safe to compile in any target context.
	static inline void ibs_tail(const uint8_t *a1, const uint8_t *a2,
		const uint8_t *b1, const uint8_t *b2, size_t nb,
		uint64_t &c0, uint64_t &c2, uint64_t &cm)
	{
		for (size_t k=0; k < nb; k++)
		{
			const uint8_t g1_1=a1[k], g1_2=a2[k], g2_1=b1[k], g2_2=b2[k];
			const uint8_t mask = (uint8_t)((g1_1 | ~g1_2) & (g2_1 | ~g2_2));
			const uint8_t e1 = (uint8_t)(g1_1 ^ g2_1), e2 = (uint8_t)(g1_2 ^ g2_2);
			c0 += __builtin_popcount((uint8_t)(e1 & e2 & mask));
			c2 += __builtin_popcount((uint8_t)(mask & (uint8_t)~e1 & (uint8_t)~e2));
			cm += __builtin_popcount(mask);
		}
	}
#endif


	/// A bound set of IBS kernels plus the name of the active ISA.
	struct ibs_kernels
	{
		ibs_count_t count;
		const char *name;
	};

	/// Thread-safe, lazy, one-time dispatch (resolved on first call).
	const ibs_kernels& ibs();

	/// Name of the active IBS ISA path ("neon", "scalar", ...); for R/tests.
	const char* active_isa();
}

#endif  /* _HEADER_SNP_VEC_EXT_ */
