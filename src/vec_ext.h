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
// Modernized, C++17 replacement for the compile-time SIMD layer (dVect.h).
// Each hot kernel is published as a function pointer bound at first use to the
// best implementation for the running CPU.  Architecture-specific kernels live
// in their own translation units (vec_ext_neon.cpp, vec_ext_sse2/avx2/avx512.cpp)
// so each is compiled with its own target options without affecting the rest.
//
//   * one rebindable table (`kernels`) groups the related kernels and records
//     the active ISA name (exported to R via gnrSIMDInfo);
//   * the table is initialized once, lazily and thread-safely (C++11 magic
//     static) to the best available ISA -- no global mutable state, no
//     "remember to call vec_init_function()" footgun;
//   * for tests, set_isa()/SNPRELATE_FORCE_ISA pin a specific path so all ISAs
//     can be exercised and asserted bit-identical in one process.
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


namespace SNPvec
{
	// =====================================================================
	// Packed genotype representation (see GWAS::PackSNPGeno1b)
	//
	// For sample i, bit-plane 1 is gi[0 .. npack-1] and bit-plane 2 is
	// gi[npack .. 2*npack-1] (likewise for sample j).  All ISA paths use the
	// same bit logic, so every one must produce identical integer counts.


	// ---- IBS: count IBS0/IBS1/IBS2 (added to ibs0/ibs1/ibs2) ----
	typedef void (*ibs_count_t)(const uint8_t *gi, const uint8_t *gj,
		size_t npack, uint32_t &ibs0, uint32_t &ibs1, uint32_t &ibs2);


	// ---- KING robust: five counters (added to the fields of TKINGRobust) ----
	struct TKINGRobust
	{
		uint32_t ibs0;   ///< loci sharing no allele
		uint32_t nloci;  ///< total non-missing loci
		uint32_t sumsq;  ///< \sum (g_i - g_j)^2  ( = #het-mismatch + 4*#ibs0 )
		uint32_t n1_aa;  ///< heterozygous loci of sample i
		uint32_t n2_aa;  ///< heterozygous loci of sample j
	};
	typedef void (*king_robust_t)(const uint8_t *gi, const uint8_t *gj,
		size_t npack, TKINGRobust &out);


	// per-ISA implementations (defined in the matching translation unit)
	void ibs_count_def(const uint8_t*, const uint8_t*, size_t,
		uint32_t&, uint32_t&, uint32_t&);
	void king_robust_def(const uint8_t*, const uint8_t*, size_t, TKINGRobust&);

#ifdef COREARRAY_SIMD_NEON
	void ibs_count_neon(const uint8_t*, const uint8_t*, size_t,
		uint32_t&, uint32_t&, uint32_t&);
	void king_robust_neon(const uint8_t*, const uint8_t*, size_t, TKINGRobust&);
#endif

#ifdef SNP_VEC_X86
	void ibs_count_sse2(const uint8_t*, const uint8_t*, size_t,
		uint32_t&, uint32_t&, uint32_t&);
	void king_robust_sse2(const uint8_t*, const uint8_t*, size_t, TKINGRobust&);
#   ifdef SNP_VEC_X86_TARGET
	void ibs_count_avx2(const uint8_t*, const uint8_t*, size_t,
		uint32_t&, uint32_t&, uint32_t&);
	void king_robust_avx2(const uint8_t*, const uint8_t*, size_t, TKINGRobust&);
	void ibs_count_avx512(const uint8_t*, const uint8_t*, size_t,
		uint32_t&, uint32_t&, uint32_t&);
	void king_robust_avx512(const uint8_t*, const uint8_t*, size_t, TKINGRobust&);
#   endif
#endif


	// ---- scalar remainders shared by the SIMD kernels (the < width tail) ----
	// npack is a multiple of 16 in practice, so these usually do nothing; they
	// are pure scalar (no SIMD), safe to compile in any target context.

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

	static inline void king_tail(const uint8_t *a1, const uint8_t *a2,
		const uint8_t *b1, const uint8_t *b2, size_t nb,
		uint64_t &c_ibs0, uint64_t &c_mask, uint64_t &c_het,
		uint64_t &c_aa1, uint64_t &c_aa2)
	{
		for (size_t k=0; k < nb; k++)
		{
			const uint8_t g1_1=a1[k], g1_2=a2[k], g2_1=b1[k], g2_2=b2[k];
			const uint8_t mask = (uint8_t)((g1_1 | ~g1_2) & (g2_1 | ~g2_2));
			const uint8_t e1 = (uint8_t)(g1_1 ^ g2_1), e2 = (uint8_t)(g1_2 ^ g2_2);
			const uint8_t het = (uint8_t)(((g1_1 ^ g1_2) ^ (g2_1 ^ g2_2)) & mask);
			c_ibs0 += __builtin_popcount((uint8_t)(e1 & e2 & mask));
			c_mask += __builtin_popcount(mask);
			c_het  += __builtin_popcount(het);
			c_aa1  += __builtin_popcount((uint8_t)(g1_1 & (uint8_t)~g1_2 & mask));
			c_aa2  += __builtin_popcount((uint8_t)(g2_1 & (uint8_t)~g2_2 & mask));
		}
	}


	// =====================================================================
	// Dispatch table

	/// A bound set of kernels plus the name of the active ISA.
	struct kernels
	{
		ibs_count_t    ibs;
		king_robust_t  king_robust;
		const char    *name;
	};

	/// The active kernel table (resolved once, lazily, to the best ISA).
	const kernels& cpu();

	/// Name of the active ISA path ("avx2", "neon", "scalar", ...); for R/tests.
	const char* active_isa();

	/// Test hook: rebind the table to a specific ISA by name; returns the ISA
	/// actually activated (may differ if unsupported on this CPU). Must not be
	/// called concurrently with a running computation.
	const char* set_isa(const char *name);
}

#endif  /* _HEADER_SNP_VEC_EXT_ */
