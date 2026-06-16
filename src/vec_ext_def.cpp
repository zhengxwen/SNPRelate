// ===========================================================
//
// vec_ext_def.cpp: scalar reference kernels + runtime dispatcher
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

#include "vec_ext.h"
#include <cstdlib>
#include <cstring>


namespace SNPvec
{

// ---------------------------------------------------------------------------
// Scalar reference: 64-bit bitwise + popcount, identical logic to the original
// genIBS.cpp inner loop.  This is both the portable fallback and the
// correctness oracle the SIMD paths are validated against.

void ibs_count_def(const uint8_t *gi, const uint8_t *gj, size_t npack,
	uint32_t &ibs0, uint32_t &ibs1, uint32_t &ibs2)
{
	const uint8_t *p1a = gi, *p1b = gi + npack;   // sample i, bit-plane 1 / 2
	const uint8_t *p2a = gj, *p2b = gj + npack;   // sample j, bit-plane 1 / 2
	uint64_t c0 = 0, c2 = 0, cm = 0;
	size_t m = npack;

	for (; m >= 8; m -= 8)
	{
		uint64_t g1_1, g1_2, g2_1, g2_2;
		__builtin_memcpy(&g1_1, p1a, 8);  __builtin_memcpy(&g1_2, p1b, 8);
		__builtin_memcpy(&g2_1, p2a, 8);  __builtin_memcpy(&g2_2, p2b, 8);
		p1a += 8; p1b += 8; p2a += 8; p2b += 8;

		const uint64_t mask = (g1_1 | ~g1_2) & (g2_1 | ~g2_2);
		const uint64_t v0 = (~((g1_1 ^ ~g2_1) | (g1_2 ^ ~g2_2))) & mask;
		const uint64_t v2 = (~((g1_1 ^  g2_1) | (g1_2 ^  g2_2))) & mask;

		c0 += __builtin_popcountll(v0);
		c2 += __builtin_popcountll(v2);
		cm += __builtin_popcountll(mask);
	}
	// tail (< 8 bytes): defensive -- npack is a multiple of 16 in practice
	for (; m > 0; m--)
	{
		const uint8_t g1_1 = *p1a++, g1_2 = *p1b++, g2_1 = *p2a++, g2_2 = *p2b++;
		const uint8_t mask = (uint8_t)((g1_1 | ~g1_2) & (g2_1 | ~g2_2));
		const uint8_t v0 = (uint8_t)((~((g1_1 ^ ~g2_1) | (g1_2 ^ ~g2_2))) & mask);
		const uint8_t v2 = (uint8_t)((~((g1_1 ^  g2_1) | (g1_2 ^  g2_2))) & mask);
		c0 += __builtin_popcount(v0);
		c2 += __builtin_popcount(v2);
		cm += __builtin_popcount(mask);
	}

	ibs0 += (uint32_t)c0;
	ibs2 += (uint32_t)c2;
	ibs1 += (uint32_t)(cm - c0 - c2);
}


// ---------------------------------------------------------------------------
// Runtime dispatcher

static inline bool isa_eq(const char *f, const char *name)
{
	return f && !std::strcmp(f, name);
}

static ibs_kernels select_ibs()
{
	// default: portable scalar
	ibs_kernels k = { &ibs_count_def, "scalar" };

	// optional override for testing / benchmarking (SNPRELATE_FORCE_ISA)
	const char *force = std::getenv("SNPRELATE_FORCE_ISA");
	if (isa_eq(force, "def") || isa_eq(force, "scalar"))
		return k;

#if defined(COREARRAY_SIMD_NEON)
	// AArch64: NEON is a mandatory baseline feature, so no runtime probe.
	if (!force || isa_eq(force, "neon"))
	{
		k.count = &ibs_count_neon;
		k.name  = "neon";
	}

#elif defined(SNP_VEC_X86)
	// x86: pick the best instruction set the running CPU actually supports.
	__builtin_cpu_init();
#  ifdef SNP_VEC_X86_TARGET
	if ((!force || isa_eq(force, "avx512")) &&
		__builtin_cpu_supports("avx512f") && __builtin_cpu_supports("avx512bw"))
	{
		k.count = &ibs_count_avx512; k.name = "avx512";
	}
	else if ((!force || isa_eq(force, "avx2")) && __builtin_cpu_supports("avx2"))
	{
		k.count = &ibs_count_avx2; k.name = "avx2";
	}
	else
#  endif
	if ((!force || isa_eq(force, "sse2")) && __builtin_cpu_supports("sse2"))
	{
		k.count = &ibs_count_sse2; k.name = "sse2";
	}

#else
	(void)force;
#endif

	return k;
}

const ibs_kernels& ibs()
{
	// C++11 magic static: thread-safe, initialized exactly once on first use.
	static const ibs_kernels K = select_ibs();
	return K;
}

const char* active_isa()
{
	return ibs().name;
}

}  // namespace SNPvec
