// ===========================================================
//
// vec_ext_neon.cpp: AArch64 NEON kernels (SNPRelate v2)
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

// Advanced SIMD is a mandatory baseline on AArch64, so no target attribute or
// runtime feature probe is required; the intrinsics compile with default flags.
// COREARRAY_SIMD_NEON (from CoreDEF.h) == __aarch64__ && __ARM_NEON, and is
// cleared by the global COREARRAY_NO_SIMD switch.
#ifdef COREARRAY_SIMD_NEON

#include <arm_neon.h>

namespace SNPvec
{

// Count IBS0/IBS1/IBS2 for one sample pair over the two bit-planes.
// Processes 16 bytes/plane per iteration.  The bitwise expressions mirror the
// scalar reference exactly, so the counts are bit-identical.
//
// NEON op map (vs. the scalar bit ops):
//   a | ~b           -> vornq_u8(a, b)
//   ~x & m  (== m AND NOT x) -> vbicq_u8(m, x)
//   a ^ ~b  (== ~(a ^ b))    -> vmvnq_u8(veorq_u8(a, b))
//   popcount(byte)   -> vcntq_u8   (then widen-accumulate via vpadalq_u8)
void ibs_count_neon(const uint8_t *gi, const uint8_t *gj, size_t npack,
	uint32_t &ibs0, uint32_t &ibs1, uint32_t &ibs2)
{
	const uint8_t *p1a = gi, *p1b = gi + npack;
	const uint8_t *p2a = gj, *p2b = gj + npack;

	// u32 lane accumulators (overflow-proof for any block size); the inner
	// loop accumulates raw byte popcounts in cheap u8 lanes and widens in
	// bulk -- this keeps the loop-carried chain short (vaddq_u8, not vpadalq).
	uint32x4_t A0 = vdupq_n_u32(0), A2 = vdupq_n_u32(0), AM = vdupq_n_u32(0);
	size_t m = npack;

	while (m >= 16)
	{
		uint8x16_t s0 = vdupq_n_u8(0), s2 = vdupq_n_u8(0), sm = vdupq_n_u8(0);
		// up to 24 iters: each byte lane <= 24*8 = 192 < 256, no u8 overflow
		for (int k = 0; m >= 16 && k < 24; m -= 16, k++)
		{
			const uint8x16_t g1_1 = vld1q_u8(p1a), g1_2 = vld1q_u8(p1b);
			const uint8x16_t g2_1 = vld1q_u8(p2a), g2_2 = vld1q_u8(p2b);
			p1a += 16; p1b += 16; p2a += 16; p2b += 16;

			// e1/e2: per-bit-plane allele disagreement between the two samples
			const uint8x16_t e1 = veorq_u8(g1_1, g2_1);
			const uint8x16_t e2 = veorq_u8(g1_2, g2_2);
			// mask = (g1_1 | ~g1_2) & (g2_1 | ~g2_2)  (both genotypes non-missing)
			const uint8x16_t mask =
				vandq_u8(vornq_u8(g1_1, g1_2), vornq_u8(g2_1, g2_2));
			// ibs0 = (e1 & e2) & mask    [De Morgan of the scalar form]
			const uint8x16_t v0 = vandq_u8(vandq_u8(e1, e2), mask);
			// ibs2 = ~e1 & ~e2 & mask
			const uint8x16_t v2 = vbicq_u8(vbicq_u8(mask, e1), e2);

			s0 = vaddq_u8(s0, vcntq_u8(v0));
			s2 = vaddq_u8(s2, vcntq_u8(v2));
			sm = vaddq_u8(sm, vcntq_u8(mask));
		}
		// widen the u8 partials: u8 ->(pairwise) u16 ->(pairwise-accum) u32
		A0 = vpadalq_u16(A0, vpaddlq_u8(s0));
		A2 = vpadalq_u16(A2, vpaddlq_u8(s2));
		AM = vpadalq_u16(AM, vpaddlq_u8(sm));
	}
	uint64_t c0 = vaddvq_u32(A0), c2 = vaddvq_u32(A2), cm = vaddvq_u32(AM);

	// tail (< 16 bytes): defensive -- npack is a multiple of 16 in practice
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

}  // namespace SNPvec

#endif  /* COREARRAY_SIMD_NEON */
