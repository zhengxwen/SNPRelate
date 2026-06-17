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

// NEON op map (vs. the scalar bit ops):
//   a | ~b    -> vornq_u8(a, b)        a & ~b -> vbicq_u8(a, b)
//   ~x & m    -> vbicq_u8(m, x)        popcount(byte) -> vcntq_u8
// All bit expressions mirror the scalar reference, so counts are bit-identical.

// ---- IBS: 16 bytes/plane per iteration ----
void ibs_count_neon(const uint8_t *gi, const uint8_t *gj, size_t npack,
	uint32_t &ibs0, uint32_t &ibs1, uint32_t &ibs2)
{
	const uint8_t *p1a=gi, *p1b=gi+npack, *p2a=gj, *p2b=gj+npack;
	uint16x8_t a0=vdupq_n_u16(0), a2=vdupq_n_u16(0), am=vdupq_n_u16(0);
	uint64_t c0=0, c2=0, cm=0;
	size_t m=npack, iter=0;
	for (; m >= 16; m -= 16)
	{
		const uint8x16_t g1_1=vld1q_u8(p1a), g1_2=vld1q_u8(p1b);
		const uint8x16_t g2_1=vld1q_u8(p2a), g2_2=vld1q_u8(p2b);
		p1a+=16; p1b+=16; p2a+=16; p2b+=16;
		const uint8x16_t e1=veorq_u8(g1_1,g2_1), e2=veorq_u8(g1_2,g2_2);
		const uint8x16_t mask=vandq_u8(vornq_u8(g1_1,g1_2), vornq_u8(g2_1,g2_2));
		const uint8x16_t v0=vandq_u8(vandq_u8(e1,e2), mask);   // ibs0
		const uint8x16_t v2=vbicq_u8(vbicq_u8(mask,e1), e2);   // ibs2 = mask&~e1&~e2
		a0=vpadalq_u8(a0, vcntq_u8(v0));
		a2=vpadalq_u8(a2, vcntq_u8(v2));
		am=vpadalq_u8(am, vcntq_u8(mask));
		if (++iter == 4000)
		{
			c0+=vaddlvq_u16(a0); c2+=vaddlvq_u16(a2); cm+=vaddlvq_u16(am);
			a0=vdupq_n_u16(0); a2=vdupq_n_u16(0); am=vdupq_n_u16(0); iter=0;
		}
	}
	c0+=vaddlvq_u16(a0); c2+=vaddlvq_u16(a2); cm+=vaddlvq_u16(am);
	ibs_tail(p1a, p1b, p2a, p2b, m, c0, c2, cm);
	ibs0 += (uint32_t)c0; ibs2 += (uint32_t)c2; ibs1 += (uint32_t)(cm-c0-c2);
}

// ---- KING robust: five counters, 16 bytes/plane per iteration ----
void king_robust_neon(const uint8_t *gi, const uint8_t *gj, size_t npack,
	TKINGRobust &out)
{
	const uint8_t *p1a=gi, *p1b=gi+npack, *p2a=gj, *p2b=gj+npack;
	uint16x8_t Ai=vdupq_n_u16(0), Am=vdupq_n_u16(0), Ah=vdupq_n_u16(0),
		A1=vdupq_n_u16(0), A2=vdupq_n_u16(0);
	uint64_t c_ibs0=0, c_mask=0, c_het=0, c_aa1=0, c_aa2=0;
	size_t m=npack, iter=0;
	for (; m >= 16; m -= 16)
	{
		const uint8x16_t g1_1=vld1q_u8(p1a), g1_2=vld1q_u8(p1b);
		const uint8x16_t g2_1=vld1q_u8(p2a), g2_2=vld1q_u8(p2b);
		p1a+=16; p1b+=16; p2a+=16; p2b+=16;
		const uint8x16_t e1=veorq_u8(g1_1,g2_1), e2=veorq_u8(g1_2,g2_2);
		const uint8x16_t mask=vandq_u8(vornq_u8(g1_1,g1_2), vornq_u8(g2_1,g2_2));
		const uint8x16_t ibs0=vandq_u8(vandq_u8(e1,e2), mask);
		const uint8x16_t het=vandq_u8(veorq_u8(e1,e2), mask);  // (e1^e2)&mask
		const uint8x16_t aa1=vandq_u8(vbicq_u8(g1_1,g1_2), mask);
		const uint8x16_t aa2=vandq_u8(vbicq_u8(g2_1,g2_2), mask);
		Ai=vpadalq_u8(Ai, vcntq_u8(ibs0));
		Am=vpadalq_u8(Am, vcntq_u8(mask));
		Ah=vpadalq_u8(Ah, vcntq_u8(het));
		A1=vpadalq_u8(A1, vcntq_u8(aa1));
		A2=vpadalq_u8(A2, vcntq_u8(aa2));
		if (++iter == 4000)
		{
			c_ibs0+=vaddlvq_u16(Ai); c_mask+=vaddlvq_u16(Am); c_het+=vaddlvq_u16(Ah);
			c_aa1+=vaddlvq_u16(A1); c_aa2+=vaddlvq_u16(A2);
			Ai=vdupq_n_u16(0); Am=vdupq_n_u16(0); Ah=vdupq_n_u16(0);
			A1=vdupq_n_u16(0); A2=vdupq_n_u16(0); iter=0;
		}
	}
	c_ibs0+=vaddlvq_u16(Ai); c_mask+=vaddlvq_u16(Am); c_het+=vaddlvq_u16(Ah);
	c_aa1+=vaddlvq_u16(A1); c_aa2+=vaddlvq_u16(A2);
	king_tail(p1a, p1b, p2a, p2b, m, c_ibs0, c_mask, c_het, c_aa1, c_aa2);
	out.ibs0  += (uint32_t)c_ibs0;
	out.nloci += (uint32_t)c_mask;
	out.sumsq += (uint32_t)(c_het + 4*c_ibs0);
	out.n1_aa += (uint32_t)c_aa1;
	out.n2_aa += (uint32_t)c_aa2;
}

// ---- PCA/GRM dot product: 2x f64 FMA lanes, 2 accumulators for ILP ----
double dot_f64_neon(const double *a, const double *b, size_t n)
{
	float64x2_t s0=vdupq_n_f64(0), s1=vdupq_n_f64(0);
	size_t k=0;
	for (; k+4 <= n; k += 4)
	{
		s0 = vfmaq_f64(s0, vld1q_f64(a+k),   vld1q_f64(b+k));
		s1 = vfmaq_f64(s1, vld1q_f64(a+k+2), vld1q_f64(b+k+2));
	}
	double s = vaddvq_f64(vaddq_f64(s0, s1));
	for (; k < n; k++) s += a[k]*b[k];
	return s;
}

}  // namespace SNPvec

#endif  /* COREARRAY_SIMD_NEON */
