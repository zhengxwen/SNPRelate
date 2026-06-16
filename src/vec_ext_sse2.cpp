// ===========================================================
//
// vec_ext_sse2.cpp: x86 SSE2 IBS kernel (SNPRelate v2)
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

#ifdef SNP_VEC_X86

// GCC: set the target BEFORE <immintrin.h> so the SSE2 intrinsics are declared
// and code is generated for this whole TU (a function-only target attribute is
// not enough on GCC, whose <immintrin.h> guards intrinsics by -m macros).
// Clang declares all intrinsics unconditionally and keys codegen off the
// per-function attribute below, so it does not need the pragma.
#if defined(__GNUC__) && !defined(__clang__)
#   pragma GCC push_options
#   pragma GCC target("sse2")
#endif
#include <immintrin.h>

#ifdef SNP_VEC_X86_TARGET
#   define KTGT  __attribute__((target("sse2")))
#else
#   define KTGT
#endif

namespace SNPvec
{

// byte-wise popcount (each byte -> 0..8), SWAR; reduced later via psadbw
KTGT static inline __m128i bytecnt(__m128i v)
{
	const __m128i m1=_mm_set1_epi8((char)0x55), m2=_mm_set1_epi8((char)0x33),
		m4=_mm_set1_epi8((char)0x0F);
	v = _mm_sub_epi8(v, _mm_and_si128(_mm_srli_epi16(v, 1), m1));
	v = _mm_add_epi8(_mm_and_si128(v, m2), _mm_and_si128(_mm_srli_epi16(v, 2), m2));
	return _mm_and_si128(_mm_add_epi8(v, _mm_srli_epi16(v, 4)), m4);
}
KTGT static inline uint64_t hsum64(__m128i v)
{
	return (uint64_t)_mm_cvtsi128_si64(_mm_add_epi64(v, _mm_srli_si128(v, 8)));
}

KTGT
void ibs_count_sse2(const uint8_t *gi, const uint8_t *gj, size_t npack,
	uint32_t &ibs0, uint32_t &ibs1, uint32_t &ibs2)
{
	const __m128i ONES = _mm_set1_epi32(-1), ZERO = _mm_setzero_si128();
	__m128i a0=ZERO, a2=ZERO, am=ZERO;
	const uint8_t *p1=gi, *p2=gj;
	size_t m = npack;
	for (; m >= 16; m -= 16)
	{
		const __m128i g1_1=_mm_loadu_si128((const __m128i*)p1);
		const __m128i g1_2=_mm_loadu_si128((const __m128i*)(p1+npack));
		const __m128i g2_1=_mm_loadu_si128((const __m128i*)p2);
		const __m128i g2_2=_mm_loadu_si128((const __m128i*)(p2+npack));
		p1 += 16; p2 += 16;
		const __m128i e1=_mm_xor_si128(g1_1,g2_1), e2=_mm_xor_si128(g1_2,g2_2);
		const __m128i mask=_mm_and_si128(
			_mm_or_si128(g1_1,_mm_andnot_si128(g1_2,ONES)),
			_mm_or_si128(g2_1,_mm_andnot_si128(g2_2,ONES)));
		const __m128i v0=_mm_and_si128(_mm_and_si128(e1,e2),mask);
		const __m128i v2=_mm_andnot_si128(e1,_mm_andnot_si128(e2,mask));
		a0=_mm_add_epi64(a0,_mm_sad_epu8(bytecnt(v0),ZERO));
		a2=_mm_add_epi64(a2,_mm_sad_epu8(bytecnt(v2),ZERO));
		am=_mm_add_epi64(am,_mm_sad_epu8(bytecnt(mask),ZERO));
	}
	uint64_t c0=hsum64(a0), c2=hsum64(a2), cm=hsum64(am);
	ibs_tail(p1, p1+npack, p2, p2+npack, m, c0, c2, cm);
	ibs0 += (uint32_t)c0; ibs2 += (uint32_t)c2; ibs1 += (uint32_t)(cm-c0-c2);
}

}  // namespace SNPvec

#if defined(__GNUC__) && !defined(__clang__)
#   pragma GCC pop_options
#endif

#endif  /* SNP_VEC_X86 */
