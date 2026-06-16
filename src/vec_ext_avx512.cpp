// ===========================================================
//
// vec_ext_avx512.cpp: x86 AVX-512 (F+BW) IBS kernel (SNPRelate v2)
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

#if defined(SNP_VEC_X86) && defined(SNP_VEC_X86_TARGET)

// GCC: enable AVX-512F/BW for the whole TU before <immintrin.h>.
#if defined(__GNUC__) && !defined(__clang__)
#   pragma GCC push_options
#   pragma GCC target("avx512f,avx512bw")
#endif
#include <immintrin.h>

#define KTGT  __attribute__((target("avx512f,avx512bw")))

namespace SNPvec
{

KTGT static inline __m512i bytecnt(__m512i v)
{
	const __m512i m1=_mm512_set1_epi8((char)0x55), m2=_mm512_set1_epi8((char)0x33),
		m4=_mm512_set1_epi8((char)0x0F);
	v = _mm512_sub_epi8(v, _mm512_and_si512(_mm512_srli_epi16(v,1), m1));
	v = _mm512_add_epi8(_mm512_and_si512(v,m2),
		_mm512_and_si512(_mm512_srli_epi16(v,2), m2));
	return _mm512_and_si512(_mm512_add_epi8(v,_mm512_srli_epi16(v,4)), m4);
}
// horizontal sum of 8 x u64 lanes (manual; avoids _mm512_reduce_* portability)
KTGT static inline uint64_t hsum64(__m512i v)
{
	__m256i s = _mm256_add_epi64(_mm512_castsi512_si256(v),
		_mm512_extracti64x4_epi64(v, 1));
	__m128i t = _mm_add_epi64(_mm256_castsi256_si128(s),
		_mm256_extracti128_si256(s, 1));
	return (uint64_t)_mm_cvtsi128_si64(_mm_add_epi64(t, _mm_srli_si128(t, 8)));
}

KTGT
void ibs_count_avx512(const uint8_t *gi, const uint8_t *gj, size_t npack,
	uint32_t &ibs0, uint32_t &ibs1, uint32_t &ibs2)
{
	const __m512i ONES = _mm512_set1_epi32(-1), ZERO = _mm512_setzero_si512();
	__m512i a0=ZERO, a2=ZERO, am=ZERO;
	const uint8_t *p1=gi, *p2=gj;
	size_t m = npack;
	for (; m >= 64; m -= 64)
	{
		const __m512i g1_1=_mm512_loadu_si512((const void*)p1);
		const __m512i g1_2=_mm512_loadu_si512((const void*)(p1+npack));
		const __m512i g2_1=_mm512_loadu_si512((const void*)p2);
		const __m512i g2_2=_mm512_loadu_si512((const void*)(p2+npack));
		p1 += 64; p2 += 64;
		const __m512i e1=_mm512_xor_si512(g1_1,g2_1), e2=_mm512_xor_si512(g1_2,g2_2);
		const __m512i mask=_mm512_and_si512(
			_mm512_or_si512(g1_1,_mm512_andnot_si512(g1_2,ONES)),
			_mm512_or_si512(g2_1,_mm512_andnot_si512(g2_2,ONES)));
		const __m512i v0=_mm512_and_si512(_mm512_and_si512(e1,e2),mask);
		const __m512i v2=_mm512_andnot_si512(e1,_mm512_andnot_si512(e2,mask));
		a0=_mm512_add_epi64(a0,_mm512_sad_epu8(bytecnt(v0),ZERO));
		a2=_mm512_add_epi64(a2,_mm512_sad_epu8(bytecnt(v2),ZERO));
		am=_mm512_add_epi64(am,_mm512_sad_epu8(bytecnt(mask),ZERO));
	}
	uint64_t c0=hsum64(a0), c2=hsum64(a2), cm=hsum64(am);
	ibs_tail(p1, p1+npack, p2, p2+npack, m, c0, c2, cm);
	ibs0 += (uint32_t)c0; ibs2 += (uint32_t)c2; ibs1 += (uint32_t)(cm-c0-c2);
}

}  // namespace SNPvec

#if defined(__GNUC__) && !defined(__clang__)
#   pragma GCC pop_options
#endif

#endif  /* SNP_VEC_X86 && SNP_VEC_X86_TARGET */
