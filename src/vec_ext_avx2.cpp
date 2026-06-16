// ===========================================================
//
// vec_ext_avx2.cpp: x86 AVX2 IBS kernel (SNPRelate v2)
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

// GCC: enable AVX2 (intrinsics + codegen) for the whole TU before <immintrin.h>.
#if defined(__GNUC__) && !defined(__clang__)
#   pragma GCC push_options
#   pragma GCC target("avx2")
#endif
#include <immintrin.h>

#define KTGT  __attribute__((target("avx2")))

namespace SNPvec
{

KTGT static inline __m256i bytecnt(__m256i v)
{
	const __m256i m1=_mm256_set1_epi8((char)0x55), m2=_mm256_set1_epi8((char)0x33),
		m4=_mm256_set1_epi8((char)0x0F);
	v = _mm256_sub_epi8(v, _mm256_and_si256(_mm256_srli_epi16(v,1), m1));
	v = _mm256_add_epi8(_mm256_and_si256(v,m2),
		_mm256_and_si256(_mm256_srli_epi16(v,2), m2));
	return _mm256_and_si256(_mm256_add_epi8(v,_mm256_srli_epi16(v,4)), m4);
}
KTGT static inline uint64_t hsum64(__m256i v)
{
	const __m128i s = _mm_add_epi64(_mm256_castsi256_si128(v),
		_mm256_extracti128_si256(v, 1));
	return (uint64_t)_mm_cvtsi128_si64(_mm_add_epi64(s, _mm_srli_si128(s, 8)));
}

KTGT
void ibs_count_avx2(const uint8_t *gi, const uint8_t *gj, size_t npack,
	uint32_t &ibs0, uint32_t &ibs1, uint32_t &ibs2)
{
	const __m256i ONES = _mm256_set1_epi32(-1), ZERO = _mm256_setzero_si256();
	__m256i a0=ZERO, a2=ZERO, am=ZERO;
	const uint8_t *p1=gi, *p2=gj;
	size_t m = npack;
	for (; m >= 32; m -= 32)
	{
		const __m256i g1_1=_mm256_loadu_si256((const __m256i*)p1);
		const __m256i g1_2=_mm256_loadu_si256((const __m256i*)(p1+npack));
		const __m256i g2_1=_mm256_loadu_si256((const __m256i*)p2);
		const __m256i g2_2=_mm256_loadu_si256((const __m256i*)(p2+npack));
		p1 += 32; p2 += 32;
		const __m256i e1=_mm256_xor_si256(g1_1,g2_1), e2=_mm256_xor_si256(g1_2,g2_2);
		const __m256i mask=_mm256_and_si256(
			_mm256_or_si256(g1_1,_mm256_andnot_si256(g1_2,ONES)),
			_mm256_or_si256(g2_1,_mm256_andnot_si256(g2_2,ONES)));
		const __m256i v0=_mm256_and_si256(_mm256_and_si256(e1,e2),mask);
		const __m256i v2=_mm256_andnot_si256(e1,_mm256_andnot_si256(e2,mask));
		a0=_mm256_add_epi64(a0,_mm256_sad_epu8(bytecnt(v0),ZERO));
		a2=_mm256_add_epi64(a2,_mm256_sad_epu8(bytecnt(v2),ZERO));
		am=_mm256_add_epi64(am,_mm256_sad_epu8(bytecnt(mask),ZERO));
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
