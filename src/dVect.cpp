// ===========================================================
//
// dVect.cpp: Classess and functions for vectorization
//
// Copyright (C) 2007-2017    Xiuwen Zheng
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
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with SNPRelate.
// If not, see <http://www.gnu.org/licenses/>.


#include "dVect.h"


namespace Vectorization
{

// count genotype sum and number of calls, not requiring 16-aligned p
COREARRAY_DLL_DEFAULT C_UInt8* vec_u8_geno_count(C_UInt8 *p,
	size_t n, C_Int32 &out_sum, C_Int32 &out_num)
{
	C_Int32 sum=0, num=0;

#if defined(COREARRAY_SIMD_AVX2)

	const __m256i three = _mm256_set1_epi8(3);
	const __m256i zero = _mm256_setzero_si256();
	__m256i sum32 = zero, num32 = zero;
	size_t limit_by_U8 = 0;

	for (; n >= 32; )
	{
		__m256i v = _mm256_loadu_si256((__m256i const*)p);
		p += 32;
		__m256i m = _mm256_cmpgt_epi8(three, _mm256_min_epu8(v, three));
		sum32 = _mm256_add_epi8(sum32, _mm256_and_si256(v, m));
		num32 = _mm256_sub_epi8(num32, m);
		n -= 32;
		limit_by_U8 ++;
		if ((limit_by_U8 >= 127) || (n < 32))
		{
			// add to sum
			sum32 = _mm256_sad_epu8(sum32, zero);
			sum32 = _mm256_add_epi32(sum32,
				_mm256_permute4x64_epi64(sum32, _MM_SHUFFLE(1,0,3,2)));
			sum32 = _mm256_add_epi32(sum32,
				_mm256_permute4x64_epi64(sum32, _MM_SHUFFLE(0,0,0,1)));
			sum += _mm_cvtsi128_si32(_mm256_castsi256_si128(sum32));
			// add to num
			num32 = _mm256_sad_epu8(num32, zero);
			num32 = _mm256_add_epi32(num32,
				_mm256_permute4x64_epi64(num32, _MM_SHUFFLE(1,0,3,2)));
			num32 = _mm256_add_epi32(num32,
				_mm256_permute4x64_epi64(num32, _MM_SHUFFLE(0,0,0,1)));
			num += _mm_cvtsi128_si32(_mm256_castsi256_si128(num32));
			// reset
			sum32 = num32 = zero;
			limit_by_U8 = 0;
		}
	}

#elif defined(COREARRAY_SIMD_SSE2)

	// header, 16-byte aligned
	size_t h = (16 - ((size_t)p & 0x0F)) & 0x0F;
	for (; (n > 0) && (h > 0); n--, h--, p++)
		if (*p <= 2) { sum += *p; num++; }

	const __m128i three = _mm_set1_epi8(3);
	const __m128i zero = _mm_setzero_si128();
	__m128i sum16=zero, num16=zero;
	size_t limit_by_U8 = 0;

	for (; n >= 16; )
	{
		__m128i v = _mm_load_si128((__m128i const*)p);
		p += 16;
		__m128i m = _mm_cmpgt_epi8(three, _mm_min_epu8(v, three));
		sum16 = _mm_add_epi8(sum16, v & m);
		num16 = _mm_sub_epi8(num16, m);
		n -= 16;
		limit_by_U8 ++;
		if ((limit_by_U8 >= 127) || (n < 16))
		{
			// add to sum
			sum16 = _mm_sad_epu8(sum16, zero);
			sum += _mm_cvtsi128_si32(sum16);
			sum += _mm_cvtsi128_si32(_mm_shuffle_epi32(sum16, 2));
			// add to num
			num16 = _mm_sad_epu8(num16, zero);
			num += _mm_cvtsi128_si32(num16);
			num += _mm_cvtsi128_si32(_mm_shuffle_epi32(num16, 2));
			// reset
			sum16 = num16 = zero;
			limit_by_U8 = 0;
		}
	}

#endif

	for (; n > 0; n--, p++)
		if (*p <= 2) { sum += *p; num++; }
	out_sum = sum;
	out_num = num;
	return p;
}


// any (*p > 2) is set to be 3
COREARRAY_DLL_DEFAULT void vec_u8_geno_valid(C_UInt8 *p, size_t n)
{
#if defined(COREARRAY_SIMD_SSE2)

	// header 1, 16-byte aligned
	size_t h = (16 - ((size_t)p & 0x0F)) & 0x0F;
	for (; (n > 0) && (h > 0); n--, h--, p++)
		if (*p > 3) *p = 3;

	const __m128i zero  = _mm_setzero_si128();
	const __m128i three = _mm_set1_epi8(3);
	for (; n >= 16; n-=16, p+=16)
	{
		__m128i v = _mm_load_si128((__m128i*)p);
		__m128i mask = _mm_or_si128(_mm_cmplt_epi8(v, zero),
			_mm_cmplt_epi8(three, v));
		if (_mm_movemask_epi8(mask) > 0)
			_mm_store_si128((__m128i*)p, _mm_min_epu8(v, three));
	}

#endif

	for (; n > 0; n--, p++) if (*p > 3) *p = 3;
}


// add *p by v and applied to all n
COREARRAY_DLL_DEFAULT void vec_i32_add(C_Int32 *p, size_t n, C_Int32 v)
{
	for (; n > 0; n--) *p++ += v;
}


// add *p by v and applied to all n
COREARRAY_DLL_DEFAULT void vec_f64_add(double *p, size_t n, double v)
{
#if defined(COREARRAY_SIMD_AVX)

	__m256d v4 = _mm256_set1_pd(v);
	switch ((size_t)p & 0x1F)
	{
	case 0x08:
		if (n > 0) { (*p++) += v; n--; }
	case 0x10:
		if (n > 0) { (*p++) += v; n--; }
	case 0x18:
		if (n > 0) { (*p++) += v; n--; }
	case 0x00:
		for (; n >= 4; n-=4)
		{
			_mm256_store_pd(p, _mm256_add_pd(_mm256_load_pd(p), v4));
			p += 4;
		}
		if (n >= 2)
		{
			_mm_store_pd(p, _mm_add_pd(_mm_load_pd(p), _mm256_castpd256_pd128(v4)));
			p += 2; n -= 2;
		}
		break;
	default:
		for (; n >= 4; n-=4)
		{
			_mm256_storeu_pd(p, _mm256_add_pd(_mm256_loadu_pd(p), v4));
			p += 4;
		}
		if (n >= 2)
		{
			_mm_storeu_pd(p, _mm_add_pd(_mm_loadu_pd(p), _mm256_castpd256_pd128(v4)));
			p += 2; n -= 2;
		}
	}

#elif defined(COREARRAY_SIMD_SSE2)

	__m128d v2 = _mm_set1_pd(v);
	switch ((size_t)p & 0x0F)
	{
	case 0x08:
		if (n > 0) { (*p++) += v; n--; }
	case 0x00:
		for (; n >= 2; n-=2)
		{
			_mm_store_pd(p, _mm_add_pd(_mm_load_pd(p), v2));
			p += 2;
		}
		break;
	default:
		for (; n >= 2; n-=2)
		{
			_mm_storeu_pd(p, _mm_add_pd(_mm_loadu_pd(p), v2));
			p += 2;
		}
	}

#endif

	for (; n > 0; n--) (*p++) += v;
}


// add *p by *s and applied to all n
COREARRAY_DLL_DEFAULT void vec_f64_add(double *p, const double *s, size_t n)
{
#if defined(COREARRAY_SIMD_AVX)

	switch ((size_t)p & 0x1F)
	{
	case 0x08:
		if (n > 0) { (*p++) += (*s++); n--; }
	case 0x10:
		if (n > 0) { (*p++) += (*s++); n--; }
	case 0x18:
		if (n > 0) { (*p++) += (*s++); n--; }
	case 0x00:
		for (; n >= 4; n-=4)
		{
			_mm256_store_pd(p, _mm256_add_pd(_mm256_load_pd(p), _mm256_loadu_pd(s)));
			p += 4; s += 4;
		}
		if (n >= 2)
		{
			_mm_store_pd(p, _mm_add_pd(_mm_load_pd(p), _mm_loadu_pd(s)));
			p += 2; s += 2; n -= 2;
		}
		break;
	default:
		for (; n >= 4; n-=4)
		{
			_mm256_storeu_pd(p, _mm256_add_pd(_mm256_loadu_pd(p), _mm256_loadu_pd(s)));
			p += 4; s += 4;
		}
		if (n >= 2)
		{
			_mm_storeu_pd(p, _mm_add_pd(_mm_loadu_pd(p), _mm_loadu_pd(s)));
			p += 2; s += 2; n -= 2;
		}
	}

#elif defined(COREARRAY_SIMD_SSE2)

	switch ((size_t)p & 0x0F)
	{
	case 0x08:
		if (n > 0) { (*p++) += (*s++); n--; }
	case 0x00:
		for (; n >= 2; n-=2)
		{
			_mm_store_pd(p, _mm_add_pd(_mm_load_pd(p), _mm_loadu_pd(s)));
			p += 2; s += 2;
		}
		break;
	default:
		for (; n >= 2; n-=2)
		{
			_mm_storeu_pd(p, _mm_add_pd(_mm_loadu_pd(p), _mm_loadu_pd(s)));
			p += 2; s += 2;
		}
	}

#endif

	for (; n > 0; n--) (*p++) += (*s++);
}


// subtract *p by v and applied to all n
COREARRAY_DLL_DEFAULT void vec_f64_sub(double *p, size_t n, double v)
{
	vec_f64_add(p, n, -v);
}


// *p = v - *p and applied to all n
COREARRAY_DLL_DEFAULT void vec_f64_sub2(double *p, size_t n, double v)
{
#if defined(COREARRAY_SIMD_AVX)

	__m256d v4 = _mm256_set1_pd(v);
	switch ((size_t)p & 0x1F)
	{
	case 0x08:
		if (n > 0) { *p = v - *p; p++; n--; }
	case 0x10:
		if (n > 0) { *p = v - *p; p++; n--; }
	case 0x18:
		if (n > 0) { *p = v - *p; p++; n--; }
	case 0x00:
		for (; n >= 4; n-=4)
		{
			_mm256_store_pd(p, _mm256_sub_pd(v4, _mm256_load_pd(p)));
			p += 4;
		}
		if (n >= 2)
		{
			_mm_store_pd(p, _mm_sub_pd(_mm256_castpd256_pd128(v4), _mm_load_pd(p)));
			p += 2; n -= 2;
		}
		break;
	default:
		for (; n >= 4; n-=4)
		{
			_mm256_storeu_pd(p, _mm256_sub_pd(v4, _mm256_loadu_pd(p)));
			p += 4;
		}
		if (n >= 2)
		{
			_mm_storeu_pd(p, _mm_sub_pd(_mm256_castpd256_pd128(v4), _mm_loadu_pd(p)));
			p += 2; n -= 2;
		}
	}

#elif defined(COREARRAY_SIMD_SSE2)

	__m128d v2 = _mm_set1_pd(v);
	switch ((size_t)p & 0x0F)
	{
	case 0x08:
		if (n > 0) { *p = v - *p; p++; n--; }
	case 0x00:
		for (; n >= 2; n-=2)
		{
			_mm_store_pd(p, _mm_sub_pd(v2, _mm_load_pd(p)));
			p += 2;
		}
		break;
	default:
		for (; n >= 2; n-=2)
		{
			_mm_storeu_pd(p, _mm_sub_pd(v2, _mm_loadu_pd(p)));
			p += 2;
		}
	}

#endif

	for (; n > 0; n--, p++) *p = v - *p;
}


// multiply *p by v and applied to all n
COREARRAY_DLL_DEFAULT void vec_f64_mul(double *p, size_t n, double v)
{
#if defined(COREARRAY_SIMD_AVX)

	const __m256d v4 = _mm256_set1_pd(v);

	switch ((size_t)p & 0x1F)
	{
	case 0x08:
		if (n > 0) { (*p++) *= v; n--; }
	case 0x10:
		if (n > 0) { (*p++) *= v; n--; }
	case 0x18:
		if (n > 0) { (*p++) *= v; n--; }
	case 0x00:
		for (; n >= 4; n-=4)
		{
			_mm256_store_pd(p, _mm256_mul_pd(_mm256_load_pd(p), v4));
			p += 4;
		}
		if (n >= 2)
		{
			_mm_store_pd(p, _mm_mul_pd(_mm_load_pd(p), _mm256_castpd256_pd128(v4)));
			p += 2; n -= 2;
		}
		break;
	default:
		for (; n >= 4; n-=4)
		{
			_mm256_storeu_pd(p, _mm256_mul_pd(_mm256_loadu_pd(p), v4));
			p += 4;
		}
		if (n >= 2)
		{
			_mm_storeu_pd(p, _mm_mul_pd(_mm_loadu_pd(p), _mm256_castpd256_pd128(v4)));
			p += 2; n -= 2;
		}
	}

#elif defined(COREARRAY_SIMD_SSE2)

	const __m128d v2 = _mm_set1_pd(v);

	switch ((size_t)p & 0x0F)
	{
	case 0x08:
		if (n > 0) { (*p++) *= v; n--; }
	case 0x00:
		for (; n >= 2; n-=2, p+=2)
			_mm_store_pd(p, _mm_mul_pd(_mm_load_pd(p), v2));
		break;
	default:
		for (; n >= 2; n-=2, p+=2)
			_mm_storeu_pd(p, _mm_mul_pd(_mm_loadu_pd(p), v2));
	}

#endif

	for (; n > 0; n--) (*p++) *= v;
}


// *p += (*s) * v
COREARRAY_DLL_DEFAULT double *vec_f64_addmul(double *p, const double *s,
	size_t n, double v)
{
#if defined(COREARRAY_SIMD_SSE2)

	const __m128d v2 = _mm_set1_pd(v);

	switch ((size_t)p & 0x0F)
	{
	case 0x08:
		if (n > 0) { (*p++) += (*s++) * v; n--; }
	case 0x00:
		for (; n >= 2; n -= 2)
		{
			_mm_store_pd(p, _mm_add_pd(_mm_load_pd(p),
				_mm_mul_pd(_mm_loadu_pd(s), v2)));
			p += 2; s += 2;
		}
		break;
	default:
		for (; n >= 2; n-=2)
		{
			_mm_storeu_pd(p, _mm_add_pd(_mm_loadu_pd(p),
				_mm_mul_pd(_mm_loadu_pd(s), v2)));
			p += 2; s += 2;
		}
	}

#endif

	for (; n > 0; n--) (*p++) += (*s++) * v;
	return p;
}

}
