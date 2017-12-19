// ===========================================================
//
// dVect.h: Classess and functions for vectorization
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
// You should have received a copy of the GNU General Public
// License along with SNPRelate.
// If not, see <http://www.gnu.org/licenses/>.

#ifndef _HEADER_VECTORIZATION_
#define _HEADER_VECTORIZATION_

#include <dType.h>

#ifdef COREARRAY_SIMD_SSE
#   include <xmmintrin.h>
#endif
#ifdef COREARRAY_SIMD_SSE2
#   include <emmintrin.h>
#endif
#ifdef COREARRAY_SIMD_SSE4_1
#   include <smmintrin.h>
#endif
#if defined(COREARRAY_SIMD_SSE4_2) || defined(__POPCNT__)
#   define VECT_HARDWARE_POPCNT
#   include <nmmintrin.h>  // SSE4_2, for POPCNT
#endif
#ifdef COREARRAY_SIMD_AVX
#   include <immintrin.h>
#endif


namespace Vectorization
{

#if defined(COREARRAY_SIMD_AVX)
	#define VEC_SIMD_ALIGN_BYTE    32u
#elif defined(COREARRAY_SIMD_SSE)
	#define VEC_SIMD_ALIGN_BYTE    16u
#else
	#define VEC_SIMD_ALIGN_BYTE    1u
#endif

	template<typename TYPE, size_t vAlign=VEC_SIMD_ALIGN_BYTE>
		struct COREARRAY_DLL_DEFAULT VEC_AUTO_PTR
	{
	public:
		static const size_t Align = vAlign;

		VEC_AUTO_PTR()
			{ alloc = NULL; base = NULL; vlen = 0; }
		VEC_AUTO_PTR(size_t n)
			{ alloc = NULL; base = NULL; vlen = 0; Reset(n); }
		~VEC_AUTO_PTR()
			{ Clear(); }

		void Reset(size_t n)
		{
			if (n != vlen)
			{
				if (alloc) delete []alloc;
				if (n > 0)
				{
					size_t m = n*sizeof(TYPE) + Align - 1;
					alloc = new char[m];
					size_t r = ((size_t)alloc) % Align;
					base = (TYPE*)(r ? (alloc+Align-r) : alloc);
					vlen = n;
				} else {
					alloc = NULL; base = NULL;
					vlen = 0;
				}
			}
		}

		void Clear()
		{
			if (alloc)
				{ delete []alloc; alloc = NULL; }
			base = NULL;
			vlen = 0;
		}

		COREARRAY_INLINE TYPE *Get() { return base; };
		COREARRAY_INLINE size_t Length() { return vlen; };

		COREARRAY_INLINE TYPE &operator[] (size_t i) { return base[i]; }

	private:
		char *alloc;
		TYPE *base;
		size_t vlen;
	};



	// ===========================================================

#ifdef COREARRAY_SIMD_SSE2
#   ifdef COREARRAY_SIMD_SSE4_1
#       define MM_MUL_LO_EPI32(dst, a, b)  dst = _mm_mullo_epi32(a, b)
#   else
#       define MM_MUL_LO_EPI32(dst, a, b)  \
        {  \
            __m128i i1 = _mm_mul_epu32(a, b);  \
            __m128i i2 = _mm_mul_epu32(_mm_srli_si128(a, 4), _mm_srli_si128(b, 4));  \
            dst = _mm_unpacklo_epi32(_mm_shuffle_epi32(i1, _MM_SHUFFLE (0,0,2,0)),  \
                _mm_shuffle_epi32(i2, _MM_SHUFFLE (0,0,2,0)));  \
        }
#   endif
#endif


	// ===========================================================

	inline static int POPCNT_U32(C_UInt32 x)
	{
		x = x - ((x >> 1) & 0x55555555);
		x = (x & 0x33333333) + ((x >> 2) & 0x33333333);
		return (((x + (x >> 4)) & 0xF0F0F0F) * 0x1010101) >> 24;
	}

	inline static int POPCNT_U64(C_UInt64 x)
	{
		x -= ((x >> 1) & 0x5555555555555555LLU);
		x = (x & 0x3333333333333333LLU) + ((x >> 2) & 0x3333333333333333LLU);
		x = (x + (x >> 4)) & 0x0F0F0F0F0F0F0F0FLLU;
		return (x * 0x0101010101010101LLU) >> 56;
	}

#ifdef VECT_HARDWARE_POPCNT

	inline static int POPCNT_M128(__m128i x)
	{
	#ifdef __LP64__
		return _mm_popcnt_u64(_mm_cvtsi128_si64(x)) + 
			_mm_popcnt_u64(_mm_cvtsi128_si64(_mm_srli_si128(x, 8)));
	#else
		return _mm_popcnt_u32(_mm_cvtsi128_si32(x)) + 
			_mm_popcnt_u32(_mm_cvtsi128_si32(_mm_srli_si128(x, 4))) +
			_mm_popcnt_u32(_mm_cvtsi128_si32(_mm_srli_si128(x, 8))) +
			_mm_popcnt_u32(_mm_cvtsi128_si32(_mm_srli_si128(x, 12)));
	#endif
	}

	#ifdef COREARRAY_SIMD_AVX2
	inline static int POPCNT_M256(__m256i x)
	{
	#ifdef __LP64__
		__m128i v1 = _mm256_extractf128_si256(x, 0);
		__m128i v2 = _mm256_extractf128_si256(x, 1);
		return _mm_popcnt_u64(_mm_cvtsi128_si64(v1)) + 
			_mm_popcnt_u64(_mm_cvtsi128_si64(_mm_srli_si128(v1, 8))) +
			_mm_popcnt_u64(_mm_cvtsi128_si64(v2)) + 
			_mm_popcnt_u64(_mm_cvtsi128_si64(_mm_srli_si128(v2, 8)));
	#else
		__m128i v1 = _mm256_extractf128_si256(x, 0);
		__m128i v2 = _mm256_extractf128_si256(x, 1);
		return _mm_popcnt_u32(_mm_cvtsi128_si32(v1)) + 
			_mm_popcnt_u32(_mm_cvtsi128_si32(_mm_srli_si128(v1, 4))) +
			_mm_popcnt_u32(_mm_cvtsi128_si32(_mm_srli_si128(v1, 8))) +
			_mm_popcnt_u32(_mm_cvtsi128_si32(_mm_srli_si128(v1, 12))) +
			_mm_popcnt_u32(_mm_cvtsi128_si32(v2)) + 
			_mm_popcnt_u32(_mm_cvtsi128_si32(_mm_srli_si128(v2, 4))) +
			_mm_popcnt_u32(_mm_cvtsi128_si32(_mm_srli_si128(v2, 8))) +
			_mm_popcnt_u32(_mm_cvtsi128_si32(_mm_srli_si128(v2, 12)));
	#endif
	}
	#endif

#endif


#ifdef COREARRAY_SIMD_SSE2

	#define POPCNT_SSE2_HEAD    \
		const __m128i pop_sse_5 = _mm_set1_epi32(0x55555555);  \
		const __m128i pop_sse_3 = _mm_set1_epi32(0x33333333);  \
		const __m128i pop_sse_f = _mm_set1_epi32(0x0F0F0F0F);  \
		const __m128i pop_sse_1 = _mm_set1_epi32(0x01010101);

	#define POPCNT_SSE2_RUN(x)    \
		x = _mm_sub_epi32(x, _mm_and_si128(_mm_srli_epi32(x, 1), pop_sse_5));  \
		x = _mm_add_epi32(_mm_and_si128(x, pop_sse_3),  \
			_mm_and_si128(_mm_srli_epi32(x, 2), pop_sse_3));  \
		x = _mm_and_si128(_mm_add_epi32(x, _mm_srli_epi32(x, 4)), pop_sse_f);  \
		MM_MUL_LO_EPI32(x, x, pop_sse_1);  \
		x = _mm_srli_epi32(x, 24);

#endif

#ifdef COREARRAY_SIMD_AVX2

	#define POPCNT_AVX2_HEAD    \
		const __m256i pop_sse_5 = _mm256_set1_epi32(0x55555555);  \
		const __m256i pop_sse_3 = _mm256_set1_epi32(0x33333333);  \
		const __m256i pop_sse_f = _mm256_set1_epi32(0x0F0F0F0F);  \
		const __m256i pop_sse_1 = _mm256_set1_epi32(0x01010101);

	#define POPCNT_AVX2_RUN(x)    \
		x = _mm256_sub_epi32(x, _mm256_and_si256(_mm256_srli_epi32(x, 1), pop_sse_5));  \
		x = _mm256_add_epi32(_mm256_and_si256(x, pop_sse_3),  \
			_mm256_and_si256(_mm256_srli_epi32(x, 2), pop_sse_3));  \
		x = _mm256_and_si256(_mm256_add_epi32(x, _mm256_srli_epi32(x, 4)), pop_sse_f);  \
		x = _mm256_mullo_epi32(x, pop_sse_1);  \
		x = _mm256_srli_epi32(x, 24);

#endif


	// ===========================================================
	// Sum all elements in a SIMD register

#ifdef COREARRAY_SIMD_SSE2
	inline static double vec_sum_f64(__m128d s)
	{
		return _mm_cvtsd_f64(_mm_add_pd(s, _mm_shuffle_pd(s, s, 1)));
	}
	inline static int vec_sum_i32(__m128i s)
	{
		s = _mm_add_epi32(s, _mm_shuffle_epi32(s, _MM_SHUFFLE(1,0,3,2)));
		s = _mm_add_epi32(s, _mm_shuffle_epi32(s, _MM_SHUFFLE(0,0,0,1)));
		return _mm_cvtsi128_si32(s);
	}
#endif

#ifdef COREARRAY_SIMD_AVX
	inline static double vec_avx_sum_f64(__m256d s)
	{
		s = _mm256_add_pd(_mm256_permute_pd(s, 5), s);
		s = _mm256_add_pd(s, _mm256_permute2f128_pd(s, s, 1));
		return _mm_cvtsd_f64(_mm256_castpd256_pd128(s));
	}
#endif

#ifdef COREARRAY_SIMD_AVX2
	inline static int vec_avx_sum_i32(__m256i s)
	{
		s = _mm256_hadd_epi32(s, s);
		s = _mm256_add_epi32(s, _mm256_permute4x64_epi64(s, _MM_SHUFFLE(1,0,3,2)));
		__m128i a = _mm256_castsi256_si128(s);
		a = _mm_add_epi32(a, _mm_shuffle_epi32(a, _MM_SHUFFLE(0,0,0,1)));
		return _mm_cvtsi128_si32(a);
	}
#endif


	// ===========================================================
	// NOT bitwise operator for compiler compatibility

#ifdef COREARRAY_SIMD_SSE2
	#define SIMD128_NOT_HEAD  \
		const __m128i ONE128 = \
			_mm_cmpeq_epi32(_mm_setzero_si128(), _mm_setzero_si128());
	#define SIMD128_NOT(x)    _mm_andnot_si128(x, ONE128)
#endif

#ifdef COREARRAY_SIMD_AVX2
	#define SIMD256_NOT_HEAD  \
		const __m256i ONE256 = \
			_mm256_cmpeq_epi32(_mm256_setzero_si256(), _mm256_setzero_si256());
	#define SIMD256_NOT(x)    _mm256_andnot_si256(x, ONE256)
#endif


	// ===========================================================

	/// count genotype sum and number of calls, not requiring 16-aligned p
	COREARRAY_DLL_DEFAULT C_UInt8* vec_u8_geno_count(C_UInt8 *p,
		size_t n, C_Int32 &out_sum, C_Int32 &out_num);

	/// any (*p > 3) is set to be 3
	COREARRAY_DLL_DEFAULT void vec_u8_geno_valid(C_UInt8 *p, size_t n);

	/// add *p by v and applied to all n
	COREARRAY_DLL_DEFAULT void vec_i32_add(C_Int32 *p, size_t n, C_Int32 v);

	/// add *p by v and applied to all n
	COREARRAY_DLL_DEFAULT void vec_f64_add(double *p, size_t n, double v);

	/// subtract *p by v and applied to all n
	COREARRAY_DLL_DEFAULT void vec_f64_sub(double *p, size_t n, double v);

	/// *p = v - *p and applied to all n
	COREARRAY_DLL_DEFAULT void vec_f64_sub2(double *p, size_t n, double v);

	/// add *p by *s and applied to all n
	COREARRAY_DLL_DEFAULT void vec_f64_add(double *p, const double *s, size_t n);

	/// multiply *p by v and applied to all n
	COREARRAY_DLL_DEFAULT void vec_f64_mul(double *p, size_t n, double v);

	/// *p += (*s) * v
	COREARRAY_DLL_DEFAULT double *vec_f64_addmul(double *p, const double *s,
		size_t n, double v);
}

#endif /* _HEADER_VECTORIZATION_ */
