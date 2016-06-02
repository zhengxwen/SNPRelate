// ===========================================================
//
// dVect.h: Classess and functions for vectorization
//
// Copyright (C) 2007-2016    Xiuwen Zheng
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

/**
 *	\file     dVect.h
 *	\author   Xiuwen Zheng [zhengxwen@gmail.com]
 *	\version  1.0
 *	\date     2007 - 2015
 *	\brief    Classess and functions for vectorization
 *	\details
**/

#ifndef _HEADER_VECTORIZATION_
#define _HEADER_VECTORIZATION_

#include <dType.h>

#ifdef COREARRAY_SIMD_SSE
#include <xmmintrin.h>
#endif
#ifdef COREARRAY_SIMD_SSE2
#include <emmintrin.h>
#endif
#ifdef COREARRAY_SIMD_AVX
#include <immintrin.h>
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


	// Vectorization Functions

	enum TFlagVectorization
	{
		vtFPU,
		vtSSE,     //< Streaming SIMD Extensions
		vtSSE2,    //< Streaming SIMD Extensions 2
		vtSSE3,    //< Streaming SIMD Extensions 3
		vtSSE4,    //< Streaming SIMD Extensions 4
		vtAVX,     //< Advanced Vector Extensions
		vtAVX2     //< Advanced Vector Extensions 2
	};

	enum TAlignVectorization { avNormal, av16Align };


#if defined(COREARRAY_SIMD_AVX2)
	const TFlagVectorization FLAG_VECTORIZATION = vtAVX2;
#elif defined(COREARRAY_SIMD_AVX)
	const TFlagVectorization FLAG_VECTORIZATION = vtAVX;
#elif defined(COREARRAY_SIMD_SSE4)
	const TFlagVectorization FLAG_VECTORIZATION = vtSSE4;
#elif defined(COREARRAY_SIMD_SSE3)
	const TFlagVectorization FLAG_VECTORIZATION = vtSSE3;
#elif defined(COREARRAY_SIMD_SSE2)
	const TFlagVectorization FLAG_VECTORIZATION = vtSSE2;
#elif defined(COREARRAY_SIMD_SSE)
	const TFlagVectorization FLAG_VECTORIZATION = vtSSE;
#else
	const TFlagVectorization FLAG_VECTORIZATION = vtFPU;
#endif

	template<typename Tx,
		TAlignVectorization fAlign = avNormal,
		bool SSE = (FLAG_VECTORIZATION >= vtSSE),
		bool SSE2 = (FLAG_VECTORIZATION >= vtSSE2),
		bool SSE3 = (FLAG_VECTORIZATION >= vtSSE3),
		bool SSE4 = (FLAG_VECTORIZATION >= vtSSE4)
	>
	struct vt
	{
		static const TAlignVectorization Align = fAlign;
		static const bool useSSE = SSE;
		static const bool useSSE2 = SSE2;
		static const bool useSSE3 = SSE3;
		static const bool useSSE4 = SSE4;

		// Add
		static void Add(Tx *d, const Tx *s1, const Tx *s2, size_t n) // d := s1 + s2
			{ while (n--) *d++ = (*s1++) + (*s2++); };
		static void Add(Tx *d, const Tx *s, const Tx v, size_t n) // d := s + v
			{ while (n--) *d++ = (*s++) + v; };
		// Sub
		static void Sub(Tx *d, const Tx *s1, const Tx *s2, size_t n) // d := s1 - s2
			{ while (n--) *d++ = (*s1++) - (*s2++); };
		static void Sub(Tx *d, const Tx *s, const Tx v, size_t n) // d := s - v
			{ while (n--) *d++ = (*s++) - v; };
		static void Sub(Tx *d, const Tx v, const Tx *s, size_t n) // d := v - s
			{ while (n--) *d++ = v - (*s++); };
		// Mul
		static void Mul(Tx *d, const Tx *s1, const Tx *s2, size_t n) // d := s1 * s2
			{ while (n--) *d++ = (*s1++) * (*s2++); };
		static void Mul(Tx *d, const Tx *s, const Tx scale, size_t n) // d := s * scale
			{ while (n--) *d++ = (*s++) * scale; };
		// Div
		static void Div(Tx *d, const Tx *s1, const Tx *s2, size_t n) // d := s1 / s2
			{ while (n--) *d++ = (*s1++) / (*s2++); };
		static void Div(Tx *d, const Tx *s, const Tx v, size_t n) // d := s / v
			{ while (n--) *d++ = (*s++) / v; };
		static void Div(Tx *d, const Tx v, const Tx *s, size_t n) // d := v / s
			{ while (n--) *d++ = v / (*s++); };

		// Sum
		static Tx Sum(const Tx *s1, size_t n)
		{
			Tx sum = 0;
			while (n--) sum += *s1++;
			return sum;
		}

		// Dot product
		static Tx DotProd(const Tx *s1, const Tx *s2, size_t n)
		{
			Tx sum = 0;
			while (n--) sum += (*s1++) * (*s2++);
			return sum;
		}
	};



// SSE/SSE2 optimization

#ifdef COREARRAY_PREDEFINED_SIMD

	#if defined(COREARRAY_CC_BORLAND) || defined(COREARRAY_CC_MSC)
		#define CORESSECALL __fastcall
	#else
		#define CORESSECALL
	#endif

#ifdef COREARRAY_SIMD_SSE
	// SSE
	// general functions

	// Add
	void CORESSECALL _SSE_Add(float *d, const float *s1, const float *s2, size_t n); // d := s1 + s2
	void CORESSECALL _SSE_Add(float *d, const float *s, const float v, size_t n); // d := s + v

	// Sub
	void CORESSECALL _SSE_Sub(float *d, const float *s1, const float *s2, size_t n); // d := s1 - s2
	void CORESSECALL _SSE_Sub(float *d, const float *s, const float v, size_t n); // d := s - v
	void CORESSECALL _SSE_Sub(float *d, const float v, const float *s, size_t n); // d := v - s

	// Mul
	void CORESSECALL _SSE_Mul(float *d, const float *s1, const float *s2, size_t n); // d := s1 * s2
	void CORESSECALL _SSE_Mul(float *d, const float *s, const float scale, size_t n); // d := s * scale

	// Div
	void CORESSECALL _SSE_Div(float *d, const float *s1, const float *s2, size_t n); // d := s1 / s2
	void CORESSECALL _SSE_Div(float *d, const float *s, const float v, size_t n); // d := s / v
	void CORESSECALL _SSE_Div(float *d, const float v, const float *s, size_t n); // d := v / s

	// Dot product
	float CORESSECALL _SSE_DotProd(const float *x, const float *y, size_t n);

	// 16-align functions

	// Add
	void CORESSECALL _SSE_Add_16(float *d, const float *s1, const float *s2, size_t n); // d := s1 + s2
	void CORESSECALL _SSE_Add_16(float *d, const float *s, const float v, size_t n); // d := s + v

	// Sub
	void CORESSECALL _SSE_Sub_16(float *d, const float *s1, const float *s2, size_t n); // d := s1 - s2
	void CORESSECALL _SSE_Sub_16(float *d, const float *s, const float v, size_t n); // d := s - v
	void CORESSECALL _SSE_Sub_16(float *d, const float v, const float *s, size_t n); // d := v - s

	// Mul
	void CORESSECALL _SSE_Mul_16(float *d, const float *s1, const float *s2, size_t n); // d := s1 * s2
	void CORESSECALL _SSE_Mul_16(float *d, const float *s, const float scale, size_t n); // d := s * scale

	// Div
	void CORESSECALL _SSE_Div_16(float *d, const float *s1, const float *s2, size_t n); // d := s1 / s2
	void CORESSECALL _SSE_Div_16(float *d, const float *s, const float v, size_t n); // d := s / v
	void CORESSECALL _SSE_Div_16(float *d, const float v, const float *s, size_t n); // d := v / s

	// Dot product
	float CORESSECALL _SSE_DotProd_16(const float *x, const float *y, size_t n);
#endif

#ifdef COREARRAY_SIMD_SSE2
	// SSE2
	// 16-align functions

	void CORESSECALL _SSE2_Add_16(double *d, const double *s1, const double *s2, size_t n); // d := s1 + s2
	void CORESSECALL _SSE2_Add_16(double *d, const double *s, const double v, size_t n); // d := s + v

	// Sub
	void CORESSECALL _SSE2_Sub_16(double *d, const double *s1, const double *s2, size_t n); // d := s1 - s2
	void CORESSECALL _SSE2_Sub_16(double *d, const double *s, const double v, size_t n); // d := s - v
	void CORESSECALL _SSE2_Sub_16(double *d, const double v, const double *s, size_t n); // d := v - s

	// Mul
	void CORESSECALL _SSE2_Mul_16(double *d, const double *s1, const double *s2, size_t n); // d := s1 * s2
	void CORESSECALL _SSE2_Mul_16(double *d, const double *s, const double scale, size_t n); // d := s * scale

	// Div
	void CORESSECALL _SSE2_Div_16(double *d, const double *s1, const double *s2, size_t n); // d := s1 / s2
	void CORESSECALL _SSE2_Div_16(double *d, const double *s, const double v, size_t n); // d := s / v
	void CORESSECALL _SSE2_Div_16(double *d, const double v, const double *s, size_t n); // d := v / s

	// Dot product
	double CORESSECALL _SSE2_DotProd_16(const double *x, const double *y, size_t n);
#endif



	// Template

#ifdef COREARRAY_SIMD_SSE
	template<bool SSE2, bool SSE3, bool SSE4>
		struct vt<float, avNormal, true, SSE2, SSE3, SSE4>
	{
		static const TAlignVectorization Align = avNormal;
		static const bool useSSE = true;
		static const bool useSSE2 = SSE2;
		static const bool useSSE3 = SSE3;
		static const bool useSSE4 = SSE4;

		// Add
		COREARRAY_INLINE static void Add(float *d, const float *s1, const float *s2, size_t n) // d := s1 + s2
			{ _SSE_Add(d, s1, s2, n); }
		COREARRAY_INLINE static void Add(float *d, const float *s, const float v, size_t n) // d := s + v
			{ _SSE_Add(d, s, v, n); }

		// Sub
		COREARRAY_INLINE static void Sub(float *d, const float *s1, const float *s2, size_t n) // d := s1 - s2
			{ _SSE_Sub(d, s1, s2, n); }
		COREARRAY_INLINE static void Sub(float *d, const float *s, const float v, size_t n) // d := s - v
			{ _SSE_Sub(d, s, v, n); }
		COREARRAY_INLINE static void Sub(float *d, const float v, const float *s, size_t n) // d := v - s
			{ _SSE_Sub(d, v, s, n); }

		// Mul
		COREARRAY_INLINE static void Mul(float *d, const float *s1, const float *s2, size_t n) // d := s1 * s2
			{ _SSE_Mul(d, s1, s2, n); }
		COREARRAY_INLINE static void Mul(float *d, const float *s, const float scale, size_t n) // d := s * scale
			{ _SSE_Mul(d, s, scale, n); }

		// Div
		COREARRAY_INLINE static void Div(float *d, const float *s1, const float *s2, size_t n) // d := s1 / s2
			{ _SSE_Div(d, s1, s2, n); }
		COREARRAY_INLINE static void Div(float *d, const float *s, const float v, size_t n) // d := s / v
			{ _SSE_Div(d, s, v, n); }
		COREARRAY_INLINE static void Div(float *d, const float v, const float *s, size_t n) // d := v / s
			{ _SSE_Div(d, v, s, n); }
	};

	template<bool SSE2, bool SSE3, bool SSE4>
		struct vt<float, av16Align, true, SSE2, SSE3, SSE4>
	{
		static const TAlignVectorization Align = av16Align;
		static const bool useSSE = true;
		static const bool useSSE2 = SSE2;
		static const bool useSSE3 = SSE3;
		static const bool useSSE4 = SSE4;

		// Add
		COREARRAY_INLINE static void Add(float *d, const float *s1, const float *s2, size_t n) // d := s1 + s2
			{ _SSE_Add_16(d, s1, s2, n); }
		COREARRAY_INLINE static void Add(float *d, const float *s, const float v, size_t n) // d := s + v
			{ _SSE_Add_16(d, s, v, n); }

		// Sub
		COREARRAY_INLINE static void Sub(float *d, const float *s1, const float *s2, size_t n) // d := s1 - s2
			{ _SSE_Sub_16(d, s1, s2, n); }
		COREARRAY_INLINE static void Sub(float *d, const float *s, const float v, size_t n) // d := s - v
			{ _SSE_Sub_16(d, s, v, n); }
		COREARRAY_INLINE static void Sub(float *d, const float v, const float *s, size_t n) // d := v - s
			{ _SSE_Sub_16(d, v, s, n); }

		// Mul
		COREARRAY_INLINE static void Mul(float *d, const float *s1, const float *s2, size_t n) // d := s1 * s2
			{ _SSE_Mul_16(d, s1, s2, n); }
		COREARRAY_INLINE static void Mul(float *d, const float *s, const float scale, size_t n) // d := s * scale
			{ _SSE_Mul_16(d, s, scale, n); }

		// Div
		COREARRAY_INLINE static void Div(float *d, const float *s1, const float *s2, size_t n) // d := s1 / s2
			{ _SSE_Div_16(d, s1, s2, n); }
		COREARRAY_INLINE static void Div(float *d, const float *s, const float v, size_t n) // d := s / v
			{ _SSE_Div_16(d, s, v, n); }
		COREARRAY_INLINE static void Div(float *d, const float v, const float *s, size_t n) // d := v / s
			{ _SSE_Div_16(d, v, s, n); }

		// Dot product
		COREARRAY_INLINE static float DotProd(const float *s1, const float *s2, size_t n)
			{ return _SSE_DotProd_16(s1, s2, n); }
	};
#endif

#ifdef COREARRAY_SIMD_SSE2
	template<bool SSE3, bool SSE4>
		struct vt<double, av16Align, true, true, SSE3, SSE4>
	{
		static const TAlignVectorization Align = av16Align;
		static const bool useSSE = true;
		static const bool useSSE2 = true;
		static const bool useSSE3 = SSE3;
		static const bool useSSE4 = SSE4;

		// Add
		COREARRAY_INLINE static void Add(double *d, const double *s1, const double *s2, size_t n) // d := s1 + s2
			{ _SSE2_Add_16(d, s1, s2, n); }
		COREARRAY_INLINE static void Add(double *d, const double *s, const double v, size_t n) // d := s + v
			{ _SSE2_Add_16(d, s, v, n); }

		// Sub
		COREARRAY_INLINE static void Sub(double *d, const double *s1, const double *s2, size_t n) // d := s1 - s2
			{ _SSE2_Sub_16(d, s1, s2, n); }
		COREARRAY_INLINE static void Sub(double *d, const double *s, const double v, size_t n) // d := s - v
			{ _SSE2_Sub_16(d, s, v, n); }
		COREARRAY_INLINE static void Sub(double *d, const double v, const double *s, size_t n) // d := v - s
			{ _SSE2_Sub_16(d, v, s, n); }

		// Mul
		COREARRAY_INLINE static void Mul(double *d, const double *s1, const double *s2, size_t n) // d := s1 * s2
			{ _SSE2_Mul_16(d, s1, s2, n); }
		COREARRAY_INLINE static void Mul(double *d, const double *s, const double scale, size_t n) // d := s * scale
			{ _SSE2_Mul_16(d, s, scale, n); }

		// Div
		COREARRAY_INLINE static void Div(double *d, const double *s1, const double *s2, size_t n) // d := s1 / s2
			{ _SSE2_Div_16(d, s1, s2, n); }
		COREARRAY_INLINE static void Div(double *d, const double *s, const double v, size_t n) // d := s / v
			{ _SSE2_Div_16(d, s, v, n); }
		COREARRAY_INLINE static void Div(double *d, const double v, const double *s, size_t n) // d := v / s
			{ _SSE2_Div_16(d, v, s, n); }

		// Dot product
		COREARRAY_INLINE static double DotProd(const double *s1, const double *s2, size_t n)
			{ return _SSE2_DotProd_16(s1, s2, n); }
	};
#endif

#endif



	// ===========================================================

#ifdef COREARRAY_SIMD_SSE2
#   ifdef __SSE4_1__
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

	/// count genotype sum and number of calls, not requiring 16-aligned p
	COREARRAY_DLL_DEFAULT C_UInt8* vec_u8_geno_count(C_UInt8 *p,
		size_t n, C_Int32 &out_sum, C_Int32 &out_num);

	/// any (*p > 3) is set to be 3
	COREARRAY_DLL_DEFAULT void vec_u8_geno_valid(C_UInt8 *p, size_t n);

	/// add *p by *s and applied to all n
	COREARRAY_DLL_DEFAULT void vec_f64_add(double *p, const double *s, size_t n);

	/// multiply *p by v and applied to all n
	COREARRAY_DLL_DEFAULT void vec_f64_mul(double *p, size_t n, double v);

	/// *p += (*s) * v
	COREARRAY_DLL_DEFAULT double *vec_f64_addmul(double *p, const double *s,
		size_t n, double v);
}

#endif /* _HEADER_VECTORIZATION_ */
