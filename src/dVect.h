// ===========================================================
//
// dVect.h: Classess and functions for vectorization
//
// Copyright (C) 2007-2015    Xiuwen Zheng
//
// This file is part of CoreArray.
//
// CoreArray is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License Version 3 as
// published by the Free Software Foundation.
//
// CoreArray is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with CoreArray.
// If not, see <http://www.gnu.org/licenses/>.

/**
 *	\file     dVect.h
 *	\author   Xiuwen Zheng [zhengxwen@gmail.com]
 *	\version  1.0
 *	\date     2007 - 2015
 *	\brief    Classess and functions for vectorization
 *	\details
**/

#ifndef _HEADER_COREARRAY_VECT_
#define _HEADER_COREARRAY_VECT_

#include <dType.h>

namespace CoreArray
{
	namespace Vectorization
	{

	#if defined(COREARRAY_SIMD_AVX)
		const size_t _SIMD_ALIGN_ = 32u;
	#else
		const size_t _SIMD_ALIGN_ = 16u;
	#endif

		template<typename Tx, size_t vAlign=_SIMD_ALIGN_> struct TdAlignPtr
		{
		public:
			static const size_t Align = vAlign;

			TdAlignPtr()
				{ alloc = NULL; base = NULL; vlen = 0; }
			TdAlignPtr(size_t n)
				{ alloc = NULL; base = NULL; vlen = 0; Reset(n); }
			~TdAlignPtr()
				{ if (alloc) delete []alloc; }

			void Reset(size_t n)
			{
				if (n != vlen)
				{
					if (alloc) delete []alloc;
					if (n > 0)
					{
						alloc = new char[n*sizeof(Tx) + Align - 1];
						size_t r = ((size_t)alloc) % Align;
						base = (Tx*)(r ? (alloc+Align-r) : alloc);
						vlen = n;
					} else {
						alloc = NULL; base = NULL;
						vlen = 0;
					}
				}
			}

			COREARRAY_INLINE Tx *get() { return base; };
			COREARRAY_INLINE size_t len() { return vlen; };
		private:
			char *alloc;
			Tx *base;
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
	}
}

#endif /* _HEADER_COREARRAY_VECT_ */
