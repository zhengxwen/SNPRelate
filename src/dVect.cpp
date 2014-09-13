// ===========================================================
//     _/_/_/   _/_/_/  _/_/_/_/    _/_/_/_/  _/_/_/   _/_/_/
//      _/    _/       _/             _/    _/    _/   _/   _/
//     _/    _/       _/_/_/_/       _/    _/    _/   _/_/_/
//    _/    _/       _/             _/    _/    _/   _/
// _/_/_/   _/_/_/  _/_/_/_/_/     _/     _/_/_/   _/_/
// ===========================================================
//
// dVect.cpp: Classess and functions for vectorization
//
// Copyright (C) 2007 - 2014	Xiuwen Zheng
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


#include <dVect.h>

#ifdef COREARRAY_PREDEFINED_SIMD
#include <xmmintrin.h>
#include <emmintrin.h>
#endif


#ifdef COREARRAY_PREDEFINED_SIMD

// SSE intrinsics

COREARRAY_INLINE static float _add(const float v1, const float v2) { return v1 + v2; }
COREARRAY_INLINE static float _sub(const float v1, const float v2) { return v1 - v2; }
COREARRAY_INLINE static float _mul(const float v1, const float v2) { return v1 * v2; }
COREARRAY_INLINE static float _div(const float v1, const float v2) { return v1 / v2; }

// Unroll loop

#define SSE_L3_F32(Fd, Fs1, Fs2, Op, OpR) \
	{ \
		while (n >= 16) { \
			Fd(d, Op(Fs1(s1), Fs2(s2))); d += 4; s1 += 4; s2 += 4; \
			Fd(d, Op(Fs1(s1), Fs2(s2))); d += 4; s1 += 4; s2 += 4; \
			Fd(d, Op(Fs1(s1), Fs2(s2))); d += 4; s1 += 4; s2 += 4; \
			Fd(d, Op(Fs1(s1), Fs2(s2))); d += 4; s1 += 4; s2 += 4; \
			n -= 16; \
		} \
		while (n >= 4) { \
			Fd(d, Op(Fs1(s1), Fs2(s2))); d += 4; s1 += 4; s2 += 4; \
			n -= 4; \
		} \
		while (n > 0) { \
			*d++ = OpR(*s1++, *s2++); n--; \
		} \
		return; \
	}

#define OpSSE_L3_F32(Op, OpR) \
	{ \
		bool pd = ((size_t)d) & 0x0F; \
		bool ps1 = ((size_t)s1) & 0x0F; \
		bool ps2 = ((size_t)s2) & 0x0F; \
		if (pd) { \
			if (ps1) { \
				if (ps2) { \
					SSE_L3_F32(_mm_storeu_ps, _mm_loadu_ps, _mm_loadu_ps, Op, OpR); \
				} else { \
					SSE_L3_F32(_mm_storeu_ps, _mm_loadu_ps, _mm_load_ps, Op, OpR); \
				} \
			} else { \
				if (ps2) { \
					SSE_L3_F32(_mm_storeu_ps, _mm_load_ps, _mm_loadu_ps, Op, OpR); \
				} else { \
					SSE_L3_F32(_mm_storeu_ps, _mm_load_ps, _mm_load_ps, Op, OpR); \
				} \
			} \
		} else { \
			if (ps1) { \
				if (ps2) { \
					SSE_L3_F32(_mm_store_ps, _mm_loadu_ps, _mm_loadu_ps, Op, OpR); \
				} else { \
					SSE_L3_F32(_mm_store_ps, _mm_loadu_ps, _mm_load_ps, Op, OpR); \
				} \
			} else { \
				if (ps2) { \
					SSE_L3_F32(_mm_store_ps, _mm_load_ps, _mm_loadu_ps, Op, OpR); \
				} else { \
					SSE_L3_F32(_mm_store_ps, _mm_load_ps, _mm_load_ps, Op, OpR); \
				} \
			} \
		} \
	}

#define OpSSE_L3_F32_A(Op, OpR) \
	{ \
		SSE_L3_F32(_mm_store_ps, _mm_load_ps, _mm_load_ps, Op, OpR); \
	}


#define SSE_L2_F32(Fd, Fs, Op, OpR) \
	{ \
		while (n >= 16) { \
			Fd(d, Op(Fs(s), rv128)); d += 4; s += 4; \
			Fd(d, Op(Fs(s), rv128)); d += 4; s += 4; \
			Fd(d, Op(Fs(s), rv128)); d += 4; s += 4; \
			Fd(d, Op(Fs(s), rv128)); d += 4; s += 4; \
			n -= 16; \
		} \
		while (n >= 4) { \
			Fd(d, Op(Fs(s), rv128)); d += 4; s += 4; \
			n -= 4; \
		} \
		while (n > 0) { \
			*d++ = OpR(*s++, v); n--; \
		} \
		return; \
	}

#define OpSSE_L2_F32(Op, OpR) \
	{ \
		float rv4[4] __attribute__((aligned(16))) = {v, v, v, v}; \
		__m128 rv128 = _mm_load_ps(rv4); \
		bool pd = ((size_t)d) & 0x0F; \
		bool ps = ((size_t)s) & 0x0F; \
		if (pd) { \
			if (ps) { \
				SSE_L2_F32(_mm_storeu_ps, _mm_loadu_ps, Op, OpR); \
			} else { \
				SSE_L2_F32(_mm_storeu_ps, _mm_load_ps, Op, OpR); \
			} \
		} else { \
			if (ps) { \
				SSE_L2_F32(_mm_store_ps, _mm_loadu_ps, Op, OpR); \
			} else { \
				SSE_L2_F32(_mm_store_ps, _mm_load_ps, Op, OpR); \
			} \
		} \
	}

#define OpSSE_L2_F32_A(Op, OpR) \
	{ \
		float rv4[4] __attribute__((aligned(16))) = {v, v, v, v}; \
		__m128 rv128 = _mm_load_ps(rv4); \
		SSE_L2_F32(_mm_store_ps, _mm_load_ps, Op, OpR); \
	}


#define SSE_L2v_F32(Fd, Fs, Op, OpR) \
	{ \
		while (n >= 16) { \
			Fd(d, Op(rv128, Fs(s))); d += 4; s += 4; \
			Fd(d, Op(rv128, Fs(s))); d += 4; s += 4; \
			Fd(d, Op(rv128, Fs(s))); d += 4; s += 4; \
			Fd(d, Op(rv128, Fs(s))); d += 4; s += 4; \
			n -= 16; \
		} \
		while (n >= 4) { \
			Fd(d, Op(rv128, Fs(s))); d += 4; s += 4; \
			n -= 4; \
		} \
		while (n > 0) { \
			*d++ = OpR(v, *s++); n--; \
		} \
		return; \
	}

#define OpSSE_L2v_F32(Op, OpR) \
	{ \
		float rv4[4] __attribute__((aligned(16))) = {v, v, v, v}; \
		__m128 rv128 = _mm_load_ps(rv4); \
		bool pd = ((size_t)d) & 0x0F; \
		bool ps = ((size_t)s) & 0x0F; \
		if (pd) { \
			if (ps) { \
				SSE_L2v_F32(_mm_storeu_ps, _mm_loadu_ps, Op, OpR); \
			} else { \
				SSE_L2v_F32(_mm_storeu_ps, _mm_load_ps, Op, OpR); \
			} \
		} else { \
			if (ps) { \
				SSE_L2v_F32(_mm_store_ps, _mm_loadu_ps, Op, OpR); \
			} else { \
				SSE_L2v_F32(_mm_store_ps, _mm_load_ps, Op, OpR); \
			} \
		} \
	}

#define OpSSE_L2v_F32_A(Op, OpR) \
	{ \
		float rv4[4] __attribute__((aligned(16))) = {v, v, v, v}; \
		__m128 rv128 = _mm_load_ps(rv4); \
		SSE_L2v_F32(_mm_store_ps, _mm_load_ps, Op, OpR); \
	}


// SSE general functions

// Add
void CORESSECALL CoreArray::Vectorization::_SSE_Add(float *d, const float *s1,
	const float *s2, size_t n) // d := s1 + s2
{
	OpSSE_L3_F32(_mm_add_ps, _add);
}

void CORESSECALL CoreArray::Vectorization::_SSE_Add(float *d, const float *s,
	const float v, size_t n) // d := s + v
{
	OpSSE_L2_F32(_mm_add_ps, _add);
}

// Sub
void CORESSECALL CoreArray::Vectorization::_SSE_Sub(float *d, const float *s1,
	const float *s2, size_t n) // d := s1 - s2
{
	OpSSE_L3_F32(_mm_sub_ps, _sub);
}

void CORESSECALL CoreArray::Vectorization::_SSE_Sub(float *d, const float *s,
	const float v, size_t n) // d := s - v
{
	OpSSE_L2_F32(_mm_sub_ps, _sub);
}

void CORESSECALL CoreArray::Vectorization::_SSE_Sub(float *d, const float v,
	const float *s, size_t n) // d := v - s
{
	OpSSE_L2v_F32(_mm_sub_ps, _sub);
}

// Mul
void CORESSECALL CoreArray::Vectorization::_SSE_Mul(float *d, const float *s1,
	const float *s2, size_t n) // d := s1 * s2
{
	OpSSE_L3_F32(_mm_mul_ps, _mul);
}

void CORESSECALL CoreArray::Vectorization::_SSE_Mul(float *d, const float *s,
	const float v, size_t n) // d := s * scale
{
	OpSSE_L2_F32(_mm_mul_ps, _mul);
}

// Div
void CORESSECALL CoreArray::Vectorization::_SSE_Div(float *d, const float *s1,
	const float *s2, size_t n) // d := s1 / s2
{
	OpSSE_L3_F32(_mm_div_ps, _div);
}

void CORESSECALL CoreArray::Vectorization::_SSE_Div(float *d, const float *s,
	const float v, size_t n) // d := s / v
{
	OpSSE_L2_F32(_mm_div_ps, _div);
}

void CORESSECALL CoreArray::Vectorization::_SSE_Div(float *d, const float v,
	const float *s, size_t n) // d := v / s
{
	OpSSE_L2v_F32(_mm_div_ps, _div);
}

// Dot
#define DOT_L2_F32(Fx, Fy) \
	{ \
		while (n >= 16) { \
			rv128 = _mm_add_ps(rv128, _mm_mul_ps(Fx(x), Fy(y))); x += 4; y += 4; \
			rv128 = _mm_add_ps(rv128, _mm_mul_ps(Fx(x), Fy(y))); x += 4; y += 4; \
			rv128 = _mm_add_ps(rv128, _mm_mul_ps(Fx(x), Fy(y))); x += 4; y += 4; \
			rv128 = _mm_add_ps(rv128, _mm_mul_ps(Fx(x), Fy(y))); x += 4; y += 4; \
			n -= 16; \
		} \
		while (n >= 4) { \
			rv128 = _mm_add_ps(rv128, _mm_mul_ps(Fx(x), Fy(y))); x += 4; y += 4; \
			n -= 4; \
		} \
		_mm_store_ps(rv4, rv128); \
		rv4[0] += rv4[1] + rv4[2] + rv4[3]; \
		while (n > 0) { \
			rv4[0] += (*x++) * (*y++); n--; \
		} \
		return rv4[0]; \
	}

// SSE 16-align functions

// Add
void CORESSECALL CoreArray::Vectorization::_SSE_Add_16(float *d, const float *s1,
	const float *s2, size_t n) // d := s1 + s2
{
	OpSSE_L3_F32_A(_mm_add_ps, _add);
}

void CORESSECALL CoreArray::Vectorization::_SSE_Add_16(float *d, const float *s,
	const float v, size_t n) // d := s + v
{
	OpSSE_L2_F32_A(_mm_add_ps, _add);
}

// Sub
void CORESSECALL CoreArray::Vectorization::_SSE_Sub_16(float *d, const float *s1,
	const float *s2, size_t n) // d := s1 - s2
{
	OpSSE_L3_F32_A(_mm_sub_ps, _sub);
}

void CORESSECALL CoreArray::Vectorization::_SSE_Sub_16(float *d, const float *s,
	const float v, size_t n) // d := s - v
{
	OpSSE_L2_F32_A(_mm_sub_ps, _sub);
}

void CORESSECALL CoreArray::Vectorization::_SSE_Sub_16(float *d, const float v,
	const float *s, size_t n) // d := v - s
{
	OpSSE_L2v_F32_A(_mm_sub_ps, _sub);
}

// Mul
void CORESSECALL CoreArray::Vectorization::_SSE_Mul_16(float *d, const float *s1,
	const float *s2, size_t n) // d := s1 * s2
{
	OpSSE_L3_F32_A(_mm_mul_ps, _mul);
}

void CORESSECALL CoreArray::Vectorization::_SSE_Mul_16(float *d, const float *s,
	const float v, size_t n) // d := s * scale
{
	OpSSE_L2_F32_A(_mm_mul_ps, _mul);
}

// Div
void CORESSECALL CoreArray::Vectorization::_SSE_Div_16(float *d, const float *s1,
	const float *s2, size_t n) // d := s1 / s2
{
	OpSSE_L3_F32_A(_mm_div_ps, _div);
}

void CORESSECALL CoreArray::Vectorization::_SSE_Div_16(float *d, const float *s,
	const float v, size_t n) // d := s / v
{
	OpSSE_L2_F32_A(_mm_div_ps, _div);
}

void CORESSECALL CoreArray::Vectorization::_SSE_Div_16(float *d, const float v,
	const float *s, size_t n) // d := v / s
{
	OpSSE_L2v_F32_A(_mm_div_ps, _div);
}

#define DOT_L2_F32(Fx, Fy) \
	{ \
		while (n >= 16) { \
			rv128 = _mm_add_ps(rv128, _mm_mul_ps(Fx(x), Fy(y))); x += 4; y += 4; \
			rv128 = _mm_add_ps(rv128, _mm_mul_ps(Fx(x), Fy(y))); x += 4; y += 4; \
			rv128 = _mm_add_ps(rv128, _mm_mul_ps(Fx(x), Fy(y))); x += 4; y += 4; \
			rv128 = _mm_add_ps(rv128, _mm_mul_ps(Fx(x), Fy(y))); x += 4; y += 4; \
			n -= 16; \
		} \
		while (n >= 4) { \
			rv128 = _mm_add_ps(rv128, _mm_mul_ps(Fx(x), Fy(y))); x += 4; y += 4; \
			n -= 4; \
		} \
		_mm_store_ps(rv4, rv128); \
		rv4[0] += rv4[1] + rv4[2] + rv4[3]; \
		while (n > 0) { \
			rv4[0] += (*x++) * (*y++); n--; \
		} \
		return rv4[0]; \
	}

float CORESSECALL CoreArray::Vectorization::_SSE_DotProd_16(const float *x, const float *y, size_t n)
{
	float rv4[4] __attribute__((aligned(16))) = {0, 0, 0, 0};
	__m128 rv128 = _mm_load_ps(rv4);
	DOT_L2_F32(_mm_load_ps, _mm_load_ps);
}

// SSE2 intrinsics

COREARRAY_INLINE static double _add(const double v1, const double v2) { return v1 + v2; }
COREARRAY_INLINE static double _sub(const double v1, const double v2) { return v1 - v2; }
COREARRAY_INLINE static double _mul(const double v1, const double v2) { return v1 * v2; }
COREARRAY_INLINE static double _div(const double v1, const double v2) { return v1 / v2; }

#define SSE2_L3_F64(Fd, Fs1, Fs2, Op, OpR) \
	{ \
		while (n >= 8) { \
			Fd(d, Op(Fs1(s1), Fs2(s2))); d += 2; s1 += 2; s2 += 2; \
			Fd(d, Op(Fs1(s1), Fs2(s2))); d += 2; s1 += 2; s2 += 2; \
			Fd(d, Op(Fs1(s1), Fs2(s2))); d += 2; s1 += 2; s2 += 2; \
			Fd(d, Op(Fs1(s1), Fs2(s2))); d += 2; s1 += 2; s2 += 2; \
			n -= 8; \
		} \
		while (n >= 2) { \
			Fd(d, Op(Fs1(s1), Fs2(s2))); d += 2; s1 += 2; s2 += 2; \
			n -= 2; \
		} \
		if (n > 0) *d = OpR(*s1, *s2); \
		return; \
	}

#define OpSSE2_L3_F64(Op, OpR) \
	{ \
		bool pd = ((size_t)d) & 0x0F; \
		bool ps1 = ((size_t)s1) & 0x0F; \
		bool ps2 = ((size_t)s2) & 0x0F; \
		if (pd) { \
			if (ps1) { \
				if (ps2) { \
					SSE2_L3_F64(_mm_storeu_pd, _mm_loadu_pd, _mm_loadu_pd, Op, OpR); \
				} else { \
					SSE2_L3_F64(_mm_storeu_pd, _mm_loadu_pd, _mm_load_pd, Op, OpR); \
				} \
			} else { \
				if (ps2) { \
					SSE2_L3_F64(_mm_storeu_pd, _mm_load_pd, _mm_loadu_pd, Op, OpR); \
				} else { \
					SSE2_L3_F64(_mm_storeu_pd, _mm_load_pd, _mm_load_pd, Op, OpR); \
				} \
			} \
		} else { \
			if (ps1) { \
				if (ps2) { \
					SSE2_L3_F64(_mm_store_pd, _mm_loadu_pd, _mm_loadu_pd, Op, OpR); \
				} else { \
					SSE2_L3_F64(_mm_store_pd, _mm_loadu_pd, _mm_load_pd, Op, OpR); \
				} \
			} else { \
				if (ps2) { \
					SSE2_L3_F64(_mm_store_pd, _mm_load_pd, _mm_loadu_pd, Op, OpR); \
				} else { \
					SSE2_L3_F64(_mm_store_pd, _mm_load_pd, _mm_load_pd, Op, OpR); \
				} \
			} \
		} \
	}

#define OpSSE2_L3_F64_A(Op, OpR) \
	{ \
		SSE2_L3_F64(_mm_store_pd, _mm_load_pd, _mm_load_pd, Op, OpR); \
	}


#define SSE2_L2_F64(Fd, Fs, Op, OpR) \
	{ \
		while (n >= 8) { \
			Fd(d, Op(Fs(s), rv128)); d += 2; s += 2; \
			Fd(d, Op(Fs(s), rv128)); d += 2; s += 2; \
			Fd(d, Op(Fs(s), rv128)); d += 2; s += 2; \
			Fd(d, Op(Fs(s), rv128)); d += 2; s += 2; \
			n -= 8; \
		} \
		while (n >= 2) { \
			Fd(d, Op(Fs(s), rv128)); d += 2; s += 2; \
			n -= 2; \
		} \
		if (n > 0) *d = OpR(*s, v); \
		return; \
	}

#define OpSSE2_L2_F64(Op, OpR) \
	{ \
		double rv2[2] __attribute__((aligned(16))) = {v, v}; \
		__m128d rv128 = _mm_load_pd(rv2); \
		bool pd = ((size_t)d) & 0x0F; \
		bool ps = ((size_t)s) & 0x0F; \
		if (pd) { \
			if (ps) { \
				SSE2_L2_F64(_mm_storeu_pd, _mm_loadu_pd, Op, OpR); \
			} else { \
				SSE2_L2_F64(_mm_storeu_pd, _mm_load_pd, Op, OpR); \
			} \
		} else { \
			if (ps) { \
				SSE2_L2_F64(_mm_store_pd, _mm_loadu_pd, Op, OpR); \
			} else { \
				SSE2_L2_F64(_mm_store_pd, _mm_load_pd, Op, OpR); \
			} \
		} \
	}

#define OpSSE2_L2_F64_A(Op, OpR) \
	{ \
		double rv2[2] __attribute__((aligned(16))) = {v, v}; \
		__m128d rv128 = _mm_load_pd(rv2); \
		SSE2_L2_F64(_mm_store_pd, _mm_load_pd, Op, OpR); \
	}


#define SSE2_L2v_F64(Fd, Fs, Op, OpR) \
	{ \
		while (n >= 8) { \
			Fd(d, Op(rv128, Fs(s))); d += 2; s += 2; \
			Fd(d, Op(rv128, Fs(s))); d += 2; s += 2; \
			Fd(d, Op(rv128, Fs(s))); d += 2; s += 2; \
			Fd(d, Op(rv128, Fs(s))); d += 2; s += 2; \
			n -= 8; \
		} \
		while (n >= 2) { \
			Fd(d, Op(rv128, Fs(s))); d += 2; s += 2; \
			n -= 2; \
		} \
		if (n > 0) *d = OpR(v, *s); \
		return; \
	}

#define OpSSE2_L2v_F64(Op, OpR) \
	{ \
		double rv2[2] __attribute__((aligned(16))) = {v, v}; \
		__m128d rv128 = _mm_load_pd(rv2); \
		bool pd = ((size_t)d) & 0x0F; \
		bool ps = ((size_t)s) & 0x0F; \
		if (pd) { \
			if (ps) { \
				SSE2_L2v_F64(_mm_storeu_pd, _mm_loadu_pd, Op, OpR); \
			} else { \
				SSE2_L2v_F64(_mm_storeu_pd, _mm_load_pd, Op, OpR); \
			} \
		} else { \
			if (ps) { \
				SSE2_L2v_F64(_mm_store_pd, _mm_loadu_pd, Op, OpR); \
			} else { \
				SSE2_L2v_F64(_mm_store_pd, _mm_load_pd, Op, OpR); \
			} \
		} \
	}

#define OpSSE2_L2v_F64_A(Op, OpR) \
	{ \
		double rv2[2] __attribute__((aligned(16))) = {v, v}; \
		__m128d rv128 = _mm_load_pd(rv2); \
		SSE2_L2v_F64(_mm_store_pd, _mm_load_pd, Op, OpR); \
	}

// SSE2 16-align functions

// Add

void CORESSECALL CoreArray::Vectorization::_SSE2_Add_16(double *d,
	const double *s1, const double *s2, size_t n) // d := s1 + s2
{
	OpSSE2_L3_F64_A(_mm_add_pd, _add);
}

void CORESSECALL CoreArray::Vectorization::_SSE2_Add_16(double *d,
	const double *s, const double v, size_t n) // d := s + v
{
	OpSSE2_L2_F64_A(_mm_add_pd, _add);
}

// Sub

void CORESSECALL CoreArray::Vectorization::_SSE2_Sub_16(double *d,
	const double *s1, const double *s2, size_t n) // d := s1 - s2
{
	OpSSE2_L3_F64_A(_mm_sub_pd, _sub);
}

void CORESSECALL CoreArray::Vectorization::_SSE2_Sub_16(double *d,
	const double *s, const double v, size_t n) // d := s - v
{
	OpSSE2_L2_F64_A(_mm_sub_pd, _sub);
}

void CORESSECALL CoreArray::Vectorization::_SSE2_Sub_16(double *d,
	const double v, const double *s, size_t n) // d := v - s
{
	OpSSE2_L2v_F64_A(_mm_sub_pd, _sub);
}

// Mul

void CORESSECALL CoreArray::Vectorization::_SSE2_Mul_16(double *d,
	const double *s1, const double *s2, size_t n) // d := s1 * s2
{
	OpSSE2_L3_F64_A(_mm_mul_pd, _mul);
}

void CORESSECALL CoreArray::Vectorization::_SSE2_Mul_16(double *d,
	const double *s, const double v, size_t n) // d := s * scale
{
	OpSSE2_L2_F64_A(_mm_mul_pd, _mul);
}

// Div

void CORESSECALL CoreArray::Vectorization::_SSE2_Div_16(double *d,
	const double *s1, const double *s2, size_t n) // d := s1 / s2
{
	OpSSE2_L3_F64_A(_mm_div_pd, _div);
}

void CORESSECALL CoreArray::Vectorization::_SSE2_Div_16(double *d,
	const double *s, const double v, size_t n) // d := s / v
{
	OpSSE2_L2_F64_A(_mm_div_pd, _div);
}

void CORESSECALL CoreArray::Vectorization::_SSE2_Div_16(double *d,
	const double v, const double *s, size_t n) // d := v / s
{
	OpSSE2_L2v_F64_A(_mm_div_pd, _div);
}

#define DOT_L2_F64(Fx, Fy) \
	{ \
		while (n >= 8) { \
			rv128 = _mm_add_pd(rv128, _mm_mul_pd(Fx(x), Fy(y))); x += 2; y += 2; \
			rv128 = _mm_add_pd(rv128, _mm_mul_pd(Fx(x), Fy(y))); x += 2; y += 2; \
			rv128 = _mm_add_pd(rv128, _mm_mul_pd(Fx(x), Fy(y))); x += 2; y += 2; \
			rv128 = _mm_add_pd(rv128, _mm_mul_pd(Fx(x), Fy(y))); x += 2; y += 2; \
			n -= 8; \
		} \
		while (n >= 2) { \
			rv128 = _mm_add_pd(rv128, _mm_mul_pd(Fx(x), Fy(y))); x += 2; y += 2; \
			n -= 2; \
		} \
		_mm_store_pd(rv2, rv128); \
		rv2[0] += rv2[1]; \
		if (n > 0) { \
			rv2[0] += (*x) * (*y); \
		} \
		return rv2[0]; \
	}

double CORESSECALL CoreArray::Vectorization::_SSE2_DotProd_16(const double *x,
	const double *y, size_t n)
{
	double rv2[2] __attribute__((aligned(16))) = {0, 0};
	__m128d rv128 = _mm_load_pd(rv2);
	DOT_L2_F64(_mm_load_pd, _mm_load_pd);
}

#endif
