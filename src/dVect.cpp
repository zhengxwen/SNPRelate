// ===========================================================
//
// dVect.cpp: Classess and functions for vectorization
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
// You should have received a copy of the GNU Lesser General Public
// License along with SNPRelate.
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
void CORESSECALL Vectorization::_SSE_Add(float *d, const float *s1,
	const float *s2, size_t n) // d := s1 + s2
{
	OpSSE_L3_F32(_mm_add_ps, _add);
}

void CORESSECALL Vectorization::_SSE_Add(float *d, const float *s,
	const float v, size_t n) // d := s + v
{
	OpSSE_L2_F32(_mm_add_ps, _add);
}

// Sub
void CORESSECALL Vectorization::_SSE_Sub(float *d, const float *s1,
	const float *s2, size_t n) // d := s1 - s2
{
	OpSSE_L3_F32(_mm_sub_ps, _sub);
}

void CORESSECALL Vectorization::_SSE_Sub(float *d, const float *s,
	const float v, size_t n) // d := s - v
{
	OpSSE_L2_F32(_mm_sub_ps, _sub);
}

void CORESSECALL Vectorization::_SSE_Sub(float *d, const float v,
	const float *s, size_t n) // d := v - s
{
	OpSSE_L2v_F32(_mm_sub_ps, _sub);
}

// Mul
void CORESSECALL Vectorization::_SSE_Mul(float *d, const float *s1,
	const float *s2, size_t n) // d := s1 * s2
{
	OpSSE_L3_F32(_mm_mul_ps, _mul);
}

void CORESSECALL Vectorization::_SSE_Mul(float *d, const float *s,
	const float v, size_t n) // d := s * scale
{
	OpSSE_L2_F32(_mm_mul_ps, _mul);
}

// Div
void CORESSECALL Vectorization::_SSE_Div(float *d, const float *s1,
	const float *s2, size_t n) // d := s1 / s2
{
	OpSSE_L3_F32(_mm_div_ps, _div);
}

void CORESSECALL Vectorization::_SSE_Div(float *d, const float *s,
	const float v, size_t n) // d := s / v
{
	OpSSE_L2_F32(_mm_div_ps, _div);
}

void CORESSECALL Vectorization::_SSE_Div(float *d, const float v,
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
void CORESSECALL Vectorization::_SSE_Add_16(float *d, const float *s1,
	const float *s2, size_t n) // d := s1 + s2
{
	OpSSE_L3_F32_A(_mm_add_ps, _add);
}

void CORESSECALL Vectorization::_SSE_Add_16(float *d, const float *s,
	const float v, size_t n) // d := s + v
{
	OpSSE_L2_F32_A(_mm_add_ps, _add);
}

// Sub
void CORESSECALL Vectorization::_SSE_Sub_16(float *d, const float *s1,
	const float *s2, size_t n) // d := s1 - s2
{
	OpSSE_L3_F32_A(_mm_sub_ps, _sub);
}

void CORESSECALL Vectorization::_SSE_Sub_16(float *d, const float *s,
	const float v, size_t n) // d := s - v
{
	OpSSE_L2_F32_A(_mm_sub_ps, _sub);
}

void CORESSECALL Vectorization::_SSE_Sub_16(float *d, const float v,
	const float *s, size_t n) // d := v - s
{
	OpSSE_L2v_F32_A(_mm_sub_ps, _sub);
}

// Mul
void CORESSECALL Vectorization::_SSE_Mul_16(float *d, const float *s1,
	const float *s2, size_t n) // d := s1 * s2
{
	OpSSE_L3_F32_A(_mm_mul_ps, _mul);
}

void CORESSECALL Vectorization::_SSE_Mul_16(float *d, const float *s,
	const float v, size_t n) // d := s * scale
{
	OpSSE_L2_F32_A(_mm_mul_ps, _mul);
}

// Div
void CORESSECALL Vectorization::_SSE_Div_16(float *d, const float *s1,
	const float *s2, size_t n) // d := s1 / s2
{
	OpSSE_L3_F32_A(_mm_div_ps, _div);
}

void CORESSECALL Vectorization::_SSE_Div_16(float *d, const float *s,
	const float v, size_t n) // d := s / v
{
	OpSSE_L2_F32_A(_mm_div_ps, _div);
}

void CORESSECALL Vectorization::_SSE_Div_16(float *d, const float v,
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

float CORESSECALL Vectorization::_SSE_DotProd_16(const float *x, const float *y, size_t n)
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

void CORESSECALL Vectorization::_SSE2_Add_16(double *d,
	const double *s1, const double *s2, size_t n) // d := s1 + s2
{
	OpSSE2_L3_F64_A(_mm_add_pd, _add);
}

void CORESSECALL Vectorization::_SSE2_Add_16(double *d,
	const double *s, const double v, size_t n) // d := s + v
{
	OpSSE2_L2_F64_A(_mm_add_pd, _add);
}

// Sub

void CORESSECALL Vectorization::_SSE2_Sub_16(double *d,
	const double *s1, const double *s2, size_t n) // d := s1 - s2
{
	OpSSE2_L3_F64_A(_mm_sub_pd, _sub);
}

void CORESSECALL Vectorization::_SSE2_Sub_16(double *d,
	const double *s, const double v, size_t n) // d := s - v
{
	OpSSE2_L2_F64_A(_mm_sub_pd, _sub);
}

void CORESSECALL Vectorization::_SSE2_Sub_16(double *d,
	const double v, const double *s, size_t n) // d := v - s
{
	OpSSE2_L2v_F64_A(_mm_sub_pd, _sub);
}

// Mul

void CORESSECALL Vectorization::_SSE2_Mul_16(double *d,
	const double *s1, const double *s2, size_t n) // d := s1 * s2
{
	OpSSE2_L3_F64_A(_mm_mul_pd, _mul);
}

void CORESSECALL Vectorization::_SSE2_Mul_16(double *d,
	const double *s, const double v, size_t n) // d := s * scale
{
	OpSSE2_L2_F64_A(_mm_mul_pd, _mul);
}

// Div

void CORESSECALL Vectorization::_SSE2_Div_16(double *d,
	const double *s1, const double *s2, size_t n) // d := s1 / s2
{
	OpSSE2_L3_F64_A(_mm_div_pd, _div);
}

void CORESSECALL Vectorization::_SSE2_Div_16(double *d,
	const double *s, const double v, size_t n) // d := s / v
{
	OpSSE2_L2_F64_A(_mm_div_pd, _div);
}

void CORESSECALL Vectorization::_SSE2_Div_16(double *d,
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

double CORESSECALL Vectorization::_SSE2_DotProd_16(const double *x,
	const double *y, size_t n)
{
	double rv2[2] __attribute__((aligned(16))) = {0, 0};
	__m128d rv128 = _mm_load_pd(rv2);
	DOT_L2_F64(_mm_load_pd, _mm_load_pd);
}

#endif



// =========================================================================

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


/// any (*p > 2) is set to be 3
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
