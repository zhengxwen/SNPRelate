// ===========================================================
//
// genPCA.cpp: Principal Component Analysis on GWAS
//
// Copyright (C) 2011-2016    Xiuwen Zheng
//
// This file is part of SNPRelate.
//
// SNPRelate is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License Version 3 as published
// by the Free Software Foundation.
//
// SNPRelate is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with SNPRelate.
// If not, see <http://www.gnu.org/licenses/>.


#ifndef _HEADER_PCA_
#define _HEADER_PCA_

// CoreArray library header
#include "dGenGWAS.h"
#include "dVect.h"
#include "ThreadPool.h"

#include <R_ext/Lapack.h>
#include <cmath>
#include <memory>
#include <algorithm>


namespace PCA
{

// using namespace
using namespace std;
using namespace CoreArray;
using namespace Vectorization;
using namespace GWAS;


// ---------------------------------------------------------------------
// PCA parameters

class CPCAMat_Alg1;

/// whether use Bayesian normalization
bool BayesianNormal = false;
/// The number of eigenvectors output
long OutputEigenDim = 32;
/// the pointer to the output buffer for SNP correlation and SNP/sample loadings
double *Out_Buffer = NULL;
/// the pointer to eigenvectors, snp loadings
double *In_EigenVect = NULL;
/// the pointer to the object storing the eigenvectors, snp loadings
CPCAMat_Alg1 *_EigenVectBuf = NULL;
/// the pointer to allele frequency
double *In_AveFreq = NULL;

/// init mutex objects
void InitMutexObject()
{
	_Mutex = GDS_Parallel_InitMutex();
}
/// destroy mutex objects
void DoneMutexObject()
{
	GDS_Parallel_DoneMutex(_Mutex); _Mutex = NULL;
}



// ---------------------------------------------------------------------
// Vectorization Computing, algorithm 1 (Arithmetic)

class COREARRAY_DLL_LOCAL CPCAMat_Base
{
public:
	VEC_AUTO_PTR<C_Int32> PCA_GenoSum, PCA_GenoNum;
	VEC_AUTO_PTR<double> tmp_var;  // avg geno: 2*\bar{p}, or ...

	CPCAMat_Base() { }

	inline void ZeroFill()
	{
		memset(PCA_GenoSum.Get(), 0, sizeof(C_Int32)*fM);
		memset(PCA_GenoNum.Get(), 0, sizeof(C_Int32)*fM);
	}

	/// calculate average genotypes (save to PCA_GenoSum and PCA_GenoNum)
	void SummarizeGeno_SampxSNP(C_UInt8 *pGeno, size_t nSNP)
	{
		C_Int32 *pS = PCA_GenoSum.Get();
		C_Int32 *pN = PCA_GenoNum.Get();
		for (size_t i=0; i < nSNP; i++)
		{
			pGeno = vec_u8_geno_count(pGeno, fN, *pS, *pN);
			pS++; pN++;
		}
		for (; nSNP < fM; nSNP++)
			*pS++ = *pN++ = 0;
	}

	/// divide PCA_GenoSum by PCA_GenoNum (PCA_GenoSum / PCA_GenoNum)
	void DivideGeno()
	{
		double *p = tmp_var.Get();
		C_Int32 *pS = PCA_GenoSum.Get();
		C_Int32 *pN = PCA_GenoNum.Get();
		size_t n = fM;

	#if defined(COREARRAY_SIMD_AVX)
		const __m256d zero = _mm256_setzero_pd();
		for (; n >= 4; n-=4)
		{
			__m256d SD = _mm256_cvtepi32_pd(_mm_load_si128((__m128i const*)pS));
			pS += 4;
			__m256d ND = _mm256_cvtepi32_pd(_mm_load_si128((__m128i const*)pN));
			pN += 4;
			_mm256_store_pd(p, _mm256_and_pd(_mm256_div_pd(SD, ND),
				_mm256_cmp_pd(zero, ND, _CMP_LT_OQ)));
			p += 4;
		}
	#elif defined(COREARRAY_SIMD_SSE2)
		const __m128d zero = _mm_setzero_pd();
		for (; n >= 4; n-=4)
		{
			__m128i S = _mm_load_si128((__m128i const*)pS);
			pS += 4;
			__m128i N = _mm_load_si128((__m128i const*)pN);
			pN += 4;

			__m128d SD = _mm_cvtepi32_pd(S);
			__m128d ND = _mm_cvtepi32_pd(N);
			_mm_store_pd(p, _mm_and_pd(_mm_div_pd(SD, ND), _mm_cmplt_pd(zero, ND)));
			p += 2;

			SD = _mm_cvtepi32_pd(_mm_shuffle_epi32(S, 0x0E));
			ND = _mm_cvtepi32_pd(_mm_shuffle_epi32(N, 0x0E));
			_mm_store_pd(p, _mm_and_pd(_mm_div_pd(SD, ND), _mm_cmplt_pd(zero, ND)));
			p += 2;
		}
	#endif
		for (; n > 0; n--)
		{
			*p++ = (*pN > 0) ? ((double)*pS / *pN) : 0;
			pS ++; pN ++;
		}
	}

	/// computing scalar vector
	void rsqrt_prod()
	{
		double *p = tmp_var.Get();
		size_t n = fM;

	#if defined(COREARRAY_SIMD_AVX)
		const __m256d half = _mm256_set1_pd(0.5);
		const __m256d zero = _mm256_set1_pd(0.0);
		const __m256d one  = _mm256_set1_pd(1.0);
		for (; n >= 4; n-=4, p+=4)
		{
			__m256d s = _mm256_mul_pd(_mm256_load_pd(p), half);
			__m256d m = _mm256_and_pd(_mm256_cmp_pd(zero, s, _CMP_LT_OQ),
				_mm256_cmp_pd(s, one, _CMP_LT_OQ));
			s = _mm256_mul_pd(s, _mm256_sub_pd(one, s));
			s = _mm256_div_pd(one, _mm256_sqrt_pd(s));
			_mm256_store_pd(p, _mm256_and_pd(s, m));
		}
	#elif defined(COREARRAY_SIMD_SSE2)
		const __m128d half = _mm_set1_pd(0.5);
		const __m128d zero = _mm_set1_pd(0.0);
		const __m128d one  = _mm_set1_pd(1.0);
		for (; n >= 2; n-=2, p+=2)
		{
			__m128d s = _mm_mul_pd(_mm_load_pd(p), half);
			__m128d m = _mm_and_pd(_mm_cmplt_pd(zero, s), _mm_cmplt_pd(s, one));
			s = _mm_mul_pd(s, _mm_sub_pd(one, s));
			s = _mm_div_pd(one, _mm_sqrt_pd(s));
			_mm_store_pd(p, _mm_and_pd(s, m));
		}
	#endif
		for (; n > 0; n--)
		{
			double s = (*p) * 0.5;
			*p++ = (0<s && s<1) ? (1.0 / sqrt(s*(1-s))) : 0;
		}
	}

	inline size_t N() const { return fN; }
	inline size_t M() const { return fM; }

protected:
	size_t fN, fM;

	void _Reset()
	{
		PCA_GenoSum.Reset(fM);
		PCA_GenoNum.Reset(fM);
		tmp_var.Reset(fM);
	}

	void _Clear()
	{
		PCA_GenoSum.Clear();
		PCA_GenoNum.Clear();
		tmp_var.Clear();
	}
};


// Mean-adjusted genotype matrix (AVX, N: # of samples, M: # of SNPs)
class COREARRAY_DLL_LOCAL CPCAMat_Alg1: public CPCAMat_Base
{
public:
	CPCAMat_Alg1() { fN = fM = 0; }
	CPCAMat_Alg1(size_t n, size_t m) { Reset(n, m); }

	void Reset(size_t n, size_t m)
	{
	#if defined(COREARRAY_SIMD_AVX)
		// 4-aligned array
		if (m & 0x03) m += 4 - (m & 0x03);
	#elif defined(COREARRAY_SIMD_SSE2)
		// 2-aligned array
		if (m & 0x01) m += 2 - (m & 0x01);
	#endif
		fGenotype.Reset(n * m);
		fN = n; fM = m;
		_Reset();
	}

	void Clear()
	{
		_Clear();
		fGenotype.Clear();
	}

	/// detect the effective value for BlockNumSNP
	static void PCA_Detect_BlockNumSNP(int nSamp)
	{
		size_t Cache = GetOptimzedCache();
		BlockNumSNP = Cache / (sizeof(double)*nSamp);
		BlockNumSNP = (BlockNumSNP / 4) * 4;
		if (BlockNumSNP < 16) BlockNumSNP = 16;
	}

	// time-consuming function
	void COREARRAY_CALL_ALIGN MulAdd(IdMatTri &Idx, size_t IdxCnt, double *pOut)
	{
		for (; IdxCnt > 0; IdxCnt--, ++Idx)
		{
			double *p1 = base() + Idx.Row() * fM;
			double *p2 = base() + Idx.Column() * fM;
			size_t n = fM;

		#if defined(COREARRAY_SIMD_AVX)

			// fM can be divided by 4
			__m256d rv4_1 = _mm256_setzero_pd();
			__m256d rv4_2 = _mm256_setzero_pd();
			for (; n >= 16; n -= 16)
			{
				rv4_1 = _mm256_add_pd(rv4_1, _mm256_mul_pd(
					_mm256_load_pd(p1), _mm256_load_pd(p2)));
				rv4_2 = _mm256_add_pd(rv4_2, _mm256_mul_pd(
					_mm256_load_pd(p1 + 4), _mm256_load_pd(p2 + 4)));
				p1 += 8; p2 += 8;
				rv4_1 = _mm256_add_pd(rv4_1, _mm256_mul_pd(
					_mm256_load_pd(p1), _mm256_load_pd(p2)));
				rv4_2 = _mm256_add_pd(rv4_2, _mm256_mul_pd(
					_mm256_load_pd(p1 + 4), _mm256_load_pd(p2 + 4)));
				p1 += 8; p2 += 8;
			}

			__m256d rv4 = _mm256_add_pd(rv4_1, rv4_2);
			for (; n >= 4; n -= 4)
			{
				rv4 = _mm256_add_pd(rv4, _mm256_mul_pd(
					_mm256_load_pd(p1), _mm256_load_pd(p2)));
				p1 += 4; p2 += 4;
			}

			(*pOut++) += vec_avx_sum_f64(rv4);

		#elif defined(COREARRAY_SIMD_SSE2)

			// unroll loop
			__m128d rv2_1 = _mm_setzero_pd();
			__m128d rv2_2 = _mm_setzero_pd();
			for (; n >= 8; n -= 8)
			{
				rv2_1 = _mm_add_pd(rv2_1, _mm_mul_pd(_mm_load_pd(p1),
					_mm_load_pd(p2)));
				rv2_2 = _mm_add_pd(rv2_2, _mm_mul_pd(_mm_load_pd(p1+2),
					_mm_load_pd(p2+2)));
				p1 += 4; p2 += 4;
				rv2_1 = _mm_add_pd(rv2_1, _mm_mul_pd(_mm_load_pd(p1),
					_mm_load_pd(p2)));
				rv2_2 = _mm_add_pd(rv2_2, _mm_mul_pd(_mm_load_pd(p1+2),
					_mm_load_pd(p2+2)));
				p1 += 4; p2 += 4;
			}

			__m128d rv2 = _mm_add_pd(rv2_1, rv2_2);
			for (; n >= 2; n -= 2)
			{
				rv2 = _mm_add_pd(rv2, _mm_mul_pd(_mm_load_pd(p1),
					_mm_load_pd(p2)));
				p1 += 2; p2 += 2;
			}

			(*pOut++) += vec_sum_f64(rv2);

		#else

			double rv = 0;
			// unroll loop
			for (; n >= 4; n -= 4)
			{
				rv += p1[0] * p2[0]; rv += p1[1] * p2[1];
				rv += p1[2] * p2[2]; rv += p1[3] * p2[3];
				p1 += 4; p2 += 4;
			}
			for (; n > 0; n--)
				rv += (*p1++) * (*p2++);
			(*pOut++) += rv;

		#endif
		}
	}

	// time-consuming function
	void COREARRAY_CALL_ALIGN MulAdd2(IdMatTri &Idx, size_t IdxCnt,
		size_t Length, double *pOut)
	{
		for (; IdxCnt > 0; IdxCnt--, ++Idx)
		{
			double *p1 = base() + Idx.Row() * fM;
			double *p2 = base() + Idx.Column() * fM;
			size_t n = Length;

		#if defined(COREARRAY_SIMD_AVX)

			// unroll loop
			__m256d rv4_1 = _mm256_setzero_pd();
			__m256d rv4_2 = _mm256_setzero_pd();
			for (; n >= 8; n -= 8)
			{
				rv4_1 = _mm256_add_pd(rv4_1, _mm256_mul_pd(
					_mm256_load_pd(p1), _mm256_load_pd(p2)));
				rv4_2 = _mm256_add_pd(rv4_2, _mm256_mul_pd(
					_mm256_load_pd(p1 + 4), _mm256_load_pd(p2 + 4)));
				p1 += 8; p2 += 8;
			}
			__m256d rv4 = _mm256_add_pd(rv4_1, rv4_2);
			if (n >= 4)
			{
				rv4 = _mm256_add_pd(rv4, _mm256_mul_pd(
					_mm256_load_pd(p1), _mm256_load_pd(p2)));
				p1 += 4; p2 += 4; n -= 4;
			}

			double x[4] __attribute__((aligned(32)));
			_mm256_store_pd(x, rv4);
			double rv = x[0] + x[1] + x[2] + x[3];

			switch (n)
			{
				case 3:
					rv += p1[0]*p2[0] + p1[1]*p2[1] + p1[2]*p2[2]; break;
				case 2:
					rv += p1[0]*p2[0] + p1[1]*p2[1]; break;
				case 1:
					rv += p1[0]*p2[0]; break;
			}

			(*pOut++) += rv;

		#elif defined(COREARRAY_SIMD_SSE2)

			// unroll loop
			__m128d rv2 = _mm_setzero_pd();
			for (; n >= 8; n -= 8)
			{
				rv2 = _mm_add_pd(rv2, _mm_mul_pd(_mm_load_pd(&p1[0]),
					_mm_load_pd(&p2[0])));
				rv2 = _mm_add_pd(rv2, _mm_mul_pd(_mm_load_pd(&p1[2]),
					_mm_load_pd(&p2[2])));
				rv2 = _mm_add_pd(rv2, _mm_mul_pd(_mm_load_pd(&p1[4]),
					_mm_load_pd(&p2[4])));
				rv2 = _mm_add_pd(rv2, _mm_mul_pd(_mm_load_pd(&p1[6]),
					_mm_load_pd(&p2[6])));
				p1 += 8; p2 += 8;
			}
			for (; n >= 2; n -= 2)
			{
				rv2 = _mm_add_pd(rv2, _mm_mul_pd(_mm_load_pd(p1),
					_mm_load_pd(p2)));
				p1 += 2; p2 += 2;
			}

			double x[2] __attribute__((aligned(16)));
			_mm_store_pd(x, rv2);
			(*pOut++) += (n <= 0) ? (x[0] + x[1]) : (x[0] + x[1] + (*p1) * (*p2));

		#else

			double rv = 0;
			// unroll loop
			for (; n >= 4; n -= 4)
			{
				rv += p1[0] * p2[0]; rv += p1[1] * p2[1];
				rv += p1[2] * p2[2]; rv += p1[3] * p2[3];
				p1 += 4; p2 += 4;
			}
			for (; n > 0; n--)
				rv += (*p1++) * (*p2++);
			(*pOut++) += rv;

		#endif
		}
	}


	inline double *base() { return fGenotype.Get(); }

	/// mean-adjusted genotypes (fGenotype - tmp_var)
	void COREARRAY_CALL_ALIGN GenoSub()
	{
		double *pGeno = fGenotype.Get();
		for (size_t num=fN; num > 0; num--, pGeno+=fM)
		{
			size_t n = fM;
			double *p = pGeno, *s = tmp_var.Get();
		#if defined(COREARRAY_SIMD_AVX)
			for (; n >= 4; n -= 4)
			{
				__m256d v = _mm256_sub_pd(_mm256_load_pd(p), _mm256_load_pd(s));
				_mm256_store_pd(p, v);
				p += 4; s += 4;
			}
		#endif
		#if defined(COREARRAY_SIMD_SSE2)
			for (; n >= 2; n -= 2)
			{
				__m128d v = _mm_sub_pd(_mm_load_pd(p), _mm_load_pd(s));
				_mm_store_pd(p, v);
				p += 2; s += 2;
			}
		#endif
			for (; n > 0; n--) *p++ -= (*s++);
		}
	}

	/// variance-adjusted genotypes (fGenotype * tmp_var)
	void COREARRAY_CALL_ALIGN GenoMul()
	{
		double *pGeno = fGenotype.Get();
		for (size_t num=fN; num > 0; num--, pGeno+=fM)
		{
			size_t n = fM;
			double *p = pGeno, *s = tmp_var.Get();
		#if defined(COREARRAY_SIMD_AVX)
			for (; n >= 4; n -= 4)
			{
				__m256d v = _mm256_mul_pd(_mm256_load_pd(p), _mm256_load_pd(s));
				_mm256_store_pd(p, v);
				p += 4; s += 4;
			}
		#endif
		#if defined(COREARRAY_SIMD_SSE2)
			for (; n >= 2; n -= 2)
			{
				__m128d v = _mm_mul_pd(_mm_load_pd(p), _mm_load_pd(s));
				_mm_store_pd(p, v);
				p += 2; s += 2;
			}
		#endif
			for (; n > 0; n--) *p++ *= (*s++);
		}
	}

private:
	VEC_AUTO_PTR<double> fGenotype;
};


/// The mean-adjusted genotype buffers
static CPCAMat_Alg1 PCA_Mat1;


// ---------------------------------------------------------------------
// Vectorization Computing, algorithm 2 (bit manipulation)

// Mean-adjusted genotype matrix (AVX, N: # of samples, M: # of SNPs)
class COREARRAY_DLL_LOCAL CPCAMat_Alg2: public CPCAMat_Base
{
public:
	CPCAMat_Alg2() { fN = fM = 0; }
	CPCAMat_Alg2(size_t n, size_t m) { Reset(n, m); }

	void Reset(size_t n, size_t m)
	{
		fN = n; fM = m;
		fGenotype.Reset(n * m / 2);
		fLookup.Reset(m * 256);
		_Reset();
	}

	/// detect the effective value for BlockNumSNP
	static void PCA_Detect_BlockNumSNP(int nSamp)
	{
		C_UInt64 L1Cache = GDS_Mach_GetCPULevelCache(1);
		if (L1Cache <= 0) L1Cache = 32*1024;
		C_UInt64 L2Cache = GDS_Mach_GetCPULevelCache(2);
		C_UInt64 L3Cache = GDS_Mach_GetCPULevelCache(3);
		C_UInt64 Cache = (L2Cache > L3Cache) ? L2Cache : L3Cache;
		if ((C_Int64)Cache <= 0) Cache = 1024*1024; // 1M
		BlockNumSNP = (Cache - 4*L1Cache) / nSamp;
		BlockNumSNP = (BlockNumSNP / 4) * 4;
		if (BlockNumSNP > 256) BlockNumSNP = 256;
		if (BlockNumSNP < 16) BlockNumSNP = 16;
	}

	void Clear()
	{
		_Clear();
		fGenotype.Clear();
		fLookup.Clear();
	}

	// time-consuming function
	void MulAdd(IdMatTri &Idx, size_t IdxCnt, size_t Length, double *pOut)
	{
		if (Length & 0x01) Length++;
		size_t M = fM >> 1;

		for (; IdxCnt > 0; IdxCnt--, ++Idx)
		{
			C_UInt8 *p1 = base() + Idx.Row() * M;
			C_UInt8 *p2 = base() + Idx.Column() * M;
			double *p = lookup();
			double sum = 0;

			for (size_t n=Length; n > 0; n -= 2)
			{
				sum += p[(*p1++) | ((*p2++) << 2)];
				p += 256;
			}

			(*pOut++) += sum;
		}
	}

	inline C_UInt8 *base() { return fGenotype.Get(); }
	inline double *lookup() { return fLookup.Get(); }

private:
	VEC_AUTO_PTR<C_UInt8> fGenotype;
	VEC_AUTO_PTR<double> fLookup;
};


/// The mean-adjusted genotype buffers
static CPCAMat_Alg2 PCA_Mat2;

static void _Fill_Lookup_16(double v[], double avg, double scale)
{
	//             0x03        0x0C
	v[0]        = (0 - avg) * (0 - avg) * scale;  // 0, 0
	v[1] = v[4] = (1 - avg) * (0 - avg) * scale;  // 1, 0
	v[2] = v[8] = (2 - avg) * (0 - avg) * scale;  // 2, 0
	v[5]        = (1 - avg) * (1 - avg) * scale;  // 1, 1
	v[6] = v[9] = (2 - avg) * (1 - avg) * scale;  // 2, 1
	v[10]       = (2 - avg) * (2 - avg) * scale;  // 2, 2
	v[3] = v[7] = v[11] = v[12] = v[13] = v[14] = v[15] = 0;
}

/// Convert the raw genotypes to the mean-adjusted genotypes (GenoBuf: sample x SNP)
static void _Do_PCA_Read_SampxSNP_Alg2(C_UInt8 *pGeno, long Start,
	long SNP_Cnt, void *Param)
{
	const size_t n = MCWorkingGeno.Space().SampleNum();

	// get the genotype sum and number
	PCA_Mat2.SummarizeGeno_SampxSNP(pGeno, SNP_Cnt);
	// get 2 * \bar{p}_l, saved in PCA_Mat2.tmp_var
	PCA_Mat2.DivideGeno();

	// transpose genotypes in PCA_Mat2, pack genotypes 0x0F and 0xF0
	C_UInt8 *pG = PCA_Mat2.base();
	for (size_t iSamp=0; iSamp < n; iSamp++)
	{
		C_UInt8 *pp = pG; pG += PCA_Mat2.M() >> 1;
		size_t m = SNP_Cnt, i = 0;
		for (; m >= 2; m-=2)
		{
			C_UInt8 g1 = pGeno[n*i     + iSamp];
			C_UInt8 g2 = pGeno[n*(i+1) + iSamp];
			i += 2;
			*pp ++ = ((g1 <= 2) ? g1 : 3) | (((g2 <= 2) ? g2 : 3) << 4);
		}
		if (m > 0)
		{
			C_UInt8 g = pGeno[n*i + iSamp];
			*pp ++ = ((g <= 2) ? g : 3) | 0x30;  // 0x30 is missing
		}
	}

	// set lookup table
	double v1[16], v2[16];
	double *p = PCA_Mat2.tmp_var.Get(); // 2\bar{p}_l
	size_t m = SNP_Cnt;
	if (m & 0x01)
	{
		m ++;
		PCA_Mat2.PCA_GenoSum[m] = PCA_Mat2.PCA_GenoNum[m] = 0;
		PCA_Mat2.tmp_var[m] = 0;
	}

	if (BayesianNormal)
	{
		C_Int32 *pSum = PCA_Mat2.PCA_GenoSum.Get();
		C_Int32 *pNum = PCA_Mat2.PCA_GenoNum.Get();
		double *pp = PCA_Mat2.lookup(), s;
		for (; m > 0; m-=2)
		{
			s = ((*pSum++) + 1.0) / (2 * (*pNum++) + 2);
			_Fill_Lookup_16(v1, *p++, 1.0 / (s * (1 - s)));
			s = ((*pSum++) + 1.0) / (2 * (*pNum++) + 2);
			_Fill_Lookup_16(v2, *p++, 1.0 / (s * (1 - s)));
			for (size_t i=0; i < 16; i++)
				for (size_t j=0; j < 16; j++)
					*pp++ = v1[j] + v2[i];
		}
	} else {
		double *pp = PCA_Mat2.lookup(), f, a;
		for (; m > 0; m-=2)
		{
			f = *p++; a = f * 0.5;
			_Fill_Lookup_16(v1, f, (0<a && a<1) ? (1.0 / (a*(1-a))) : 0);
			f = *p++; a = f * 0.5;
			_Fill_Lookup_16(v2, f, (0<a && a<1) ? (1.0 / (a*(1-a))) : 0);
			for (size_t i=0; i < 16; i++)
				for (size_t j=0; j < 16; j++)
					*pp++ = v1[j] + v2[i];
		}
	}
}

/// Compute the covariance matrix
static void _Do_PCA_ComputeCov_Alg2(int IdxThread, long Start,
	long SNP_Cnt, void* Param)
{
	double *base = (double*)Param;
	IdMatTri I = Array_Thread_MatIdx[IdxThread];
	PCA_Mat2.MulAdd(I, Array_Thread_MatCnt[IdxThread], SNP_Cnt,
		base + I.Offset());
}



// ================== PCA covariance matrix ==================

// -----------------------------------------------------------
// Exact PCA

class COREARRAY_DLL_LOCAL CExactPCA: protected CPCAMat_Alg1
{
private:
	CdBaseWorkSpace &Space;
	double *ptrCov;

	void thread_cov_outer(size_t i, size_t n)
	{
		IdMatTri I = Array_Thread_MatIdx[i];
		MulAdd(I, Array_Thread_MatCnt[i], ptrCov + I.Offset());
	}

public:
	CExactPCA(CdBaseWorkSpace &space): Space(space) { }

	/// run the algorithm
	void Run(CdMatTri<double> &Cov, int NumThread, bool Bayesian,
		bool verbose)
	{
		if (NumThread < 1) NumThread = 1;
		const size_t nSamp = Space.SampleNum();

		// detect the appropriate block size
		PCA_Detect_BlockNumSNP(nSamp);
		if (verbose)
		{
			Rprintf("%s    (internal increment: %d)\n", TimeToStr(),
				(int)BlockNumSNP);
		}

		// initialize
		Reset(nSamp, BlockNumSNP);
		ptrCov = Cov.Get();
		memset(ptrCov, 0, sizeof(double)*Cov.Size());

		// thread pool
		CThreadPoolEx<CExactPCA> thpool(NumThread);
		Array_SplitJobs(NumThread, nSamp, Array_Thread_MatIdx,
			Array_Thread_MatCnt);

		// genotypes (0, 1, 2 and NA)
		VEC_AUTO_PTR<C_UInt8> Geno(nSamp * BlockNumSNP);
		C_UInt8 *pGeno = Geno.Get();

		// genotype buffer, false for no memory buffer
		CGenoReadBySNP WS(NumThread, Space, BlockNumSNP, verbose ? -1 : 0, false);

		// for-loop
		WS.Init();
		while (WS.Read(Geno.Get()))
		{
			// get the genotype sum and number
			SummarizeGeno_SampxSNP(pGeno, WS.Count());
			// get 2 * \bar{p}_l, saved in PCA_Mat.tmp_var
			DivideGeno();

			// transpose genotypes in PCA_Mat, set missing values to the average
			double *p = base();
			for (size_t i=0; i < nSamp; i++)
			{
				double *pp = p; p += M();
				size_t m = WS.Count();
				for (size_t j=0; j < m; j++)
				{
					C_UInt8 g = pGeno[nSamp*j + i];
					*pp ++ = (g <= 2) ? g : tmp_var[j];
				}
				for (; m < M(); m++)
					*pp ++ = 0;
			}

			// G@ij - 2\bar{p}_l
			GenoSub();

			// 1 / sqrt(p@j*(1-p@j))
			if (!Bayesian)
			{
				rsqrt_prod();
			} else {
				p = tmp_var.Get();
				C_Int32 *pSum = PCA_GenoSum.Get();
				C_Int32 *pNum = PCA_GenoNum.Get();
				for (size_t m=WS.Count(); m > 0; m--)
				{
					double s = ((*pSum++) + 1.0) / (2 * (*pNum++) + 2);
					*p++ = 1.0 / sqrt(s * (1 - s));
				}
			}

			// (G@ij - 2p@j) / sqrt(p@j*(1-p@j)), get the normalized genotypes
			GenoMul();

			// outer product, using thread thpool
			thpool.BatchWork(this, &CExactPCA::thread_cov_outer, NumThread);

			// update
			WS.ProgressForward(WS.Count());
		}
	}
};


    /// Calculate the genetic covariance
	void DoCovCalculate_Alg2(CdMatTri<double> &PublicCov, int NumThread,
		const char *Info, bool verbose)
	{
		// initialize
		PCA_Mat2.Reset(PublicCov.N(), BlockNumSNP);
		memset(PublicCov.Get(), 0, sizeof(double)*PublicCov.Size());

		MCWorkingGeno.Progress.Info = Info;
		MCWorkingGeno.Progress.Show() = verbose;
		Array_SplitJobs(NumThread, PublicCov.N(), Array_Thread_MatIdx,
			Array_Thread_MatCnt);

		MCWorkingGeno.InitParam(true, RDim_Sample_X_SNP, BlockNumSNP);
		MCWorkingGeno.Run(NumThread, &_Do_PCA_Read_SampxSNP_Alg2,
			&_Do_PCA_ComputeCov_Alg2, PublicCov.Get());

		PCA_Mat2.Clear();
	}


// -----------------------------------------------------------
// Fast randomized PCA

class COREARRAY_DLL_LOCAL CRandomPCA
{
private:
	CdBaseWorkSpace &Space;

	size_t nSamp;  /// the number of selected samples
	size_t nSNP;   /// the number of selected SNPs
	double *AuxMat;  /// auxiliary matrix (G_i = nSamp X AuxDim)
	size_t AuxDim;  /// the number of columns
	int IterNum;  /// the number of iterations

	size_t hsize;  /// the number of columns
	VEC_AUTO_PTR<double> MatH;
	VEC_AUTO_PTR<double> LookupY;  /// normalized genotypes, Y

	VEC_AUTO_PTR<C_UInt8> Geno;  /// the genotype buffer (0, 1, 2 and NA)
	VEC_AUTO_PTR<double> Y_mc;  /// normalized genotypes at a site, with multiple cores
	VEC_AUTO_PTR<double> AuxMat_mc; /// Aux matrices, with multiple cores
	VEC_AUTO_PTR<double> MatT;  /// T = U_H^T * Y

	size_t iSNP;  /// the starting SNP index
	int iteration;  /// the iteration

	// auxiliary thread variables
	vector<size_t> thread_start, thread_length;

	void thread_lookup_y(size_t i, size_t n)
	{
		C_UInt8 *g = &Geno[i * nSamp];
		double *Y = &LookupY[(iSNP + i) * 4];
		for (; n > 0; n--)
		{
			C_Int32 sum, num;
			vec_u8_geno_count(g, nSamp, sum, num);
			g += nSamp;

			double avg = (num>0) ? (double)sum / num : 0;
			double p = avg * 0.5;
			double s = (0<p && p<1) ? 1.0 / sqrt(2*p*(1-p)) : 0;

			Y[0] = (0 - avg) * s; Y[1] = (1 - avg) * s;
			Y[2] = (2 - avg) * s; Y[3] = 0;
			Y += 4;
		}
	}

	void thread_Y_x_G_i(size_t i, size_t num)
	{
		C_UInt8 *pGeno = &Geno[i * nSamp];
		i += iSNP;
		for (; num > 0; num--, i++)
		{
			double *pH = &MatH[(AuxDim * iteration) + (i * hsize)];
			double *pY = &LookupY[i * 4];
			double *pA = AuxMat;

			for (size_t j=0; j < AuxDim; j++)
			{
				C_UInt8 *pG = pGeno;
				size_t n = nSamp;
				double sum = 0;

			#if defined(COREARRAY_SIMD_AVX)

				__m256d sum4 = _mm256_setzero_pd();
				for (; n >= 4; n-=4)
				{
					__m256d a = _mm256_loadu_pd(pA); pA += 4;

				#if defined(COREARRAY_SIMD_AVX2)
					__m256i k4 = _mm256_set_epi64x(pG[3], pG[2], pG[1], pG[0]);
					__m256d y = _mm256_i64gather_pd(pY, k4, sizeof(double));
				#else
					__m256d y = _mm256_set_pd(pY[pG[3]], pY[pG[2]],
						pY[pG[1]], pY[pG[0]]);
				#endif

					pG += 4;
					sum4 = _mm256_add_pd(sum4, _mm256_mul_pd(a, y));
				}
				sum = vec_avx_sum_f64(sum4);

			#elif defined(COREARRAY_SIMD_SSE2)

				__m128d sum2 = _mm_setzero_pd();
				for (; n >= 2; n-=2)
				{
					__m128d a = _mm_loadu_pd(pA); pA += 2;
					__m128d y = _mm_set_pd(pY[pG[1]], pY[pG[0]]); pG += 2;
					sum2 = _mm_add_pd(sum2, _mm_mul_pd(a, y));
				}
				sum = vec_sum_f64(sum2);

			#endif

				for (; n > 0; n--)
					sum += pY[*pG++] * (*pA++);
				*pH++ = sum;
			}

			pGeno += nSamp;
		}
	}

	void thread_YT_x_H_i(size_t i, size_t num)
	{
		size_t ii = thread_start[i];
		size_t n  = thread_length[i];
		double *pH = &MatH[AuxDim*iteration + (iSNP + ii)*hsize];
		double *pY = &LookupY[(iSNP + ii) * 4];

		for (; n > 0; n--, ii++)
		{
			// set Y_mc
			C_UInt8 *pG = &Geno[ii * nSamp];
			double *Y = &Y_mc[i * nSamp];
			for (size_t j=0; j < nSamp; j++)
			{
				Y[j] = pY[(*pG < 4) ? *pG : 3];
				pG ++;
			}

			double *pA = (i == 0) ? AuxMat : &AuxMat_mc[(i-1)*nSamp*AuxDim];
			for (size_t j=0; j < AuxDim; j++)
			{
				// pA += Y * pH[j]
				pA = vec_f64_addmul(pA, Y, nSamp, pH[j]);
			}

			pY += 4;
			pH += hsize;
		}
	}

	void thread_U_H_x_Y(size_t i, size_t num)
	{
		size_t ii = thread_start[i];
		size_t n  = thread_length[i];
		double *pH = &MatH[(iSNP + ii) * hsize];
		double *pY = &LookupY[(iSNP + ii) * 4];

		for (; n > 0; n--, ii++)
		{
			C_UInt8 *pG = &Geno[ii * nSamp];
			double *pT = &MatT[i * hsize * nSamp];

			for (size_t j=0; j < nSamp; j++)
			{
				double y = pY[(*pG < 4) ? *pG : 3];
				pG ++;
				pT = vec_f64_addmul(pT, pH, hsize, y);
			}

			pY += 4;
			pH += hsize;
		}
	}

	static void svd_vt(double *a, int m, int n, double *d)
	{
		int info = 0;
		double u=0, vt=0, w=0;
		vector<double> ds;
		if (!d)
			{ ds.resize(min(m, n)); d = &ds[0]; }

		int lwork = -1;
		F77_NAME(dgesvd)("N", "O", &m, &n, a, &m, d, &u, &m, &vt, &n,
			&w, &lwork, &info);
		if (info != 0)
			throw ErrCoreArray("LAPACK::DGESVD error (%d).", info);

		lwork = (int)w;
		vector<double> work(lwork);
		F77_NAME(dgesvd)("N", "O", &m, &n, a, &m, d, &u, &m, &vt, &n,
			&work[0], &lwork, &info);
		if (info != 0)
			throw ErrCoreArray("LAPACK::DGESVD error (%d).", info);
	}

public:
	/// constructor
	CRandomPCA(CdBaseWorkSpace &space, double *aux_mat, size_t aux_dim,
		int iter_num): Space(space)
	{
		nSamp = Space.SampleNum();
		nSNP  = Space.SNPNum();
		AuxMat = aux_mat;
		AuxDim = aux_dim;
		IterNum = iter_num;
		hsize = AuxDim * (IterNum + 1);
		MatH.Reset(nSNP * hsize);
		LookupY.Reset(nSNP * 4);
		iSNP = 0;
		iteration = 0;
	}

	/// run the algorithm
	SEXP Run(int NumThread, bool verbose)
	{
		if (NumThread < 1) NumThread = 1;
		size_t IncSNP = (256 / NumThread) * NumThread;
		if (IncSNP < 16) IncSNP = 16;
		if (verbose)
			Rprintf("%s\n", TimeToStr());

		Geno.Reset(nSamp * IncSNP);
		Y_mc.Reset(nSamp * NumThread);
		AuxMat_mc.Reset(nSamp * AuxDim * (NumThread - 1));
		thread_start.resize(NumThread);
		thread_length.resize(NumThread);

		// thread thpool
		CThreadPoolEx<CRandomPCA> thpool(NumThread);

		// genotype buffer, false for no memory buffer
		CGenoReadBySNP WS(NumThread, Space, IncSNP,
			verbose ? C_Int64(nSNP)*(2*IterNum + 1) : 0, false);

		// for-loop
		for (iteration=0; iteration <= IterNum; iteration++)
		{
			WS.Init();
			while (WS.Read(Geno.Get(), iSNP))
			{
				// need to initialize the variable LookupY
				if (iteration == 0)
					thpool.BatchWork(this, &CRandomPCA::thread_lookup_y, WS.Count());

				// update H_i = Y * G_i (G_i stored in AuxMat)
				thpool.BatchWork(this, &CRandomPCA::thread_Y_x_G_i, WS.Count());

				// update
				WS.ProgressForward(WS.Count());
			}

			// update G_{i+1} = Y^T * H_i / nSNP
			if (iteration < IterNum)
			{
				memset(AuxMat, 0, sizeof(double)*AuxDim*nSamp);
				memset(AuxMat_mc.Get(), 0, sizeof(double)*AuxMat_mc.Length());

				WS.Init();
				while (WS.Read(Geno.Get(), iSNP))
				{
					// update G_{i+1} = Y^T * H_i
					thpool.Split(NumThread, WS.Count(), &thread_start[0], &thread_length[0]);
					thpool.BatchWork(this, &CRandomPCA::thread_YT_x_H_i, NumThread);

					if (NumThread > 1)
					{
						// update G_{i+1} from AuxMat_mc
						size_t n = nSamp * AuxDim;
						for (size_t i=0; i < (size_t)NumThread-1; i++)
							vec_f64_add(AuxMat, &AuxMat_mc[i * n], n);
					}

					// update
					WS.ProgressForward(WS.Count());
				}

				// divide AuxMat by nSNP
				vec_f64_mul(AuxMat, AuxDim*nSamp, 1.0/nSNP);
			}
		}

		if (verbose)
			Rprintf("%s    Begin projecting genotypes and SVD\n", TimeToStr());

		// SVD MatH, get vt stored in MatH
		svd_vt(&MatH[0], hsize, nSNP, NULL);

		// T = U_H^T * Y
		MatT.Reset(hsize*nSamp*NumThread);
		memset(MatT.Get(), 0, sizeof(double)*MatT.Length());

		// for-loop
		WS.Init();
		while (WS.Read(Geno.Get(), iSNP))
		{
//			thpool.Split(1, inc_snp, &thread_start[0], &thread_length[0]);
//			CRandomPCA::thread_U_H_x_Y(0, 1);

			// update T = U_H^T * Y
			thpool.Split(1, WS.Count(), &thread_start[0], &thread_length[0]);
			thpool.BatchWork(this, &CRandomPCA::thread_U_H_x_Y, 1);

/*			if (NumThread > 1)
			{
				// update T = U_H^T * Y from MatT
				size_t n = hsize*nSamp;
				for (size_t i=1; i < (size_t)NumThread; i++)
					vec_f64_add(&MatT[0], &MatT[i * n], n);
			}
*/
		}

		vector<double> sigma(nSamp);
		svd_vt(&MatT[0], hsize, nSamp, &sigma[0]);

		SEXP rv_ans = PROTECT(NEW_LIST(2));
		{
			SEXP d = NEW_NUMERIC(nSamp);
			memcpy(REAL(d), &sigma[0], sizeof(double) * nSamp);
			SET_ELEMENT(rv_ans, 0, d);

			SEXP h = Rf_allocMatrix(REALSXP, hsize, nSamp);
			memcpy(REAL(h), &MatT[0], sizeof(double) * nSamp * hsize);
			SET_ELEMENT(rv_ans, 1, h);
		}

		UNPROTECT(1);

		if (verbose)
			Rprintf("%s    End\n", TimeToStr());

		return rv_ans;
	}
};



// ==================== SNP Correlations =====================

class COREARRAY_DLL_LOCAL CPCA_SNPCorr
{
private:
	CdBaseWorkSpace &Space;  ///< working genotypes
	VEC_AUTO_PTR<C_UInt8> Geno;  ///< genotypes (0, 1, 2 and NA)
	size_t nSamp;      ///< the number of samples
	size_t NumEigVal;  ///< the number of eigenvalues
	double *pEigVect;  ///< the pointer to eigenvectors
	double *pCorr;     ///< the pointer to correlation

	// Correlation
	static double SNP_PC_Corr(double *pX, C_UInt8 *pY, size_t n)
	{
		size_t m=0;
		double XY=0, X=0, XX=0, Y=0, YY=0, ans=R_NaN;
		while (n > 0)
		{
			if (*pY < 3)
			{
				XY += (*pX) * (*pY);
				X += *pX; XX += (*pX) * (*pX);
				Y += *pY; YY += (*pY) * (*pY);
				m ++;
			}
			pX++; pY++; n--;
		}
		if (m > 1)
		{
			double c1 = XX - X*X/m, c2 = YY - Y*Y/m, val = c1*c2;
			if (val > 0)
				ans = (XY - X*Y/m) / sqrt(val);
		}
		return ans;
	}

	void thread_corr(size_t i, size_t num)
	{
		C_UInt8 *pGeno = Geno.Get() + nSamp * i;
		double *pOut = pCorr + NumEigVal * i;
		for (; num > 0; num--, i++)
		{
			double *pEig = pEigVect;
			for (size_t j=0; j < NumEigVal; j++)
			{
				*pOut++ = SNP_PC_Corr(pEig, pGeno, nSamp);
				pEig += nSamp;
			}
			pGeno += nSamp;
		}
	}

public:
	/// constructor
	CPCA_SNPCorr(CdBaseWorkSpace &space): Space(space) { }

	/// run the algorithm
	void Run(double OutSNPCorr[], size_t NumEig, double EigVect[],
		int NumThread, bool verbose)
	{
		if (NumThread < 1) NumThread = 1;
		nSamp = Space.SampleNum();
		NumEigVal = NumEig; pEigVect = EigVect;

		// detect the appropriate block size
		size_t Cache = GetOptimzedCache();
		size_t nBlock = Cache / nSamp;
		nBlock = (nBlock / 4) * 4;
		if (nBlock < 128) nBlock = 128;
		if (nBlock > 65536) nBlock = 65536;
		if (verbose)
			Rprintf("%s    (internal increment: %d)\n", TimeToStr(), (int)nBlock);

		// thread thpool
		CThreadPoolEx<CPCA_SNPCorr> thpool(NumThread);
		// genotypes (0, 1, 2 and NA)
		Geno.Reset(nSamp * nBlock);
		// genotype buffer, false for no memory buffer
		CGenoReadBySNP WS(NumThread, Space, nBlock, verbose ? -1 : 0, false);

		// for-loop
		WS.Init();
		while (WS.Read(Geno.Get()))
		{
			pCorr = OutSNPCorr + WS.Index() * NumEig;
			// using thread thpool
			thpool.BatchWork(this, &CPCA_SNPCorr::thread_corr, WS.Count());
			// update
			WS.ProgressForward(WS.Count());
		}
	}
};



// ================== SNP Loadings ==================

class COREARRAY_DLL_LOCAL CPCA_SNPLoad
{
private:
	CdBaseWorkSpace &Space;  ///< working genotypes
	VEC_AUTO_PTR<C_UInt8> Geno;  ///< genotypes (0, 1, 2 and NA)
	size_t nSamp;      ///< the number of samples
	size_t NumEigVal;  ///< the number of eigenvalues
	double *pEigVect;  ///< the pointer to eigenvectors
	double *pLoading;  ///< the pointer to correlation
	double *pAFreq;    ///< the pointer to allele frequency
	double *pScale;    ///< the pointer to scale factor

	void thread_loading(size_t i, size_t num)
	{
		C_UInt8 *pGeno = Geno.Get() + nSamp * i;
		double *pOut = pLoading + NumEigVal * i;
		for (; num > 0; num--, i++)
		{
			double avg, scale;
			C_Int32 gSum, gNum;
			vec_u8_geno_count(pGeno, nSamp, gSum, gNum);
			if (gNum > 0)
			{
				avg = double(gSum) / gNum;
				if (!BayesianNormal)
				{
					scale = avg * 0.5;
					scale = ((0.0 < scale) && (scale < 1.0)) ?
						(1.0 / sqrt(scale*(1.0-scale))) : 0.0;
				} else {
					scale = double(gSum + 1) / (2*gNum + 2);
					scale = 1.0 / sqrt(scale*(1.0-scale));
				}
			} else {
				avg = scale = 0;
			}
			pAFreq[i] = avg;
			pScale[i] = scale;

			// zero filling
			memset(pOut, 0, sizeof(double)*NumEigVal);
			// dot product
			for (size_t j=0; j < nSamp; j++)
			{
				double g = (*pGeno < 3) ? ((*pGeno - avg) * scale) : 0.0;
				pGeno ++;
				double *pEig = pEigVect + j;
				for (size_t k=0; k < NumEigVal; k++)
				{
					pOut[k] += g * (*pEig);
					pEig += nSamp;
				}
			}
			pOut += NumEigVal;
		}
	}

public:
	/// constructor
	CPCA_SNPLoad(CdBaseWorkSpace &space): Space(space) { }

	/// run the algorithm
	void Run(double OutSNPLoading[], double OutAFreq[], double OutScale[],
		size_t NumEig, double EigVect[], int NumThread, bool verbose)
	{
		if (NumThread < 1) NumThread = 1;
		nSamp = Space.SampleNum();
		NumEigVal = NumEig; pEigVect = EigVect;

		// detect the appropriate block size
		size_t Cache = GetOptimzedCache();
		size_t nBlock = Cache / nSamp;
		nBlock = (nBlock / 4) * 4;
		if (nBlock < 128) nBlock = 128;
		if (nBlock > 65536) nBlock = 65536;
		if (verbose)
			Rprintf("%s    (internal increment: %d)\n", TimeToStr(), (int)nBlock);

		// thread thpool
		CThreadPoolEx<CPCA_SNPLoad> thpool(NumThread);
		// genotypes (0, 1, 2 and NA)
		Geno.Reset(nSamp * nBlock);
		// genotype buffer, false for no memory buffer
		CGenoReadBySNP WS(NumThread, Space, nBlock, verbose ? -1 : 0, false);

		// for-loop
		WS.Init();
		while (WS.Read(Geno.Get()))
		{
			pLoading = OutSNPLoading + WS.Index() * NumEig;
			pAFreq = OutAFreq + WS.Index();
			pScale = OutScale + WS.Index();
			// using thread thpool
			thpool.BatchWork(this, &CPCA_SNPLoad::thread_loading, WS.Count());
			// update
			WS.ProgressForward(WS.Count());
		}
	}
};



// ================== Sample Loadings ==================

class COREARRAY_DLL_LOCAL CPCA_SampleLoad
{
private:
	CdBaseWorkSpace &Space;  ///< working genotypes
	VEC_AUTO_PTR<C_UInt8> Geno;  ///< genotypes (0, 1, 2 and NA)
	size_t nSamp;      ///< the number of samples
	size_t nEigVal;    ///< the number of eigenvalues
	size_t nSubSNP;    ///< the number of genotypes in a block
	double *pLoading;  ///< the pointer to SNP eigenvectors
	double *pAFreq;    ///< the pointer to allele frequency
	double *pScale;    ///< the pointer to scale factor
	double *pOutEig;   ///< the pointer to sample eigenvectors

	void thread_loading(size_t i, size_t num)
	{
		// for-loop each individual i
		for (; num > 0; num--, i++)
		{
			C_UInt8 *pGeno = Geno.Get() + i;
			double *pLoad = pLoading;

			for (size_t j=0; j < nSubSNP; j++)
			{
				double g = (*pGeno < 3) ? (*pGeno - pAFreq[j]) * pScale[j] : 0.0;
				pGeno += nSamp;

				double *p = pOutEig + i;
				for (size_t k=0; k < nEigVal; k++)
				{
					*p += g * (*pLoad++);
					p += nSamp;
				}
			}
		}
	}

public:
	/// constructor
	CPCA_SampleLoad(CdBaseWorkSpace &space): Space(space) { }

	/// run the algorithm
	void Run(double OutSampLoad[], size_t NumEig, double SNPLoading[],
		double AFreq[], double Scale[], int NumThread, bool verbose)
	{
		if (NumThread < 1) NumThread = 1;
		nSamp = Space.SampleNum();
		nEigVal = NumEig;
		pOutEig = OutSampLoad;

		// detect the appropriate block size
		size_t Cache = GetOptimzedCache();
		size_t nBlock = Cache / nSamp;
		nBlock = (nBlock / 4) * 4;
		if (nBlock < 128) nBlock = 128;
		if (nBlock > 65536) nBlock = 65536;
		if (verbose)
			Rprintf("%s    (internal increment: %d)\n", TimeToStr(), (int)nBlock);

		// thread thpool
		CThreadPoolEx<CPCA_SampleLoad> thpool(NumThread);
		// genotypes (0, 1, 2 and NA)
		Geno.Reset(nSamp * nBlock);
		// genotype buffer, false for no memory buffer
		CGenoReadBySNP WS(NumThread, Space, nBlock, verbose ? -1 : 0, false);

		// zero filling
		memset(OutSampLoad, 0, sizeof(double)*NumEig*nSamp);

		// for-loop
		WS.Init();
		while (WS.Read(Geno.Get()))
		{
			pLoading = SNPLoading + WS.Index() * NumEig;
			pAFreq = AFreq + WS.Index();
			pScale = Scale + WS.Index();
			nSubSNP = WS.Count();
			// using thread thpool
			thpool.BatchWork(this, &CPCA_SampleLoad::thread_loading, nSamp);
			// update
			WS.ProgressForward(WS.Count());
		}
	}
};




	// ---------------------------------------------------------------------
	// Admixture Analyses
	// ---------------------------------------------------------------------

	// ================== Relative IBD matrix ==================

	/// Missing genotypic flags:
	//      0 -- missing value, other value -- valid genotype
	static vector<C_UInt8> Admix_Missing_Flag;

	/// Adjusted genotypic values for diagonal
	static vector<double> Admix_Adj_Geno;
	/// the product of allele frequencies: 4 * p_l * (1 - p_l)
	static vector<double> Admix_P_ONE_P;

	/// Convert the raw genotypes to the mean-adjusted genotypes
	static void _Do_Admix_RatioOfAvg_ReadBlock(C_UInt8 *GenoBuf,
		long Start, long SNP_Cnt, void* Param)
	{
		// init ...
		const int n = MCWorkingGeno.Space().SampleNum();
		PCA_Mat1.ZeroFill();

		C_UInt8 *p, *pMissing;
		double *pf, *pp;
		double *pAdjGeno = &Admix_Adj_Geno[0];

		// calculate the averages of genotypes for each SelSNP
		// store the values in PCA_GenoSum
		p = GenoBuf;
		pMissing = &Admix_Missing_Flag[0];
		pp = PCA_Mat1.base();
		for (long iSample=0; iSample < n; iSample++)
		{
			pf = pp; pp += PCA_Mat1.M();
			for (long iSNP=0; iSNP < SNP_Cnt; iSNP++)
			{
				if (*p < 3)
				{
					PCA_Mat1.PCA_GenoSum[iSNP] += *p;
					PCA_Mat1.PCA_GenoNum[iSNP]++;
					*pAdjGeno += (*p) * (2 - *p);
					*pMissing++ = true;
				} else {
					*pMissing++ = false;
				}
				*pf++ = *p++;
			}
			pAdjGeno ++;
		}

		// PCA_GenoSum = 2 * \bar{p}_l
		for (long iSNP=0; iSNP < SNP_Cnt; iSNP++)
		{
			double &fv = PCA_Mat1.tmp_var.Get()[iSNP];
			long gN = PCA_Mat1.PCA_GenoNum[iSNP];
			if (gN > 0)
				fv = (double)PCA_Mat1.PCA_GenoSum[iSNP] / gN;
			else
				fv = 0;
		}

		// Averaging and set missing values to ZERO
		pp = PCA_Mat1.base(); // G@ij - 2\bar{p}_l
		for (long iSample=0; iSample < n; iSample++)
		{
			vt<double, av16Align>::Sub(pp, pp, PCA_Mat1.tmp_var.Get(), SNP_Cnt);
			pp += PCA_Mat1.M();
		}

		// missing values
		p = GenoBuf; pp = PCA_Mat1.base();
		for (long iSample=0; iSample < n; iSample++)
		{
			pf = pp; pp += PCA_Mat1.M();
			for (long iSNP=0; iSNP < SNP_Cnt; iSNP++)
			{
				if (*p > 2) *pf = 0;
				p++; pf++;
			}
		}

		// 4 * p_l * (1 - p_l) -> tmp_var
		for (long iSNP=0; iSNP < SNP_Cnt; iSNP++)
		{
			double &f = PCA_Mat1.tmp_var.Get()[iSNP];
			f *= 0.5;
			f = 4 * f * (1 - f);
		}
	}

	/// Compute the relative IBD matrix
	static void _Do_Admix_RatioOfAvg_Compute(int IdxThread,
		long Start, long SNP_Cnt, void* Param)
	{
		double **Ptr = (double**)Param;

		// numerator
		IdMatTri I = Array_Thread_MatIdx[IdxThread];
		// TODO
		PCA_Mat1.MulAdd2(I, Array_Thread_MatCnt[IdxThread], SNP_Cnt, Ptr[0] + I.Offset());

		// denominator
		I = Array_Thread_MatIdx[IdxThread];
		double *pAFreq = Ptr[1] + I.Offset();
		for (C_Int64 L = Array_Thread_MatCnt[IdxThread]; L > 0; L--)
		{
			C_UInt8 *p1 = &(Admix_Missing_Flag[0]) + SNP_Cnt*I.Row();
			C_UInt8 *p2 = &(Admix_Missing_Flag[0]) + SNP_Cnt*I.Column();
			for (long i=0; i< SNP_Cnt; i++)
			{
				if (p1[i] && p2[i])
					*pAFreq += PCA_Mat1.tmp_var.Get()[i];
			}
			pAFreq ++;
			++ I;
		}
	}

	/// Calculate the genetic covariace by ratio of averages (or ratio of sums)
	void DoAdmixCalc_RatioOfAvg(CdMatTri<double> &OutIBD, bool DiagAdj,
		int NumThread, bool verbose)
	{
		// initialize ...
		const long n = OutIBD.N();  // the number of individuals

		PCA_Mat1.Reset(n, BlockNumSNP);
		Admix_Missing_Flag.resize(BlockNumSNP * n);
		Admix_Adj_Geno.resize(n);

		memset(OutIBD.Get(), 0, sizeof(double)*OutIBD.Size());
		memset(&Admix_Adj_Geno[0], 0, sizeof(double)*n);

		MCWorkingGeno.Progress.Info = "Eigen-analysis:";
		MCWorkingGeno.Progress.Show() = verbose;
		MCWorkingGeno.InitParam(true, RDim_SNP_X_Sample, BlockNumSNP);

		CdMatTri<double> AFreqProd(n);
		memset(AFreqProd.Get(), 0, sizeof(double)*AFreqProd.Size());

		double *Ptr[2] = { OutIBD.Get(), AFreqProd.Get() };

		Array_SplitJobs(NumThread, n, Array_Thread_MatIdx,
			Array_Thread_MatCnt);
		MCWorkingGeno.Run(NumThread, &_Do_Admix_RatioOfAvg_ReadBlock,
			&_Do_Admix_RatioOfAvg_Compute, (void*)Ptr);

		// output
		double *p = OutIBD.Get();
		double *s = AFreqProd.Get();
		for (long i=0; i < n; i++)
		{
			if (DiagAdj)
				(*p) -= Admix_Adj_Geno[i];
			for (long j=i; j < n; j++)
				{ *p /= *s; p++; s++; }
		}
	}


	// ==================================== //

	/// Convert the raw genotypes to the mean-adjusted genotypes
	static void _Do_Admix_AvgOfRatio_ReadBlock(C_UInt8 *GenoBuf,
		long Start, long SNP_Cnt, void* Param)
	{
		// init ...
		const int n = MCWorkingGeno.Space().SampleNum();
		PCA_Mat1.ZeroFill();

		C_UInt8 *p, *pMissing;
		double *pf, *pp;

		// calculate the averages of genotypes for each SelSNP
		// store the values in PCA_GenoSum
		p = GenoBuf;
		pMissing = &Admix_Missing_Flag[0];
		pp = PCA_Mat1.base();
		for (long iSample=0; iSample < n; iSample++)
		{
			pf = pp; pp += PCA_Mat1.M();
			for (long iSNP=0; iSNP < SNP_Cnt; iSNP++)
			{
				if (*p < 3)
				{
					PCA_Mat1.PCA_GenoSum[iSNP] += *p;
					PCA_Mat1.PCA_GenoNum[iSNP]++;
					*pMissing++ = true;
				} else {
					*pMissing++ = false;
				}
				*pf++ = *p++;
			}
		}

		// PCA_GenoSum = 2 * \bar{p}_l
		for (long iSNP=0; iSNP < SNP_Cnt; iSNP++)
		{
			double &fv = PCA_Mat1.tmp_var.Get()[iSNP];
			long gN = PCA_Mat1.PCA_GenoNum[iSNP];
			if (gN > 0)
				fv = (double)PCA_Mat1.PCA_GenoSum[iSNP] / gN;
			else
				fv = 0;
		}

		// Averaging and set missing values to ZERO
		pp = PCA_Mat1.base(); // G@ij - 2\bar{p}_l
		for (long iSample=0; iSample < n; iSample++)
		{
			vt<double, av16Align>::Sub(pp, pp, PCA_Mat1.tmp_var.Get(), SNP_Cnt);
			pp += PCA_Mat1.M();
		}

		// 1 / (2*sqrt(p@j*(1-p@j)))
		pf = PCA_Mat1.tmp_var.Get();
		for (long iSNP=0; iSNP < SNP_Cnt; iSNP++, pf++)
		{
			double scale = (*pf) * 0.5;
			if ((0.0 < scale) && (scale < 1.0))
				*pf = 0.5 / sqrt(scale * (1.0-scale));
			else
				*pf = 0;
		}

		// G * (2 - G) / (4*p*(1-p))
		double *pAdjGeno = &Admix_Adj_Geno[0];
		p = GenoBuf;
		for (long iSample=0; iSample < n; iSample++)
		{
			pf = PCA_Mat1.tmp_var.Get();
			for (long iSNP=0; iSNP < SNP_Cnt; iSNP++)
			{
				if (*p < 3)
					*pAdjGeno += (*p) * (2 - *p) * (*pf) * (*pf);
				p++; pf++;
			}
			pAdjGeno ++;
		}

		// (G@ij - 2p@j) / (2 * sqrt(p@j*(1-p@j)))
		pp = PCA_Mat1.base();
		for (long iSample=0; iSample < n; iSample++)
		{
			vt<double, av16Align>::Mul(pp, pp, PCA_Mat1.tmp_var.Get(), SNP_Cnt);
			pp += PCA_Mat1.M();
		}

		// missing values
		p = GenoBuf; pp = PCA_Mat1.base();
		for (long iSample=0; iSample < n; iSample++)
		{
			pf = pp; pp += PCA_Mat1.M();
			for (long iSNP=0; iSNP < SNP_Cnt; iSNP++)
			{
				if (*p > 2) *pf = 0;
				p++; pf++;
			}
		}
	}

	/// Compute the relative IBD matrix
	static void _Do_Admix_AvgOfRatio_Compute(int IdxThread,
		long Start, long SNP_Cnt, void* Param)
	{
		void **Ptr = (void**)Param;
		double *Ptr1 = (double*)(Ptr[0]);

		// numerator
		IdMatTri I = Array_Thread_MatIdx[IdxThread];
		// TODO
		PCA_Mat1.MulAdd2(I, Array_Thread_MatCnt[IdxThread], SNP_Cnt,
			Ptr1 + I.Offset());

		// denominator
		I = Array_Thread_MatIdx[IdxThread];
		int *Ptr2 = (int*)(Ptr[1]) + I.Offset();
		for (C_Int64 L = Array_Thread_MatCnt[IdxThread]; L > 0; L--)
		{
			C_UInt8 *p1 = &Admix_Missing_Flag[0] + SNP_Cnt*I.Row();
			C_UInt8 *p2 = &Admix_Missing_Flag[0] + SNP_Cnt*I.Column();
			for (long i=0; i< SNP_Cnt; i++)
			{
				if ((*p1) && (*p2)) (*Ptr2) ++;
				p1 ++; p2 ++;
			}
			Ptr2 ++;
			++ I;
		}
	}

	/// Calculate the genetic covariace by averaging ratios
	void DoAdmixCalc_AvgOfRatios(CdMatTri<double> &OutIBD, bool DiagAdj,
		int NumThread, bool verbose)
	{
		// initialize ...
		const long n = OutIBD.N();  // the number of individuals

		PCA_Mat1.Reset(n, BlockNumSNP);
		Admix_Missing_Flag.resize(BlockNumSNP * n);
		Admix_Adj_Geno.resize(n);

		memset(OutIBD.Get(), 0, sizeof(double)*OutIBD.Size());
		memset(&Admix_Adj_Geno[0], 0, sizeof(double)*n);

		MCWorkingGeno.Progress.Info = "Eigen-analysis:";
		MCWorkingGeno.Progress.Show() = verbose;
		MCWorkingGeno.InitParam(true, RDim_SNP_X_Sample, BlockNumSNP);

		CdMatTri<int> NumValid(n);
		memset(NumValid.Get(), 0, sizeof(int)*NumValid.Size());

		void *Ptr[2] = { OutIBD.Get(), NumValid.Get() };

		Array_SplitJobs(NumThread, n, Array_Thread_MatIdx,
			Array_Thread_MatCnt);
		MCWorkingGeno.Run(NumThread, &_Do_Admix_AvgOfRatio_ReadBlock,
			&_Do_Admix_AvgOfRatio_Compute, (void*)Ptr);

		// output
		double *p = OutIBD.Get();
		int *s = NumValid.Get();
		for (long i=0; i < n; i++)
		{
			if (DiagAdj)
				(*p) -= Admix_Adj_Geno[i];
			for (long j=i; j < n; j++)
				{ *p /= *s; p++; s++; }
		}
	}


	// ==================================== //

	/// Convert the raw genotypes to the mean-adjusted genotypes
	static void _Do_GRM_AvgOfRatio_ReadBlock(C_UInt8 *GenoBuf,
		long Start, long SNP_Cnt, void* Param)
	{
		const static double SCALE = 1.0 / sqrt(2.0);

		// init ...
		const int n = MCWorkingGeno.Space().SampleNum();
		PCA_Mat1.ZeroFill();

		C_UInt8 *p, *pMissing;
		double *pf, *pp;

		// calculate the averages of genotypes for each SelSNP
		// store the values in PCA_GenoSum
		p = GenoBuf;
		pMissing = &Admix_Missing_Flag[0];
		pp = PCA_Mat1.base();
		for (long iSample=0; iSample < n; iSample++)
		{
			pf = pp; pp += PCA_Mat1.M();
			for (long iSNP=0; iSNP < SNP_Cnt; iSNP++)
			{
				if (*p < 3)
				{
					PCA_Mat1.PCA_GenoSum[iSNP] += *p;
					PCA_Mat1.PCA_GenoNum[iSNP]++;
					*pMissing++ = true;
				} else {
					*pMissing++ = false;
				}
				*pf++ = *p++;
			}
		}

		// PCA_GenoSum = 2 * \bar{p}_l
		for (long iSNP=0; iSNP < SNP_Cnt; iSNP++)
		{
			double &fv = PCA_Mat1.tmp_var.Get()[iSNP];
			long gN = PCA_Mat1.PCA_GenoNum[iSNP];
			if (gN > 0)
				fv = (double)PCA_Mat1.PCA_GenoSum[iSNP] / gN;
			else
				fv = 0;
		}

		// Averaging and set missing values to ZERO
		pp = PCA_Mat1.base(); // G@ij - 2\bar{p}_l
		for (long iSample=0; iSample < n; iSample++)
		{
			vt<double, av16Align>::Sub(pp, pp, PCA_Mat1.tmp_var.Get(), SNP_Cnt);
			pp += PCA_Mat1.M();
		}

		// (1 - 2p) * (G - 2p) / (2*p*(1-p))
		double *pAdjGeno = &Admix_Adj_Geno[0];
		p = GenoBuf;
		for (long iSample=0; iSample < n; iSample++)
		{
			pf = PCA_Mat1.tmp_var.Get();
			for (long iSNP=0; iSNP < SNP_Cnt; iSNP++)
			{
				double scale = (*pf) * 0.5;
				if ((0.0 < scale) && (scale < 1.0))
					scale = 1 / (2 * scale * (1 - scale));
				else
					scale = 0;
				if (*p < 3)
					*pAdjGeno += (1 - *pf) * (*p - *pf) * scale;
				p++; pf++;
			}
			pAdjGeno ++;
		}

		// 1 / (sqrt(2)*sqrt(p@j*(1-p@j)))
		pf = PCA_Mat1.tmp_var.Get();
		for (long iSNP=0; iSNP < SNP_Cnt; iSNP++, pf++)
		{
			double scale = (*pf) * 0.5;
			if ((0.0 < scale) && (scale < 1.0))
				*pf = SCALE / sqrt(scale * (1 - scale));
			else
				*pf = 0;
		}

		// (G@ij - 2p@j) / (sqrt(2) * sqrt(p@j*(1-p@j)))
		pp = PCA_Mat1.base();
		for (long iSample=0; iSample < n; iSample++)
		{
			vt<double, av16Align>::Mul(pp, pp, PCA_Mat1.tmp_var.Get(), SNP_Cnt);
			pp += PCA_Mat1.M();
		}

		// missing values
		p = GenoBuf; pp = PCA_Mat1.base();
		for (long iSample=0; iSample < n; iSample++)
		{
			pf = pp; pp += PCA_Mat1.M();
			for (long iSNP=0; iSNP < SNP_Cnt; iSNP++)
			{
				if (*p > 2) *pf = 0;
				p++; pf++;
			}
		}
	}

	/// Calculate the genetic relationship matrix (GRM)
	void DoGRMCalc(CdMatTri<double> &OutIBD, int NumThread,
		bool verbose)
	{
		// initialize ...
		const long n = OutIBD.N();  // the number of individuals

		PCA_Mat1.Reset(n, BlockNumSNP);
		Admix_Missing_Flag.resize(BlockNumSNP * n);
		Admix_Adj_Geno.resize(n);

		memset(OutIBD.Get(), 0, sizeof(double)*OutIBD.Size());
		memset(&Admix_Adj_Geno[0], 0, sizeof(double)*n);

		MCWorkingGeno.Progress.Info = "GRM-analysis:";
		MCWorkingGeno.Progress.Show() = verbose;
		MCWorkingGeno.InitParam(true, RDim_SNP_X_Sample, BlockNumSNP);

		CdMatTri<int> NumValid(n);
		memset(NumValid.Get(), 0, sizeof(int)*NumValid.Size());

		void *Ptr[2] = { OutIBD.Get(), NumValid.Get() };

		Array_SplitJobs(NumThread, n, Array_Thread_MatIdx,
			Array_Thread_MatCnt);
		MCWorkingGeno.Run(NumThread, &_Do_GRM_AvgOfRatio_ReadBlock,
			&_Do_Admix_AvgOfRatio_Compute, (void*)Ptr);

		// output
		double *p = OutIBD.Get();
		int *s = NumValid.Get();
		for (long i=0; i < n; i++)
		{
			(*p) -= Admix_Adj_Geno[i];
			for (long j=i; j < n; j++)
				{ *p /= *s; p++; s++; }
		}
	}
}


using namespace PCA;

extern "C"
{

/// get the eigenvalues and eigenvectors, return 'nProtect'
static int GetEigen(double *pMat, int n, int nEig, const char *EigMethod,
	SEXP &EigVal, SEXP &EigVect)
{
	int nProtected = 0;

	// the method to compute eigenvalues and eigenvectors
	if (strcmp(EigMethod, "DSPEV") == 0)
	{
		vector<double> tmp_Work(n*3);
		vector<double> tmp_EigenVec(n*n);

		EigVal = PROTECT(NEW_NUMERIC(n));
		nProtected ++;

		// DSPEV -- computes all the eigenvalues and, optionally,
		// eigenvectors of a real symmetric matrix A in packed storage

		int info = 0;
		F77_NAME(dspev)("V", "L", &n, pMat, REAL(EigVal),
			&tmp_EigenVec[0], &n, &tmp_Work[0], &info);
		if (info != 0)
		{
			throw ErrCoreArray(
				"LAPACK::DSPEV error (%d), infinite or missing values in the genetic covariance matrix!",
				info);
		}

		// output eigenvalues
		vt<double>::Sub(REAL(EigVal), 0.0, REAL(EigVal), n);

		// output eigenvectors
		EigVect = PROTECT(Rf_allocMatrix(REALSXP, n, nEig));
		nProtected ++;

		for (size_t i=0; i < (size_t)nEig; i++)
		{
			memmove(REAL(EigVect) + n*i,
				&tmp_EigenVec[0] + n*i, sizeof(double)*n);
		}

	} else if (strcmp(EigMethod, "DSPEVX") == 0)
	{
		vector<double> tmp_Work(n*8);
		vector<int>    tmp_IWork(n*5);

		EigVal = PROTECT(NEW_NUMERIC(n));
		EigVect = PROTECT(Rf_allocMatrix(REALSXP, n, nEig));
		nProtected += 2;

		// DSPEVX -- compute selected eigenvalues and, optionally,
		// eigenvectors of a real symmetric matrix A in packed storage

		int IL = 1, IU = nEig, LDZ = n;
		double VL = 0, VU = 0;
		int M = 0;
		// it is suggested: ABSTOL is set to twice the underflow threshold, not zero
		double ABSTOL = 2 * F77_NAME(dlamch)("S");
		vector<int> ifail(n);
		int info = 0;

		F77_NAME(dspevx)("V", "I", "L", &n, pMat, &VL, &VU, &IL, &IU, &ABSTOL,
			&M, REAL(EigVal), REAL(EigVect), &LDZ,
			&tmp_Work[0], &tmp_IWork[0], &ifail[0], &info);
		if (info != 0)
		{
			throw ErrCoreArray(
				"LAPACK::DSPEVX error (%d), infinite or missing values in the genetic covariance matrix!",
				info);
		}

		// output eigenvalues
		double *p = REAL(EigVal);
		for (int i=0; i < nEig; i++) p[i] = -p[i];
		for (int i=nEig; i < n; i++) p[i] = R_NaN;
	} else
		throw ErrCoreArray("Unknown 'eigen.method'.");

	return nProtected;
}



// ========================================================================
// Principal Component Analysis (PCA)
// ========================================================================

/// Compute the eigenvalues and eigenvectors
COREARRAY_DLL_EXPORT SEXP gnrPCA(SEXP EigenCnt, SEXP Algorithm,
	SEXP NumThread, SEXP ParamList, SEXP Verbose)
{
	const bool verbose = SEXP_Verbose(Verbose);
	int nThread = Rf_asInteger(NumThread);
	if (nThread <= 0) nThread = 1;

	COREARRAY_TRY

		// cache the genotype data
		CachingSNPData("PCA", verbose);

		// PCA algorithm: exact and fast randomized
		if (strcmp(CHAR(STRING_ELT(Algorithm, 0)), "exact") == 0)
		{
			// set parameters
			PCA::BayesianNormal =
				Rf_asLogical(RGetListElement(ParamList, "bayesian")) == TRUE;

			// the number of samples
			const R_xlen_t n = MCWorkingGeno.Space().SampleNum();

			// the upper-triangle genetic covariance matrix
			CdMatTri<double> Cov(n);

			// Calculate the genetic covariace
			const char *str = CHAR(STRING_ELT(RGetListElement(ParamList, "covalg"), 0));
			if (strcmp(str, "arith") == 0)
			{
				CExactPCA pca(MCWorkingGeno.Space());
				pca.Run(Cov, nThread, BayesianNormal, verbose);
			} else if (strcmp(str, "bitops") == 0)
			{
				PCA::CPCAMat_Alg2::PCA_Detect_BlockNumSNP(n);
				if (verbose)
					Rprintf("\tinternal increment: %d\n", (int)BlockNumSNP);
				PCA::DoCovCalculate_Alg2(Cov, nThread, "PCA:", verbose);
			} else
				throw "Invalid 'covalg'.";

			// Normalize
			double TraceXTX = Cov.Trace();
			double scale = double(n-1) / TraceXTX;
			vec_f64_mul(Cov.Get(), Cov.Size(), scale);
			double TraceVal = Cov.Trace();

			// ======== output ========

			int nProtected = 0;
			PROTECT(rv_ans = NEW_LIST(5));
			nProtected ++;

			SET_ELEMENT(rv_ans, 0, ScalarReal(TraceXTX));
			SET_ELEMENT(rv_ans, 4, ScalarReal(TraceVal));

			if (Rf_asLogical(RGetListElement(ParamList, "need.genmat")) == TRUE)
			{
				SEXP tmp;
				PROTECT(tmp = Rf_allocMatrix(REALSXP, n, n));
				nProtected ++;
				SET_ELEMENT(rv_ans, 1, tmp);
				Cov.SaveTo(REAL(tmp));
			}

			// ======== eigenvectors and eigenvalues ========

			if (Rf_asLogical(RGetListElement(ParamList, "genmat.only")) != TRUE)
			{
				if (verbose)
				{
					Rprintf("%s    Begin (eigenvalues and eigenvectors)\n",
						TimeToStr());
				}

				vt<double>::Sub(Cov.Get(), 0.0, Cov.Get(), Cov.Size());

				int nEig = Rf_asInteger(EigenCnt);
				if (nEig <= 0)
					throw ErrCoreArray("Invalid 'eigen.cnt'.");
				if (nEig > n) nEig = n;

				SEXP EigVal  = R_NilValue;
				SEXP EigVect = R_NilValue;
				nProtected += GetEigen(Cov.Get(), n, nEig,
					CHAR(STRING_ELT(RGetListElement(ParamList, "eigen.method"), 0)),
					EigVal, EigVect);
				SET_ELEMENT(rv_ans, 2, EigVal);
				SET_ELEMENT(rv_ans, 3, EigVect);
			}

			UNPROTECT(nProtected);

		} else if (strcmp(CHAR(STRING_ELT(Algorithm, 0)), "randomized") == 0)
		{
			CRandomPCA pca(MCWorkingGeno.Space(),
				REAL(RGetListElement(ParamList, "aux.mat")),
				Rf_asInteger(RGetListElement(ParamList, "aux.dim")),
				Rf_asInteger(RGetListElement(ParamList, "iter.num")));
			rv_ans = pca.Run(nThread, verbose);
		} else
			throw "Invalid 'algorithm'.";

		if (verbose)
			Rprintf("%s    Done.\n", TimeToStr());

	COREARRAY_CATCH
}


/// Calculate the SNP correlations
COREARRAY_DLL_EXPORT SEXP gnrPCACorr(SEXP LenEig, SEXP EigenVect,
	SEXP NumThread, SEXP _Verbose)
{
	const bool verbose = SEXP_Verbose(_Verbose);
	COREARRAY_TRY

		// cache the genotype data
		CachingSNPData("SNP Correlation", verbose);

		// output variable
		PROTECT(rv_ans = Rf_allocMatrix(REALSXP, Rf_asInteger(LenEig),
			MCWorkingGeno.Space().SNPNum()));
		{
			CPCA_SNPCorr Work(MCWorkingGeno.Space());
			Work.Run(REAL(rv_ans), Rf_asInteger(LenEig), REAL(EigenVect),
				Rf_asInteger(NumThread), verbose);
		}
		if (verbose)
			Rprintf("%s    Done.\n", TimeToStr());
		UNPROTECT(1);

	COREARRAY_CATCH
}


/// Calculate the SNP loadings
COREARRAY_DLL_EXPORT SEXP gnrPCASNPLoading(SEXP EigenVal, SEXP LenEig,
	SEXP EigenVect, SEXP TraceXTX, SEXP NumThread, SEXP Bayesian,
	SEXP _Verbose)
{
	const bool verbose = SEXP_Verbose(_Verbose);
	COREARRAY_TRY

		// cache the genotype data
		CachingSNPData("SNP Loading", verbose);

		// scale eigenvectors with eigenvalues
		SEXP EigVect = PROTECT(duplicate(EigenVect));
		{
			const size_t n = MCWorkingGeno.Space().SampleNum();
			const double Scale = double(n - 1) / Rf_asReal(TraceXTX);
			for (int i=0; i < Rf_asInteger(LenEig); i++)
			{
				vec_f64_mul(REAL(EigVect) + i*n, n,
					sqrt(Scale / REAL(EigenVal)[i]));
			}
		}

		PCA::BayesianNormal = (Rf_asLogical(Bayesian) == TRUE);
		const size_t n = MCWorkingGeno.Space().SNPNum();
		PROTECT(rv_ans = NEW_LIST(3));

		SEXP loading, afreq, scale;
		PROTECT(loading = Rf_allocMatrix(REALSXP, Rf_asInteger(LenEig), n));
		SET_ELEMENT(rv_ans, 0, loading);

		PROTECT(afreq = NEW_NUMERIC(n));
		SET_ELEMENT(rv_ans, 1, afreq);

		PROTECT(scale = NEW_NUMERIC(n));
		SET_ELEMENT(rv_ans, 2, scale);

		{
			CPCA_SNPLoad Work(MCWorkingGeno.Space());
			Work.Run(REAL(loading), REAL(afreq), REAL(scale),
				Rf_asInteger(LenEig), REAL(EigVect),
				Rf_asInteger(NumThread), verbose);
		}
		if (verbose)
			Rprintf("%s    Done.\n", TimeToStr());

		UNPROTECT(5);

	COREARRAY_CATCH
}


/// Calculate the sample loadings from SNP loadings
COREARRAY_DLL_EXPORT SEXP gnrPCASampLoading(SEXP EigenCnt, SEXP SNPLoadings,
	SEXP AveFreq, SEXP Scale, SEXP NumThread, SEXP _Verbose)
{
	const bool verbose = SEXP_Verbose(_Verbose);
	COREARRAY_TRY

		// cache the genotype data
		CachingSNPData("Sample Loading", verbose);

		PROTECT(rv_ans = Rf_allocMatrix(REALSXP,
			MCWorkingGeno.Space().SampleNum(), Rf_asInteger(EigenCnt)));
		{
			CPCA_SampleLoad Work(MCWorkingGeno.Space());
			Work.Run(REAL(rv_ans), Rf_asInteger(EigenCnt), REAL(SNPLoadings),
				REAL(AveFreq), REAL(Scale), Rf_asInteger(NumThread), verbose);
		}
		if (verbose)
			Rprintf("%s    Done.\n", TimeToStr());
		UNPROTECT(1);

	COREARRAY_CATCH
}



// =======================================================================
// Genetic relationship matrix
// =======================================================================

/// Compute the eigenvalues and eigenvectors
COREARRAY_DLL_EXPORT SEXP gnrGRM(SEXP _NumThread, SEXP _Method, SEXP _Verbose)
{
	const int nThread  = Rf_asInteger(_NumThread);
	const char *Method = CHAR(STRING_ELT(_Method, 0));
	const bool verbose = SEXP_Verbose(_Verbose);

	COREARRAY_TRY

		// ======== To cache the genotype data ========
		CachingSNPData("GRM Calculation", verbose);

		// ======== The calculation of genetic covariance matrix ========

		// the number of samples
		const R_xlen_t n = MCWorkingGeno.Space().SampleNum();

		// set internal parameters
		PCA::CPCAMat_Alg1::PCA_Detect_BlockNumSNP(n);

		// the upper-triangle IBD matrix
		CdMatTri<double> IBD(n);

		if (strcmp(Method, "Eigenstrat") == 0)
		{
			CExactPCA pca(MCWorkingGeno.Space());
			pca.Run(IBD, nThread, false, verbose);

			// normalize
			double TraceXTX = IBD.Trace();
			double scale = double(n-1) / TraceXTX;
			vt<double, av16Align>::Mul(IBD.Get(), IBD.Get(), scale, IBD.Size());

		} else if (strcmp(Method, "GCTA") == 0)
		{
			// Calculate Visscher's GRM
			PCA::DoGRMCalc(IBD, nThread, verbose);

		} else if (strcmp(Method, "EIGMIX") == 0)
		{
			// Calculate Zheng's coancestry matrix (EIGMIX)
			PCA::DoAdmixCalc_RatioOfAvg(IBD, true, nThread, verbose);

		} else if (strcmp(Method, "W&Z15") == 0)
		{
			const int nSNP = MCWorkingGeno.Space().SNPNum();
			CdBufSpace BufSNP(MCWorkingGeno.Space(), true, CdBufSpace::acInc);

			IBD.Clear(0);
			CdMatTri<double> Denom(n, 0);
			CdMatTri<C_Int8> M(n);

			CdProgression Progress(0, verbose);
			Progress.Init(nSNP, true);

			for (int iSNP=0; iSNP < nSNP; iSNP++)
			{
				const C_UInt8 *pG = BufSNP.ReadGeno(iSNP);
				C_Int8 *pM = M.Get();
				C_Int64 Sum = 0;
				int nSum = 0;
				for (R_xlen_t i=0; i < n; i++)
				{
					for (R_xlen_t j=i; j < n; j++, pM++)
					{
						if ((pG[i] < 3) && (pG[j] < 3))
						{
							if (j != i)
							{
								*pM = 1 + (1 - (int)pG[i]) * (1 - (int)pG[j]);
								Sum += (*pM); nSum ++;
							} else
								*pM = 2 * (1 - (int)pG[i]) * (1 - (int)pG[j]);
						} else
							*pM = -1;
					}
				}

				if (nSum > 0)
				{
					double Mb    = (double)Sum / (2*nSum);
					double OneMb = 1 - Mb;
					double *pI = IBD.Get(), *pD = Denom.Get();
					pM = M.Get();
					for (size_t i=IBD.Size(); i > 0; i--)
					{
						if (*pM >= 0)
						{
							*pI += (*pM) * 0.5 - Mb;
							*pD += OneMb;
						}
						pI ++; pD ++; pM ++;
					}
				}

				Progress.Forward(1);
			}

			// the ratio
			vt<double>::Div(IBD.Get(), IBD.Get(), Denom.Get(), IBD.Size());

		} else
			throw ErrCoreArray("Invalid 'method'!");

		// Output
		PROTECT(rv_ans = Rf_allocMatrix(REALSXP, n, n));
		double *base = REAL(rv_ans);
		double *p = IBD.Get();
		for (R_xlen_t i=0; i < n; i++)
		{
			for (R_xlen_t j=i; j < n; j++)
			{
				base[i*n + j] = base[j*n + i] = *p;
				p ++;
			}
		}
		UNPROTECT(1);

	COREARRAY_CATCH
}


// =======================================================================
// the functions for eigen-analysis for admixtures
//

struct TEigPair
{
	double EigVal;
	int Index;
	TEigPair(double e, int i) { EigVal = e; Index = i; }
};

static bool _EigComp(const TEigPair &i, const TEigPair &j)
{
	return (fabs(i.EigVal) >= fabs(j.EigVal));
}

/// to compute the eigenvalues and eigenvectors
COREARRAY_DLL_EXPORT SEXP gnrEIGMIX(SEXP _EigenCnt, SEXP _NumThread,
	SEXP _NeedIBDMat, SEXP _IBDMatOnly, SEXP _Verbose)
{
	const bool verbose = SEXP_Verbose(_Verbose);

	COREARRAY_TRY

		// ======== To cache the genotype data ========
		CachingSNPData("Eigen-analysis", verbose);

		// ======== The calculation of genetic covariance matrix ========

		// the number of samples
		const R_xlen_t n = MCWorkingGeno.Space().SampleNum();

		// set internal parameters
		PCA::CPCAMat_Alg1::PCA_Detect_BlockNumSNP(n);

		// the upper-triangle IBD matrix
		CdMatTri<double> IBD(n);

		// calculate the EIGMIX coancestry matrix
		PCA::DoAdmixCalc_RatioOfAvg(IBD, true, Rf_asInteger(_NumThread),
			verbose);

		// ======== The calculation of eigenvectors and eigenvalues ========

		int nProtected = 0;
		SEXP EigenVal=NULL, EigenVec=NULL, IBDMat=NULL;

		if (Rf_asLogical(_NeedIBDMat) == TRUE)
		{
			PROTECT(IBDMat = Rf_allocMatrix(REALSXP, n, n));
			nProtected ++;

			double *base = REAL(IBDMat);
			double *p = IBD.Get();
			for (R_xlen_t i=0; i < n; i++)
			{
				for (R_xlen_t j=i; j < n; j++)
				{
					base[i*n + j] = base[j*n + i] = *p;
					p ++;
				}
			}
		}

		if (Rf_asLogical(_IBDMatOnly) != TRUE)
		{
			const size_t NN = n;
			vector<double> tmp_Work(NN*3);
			vector<double> tmp_EigenVec(NN*NN);

			vt<double>::Sub(IBD.Get(), 0.0, IBD.Get(), IBD.Size());
			if (verbose)
			{
				Rprintf("Eigen-analysis:\t%s\tBegin (eigenvalues and eigenvectors)\n",
					NowDateToStr().c_str());
			}

			int eigencnt = Rf_asInteger(_EigenCnt);
			if (eigencnt > n) eigencnt = n;

			PROTECT(EigenVal = NEW_NUMERIC(n));
			PROTECT(EigenVec = Rf_allocMatrix(REALSXP, n, eigencnt));
			nProtected += 2;

			{
				int info = 0;
				int _n = n;
				F77_NAME(dspev)("V", "L", &_n, IBD.Get(), REAL(EigenVal),
					&tmp_EigenVec[0], &_n, &tmp_Work[0], &info);
				if (info != 0)
					throw "LAPACK::DSPEV error!";
				double *p = REAL(EigenVal);
				for (int i=0; i < n; i++)
				{
					if (!R_FINITE(p[i]))
						throw "LAPACK::DSPEV returns non-finite eigenvalues!";
				}
			}

			// sort in a decreasing order
			vector<TEigPair> lst;
			for (int i=0; i < n; i++)
				lst.push_back(TEigPair(REAL(EigenVal)[i], i));
			sort(lst.begin(), lst.end(), _EigComp);

			// output eigenvalues
			for (int i=0; i < n; i++)
				REAL(EigenVal)[i] = - lst[i].EigVal;

			// output eigenvectors
			for (int i=0; i < eigencnt; i++)
			{
				double *p = &tmp_EigenVec[0] + lst[i].Index * n;
				memmove(&(REAL(EigenVec)[i*n]), p, sizeof(double)*n);
			}

			if (verbose)
			{
				Rprintf("Eigen-analysis:\t%s\tEnd (eigenvalues and eigenvectors)\n",
					NowDateToStr().c_str());
			}
		}

		PROTECT(rv_ans = NEW_LIST(3));
		nProtected ++;

		if (EigenVal != NULL)
			SET_ELEMENT(rv_ans, 0, EigenVal);
		if (EigenVec != NULL)
			SET_ELEMENT(rv_ans, 1, EigenVec);
		if (IBDMat != NULL)
			SET_ELEMENT(rv_ans, 2, IBDMat);

		SEXP tmp;
		PROTECT(tmp = NEW_CHARACTER(3));
		nProtected ++;
			SET_STRING_ELT(tmp, 0, mkChar("eigenval"));
			SET_STRING_ELT(tmp, 1, mkChar("eigenvect"));
			SET_STRING_ELT(tmp, 2, mkChar("ibdmat"));
		SET_NAMES(rv_ans, tmp);

		UNPROTECT(nProtected);

	COREARRAY_CATCH
}

}

#endif  /* _HEADER_PCA_ */
