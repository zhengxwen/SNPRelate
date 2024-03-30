// ===========================================================
//
// genKING.cpp: KINK Identity-by-descent (IBD) Analysis on GWAS
//
// Copyright (C) 2011-2024    Xiuwen Zheng
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
// with SNPRelate.  If not, see <http://www.gnu.org/licenses/>.


#ifndef _HEADER_IBD_KING_
#define _HEADER_IBD_KING_

// Standard library header
#include <cmath>
#include <cfloat>
#include <memory>
#include <algorithm>

// CoreArray library header
#include <dGenGWAS.h>
#include <dVect.h>
#include "ThreadPool.h"


namespace IBD_KING
{

using namespace std;
using namespace CoreArray;
using namespace Vectorization;
using namespace GWAS;


// ---------------------------------------------------------------------
// Counting IBS variables for KING homo method

/// The structure of KING IBD Homo estimator
struct TS_KINGHomo
{
	C_UInt32 IBS0;     ///< the number of loci sharing no allele
	C_UInt32 SumSq;    ///< \sum_m (g_m^{(i)} - g_m^{(j)})^2
	double SumAFreq;   ///< \sum_m p_m (1 - p_m)
	double SumAFreq2;  ///< \sum_m p_m^2 (1 - p_m)^2
};

class COREARRAY_DLL_LOCAL CKINGHomo
{
private:
	CdBaseWorkSpace &Space;

	size_t nBlock; /// the number of SNPs in a block, a multiple of 128
	VEC_AUTO_PTR<C_UInt8> Geno1b;  /// the genotype 1b representation
	VEC_AUTO_PTR<double> AF_1_AF;    /// p * (1 - p)
	VEC_AUTO_PTR<double> AF_1_AF_2;  /// (p * (1 - p))^2
	TS_KINGHomo *ptrKING;

	void thread_ibs_num(size_t i, size_t n)
	{
		const size_t npack  = nBlock >> 3;
		const size_t npack2 = npack * 2;

		C_UInt8 *Base = Geno1b.Get();
		IdMatTri I = Array_Thread_MatIdx[i];
		C_Int64 N = Array_Thread_MatCnt[i];
		TS_KINGHomo *p = ptrKING + I.Offset();

		for (; N > 0; N--, ++I, p++)
		{
			C_UInt8 *p1 = Base + I.Row() * npack2;
			C_UInt8 *p2 = Base + I.Column() * npack2;
			double *pAF  = AF_1_AF.Get();
			double *pAF2 = AF_1_AF_2.Get();
			ssize_t m = npack;

		#if defined(COREARRAY_SIMD_SSE2)
		{
			POPCNT_SSE2_HEAD
			SIMD128_NOT_HEAD
			__m128i ibs0_sum, sumsq_sum;
			ibs0_sum = sumsq_sum = _mm_setzero_si128();
			__m128d sq_sum, sq_sum2;
			sq_sum = sq_sum2 = _mm_setzero_pd();

			for (; m > 0; m-=16)
			{
				__m128i g1_1 = _mm_load_si128((__m128i*)p1);
				__m128i g1_2 = _mm_load_si128((__m128i*)(p1 + npack));
				__m128i g2_1 = _mm_load_si128((__m128i*)p2);
				__m128i g2_2 = _mm_load_si128((__m128i*)(p2 + npack));
				p1 += 16; p2 += 16;

				__m128i mask = (g1_1 | SIMD128_NOT(g1_2)) & (g2_1 | SIMD128_NOT(g2_2));
				__m128i ibs0 = SIMD128_NOT((g1_1 ^ SIMD128_NOT(g2_1)) | (g1_2 ^ SIMD128_NOT(g2_2))) & mask;
				__m128i het  = ((g1_1 ^ g1_2) ^ (g2_1 ^ g2_2)) & mask;

				POPCNT_SSE2_RUN(ibs0)
				ibs0_sum = _mm_add_epi32(ibs0_sum, ibs0);

				POPCNT_SSE2_RUN(het)
				sumsq_sum = _mm_add_epi32(_mm_add_epi32(sumsq_sum, het),
					_mm_slli_epi32(ibs0, 2));

				C_UInt64 m[2] COREARRAY_SIMD_ATTR_ALIGN;
				_mm_store_si128((__m128i*)m, mask);
				for (size_t k=32; k > 0; k--)
				{
					switch (m[0] & 0x03)
					{
					case 3:
						sq_sum = _mm_add_pd(sq_sum, _mm_load_pd(pAF));
						sq_sum2 = _mm_add_pd(sq_sum2, _mm_load_pd(pAF2));
						break;
					case 1:
						sq_sum = _mm_add_pd(sq_sum, _mm_set_pd(0, pAF[0]));
						sq_sum2 = _mm_add_pd(sq_sum2, _mm_set_pd(0, pAF2[0]));
						break;
					case 2:
						sq_sum = _mm_add_pd(sq_sum, _mm_set_pd(pAF[1], 0));
						sq_sum2 = _mm_add_pd(sq_sum2, _mm_set_pd(pAF2[1], 0));
						break;
					}
					pAF += 2; pAF2 += 2; m[0] >>= 2;
				}
				for (size_t k=32; k > 0; k--)
				{
					switch (m[1] & 0x03)
					{
					case 3:
						sq_sum = _mm_add_pd(sq_sum, _mm_load_pd(pAF));
						sq_sum2 = _mm_add_pd(sq_sum2, _mm_load_pd(pAF2));
						break;
					case 1:
						sq_sum = _mm_add_pd(sq_sum, _mm_set_pd(0, pAF[0]));
						sq_sum2 = _mm_add_pd(sq_sum2, _mm_set_pd(0, pAF2[0]));
						break;
					case 2:
						sq_sum = _mm_add_pd(sq_sum, _mm_set_pd(pAF[1], 0));
						sq_sum2 = _mm_add_pd(sq_sum2, _mm_set_pd(pAF2[1], 0));
						break;
					}
					pAF += 2; pAF2 += 2; m[1] >>= 2;
				}
			}

			p->IBS0  += vec_sum_i32(ibs0_sum);
			p->SumSq += vec_sum_i32(sumsq_sum);
			p->SumAFreq += vec_sum_f64(sq_sum);
			p->SumAFreq2 += vec_sum_f64(sq_sum2);
		}
		#else
			for (; m > 0; m-=8)
			{
				C_UInt64 g1_1 = *((C_UInt64*)p1);
				C_UInt64 g1_2 = *((C_UInt64*)(p1 + npack));
				C_UInt64 g2_1 = *((C_UInt64*)p2);
				C_UInt64 g2_2 = *((C_UInt64*)(p2 + npack));
				p1 += 8; p2 += 8;

				C_UInt64 mask = (g1_1 | ~g1_2) & (g2_1 | ~g2_2);
				C_UInt64 ibs0 = (~((g1_1 ^ ~g2_1) | (g1_2 ^ ~g2_2))) & mask;
				C_UInt64 het  = ((g1_1 ^ g1_2) ^ (g2_1 ^ g2_2)) & mask;

				p->IBS0  += POPCNT_U64(ibs0);
				p->SumSq += POPCNT_U64(het) + POPCNT_U64(ibs0)*4;

				double sum=0, sum2=0;
				for (size_t k=64; k > 0; k--)
				{
					if (mask & 0x01)
						{ sum += (*pAF); sum2 += (*pAF2); }
					pAF ++; pAF2 ++;
					mask >>= 1;
				}
				p->SumAFreq += sum;
				p->SumAFreq2 += sum2;
			}
		#endif
		}
	}

public:
	/// constructor
	CKINGHomo(CdBaseWorkSpace &space): Space(space) { }

	/// run the algorithm
	void Run(CdMatTri<TS_KINGHomo> &IBS, int NumThread, bool verbose)
	{
		if (NumThread < 1) NumThread = 1;
		const size_t nSamp = Space.SampleNum();

		// detect the appropriate block size
		nBlock = 4 * GetOptimzedCache() / nSamp;
		nBlock = (nBlock / 128) * 128;
		if (nBlock < 256) nBlock = 256;
		if (nBlock > 65536) nBlock = 65536;
		const size_t nPack = nBlock / 8;
		if (verbose)
			Rprintf("%s    (internal increment: %d)\n", TimeToStr(), (int)nBlock);

		// initialize
		ptrKING = IBS.Get();
		memset(ptrKING, 0, sizeof(TS_KINGHomo)*IBS.Size());

		// thread pool
		CThreadPoolEx<CKINGHomo> thpool(NumThread);
		Array_SplitJobs(NumThread, nSamp, Array_Thread_MatIdx,
			Array_Thread_MatCnt);

		// genotypes
		Geno1b.Reset(nSamp * nBlock / 4);
		VEC_AUTO_PTR<C_UInt8> Geno(nSamp * nBlock);

		AF_1_AF.Reset(nBlock);
		AF_1_AF_2.Reset(nBlock);

		// genotype buffer, false for no memory buffer
		CGenoReadBySNP WS(NumThread, Space, nBlock, verbose ? -1 : 0, false);

		// for-loop
		WS.Init();
		while (WS.Read(Geno.Get()))
		{
			// computing "p * (1 - p)" and "(p * (1 - p))^2"
			double *pAF = AF_1_AF.Get(), *pAF2 = AF_1_AF_2.Get();
			C_UInt8 *pG = Geno.Get();
			for (size_t m=WS.Count(); m > 0; m--)
			{
				C_Int32 sum, num;
				vec_u8_geno_count(pG, nSamp, sum, num);
				pG += nSamp;
				double p = (num > 0) ? 0.5*sum/num : 0;
				double s = p * (1 - p);
				*pAF++ = s; *pAF2++ = s * s;
			}
			for (size_t m=WS.Count(); m < nBlock; m++)
				*pAF++ = *pAF2++ = 0;

			// pack genotypes
			pG = Geno.Get();
			C_UInt8 *pB = Geno1b.Get();
			for (size_t m=nSamp; m > 0; m--)
			{
				PackSNPGeno1b(pB, pB + nPack, pG, WS.Count(), nSamp, nBlock);
				pB += (nPack << 1);
				pG ++;
			}

			// using thread thpool
			thpool.BatchWork(this, &CKINGHomo::thread_ibs_num, NumThread);
			// update
			WS.ProgressForward(WS.Count());
		}
	}
};



// ---------------------------------------------------------------------
// Counting IBS variables for KING robust method

/// The structure of KING IBD Robust Estimator
struct TS_KINGRobust
{
	C_UInt32 IBS0;   ///< the number of loci sharing no allele
	C_UInt32 nLoci;  ///< the total number of loci
	C_UInt32 SumSq;  ///< \sum_m (g_m^{(i)} - g_m^{(j)})^2
	C_UInt32 N1_Aa;  ///< the number of hetet loci for the first individual
	C_UInt32 N2_Aa;  ///< the number of hetet loci for the second individual
};

class COREARRAY_DLL_LOCAL CKINGRobust
{
private:
	CdBaseWorkSpace &Space;

	size_t nBlock; /// the number of SNPs in a block, a multiple of 128
	VEC_AUTO_PTR<C_UInt8> Geno1b;  /// the genotype 1b representation
	TS_KINGRobust *ptrKING;

	void thread_ibs_num(size_t i, size_t n)
	{
		const size_t npack  = nBlock >> 3;
		const size_t npack2 = npack * 2;

		C_UInt8 *Base = Geno1b.Get();
		IdMatTri I = Array_Thread_MatIdx[i];
		C_Int64 N = Array_Thread_MatCnt[i];
		TS_KINGRobust *p = ptrKING + I.Offset();

		for (; N > 0; N--, ++I, p++)
		{
			C_UInt8 *p1 = Base + I.Row() * npack2;
			C_UInt8 *p2 = Base + I.Column() * npack2;
			ssize_t m = npack;

		#if defined(COREARRAY_SIMD_AVX2)
			POPCNT_AVX2_HEAD
			SIMD256_NOT_HEAD
			__m256i ibs0_sum, nloci_sum, sumsq_sum, n1_Aa, n2_Aa;
			ibs0_sum = nloci_sum = sumsq_sum = n1_Aa = n2_Aa = _mm256_setzero_si256();

			for (; m >= 32; m-=32)
			{
				__m256i g1_1 = _mm256_load_si256((__m256i*)p1);
				__m256i g1_2 = _mm256_load_si256((__m256i*)(p1 + npack));
				__m256i g2_1 = _mm256_load_si256((__m256i*)p2);
				__m256i g2_2 = _mm256_load_si256((__m256i*)(p2 + npack));
				p1 += 32; p2 += 32;
				__m256i r_g1_2 = SIMD256_NOT(g1_2);
				__m256i r_g2_1 = SIMD256_NOT(g2_1);
				__m256i r_g2_2 = SIMD256_NOT(g2_2);

				__m256i mask = (g1_1 | r_g1_2) & (g2_1 | r_g2_2);
				__m256i ibs0 = SIMD256_NOT((g1_1 ^ r_g2_1) | (g1_2 ^ r_g2_2)) & mask;
				__m256i het  = ((g1_1 ^ g1_2) ^ (g2_1 ^ g2_2)) & mask;
				__m256i Aa1  = g1_1 & r_g1_2 & mask;
				__m256i Aa2  = g2_1 & r_g2_2 & mask;

				POPCNT_AVX2_RUN(ibs0)
				ibs0_sum = _mm256_add_epi32(ibs0_sum, ibs0);

				POPCNT_AVX2_RUN(mask)
				nloci_sum = _mm256_add_epi32(nloci_sum, mask);

				POPCNT_AVX2_RUN(het)
				sumsq_sum = _mm256_add_epi32(_mm256_add_epi32(sumsq_sum, het),
					_mm256_slli_epi32(ibs0, 2));

				POPCNT_AVX2_RUN(Aa1)
				n1_Aa = _mm256_add_epi32(n1_Aa, Aa1);

				POPCNT_AVX2_RUN(Aa2)
				n2_Aa = _mm256_add_epi32(n2_Aa, Aa2);
			}

			p->IBS0  += vec_avx_sum_i32(ibs0_sum);
			p->nLoci += vec_avx_sum_i32(nloci_sum);
			p->SumSq += vec_avx_sum_i32(sumsq_sum);
			p->N1_Aa += vec_avx_sum_i32(n1_Aa);
			p->N2_Aa += vec_avx_sum_i32(n2_Aa);

			if (m >= 16)
		#endif
		#if defined(COREARRAY_SIMD_SSE2)
		{
			POPCNT_SSE2_HEAD
			SIMD128_NOT_HEAD
			__m128i ibs0_sum, nloci_sum, sumsq_sum, n1_Aa, n2_Aa;
			ibs0_sum = nloci_sum = sumsq_sum = n1_Aa = n2_Aa = _mm_setzero_si128();

			for (; m > 0; m-=16)
			{
				__m128i g1_1 = _mm_load_si128((__m128i*)p1);
				__m128i g1_2 = _mm_load_si128((__m128i*)(p1 + npack));
				__m128i g2_1 = _mm_load_si128((__m128i*)p2);
				__m128i g2_2 = _mm_load_si128((__m128i*)(p2 + npack));
				p1 += 16; p2 += 16;
				__m128i r_g1_2 = SIMD128_NOT(g1_2);
				__m128i r_g2_1 = SIMD128_NOT(g2_1);
				__m128i r_g2_2 = SIMD128_NOT(g2_2);

				__m128i mask = (g1_1 | r_g1_2) & (g2_1 | r_g2_2);
				__m128i ibs0 = SIMD128_NOT((g1_1 ^ r_g2_1) | (g1_2 ^ r_g2_2)) & mask;
				__m128i het  = ((g1_1 ^ g1_2) ^ (g2_1 ^ g2_2)) & mask;
				__m128i Aa1  = g1_1 & r_g1_2 & mask;
				__m128i Aa2  = g2_1 & r_g2_2 & mask;

				POPCNT_SSE2_RUN(ibs0)
				ibs0_sum = _mm_add_epi32(ibs0_sum, ibs0);

				POPCNT_SSE2_RUN(mask)
				nloci_sum = _mm_add_epi32(nloci_sum, mask);

				POPCNT_SSE2_RUN(het)
				sumsq_sum = _mm_add_epi32(_mm_add_epi32(sumsq_sum, het),
					_mm_slli_epi32(ibs0, 2));

				POPCNT_SSE2_RUN(Aa1)
				n1_Aa = _mm_add_epi32(n1_Aa, Aa1);

				POPCNT_SSE2_RUN(Aa2)
				n2_Aa = _mm_add_epi32(n2_Aa, Aa2);
			}

			p->IBS0  += vec_sum_i32(ibs0_sum);
			p->nLoci += vec_sum_i32(nloci_sum);
			p->SumSq += vec_sum_i32(sumsq_sum);
			p->N1_Aa += vec_sum_i32(n1_Aa);
			p->N2_Aa += vec_sum_i32(n2_Aa);
		}
		#else
			for (; m > 0; m-=8)
			{
				C_UInt64 g1_1 = *((C_UInt64*)p1);
				C_UInt64 g1_2 = *((C_UInt64*)(p1 + npack));
				C_UInt64 g2_1 = *((C_UInt64*)p2);
				C_UInt64 g2_2 = *((C_UInt64*)(p2 + npack));
				p1 += 8; p2 += 8;

				C_UInt64 mask = (g1_1 | ~g1_2) & (g2_1 | ~g2_2);
				C_UInt64 ibs0 = (~((g1_1 ^ ~g2_1) | (g1_2 ^ ~g2_2))) & mask;
				C_UInt64 het  = ((g1_1 ^ g1_2) ^ (g2_1 ^ g2_2)) & mask;
				C_UInt64 Aa1  = g1_1 & ~g1_2 & mask;
				C_UInt64 Aa2  = g2_1 & ~g2_2 & mask;

				p->IBS0  += POPCNT_U64(ibs0);
				p->nLoci += POPCNT_U64(mask);
				p->SumSq += POPCNT_U64(het) + POPCNT_U64(ibs0)*4;
				p->N1_Aa += POPCNT_U64(Aa1);
				p->N2_Aa += POPCNT_U64(Aa2);
			}
		#endif
		}
	}

public:
	/// constructor
	CKINGRobust(CdBaseWorkSpace &space): Space(space) { }

	/// run the algorithm
	void Run(CdMatTri<TS_KINGRobust> &IBS, int NumThread, bool verbose)
	{
		if (NumThread < 1) NumThread = 1;
		const size_t nSamp = Space.SampleNum();

		// detect the appropriate block size
		nBlock = 4 * GetOptimzedCache() / nSamp;
		nBlock = (nBlock / 128) * 128;
		if (nBlock < 256) nBlock = 256;
		if (nBlock > 65536) nBlock = 65536;
		const size_t nPack = nBlock / 8;
		if (verbose)
			Rprintf("%s    (internal increment: %d)\n", TimeToStr(), (int)nBlock);

		// initialize
		ptrKING = IBS.Get();
		memset(ptrKING, 0, sizeof(TS_KINGRobust)*IBS.Size());

		// thread pool
		CThreadPoolEx<CKINGRobust> thpool(NumThread);
		Array_SplitJobs(NumThread, nSamp, Array_Thread_MatIdx,
			Array_Thread_MatCnt);

		// genotypes
		Geno1b.Reset(nSamp * nBlock / 4);
		VEC_AUTO_PTR<C_UInt8> Geno(nSamp * nBlock);

		// genotype buffer, false for no memory buffer
		CGenoReadBySNP WS(NumThread, Space, nBlock, verbose ? -1 : 0, false);

		// for-loop
		WS.Init();
		while (WS.Read(Geno.Get()))
		{
			C_UInt8 *pG = Geno.Get();
			C_UInt8 *pB = Geno1b.Get();
			for (size_t m=nSamp; m > 0; m--)
			{
				PackSNPGeno1b(pB, pB + nPack, pG, WS.Count(), nSamp, nBlock);
				pB += (nPack << 1);
				pG ++;
			}

			// using thread thpool
			thpool.BatchWork(this, &CKINGRobust::thread_ibs_num, NumThread);
			// update
			WS.ProgressForward(WS.Count());
		}
	}
};

}


extern "C"
{

using namespace IBD_KING;

/// Compute the IBD coefficients by KING method of moment (KING-homo)
COREARRAY_DLL_EXPORT SEXP gnrIBD_KING_Homo(SEXP NumThread, SEXP useMatrix,
	SEXP _Verbose)
{
	bool verbose = SEXP_Verbose(_Verbose);

	COREARRAY_TRY

		// cache the genotype data
		CachingSNPData("KING IBD", verbose);

		// the number of samples
		const size_t n = MCWorkingGeno.Space().SampleNum();

		// the upper-triangle IBS matrix
		CdMatTri<TS_KINGHomo> IBS(n);
		{
			CKINGHomo Work(MCWorkingGeno.Space());
			Work.Run(IBS, Rf_asInteger(NumThread), verbose);
		}

		// output variables
		SEXP K0, K1;
		if (Rf_asLogical(useMatrix) != TRUE)
		{
			K0 = PROTECT(Rf_allocMatrix(REALSXP, n, n));
			K1 = PROTECT(Rf_allocMatrix(REALSXP, n, n));
			double *pK0 = REAL(K0);
			double *pK1 = REAL(K1);
			// output
			TS_KINGHomo *p = IBS.Get();
			for (size_t i=0; i < n; i++)
			{
				pK0[i*n + i] = pK1[i*n + i] = 0;
				p ++;
				for (size_t j=i+1; j < n; j++, p++)
				{
					double theta = 0.5 - p->SumSq / (8 * p->SumAFreq);
					double k0 = p->IBS0 / (2 * p->SumAFreq2);
					double k1 = 2 - 2*k0 - 4*theta;
					pK0[i*n + j] = pK0[j*n + i] = R_FINITE(k0) ? k0 : R_NaN;
					pK1[i*n + j] = pK1[j*n + i] = R_FINITE(k1) ? k1 : R_NaN;
				}
			}
		} else {
			// triangle matrix
			const size_t ns = n*(n+1)/2;
			K0 = PROTECT(NEW_NUMERIC(ns));
			K1 = PROTECT(NEW_NUMERIC(ns));
			double *pK0 = REAL(K0);
			double *pK1 = REAL(K1);
			// output
			TS_KINGHomo *p = IBS.Get();
			for (size_t i=0; i < n; i++)
			{
				*pK0++ = *pK1++ = 0;
				p ++;
				for (size_t j=i+1; j < n; j++, p++)
				{
					double theta = 0.5 - p->SumSq / (8 * p->SumAFreq);
					double k0 = p->IBS0 / (2 * p->SumAFreq2);
					double k1 = 2 - 2*k0 - 4*theta;
					*pK0++ = R_FINITE(k0) ? k0 : R_NaN;
					*pK1++ = R_FINITE(k1) ? k1 : R_NaN;
				}
			}
		}

		// output
		PROTECT(rv_ans = NEW_LIST(2));
		SET_ELEMENT(rv_ans, 0, K0);
		SET_ELEMENT(rv_ans, 1, K1);

		if (verbose)
			Rprintf("%s    Done.\n", TimeToStr());
		UNPROTECT(3);

	COREARRAY_CATCH
}




/// Compute the IBD coefficients by KING method of moment (KING-robust)
COREARRAY_DLL_EXPORT SEXP gnrIBD_KING_Robust(SEXP FamilyID, SEXP NumThread,
	SEXP useMatrix, SEXP _Verbose)
{
	bool verbose = SEXP_Verbose(_Verbose);

	COREARRAY_TRY

		// cache the genotype data
		CachingSNPData("KING IBD", verbose);
		if (verbose)
		{
			Rprintf("CPU capabilities:");
		#ifdef COREARRAY_SIMD_SSE2
			Rprintf(" Double-Precision SSE2");
		#endif
		#ifdef COREARRAY_SIMD_AVX2
			Rprintf(" AVX2");
		#endif
			Rprintf("\n");
		}

		// check the number of SNPs, 2^32 / 4
		if (MCWorkingGeno.Space().SNPNum() >= 1073741824)
		{
			throw ErrCoreArray(
				"The number of SNPs should be less than 1,073,741,824.");
		}

		// the number of samples
		const size_t n = MCWorkingGeno.Space().SampleNum();

		// the upper-triangle IBS matrix
		CdMatTri<TS_KINGRobust> IBS(n);
		{
			CKINGRobust Work(MCWorkingGeno.Space());
			Work.Run(IBS, Rf_asInteger(NumThread), verbose);
		}

		// output variables
		SEXP IBS0, Kinship;
		int *FamID_Ptr = INTEGER(FamilyID);

		if (Rf_asLogical(useMatrix) != TRUE)
		{
			IBS0 = PROTECT(Rf_allocMatrix(REALSXP, n, n));
			Kinship = PROTECT(Rf_allocMatrix(REALSXP, n, n));
			double *pIBS0 = REAL(IBS0);
			double *pKinship = REAL(Kinship);
			// output
			TS_KINGRobust *p = IBS.Get();
			for (size_t i=0; i < n; i++)
			{
				pIBS0[i*n + i] = 0; pKinship[i*n + i] = 0.5;
				p ++;
				for (size_t j=i+1; j < n; j++, p++)
				{
					pIBS0[i*n + j] = pIBS0[j*n + i] = (p->nLoci > 0) ?
						(double(p->IBS0) / p->nLoci) : R_NaN;
					int f1 = FamID_Ptr[i], f2 = FamID_Ptr[j];
					double v = (f1==f2 && f1!=NA_INTEGER) ?
						(0.5 - p->SumSq / (2.0 *(p->N1_Aa + p->N2_Aa))) :
						(0.5 - p->SumSq / (4.0 * min(p->N1_Aa, p->N2_Aa)));
					if (!R_FINITE(v)) v = R_NaN;
					pKinship[i*n + j] = pKinship[j*n + i] = v;
				}
			}
		} else {
			// triangle matrix
			const size_t ns = n*(n+1)/2;
			IBS0 = PROTECT(NEW_NUMERIC(ns));
			Kinship = PROTECT(NEW_NUMERIC(ns));
			double *pIBS0 = REAL(IBS0);
			double *pKinship = REAL(Kinship);
			// output
			TS_KINGRobust *p = IBS.Get();
			for (size_t i=0; i < n; i++)
			{
				*pIBS0++ = 0; *pKinship++ = 0.5;
				p ++;
				for (size_t j=i+1; j < n; j++, p++)
				{
					*pIBS0++ = (p->nLoci > 0) ?
						(double(p->IBS0) / p->nLoci) : R_NaN;
					int f1 = FamID_Ptr[i], f2 = FamID_Ptr[j];
					double v = (f1==f2 && f1!=NA_INTEGER) ?
						(0.5 - p->SumSq / (2.0 *(p->N1_Aa + p->N2_Aa))) :
						(0.5 - p->SumSq / (4.0 * min(p->N1_Aa, p->N2_Aa)));
					if (!R_FINITE(v)) v = R_NaN;
					*pKinship++ = v;
				}
			}
		}

		// output
		PROTECT(rv_ans = NEW_LIST(2));
		SET_ELEMENT(rv_ans, 0, IBS0);
		SET_ELEMENT(rv_ans, 1, Kinship);

		if (verbose)
			Rprintf("%s    Done.\n", TimeToStr());
		UNPROTECT(3);

	COREARRAY_CATCH
}

}

#endif  /* _HEADER_IBD_KING_ */
