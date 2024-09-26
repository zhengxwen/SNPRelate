// ===========================================================
//
// genBeta.cpp: Individual Inbreeding and Relatedness (Beta) on GWAS
//
// Copyright (C) 2016-2024    Xiuwen Zheng
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


#ifndef _HEADER_IBD_BETA_
#define _HEADER_IBD_BETA_

// Standard library header
#include <cmath>
#include <cfloat>
#include <memory>
#include <algorithm>

// CoreArray library header
#include <dGenGWAS.h>
#include <dVect.h>
#include "ThreadPool.h"


namespace IBD_BETA
{

using namespace std;
using namespace CoreArray;
using namespace Vectorization;
using namespace GWAS;


// ---------------------------------------------------------------------
// Counting IBS variables for individual beta method

/// The structure of Individual Beta Estimator
struct TS_Beta
{
	C_UInt32 ibscnt;  ///< the number shared states defined in beta estimator
	C_UInt32 num;     ///< the total number of valid loci
};

class COREARRAY_DLL_LOCAL CIndivBeta
{
private:
	CdBaseWorkSpace &Space;

	size_t nBlock; /// the number of SNPs in a block, a multiple of 128
	VEC_AUTO_PTR<C_UInt8> Geno1b;  /// the genotype 1b representation
	TS_Beta *ptrBeta;

	void thread_ibs_num(size_t i, size_t n)
	{
		const size_t npack  = nBlock >> 3;
		const size_t npack2 = npack * 2;

		C_UInt8 *Base = Geno1b.Get();
		IdMatTri I = Array_Thread_MatIdx[i];
		C_Int64 N = Array_Thread_MatCnt[i];
		TS_Beta *p = ptrBeta + I.Offset();

		for (; N > 0; N--, ++I, p++)
		{
			C_UInt8 *p1 = Base + I.Row() * npack2;
			C_UInt8 *p2 = Base + I.Column() * npack2;
			ssize_t m = npack;

			#if defined(COREARRAY_SIMD_AVX2)
			{
				SIMD256_NOT_HEAD
				for (; m >= 32; m-=32)
				{
					__m256i g1_1 = _mm256_loadu_si256((__m256i*)p1);
					__m256i g1_2 = _mm256_loadu_si256((__m256i*)(p1 + npack));
					__m256i g2_1 = _mm256_loadu_si256((__m256i*)p2);
					__m256i g2_2 = _mm256_loadu_si256((__m256i*)(p2 + npack));
					p1 += 32; p2 += 32;

					__m256i mask = (g1_1 | SIMD256_NOT(g1_2)) & (g2_1 | SIMD256_NOT(g2_2));
					__m256i het  = (g1_1 ^ g1_2) | (g2_1 ^ g2_2);
					__m256i ibs2 = SIMD256_NOT(het | (g1_1 ^ g2_1));
					het &= mask;
					ibs2 &= mask;

					p->ibscnt += POPCNT_M256(het) + (POPCNT_M256(ibs2) << 1);
					p->num    += POPCNT_M256(mask);
				}
			}
			#endif

			#if defined(VECT_HARDWARE_POPCNT)
			{
				SIMD128_NOT_HEAD
				for (; m > 0; m-=16)
				{
					__m128i g1_1 = _mm_load_si128((__m128i*)p1);
					__m128i g1_2 = _mm_load_si128((__m128i*)(p1 + npack));
					__m128i g2_1 = _mm_load_si128((__m128i*)p2);
					__m128i g2_2 = _mm_load_si128((__m128i*)(p2 + npack));
					p1 += 16; p2 += 16;

					__m128i mask = (g1_1 | SIMD128_NOT(g1_2)) & (g2_1 | SIMD128_NOT(g2_2));
					__m128i het  = (g1_1 ^ g1_2) | (g2_1 ^ g2_2);
					__m128i ibs2 = SIMD128_NOT(het | (g1_1 ^ g2_1));
					het &= mask;
					ibs2 &= mask;

					p->ibscnt += POPCNT_M128(het) + (POPCNT_M128(ibs2) << 1);
					p->num    += POPCNT_M128(mask);
				}
			}
			#elif defined(COREARRAY_SIMD_SSE2)
			{
				POPCNT_SSE2_HEAD
				SIMD128_NOT_HEAD
				__m128i ibscnt4i, num4i;
				ibscnt4i = num4i = _mm_setzero_si128();

				for (; m > 0; m-=16)
				{
					__m128i g1_1 = _mm_load_si128((__m128i*)p1);
					__m128i g1_2 = _mm_load_si128((__m128i*)(p1 + npack));
					__m128i g2_1 = _mm_load_si128((__m128i*)p2);
					__m128i g2_2 = _mm_load_si128((__m128i*)(p2 + npack));
					p1 += 16; p2 += 16;

					__m128i mask = (g1_1 | SIMD128_NOT(g1_2)) & (g2_1 | SIMD128_NOT(g2_2));
					__m128i het  = (g1_1 ^ g1_2) | (g2_1 ^ g2_2);
					__m128i ibs2 = SIMD128_NOT(het | (g1_1 ^ g2_1));
					het &= mask;
					ibs2 &= mask;

					POPCNT_SSE2_RUN(het)
					ibscnt4i = _mm_add_epi32(ibscnt4i, het);

					POPCNT_SSE2_RUN(ibs2)
					ibscnt4i = _mm_add_epi32(ibscnt4i, ibs2);
					ibscnt4i = _mm_add_epi32(ibscnt4i, ibs2);

					POPCNT_SSE2_RUN(mask)
					num4i = _mm_add_epi32(num4i, mask);
				}

				p->ibscnt += vec_sum_i32(ibscnt4i);
				p->num    += vec_sum_i32(num4i);
			}
			#else
				// No SIMD
				for (; m > 0; m-=8)
				{
					C_UInt64 g1_1 = *((C_UInt64*)p1);
					C_UInt64 g1_2 = *((C_UInt64*)(p1 + npack));
					C_UInt64 g2_1 = *((C_UInt64*)p2);
					C_UInt64 g2_2 = *((C_UInt64*)(p2 + npack));
					p1 += 8; p2 += 8;

					C_UInt64 mask = (g1_1 | ~g1_2) & (g2_1 | ~g2_2);
					C_UInt64 het  = (g1_1 ^ g1_2) | (g2_1 ^ g2_2);
					C_UInt64 ibs2 = ~(het | (g1_1 ^ g2_1));
					p->ibscnt += POPCNT_U64(het & mask) + 2*POPCNT_U64(ibs2 & mask);
					p->num    += POPCNT_U64(mask);
				}
			#endif
		}
	}

public:
	/// constructor
	CIndivBeta(CdBaseWorkSpace &space): Space(space) { }

	/// run the algorithm
	void Run(CdMatTri<TS_Beta> &IBS, int NumThread, bool verbose)
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
		ptrBeta = IBS.Get();
		memset(ptrBeta, 0, sizeof(TS_Beta)*IBS.Size());

		// thread pool
		CThreadPoolEx<CIndivBeta> thpool(NumThread);
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
			thpool.BatchWork(this, &CIndivBeta::thread_ibs_num, NumThread);
			// update
			WS.ProgressForward(WS.Count());
		}
	}
};

}


extern "C"
{

using namespace IBD_BETA;

static void CPU_Flag()
{
	Rprintf("CPU capabilities:");
	#ifdef COREARRAY_SIMD_SSE2
		Rprintf(" Double-Precision SSE2");
	#endif
	#ifdef VECT_HARDWARE_POPCNT
		Rprintf(" POPCNT");
	#endif
	#ifdef COREARRAY_SIMD_AVX2
		Rprintf(" AVX2");
	#endif
	Rprintf("\n");
}


extern double grm_avg_value;

/// Compute individual beta - based GRM
COREARRAY_DLL_EXPORT void CalcIndivBetaGRM_Mat(CdMatTri<double> &beta,
	int NumThread, bool Verbose)
{
	if (Verbose) CPU_Flag();

	const size_t n = MCWorkingGeno.Space().SampleNum();
	CdMatTri<TS_Beta> IBS(n);
	CIndivBeta Work(MCWorkingGeno.Space());
	Work.Run(IBS, NumThread, Verbose);

	// output variables
	double *pBeta = beta.Get();
	TS_Beta *p = IBS.Get();
	double avg = 0;
	double min = double(p->ibscnt) / p->num - 1, r;

	// for-loop, average
	for (size_t i=0; i < n; i++)
	{
		*pBeta++ = r = double(p->ibscnt) / p->num - 1; p++;
		if (min > r) min = r;
		for (size_t j=i+1; j < n; j++)
		{
			*pBeta++ = r = (0.5 * p->ibscnt) / p->num; p++;
			avg += r;
			if (min > r) min = r;
		}
	}
	avg /= C_Int64(n) * (n-1) / 2;
	grm_avg_value = avg;

	// transformation
	pBeta = beta.Get();
	double scale = 2.0 / (1 - min);
	for (size_t i=0; i < n; i++)
	{
		*pBeta = (*pBeta - min) * scale * 0.5 + 1;
		pBeta ++;
		for (size_t j=i+1; j < n; j++, pBeta++)
			*pBeta = (*pBeta - min) * scale;
	}
}


/// Compute individual beta-based GRM, return SEXP
COREARRAY_DLL_EXPORT SEXP CalcIndivBetaGRM(int NumThread, bool Verbose)
{
	if (Verbose) CPU_Flag();

	const size_t n = MCWorkingGeno.Space().SampleNum();
	CdMatTri<TS_Beta> IBS(n);
	CIndivBeta Work(MCWorkingGeno.Space());
	Work.Run(IBS, NumThread, Verbose);

	// output variables
	SEXP rv_ans = PROTECT(Rf_allocMatrix(REALSXP, n, n));
	double *pBeta = REAL(rv_ans);
	TS_Beta *p = IBS.Get();
	double avg = 0;

	// for-loop, average
	for (size_t i=0; i < n; i++)
	{
		pBeta[i*n + i] = double(p->ibscnt) / p->num - 1; p++;
		for (size_t j=i+1; j < n; j++)
		{
			double s = (0.5 * p->ibscnt) / p->num; p++;
			pBeta[i*n + j] = s;
			avg += s;
		}
	}
	avg /= C_Int64(n) * (n-1) / 2;
	grm_avg_value = avg;

	// find the minimum of pairwise beta
	double min = pBeta[0];
	for (size_t i=0; i < n; i++)
	{
		double *p = pBeta + i*n;
		for (size_t j=i; j < n; j++)
			if (min > p[j]) min = p[j];
	}

	// transformation
	double scale = 2.0 / (1 - min);
	for (size_t i=0; i < n; i++)
	{
		pBeta[i*n + i] = (pBeta[i*n + i] - min) * scale * 0.5 + 1;
		for (size_t j=i+1; j < n; j++)
			pBeta[i*n + j] = pBeta[j*n + i] = (pBeta[i*n + j] - min) * scale;
	}

	UNPROTECT(1);
	return rv_ans;
}


/// Compute the IBD coefficients by individual relatedness beta
COREARRAY_DLL_EXPORT SEXP gnrIBD_Beta(SEXP Inbreeding, SEXP NumThread,
	SEXP useMatrix, SEXP Verbose)
{
	int inbreeding = Rf_asLogical(Inbreeding);
	if (inbreeding == NA_LOGICAL)
		Rf_error("'inbreeding' must be TRUE or FALSE.");
	bool verbose = SEXP_Verbose(Verbose);

	COREARRAY_TRY

		// cache the genotype data
		CachingSNPData("Individual Beta", verbose);
		if (verbose) CPU_Flag();

		// the number of samples
		const size_t n = MCWorkingGeno.Space().SampleNum();
		// the upper-triangle IBS matrix
		CdMatTri<TS_Beta> IBS(n);
		{
			CIndivBeta Work(MCWorkingGeno.Space());
			Work.Run(IBS, Rf_asInteger(NumThread), verbose);
		}

		// output variables
		if (Rf_asLogical(useMatrix) != TRUE)
		{
			rv_ans = PROTECT(Rf_allocMatrix(REALSXP, n, n));
			double *pBeta = REAL(rv_ans);
			TS_Beta *p = IBS.Get();
			double avg = 0;
			// for-loop, average
			for (size_t i=0; i < n; i++)
			{
				if (inbreeding)
					pBeta[i*n + i] = double(p->ibscnt) / p->num - 1;
				else
					pBeta[i*n + i] = (0.5 * p->ibscnt) / p->num;
				p ++;
				for (size_t j=i+1; j < n; j++, p++)
				{
					double s = (0.5 * p->ibscnt) / p->num;
					pBeta[i*n + j] = s;
					avg += s;
				}
			}
			avg /= C_Int64(n) * (n-1) / 2;
			grm_avg_value = avg;
			// for-loop, final update
			double bt = 1.0 / (1 - avg);
			for (size_t i=0; i < n; i++)
			{
				pBeta[i*n + i] = (pBeta[i*n + i] - avg) * bt;
				for (size_t j=i+1; j < n; j++)
					pBeta[i*n + j] = pBeta[j*n + i] = (pBeta[i*n + j] - avg) * bt;
			}
		} else {
			// triangle matrix
			const size_t ns = n*(n+1)/2;
			rv_ans = PROTECT(NEW_NUMERIC(ns));
			double *pBeta = REAL(rv_ans);
			TS_Beta *p = IBS.Get();
			double avg = 0;
			// for-loop, average
			for (size_t i=0; i < n; i++)
			{
				if (inbreeding)
					*pBeta++ = double(p->ibscnt) / p->num - 1;
				else
					*pBeta++ = (0.5 * p->ibscnt) / p->num;
				p ++;
				for (size_t j=i+1; j < n; j++, p++)
				{
					double s = (0.5 * p->ibscnt) / p->num;
					*pBeta++ = s;
					avg += s;
				}
			}
			avg /= C_Int64(n) * (n-1) / 2;
			grm_avg_value = avg;
			// for-loop, final update
			double bt = 1.0 / (1 - avg);
			pBeta = REAL(rv_ans);
			for (size_t i=0; i < n; i++)
			{
				*pBeta = (*pBeta - avg) * bt;
				pBeta++;
				for (size_t j=i+1; j < n; j++)
				{
					*pBeta = (*pBeta - avg) * bt;
					pBeta++;
				}
			}
		}

		if (verbose)
			Rprintf("%s    Done.\n", TimeToStr());
		UNPROTECT(1);

	COREARRAY_CATCH
}

}

#endif  /* _HEADER_IBD_BETA_ */
