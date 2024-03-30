// ===========================================================
//
// genIBS.cpp: Identity by state (IBS) Analysis on GWAS
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
// with SNPRelate.
// If not, see <http://www.gnu.org/licenses/>.


#ifndef _HEADER_IBS_
#define _HEADER_IBS_

// Standard library header
#include <cmath>
#include <cfloat>
#include <memory>
#include <algorithm>

// CoreArray library header
#include "dGenGWAS.h"
#include "dVect.h"
#include "ThreadPool.h"


namespace IBS
{
	// using namespace
	using namespace std;
	using namespace CoreArray;
	using namespace Vectorization;
	using namespace GWAS;


	/// Packed size
	static const long PACKED_SIZE = 256*256;

	/// IBS
	/// The number of IBS 0 in the packed genotype
	C_UInt8 IBS0_Num_SNP[PACKED_SIZE];
	/// The number of IBS 1 in the packed genotype
	C_UInt8 IBS1_Num_SNP[PACKED_SIZE];
	/// The number of IBS 2 in the packed genotype
	C_UInt8 IBS2_Num_SNP[PACKED_SIZE];

	/// Genetic dissimilarity
	/// The dissimilarity in the packed genotype
	C_UInt8 Gen_Diss_SNP[PACKED_SIZE];
	/// The flag of use of allele frequencies
	C_UInt8 Gen_Both_Valid[PACKED_SIZE];


	/// The packed genotype buffer
	vector<C_UInt8> GenoPacked;
	/// The allele frequencies
	vector<double> GenoAlleleFreq;


	/// The structure of IBS states
	struct TIBS
	{
		C_UInt32 IBS0;  //< the number of loci sharing no allele
		C_UInt32 IBS1;  //< the number of loci sharing only one allele
		C_UInt32 IBS2;  //< the number of loci sharing two alleles
	};


	/// The pointer to the variable 'PublicDiss' in the function "DoDissCalculate"
	/// The structure of genetic distance
	struct TS_Dissimilarity
	{
		C_Int64 SumGeno;
		double SumAFreq;
		TS_Dissimilarity() { SumGeno = 0; SumAFreq = 0; }
	};



	// TInit object
	class TInit
	{
	public:
		TInit()
		{
			#define PACKED_COND(cond, var, op)	\
				for (int s=0; s < PACKED_SIZE; s++)	\
				{	\
					int g1 = s/256, g2 = s%256;	\
					int sum = 0;	\
					for (int i=0; i < 4; i++)	\
					{	\
						int b1 = g1 & 0x03, b2 = g2 & 0x03;	\
						if (cond) op;	\
						g1 >>= 2; g2 >>= 2;	\
					}	\
					var[s] = sum;	\
				}

			/// The number of IBS 0 in the packed genotype
			PACKED_COND((b1 < 3) && (b2 < 3) && (abs(b1-b2)==2), IBS0_Num_SNP, sum++);
			/// The number of IBS 1 in the packed genotype
			PACKED_COND((b1 < 3) && (b2 < 3) && (abs(b1-b2)==1), IBS1_Num_SNP, sum++);
			/// The number of IBS 2 in the packed genotype
			PACKED_COND((b1 < 3) && (b2 < 3) && (abs(b1-b2)==0), IBS2_Num_SNP, sum++);
			
			/// The distance in the packed genotype
			PACKED_COND((b1 < 3) && (b2 < 3), Gen_Diss_SNP, sum += b1*(2-b2) + (2-b1)*b2);
			PACKED_COND((b1 < 3) && (b2 < 3), Gen_Both_Valid, sum |= (1 << i));
		}
	} InitObj;


	/// detect the effective value for BlockNumSNP
	static void AutoDetectSNPBlockSize(int nSamp, bool Detect=true)
	{
		if (Detect)
		{
			C_UInt64 L2Cache = GDS_Mach_GetCPULevelCache(2);
			C_UInt64 L3Cache = GDS_Mach_GetCPULevelCache(3);
			C_UInt64 Cache = (L2Cache > L3Cache) ? L2Cache : L3Cache;
			if ((C_Int64)Cache <= 0)
				Cache = 1024*1024; // 1MiB
			BlockNumSNP = (Cache - 3*256*256 - 8*1024) / nSamp * 4;
		}
		BlockNumSNP = (BlockNumSNP / 4) * 4;
		if (BlockNumSNP < 16) BlockNumSNP = 16;
	}


// ---------------------------------------------------------------------
// Counting IBS0, IBS1 and IBS2

class COREARRAY_DLL_LOCAL CIBSCount
{
private:
	CdBaseWorkSpace &Space;

	size_t nBlock; /// the number of SNPs in a block, a multiple of 128
	VEC_AUTO_PTR<C_UInt8> Geno1b;  /// the genotype 1b representation
	TIBS *ptrIBS;

	void thread_ibs_num(size_t i, size_t n)
	{
		const size_t npack  = nBlock >> 3;
		const size_t npack2 = npack * 2;

		C_UInt8 *Base = Geno1b.Get();
		IdMatTri I = Array_Thread_MatIdx[i];
		C_Int64 N = Array_Thread_MatCnt[i];
		TIBS *p = ptrIBS + I.Offset();

		for (; N > 0; N--, ++I, p++)
		{
			C_UInt8 *p1 = Base + I.Row() * npack2;
			C_UInt8 *p2 = Base + I.Column() * npack2;
			ssize_t m = npack;

		#if defined(COREARRAY_SIMD_AVX2_XXX)
		{	// disable
			const __m256i ones = _mm256_set1_epi32(-1);
			POPCNT_AVX2_HEAD
			__m256i ibs0_sum, ibs1_sum, ibs2_sum;
			ibs0_sum = ibs1_sum = ibs2_sum = _mm256_setzero_si256();

			for (; m > 0; m-=32)
			{
				__m256i g1_1 = _mm256_load_si256((__m256i*)p1);
				__m256i g1_2 = _mm256_load_si256((__m256i*)(p1 + npack));
				__m256i g2_1 = _mm256_load_si256((__m256i*)p2);
				__m256i g2_2 = _mm256_load_si256((__m256i*)(p2 + npack));
				p1 += 32; p2 += 32;

				__m256i g1_2_x = _mm256_andnot_si256(g1_2, ones);
				__m256i g2_1_x = _mm256_andnot_si256(g2_1, ones);
				__m256i g2_2_x = _mm256_andnot_si256(g2_2, ones);

				// (g1_1 | ~g1_2) & (g2_1 | ~g2_2)
				__m256i mask = _mm256_and_si256(
					_mm256_or_si256(g1_1, g1_2_x), _mm256_or_si256(g2_1, g2_2_x));
				// (~((g1_1 ^ ~g2_1) | (g1_2 ^ ~g2_2))) & mask
				__m256i ibs0 = _mm256_andnot_si256(_mm256_or_si256(
					_mm256_xor_si256(g1_1, g2_1_x), _mm256_xor_si256(g1_2, g2_2_x)), mask);
				// (~((g1_1 ^ g2_1) | (g1_2 ^ g2_2))) & mask
				__m256i ibs2 = _mm256_andnot_si256(_mm256_or_si256(
					_mm256_xor_si256(g1_1, g2_1), _mm256_xor_si256(g1_2, g2_2)), mask);

				POPCNT_AVX2_RUN(ibs0)
				ibs0_sum = _mm256_add_epi32(ibs0_sum, ibs0);

				POPCNT_AVX2_RUN(ibs2)
				ibs2_sum = _mm256_add_epi32(ibs2_sum, ibs2);

				POPCNT_AVX2_RUN(mask)
				mask = _mm256_sub_epi32(_mm256_sub_epi32(mask, ibs0), ibs2);
				ibs1_sum = _mm256_add_epi32(ibs1_sum, mask);
			}

			p->IBS0 += vec_avx_sum_i32(ibs0_sum);
			p->IBS1 += vec_avx_sum_i32(ibs1_sum);
			p->IBS2 += vec_avx_sum_i32(ibs2_sum);
		}
			if (m > 0)
		#endif
		#if defined(COREARRAY_SIMD_SSE2)
		{
			POPCNT_SSE2_HEAD
			SIMD128_NOT_HEAD
			__m128i ibs0_sum, ibs1_sum, ibs2_sum;
			ibs0_sum = ibs1_sum = ibs2_sum = _mm_setzero_si128();

			for (; m > 0; m-=16)
			{
				__m128i g1_1 = _mm_load_si128((__m128i*)p1);
				__m128i g1_2 = _mm_load_si128((__m128i*)(p1 + npack));
				__m128i g2_1 = _mm_load_si128((__m128i*)p2);
				__m128i g2_2 = _mm_load_si128((__m128i*)(p2 + npack));
				p1 += 16; p2 += 16;

				__m128i mask = (g1_1 | SIMD128_NOT(g1_2)) & (g2_1 | SIMD128_NOT(g2_2));
				__m128i ibs0 = SIMD128_NOT((g1_1 ^ SIMD128_NOT(g2_1)) | (g1_2 ^ SIMD128_NOT(g2_2))) & mask;
				__m128i ibs2 = SIMD128_NOT((g1_1 ^ g2_1) | (g1_2 ^ g2_2)) & mask;

				POPCNT_SSE2_RUN(ibs0)
				ibs0_sum = _mm_add_epi32(ibs0_sum, ibs0);

				POPCNT_SSE2_RUN(ibs2)
				ibs2_sum = _mm_add_epi32(ibs2_sum, ibs2);

				POPCNT_SSE2_RUN(mask)
				mask = _mm_sub_epi32(_mm_sub_epi32(mask, ibs0), ibs2);
				ibs1_sum = _mm_add_epi32(ibs1_sum, mask);
			}

			p->IBS0 += vec_sum_i32(ibs0_sum);
			p->IBS1 += vec_sum_i32(ibs1_sum);
			p->IBS2 += vec_sum_i32(ibs2_sum);
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
				C_UInt64 ibs2 = (~((g1_1 ^ g2_1) | (g1_2 ^ g2_2))) & mask;

				size_t i0 = POPCNT_U64(ibs0);
				size_t i2 = POPCNT_U64(ibs2);
				size_t i1 = POPCNT_U64(mask) - i0 - i2;

				p->IBS0 += i0;
				p->IBS1 += i1;
				p->IBS2 += i2;
			}
		#endif
		}
	}

public:
	/// constructor
	CIBSCount(CdBaseWorkSpace &space): Space(space) { }

	/// run the algorithm
	void Run(CdMatTri<TIBS> &IBS, int NumThread, bool verbose)
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
		ptrIBS = IBS.Get();
		memset(ptrIBS, 0, sizeof(TIBS)*IBS.Size());

		// thread pool
		CThreadPoolEx<CIBSCount> thpool(NumThread);
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
			thpool.BatchWork(this, &CIBSCount::thread_ibs_num, NumThread);
			// update
			WS.ProgressForward(WS.Count());
		}
	}
};



	/// ======================================================================
	/// Individual Dissimilarity
	/// ======================================================================

	/// Convert the raw genotypes
	static void _Do_Diss_ReadBlock(C_UInt8 *GenoBuf, long Start,
		long SNP_Cnt, void* Param)
	{
		// initialize
		const int nSamp = MCWorkingGeno.Space().SampleNum();
		C_UInt8 *pG = GenoBuf;
		C_UInt8 *pPack = &GenoPacked[0];

		// pack genotypes
		for (long iSamp=0; iSamp < nSamp; iSamp++)
		{
			pPack = PackGeno2b(pG, SNP_Cnt, pPack);
			pG += SNP_Cnt;
		}
		// calculate the allele frequencies
		for (long iSNP=0; iSNP < SNP_Cnt; iSNP++)
		{
			C_UInt8 *p = GenoBuf + iSNP;
			double &Freq = GenoAlleleFreq[iSNP];
			int n = 0; Freq = 0;
			for (long iSamp=0; iSamp < nSamp; iSamp++)
			{
				if (*p < 3) { Freq += *p; n += 2; }
				p += SNP_Cnt;
			}
			Freq = (n > 0) ? Freq/n : 0;
			Freq = 8 * Freq * (1 - Freq);
		}
	}

	/// Compute the covariate matrix
	static void _Do_Diss_Compute(int ThreadIndex, long Start,
		long SNP_Cnt, void* Param)
	{
		long Cnt = Array_Thread_MatCnt[ThreadIndex];
		IdMatTri I = Array_Thread_MatIdx[ThreadIndex];
		TS_Dissimilarity *p = ((TS_Dissimilarity*)Param) + I.Offset();
		long _PackSNPLen = (SNP_Cnt / 4) + (SNP_Cnt % 4 ? 1 : 0);

		for (; Cnt > 0; Cnt--, ++I, p++)
		{
			C_UInt8 *p1 = &GenoPacked[0] + I.Row()*_PackSNPLen;
			C_UInt8 *p2 = &GenoPacked[0] + I.Column()*_PackSNPLen;
			for (long k=0; k < _PackSNPLen; k++, p1++, p2++)
			{
				size_t t = (size_t(*p1) << 8) | (*p2);
				p->SumGeno += Gen_Diss_SNP[t];

				C_UInt8 flag = Gen_Both_Valid[t];
				if (flag & 0x01)
					p->SumAFreq += GenoAlleleFreq[4*k];
				if (flag & 0x02)
					p->SumAFreq += GenoAlleleFreq[4*k+1];
				if (flag & 0x04)
					p->SumAFreq += GenoAlleleFreq[4*k+2];
				if (flag & 0x08)
					p->SumAFreq += GenoAlleleFreq[4*k+3];
			}
		}
	}

	/// Calculate the genetic distance matrix
	void DoDissCalculate(CdMatTri<TS_Dissimilarity> &PublicDiss, int NumThread,
		const char *Info, bool verbose)
	{
		// Initialize ...
		GenoPacked.resize(BlockNumSNP * PublicDiss.N());
		void *p = (void*)PublicDiss.Get();
		size_t n = sizeof(TS_Dissimilarity)*PublicDiss.Size();
		memset(p, 0, n);
		GenoAlleleFreq.resize(BlockNumSNP);

		MCWorkingGeno.Progress.Info = Info;
		MCWorkingGeno.Progress.Show() = verbose;
		MCWorkingGeno.InitParam(true, RDim_SNP_X_Sample, BlockNumSNP);

		Array_SplitJobs(NumThread, PublicDiss.N(),
			Array_Thread_MatIdx, Array_Thread_MatCnt);
		MCWorkingGeno.Run(NumThread, &_Do_Diss_ReadBlock,
			&_Do_Diss_Compute, PublicDiss.Get());
	}
}

namespace IBD
{
	// PLINK method of moment
	COREARRAY_DLL_LOCAL void Init_EPrIBD_IBS(const double in_afreq[],
		double out_afreq[], bool CorrectFactor, long nSNP = -1);

	COREARRAY_DLL_LOCAL void Est_PLINK_Kinship(int IBS0, int IBS1, int IBS2,
		double &k0, double &k1, bool KinshipConstrict);
}


using namespace IBS;

extern "C"
{
// =======================================================================
// Identity-By-State (IBS)
//

/// Compute the average IBS
COREARRAY_DLL_EXPORT SEXP gnrIBSAve(SEXP NumThread, SEXP useMatrix,
	SEXP _Verbose)
{
	const bool verbose = SEXP_Verbose(_Verbose);

	COREARRAY_TRY

		// cache the genotype data
		CachingSNPData("IBS", verbose);

		// the number of samples
		const size_t n = MCWorkingGeno.Space().SampleNum();

		// the upper-triangle IBS matrix
		CdMatTri<IBS::TIBS> IBS(n);
		{
			// Calculate the IBS matrix
			CIBSCount Work(MCWorkingGeno.Space());
			Work.Run(IBS, Rf_asInteger(NumThread), verbose);
		}

		// output variables
		if (Rf_asLogical(useMatrix) != TRUE)
		{
			rv_ans = PROTECT(Rf_allocMatrix(REALSXP, n, n));
			double *pIBS = REAL(rv_ans);
			IBS::TIBS *p = IBS.Get();
			for (size_t i=0; i < n; i++)
			{
				for (size_t j=i; j < n; j++, p++)
				{
					pIBS[i*n + j] = pIBS[j*n + i] =
						(0.5 * p->IBS1 + p->IBS2) / (p->IBS0 + p->IBS1 + p->IBS2);
				}
			}
		} else {
			// triangle matrix
			const size_t ns = n*(n+1)/2;
			rv_ans = PROTECT(NEW_NUMERIC(ns));
			double *pIBS = REAL(rv_ans);
			IBS::TIBS *p = IBS.Get();
			for (size_t i=0; i < n; i++)
			{
				for (size_t j=i; j < n; j++, p++)
				{
					*pIBS++ = (0.5 * p->IBS1 + p->IBS2) /
						(p->IBS0 + p->IBS1 + p->IBS2);
				}
			}
		}

		if (verbose)
			Rprintf("%s    Done.\n", TimeToStr());
		UNPROTECT(1);

	COREARRAY_CATCH
}

/// Compute the counts of IBS0, IBS1, IBS2
COREARRAY_DLL_EXPORT SEXP gnrIBSNum(SEXP NumThread, SEXP _Verbose)
{
	const bool verbose = SEXP_Verbose(_Verbose);

	COREARRAY_TRY

		// cache the genotype data
		CachingSNPData("IBS", verbose);

		// the number of samples
		const size_t n = MCWorkingGeno.Space().SampleNum();

		// the upper-triangle IBS matrix
		CdMatTri<IBS::TIBS> IBS(n);
		{
			// Calculate the IBS matrix
			CIBSCount Work(MCWorkingGeno.Space());
			Work.Run(IBS, Rf_asInteger(NumThread), verbose);
		}

		// output variables
		SEXP IBS0 = PROTECT(Rf_allocMatrix(INTSXP, n, n));
		SEXP IBS1 = PROTECT(Rf_allocMatrix(INTSXP, n, n));
		SEXP IBS2 = PROTECT(Rf_allocMatrix(INTSXP, n, n));

		PROTECT(rv_ans = NEW_LIST(3));
		SET_ELEMENT(rv_ans, 0, IBS0);
		SET_ELEMENT(rv_ans, 1, IBS1);
		SET_ELEMENT(rv_ans, 2, IBS2);

		int *pIBS0 = INTEGER(IBS0);
		int *pIBS1 = INTEGER(IBS1);
		int *pIBS2 = INTEGER(IBS2);

		IBS::TIBS *p = IBS.Get();
		for (size_t i=0; i < n; i++)
		{
			for (size_t j=i; j < n; j++, p++)
			{
				pIBS0[i*n + j] = pIBS0[j*n + i] = p->IBS0;
				pIBS1[i*n + j] = pIBS1[j*n + i] = p->IBS1;
				pIBS2[i*n + j] = pIBS2[j*n + i] = p->IBS2;
			}
		}

		if (verbose)
			Rprintf("%s    Done.\n", TimeToStr());
		UNPROTECT(4);

	COREARRAY_CATCH
}


// =======================================================================
// PLINK Identity-By-Descent (IBD)
//

/// Compute the IBD coefficients by PLINK method of moment
COREARRAY_DLL_EXPORT SEXP gnrIBD_PLINK(SEXP NumThread, SEXP AlleleFreq,
	SEXP UseSpecificAFreq, SEXP KinshipConstrict, SEXP useMatrix, SEXP _Verbose)
{
	const bool kc = Rf_asLogical(KinshipConstrict)==TRUE;
	const bool verbose = SEXP_Verbose(_Verbose);

	COREARRAY_TRY

		// cache the genotype data
		CachingSNPData("PLINK IBD", verbose);

		// the number of individuals
		const size_t n = MCWorkingGeno.Space().SampleNum();

		// the upper-triangle IBS matrix
		CdMatTri<IBS::TIBS> IBS(n);
		{
			// Calculate the IBS matrix
			CIBSCount Work(MCWorkingGeno.Space());
			Work.Run(IBS, Rf_asInteger(NumThread), verbose);
		}

		// output variables
		SEXP afreq = PROTECT(NEW_NUMERIC(MCWorkingGeno.Space().SNPNum()));
		double *out_afreq = REAL(afreq);

		// initialize the internal matrix
		IBD::Init_EPrIBD_IBS(
			Rf_asLogical(UseSpecificAFreq)==TRUE ? REAL(AlleleFreq) : NULL,
			out_afreq, Rf_asLogical(UseSpecificAFreq)!=TRUE);

		SEXP k0, k1;
		if (Rf_asLogical(useMatrix) != TRUE)
		{
			k0 = PROTECT(Rf_allocMatrix(REALSXP, n, n));
			k1 = PROTECT(Rf_allocMatrix(REALSXP, n, n));
			double *out_k0=REAL(k0), *out_k1=REAL(k1);
			IBS::TIBS *p = IBS.Get();
			for (size_t i=0; i < n; i++)
			{
				out_k0[i*n + i] = out_k1[i*n + i] = 0;
				p++;
				for (size_t j=i+1; j < n; j++, p++)
				{
					double k0, k1;
					IBD::Est_PLINK_Kinship(p->IBS0, p->IBS1, p->IBS2, k0, k1, kc);
					out_k0[i*n + j] = out_k0[j*n + i] = k0;
					out_k1[i*n + j] = out_k1[j*n + i] = k1;
				}
			}
		} else {
			// triangle matrix
			const size_t ns = n*(n+1)/2;
			k0 = PROTECT(NEW_NUMERIC(ns));
			k1 = PROTECT(NEW_NUMERIC(ns));
			double *out_k0=REAL(k0), *out_k1=REAL(k1);
			IBS::TIBS *p = IBS.Get();
			for (size_t i=0; i < n; i++)
			{
				*out_k0++ = *out_k1++ = 0;
				p++;
				for (size_t j=i+1; j < n; j++, p++)
				{
					double k0, k1;
					IBD::Est_PLINK_Kinship(p->IBS0, p->IBS1, p->IBS2, k0, k1, kc);
					*out_k0++ = k0; *out_k1++ = k1;
				}
			}
		}

		// output
		PROTECT(rv_ans = NEW_LIST(3));
		SET_ELEMENT(rv_ans, 0, k0);
		SET_ELEMENT(rv_ans, 1, k1);
		SET_ELEMENT(rv_ans, 2, afreq);

		if (verbose)
			Rprintf("%s    Done.\n", TimeToStr());
		UNPROTECT(4);

	COREARRAY_CATCH
}



// =========================================================================
// Individual dissimilarity

/// Compute the individual dissimilarity
COREARRAY_DLL_EXPORT SEXP gnrDiss(SEXP NumThread, SEXP _Verbose)
{
	const bool verbose = SEXP_Verbose(_Verbose);

	COREARRAY_TRY

		// cache the genotype data
		CachingSNPData("Dissimilarity", verbose);

		// the number of samples
		const R_xlen_t n = MCWorkingGeno.Space().SampleNum();
		// to detect the block size
		IBS::AutoDetectSNPBlockSize(n);
		// the upper-triangle genetic covariance matrix
		CdMatTri<IBS::TS_Dissimilarity> Dist(n);

		// Calculate the genetic distance matrix
		IBS::DoDissCalculate(Dist, INTEGER(NumThread)[0], "Dissimilarity:",
			verbose);

		// output
		rv_ans = PROTECT(Rf_allocMatrix(REALSXP, n, n));

		IBS::TS_Dissimilarity *p = Dist.Get();
		double *out_Diss = REAL(rv_ans);
		for (R_xlen_t i=0; i < n; i++)
		{
			out_Diss[i*n + i] = 2 * (p->SumGeno / p->SumAFreq);
			p ++;
			for (R_xlen_t j=i+1; j < n; j++, p++)
				out_Diss[i*n + j] = out_Diss[j*n + i] = (p->SumGeno / p->SumAFreq);
		}

		UNPROTECT(1);

	COREARRAY_CATCH
}



// =========================================================================
// the functions for genotype score

static void flap_allele(C_UInt8 *pGeno, const int *c_1, const int *c_2,
	int nSamp, int nPair)
{
	int n=0, gsum=0;
	for (int j=0; j < nPair; j++)
	{
		C_UInt8 g1 = pGeno[c_1[j]], g2 = pGeno[c_2[j]];
		if (g1 < 3) { n++; gsum+=g1; }
		if (g2 < 3) { n++; gsum+=g2; }
	}
	if (gsum < n)
	{
		for (int j=0; j < nSamp; j++)
		{
			C_UInt8 &g = pGeno[j];
			if (g < 3) g = 2 - g;
		}
	}
}

/// Compute the individual dissimilarity
COREARRAY_DLL_EXPORT SEXP gnrPairScore(SEXP SampIdx1, SEXP SampIdx2,
	SEXP Method, SEXP Type, SEXP Dosage, SEXP GDSNode, SEXP _Verbose)
{
	typedef int TMat[4][4];

	const int M = -1;
	static const TMat IBS =
		{ { 2, 1, 0, M }, { 1, 2, 1, M }, { 0, 1, 2, M }, { M, M, M, M } };
	static const TMat IBS_1 =
		{ { 1, 1, 0, M }, { 1, 1, 1, M }, { 0, 1, 1, M }, { M, M, M, M } };
	static const TMat GVH =
		{ { 0, 0, 2, M }, { 1, 0, 1, M }, { 2, 0, 0, M }, { M, M, M, M } };
	static const TMat GVH_1 =
		{ { 0, 0, 1, M }, { 1, 0, 1, M }, { 1, 0, 0, M }, { M, M, M, M } };
	static const TMat HVG =
		{ { 0, 1, 2, M }, { 0, 0, 0, M }, { 2, 1, 0, M }, { M, M, M, M } };
	static const TMat HVG_1 =
		{ { 0, 1, 1, M }, { 0, 0, 0, M }, { 1, 1, 0, M }, { M, M, M, M } };
	static const TMat GVH_major =
		{ { 0, 0, 0, M }, { 1, 0, 0, M }, { 1, 0, 0, M }, { M, M, M, M } };
	static const TMat GVH_minor =
		{ { 0, 0, 1, M }, { 0, 0, 1, M }, { 0, 0, 0, M }, { M, M, M, M } };
	static const TMat GVH_major_only =
		{ { 0, 0, M, M }, { 1, 0, M, M }, { 1, 0, 0, M }, { M, M, M, M } };
	static const TMat GVH_minor_only =
		{ { 0, 0, 1, M }, { M, 0, 1, M }, { M, 0, 0, M }, { M, M, M, M } };

	const int *c_1 = INTEGER(SampIdx1);
	const int *c_2 = INTEGER(SampIdx2);
	const int nPair = Rf_length(SampIdx1);
	const char *c_Method = CHAR(STRING_ELT(Method, 0));
	const char *c_Type = CHAR(STRING_ELT(Type, 0));
	const bool verbose = SEXP_Verbose(_Verbose);

	const int dosage = Rf_asLogical(Dosage);
	if (dosage == NA_LOGICAL)
		error("'dosage' must be TRUE or FALSE.");

	COREARRAY_TRY

		const TMat *m;
		bool need_major = false;
		if (strcmp(c_Method, "IBS") == 0)
		{
			m = dosage ? &IBS : &IBS_1;
		} else if (strcmp(c_Method, "GVH") == 0)
		{
			m = dosage ? &GVH : &GVH_1;
		} else if (strcmp(c_Method, "HVG") == 0)
		{
			m = dosage ? &HVG : &HVG_1;
		} else if (strcmp(c_Method, "GVH.major") == 0)
		{
			m = &GVH_major; need_major = true;
		} else if (strcmp(c_Method, "GVH.minor") == 0)
		{
			m = &GVH_minor; need_major = true;
		} else if (strcmp(c_Method, "GVH.major.only") == 0)
		{
			m = &GVH_major_only; need_major = true;
		} else if (strcmp(c_Method, "GVH.minor.only") == 0)
		{
			m = &GVH_minor_only; need_major = true;
		} else
			throw ErrCoreArray("Invalid 'method'.");
		const TMat &map = *m;

		// ======== To cache the genotype data ========
		CachingSNPData("Genotype Score", verbose);

		// ======== Genotype score for individual pairs ========
		const int nSNP = MCWorkingGeno.Space().SNPNum();
		const int nSamp = MCWorkingGeno.Space().SampleNum();
		CdBufSpace BufSNP(MCWorkingGeno.Space(), true, CdBufSpace::acInc);

		if (strcmp(c_Type, "per.pair") == 0)
		{
			vector<CSummary_AvgSD> Sums(nPair);
			// for-loop
			for (int i=0; i < nSNP; i++)
			{
				C_UInt8 *pGeno = BufSNP.ReadGeno(i);
				if (need_major)
					flap_allele(pGeno, c_1, c_2, nSamp, nPair);
				for (int j=0; j < nPair; j++)
				{
					C_UInt8 g1 = pGeno[c_1[j]];
					C_UInt8 g2 = pGeno[c_2[j]];
					if ((g1 < 3) && (g2 < 3))
						Sums[j].Add(map[g1][g2]);
				}
			}

			rv_ans = Rf_allocMatrix(REALSXP, nPair, 3);
			double *Out = REAL(rv_ans);
			for (int i=0; i < nPair; i++)
			{
				Sums[i].CalcAvgSD();
				double *p = Out + i;
				p[0      ] = Sums[i].Avg;
				p[nPair  ] = Sums[i].SD;
				p[nPair*2] = Sums[i].Num;
			}
		} else if (strcmp(c_Type, "per.snp") == 0)
		{
			vector<double> Buffer(nPair);
			rv_ans = Rf_allocMatrix(REALSXP, 3, nSNP);
			double *Out = REAL(rv_ans);

			// for-loop
			for (int i=0; i < nSNP; i++)
			{
				C_UInt8 *pGeno = BufSNP.ReadGeno(i);
				if (need_major)
					flap_allele(pGeno, c_1, c_2, nSamp, nPair);
				for (int j=0; j < nPair; j++)
				{
					C_UInt8 g1 = pGeno[c_1[j]];
					C_UInt8 g2 = pGeno[c_2[j]];
					if ((g1 < 3) && (g2 < 3))
						Buffer[j] = map[g1][g2];
					else
						Buffer[j] = R_NaN;
				}

				CSummary_AvgSD Sum;
				Sum.Add(&Buffer[0], nPair);
				Sum.CalcAvgSD();
				Out[0] = Sum.Avg; Out[1] = Sum.SD; Out[2] = Sum.Num;
				Out += 3;
			}
		} else if (strcmp(c_Type, "matrix") == 0)
		{
			rv_ans = Rf_allocMatrix(INTSXP, nPair, nSNP);
			int *Out = INTEGER(rv_ans);

			// for-loop
			for (int i=0; i < nSNP; i++)
			{
				C_UInt8 *pGeno = BufSNP.ReadGeno(i);
				if (need_major)
					flap_allele(pGeno, c_1, c_2, nSamp, nPair);
				for (int j=0; j < nPair; j++, Out++)
				{
					C_UInt8 g1 = pGeno[c_1[j]];
					C_UInt8 g2 = pGeno[c_2[j]];
					if ((g1 < 3) && (g2 < 3))
						*Out = map[g1][g2];
					else
						*Out = NA_INTEGER;
				}
			}
		} else if (strcmp(c_Type, "gds.file") == 0)
		{
			PdAbstractArray Obj = GDS_R_SEXP2Obj(GDSNode, FALSE);
			vector<C_UInt8> Buffer(nPair);

			// for-loop
			for (int i=0; i < nSNP; i++)
			{
				C_UInt8 *pGeno = BufSNP.ReadGeno(i);
				if (need_major)
					flap_allele(pGeno, c_1, c_2, nSamp, nPair);
				C_UInt8 *Out = &Buffer[0];
				for (int j=0; j < nPair; j++, Out++)
				{
					C_UInt8 g1 = pGeno[c_1[j]];
					C_UInt8 g2 = pGeno[c_2[j]];
					if ((g1 < 3) && (g2 < 3))
						*Out = map[g1][g2];
					else
						*Out = 3; // NA
				}

				GDS_Array_AppendData(Obj, nPair, &Buffer[0], svUInt8);
			}
		} else
			throw ErrCoreArray("Invalid 'type'.");

	COREARRAY_CATCH
}

}

#endif  /* _HEADER_IBS_ */
