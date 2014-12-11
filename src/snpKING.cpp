// ===========================================================
//     _/_/_/   _/_/_/  _/_/_/_/    _/_/_/_/  _/_/_/   _/_/_/
//      _/    _/       _/             _/    _/    _/   _/   _/
//     _/    _/       _/_/_/_/       _/    _/    _/   _/_/_/
//    _/    _/       _/             _/    _/    _/   _/
// _/_/_/   _/_/_/  _/_/_/_/_/     _/     _/_/_/   _/_/
// ===========================================================
//
// snpKING.cpp: KINK Identity-by-descent (IBD) analysis on GWAS
//
// Copyright (C) 2011 - 2015	Xiuwen Zheng [zhengxwen@gmail.com]
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


#ifndef _H_IBD_KING_
#define _H_IBD_KING_

// CoreArray library header
#include <dGenGWAS.h>
#include <dVect.h>

// Standard library header
#include <cmath>
#include <cfloat>
#include <memory>
#include <algorithm>


#ifdef COREARRAY_SIMD_SSE
#include <xmmintrin.h>
#endif
#ifdef COREARRAY_SIMD_SSE2
#include <emmintrin.h>
#endif


namespace KING_IBD
{
	// using namespace
	using namespace std;
	using namespace CoreArray;
	using namespace CoreArray::Vectorization;
	using namespace GWAS;


	/// Packed size
	static const long PACKED_SIZE = 256*256;

	//================    KING robust estimator    ================//

	/// The number of IBS 0 in the packed genotype
	C_UInt8 IBS0_Num_SNP[PACKED_SIZE];

	/// The square value of genotype difference, (X_m^{(i)} - X_m^{(j)})^2
	C_UInt8 Gen_KING_SqDiff[PACKED_SIZE];
	/// N1_Aa requiring both genotypes are available
	C_UInt8 Gen_KING_N1_Aa[PACKED_SIZE];
	/// N2_Aa requiring both genotypes are available
	C_UInt8 Gen_KING_N2_Aa[PACKED_SIZE];
	/// the number of valid loci
	C_UInt8 Gen_KING_Num_Loci[PACKED_SIZE];
	/// The flag of use of allele frequencies
	C_UInt8 Gen_Both_Valid[PACKED_SIZE];


	/// The packed genotype buffer
	vector<C_UInt8> GenoPacked;
	/// The allele frequencies
	vector<double> GenoAlleleFreq;


	/// Thread variables
	const int N_MAX_THREAD = 256;

	// IBS, KING IBD, Individual Similarity
	IdMatTri IBS_Thread_MatIdx[N_MAX_THREAD];
	C_Int64 IBS_Thread_MatCnt[N_MAX_THREAD];


	/// The pointer to the variable 'PublicKING' in the function "DoKINGCalculate"
	/// The structure of KING IBD estimator
	struct TS_KINGHomo
	{
		C_UInt32 IBS0;       //< the number of loci sharing no allele
		C_UInt32 SumSq;      //< \sum_m (X_m^{(i)} - X_m^{(j)})^2
		double SumAFreq;   //< \sum_m p_m (1 - p_m)
		double SumAFreq2;  //< \sum_m p_m^2 (1 - p_m)^2
		TS_KINGHomo() { IBS0 = SumSq = 0; SumAFreq = SumAFreq2 = 0; }
	};

	struct TS_KINGRobust
	{
		C_UInt32 IBS0;       //< the number of loci sharing no allele
		C_UInt32 nLoci;      //< the total number of loci
		C_UInt32 SumSq;      //< \sum_m (X_m^{(i)} - X_m^{(j)})^2
		C_UInt32 N1_Aa;      //< the number of hetet loci for the first individual
		C_UInt32 N2_Aa;      //< the number of hetet loci for the second individual
		TS_KINGRobust() { IBS0 = nLoci = SumSq = N1_Aa = N2_Aa = 0; }
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

			/// \sum_m (X_m^{(i)} - X_m^{(j)})^2
			PACKED_COND((b1 < 3) && (b2 < 3) && (abs(b1-b2)==2), IBS0_Num_SNP, sum++);
			PACKED_COND((b1 < 3) && (b2 < 3), Gen_KING_SqDiff, sum += (b1-b2)*(b1-b2));
			PACKED_COND((b1 < 3) && (b2 < 3), Gen_KING_N1_Aa, sum += (b1==1) ? 1:0);
			PACKED_COND((b1 < 3) && (b2 < 3), Gen_KING_N2_Aa, sum += (b2==1) ? 1:0);
			PACKED_COND((b1 < 3) && (b2 < 3), Gen_KING_Num_Loci, sum ++);
			PACKED_COND((b1 < 3) && (b2 < 3), Gen_Both_Valid, sum |= (1 << i));
		}
	} InitObj;


	/// detect the effective value for BlockSNP
	void AutoDetectSNPBlockSize(int nSamp, bool Detect=true)
	{
		if (Detect)
		{
			C_UInt64 L2Cache = GDS_Mach_GetCPULevelCache(2);
			C_UInt64 L3Cache = GDS_Mach_GetCPULevelCache(3);
			C_UInt64 Cache = (L2Cache > L3Cache) ? L2Cache : L3Cache;
			if ((C_Int64)Cache <= 0) Cache = 1024*1024; // 1M
			BlockSNP = (Cache - 3*256*256 - 8*1024) / nSamp * 4;
		}
		BlockSNP = (BlockSNP / 4) * 4;
		if (BlockSNP < 16) BlockSNP = 16;
	}



	/// ======================================================================
	/// KING Robust Estimator
	/// ======================================================================

	/// Convert the raw genotypes
	static void _Do_KING_ReadBlock(C_UInt8 *GenoBuf, long Start,
		long SNP_Cnt, void *Param)
	{
		// initialize
		const long nSamp = MCWorkingGeno.Space.SampleNum();
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
			Freq = 0;
			long n = 0;
			for (long iSamp=0; iSamp < nSamp; iSamp++)
			{
				if (*p < 3) { Freq += *p; n += 2; }
				p += SNP_Cnt;
			}
			Freq = (n > 0) ? Freq/n : 0;
			Freq = Freq * (1 - Freq);
		}
	}

	/// Compute IBD estimator in KING-homo
	static void _Do_KING_Homo_Compute(int ThreadIndex, long Start,
		long SNP_Cnt, void *Param)
	{
		long Cnt = IBS_Thread_MatCnt[ThreadIndex];
		IdMatTri I = IBS_Thread_MatIdx[ThreadIndex];
		TS_KINGHomo *p = ((TS_KINGHomo*)Param) + I.Offset();
		long _PackSNPLen = (SNP_Cnt / 4) + (SNP_Cnt % 4 ? 1 : 0);

		for (; Cnt > 0; Cnt--, ++I, p++)
		{
			C_UInt8 *p1 = &GenoPacked[0] + I.Row()*_PackSNPLen;
			C_UInt8 *p2 = &GenoPacked[0] + I.Column()*_PackSNPLen;
			for (long k=0; k < _PackSNPLen; k++, p1++, p2++)
			{
				size_t t = (size_t(*p1) << 8) | (*p2);

				p->IBS0 += IBS0_Num_SNP[t];
				p->SumSq += Gen_KING_SqDiff[t];

				C_UInt8 flag = Gen_Both_Valid[t];
				if (flag & 0x01)
				{
					double f = GenoAlleleFreq[4*k + 0];
					p->SumAFreq += f; p->SumAFreq2 += f*f;
				}
				if (flag & 0x02)
				{
					double f = GenoAlleleFreq[4*k + 1];
					p->SumAFreq += f; p->SumAFreq2 += f*f;
				}
				if (flag & 0x04)
				{
					double f = GenoAlleleFreq[4*k + 2];
					p->SumAFreq += f; p->SumAFreq2 += f*f;
				}
				if (flag & 0x08)
				{
					double f = GenoAlleleFreq[4*k + 3];
					p->SumAFreq += f; p->SumAFreq2 += f*f;
				}
			}
		}
	}

	/// Compute IBD estimator in KING-robust
	static void _Do_KING_Robust_Compute(int ThreadIndex, long Start,
		long SNP_Cnt, void* Param)
	{
		long Cnt = IBS_Thread_MatCnt[ThreadIndex];
		IdMatTri I = IBS_Thread_MatIdx[ThreadIndex];
		TS_KINGRobust *p = ((TS_KINGRobust*)Param) + I.Offset();
		long _PackSNPLen = (SNP_Cnt / 4) + (SNP_Cnt % 4 ? 1 : 0);

		for (; Cnt > 0; Cnt--, ++I, p++)
		{
			C_UInt8 *p1 = &GenoPacked[0] + I.Row()*_PackSNPLen;
			C_UInt8 *p2 = &GenoPacked[0] + I.Column()*_PackSNPLen;
			for (long k=0; k < _PackSNPLen; k++, p1++, p2++)
			{
				size_t t = (size_t(*p1) << 8) | (*p2);
				p->IBS0 += IBS0_Num_SNP[t];
				p->nLoci += Gen_KING_Num_Loci[t];
				p->SumSq += Gen_KING_SqDiff[t];
				p->N1_Aa += Gen_KING_N1_Aa[t];
				p->N2_Aa += Gen_KING_N2_Aa[t];
			}
		}
	}

	/// Calculate KING IBD Homo Estimator
	void DoKINGCalculate(CdMatTri<TS_KINGHomo> &PublicKING, int NumThread,
		const char *Info, bool verbose)
	{
		// Initialize ...
		GenoPacked.resize(BlockSNP * PublicKING.N());
		memset(PublicKING.get(), 0, sizeof(TS_KINGHomo)*PublicKING.Size());
		GenoAlleleFreq.resize(BlockSNP);

		MCWorkingGeno.Progress.Info = Info;
		MCWorkingGeno.Progress.Show() = verbose;
		MCWorkingGeno.InitParam(true, true, BlockSNP);

		MCWorkingGeno.SplitJobs(NumThread, PublicKING.N(),
			IBS_Thread_MatIdx, IBS_Thread_MatCnt);
		MCWorkingGeno.Run(NumThread, &_Do_KING_ReadBlock,
			&_Do_KING_Homo_Compute, PublicKING.get());
	}

	/// Calculate KING IBD Robust Estimator
	void DoKINGCalculate(CdMatTri<TS_KINGRobust> &PublicKING, int NumThread,
		const char *Info, bool verbose)
	{
		// Initialize ...
		GenoPacked.resize(BlockSNP * PublicKING.N());
		memset(PublicKING.get(), 0, sizeof(TS_KINGRobust)*PublicKING.Size());
		GenoAlleleFreq.resize(BlockSNP);

		MCWorkingGeno.Progress.Info = Info;
		MCWorkingGeno.Progress.Show() = verbose;
		MCWorkingGeno.InitParam(true, true, BlockSNP);

		MCWorkingGeno.SplitJobs(NumThread, PublicKING.N(),
			IBS_Thread_MatIdx, IBS_Thread_MatCnt);
		MCWorkingGeno.Run(NumThread, &_Do_KING_ReadBlock,
			&_Do_KING_Robust_Compute, PublicKING.get());
	}
}

using namespace KING_IBD;

extern "C"
{
/// to compute the IBD coefficients by KING method of moment (KING-homo)
COREARRAY_DLL_EXPORT SEXP gnrIBD_KING_Homo(SEXP NumThread, SEXP _Verbose)
{
	bool verbose = SEXP_Verbose(_Verbose);

	COREARRAY_TRY

		// =======* To cache the genotype data =======*
		CachingSNPData("KING IBD", verbose);

		// =======* The calculation of genetic covariance matrix =======*

		// the number of samples
		const R_xlen_t n = MCWorkingGeno.Space.SampleNum();
		// to detect the block size
		AutoDetectSNPBlockSize(n);

		// the upper-triangle IBD matrix
		CdMatTri<TS_KINGHomo> IBD(n);
		// Calculate the IBD matrix
		DoKINGCalculate(IBD, INTEGER(NumThread)[0], "KING IBD:", verbose);

		// initialize output variables
		SEXP dim;
		PROTECT(dim = NEW_INTEGER(2));
		INTEGER(dim)[0] = n; INTEGER(dim)[1] = n;

		SEXP k0    = NEW_NUMERIC(n*n);  PROTECT(k0);
		setAttrib(k0, R_DimSymbol, dim);
		double *out_k0 = REAL(k0);

		SEXP k1    = NEW_NUMERIC(n*n);  PROTECT(k1);
		setAttrib(k1, R_DimSymbol, dim);
		double *out_k1 = REAL(k1);

		// output
		TS_KINGHomo *p = IBD.get();
		for (int i=0; i < n; i++)
		{
			out_k0[i*n + i] = out_k1[i*n + i] = 0;
			p ++;
			for (int j=i+1; j < n; j++, p++)
			{
				double theta = 0.5 - p->SumSq / (8 * p->SumAFreq);
				double k0 = p->IBS0 / (2 * p->SumAFreq2);
				double k1 = 2 - 2*k0 - 4*theta;
				out_k0[i*n + j] = out_k0[j*n + i] = k0;
				out_k1[i*n + j] = out_k1[j*n + i] = k1;
			}
		}

		// output
		PROTECT(rv_ans = NEW_LIST(2));
		SET_ELEMENT(rv_ans, 0, k0);
		SET_ELEMENT(rv_ans, 1, k1);
		UNPROTECT(4);

	COREARRAY_CATCH
}

/// to compute the IBD coefficients by KING method of moment (KING-robust)
COREARRAY_DLL_EXPORT SEXP gnrIBD_KING_Robust(SEXP FamilyID, SEXP NumThread,
	SEXP _Verbose)
{
	bool verbose = SEXP_Verbose(_Verbose);

	COREARRAY_TRY

		// =======* To cache the genotype data =======*
		CachingSNPData("KING IBD", verbose);

		// =======* The calculation of genetic covariance matrix =======*

		// the number of samples
		const R_xlen_t n = MCWorkingGeno.Space.SampleNum();
		// to detect the block size
		AutoDetectSNPBlockSize(n);

		// IBD KING robust IBD estimator across family
		// the upper-triangle IBD matrix
		CdMatTri<TS_KINGRobust> IBD(n);
		// Calculate the IBD matrix
		DoKINGCalculate(IBD, INTEGER(NumThread)[0], "KING IBD:", verbose);

		// initialize output variables
		SEXP dim;
		PROTECT(dim = NEW_INTEGER(2));
		INTEGER(dim)[0] = n; INTEGER(dim)[1] = n;

		SEXP IBS0    = NEW_NUMERIC(n*n);  PROTECT(IBS0);
		setAttrib(IBS0, R_DimSymbol, dim);
		double *out_IBS0 = REAL(IBS0);

		SEXP Kinship = NEW_NUMERIC(n*n);  PROTECT(Kinship);
		setAttrib(Kinship, R_DimSymbol, dim);
		double *out_kinship = REAL(Kinship);

		int *FamID_Ptr = INTEGER(FamilyID);

		// output
		TS_KINGRobust *p = IBD.get();
		for (int i=0; i < n; i++)
		{
			out_IBS0[i*n + i] = 0; out_kinship[i*n + i] = 0.5;
			p ++;
			for (int j=i+1; j < n; j++, p++)
			{
				out_IBS0[i*n + j] = out_IBS0[j*n + i] = double(p->IBS0) / p->nLoci;

				int f1 = FamID_Ptr[i];
				int f2 = FamID_Ptr[j];
				out_kinship[i*n + j] = out_kinship[j*n + i] =
						((f1 == f2) && (f1 != NA_INTEGER)) ?
					(0.5 - p->SumSq / (2.0 *(p->N1_Aa + p->N2_Aa))) :
					(0.5 - p->SumSq / (4.0 * min(p->N1_Aa, p->N2_Aa)));
			}
		}

		// output
		PROTECT(rv_ans = NEW_LIST(2));
		SET_ELEMENT(rv_ans, 0, IBS0);
		SET_ELEMENT(rv_ans, 1, Kinship);
		UNPROTECT(4);

	COREARRAY_CATCH
}
}

#endif  /* _H_IBD_KING_ */
