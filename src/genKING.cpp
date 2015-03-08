// ===========================================================
//
// genKING.cpp: KINK Identity-by-descent (IBD) Analysis on GWAS
//
// Copyright (C) 2011-2015    Xiuwen Zheng
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


#ifndef _HEADER_IBD_KING_
#define _HEADER_IBD_KING_

// CoreArray library header
#include <dGenGWAS.h>
#include <dVect.h>

// Standard library header
#include <cmath>
#include <cfloat>
#include <memory>
#include <algorithm>


namespace KING_IBD
{
	// using namespace
	using namespace std;
	using namespace GWAS;


	/// the size of a lookup table
	#define SIZE_LOOKUP_TABLE    (256*256)

	/// The number of IBS 0 in the packed genotype
	static C_UInt8 LT_nIBS0[SIZE_LOOKUP_TABLE];

	/// The square value of genotype difference, (X_m^{(i)} - X_m^{(j)})^2
	static C_UInt8 LT_SqDiff[SIZE_LOOKUP_TABLE];

	/// N1_Aa requiring both genotypes are available
	static C_UInt8 LT_N1_Aa[SIZE_LOOKUP_TABLE];

	/// N2_Aa requiring both genotypes are available
	static C_UInt8 LT_N2_Aa[SIZE_LOOKUP_TABLE];

	/// the number of valid loci
	static C_UInt8 LT_nLoci[SIZE_LOOKUP_TABLE];

	/// The flag of use of allele frequencies
	static C_UInt8 LT_ValidFlag[SIZE_LOOKUP_TABLE];

	/// It is used in 'DetectOptimizedNumOfSNP'
	static const size_t LT_MEM_SIZE =
		sizeof(LT_nIBS0) + sizeof(LT_SqDiff) + sizeof(LT_N1_Aa) +
		sizeof(LT_N2_Aa) + sizeof(LT_nLoci) + sizeof(LT_ValidFlag);


	#ifdef __GNUC__
	#   define COREARRAY_PACKED    __attribute__((packed))
	#else
	#   define COREARRAY_PACKED
	#   ifdef __IBMC__
	#       pragma pack(1)
	#   else
	#       pragma pack(push, 1)
	#   endif
	#endif

	/// The structure of KING IBD Homo estimator
	struct COREARRAY_PACKED TS_KINGHomo
	{
		C_UInt32 IBS0;     //< the number of loci sharing no allele
		C_UInt32 SumSq;    //< \sum_m (X_m^{(i)} - X_m^{(j)})^2
		double SumAFreq;   //< \sum_m p_m (1 - p_m)
		double SumAFreq2;  //< \sum_m p_m^2 (1 - p_m)^2
		TS_KINGHomo()
			{ IBS0 = SumSq = 0; SumAFreq = SumAFreq2 = 0; }
	};

	/// The structure of KING IBD Robust Estimator
	struct COREARRAY_PACKED TS_KINGRobust
	{
		C_UInt32 IBS0;   //< the number of loci sharing no allele
		C_UInt32 nLoci;  //< the total number of loci
		C_UInt32 SumSq;  //< \sum_m (X_m^{(i)} - X_m^{(j)})^2
		C_UInt32 N1_Aa;  //< the number of hetet loci for the first individual
		C_UInt32 N2_Aa;  //< the number of hetet loci for the second individual
		TS_KINGRobust()
			{ IBS0 = nLoci = SumSq = N1_Aa = N2_Aa = 0; }
	};

	#ifndef __GNUC__
	#   pragma pack(pop)
	#endif



	// TInit object
	class COREARRAY_DLL_LOCAL TInit
	{
	public:
		TInit()
		{
			#define LOOKUP_COND(cond, var, op)	\
				for (int s=0; s < SIZE_LOOKUP_TABLE; s++)	\
				{	\
					int g1 = s >> 8, g2 = s & 0xFF;	\
					int sum = 0;	\
					for (int i=0; i < 4; i++)	\
					{	\
						int b1 = g1 & 0x03, b2 = g2 & 0x03;	\
						if (cond) op;	\
						g1 >>= 2; g2 >>= 2;	\
					}	\
					var[s] = sum;	\
				}

			LOOKUP_COND((b1 < 3) && (b2 < 3) && (abs(b1-b2)==2), LT_nIBS0, sum++);
			LOOKUP_COND((b1 < 3) && (b2 < 3), LT_SqDiff, sum += (b1-b2)*(b1-b2));
			LOOKUP_COND((b1 < 3) && (b2 < 3), LT_N1_Aa, sum += (b1==1) ? 1:0);
			LOOKUP_COND((b1 < 3) && (b2 < 3), LT_N2_Aa, sum += (b2==1) ? 1:0);
			LOOKUP_COND((b1 < 3) && (b2 < 3), LT_nLoci, sum ++);
			LOOKUP_COND((b1 < 3) && (b2 < 3), LT_ValidFlag, sum |= (1 << i));
		}
	} InitObj;



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
		C_UInt8 *pPack = &Array_PackedGeno[0];

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
			double &Freq = Array_AlleleFreq[iSNP];
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
		long Cnt = Array_Thread_MatCnt[ThreadIndex];
		IdMatTri I = Array_Thread_MatIdx[ThreadIndex];
		TS_KINGHomo *p = ((TS_KINGHomo*)Param) + I.Offset();
		long _PackSNPLen = (SNP_Cnt / 4) + (SNP_Cnt % 4 ? 1 : 0);

		for (; Cnt > 0; Cnt--, ++I, p++)
		{
			C_UInt8 *p1 = &Array_PackedGeno[0] + I.Row()*_PackSNPLen;
			C_UInt8 *p2 = &Array_PackedGeno[0] + I.Column()*_PackSNPLen;
			for (long k=0; k < _PackSNPLen; k++, p1++, p2++)
			{
				size_t t = (size_t(*p1) << 8) | (*p2);

				p->IBS0  += LT_nIBS0[t];
				p->SumSq += LT_SqDiff[t];

				C_UInt8 flag = LT_ValidFlag[t];
				if (flag & 0x01)
				{
					double f = Array_AlleleFreq[4*k + 0];
					p->SumAFreq += f; p->SumAFreq2 += f*f;
				}
				if (flag & 0x02)
				{
					double f = Array_AlleleFreq[4*k + 1];
					p->SumAFreq += f; p->SumAFreq2 += f*f;
				}
				if (flag & 0x04)
				{
					double f = Array_AlleleFreq[4*k + 2];
					p->SumAFreq += f; p->SumAFreq2 += f*f;
				}
				if (flag & 0x08)
				{
					double f = Array_AlleleFreq[4*k + 3];
					p->SumAFreq += f; p->SumAFreq2 += f*f;
				}
			}
		}
	}

	/// Compute IBD estimator in KING-robust
	static void _Do_KING_Robust_Compute(int ThreadIndex, long Start,
		long SNP_Cnt, void* Param)
	{
		long Cnt = Array_Thread_MatCnt[ThreadIndex];
		IdMatTri I = Array_Thread_MatIdx[ThreadIndex];
		TS_KINGRobust *p = ((TS_KINGRobust*)Param) + I.Offset();
		long _PackSNPLen = (SNP_Cnt / 4) + (SNP_Cnt % 4 ? 1 : 0);

		for (; Cnt > 0; Cnt--, ++I, p++)
		{
			C_UInt8 *p1 = &Array_PackedGeno[0] + I.Row()*_PackSNPLen;
			C_UInt8 *p2 = &Array_PackedGeno[0] + I.Column()*_PackSNPLen;
			for (long k=0; k < _PackSNPLen; k++, p1++, p2++)
			{
				size_t t = (size_t(*p1) << 8) | (*p2);

				p->IBS0  += LT_nIBS0[t];
				p->nLoci += LT_nLoci[t];
				p->SumSq += LT_SqDiff[t];
				p->N1_Aa += LT_N1_Aa[t];
				p->N2_Aa += LT_N2_Aa[t];
			}
		}
	}

	/// Calculate KING IBD Homo Estimator
	COREARRAY_DLL_LOCAL void DoKINGCalculate(CdMatTri<TS_KINGHomo> &PublicKING,
		int NumThread, const char *Info, bool verbose)
	{
		// Initialize ...
		Array_PackedGeno.resize(BlockNumSNP * PublicKING.N());
		memset(PublicKING.get(), 0, sizeof(TS_KINGHomo)*PublicKING.Size());
		Array_AlleleFreq.resize(BlockNumSNP);

		MCWorkingGeno.Progress.Info = Info;
		MCWorkingGeno.Progress.Show() = verbose;
		MCWorkingGeno.InitParam(true, true, BlockNumSNP);

		MCWorkingGeno.SplitJobs(NumThread, PublicKING.N(),
			Array_Thread_MatIdx, Array_Thread_MatCnt);
		MCWorkingGeno.Run(NumThread, &_Do_KING_ReadBlock,
			&_Do_KING_Homo_Compute, PublicKING.get());
	}

	/// Calculate KING IBD Robust Estimator
	COREARRAY_DLL_LOCAL void DoKINGCalculate(CdMatTri<TS_KINGRobust> &PublicKING,
		int NumThread, const char *Info, bool verbose)
	{
		// Initialize ...
		Array_PackedGeno.resize(BlockNumSNP * PublicKING.N());
		memset(PublicKING.get(), 0, sizeof(TS_KINGRobust)*PublicKING.Size());
		Array_AlleleFreq.resize(BlockNumSNP);

		MCWorkingGeno.Progress.Info = Info;
		MCWorkingGeno.Progress.Show() = verbose;
		MCWorkingGeno.InitParam(true, true, BlockNumSNP);

		MCWorkingGeno.SplitJobs(NumThread, PublicKING.N(),
			Array_Thread_MatIdx, Array_Thread_MatCnt);
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

		// ======= To cache the genotype data =======
		CachingSNPData("KING IBD", verbose);

		// ======= The calculation of genetic covariance matrix =======

		// the number of samples
		const R_xlen_t n = MCWorkingGeno.Space.SampleNum();

		// to detect the block size
		DetectOptimizedNumOfSNP(n, LT_MEM_SIZE);

		// the upper-triangle IBD matrix
		CdMatTri<TS_KINGHomo> IBD(n);
		// Calculate the IBD matrix
		DoKINGCalculate(IBD, INTEGER(NumThread)[0], "KING IBD:", verbose);

		// initialize output variables
		SEXP K0 = PROTECT(allocMatrix(REALSXP, n, n));
		SEXP K1 = PROTECT(allocMatrix(REALSXP, n, n));

		double *pK0 = REAL(K0);
		double *pK1 = REAL(K1);

		// output
		TS_KINGHomo *p = IBD.get();
		for (R_xlen_t i=0; i < n; i++)
		{
			pK0[i*n + i] = pK1[i*n + i] = 0;
			p ++;
			for (R_xlen_t j=i+1; j < n; j++, p++)
			{
				double theta = 0.5 - p->SumSq / (8 * p->SumAFreq);
				double k0 = p->IBS0 / (2 * p->SumAFreq2);
				double k1 = 2 - 2*k0 - 4*theta;
				pK0[i*n + j] = pK0[j*n + i] = k0;
				pK1[i*n + j] = pK1[j*n + i] = k1;
			}
		}

		// output
		PROTECT(rv_ans = NEW_LIST(2));
		SET_ELEMENT(rv_ans, 0, K0);
		SET_ELEMENT(rv_ans, 1, K1);
		UNPROTECT(3);

	COREARRAY_CATCH
}

/// to compute the IBD coefficients by KING method of moment (KING-robust)
COREARRAY_DLL_EXPORT SEXP gnrIBD_KING_Robust(SEXP FamilyID, SEXP NumThread,
	SEXP _Verbose)
{
	bool verbose = SEXP_Verbose(_Verbose);

	COREARRAY_TRY

		// ======= To cache the genotype data =======
		CachingSNPData("KING IBD", verbose);

		// ======= The calculation of genetic covariance matrix =======

		// the number of samples
		const R_xlen_t n = MCWorkingGeno.Space.SampleNum();

		// to detect the block size
		DetectOptimizedNumOfSNP(n, LT_MEM_SIZE);

		// IBD KING robust IBD estimator across family

		// the upper-triangle IBD matrix
		CdMatTri<TS_KINGRobust> IBD(n);
		// Calculate the IBD matrix
		DoKINGCalculate(IBD, INTEGER(NumThread)[0], "KING IBD:", verbose);

		// initialize output variables
		SEXP IBS0 = PROTECT(allocMatrix(REALSXP, n, n));
		SEXP Kinship = PROTECT(allocMatrix(REALSXP, n, n));

		double *pIBS0 = REAL(IBS0);
		double *pKinship = REAL(Kinship);
		int *FamID_Ptr = INTEGER(FamilyID);

		// output
		TS_KINGRobust *p = IBD.get();
		for (R_xlen_t i=0; i < n; i++)
		{
			pIBS0[i*n + i] = 0; pKinship[i*n + i] = 0.5;
			p ++;
			for (R_xlen_t j=i+1; j < n; j++, p++)
			{
				pIBS0[i*n + j] = pIBS0[j*n + i] = double(p->IBS0) / p->nLoci;

				int f1 = FamID_Ptr[i];
				int f2 = FamID_Ptr[j];
				pKinship[i*n + j] = pKinship[j*n + i] =
					((f1 == f2) && (f1 != NA_INTEGER)) ?
					(0.5 - p->SumSq / (2.0 *(p->N1_Aa + p->N2_Aa))) :
					(0.5 - p->SumSq / (4.0 * min(p->N1_Aa, p->N2_Aa)));
			}
		}

		// output
		PROTECT(rv_ans = NEW_LIST(2));
		SET_ELEMENT(rv_ans, 0, IBS0);
		SET_ELEMENT(rv_ans, 1, Kinship);
		UNPROTECT(3);

	COREARRAY_CATCH
}

}

#endif  /* _HEADER_IBD_KING_ */
