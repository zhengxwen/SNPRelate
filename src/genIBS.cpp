// ===========================================================
//     _/_/_/   _/_/_/  _/_/_/_/    _/_/_/_/  _/_/_/   _/_/_/
//      _/    _/       _/             _/    _/    _/   _/   _/
//     _/    _/       _/_/_/_/       _/    _/    _/   _/_/_/
//    _/    _/       _/             _/    _/    _/   _/
// _/_/_/   _/_/_/  _/_/_/_/_/     _/     _/_/_/   _/_/
// ===========================================================
//
// genIBS.cpp: Identity by state (IBS) analysis on GWAS
//
// Copyright (C) 2011 - 2014	Xiuwen Zheng [zhengxwen@gmail.com]
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


namespace IBS
{
	// using namespace
	using namespace std;
	using namespace CoreArray;
	using namespace CoreArray::Vectorization;
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


	//================    KING robust estimator    ================//
	/// The square value of genotype difference, (X_m^{(i)} - X_m^{(j)})^2
	C_UInt8 Gen_KING_SqDiff[PACKED_SIZE];
	/// N1_Aa requiring both genotypes are available
	C_UInt8 Gen_KING_N1_Aa[PACKED_SIZE];
	/// N2_Aa requiring both genotypes are available
	C_UInt8 Gen_KING_N2_Aa[PACKED_SIZE];
	/// the number of valid loci
	C_UInt8 Gen_KING_Num_Loci[PACKED_SIZE];


	/// The packed genotype buffer
	vector<C_UInt8> GenoPacked;
	/// The allele frequencies
	vector<double> GenoAlleleFreq;


	/// Thread variables
	const int N_MAX_THREAD = 256;
	// PLINK -- IBS
	IdMatTriD PLINKIBS_Thread_MatIdx[N_MAX_THREAD];
	C_Int64 PLINKIBS_Thread_MatCnt[N_MAX_THREAD];

	// IBS, KING IBD, Individual Similarity
	IdMatTri IBS_Thread_MatIdx[N_MAX_THREAD];
	C_Int64 IBS_Thread_MatCnt[N_MAX_THREAD];



	/// The pointer to the variable 'PublicIBS' in the function "DoIBSCalculate"
	/// The structure of IBS states
	struct TS_IBS
	{
		C_UInt32 IBS0;  //< the number of loci sharing no allele
		C_UInt32 IBS1;  //< the number of loci sharing only one allele
		C_UInt32 IBS2;  //< the number of loci sharing two alleles
		TS_IBS() { IBS0 = IBS1 = IBS2 = 0; }
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
	void AutoDetectSNPBlockSize(int nSamp, bool Detect=true)
	{
		if (Detect)
		{
			C_UInt64 L2Cache = GDS_Mach_GetCPULevelCache(2);
			C_UInt64 L3Cache = GDS_Mach_GetCPULevelCache(3);
			C_UInt64 Cache = (L2Cache > L3Cache) ? L2Cache : L3Cache;
			if ((C_Int64)Cache <= 0) Cache = 1024*1024; // 1M
			BlockNumSNP = (Cache - 3*256*256 - 8*1024) / nSamp * 4;
		}
		BlockNumSNP = (BlockNumSNP / 4) * 4;
		if (BlockNumSNP < 16) BlockNumSNP = 16;
	}

	/// Convert the raw genotypes
	static void _Do_IBS_ReadBlock(C_UInt8 *GenoBuf, long Start,
		long SNP_Cnt, void* Param)
	{
		// initialize
		const int nSamp = MCWorkingGeno.Space.SampleNum();
		C_UInt8 *pG = GenoBuf;
		C_UInt8 *pPack = &GenoPacked[0];

		// pack genotypes
		for (long iSamp=0; iSamp < nSamp; iSamp++)
		{
			pPack = PackGeno2b(pG, SNP_Cnt, pPack);
			pG += SNP_Cnt;
		}
	}

	/// Compute the pairwise IBS matrix for PLINK
	static void _Do_PLINKIBS_Compute(int ThreadIndex, long Start,
		long SNP_Cnt, void* Param)
	{
		long Cnt = PLINKIBS_Thread_MatCnt[ThreadIndex];
		IdMatTriD I = PLINKIBS_Thread_MatIdx[ThreadIndex];
		TS_IBS *p = ((TS_IBS*)Param) + I.Offset();
		long _PackSNPLen = (SNP_Cnt / 4) + (SNP_Cnt % 4 ? 1 : 0);

		for (; Cnt > 0; Cnt--, ++I, p++)
		{
			C_UInt8 *p1 = &GenoPacked[0] + I.Row()*_PackSNPLen;
			C_UInt8 *p2 = &GenoPacked[0] + I.Column()*_PackSNPLen;
			for (long k=_PackSNPLen; k > 0; k--, p1++, p2++)
			{
				size_t t = (size_t(*p1) << 8) | (*p2);
				p->IBS0 += IBS0_Num_SNP[t];
				p->IBS1 += IBS1_Num_SNP[t];
				p->IBS2 += IBS2_Num_SNP[t];
			}
		}
	}

	/// Compute the pairwise IBS matrix
	static void _Do_IBS_Compute(int ThreadIndex, long Start,
		long SNP_Cnt, void* Param)
	{
		long Cnt = IBS_Thread_MatCnt[ThreadIndex];
		IdMatTri I = IBS_Thread_MatIdx[ThreadIndex];
		TS_IBS *p = ((TS_IBS*)Param) + I.Offset();
		long _PackSNPLen = (SNP_Cnt / 4) + (SNP_Cnt % 4 ? 1 : 0);

		for (; Cnt > 0; Cnt--, ++I, p++)
		{
			C_UInt8 *p1 = &GenoPacked[0] + I.Row()*_PackSNPLen;
			C_UInt8 *p2 = &GenoPacked[0] + I.Column()*_PackSNPLen;
			for (long k=_PackSNPLen; k > 0; k--, p1++, p2++)
			{
				size_t t = (size_t(*p1) << 8) | (*p2);
				p->IBS0 += IBS0_Num_SNP[t];
				p->IBS1 += IBS1_Num_SNP[t];
				p->IBS2 += IBS2_Num_SNP[t];
			}
		}
	}

	/// Calculate the IBS matrix for PLINK
	void DoPLINKIBSCalculate(CdMatTriDiag<TS_IBS> &PublicIBS, int NumThread,
		const char *Info, bool verbose)
	{
		// Initialize ...
		GenoPacked.resize(BlockNumSNP * PublicIBS.N());
		memset(PublicIBS.get(), 0, sizeof(TS_IBS)*PublicIBS.Size());

		MCWorkingGeno.Progress.Info = Info;
		MCWorkingGeno.Progress.Show() = verbose;
		MCWorkingGeno.InitParam(true, true, BlockNumSNP);

		MCWorkingGeno.SplitJobs(NumThread, PublicIBS.N(), PLINKIBS_Thread_MatIdx, PLINKIBS_Thread_MatCnt);
		MCWorkingGeno.Run(NumThread, &_Do_IBS_ReadBlock, &_Do_PLINKIBS_Compute, PublicIBS.get());
	}

	/// Calculate the IBS matrix
	void DoIBSCalculate(CdMatTri<TS_IBS> &PublicIBS, int NumThread,
		const char *Info, bool verbose)
	{
		// Initialize ...
		GenoPacked.resize(BlockNumSNP * PublicIBS.N());
		memset(PublicIBS.get(), 0, sizeof(TS_IBS)*PublicIBS.Size());

		MCWorkingGeno.Progress.Info = Info;
		MCWorkingGeno.Progress.Show() = verbose;
		MCWorkingGeno.InitParam(true, true, BlockNumSNP);

		MCWorkingGeno.SplitJobs(NumThread, PublicIBS.N(), IBS_Thread_MatIdx, IBS_Thread_MatCnt);
		MCWorkingGeno.Run(NumThread, &_Do_IBS_ReadBlock, &_Do_IBS_Compute, PublicIBS.get());
	}



	/// ======================================================================
	/// Individual Dissimilarity
	/// ======================================================================

	/// Convert the raw genotypes
	static void _Do_Diss_ReadBlock(C_UInt8 *GenoBuf, long Start,
		long SNP_Cnt, void* Param)
	{
		// initialize
		const int nSamp = MCWorkingGeno.Space.SampleNum();
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
		long Cnt = IBS_Thread_MatCnt[ThreadIndex];
		IdMatTri I = IBS_Thread_MatIdx[ThreadIndex];
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
		memset(PublicDiss.get(), 0, sizeof(TS_Dissimilarity)*PublicDiss.Size());
		GenoAlleleFreq.resize(BlockNumSNP);

		MCWorkingGeno.Progress.Info = Info;
		MCWorkingGeno.Progress.Show() = verbose;
		MCWorkingGeno.InitParam(true, true, BlockNumSNP);

		MCWorkingGeno.SplitJobs(NumThread, PublicDiss.N(),
			IBS_Thread_MatIdx, IBS_Thread_MatCnt);
		MCWorkingGeno.Run(NumThread, &_Do_Diss_ReadBlock,
			&_Do_Diss_Compute, PublicDiss.get());
	}
}

namespace IBD
{
	// PLINK method of moment
	void Init_EPrIBD_IBS(const double in_afreq[], double out_afreq[],
		bool CorrectFactor, long nSNP = -1);

	void Est_PLINK_Kinship(int IBS0, int IBS1, int IBS2,
		double &k0, double &k1, bool KinshipConstrict);
}


using namespace IBS;

extern "C"
{
// ======================================================================*
// the functions for identity-by-state (IBS)
//

// the functions for Identity by state (IBS)

/// to compute the average IBS
COREARRAY_DLL_EXPORT SEXP gnrIBSAve(SEXP NumThread, SEXP _Verbose)
{
	bool verbose = SEXP_Verbose(_Verbose);

	COREARRAY_TRY

		// =======* To cache the genotype data =======*
		CachingSNPData("IBS", verbose);

		// =======* The calculation of IBS matrix =======*

		// the number of samples
		const R_xlen_t n = MCWorkingGeno.Space.SampleNum();
		// to detect the block size
		IBS::AutoDetectSNPBlockSize(n);
		// the upper-triangle genetic covariance matrix
		CdMatTri<IBS::TS_IBS> IBS(n);

		// initialize output variables
		PROTECT(rv_ans = NEW_NUMERIC(n*n));
		double *out_IBSMat = REAL(rv_ans);

		SEXP dim;
		PROTECT(dim = NEW_INTEGER(2));
		INTEGER(dim)[0] = n; INTEGER(dim)[1] = n;
		setAttrib(rv_ans, R_DimSymbol, dim);

		// Calculate the IBS matrix
		IBS::DoIBSCalculate(IBS, INTEGER(NumThread)[0], "IBS:", verbose);

		// output
		IBS::TS_IBS *p = IBS.get();
		for (int i=0; i < n; i++)
		{
			for (int j=i; j < n; j++, p++)
			{
				out_IBSMat[i*n + j] = out_IBSMat[j*n + i] =
					double(0.5*p->IBS1 + p->IBS2) /
					(p->IBS0 + p->IBS1 + p->IBS2);
			}
		}

		UNPROTECT(2);

	COREARRAY_CATCH
}

/// to compute the average IBS
COREARRAY_DLL_EXPORT SEXP gnrIBSNum(SEXP NumThread, SEXP _Verbose)
{
	bool verbose = SEXP_Verbose(_Verbose);

	COREARRAY_TRY

		// =======* To cache the genotype data =======*
		CachingSNPData("IBS", verbose);

		// =======* The calculation of IBS matrix =======*

		// the number of samples
		const R_xlen_t n = MCWorkingGeno.Space.SampleNum();
		// to detect the block size
		IBS::AutoDetectSNPBlockSize(n);

		// the upper-triangle IBS matrix
		CdMatTri<IBS::TS_IBS> IBS(n);

		// Calculate the IBS matrix
		IBS::DoIBSCalculate(IBS, INTEGER(NumThread)[0], "IBS:", verbose);

		// initialize output variables
		SEXP dim;
		PROTECT(dim = NEW_INTEGER(2));
		INTEGER(dim)[0] = n; INTEGER(dim)[1] = n;

		SEXP IBS0    = NEW_NUMERIC(n*n);  PROTECT(IBS0);
		setAttrib(IBS0, R_DimSymbol, dim);
		SEXP IBS1    = NEW_NUMERIC(n*n);  PROTECT(IBS1);
		setAttrib(IBS1, R_DimSymbol, dim);
		SEXP IBS2    = NEW_NUMERIC(n*n);  PROTECT(IBS2);
		setAttrib(IBS2, R_DimSymbol, dim);

		double *out_IBS0    = REAL(IBS0);
		double *out_IBS1    = REAL(IBS1);
		double *out_IBS2    = REAL(IBS2);

		// output
		IBS::TS_IBS *p = IBS.get();
		for (int i=0; i < n; i++)
		{
			for (int j=i; j < n; j++, p++)
			{
				out_IBS0[i*n + j] = out_IBS0[j*n + i] = p->IBS0;
				out_IBS1[i*n + j] = out_IBS1[j*n + i] = p->IBS1;
				out_IBS2[i*n + j] = out_IBS2[j*n + i] = p->IBS2;
			}
		}

		// output
		PROTECT(rv_ans = NEW_LIST(3));
		SET_ELEMENT(rv_ans, 0, IBS0);
		SET_ELEMENT(rv_ans, 1, IBS1);
		SET_ELEMENT(rv_ans, 2, IBS2);
		UNPROTECT(5);

	COREARRAY_CATCH
}


// ======================================================================*
// the functions for identity-by-descent (IBD)
//

/// to compute the IBD coefficients by PLINK method of moment
COREARRAY_DLL_EXPORT SEXP gnrIBD_PLINK(SEXP NumThread, SEXP AlleleFreq,
	SEXP UseSpecificAFreq, SEXP KinshipConstrict, SEXP _Verbose)
{
	bool verbose = SEXP_Verbose(_Verbose);

	COREARRAY_TRY

		// =======* To cache the genotype data =======*
		CachingSNPData("PLINK IBD", verbose);

		// the number of individuals
		const R_xlen_t n = MCWorkingGeno.Space.SampleNum();
		// the number of SNPs
		const R_xlen_t m = MCWorkingGeno.Space.SNPNum();
		// to detect the block size
		IBS::AutoDetectSNPBlockSize(n);

		// =======* PLINK method of moment =======*

		// initialize output variables
		SEXP dim;
		PROTECT(dim = NEW_INTEGER(2));
		INTEGER(dim)[0] = n; INTEGER(dim)[1] = n;

		SEXP k0    = NEW_NUMERIC(n*n);  PROTECT(k0);
		setAttrib(k0, R_DimSymbol, dim);

		SEXP k1    = NEW_NUMERIC(n*n);  PROTECT(k1);
		setAttrib(k1, R_DimSymbol, dim);

		SEXP afreq = NEW_NUMERIC(m);    PROTECT(afreq);
		double *out_k0    = REAL(k0);
		double *out_k1    = REAL(k1);
		double *out_afreq = REAL(afreq);

		// initialize the internal matrix
		IBD::Init_EPrIBD_IBS(LOGICAL(UseSpecificAFreq)[0] ? REAL(AlleleFreq) : NULL,
			out_afreq, !(LOGICAL(UseSpecificAFreq)[0]));
		// the upper-triangle genetic covariance matrix
		CdMatTriDiag<IBS::TS_IBS> IBS(IBS::TS_IBS(), n);
		// Calculate the IBS matrix
		IBS::DoPLINKIBSCalculate(IBS, INTEGER(NumThread)[0],
			"PLINK IBD:", verbose);

		// output
		bool kc = LOGICAL(KinshipConstrict)[0] == TRUE;
		IBS::TS_IBS *p = IBS.get();
		for (int i=0; i < n; i++)
		{
			out_k0[i*n + i] = out_k1[i*n + i] = 0;
			for (int j=i+1; j < n; j++, p++)
			{
				double k0, k1;
				IBD::Est_PLINK_Kinship(p->IBS0, p->IBS1, p->IBS2, k0, k1, kc);
				out_k0[i*n + j] = out_k0[j*n + i] = k0;
				out_k1[i*n + j] = out_k1[j*n + i] = k1;
			}
		}

		// output
		PROTECT(rv_ans = NEW_LIST(3));
		SET_ELEMENT(rv_ans, 0, k0);
		SET_ELEMENT(rv_ans, 1, k1);
		SET_ELEMENT(rv_ans, 2, afreq);
		UNPROTECT(5);

	COREARRAY_CATCH
}



// =========================================================================
// the functions for individual dissimilarity

/// to compute the individual dissimilarity
COREARRAY_DLL_EXPORT SEXP gnrDiss(SEXP NumThread, SEXP _Verbose)
{
	bool verbose = SEXP_Verbose(_Verbose);

	COREARRAY_TRY

		// =======* To cache the genotype data =======*
		CachingSNPData("Dissimilarity", verbose);

		// =======* The calculation of genetic covariance matrix =======*

		// the number of samples
		const R_xlen_t n = MCWorkingGeno.Space.SampleNum();
		// to detect the block size
		IBS::AutoDetectSNPBlockSize(n);
		// the upper-triangle genetic covariance matrix
		CdMatTri<IBS::TS_Dissimilarity> Dist(n);

		// Calculate the genetic distance matrix
		IBS::DoDissCalculate(Dist, INTEGER(NumThread)[0], "Dissimilarity:",
			verbose);

		// output
		PROTECT(rv_ans = allocMatrix(REALSXP, n, n));

		IBS::TS_Dissimilarity *p = Dist.get();
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

}

#endif  /* _HEADER_IBS_ */
