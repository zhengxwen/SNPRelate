// ===========================================================
//     _/_/_/   _/_/_/  _/_/_/_/    _/_/_/_/  _/_/_/   _/_/_/
//      _/    _/       _/             _/    _/    _/   _/   _/
//     _/    _/       _/_/_/_/       _/    _/    _/   _/_/_/
//    _/    _/       _/             _/    _/    _/   _/
// _/_/_/   _/_/_/  _/_/_/_/_/     _/     _/_/_/   _/_/
// ===========================================================
//
// genIBS.cpp: Identity by state (IBS) analysis on genome-wide association studies
//
// Copyright (C) 2013	Xiuwen Zheng
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


// CoreArray library header
#include <dType.h>
#include <dVect.h>
#include <CoreGDSLink.h>
#include <dGenGWAS.h>

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


#ifndef _FuncIBS_H_
#define _FuncIBS_H_

namespace IBS
{
	// using namespace
	using namespace std;
	using namespace CoreArray;
	using namespace CoreArray::Vectorization;
	using namespace GWAS;


	/// Packed size
	static const long _SIZE_ = 256*256;

	/// IBS
	/// The number of IBS 0 in the packed genotype
	UInt8 IBS0_Num_SNP[_SIZE_];
	/// The number of IBS 1 in the packed genotype
	UInt8 IBS1_Num_SNP[_SIZE_];
	/// The number of IBS 2 in the packed genotype
	UInt8 IBS2_Num_SNP[_SIZE_];

	/// Genetic Distance
	/// The distance in the packed genotype
	UInt8 Gen_Dist_SNP[_SIZE_];
	/// The flag of use of allele frequencies
	UInt8 Gen_Both_Valid[_SIZE_];

	/// KING robust estimator
	/// The square value of genotype difference, (X_m^{(i)} - X_m^{(j)})^2
	UInt8 Gen_KING_SqDiff[_SIZE_];
	/// N1_Aa requiring both genotypes are available
	UInt8 Gen_KING_N1_Aa[_SIZE_];
	/// N2_Aa requiring both genotypes are available
	UInt8 Gen_KING_N2_Aa[_SIZE_];
	/// the number of valid loci
	UInt8 Gen_KING_Num_Loci[_SIZE_];


	/// The packed genotype buffer
	auto_ptr<UInt8> GenoPacked;
	/// The allele frequencies
	auto_ptr<double> GenoAlleleFreq;


	/// Thread variables
	const int N_MAX_THREAD = 256;
	// PLINK -- IBS
	IdMatTriD PLINKIBS_Thread_MatIdx[N_MAX_THREAD];
	Int64 PLINKIBS_Thread_MatCnt[N_MAX_THREAD];

	// IBS, KING IBD, Individual Similarity
	IdMatTri IBS_Thread_MatIdx[N_MAX_THREAD];
	Int64 IBS_Thread_MatCnt[N_MAX_THREAD];



	/// The pointer to the variable 'PublicIBS' in the function "DoIBSCalculate"
	/// The structure of IBS states
	struct TIBS_Flag
	{
		UInt32 IBS0;  //< the number of loci sharing no allele
		UInt32 IBS1;  //< the number of loci sharing only one allele
		UInt32 IBS2;  //< the number of loci sharing two alleles
		TIBS_Flag() { IBS0 = IBS1 = IBS2 = 0; }
	};


	/// The pointer to the variable 'PublicKING' in the function "DoKINGCalculate"
	/// The structure of KING IBD estimator
	struct TKINGHomoFlag
	{
		UInt32 IBS0;       //< the number of loci sharing no allele
		UInt32 SumSq;      //< \sum_m (X_m^{(i)} - X_m^{(j)})^2
		double SumAFreq;   //< \sum_m p_m (1 - p_m)
		double SumAFreq2;  //< \sum_m p_m^2 (1 - p_m)^2
		TKINGHomoFlag() { IBS0 = SumSq = 0; SumAFreq = SumAFreq2 = 0; }
	};

	struct TKINGRobustFlag
	{
		UInt32 IBS0;       //< the number of loci sharing no allele
		UInt32 nLoci;      //< the total number of loci
		UInt32 SumSq;      //< \sum_m (X_m^{(i)} - X_m^{(j)})^2
		UInt32 N1_Aa;      //< the number of hetet loci for the first individual
		UInt32 N2_Aa;      //< the number of hetet loci for the second individual
		TKINGRobustFlag() { IBS0 = nLoci = SumSq = N1_Aa = N2_Aa = 0; }
	};


	/// The pointer to the variable 'PublicDist' in the function "DoDistCalculate"
	/// The structure of genetic distance
	struct TDissflag
	{
		Int64 SumGeno;
		double SumAFreq;
		TDissflag() { SumGeno = 0; SumAFreq = 0; }
	};



	// TInit object
	class TInit
	{
	public:
		TInit()
		{
			#define PACKED_COND(cond, var, op)	\
				for (int s=0; s < _SIZE_; s++)	\
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
			PACKED_COND((b1 < 3) && (b2 < 3), Gen_Dist_SNP, sum += b1*(2-b2) + (2-b1)*b2);
			PACKED_COND((b1 < 3) && (b2 < 3), Gen_Both_Valid, sum |= (1 << i));

			/// \sum_m (X_m^{(i)} - X_m^{(j)})^2
			PACKED_COND((b1 < 3) && (b2 < 3), Gen_KING_SqDiff, sum += (b1-b2)*(b1-b2));
			PACKED_COND((b1 < 3) && (b2 < 3), Gen_KING_N1_Aa, sum += (b1==1) ? 1:0);
			PACKED_COND((b1 < 3) && (b2 < 3), Gen_KING_N2_Aa, sum += (b2==1) ? 1:0);
			PACKED_COND((b1 < 3) && (b2 < 3), Gen_KING_Num_Loci, sum ++);
		}
	} InitObj;


	/// detect the effective value for BlockSNP
	void AutoDetectSNPBlockSize(int nSamp, bool Detect=true)
	{
		if (Detect)
		{
			long L2Cache = conf_GetL2CacheMemory();
			if (L2Cache <= 0) L2Cache = 1024*1024; // 1M
			BlockSNP = (L2Cache - 3*256*256 - 8*1024) / nSamp * 4;
		}
		BlockSNP = (BlockSNP / 4) * 4;
		if (BlockSNP < 16) BlockSNP = 16;
	}

	/// Convert the raw genotypes
	static void _Do_IBS_ReadBlock(UInt8 *GenoBuf, long Start, long SNP_Cnt, void* Param)
	{
		// initialize
		const int nSamp = MCWorkingGeno.Space.SampleNum();
		UInt8 *pG = GenoBuf;
		UInt8 *pPack = GenoPacked.get();

		// pack genotypes
		for (long iSamp=0; iSamp < nSamp; iSamp++)
		{
			pPack = PackGenotypes(pG, SNP_Cnt, pPack);
			pG += SNP_Cnt;
		}
	}

	/// Compute the pairwise IBS matrix for PLINK
	static void _Do_PLINKIBS_Compute(int ThreadIndex, long Start, long SNP_Cnt, void* Param)
	{
		long Cnt = PLINKIBS_Thread_MatCnt[ThreadIndex];
		IdMatTriD I = PLINKIBS_Thread_MatIdx[ThreadIndex];
		TIBS_Flag *p = ((TIBS_Flag*)Param) + I.Offset();
		long _PackSNPLen = (SNP_Cnt / 4) + (SNP_Cnt % 4 ? 1 : 0);

		for (; Cnt > 0; Cnt--, ++I, p++)
		{
			UInt8 *p1 = GenoPacked.get() + I.Row()*_PackSNPLen;
			UInt8 *p2 = GenoPacked.get() + I.Column()*_PackSNPLen;
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
	static void _Do_IBS_Compute(int ThreadIndex, long Start, long SNP_Cnt, void* Param)
	{
		long Cnt = IBS_Thread_MatCnt[ThreadIndex];
		IdMatTri I = IBS_Thread_MatIdx[ThreadIndex];
		TIBS_Flag *p = ((TIBS_Flag*)Param) + I.Offset();
		long _PackSNPLen = (SNP_Cnt / 4) + (SNP_Cnt % 4 ? 1 : 0);

		for (; Cnt > 0; Cnt--, ++I, p++)
		{
			UInt8 *p1 = GenoPacked.get() + I.Row()*_PackSNPLen;
			UInt8 *p2 = GenoPacked.get() + I.Column()*_PackSNPLen;
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
	void DoPLINKIBSCalculate(CdMatTriDiag<TIBS_Flag> &PublicIBS, int NumThread,
		const char *Info, bool verbose)
	{
		// Initialize ...
		GenoPacked.reset(new UInt8[BlockSNP * PublicIBS.N()]);
		memset(PublicIBS.get(), 0, sizeof(TIBS_Flag)*PublicIBS.Size());

		MCWorkingGeno.Progress.Info = Info;
		MCWorkingGeno.Progress.Show() = verbose;
		MCWorkingGeno.InitParam(true, true, BlockSNP);

		MCWorkingGeno.SplitJobs(NumThread, PublicIBS.N(), PLINKIBS_Thread_MatIdx, PLINKIBS_Thread_MatCnt);
		MCWorkingGeno.Run(NumThread, &_Do_IBS_ReadBlock, &_Do_PLINKIBS_Compute, PublicIBS.get());
	}

	/// Calculate the IBS matrix
	void DoIBSCalculate(CdMatTri<TIBS_Flag> &PublicIBS, int NumThread,
		const char *Info, bool verbose)
	{
		// Initialize ...
		GenoPacked.reset(new UInt8[BlockSNP * PublicIBS.N()]);
		memset(PublicIBS.get(), 0, sizeof(TIBS_Flag)*PublicIBS.Size());

		MCWorkingGeno.Progress.Info = Info;
		MCWorkingGeno.Progress.Show() = verbose;
		MCWorkingGeno.InitParam(true, true, BlockSNP);

		MCWorkingGeno.SplitJobs(NumThread, PublicIBS.N(), IBS_Thread_MatIdx, IBS_Thread_MatCnt);
		MCWorkingGeno.Run(NumThread, &_Do_IBS_ReadBlock, &_Do_IBS_Compute, PublicIBS.get());
	}


	/// *********************************************************************************
	/// **  KING robust estimator  **
	/// *********************************************************************************

	/// Convert the raw genotypes
	static void _Do_KING_ReadBlock(UInt8 *GenoBuf, long Start, long SNP_Cnt, void* Param)
	{
		// initialize
		const int nSamp = MCWorkingGeno.Space.SampleNum();
		UInt8 *pG = GenoBuf;
		UInt8 *pPack = GenoPacked.get();

		// pack genotypes
		for (long iSamp=0; iSamp < nSamp; iSamp++)
		{
			pPack = PackGenotypes(pG, SNP_Cnt, pPack);
			pG += SNP_Cnt;
		}
		// calculate the allele frequencies
		for (long iSNP=0; iSNP < SNP_Cnt; iSNP++)
		{
			UInt8 *p = GenoBuf + iSNP;
			double &Freq = GenoAlleleFreq.get()[iSNP];
			int n = 0; Freq = 0;
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
	static void _Do_KING_Homo_Compute(int ThreadIndex, long Start, long SNP_Cnt, void* Param)
	{
		long Cnt = IBS_Thread_MatCnt[ThreadIndex];
		IdMatTri I = IBS_Thread_MatIdx[ThreadIndex];
		TKINGHomoFlag *p = ((TKINGHomoFlag*)Param) + I.Offset();
		long _PackSNPLen = (SNP_Cnt / 4) + (SNP_Cnt % 4 ? 1 : 0);

		for (; Cnt > 0; Cnt--, ++I, p++)
		{
			UInt8 *p1 = GenoPacked.get() + I.Row()*_PackSNPLen;
			UInt8 *p2 = GenoPacked.get() + I.Column()*_PackSNPLen;
			for (long k=0; k < _PackSNPLen; k++, p1++, p2++)
			{
				size_t t = (size_t(*p1) << 8) | (*p2);

				p->IBS0 += IBS0_Num_SNP[t];
				p->SumSq += Gen_KING_SqDiff[t];

				UInt8 flag = Gen_Both_Valid[t];
				if (flag & 0x01)
				{
					double f = GenoAlleleFreq.get()[4*k + 0];
					p->SumAFreq += f; p->SumAFreq2 += f*f;
				}
				if (flag & 0x02)
				{
					double f = GenoAlleleFreq.get()[4*k + 1];
					p->SumAFreq += f; p->SumAFreq2 += f*f;
				}
				if (flag & 0x04)
				{
					double f = GenoAlleleFreq.get()[4*k + 2];
					p->SumAFreq += f; p->SumAFreq2 += f*f;
				}
				if (flag & 0x08)
				{
					double f = GenoAlleleFreq.get()[4*k + 3];
					p->SumAFreq += f; p->SumAFreq2 += f*f;
				}
			}
		}
	}

	/// Compute IBD estimator in KING-robust
	static void _Do_KING_Robust_Compute(int ThreadIndex, long Start, long SNP_Cnt, void* Param)
	{
		long Cnt = IBS_Thread_MatCnt[ThreadIndex];
		IdMatTri I = IBS_Thread_MatIdx[ThreadIndex];
		TKINGRobustFlag *p = ((TKINGRobustFlag*)Param) + I.Offset();
		long _PackSNPLen = (SNP_Cnt / 4) + (SNP_Cnt % 4 ? 1 : 0);

		for (; Cnt > 0; Cnt--, ++I, p++)
		{
			UInt8 *p1 = GenoPacked.get() + I.Row()*_PackSNPLen;
			UInt8 *p2 = GenoPacked.get() + I.Column()*_PackSNPLen;
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

	/// Calculate KING IBD estimators
	void DoKINGCalculate(CdMatTri<TKINGHomoFlag> &PublicKING, int NumThread,
		const char *Info, bool verbose)
	{
		// Initialize ...
		GenoPacked.reset(new UInt8[BlockSNP * PublicKING.N()]);
		memset(PublicKING.get(), 0, sizeof(TKINGHomoFlag)*PublicKING.Size());
		GenoAlleleFreq.reset(new double[BlockSNP]);

		MCWorkingGeno.Progress.Info = Info;
		MCWorkingGeno.Progress.Show() = verbose;
		MCWorkingGeno.InitParam(true, true, BlockSNP);

		MCWorkingGeno.SplitJobs(NumThread, PublicKING.N(), IBS_Thread_MatIdx, IBS_Thread_MatCnt);
		MCWorkingGeno.Run(NumThread, &_Do_KING_ReadBlock, &_Do_KING_Homo_Compute, PublicKING.get());
	}

	/// Calculate KING IBD estimators
	void DoKINGCalculate(CdMatTri<TKINGRobustFlag> &PublicKING, int NumThread,
		const char *Info, bool verbose)
	{
		// Initialize ...
		GenoPacked.reset(new UInt8[BlockSNP * PublicKING.N()]);
		memset(PublicKING.get(), 0, sizeof(TKINGRobustFlag)*PublicKING.Size());
		GenoAlleleFreq.reset(new double[BlockSNP]);

		MCWorkingGeno.Progress.Info = Info;
		MCWorkingGeno.Progress.Show() = verbose;
		MCWorkingGeno.InitParam(true, true, BlockSNP);

		MCWorkingGeno.SplitJobs(NumThread, PublicKING.N(), IBS_Thread_MatIdx, IBS_Thread_MatCnt);
		MCWorkingGeno.Run(NumThread, &_Do_KING_ReadBlock, &_Do_KING_Robust_Compute, PublicKING.get());
	}



	/// *********************************************************************************
	/// **  Individual Dissimilarity  **
	/// *********************************************************************************

	/// Convert the raw genotypes
	static void _Do_Diss_ReadBlock(UInt8 *GenoBuf, long Start, long SNP_Cnt, void* Param)
	{
		// initialize
		const int nSamp = MCWorkingGeno.Space.SampleNum();
		UInt8 *pG = GenoBuf;
		UInt8 *pPack = GenoPacked.get();

		// pack genotypes
		for (long iSamp=0; iSamp < nSamp; iSamp++)
		{
			pPack = PackGenotypes(pG, SNP_Cnt, pPack);
			pG += SNP_Cnt;
		}
		// calculate the allele frequencies
		for (long iSNP=0; iSNP < SNP_Cnt; iSNP++)
		{
			UInt8 *p = GenoBuf + iSNP;
			double &Freq = GenoAlleleFreq.get()[iSNP];
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
	static void _Do_Diss_Compute(int ThreadIndex, long Start, long SNP_Cnt, void* Param)
	{
		long Cnt = IBS_Thread_MatCnt[ThreadIndex];
		IdMatTri I = IBS_Thread_MatIdx[ThreadIndex];
		TDissflag *p = ((TDissflag*)Param) + I.Offset();
		long _PackSNPLen = (SNP_Cnt / 4) + (SNP_Cnt % 4 ? 1 : 0);

		for (; Cnt > 0; Cnt--, ++I, p++)
		{
			UInt8 *p1 = GenoPacked.get() + I.Row()*_PackSNPLen;
			UInt8 *p2 = GenoPacked.get() + I.Column()*_PackSNPLen;
			for (long k=0; k < _PackSNPLen; k++, p1++, p2++)
			{
				size_t t = (size_t(*p1) << 8) | (*p2);
				p->SumGeno += Gen_Dist_SNP[t];

				UInt8 flag = Gen_Both_Valid[t];
				if (flag & 0x01) p->SumAFreq += GenoAlleleFreq.get()[4*k];
				if (flag & 0x02) p->SumAFreq += GenoAlleleFreq.get()[4*k+1];
				if (flag & 0x04) p->SumAFreq += GenoAlleleFreq.get()[4*k+2];
				if (flag & 0x08) p->SumAFreq += GenoAlleleFreq.get()[4*k+3];
			}
		}
	}

	/// Calculate the genetic distance matrix
	void DoDissCalculate(CdMatTri<TDissflag> &PublicDist, int NumThread,
		const char *Info, bool verbose)
	{
		// Initialize ...
		GenoPacked.reset(new UInt8[BlockSNP * PublicDist.N()]);
		memset(PublicDist.get(), 0, sizeof(TDissflag)*PublicDist.Size());
		GenoAlleleFreq.reset(new double[BlockSNP]);

		MCWorkingGeno.Progress.Info = Info;
		MCWorkingGeno.Progress.Show() = verbose;
		MCWorkingGeno.InitParam(true, true, BlockSNP);

		MCWorkingGeno.SplitJobs(NumThread, PublicDist.N(), IBS_Thread_MatIdx, IBS_Thread_MatCnt);
		MCWorkingGeno.Run(NumThread, &_Do_Diss_ReadBlock, &_Do_Diss_Compute, PublicDist.get());
	}
}

#endif  /* _FuncIBS_H_ */
