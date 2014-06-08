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
	static const long _size = 256*256;

	/// IBS
	/// The number of IBS 0 in the packed genotype
	UInt8 IBS0_Num_SNP[_size];
	/// The number of IBS 1 in the packed genotype
	UInt8 IBS1_Num_SNP[_size];
	/// The number of IBS 2 in the packed genotype
	UInt8 IBS2_Num_SNP[_size];

	/// Genetic Distance
	/// The distance in the packed genotype
	UInt8 Gen_Dist_SNP[_size];
	/// The flag of use of allele frequencies
	UInt8 Gen_Freq_Flag[_size];


	/// The packed genotype buffer
	auto_ptr<UInt8> GenoPacked;
	auto_ptr<double> GenoAlleleFreq;

	/// Thread variables
	const int N_MAX_THREAD = 256;
	// IBS
	IdMatTriD IBS_Thread_MatIdx[N_MAX_THREAD];
	Int64 IBS_Thread_MatCnt[N_MAX_THREAD];
	// Individual Similarity
	IdMatTri IS_Thread_MatIdx[N_MAX_THREAD];
	Int64 IS_Thread_MatCnt[N_MAX_THREAD];


	/// The pointer to the variable PublicIBS in the function "DoIBSCalculate"
	/// The structure of IBS states
	struct TIBSflag
	{
		UInt32 IBS0, IBS1, IBS2;
		TIBSflag() { IBS0 = IBS1 = IBS2 = 0; }
	};

	/// The pointer to the variable PublicIBS in the function "DoDistCalculate"
	/// The structure of genetic distance
	struct TDistflag
	{
		Int64 SumGeno;
		double SumAFreq;
		TDistflag() { SumGeno = 0; SumAFreq = 0; }
	};


	// TInit object
	class TInit
	{
	public:
		TInit()
		{
			#define PACKED_COND(cond, var, op)	\
				for (int s=0; s < _size; s++)	\
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
			PACKED_COND((b1 < 3) && (b2 < 3), Gen_Freq_Flag, sum |= (1 << i));
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
		// init ...
		const int n = MCWorkingGeno.Space.SampleNum();
		UInt8 *pG = GenoBuf;
		UInt8 *pPack = GenoPacked.get();

		// pack genotypes
		for (long iSamp=0; iSamp < n; iSamp++)
		{
			pPack = PackGenotypes(pG, SNP_Cnt, pPack);
			pG += SNP_Cnt;
		}
	}

	/// Compute the covariate matrix
	static void _Do_IBS_Compute(int ThreadIndex, long Start, long SNP_Cnt, void* Param)
	{
		long Cnt = IBS_Thread_MatCnt[ThreadIndex];
		IdMatTriD I = IBS_Thread_MatIdx[ThreadIndex];
		TIBSflag *p = ((TIBSflag*)Param) + I.Offset();
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

	/// Calculate the IBS matrix
	void DoIBSCalculate(CdMatTriDiag<TIBSflag> &PublicIBS, int NumThread,
		const char *Info, bool verbose)
	{
		// Initialize ...
		GenoPacked.reset(new UInt8[BlockSNP * PublicIBS.N()]);
		memset(PublicIBS.get(), 0, sizeof(TIBSflag)*PublicIBS.Size());

		MCWorkingGeno.Progress.Info = Info;
		MCWorkingGeno.Progress.Show() = verbose;
		MCWorkingGeno.InitParam(true, true, BlockSNP);

		MCWorkingGeno.SplitJobs(NumThread, PublicIBS.N(), IBS_Thread_MatIdx, IBS_Thread_MatCnt);
		MCWorkingGeno.Run(NumThread, &_Do_IBS_ReadBlock, &_Do_IBS_Compute, PublicIBS.get());
	}



	/// *********************************************************************************
	/// **    **
	/// *********************************************************************************

	/// Convert the raw genotypes
	static void _Do_Dist_ReadBlock(UInt8 *GenoBuf, long Start, long SNP_Cnt, void* Param)
	{
		// init ...
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
	static void _Do_Dist_Compute(int ThreadIndex, long Start, long SNP_Cnt, void* Param)
	{
		long Cnt = IS_Thread_MatCnt[ThreadIndex];
		IdMatTri I = IS_Thread_MatIdx[ThreadIndex];
		TDistflag *p = ((TDistflag*)Param) + I.Offset();
		long _PackSNPLen = (SNP_Cnt / 4) + (SNP_Cnt % 4 ? 1 : 0);

		for (; Cnt > 0; Cnt--, ++I, p++)
		{
			UInt8 *p1 = GenoPacked.get() + I.Row()*_PackSNPLen;
			UInt8 *p2 = GenoPacked.get() + I.Column()*_PackSNPLen;
			for (long k=0; k < _PackSNPLen; k++, p1++, p2++)
			{
				size_t t = (size_t(*p1) << 8) | (*p2);
				p->SumGeno += Gen_Dist_SNP[t];

				UInt8 flag = Gen_Freq_Flag[t];
				if (flag & 0x01) p->SumAFreq += GenoAlleleFreq.get()[4*k];
				if (flag & 0x02) p->SumAFreq += GenoAlleleFreq.get()[4*k+1];
				if (flag & 0x04) p->SumAFreq += GenoAlleleFreq.get()[4*k+2];
				if (flag & 0x08) p->SumAFreq += GenoAlleleFreq.get()[4*k+3];
			}
		}
	}

	/// Calculate the genetic distance matrix
	void DoDistCalculate(CdMatTri<TDistflag> &PublicDist, int NumThread,
		const char *Info, bool verbose)
	{
		// Initialize ...
		GenoPacked.reset(new UInt8[BlockSNP * PublicDist.N()]);
		memset(PublicDist.get(), 0, sizeof(TDistflag)*PublicDist.Size());
		GenoAlleleFreq.reset(new double[BlockSNP]);

		MCWorkingGeno.Progress.Info = Info;
		MCWorkingGeno.Progress.Show() = verbose;
		MCWorkingGeno.InitParam(true, true, BlockSNP);

		MCWorkingGeno.SplitJobs(NumThread, PublicDist.N(), IS_Thread_MatIdx, IS_Thread_MatCnt);
		MCWorkingGeno.Run(NumThread, &_Do_Dist_ReadBlock, &_Do_Dist_Compute, PublicDist.get());
	}
}


#endif  /* _FuncIBS_H_ */
