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
// Copyright (C) 2011	Xiuwen Zheng
//
// This file is part of CoreArray.
//
// CoreArray is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License Version 3 as
// published by the Free Software Foundation.
//
// CoreArray is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with CoreArray.
// If not, see <http://www.gnu.org/licenses/>.


// CoreArray library header
#include <dType.hpp>
#include <dVect.hpp>
#include <CoreGDSLink.hpp>
#include <dGenGWAS.hpp>

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

	/// The number of IBS 0 in the packed genotype
	UInt8 IBS0_Num_SNP[_size];
	/// The number of IBS 1 in the packed genotype
	UInt8 IBS1_Num_SNP[_size];
	/// The number of IBS 2 in the packed genotype
	UInt8 IBS2_Num_SNP[_size];

	/// The packed genotype buffer
	auto_ptr<UInt8> GenoPacked;

	/// Thread variables
	const int N_MAX_THREAD = 256;
	IdMatTriD IBS_Thread_MatIdx[N_MAX_THREAD];
	Int64 IBS_Thread_MatCnt[N_MAX_THREAD];


	/// The pointer to the variable PublicIBS in the function "DoIBSCalculate"
	/// The structure of IBS states
	struct TIBSflag
	{
		UInt32 IBS0, IBS1, IBS2;
		TIBSflag() { IBS0 = IBS1 = IBS2 = 0; }
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

	/// Convert the raw genotypes to the mean-adjusted genotypes
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

}


#endif  /* _FuncIBS_H_ */
