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


	/// The mutex object for ptrPublicCov
	TdMutex _IBSMutex = NULL;

	/// Packed size
	static const long _size = 256*256;

	/// The number of IBS 0 in the packed genotype
	UInt8 IBS0_Num_SNP[_size];
	/// The number of IBS 1 in the packed genotype
	UInt8 IBS1_Num_SNP[_size];
	/// The number of IBS 2 in the packed genotype
	UInt8 IBS2_Num_SNP[_size];

	/// The pointer to the variable PublicIBS in the function "DoIBSCalculate"
	/// The structure of IBS states
	struct TIBSflag
	{
		UInt32 IBS0, IBS1, IBS2;
		TIBSflag() { IBS0 = IBS1 = IBS2 = 0; }
	};
	TIBSflag *ptrPublicIBS;



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

	/// The thread entry for the calculation of genetic covariace matrix
	void Entry_IBSCalc(TdThread Thread, int ThreadIndex, void *Param)
	{
		// The number of working samples
		const long n = GenoSpace.SampleNum();

		// The buffer of genotype
		auto_ptr<UInt8> GenoBlock(new UInt8[n * BlockSNP]);
		auto_ptr<UInt8> GenoPacked(new UInt8[n * BlockSNP / 4]);

		CdMatTriDiag<TIBSflag> fpIBS;
		TIBSflag *ptrfpIBS = ptrPublicIBS;
		if (Thread)
		{
			fpIBS.Reset(n); fpIBS.Clear(TIBSflag());
			ptrfpIBS = fpIBS.get();
		}

		long _SNPstart, _SNPlen;
		while (RequireWork(GenoBlock.get(), _SNPstart, _SNPlen, true))
		{
			UInt8 *pG = GenoBlock.get();
			UInt8 *pPack = GenoPacked.get();

			// pack genotypes
			for (long iSamp=0; iSamp < n; iSamp++)
			{
				pPack = PackGenotypes(pG, _SNPlen, pPack);
				pG += _SNPlen;
			}

			TIBSflag *pr = ptrfpIBS;
			long _PackSNPLen = (_SNPlen / 4) + (_SNPlen % 4 ? 1 : 0);

			for (long i=0; i < n; i++)
			{
				for (long j=i+1; j < n; j++)
				{
					UInt8 *p1 = GenoPacked.get() + i*_PackSNPLen;
					UInt8 *p2 = GenoPacked.get() + j*_PackSNPLen;
					for (long k=_PackSNPLen; k > 0; k--, p1++, p2++)
					{
						size_t t = (size_t(*p1) << 8) | (*p2);
						pr->IBS0 += IBS0_Num_SNP[t];
						pr->IBS1 += IBS1_Num_SNP[t];
						pr->IBS2 += IBS2_Num_SNP[t];
					}
					pr++;
				}
			}

			// Update progress
			{
				TdAutoMutex _m(_Mutex);
				Progress.Forward(_SNPlen);
			}
		}

		// finally, update fpCov
		if (Thread)
		{
			// auto Lock and Unlock
			TdAutoMutex _m(_IBSMutex);
			TIBSflag *s = fpIBS.get(), *d = ptrPublicIBS;
			for (size_t i=fpIBS.Size(); i > 0; i--)
			{
				d->IBS0 += s->IBS0; d->IBS1 += s->IBS1;
				d->IBS2 += s->IBS2;
				d++; s++;
			}
		} else
			// Unlock the elements of the variable "IBSMutex"
			plc_UnlockMutex(_IBSMutex);
	}

	/// Calculate the IBS matrix
	void DoIBSCalculate(CdMatTriDiag<TIBSflag> &PublicIBS, int NumThread,
		const char *Info, bool verbose)
	{
		// initialize mutex objects
		_Mutex = plc_InitMutex();
		_IBSMutex = plc_InitMutex();

		PublicIBS.Clear(TIBSflag());
		ptrPublicIBS = PublicIBS.get();

		// Lock the elements of the varaible "PublicIBS"
		plc_LockMutex(_IBSMutex);

		// Initialize progress information
		Progress.Info = Info;
		Progress.Show() = verbose;
		Progress.Init(GenoSpace.SNPNum());
		SNPStart = 0;

		// Threads
		plc_DoBaseThread(Entry_IBSCalc, NULL, NumThread);

		/// destroy mutex objects
		plc_DoneMutex(_Mutex); _Mutex = NULL;
		plc_DoneMutex(_IBSMutex); _IBSMutex = NULL;
	}

}


#endif  /* _FuncIBS_H_ */
