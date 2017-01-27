// ===========================================================
//
// genPCA.h: Principal Component Analysis on GWAS
//
// Copyright (C) 2011-2017    Xiuwen Zheng
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


#ifndef _HEADER_PCA_
#define _HEADER_PCA_

// CoreArray library header
#include "dGenGWAS.h"
#include "dVect.h"
#include "ThreadPool.h"


namespace PCA
{
	// using namespace
	using namespace std;
	using namespace CoreArray;
	using namespace Vectorization;
	using namespace GWAS;


	// ---------------------------------------------------------------------
	// Vectorization Computing, algorithm 1 (Arithmetic)

	class COREARRAY_DLL_LOCAL CProdMat_Base
	{
	public:
		VEC_AUTO_PTR<C_Int32> GenoSum, GenoNum;
		VEC_AUTO_PTR<double> avg_geno;  // avg geno: 2*\bar{p}, or ...

		CProdMat_Base() { }
		/// initialize GenoSum, GenoNum
		void ZeroFill();
		/// calculate average genotypes (save to GenoSum and GenoNum)
		void SummarizeGeno_SampxSNP(C_UInt8 *pGeno, size_t nSNP);
		/// divide GenoSum by GenoNum (GenoSum / GenoNum)
		void DivideGeno();
		/// computing scalar vector
		void rsqrt_prod();

		inline size_t N() const { return fN; }
		inline size_t M() const { return fM; }

	protected:
		size_t fN, fM;
		void _Reset();
		void _Clear();
	};


	// Mean-adjusted genotype matrix (AVX, N: # of samples, M: # of SNPs)
	class COREARRAY_DLL_LOCAL CProdMat_AlgArith: public CProdMat_Base
	{
	public:
		CProdMat_AlgArith() { fN = fM = 0; }
		CProdMat_AlgArith(size_t n, size_t m) { Reset(n, m); }

		void Reset(size_t n, size_t m);
		void Clear();
		/// detect the effective value for BlockNumSNP
		static void PCA_Detect_BlockNumSNP(int nSamp);
		// time-consuming function
		void COREARRAY_CALL_ALIGN MulAdd(IdMatTri &Idx, size_t IdxCnt, double *pOut);
		/// mean-adjusted genotypes (fGenotype - avg_geno)
		void COREARRAY_CALL_ALIGN GenoSub();
		/// variance-adjusted genotypes (fGenotype * avg_geno)
		void COREARRAY_CALL_ALIGN GenoMul();

		inline double *base() { return fGenotype.Get(); }

	protected:
		VEC_AUTO_PTR<double> fGenotype;

		inline void TransposeGenotype(size_t nSamp, size_t nSNP, C_UInt8 *pGeno)
		{
			double *p = fGenotype.Get();
			for (size_t i=0; i < nSamp; i++)
			{
				double *pp = p; p += fM;
				size_t m = nSNP;
				for (size_t j=0; j < m; j++)
				{
					C_UInt8 g = pGeno[nSamp*j + i];
					*pp ++ = (g <= 2) ? g : avg_geno[j];
				}
				// zero fill the left part
				for (; m < fM; m++) *pp ++ = 0;
			}
		}
	};

}


extern "C"
{
	/// get the eigenvalues and eigenvectors, return 'nProtect'
	COREARRAY_DLL_LOCAL int CalcEigen(double *pMat, int n, int nEig,
		const char *EigMethod, SEXP &EigVal, SEXP &EigVect);
}

#endif  /* _HEADER_PCA_ */
