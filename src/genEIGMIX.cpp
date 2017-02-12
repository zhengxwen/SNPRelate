// ===========================================================
//
// genEIGMIX.cpp: Eigen-analysis with admixture on GWAS
//
// Copyright (C) 2017    Xiuwen Zheng
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


#ifndef _HEADER_EIGMIX_
#define _HEADER_EIGMIX_

#include "genPCA.h"
#include <algorithm>


namespace EIGMIX
{

// using namespace
using namespace std;
using namespace CoreArray;
using namespace Vectorization;
using namespace GWAS;


// -----------------------------------------------------------
// EIGMIX IBD Calculation with arithmetic algorithm

class COREARRAY_DLL_LOCAL CEigMix_AlgArith: protected PCA::CProdMat_AlgArith
{
private:
	CdBaseWorkSpace &Space;
	double *ptrIBD;

	void thread_cov_outer(size_t i, size_t n)
	{
		IdMatTri I = Array_Thread_MatIdx[i];
		MulAdd(I, Array_Thread_MatCnt[i], ptrIBD + I.Offset());
	}

public:
	CEigMix_AlgArith(CdBaseWorkSpace &space): CProdMat_AlgArith(), Space(space)
		{ }

	/// run the algorithm
	void Run(CdMatTri<double> &IBD, int NumThread, double AFreq[], bool DiagAdj,
		bool verbose)
	{
		if (NumThread < 1) NumThread = 1;
		const size_t nSamp = Space.SampleNum();

		// detect the appropriate block size
		PCA_Detect_BlockNumSNP(nSamp);
		if (verbose)
		{
			Rprintf("%s    (internal increment: %d)\n", TimeToStr(),
				(int)BlockNumSNP);
		}

		// initialize
		Reset(nSamp, BlockNumSNP);
		ptrIBD = IBD.Get();
		memset(ptrIBD, 0, sizeof(double)*IBD.Size());

		// denominator 4*p*(1-p)
		CdMatTri<double> Denom(nSamp);
		memset(Denom.Get(), 0, sizeof(double)*Denom.Size());
		double SumDenominator = 0;

		// diagonal
		vector<int> DiagAdjVal(nSamp, 0);

		// thread pool
		CThreadPoolEx<CEigMix_AlgArith> thpool(NumThread);
		Array_SplitJobs(NumThread, nSamp, Array_Thread_MatIdx,
			Array_Thread_MatCnt);

		// genotypes (0, 1, 2 and NA)
		VEC_AUTO_PTR<C_UInt8> Geno(nSamp * BlockNumSNP);
		C_UInt8 *pGeno = Geno.Get();

		// genotype buffer, false for no memory buffer
		CGenoReadBySNP WS(NumThread, Space, BlockNumSNP, verbose ? -1 : 0, false);
		WS.Init();

		// for-loop
		while (WS.Read(Geno.Get()))
		{
			// get the genotype sum and number
			SummarizeGeno_SampxSNP(pGeno, WS.Count());
			// get 2 * \bar{p}_l, saved in PCA_Mat.avg_geno
			DivideGeno();

			// transpose genotypes in PCA_Mat, set missing values to the average
			TransposeGenotype(nSamp, WS.Count(), pGeno);
			// G@ij - 2\bar{p}_l
			GenoSub();

			// denominator
			size_t nsnp = WS.Count();
			C_UInt8 *pG = pGeno;
			for (size_t i=0; i < nsnp; i++)
			{
				const double af = 0.5 * avg_geno[i];
				if (AFreq)
					AFreq[WS.Index() + i] = af;
				double denom = 4 * af * (1 - af);
				SumDenominator += denom;
				C_UInt8 *pGG = pG;
				for (size_t j=0; j < nSamp; j++, pG++)
				{
					if (*pG == 1)
					{
						DiagAdjVal[j] ++;
					} else if (*pG > 2)
					{
						double *pp = Denom.Get();
						vec_f64_add(pp + j*(2*nSamp-j+1)/2, nSamp-j, denom);
						for (ssize_t k=(ssize_t)j - 1; k >= 0; k--)
							if (pGG[k] <= 2)
								pp[j + k*(2*nSamp-k-1)/2] += denom;
					}
				}
			}

			// outer product, using thread thpool
			thpool.BatchWork(this, &CEigMix_AlgArith::thread_cov_outer, NumThread);
			// update
			WS.ProgressForward(WS.Count());
		}

		// diagonal adjustment 
		if (DiagAdj)
		{
			for (size_t i=0; i < nSamp; i++)
				IBD.At(i, i) -= DiagAdjVal[i];
		}
		// normalized by the denominator
		double *p = IBD.Get(), *s = Denom.Get();
		for (size_t n=IBD.Size(); n > 0; n--)
			*p++ /= (SumDenominator - *s++);
	}
};



// -----------------------------------------------------------
// EIGMIX IBD Calculation with genotype indexing

class COREARRAY_DLL_LOCAL CEigMix_AlgIndexing
{
private:
	CdBaseWorkSpace &Space;
	double *ptrNom, *ptrDenom;
	size_t nSamp;          ///< the number of samples
	vector<C_UInt8> Geno;  ///< genotypes (0, 1, 2 and NA)
	vector< vector<size_t> > GenoIndex;
	double SumNumerator, SumDenominator;
	CMutex mutex;
	vector<CMutex> mutex_list;

	void thread_cov(int thread_idx, size_t i, size_t num)
	{
		size_t n4[4];
		size_t *i4[4] = {
			&GenoIndex[thread_idx][0],       &GenoIndex[thread_idx][nSamp],
			&GenoIndex[thread_idx][2*nSamp], &GenoIndex[thread_idx][3*nSamp]
		};
		C_UInt8 *pGeno = &Geno[0] + nSamp * i;
		double SumNom=0, SumDenom=0;

		for (; num > 0; num--, pGeno += nSamp)
		{
			// classify genotypes
			PackGenoIndex(pGeno, nSamp, n4, i4[0], i4[1], i4[2], i4[3]);
			if (n4[2] < n4[0])
			{
				swap(n4[2], n4[0]);
				swap(i4[2], i4[0]);			
			}
			size_t *i4_end[4] = { i4[0]+n4[0], i4[1]+n4[1], i4[2]+n4[2], i4[3]+n4[3] };
			// allele frequency
			if (n4[3] >= nSamp) continue;
			double freq = double(2*n4[2] + n4[1]) / (nSamp - n4[3]);
			if (freq<=0 || freq>=2) continue;

			// update
			double nom;
			if (n4[2] >= n4[1])
			{
				// update numerator
				nom = (2 - freq) * (2 - freq);
				SumNom += nom;
				// 2 - 1 / 2 - 0
				double v0  = (2 - freq) * (-freq) - nom;
				double v1  = (2 - freq) * (1 - freq) - nom;
				for (size_t *p=i4[2]; p < i4_end[2]; p++)
				{
					mutex_list[*p].Lock();
					double *v = ptrNom + (*p)*nSamp;
					for (size_t *s=i4[1]; s < i4_end[1]; s++) v[*s] += v1;
					for (size_t *s=i4[0]; s < i4_end[0]; s++) v[*s] += v0;
					mutex_list[*p].Unlock();
				}
				// 1 - 1 / 1 - 0
				v0  = (1 - freq) * (-freq) - nom;
				v1  = (1 - freq) * (1 - freq) - nom;
				for (size_t *p=i4[1]; p < i4_end[1]; p++)
				{
					mutex_list[*p].Lock();
					double *v = ptrNom + (*p)*nSamp;
					v[*p] += v1 - 1;
					for (size_t *s=p+1; s < i4_end[1]; s++) v[*s] += v1;
					for (size_t *s=i4[0]; s < i4_end[0]; s++) v[*s] += v0;
					mutex_list[*p].Unlock();
				}
			} else {
				// update numerator
				nom = (1 - freq) * (1 - freq);
				SumNom += nom;
				// 1 - 0 / 1 - 2
				double v0  = (1 - freq) * (-freq) - nom;
				double v2  = (1 - freq) * (2 - freq) - nom;
				for (size_t *p=i4[1]; p < i4_end[1]; p++)
				{
					mutex_list[*p].Lock();
					double *v = ptrNom + (*p)*nSamp;
					v[*p] -= 1;  // 1 - 1
					for (size_t *s=i4[2]; s < i4_end[2]; s++) v[*s] += v2;
					for (size_t *s=i4[0]; s < i4_end[0]; s++) v[*s] += v0;
					mutex_list[*p].Unlock();
				}
				// 2 - 2 / 2 - 0
				v0  = (2 - freq) * (-freq) - nom;
				v2  = (2 - freq) * (2 - freq) - nom;
				for (size_t *p=i4[2]; p < i4_end[2]; p++)
				{
					mutex_list[*p].Lock();
					double *v = ptrNom + (*p)*nSamp;
					for (size_t *s=p; s < i4_end[2]; s++) v[*s] += v2;
					for (size_t *s=i4[0]; s < i4_end[0]; s++) v[*s] += v0;
					mutex_list[*p].Unlock();
				}
			}
			// 0 - 0
			double v0 = freq * freq - nom;
			for (size_t *p=i4[0]; p < i4_end[0]; p++)
			{
				mutex_list[*p].Lock();
				double *v = ptrNom + (*p)*nSamp;
				for (size_t *s=p; s < i4_end[0]; s++) v[*s] += v0;
				mutex_list[*p].Unlock();
			}

			// update denominator, a bug, need to be fixed later
			double den = 2 * freq * (1 - 0.5*freq);
			SumDenom += den;
			for (size_t *p=i4[3]; p < i4_end[3]; p++)
				vec_f64_add(ptrDenom + (*p) * nSamp, den, nSamp);
		}

		mutex.Lock();
		SumNumerator += SumNom;
		SumDenominator += SumDenom;
		mutex.Unlock();
	}

	void thread_cov_nosync(size_t num)
	{
		size_t n4[4];
		size_t *i4[4] = {
			&GenoIndex[0][0],       &GenoIndex[0][nSamp],
			&GenoIndex[0][2*nSamp], &GenoIndex[0][3*nSamp]
		};
		C_UInt8 *pGeno = &Geno[0];
		double SumNom=0, SumDenom=0;

		for (; num > 0; num--, pGeno += nSamp)
		{
			// classify genotypes
			PackGenoIndex(pGeno, nSamp, n4, i4[0], i4[1], i4[2], i4[3]);
			if (n4[2] < n4[0])
			{
				swap(n4[2], n4[0]);
				swap(i4[2], i4[0]);			
			}
			size_t *i4_end[4] = { i4[0]+n4[0], i4[1]+n4[1], i4[2]+n4[2], i4[3]+n4[3] };
			// allele frequency
			if (n4[3] >= nSamp) continue;
			double freq = double(2*n4[2] + n4[1]) / (nSamp - n4[3]);
			if (freq<=0 || freq>=2) continue;

			// update
			double nom;
			if (n4[2] >= n4[1])
			{
				// update numerator
				nom = (2 - freq) * (2 - freq);
				SumNom += nom;
				// 2 - 1 / 2 - 0
				double v0  = (2 - freq) * (-freq) - nom;
				double v1  = (2 - freq) * (1 - freq) - nom;
				for (size_t *p=i4[2]; p < i4_end[2]; p++)
				{
					double *v = ptrNom + (*p)*nSamp;
					for (size_t *s=i4[1]; s < i4_end[1]; s++) v[*s] += v1;
					for (size_t *s=i4[0]; s < i4_end[0]; s++) v[*s] += v0;
				}
				// 1 - 1 / 1 - 0
				v0  = (1 - freq) * (-freq) - nom;
				v1  = (1 - freq) * (1 - freq) - nom;
				for (size_t *p=i4[1]; p < i4_end[1]; p++)
				{
					double *v = ptrNom + (*p)*nSamp;
					v[*p] += v1 - 1;
					for (size_t *s=p+1; s < i4_end[1]; s++) v[*s] += v1;
					for (size_t *s=i4[0]; s < i4_end[0]; s++) v[*s] += v0;
				}
			} else {
				// update numerator
				nom = (1 - freq) * (1 - freq);
				SumNom += nom;
				// 1 - 0 / 1 - 2
				double v0  = (1 - freq) * (-freq) - nom;
				double v2  = (1 - freq) * (2 - freq) - nom;
				for (size_t *p=i4[1]; p < i4_end[1]; p++)
				{
					double *v = ptrNom + (*p)*nSamp;
					v[*p] -= 1;  // 1 - 1
					for (size_t *s=i4[2]; s < i4_end[2]; s++) v[*s] += v2;
					for (size_t *s=i4[0]; s < i4_end[0]; s++) v[*s] += v0;
				}
				// 2 - 2 / 2 - 0
				v0  = (2 - freq) * (-freq) - nom;
				v2  = (2 - freq) * (2 - freq) - nom;
				for (size_t *p=i4[2]; p < i4_end[2]; p++)
				{
					double *v = ptrNom + (*p)*nSamp;
					for (size_t *s=p; s < i4_end[2]; s++) v[*s] += v2;
					for (size_t *s=i4[0]; s < i4_end[0]; s++) v[*s] += v0;
				}
			}
			// 0 - 0
			double v0 = freq * freq - nom;
			for (size_t *p=i4[0]; p < i4_end[0]; p++)
			{
				double *v = ptrNom + (*p)*nSamp;
				for (size_t *s=p; s < i4_end[0]; s++) v[*s] += v0;
			}

			// update denominator, a bug, need to be fixed later
			double den = 2 * freq * (1 - 0.5*freq);
			SumDenom += den;
			for (size_t *p=i4[3]; p < i4_end[3]; p++)
				vec_f64_add(ptrDenom + (*p) * nSamp, den, nSamp);
		}

		SumNumerator += SumNom;
		SumDenominator += SumDenom;
	}

public:
	CEigMix_AlgIndexing(CdBaseWorkSpace &space): Space(space) { }

	/// run the algorithm
	void Run(CdMat<double> &IBD, int NumThread, bool verbose)
	{
		if (NumThread < 1) NumThread = 1;
		nSamp = Space.SampleNum();

		// initialize numerator and denominator
		ptrNom = IBD.Get();
		memset(ptrNom, 0, sizeof(double)*IBD.Size());
		CdMat<double> Denom(nSamp);
		ptrDenom = Denom.Get();
		memset(ptrDenom, 0, sizeof(double)*Denom.Size());
		SumNumerator = SumDenominator = 0;

		// genotypes (0, 1, 2 and NA)
		const size_t BlockSNPNum = 512;
		Geno.resize(nSamp * BlockSNPNum);
		GenoIndex.resize(NumThread);
		for (int i=0; i < NumThread; i++) GenoIndex[i].resize(4*nSamp);

		// thread thpool
		CThreadPoolEx<CEigMix_AlgIndexing> thpool(NumThread);
		// genotype buffer, false for no memory buffer
		CGenoReadBySNP WS(NumThread, Space, BlockSNPNum, verbose ? -1 : 0, false);
		// mutex list
		mutex_list.resize(nSamp);

		// for-loop
		WS.Init();
		while (WS.Read(&Geno[0]))
		{
			if (NumThread > 1)
			{
				// using thread thpool
				thpool.BatchWork2(this, &CEigMix_AlgIndexing::thread_cov, WS.Count());
			} else {
				thread_cov_nosync(WS.Count());
			}
			// update
			WS.ProgressForward(WS.Count());
		}

		// sum lower and upper triangles
		for (size_t i=0; i < nSamp; i++)
		{
			// diagonal
			double &v_n = ptrNom[i + nSamp*i];
			double &v_d = ptrDenom[i + nSamp*i];
			v_n = (SumNumerator + v_n) / (SumDenominator + v_d);
			// off-diagonal
			for (size_t j=i+1; j < nSamp; j++)
			{
				size_t p1 = j + nSamp*i;
				size_t p2 = i + nSamp*j;
				double &v1 = ptrNom[p1], &v2 = ptrNom[p2];
				v1 = v2 = (SumNumerator + v1 + v2) /
					(SumDenominator - ptrDenom[p1] - ptrDenom[p2]);
			}
		}
	}
};



// ================== SNP Loadings ==================

class COREARRAY_DLL_LOCAL CEigMix_SNPLoad
{
private:
	CdBaseWorkSpace &Space;  ///< working genotypes
	VEC_AUTO_PTR<C_UInt8> Geno;  ///< genotypes (0, 1, 2 and NA)
	size_t nSamp;       ///< the number of samples
	size_t NumEigVal;   ///< the number of eigenvalues
	double *pEigVect;   ///< the pointer to eigenvectors
	double *pLoading;   ///< the pointer to correlation
	double *pAFreq;     ///< the pointer to allele frequency
	double AFreqScale;  ///< the scalar of allele frequency

	void thread_loading(size_t i, size_t num)
	{
		C_UInt8 *pG = Geno.Get() + nSamp * i;
		double *pOut = pLoading + NumEigVal * i;
		for (; num > 0; num--, i++)
		{
			// zero filling
			memset(pOut, 0, sizeof(double)*NumEigVal);
			// dot product
			for (size_t j=0; j < nSamp; j++)
			{
				double g = (*pG < 3) ? ((*pG - 2*pAFreq[i]) * AFreqScale) : 0.0;
				pG ++;
				double *pEig = pEigVect + j;
				for (size_t k=0; k < NumEigVal; k++)
				{
					pOut[k] += g * (*pEig);
					pEig += nSamp;
				}
			}
			pOut += NumEigVal;
		}
	}

public:
	/// constructor
	CEigMix_SNPLoad(CdBaseWorkSpace &space): Space(space) { }

	/// run the algorithm
	void Run(double OutSNPLoading[], double AFreq[], size_t NumEig,
		double EigVect[], int NumThread, bool verbose)
	{
		if (NumThread < 1) NumThread = 1;
		nSamp = Space.SampleNum();
		NumEigVal = NumEig; pEigVect = EigVect;

		// detect the appropriate block size
		size_t Cache = GetOptimzedCache();
		size_t nBlock = Cache / nSamp;
		nBlock = (nBlock / 4) * 4;
		if (nBlock < 128) nBlock = 128;
		if (nBlock > 65536) nBlock = 65536;
		if (verbose)
			Rprintf("%s    (internal increment: %d)\n", TimeToStr(), (int)nBlock);

		// normalize
		size_t nsnp = Space.SNPNum();
		double sum = 0;
		for (size_t i=0; i < nsnp; i++)
			sum += 4 * AFreq[i] * (1 - AFreq[i]);
		AFreqScale = 1 / sqrt(sum);

		// thread thpool
		CThreadPoolEx<CEigMix_SNPLoad> thpool(NumThread);
		// genotypes (0, 1, 2 and NA)
		Geno.Reset(nSamp * nBlock);
		// genotype buffer, false for no memory buffer
		CGenoReadBySNP WS(NumThread, Space, nBlock, verbose ? -1 : 0, false);

		// for-loop
		WS.Init();
		while (WS.Read(Geno.Get()))
		{
			pLoading = OutSNPLoading + WS.Index() * NumEig;
			pAFreq = AFreq + WS.Index();
			// using thread thpool
			thpool.BatchWork(this, &CEigMix_SNPLoad::thread_loading, WS.Count());
			// update
			WS.ProgressForward(WS.Count());
		}
	}
};



// ================== Sample Loadings ==================

class COREARRAY_DLL_LOCAL CEigMix_SampleLoad
{
private:
	CdBaseWorkSpace &Space;  ///< working genotypes
	VEC_AUTO_PTR<C_UInt8> Geno;  ///< genotypes (0, 1, 2 and NA)
	size_t nSamp;    ///< the number of samples
	size_t nEigVal;  ///< the number of eigenvalues
	size_t nSubSNP;  ///< the number of genotypes in a block
	double *pLoading;   ///< the pointer to SNP eigenvectors
	double *pAFreq;     ///< the pointer to allele frequency
	double AFreqScale;  ///< the scalar of allele frequency
	double *pOutEig;    ///< the pointer to sample eigenvectors

	void thread_loading(size_t i, size_t num)
	{
		// for-loop each individual i
		for (; num > 0; num--, i++)
		{
			C_UInt8 *pGeno = Geno.Get() + i;
			double *pLoad = pLoading;
			for (size_t j=0; j < nSubSNP; j++)
			{
				double g = (*pGeno < 3) ? (*pGeno - 2*pAFreq[j]) * AFreqScale : 0.0;
				pGeno += nSamp;

				double *p = pOutEig + i;
				for (size_t k=0; k < nEigVal; k++)
				{
					*p += g * (*pLoad++);
					p += nSamp;
				}
			}
		}
	}

public:
	/// constructor
	CEigMix_SampleLoad(CdBaseWorkSpace &space): Space(space) { }

	/// run the algorithm
	void Run(double OutSampLoad[], size_t NumEig, double SNPLoading[],
		double AFreq[], int NumThread, bool verbose)
	{
		if (NumThread < 1) NumThread = 1;
		nSamp = Space.SampleNum();
		nEigVal = NumEig;
		pOutEig = OutSampLoad;

		// detect the appropriate block size
		size_t Cache = GetOptimzedCache();
		size_t nBlock = Cache / nSamp;
		nBlock = (nBlock / 4) * 4;
		if (nBlock < 128) nBlock = 128;
		if (nBlock > 65536) nBlock = 65536;
		if (verbose)
			Rprintf("%s    (internal increment: %d)\n", TimeToStr(), (int)nBlock);

		// thread thpool
		CThreadPoolEx<CEigMix_SampleLoad> thpool(NumThread);
		// genotypes (0, 1, 2 and NA)
		Geno.Reset(nSamp * nBlock);
		// genotype buffer, false for no memory buffer
		CGenoReadBySNP WS(NumThread, Space, nBlock, verbose ? -1 : 0, false);

		// zero filling
		memset(OutSampLoad, 0, sizeof(double)*NumEig*nSamp);
		// normalize
		size_t nsnp = Space.SNPNum();
		double sum = 0;
		for (size_t i=0; i < nsnp; i++)
			sum += 4 * AFreq[i] * (1 - AFreq[i]);
		AFreqScale = 1 / sqrt(sum);

		// for-loop
		WS.Init();
		while (WS.Read(Geno.Get()))
		{
			pLoading = SNPLoading + WS.Index() * NumEig;
			pAFreq = AFreq + WS.Index();
			nSubSNP = WS.Count();
			// using thread thpool
			thpool.BatchWork(this, &CEigMix_SampleLoad::thread_loading, nSamp);
			// update
			WS.ProgressForward(WS.Count());
		}
	}
};

}


extern "C"
{

using namespace EIGMIX;


static void CPU_Flag()
{
	Rprintf("CPU capabilities:");
	#ifdef COREARRAY_SIMD_SSE2
		Rprintf(" Double-Precision SSE2");
	#endif
	#ifdef COREARRAY_SIMD_AVX
		Rprintf(" AVX");
	#endif
	Rprintf("\n");
}


/// Calculate EigMix GRM matrix
COREARRAY_DLL_LOCAL void CalcEigMixGRM(CdMatTri<double> &grm, int NumThread,
	bool Verbose)
{
	if (Verbose) CPU_Flag();
	CEigMix_AlgArith eigmix(MCWorkingGeno.Space());
	eigmix.Run(grm, NumThread, NULL, false, Verbose);
	vec_f64_mul(grm.Get(), grm.Size(), 2);
}


/// Calculate the eigenvalues and eigenvectors from EIGMIX matrix
COREARRAY_DLL_EXPORT SEXP gnrEigMix(SEXP EigenCnt, SEXP NumThread,
	SEXP ParamList, SEXP Verbose)
{
	const bool verbose = SEXP_Verbose(Verbose);
	int diag_adj = Rf_asLogical(RGetListElement(ParamList, "diagadj"));
	if (diag_adj == NA_LOGICAL)
		error("'diagadj' must be TRUE or FALSE.");
	int need_ibd = Rf_asLogical(RGetListElement(ParamList, "ibdmat"));
	if (need_ibd == NA_LOGICAL)
		error("'ibdmat' must be TRUE or FALSE.");

	COREARRAY_TRY

		// cache the genotype data
		CachingSNPData("Eigen-analysis", verbose);
		if (verbose) CPU_Flag();

		// the number of samples
		const R_xlen_t n = MCWorkingGeno.Space().SampleNum();
		int nEig = Rf_asInteger(EigenCnt);
		if (nEig < 0 || nEig > n) nEig = n;

		int nProtected = 1;
		SEXP ibd_mat, EigVal, EigVect, AFreq;
		ibd_mat = EigVal = EigVect = R_NilValue;
		AFreq = PROTECT(NEW_NUMERIC(MCWorkingGeno.Space().SNPNum()));

		if (1)
		{
			// the upper-triangle IBD matrix
			CdMatTri<double> IBD(n);
			CEigMix_AlgArith eigmix(MCWorkingGeno.Space());
			eigmix.Run(IBD, Rf_asInteger(NumThread), REAL(AFreq),
				diag_adj==TRUE, verbose);
			if (need_ibd)
			{
				ibd_mat = PROTECT(Rf_allocMatrix(REALSXP, n, n));
				nProtected ++;
				IBD.SaveTo(REAL(ibd_mat));
			}
			if (verbose)
				Rprintf("%s    Begin (eigenvalues and eigenvectors)\n", TimeToStr());
			vec_f64_mul(IBD.Get(), IBD.Size(), -1);
			nProtected += CalcEigen(IBD.Get(), n, nEig, "DSPEVX", EigVal, EigVect);
		} else {
			// the IBD matrix
			CdMat<double> IBD(n);
			CEigMix_AlgIndexing eigmix(MCWorkingGeno.Space());
			eigmix.Run(IBD, Rf_asInteger(NumThread), verbose);
			if (need_ibd)
			{
				ibd_mat = PROTECT(Rf_allocMatrix(REALSXP, n, n));
				nProtected ++;
				memcpy(REAL(ibd_mat), IBD.Get(), sizeof(double)*IBD.Size());
			}
			// to upper triangle matrix
			double *p = IBD.Get();
			for (R_xlen_t i=0; i < n; i++)
			{
				size_t m = n - i;
				memmove(p, IBD.Get() + i*n + i, sizeof(double)*m);
				p += m;
			}
			vec_f64_mul(IBD.Get(), n*(n+1)/2, -1);
			nProtected += CalcEigen(IBD.Get(), n, nEig, "DSPEVX", EigVal, EigVect);
		}

		if (verbose)
			Rprintf("%s    Done.\n", TimeToStr());

		// output
		PROTECT(rv_ans = NEW_LIST(4)); nProtected++;
		SET_ELEMENT(rv_ans, 0, EigVal);
		SET_ELEMENT(rv_ans, 1, EigVect);
		SET_ELEMENT(rv_ans, 2, AFreq);
		SET_ELEMENT(rv_ans, 3, ibd_mat);
		UNPROTECT(nProtected);

	COREARRAY_CATCH
}


/// Calculate the SNP loadings
COREARRAY_DLL_EXPORT SEXP gnrEigMixSNPLoading(SEXP EigenVal, SEXP EigenVect,
	SEXP AFreq, SEXP NumThread, SEXP Verbose)
{
	const bool verbose = SEXP_Verbose(Verbose);
	int LenEig = INTEGER(GET_DIM(EigenVect))[1];

	COREARRAY_TRY

		// cache the genotype data
		CachingSNPData("SNP Loading", verbose);

		// scale eigenvectors with eigenvalues
		SEXP EigVect = PROTECT(duplicate(EigenVect));
		{
			const size_t n = MCWorkingGeno.Space().SampleNum();
			for (int i=0; i < LenEig; i++)
			{
				vec_f64_mul(REAL(EigVect) + i*n, n,
					sqrt(1 / REAL(EigenVal)[i]));
			}
		}

		rv_ans = PROTECT(Rf_allocMatrix(REALSXP, LenEig,
			MCWorkingGeno.Space().SNPNum()));
		{
			CEigMix_SNPLoad Work(MCWorkingGeno.Space());
			Work.Run(REAL(rv_ans), REAL(AFreq), LenEig, REAL(EigVect),
				Rf_asInteger(NumThread), verbose);
		}
		if (verbose)
			Rprintf("%s    Done.\n", TimeToStr());
		UNPROTECT(2);

	COREARRAY_CATCH
}


/// Calculate the sample loadings from SNP loadings
COREARRAY_DLL_EXPORT SEXP gnrEigMixSampLoading(SEXP SNPLoadings,
	SEXP AFreq, SEXP NumThread, SEXP Verbose)
{
	const bool verbose = SEXP_Verbose(Verbose);
	int EigenCnt = INTEGER(GET_DIM(SNPLoadings))[0];

	COREARRAY_TRY

		// cache the genotype data
		CachingSNPData("Sample Loading", verbose);

		rv_ans = PROTECT(Rf_allocMatrix(REALSXP,
			MCWorkingGeno.Space().SampleNum(), EigenCnt));
		{
			CEigMix_SampleLoad Work(MCWorkingGeno.Space());
			Work.Run(REAL(rv_ans), EigenCnt, REAL(SNPLoadings),
				REAL(AFreq), Rf_asInteger(NumThread), verbose);
		}
		if (verbose)
			Rprintf("%s    Done.\n", TimeToStr());
		UNPROTECT(1);

	COREARRAY_CATCH
}

}

#endif  /* _HEADER_EIGMIX_ */
