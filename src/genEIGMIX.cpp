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

class COREARRAY_DLL_LOCAL CEigMix_AlgArith: protected PCA::CPCAMat_AlgArith
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
	CEigMix_AlgArith(CdBaseWorkSpace &space): CPCAMat_AlgArith(), Space(space)
		{ }

	/// run the algorithm
	void Run(CdMatTri<double> &IBD, int NumThread, bool DiagAdj, bool verbose)
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
			double *p = base();
			for (size_t i=0; i < nSamp; i++)
			{
				double *pp = p; p += M();
				size_t m = WS.Count();
				for (size_t j=0; j < m; j++)
				{
					C_UInt8 g = pGeno[nSamp*j + i];
					*pp ++ = (g <= 2) ? g : avg_geno[j];
				}
				// zero fill the left part
				for (; m < fM; m++) *pp ++ = 0;
			}
			// G@ij - 2\bar{p}_l
			GenoSub();

			// denominator
			size_t nsnp = WS.Count();
			C_UInt8 *pG = pGeno;
			for (size_t i=0; i < nsnp; i++)
			{
				double denom = 2 * avg_geno[i] * (1 - 0.5*avg_geno[i]);
				SumDenominator += denom;
				if (GenoNum[i] < nSamp)
				{
					for (size_t j=0; j < nSamp; j++, pG++)
					{
						if (*pG == 1)
						{
							DiagAdjVal[j] ++;
						} else if (*pG > 2)
						{
							for (size_t k=j; k < nSamp; k++)
								Denom.Get()[k + j*(2*nSamp-j-1)/2] += denom;
						}
					}
				} else {
					for (size_t j=0; j < nSamp; j++)
						if (*pG++ == 1) DiagAdjVal[j] ++;
				}
			}

			// outer product, using thread thpool
			thpool.BatchWork(this, &CEigMix_AlgArith::thread_cov_outer, NumThread);
			// update
			WS.ProgressForward(WS.Count());
		}

		// finally
		if (DiagAdj)
		{
			for (size_t i=0; i < nSamp; i++)
				IBD.At(i, i) -= DiagAdjVal[i];
		}
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

			// update denominator
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

			// update denominator
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

}


extern "C"
{
using namespace EIGMIX;

/// to compute the eigenvalues and eigenvectors
COREARRAY_DLL_EXPORT SEXP gnrEIGMIX(SEXP EigenCnt, SEXP NumThread,
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

		// the number of samples
		const R_xlen_t n = MCWorkingGeno.Space().SampleNum();
		int nEig = Rf_asInteger(EigenCnt);
		if (nEig < 0 || nEig > n) nEig = n;

		int nProtected = 0;
		SEXP ibd_mat, EigVal, EigVect;
		ibd_mat = EigVal = EigVect = R_NilValue;

		if (1)
		{
			// the upper-triangle IBD matrix
			CdMatTri<double> IBD(n);
			CEigMix_AlgArith eigmix(MCWorkingGeno.Space());
			eigmix.Run(IBD, Rf_asInteger(NumThread), diag_adj==TRUE, verbose);
			if (need_ibd)
			{
				ibd_mat = PROTECT(Rf_allocMatrix(REALSXP, n, n));
				nProtected ++;
				IBD.SaveTo(REAL(ibd_mat));
			}
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

		// output
		PROTECT(rv_ans = NEW_LIST(3)); nProtected++;
		SET_ELEMENT(rv_ans, 0, EigVal);
		SET_ELEMENT(rv_ans, 1, EigVect);
		SET_ELEMENT(rv_ans, 2, ibd_mat);
		UNPROTECT(nProtected);

	COREARRAY_CATCH
}

}

#endif  /* _HEADER_EIGMIX_ */
