// ===========================================================
//
// genPCA.cpp: Principal Component Analysis on GWAS
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

#include "genPCA.h"
#include <R_ext/Lapack.h>
#include <cmath>
#include <memory>
#include <algorithm>


namespace PCA
{

// using namespace
using namespace std;
using namespace CoreArray;
using namespace Vectorization;
using namespace GWAS;


// ---------------------------------------------------------------------
// PCA parameters

class CProdMat_AlgArith;

/// whether use Bayesian normalization
bool BayesianNormal = false;
/// The number of eigenvectors output
long OutputEigenDim = 32;
/// the pointer to the output buffer for SNP correlation and SNP/sample loadings
double *Out_Buffer = NULL;
/// the pointer to eigenvectors, snp loadings
double *In_EigenVect = NULL;
/// the pointer to the object storing the eigenvectors, snp loadings
CProdMat_AlgArith *_EigenVectBuf = NULL;
/// the pointer to allele frequency
double *In_AveFreq = NULL;

/// init mutex objects
void InitMutexObject()
{
	_Mutex = GDS_Parallel_InitMutex();
}
/// destroy mutex objects
void DoneMutexObject()
{
	GDS_Parallel_DoneMutex(_Mutex); _Mutex = NULL;
}



// ---------------------------------------------------------------------
// Vectorization Computing, algorithm 1 (Arithmetic)

void CProdMat_Base::ZeroFill()
{
	memset(GenoSum.Get(), 0, sizeof(C_Int32)*fM);
	memset(GenoNum.Get(), 0, sizeof(C_Int32)*fM);
}

// calculate average genotypes (save to GenoSum and GenoNum)
void CProdMat_Base::SummarizeGeno_SampxSNP(C_UInt8 *pGeno, size_t nSNP)
{
	C_Int32 *pS = GenoSum.Get();
	C_Int32 *pN = GenoNum.Get();
	for (size_t i=0; i < nSNP; i++)
	{
		pGeno = vec_u8_geno_count(pGeno, fN, *pS, *pN);
		pS++; pN++;
	}
	for (; nSNP < fM; nSNP++)
		*pS++ = *pN++ = 0;
}

// divide GenoSum by GenoNum (GenoSum / GenoNum)
void CProdMat_Base::DivideGeno()
{
	double *p = avg_geno.Get();
	C_Int32 *pS = GenoSum.Get();
	C_Int32 *pN = GenoNum.Get();
	size_t n = fM;

#if defined(COREARRAY_SIMD_AVX)
	const __m256d zero = _mm256_setzero_pd();
	for (; n >= 4; n-=4)
	{
		__m256d SD = _mm256_cvtepi32_pd(_mm_load_si128((__m128i const*)pS));
		pS += 4;
		__m256d ND = _mm256_cvtepi32_pd(_mm_load_si128((__m128i const*)pN));
		pN += 4;
		_mm256_store_pd(p, _mm256_and_pd(_mm256_div_pd(SD, ND),
			_mm256_cmp_pd(zero, ND, _CMP_LT_OQ)));
		p += 4;
	}
#elif defined(COREARRAY_SIMD_SSE2)
	const __m128d zero = _mm_setzero_pd();
	for (; n >= 4; n-=4)
	{
		__m128i S = _mm_load_si128((__m128i const*)pS);
		pS += 4;
		__m128i N = _mm_load_si128((__m128i const*)pN);
		pN += 4;

		__m128d SD = _mm_cvtepi32_pd(S);
		__m128d ND = _mm_cvtepi32_pd(N);
		_mm_store_pd(p, _mm_and_pd(_mm_div_pd(SD, ND), _mm_cmplt_pd(zero, ND)));
		p += 2;

		SD = _mm_cvtepi32_pd(_mm_shuffle_epi32(S, 0x0E));
		ND = _mm_cvtepi32_pd(_mm_shuffle_epi32(N, 0x0E));
		_mm_store_pd(p, _mm_and_pd(_mm_div_pd(SD, ND), _mm_cmplt_pd(zero, ND)));
		p += 2;
	}
#endif
	for (; n > 0; n--)
	{
		*p++ = (*pN > 0) ? ((double)*pS / *pN) : 0;
		pS ++; pN ++;
	}
}

// computing scalar vector
void CProdMat_Base::rsqrt_prod()
{
	double *p = avg_geno.Get();
	size_t n = fM;

#if defined(COREARRAY_SIMD_AVX)
	const __m256d half = _mm256_set1_pd(0.5);
	const __m256d zero = _mm256_set1_pd(0.0);
	const __m256d one  = _mm256_set1_pd(1.0);
	for (; n >= 4; n-=4, p+=4)
	{
		__m256d s = _mm256_mul_pd(_mm256_load_pd(p), half);
		__m256d m = _mm256_and_pd(_mm256_cmp_pd(zero, s, _CMP_LT_OQ),
			_mm256_cmp_pd(s, one, _CMP_LT_OQ));
		s = _mm256_mul_pd(s, _mm256_sub_pd(one, s));
		s = _mm256_div_pd(one, _mm256_sqrt_pd(s));
		_mm256_store_pd(p, _mm256_and_pd(s, m));
	}
#elif defined(COREARRAY_SIMD_SSE2)
	const __m128d half = _mm_set1_pd(0.5);
	const __m128d zero = _mm_set1_pd(0.0);
	const __m128d one  = _mm_set1_pd(1.0);
	for (; n >= 2; n-=2, p+=2)
	{
		__m128d s = _mm_mul_pd(_mm_load_pd(p), half);
		__m128d m = _mm_and_pd(_mm_cmplt_pd(zero, s), _mm_cmplt_pd(s, one));
		s = _mm_mul_pd(s, _mm_sub_pd(one, s));
		s = _mm_div_pd(one, _mm_sqrt_pd(s));
		_mm_store_pd(p, _mm_and_pd(s, m));
	}
#endif
	for (; n > 0; n--)
	{
		double s = (*p) * 0.5;
		*p++ = (0<s && s<1) ? (1.0 / sqrt(s*(1-s))) : 0;
	}
}

void CProdMat_Base::_Reset()
{
	GenoSum.Reset(fM);
	GenoNum.Reset(fM);
	avg_geno.Reset(fM);
}

void CProdMat_Base::_Clear()
{
	GenoSum.Clear();
	GenoNum.Clear();
	avg_geno.Clear();
}



void CProdMat_AlgArith::Reset(size_t n, size_t m)
{
#if defined(COREARRAY_SIMD_AVX)
	// 4-aligned array
	if (m & 0x03) m += 4 - (m & 0x03);
#elif defined(COREARRAY_SIMD_SSE2)
	// 2-aligned array
	if (m & 0x01) m += 2 - (m & 0x01);
#endif
	fGenotype.Reset(n * m);
	fN = n; fM = m;
	_Reset();
}

void CProdMat_AlgArith::Clear()
{
	_Clear();
	fGenotype.Clear();
}

// detect the effective value for BlockNumSNP
void CProdMat_AlgArith::PCA_Detect_BlockNumSNP(int nSamp)
{
	size_t Cache = GetOptimzedCache();
	BlockNumSNP = Cache / (sizeof(double)*nSamp);
	BlockNumSNP = (BlockNumSNP / 4) * 4;
	if (BlockNumSNP < 16) BlockNumSNP = 16;
}

// time-consuming function
void COREARRAY_CALL_ALIGN CProdMat_AlgArith::MulAdd(IdMatTri &Idx,
	size_t IdxCnt, double *pOut)
{
	for (; IdxCnt > 0; IdxCnt--, ++Idx)
	{
		double *p1 = base() + Idx.Row() * fM;
		double *p2 = base() + Idx.Column() * fM;
		size_t n = fM;

	#if defined(COREARRAY_SIMD_AVX)

		// fM can be divided by 4
		__m256d rv4_1, rv4_2;
		rv4_1 = rv4_2 = _mm256_setzero_pd();
		for (; n >= 16; n -= 16)
		{
			rv4_1 = _mm256_add_pd(rv4_1, _mm256_mul_pd(
				_mm256_load_pd(p1), _mm256_load_pd(p2)));
			rv4_2 = _mm256_add_pd(rv4_2, _mm256_mul_pd(
				_mm256_load_pd(p1 + 4), _mm256_load_pd(p2 + 4)));
			p1 += 8; p2 += 8;
			rv4_1 = _mm256_add_pd(rv4_1, _mm256_mul_pd(
				_mm256_load_pd(p1), _mm256_load_pd(p2)));
			rv4_2 = _mm256_add_pd(rv4_2, _mm256_mul_pd(
				_mm256_load_pd(p1 + 4), _mm256_load_pd(p2 + 4)));
			p1 += 8; p2 += 8;
		}

		__m256d rv4 = _mm256_add_pd(rv4_1, rv4_2);
		for (; n >= 4; n -= 4)
		{
			rv4 = _mm256_add_pd(rv4, _mm256_mul_pd(
				_mm256_load_pd(p1), _mm256_load_pd(p2)));
			p1 += 4; p2 += 4;
		}

		(*pOut++) += vec_avx_sum_f64(rv4);

	#elif defined(COREARRAY_SIMD_SSE2)

		// unroll loop
		__m128d rv2_1, rv2_2;
		rv2_1 = rv2_2 = _mm_setzero_pd();
		for (; n >= 8; n -= 8)
		{
			rv2_1 = _mm_add_pd(rv2_1, _mm_mul_pd(_mm_load_pd(p1),
				_mm_load_pd(p2)));
			rv2_2 = _mm_add_pd(rv2_2, _mm_mul_pd(_mm_load_pd(p1+2),
				_mm_load_pd(p2+2)));
			p1 += 4; p2 += 4;
			rv2_1 = _mm_add_pd(rv2_1, _mm_mul_pd(_mm_load_pd(p1),
				_mm_load_pd(p2)));
			rv2_2 = _mm_add_pd(rv2_2, _mm_mul_pd(_mm_load_pd(p1+2),
				_mm_load_pd(p2+2)));
			p1 += 4; p2 += 4;
		}

		__m128d rv2 = _mm_add_pd(rv2_1, rv2_2);
		for (; n >= 2; n -= 2)
		{
			rv2 = _mm_add_pd(rv2, _mm_mul_pd(_mm_load_pd(p1),
				_mm_load_pd(p2)));
			p1 += 2; p2 += 2;
		}

		(*pOut++) += vec_sum_f64(rv2);

	#else

		double rv = 0;
		// unroll loop
		for (; n >= 4; n -= 4)
		{
			rv += p1[0] * p2[0]; rv += p1[1] * p2[1];
			rv += p1[2] * p2[2]; rv += p1[3] * p2[3];
			p1 += 4; p2 += 4;
		}
		for (; n > 0; n--)
			rv += (*p1++) * (*p2++);
		(*pOut++) += rv;

	#endif
	}
}

// mean-adjusted genotypes (fGenotype - avg_geno)
void COREARRAY_CALL_ALIGN CProdMat_AlgArith::GenoSub()
{
	double *pGeno = fGenotype.Get();
	for (size_t num=fN; num > 0; num--, pGeno+=fM)
	{
		size_t n = fM;
		double *p = pGeno, *s = avg_geno.Get();
	#if defined(COREARRAY_SIMD_AVX)
		for (; n >= 4; n -= 4)
		{
			__m256d v = _mm256_sub_pd(_mm256_load_pd(p), _mm256_load_pd(s));
			_mm256_store_pd(p, v);
			p += 4; s += 4;
		}
	#endif
	#if defined(COREARRAY_SIMD_SSE2)
		for (; n >= 2; n -= 2)
		{
			__m128d v = _mm_sub_pd(_mm_load_pd(p), _mm_load_pd(s));
			_mm_store_pd(p, v);
			p += 2; s += 2;
		}
	#endif
		for (; n > 0; n--) *p++ -= (*s++);
	}
}

// variance-adjusted genotypes (fGenotype * avg_geno)
void COREARRAY_CALL_ALIGN CProdMat_AlgArith::GenoMul()
{
	double *pGeno = fGenotype.Get();
	for (size_t num=fN; num > 0; num--, pGeno+=fM)
	{
		size_t n = fM;
		double *p = pGeno, *s = avg_geno.Get();
	#if defined(COREARRAY_SIMD_AVX)
		for (; n >= 4; n -= 4)
		{
			__m256d v = _mm256_mul_pd(_mm256_load_pd(p), _mm256_load_pd(s));
			_mm256_store_pd(p, v);
			p += 4; s += 4;
		}
	#endif
	#if defined(COREARRAY_SIMD_SSE2)
		for (; n >= 2; n -= 2)
		{
			__m128d v = _mm_mul_pd(_mm_load_pd(p), _mm_load_pd(s));
			_mm_store_pd(p, v);
			p += 2; s += 2;
		}
	#endif
		for (; n > 0; n--) *p++ *= (*s++);
	}
}




// ================== PCA covariance matrix ==================

// -----------------------------------------------------------
// Exact PCA

class COREARRAY_DLL_LOCAL CExactPCA: protected CProdMat_AlgArith
{
private:
	CdBaseWorkSpace &Space;
	double *ptrCov;

	void thread_cov_outer(size_t i, size_t n)
	{
		IdMatTri I = Array_Thread_MatIdx[i];
		MulAdd(I, Array_Thread_MatCnt[i], ptrCov + I.Offset());
	}

public:
	CExactPCA(CdBaseWorkSpace &space): CProdMat_AlgArith(), Space(space)
		{ }

	/// run the algorithm
	void Run(CdMatTri<double> &Cov, int NumThread, bool Bayesian,
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
		ptrCov = Cov.Get();
		memset(ptrCov, 0, sizeof(double)*Cov.Size());

		// thread pool
		CThreadPoolEx<CExactPCA> thpool(NumThread);
		Array_SplitJobs(NumThread, nSamp, Array_Thread_MatIdx,
			Array_Thread_MatCnt);

		// genotypes (0, 1, 2 and NA)
		VEC_AUTO_PTR<C_UInt8> Geno(nSamp * BlockNumSNP);
		C_UInt8 *pGeno = Geno.Get();

		// genotype buffer, false for no memory buffer
		CGenoReadBySNP WS(NumThread, Space, BlockNumSNP, verbose ? -1 : 0, false);

		// for-loop
		WS.Init();
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

			// 1 / sqrt(p@j*(1-p@j))
			if (!Bayesian)
			{
				rsqrt_prod();
			} else {
				double *p = avg_geno.Get();
				C_Int32 *pSum = GenoSum.Get();
				C_Int32 *pNum = GenoNum.Get();
				for (size_t m=WS.Count(); m > 0; m--)
				{
					double s = ((*pSum++) + 1.0) / (2 * (*pNum++) + 2);
					*p++ = 1.0 / sqrt(s * (1 - s));
				}
			}

			// (G@ij - 2p@j) / sqrt(p@j*(1-p@j)), get the normalized genotypes
			GenoMul();

			// outer product, using thread thpool
			thpool.BatchWork(this, &CExactPCA::thread_cov_outer, NumThread);

			// update
			WS.ProgressForward(WS.Count());
		}
	}
};



// -----------------------------------------------------------
// Fast randomized PCA

class COREARRAY_DLL_LOCAL CRandomPCA
{
private:
	CdBaseWorkSpace &Space;

	size_t nSamp;  /// the number of selected samples
	size_t nSNP;   /// the number of selected SNPs
	double *AuxMat;  /// auxiliary matrix (G_i = nSamp X AuxDim)
	size_t AuxDim;  /// the number of columns
	int IterNum;  /// the number of iterations

	size_t hsize;  /// the number of columns
	VEC_AUTO_PTR<double> MatH;
	VEC_AUTO_PTR<double> LookupY;  /// normalized genotypes, Y

	VEC_AUTO_PTR<C_UInt8> Geno;  /// the genotype buffer (0, 1, 2 and NA)
	VEC_AUTO_PTR<double> Y_mc;  /// normalized genotypes at a site, with multiple cores
	VEC_AUTO_PTR<double> AuxMat_mc; /// Aux matrices, with multiple cores
	VEC_AUTO_PTR<double> MatT;  /// T = U_H^T * Y

	size_t iSNP;  /// the starting SNP index
	int iteration;  /// the iteration

	// auxiliary thread variables
	vector<size_t> thread_start, thread_length;

	void thread_lookup_y(size_t i, size_t n)
	{
		C_UInt8 *g = &Geno[i * nSamp];
		double *Y = &LookupY[(iSNP + i) * 4];
		for (; n > 0; n--)
		{
			C_Int32 sum, num;
			vec_u8_geno_count(g, nSamp, sum, num);
			g += nSamp;

			double avg = (num>0) ? (double)sum / num : 0;
			double p = avg * 0.5;
			double s = (0<p && p<1) ? 1.0 / sqrt(2*p*(1-p)) : 0;

			Y[0] = (0 - avg) * s; Y[1] = (1 - avg) * s;
			Y[2] = (2 - avg) * s; Y[3] = 0;
			Y += 4;
		}
	}

	void thread_Y_x_G_i(size_t i, size_t num)
	{
		C_UInt8 *pGeno = &Geno[i * nSamp];
		i += iSNP;
		for (; num > 0; num--, i++)
		{
			double *pH = &MatH[(AuxDim * iteration) + (i * hsize)];
			double *pY = &LookupY[i * 4];
			double *pA = AuxMat;

			for (size_t j=0; j < AuxDim; j++)
			{
				C_UInt8 *pG = pGeno;
				size_t n = nSamp;
				double sum = 0;

			#if defined(COREARRAY_SIMD_AVX)

				__m256d sum4 = _mm256_setzero_pd();
				for (; n >= 4; n-=4)
				{
					__m256d a = _mm256_loadu_pd(pA); pA += 4;

				#if defined(COREARRAY_SIMD_AVX2)
					__m256i k4 = _mm256_set_epi64x(pG[3], pG[2], pG[1], pG[0]);
					__m256d y = _mm256_i64gather_pd(pY, k4, sizeof(double));
				#else
					__m256d y = _mm256_set_pd(pY[pG[3]], pY[pG[2]],
						pY[pG[1]], pY[pG[0]]);
				#endif

					pG += 4;
					sum4 = _mm256_add_pd(sum4, _mm256_mul_pd(a, y));
				}
				sum = vec_avx_sum_f64(sum4);

			#elif defined(COREARRAY_SIMD_SSE2)

				__m128d sum2 = _mm_setzero_pd();
				for (; n >= 2; n-=2)
				{
					__m128d a = _mm_loadu_pd(pA); pA += 2;
					__m128d y = _mm_set_pd(pY[pG[1]], pY[pG[0]]); pG += 2;
					sum2 = _mm_add_pd(sum2, _mm_mul_pd(a, y));
				}
				sum = vec_sum_f64(sum2);

			#endif

				for (; n > 0; n--)
					sum += pY[*pG++] * (*pA++);
				*pH++ = sum;
			}

			pGeno += nSamp;
		}
	}

	void thread_YT_x_H_i(size_t i, size_t num)
	{
		size_t ii = thread_start[i];
		size_t n  = thread_length[i];
		double *pH = &MatH[AuxDim*iteration + (iSNP + ii)*hsize];
		double *pY = &LookupY[(iSNP + ii) * 4];

		for (; n > 0; n--, ii++)
		{
			// set Y_mc
			C_UInt8 *pG = &Geno[ii * nSamp];
			double *Y = &Y_mc[i * nSamp];
			for (size_t j=0; j < nSamp; j++)
			{
				Y[j] = pY[(*pG < 4) ? *pG : 3];
				pG ++;
			}

			double *pA = (i == 0) ? AuxMat : &AuxMat_mc[(i-1)*nSamp*AuxDim];
			for (size_t j=0; j < AuxDim; j++)
			{
				// pA += Y * pH[j]
				pA = vec_f64_addmul(pA, Y, nSamp, pH[j]);
			}

			pY += 4;
			pH += hsize;
		}
	}

	void thread_U_H_x_Y(size_t i, size_t num)
	{
		size_t ii = thread_start[i];
		size_t n  = thread_length[i];
		double *pH = &MatH[(iSNP + ii) * hsize];
		double *pY = &LookupY[(iSNP + ii) * 4];

		for (; n > 0; n--, ii++)
		{
			C_UInt8 *pG = &Geno[ii * nSamp];
			double *pT = &MatT[i * hsize * nSamp];

			for (size_t j=0; j < nSamp; j++)
			{
				double y = pY[(*pG < 4) ? *pG : 3];
				pG ++;
				pT = vec_f64_addmul(pT, pH, hsize, y);
			}

			pY += 4;
			pH += hsize;
		}
	}

	static void svd_vt(double *a, int m, int n, double *d)
	{
		int info = 0;
		double u=0, vt=0, w=0;
		vector<double> ds;
		if (!d)
			{ ds.resize(min(m, n)); d = &ds[0]; }

		int lwork = -1;
		F77_NAME(dgesvd)("N", "O", &m, &n, a, &m, d, &u, &m, &vt, &n,
			&w, &lwork, &info);
		if (info != 0)
			throw ErrCoreArray("LAPACK::DGESVD error (%d).", info);

		lwork = (int)w;
		vector<double> work(lwork);
		F77_NAME(dgesvd)("N", "O", &m, &n, a, &m, d, &u, &m, &vt, &n,
			&work[0], &lwork, &info);
		if (info != 0)
			throw ErrCoreArray("LAPACK::DGESVD error (%d).", info);
	}

public:
	/// constructor
	CRandomPCA(CdBaseWorkSpace &space, double *aux_mat, size_t aux_dim,
		int iter_num): Space(space)
	{
		nSamp = Space.SampleNum();
		nSNP  = Space.SNPNum();
		AuxMat = aux_mat;
		AuxDim = aux_dim;
		IterNum = iter_num;
		hsize = AuxDim * (IterNum + 1);
		MatH.Reset(nSNP * hsize);
		LookupY.Reset(nSNP * 4);
		iSNP = 0;
		iteration = 0;
	}

	/// run the algorithm
	SEXP Run(int NumThread, bool verbose)
	{
		if (NumThread < 1) NumThread = 1;
		size_t IncSNP = (256 / NumThread) * NumThread;
		if (IncSNP < 16) IncSNP = 16;
		if (verbose)
			Rprintf("%s\n", TimeToStr());

		Geno.Reset(nSamp * IncSNP);
		Y_mc.Reset(nSamp * NumThread);
		AuxMat_mc.Reset(nSamp * AuxDim * (NumThread - 1));
		thread_start.resize(NumThread);
		thread_length.resize(NumThread);

		// thread thpool
		CThreadPoolEx<CRandomPCA> thpool(NumThread);

		// genotype buffer, false for no memory buffer
		CGenoReadBySNP WS(NumThread, Space, IncSNP,
			verbose ? C_Int64(nSNP)*(2*IterNum + 1) : 0, false);

		// for-loop
		for (iteration=0; iteration <= IterNum; iteration++)
		{
			WS.Init();
			while (WS.Read(Geno.Get(), iSNP))
			{
				// need to initialize the variable LookupY
				if (iteration == 0)
					thpool.BatchWork(this, &CRandomPCA::thread_lookup_y, WS.Count());

				// update H_i = Y * G_i (G_i stored in AuxMat)
				thpool.BatchWork(this, &CRandomPCA::thread_Y_x_G_i, WS.Count());

				// update
				WS.ProgressForward(WS.Count());
			}

			// update G_{i+1} = Y^T * H_i / nSNP
			if (iteration < IterNum)
			{
				memset(AuxMat, 0, sizeof(double)*AuxDim*nSamp);
				memset(AuxMat_mc.Get(), 0, sizeof(double)*AuxMat_mc.Length());

				WS.Init();
				while (WS.Read(Geno.Get(), iSNP))
				{
					// update G_{i+1} = Y^T * H_i
					thpool.Split(NumThread, WS.Count(), &thread_start[0], &thread_length[0]);
					thpool.BatchWork(this, &CRandomPCA::thread_YT_x_H_i, NumThread);

					if (NumThread > 1)
					{
						// update G_{i+1} from AuxMat_mc
						size_t n = nSamp * AuxDim;
						for (size_t i=0; i < (size_t)NumThread-1; i++)
							vec_f64_add(AuxMat, &AuxMat_mc[i * n], n);
					}

					// update
					WS.ProgressForward(WS.Count());
				}

				// divide AuxMat by nSNP
				vec_f64_mul(AuxMat, AuxDim*nSamp, 1.0/nSNP);
			}
		}

		if (verbose)
			Rprintf("%s    Begin projecting genotypes and SVD\n", TimeToStr());

		// SVD MatH, get vt stored in MatH
		svd_vt(&MatH[0], hsize, nSNP, NULL);

		// T = U_H^T * Y
		MatT.Reset(hsize*nSamp*NumThread);
		memset(MatT.Get(), 0, sizeof(double)*MatT.Length());

		// for-loop
		WS.Init();
		while (WS.Read(Geno.Get(), iSNP))
		{
//			thpool.Split(1, inc_snp, &thread_start[0], &thread_length[0]);
//			CRandomPCA::thread_U_H_x_Y(0, 1);

			// update T = U_H^T * Y
			thpool.Split(1, WS.Count(), &thread_start[0], &thread_length[0]);
			thpool.BatchWork(this, &CRandomPCA::thread_U_H_x_Y, 1);

/*			if (NumThread > 1)
			{
				// update T = U_H^T * Y from MatT
				size_t n = hsize*nSamp;
				for (size_t i=1; i < (size_t)NumThread; i++)
					vec_f64_add(&MatT[0], &MatT[i * n], n);
			}
*/
		}

		vector<double> sigma(nSamp);
		svd_vt(&MatT[0], hsize, nSamp, &sigma[0]);

		SEXP rv_ans = PROTECT(NEW_LIST(2));
		{
			SEXP d = NEW_NUMERIC(nSamp);
			memcpy(REAL(d), &sigma[0], sizeof(double) * nSamp);
			SET_ELEMENT(rv_ans, 0, d);

			SEXP h = Rf_allocMatrix(REALSXP, hsize, nSamp);
			memcpy(REAL(h), &MatT[0], sizeof(double) * nSamp * hsize);
			SET_ELEMENT(rv_ans, 1, h);
		}

		UNPROTECT(1);

		if (verbose)
			Rprintf("%s    End\n", TimeToStr());

		return rv_ans;
	}
};



// ==================== SNP Correlations =====================

class COREARRAY_DLL_LOCAL CPCA_SNPCorr
{
private:
	CdBaseWorkSpace &Space;  ///< working genotypes
	VEC_AUTO_PTR<C_UInt8> Geno;  ///< genotypes (0, 1, 2 and NA)
	size_t nSamp;      ///< the number of samples
	size_t NumEigVal;  ///< the number of eigenvalues
	double *pEigVect;  ///< the pointer to eigenvectors
	double *pCorr;     ///< the pointer to correlation

	// Correlation
	inline static double SNP_PC_Corr(double *pX, C_UInt8 *pY, size_t n)
	{
		size_t m=0;
		double XY=0, X=0, XX=0, Y=0, YY=0, ans=R_NaN;
		while (n > 0)
		{
			if (*pY < 3)
			{
				XY += (*pX) * (*pY);
				X += *pX; XX += (*pX) * (*pX);
				Y += *pY; YY += (*pY) * (*pY);
				m ++;
			}
			pX++; pY++; n--;
		}
		if (m > 1)
		{
			double c1 = XX - X*X/m, c2 = YY - Y*Y/m, val = c1*c2;
			if (val > 0)
				ans = (XY - X*Y/m) / sqrt(val);
		}
		return ans;
	}

	void thread_corr(size_t i, size_t num)
	{
		C_UInt8 *pGeno = Geno.Get() + nSamp * i;
		double *pOut = pCorr + NumEigVal * i;
		for (; num > 0; num--, i++)
		{
			double *pEig = pEigVect;
			for (size_t j=0; j < NumEigVal; j++)
			{
				*pOut++ = SNP_PC_Corr(pEig, pGeno, nSamp);
				pEig += nSamp;
			}
			pGeno += nSamp;
		}
	}

public:
	/// constructor
	CPCA_SNPCorr(CdBaseWorkSpace &space): Space(space) { }

	/// run the algorithm
	void Run(double OutSNPCorr[], size_t NumEig, double EigVect[],
		int NumThread, bool verbose)
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

		// thread thpool
		CThreadPoolEx<CPCA_SNPCorr> thpool(NumThread);
		// genotypes (0, 1, 2 and NA)
		Geno.Reset(nSamp * nBlock);
		// genotype buffer, false for no memory buffer
		CGenoReadBySNP WS(NumThread, Space, nBlock, verbose ? -1 : 0, false);

		// for-loop
		WS.Init();
		while (WS.Read(Geno.Get()))
		{
			pCorr = OutSNPCorr + WS.Index() * NumEig;
			// using thread thpool
			thpool.BatchWork(this, &CPCA_SNPCorr::thread_corr, WS.Count());
			// update
			WS.ProgressForward(WS.Count());
		}
	}
};



// ================== SNP Loadings ==================

class COREARRAY_DLL_LOCAL CPCA_SNPLoad
{
private:
	CdBaseWorkSpace &Space;  ///< working genotypes
	VEC_AUTO_PTR<C_UInt8> Geno;  ///< genotypes (0, 1, 2 and NA)
	size_t nSamp;      ///< the number of samples
	size_t NumEigVal;  ///< the number of eigenvalues
	double *pEigVect;  ///< the pointer to eigenvectors
	double *pLoading;  ///< the pointer to correlation
	double *pGFreq;    ///< the pointer to allele frequency (times 2)
	double *pScale;    ///< the pointer to scale factor

	void thread_loading(size_t i, size_t num)
	{
		C_UInt8 *pGeno = Geno.Get() + nSamp * i;
		double *pOut = pLoading + NumEigVal * i;
		for (; num > 0; num--, i++)
		{
			double avg, scale;
			C_Int32 gSum, gNum;
			vec_u8_geno_count(pGeno, nSamp, gSum, gNum);
			if (gNum > 0)
			{
				avg = double(gSum) / gNum;
				if (!BayesianNormal)
				{
					scale = avg * 0.5;
					scale = ((0.0 < scale) && (scale < 1.0)) ?
						(1.0 / sqrt(scale*(1.0-scale))) : 0.0;
				} else {
					scale = double(gSum + 1) / (2*gNum + 2);
					scale = 1.0 / sqrt(scale*(1.0-scale));
				}
			} else {
				avg = scale = 0;
			}
			pGFreq[i] = avg;
			pScale[i] = scale;

			// zero filling
			memset(pOut, 0, sizeof(double)*NumEigVal);
			// dot product
			for (size_t j=0; j < nSamp; j++)
			{
				double g = (*pGeno < 3) ? ((*pGeno - avg) * scale) : 0.0;
				pGeno ++;
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
	CPCA_SNPLoad(CdBaseWorkSpace &space): Space(space) { }

	/// run the algorithm
	void Run(double OutSNPLoading[], double OutGFreq[], double OutScale[],
		size_t NumEig, double EigVect[], int NumThread, bool verbose)
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

		// thread thpool
		CThreadPoolEx<CPCA_SNPLoad> thpool(NumThread);
		// genotypes (0, 1, 2 and NA)
		Geno.Reset(nSamp * nBlock);
		// genotype buffer, false for no memory buffer
		CGenoReadBySNP WS(NumThread, Space, nBlock, verbose ? -1 : 0, false);

		// for-loop
		WS.Init();
		while (WS.Read(Geno.Get()))
		{
			pLoading = OutSNPLoading + WS.Index() * NumEig;
			pGFreq = OutGFreq + WS.Index();
			pScale = OutScale + WS.Index();
			// using thread thpool
			thpool.BatchWork(this, &CPCA_SNPLoad::thread_loading, WS.Count());
			// update
			WS.ProgressForward(WS.Count());
		}
	}
};



// ================== Sample Loadings ==================

class COREARRAY_DLL_LOCAL CPCA_SampleLoad
{
private:
	CdBaseWorkSpace &Space;  ///< working genotypes
	VEC_AUTO_PTR<C_UInt8> Geno;  ///< genotypes (0, 1, 2 and NA)
	size_t nSamp;      ///< the number of samples
	size_t nEigVal;    ///< the number of eigenvalues
	size_t nSubSNP;    ///< the number of genotypes in a block
	double *pLoading;  ///< the pointer to SNP eigenvectors
	double *pGFreq;    ///< the pointer to allele frequency (times 2)
	double *pScale;    ///< the pointer to scale factor
	double *pOutEig;   ///< the pointer to sample eigenvectors

	void thread_loading(size_t i, size_t num)
	{
		// for-loop each individual i
		for (; num > 0; num--, i++)
		{
			C_UInt8 *pGeno = Geno.Get() + i;
			double *pLoad = pLoading;
			for (size_t j=0; j < nSubSNP; j++)
			{
				double g = (*pGeno < 3) ? (*pGeno - pGFreq[j]) * pScale[j] : 0.0;
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
	CPCA_SampleLoad(CdBaseWorkSpace &space): Space(space) { }

	/// run the algorithm
	void Run(double OutSampLoad[], size_t NumEig, double SNPLoading[],
		double GFreq[], double Scale[], int NumThread, bool verbose)
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
		CThreadPoolEx<CPCA_SampleLoad> thpool(NumThread);
		// genotypes (0, 1, 2 and NA)
		Geno.Reset(nSamp * nBlock);
		// genotype buffer, false for no memory buffer
		CGenoReadBySNP WS(NumThread, Space, nBlock, verbose ? -1 : 0, false);

		// zero filling
		memset(OutSampLoad, 0, sizeof(double)*NumEig*nSamp);

		// for-loop
		WS.Init();
		while (WS.Read(Geno.Get()))
		{
			pLoading = SNPLoading + WS.Index() * NumEig;
			pGFreq = GFreq + WS.Index();
			pScale = Scale + WS.Index();
			nSubSNP = WS.Count();
			// using thread thpool
			thpool.BatchWork(this, &CPCA_SampleLoad::thread_loading, nSamp);
			// update
			WS.ProgressForward(WS.Count());
		}
	}
};




// -----------------------------------------------------------
// GCTA GRM Calculation with arithmetic algorithm

class COREARRAY_DLL_LOCAL CGCTA_AlgArith: protected PCA::CProdMat_AlgArith
{
private:
	CdBaseWorkSpace &Space;
	double *ptrCov;

	void thread_cov_outer(size_t i, size_t n)
	{
		IdMatTri I = Array_Thread_MatIdx[i];
		MulAdd(I, Array_Thread_MatCnt[i], ptrCov + I.Offset());
	}

public:
	CGCTA_AlgArith(CdBaseWorkSpace &space): CProdMat_AlgArith(), Space(space)
		{ }

	/// run the algorithm
	void Run(CdMatTri<double> &Cov, int NumThread, bool verbose)
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
		ptrCov = Cov.Get();
		memset(ptrCov, 0, sizeof(double)*Cov.Size());

		// denominator the number of valid sites
		CdMatTri<C_Int32> Denom(nSamp);
		memset(Denom.Get(), 0, sizeof(C_Int32)*Denom.Size());
		C_Int64 nLocus = 0;

		// thread pool
		CThreadPoolEx<CGCTA_AlgArith> thpool(NumThread);
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
			// 1 / sqrt(p@j*(1-p@j))
			rsqrt_prod();
			// (G@ij - 2p@j) / sqrt(p@j*(1-p@j)), get the normalized genotypes
			GenoMul();

			// denominator, missing values
			size_t nsnp = WS.Count();
			C_UInt8 *pG = pGeno;
			for (size_t i=0; i < nsnp; i++)
			{
				if ((0 < GenoSum[i]) && (GenoSum[i] < 2*GenoNum[i]))
				{
					nLocus ++;
					C_UInt8 *pGG = pG;
					for (size_t j=0; j < nSamp; j++)
					{
						if (*pG++ > 2)
						{
							C_Int32 *pp = Denom.Get();
							vec_i32_add(pp + j*(2*nSamp-j+1)/2, nSamp-j, 1);
							for (ssize_t k=(ssize_t)j - 1; k >= 0; k--)
								if (pGG[k] <= 2)
									pp[j + k*(2*nSamp-k-1)/2] ++;
						}
					}
				} else {
					pG += nSamp;
				}
			}

			// outer product, using thread thpool
			thpool.BatchWork(this, &CGCTA_AlgArith::thread_cov_outer, NumThread);
			// update
			WS.ProgressForward(WS.Count());
		}

		// normalized by the denominator
		double  *p = Cov.Get();
		C_Int32 *s = Denom.Get();
		for (size_t n=Cov.Size(); n > 0; n--)
			*p++ /= 2*C_Int64(nLocus - *s++);
	}
};

}


using namespace PCA;

extern "C"
{

/// get the eigenvalues and eigenvectors, return 'nProtect'
COREARRAY_DLL_LOCAL int CalcEigen(double *pMat, int n, int nEig,
	const char *EigMethod, SEXP &EigVal, SEXP &EigVect)
{
	if (nEig <= 0)
	{
		EigVal = EigVect = R_NilValue;
		return 0;
	}

	int nProtected = 0;

	// the method to compute eigenvalues and eigenvectors
	if (strcmp(EigMethod, "DSPEV") == 0)
	{
		vector<double> tmp_Work(n*3);
		vector<double> tmp_EigenVec(n*n);

		EigVal = PROTECT(NEW_NUMERIC(n));
		nProtected ++;

		// DSPEV -- computes all the eigenvalues and, optionally,
		// eigenvectors of a real symmetric matrix A in packed storage

		int info = 0;
		F77_NAME(dspev)("V", "L", &n, pMat, REAL(EigVal),
			&tmp_EigenVec[0], &n, &tmp_Work[0], &info);
		if (info != 0)
		{
			throw ErrCoreArray(
				"LAPACK::DSPEV error (%d), infinite or missing values in the genetic covariance matrix!",
				info);
		}

		// output eigenvalues
		vt<double>::Sub(REAL(EigVal), 0.0, REAL(EigVal), n);

		// output eigenvectors
		EigVect = PROTECT(Rf_allocMatrix(REALSXP, n, nEig));
		nProtected ++;

		for (size_t i=0; i < (size_t)nEig; i++)
		{
			memmove(REAL(EigVect) + n*i,
				&tmp_EigenVec[0] + n*i, sizeof(double)*n);
		}

	} else if (strcmp(EigMethod, "DSPEVX") == 0)
	{
		vector<double> tmp_Work(n*8);
		vector<int>    tmp_IWork(n*5);

		EigVal = PROTECT(NEW_NUMERIC(n));
		EigVect = PROTECT(Rf_allocMatrix(REALSXP, n, nEig));
		nProtected += 2;

		// DSPEVX -- compute selected eigenvalues and, optionally,
		// eigenvectors of a real symmetric matrix A in packed storage

		int IL = 1, IU = nEig, LDZ = n;
		double VL = 0, VU = 0;
		int M = 0;
		// it is suggested: ABSTOL is set to twice the underflow threshold, not zero
		double ABSTOL = 2 * F77_NAME(dlamch)("S");
		vector<int> ifail(n);
		int info = 0;

		F77_NAME(dspevx)("V", "I", "L", &n, pMat, &VL, &VU, &IL, &IU, &ABSTOL,
			&M, REAL(EigVal), REAL(EigVect), &LDZ,
			&tmp_Work[0], &tmp_IWork[0], &ifail[0], &info);
		if (info != 0)
		{
			throw ErrCoreArray(
				"LAPACK::DSPEVX error (%d), infinite or missing values in the genetic covariance matrix!",
				info);
		}

		// output eigenvalues
		double *p = REAL(EigVal);
		for (int i=0; i < nEig; i++) p[i] = -p[i];
		for (int i=nEig; i < n; i++) p[i] = R_NaN;
	} else
		throw ErrCoreArray("Unknown 'eigen.method'.");

	return nProtected;
}



// ========================================================================
// Principal Component Analysis (PCA)
// ========================================================================

/// Compute the eigenvalues and eigenvectors
COREARRAY_DLL_EXPORT SEXP gnrPCA(SEXP EigenCnt, SEXP Algorithm,
	SEXP NumThread, SEXP ParamList, SEXP Verbose)
{
	const bool verbose = SEXP_Verbose(Verbose);
	int nThread = Rf_asInteger(NumThread);
	if (nThread <= 0) nThread = 1;

	COREARRAY_TRY

		// cache the genotype data
		CachingSNPData("PCA", verbose);

		// PCA algorithm: exact and fast randomized
		if (strcmp(CHAR(STRING_ELT(Algorithm, 0)), "exact") == 0)
		{
			// set parameters
			PCA::BayesianNormal =
				Rf_asLogical(RGetListElement(ParamList, "bayesian")) == TRUE;

			// the number of samples
			const R_xlen_t n = MCWorkingGeno.Space().SampleNum();
			// the upper-triangle genetic covariance matrix
			CdMatTri<double> Cov(n);

			// Calculate the genetic covariace
			{
				CExactPCA pca(MCWorkingGeno.Space());
				pca.Run(Cov, nThread, BayesianNormal, verbose);
			}

			// Normalize
			double TraceXTX = Cov.Trace();
			double scale = double(n-1) / TraceXTX;
			vec_f64_mul(Cov.Get(), Cov.Size(), scale);
			double TraceVal = Cov.Trace();

			// ======== output ========

			int nProtected = 0;
			PROTECT(rv_ans = NEW_LIST(5));
			nProtected ++;

			SET_ELEMENT(rv_ans, 0, ScalarReal(TraceXTX));
			SET_ELEMENT(rv_ans, 4, ScalarReal(TraceVal));

			if (Rf_asLogical(RGetListElement(ParamList, "need.genmat")) == TRUE)
			{
				SEXP tmp;
				PROTECT(tmp = Rf_allocMatrix(REALSXP, n, n));
				nProtected ++;
				SET_ELEMENT(rv_ans, 1, tmp);
				Cov.SaveTo(REAL(tmp));
			}

			// ======== eigenvectors and eigenvalues ========

			if (Rf_asLogical(RGetListElement(ParamList, "genmat.only")) != TRUE)
			{
				if (verbose)
				{
					Rprintf("%s    Begin (eigenvalues and eigenvectors)\n",
						TimeToStr());
				}

				vt<double>::Sub(Cov.Get(), 0.0, Cov.Get(), Cov.Size());

				int nEig = Rf_asInteger(EigenCnt);
				if (nEig < 0)
					throw ErrCoreArray("Invalid 'eigen.cnt'.");
				if (nEig > n) nEig = n;

				SEXP EigVal  = R_NilValue;
				SEXP EigVect = R_NilValue;
				nProtected += CalcEigen(Cov.Get(), n, nEig,
					CHAR(STRING_ELT(RGetListElement(ParamList, "eigen.method"), 0)),
					EigVal, EigVect);
				SET_ELEMENT(rv_ans, 2, EigVal);
				SET_ELEMENT(rv_ans, 3, EigVect);
			}

			UNPROTECT(nProtected);

		} else if (strcmp(CHAR(STRING_ELT(Algorithm, 0)), "randomized") == 0)
		{
			CRandomPCA pca(MCWorkingGeno.Space(),
				REAL(RGetListElement(ParamList, "aux.mat")),
				Rf_asInteger(RGetListElement(ParamList, "aux.dim")),
				Rf_asInteger(RGetListElement(ParamList, "iter.num")));
			rv_ans = pca.Run(nThread, verbose);
		} else
			throw "Invalid 'algorithm'.";

		if (verbose)
			Rprintf("%s    Done.\n", TimeToStr());

	COREARRAY_CATCH
}


/// Calculate the SNP correlations
COREARRAY_DLL_EXPORT SEXP gnrPCACorr(SEXP LenEig, SEXP EigenVect,
	SEXP NumThread, SEXP _Verbose)
{
	const bool verbose = SEXP_Verbose(_Verbose);
	COREARRAY_TRY

		// cache the genotype data
		CachingSNPData("SNP Correlation", verbose);

		// output variable
		PROTECT(rv_ans = Rf_allocMatrix(REALSXP, Rf_asInteger(LenEig),
			MCWorkingGeno.Space().SNPNum()));
		{
			CPCA_SNPCorr Work(MCWorkingGeno.Space());
			Work.Run(REAL(rv_ans), Rf_asInteger(LenEig), REAL(EigenVect),
				Rf_asInteger(NumThread), verbose);
		}
		if (verbose)
			Rprintf("%s    Done.\n", TimeToStr());
		UNPROTECT(1);

	COREARRAY_CATCH
}


/// Calculate the SNP loadings
COREARRAY_DLL_EXPORT SEXP gnrPCASNPLoading(SEXP EigenVal, SEXP EigenVect,
	SEXP TraceXTX, SEXP NumThread, SEXP Bayesian, SEXP _Verbose)
{
	const bool verbose = SEXP_Verbose(_Verbose);
	int LenEig = INTEGER(GET_DIM(EigenVect))[1];

	COREARRAY_TRY

		// cache the genotype data
		CachingSNPData("SNP Loading", verbose);

		// scale eigenvectors with eigenvalues
		SEXP EigVect = PROTECT(duplicate(EigenVect));
		{
			const size_t n = MCWorkingGeno.Space().SampleNum();
			const double Scale = double(n - 1) / Rf_asReal(TraceXTX);
			for (int i=0; i < LenEig; i++)
			{
				vec_f64_mul(REAL(EigVect) + i*n, n,
					sqrt(Scale / REAL(EigenVal)[i]));
			}
		}

		PCA::BayesianNormal = (Rf_asLogical(Bayesian) == TRUE);
		const size_t n = MCWorkingGeno.Space().SNPNum();
		PROTECT(rv_ans = NEW_LIST(3));

		SEXP loading, afreq, scale;
		PROTECT(loading = Rf_allocMatrix(REALSXP, LenEig, n));
		SET_ELEMENT(rv_ans, 0, loading);

		PROTECT(afreq = NEW_NUMERIC(n));
		SET_ELEMENT(rv_ans, 1, afreq);

		PROTECT(scale = NEW_NUMERIC(n));
		SET_ELEMENT(rv_ans, 2, scale);

		{
			CPCA_SNPLoad Work(MCWorkingGeno.Space());
			Work.Run(REAL(loading), REAL(afreq), REAL(scale),
				LenEig, REAL(EigVect),
				Rf_asInteger(NumThread), verbose);
		}
		if (verbose)
			Rprintf("%s    Done.\n", TimeToStr());

		UNPROTECT(5);

	COREARRAY_CATCH
}


/// Calculate the sample loadings from SNP loadings
COREARRAY_DLL_EXPORT SEXP gnrPCASampLoading(SEXP EigenCnt, SEXP SNPLoadings,
	SEXP AvgFreq, SEXP Scale, SEXP NumThread, SEXP _Verbose)
{
	const bool verbose = SEXP_Verbose(_Verbose);
	COREARRAY_TRY

		// cache the genotype data
		CachingSNPData("Sample Loading", verbose);

		PROTECT(rv_ans = Rf_allocMatrix(REALSXP,
			MCWorkingGeno.Space().SampleNum(), Rf_asInteger(EigenCnt)));
		{
			CPCA_SampleLoad Work(MCWorkingGeno.Space());
			Work.Run(REAL(rv_ans), Rf_asInteger(EigenCnt), REAL(SNPLoadings),
				REAL(AvgFreq), REAL(Scale), Rf_asInteger(NumThread), verbose);
		}
		if (verbose)
			Rprintf("%s    Done.\n", TimeToStr());
		UNPROTECT(1);

	COREARRAY_CATCH
}



// =======================================================================
// Genetic relationship matrix
// =======================================================================

/// Compute the eigenvalues and eigenvectors
COREARRAY_DLL_EXPORT SEXP gnrGRM(SEXP _NumThread, SEXP _Method, SEXP _Verbose)
{
	const int nThread  = Rf_asInteger(_NumThread);
	const char *Method = CHAR(STRING_ELT(_Method, 0));
	const bool verbose = SEXP_Verbose(_Verbose);

	COREARRAY_TRY

		// ======== To cache the genotype data ========
		CachingSNPData("GRM Calculation", verbose);

		// ======== The calculation of genetic covariance matrix ========

		// the number of samples
		const R_xlen_t n = MCWorkingGeno.Space().SampleNum();

		// set internal parameters
		PCA::CProdMat_AlgArith::PCA_Detect_BlockNumSNP(n);

		// the upper-triangle IBD matrix
		CdMatTri<double> IBD(n);

		if (strcmp(Method, "Eigenstrat") == 0)
		{
			CExactPCA pca(MCWorkingGeno.Space());
			pca.Run(IBD, nThread, false, verbose);
			// normalize
			double TraceXTX = IBD.Trace();
			double scale = double(n-1) / TraceXTX;
			vt<double, av16Align>::Mul(IBD.Get(), IBD.Get(), scale, IBD.Size());

		} else if (strcmp(Method, "GCTA") == 0)
		{
			CachingSNPData("GCTA GRM", verbose);
			CGCTA_AlgArith GCTA(MCWorkingGeno.Space());
			GCTA.Run(IBD, nThread, verbose);
		} else if (strcmp(Method, "EIGMIX") == 0)
		{
			extern void CalcEigMixGRM(CdMatTri<double> &grm, int NumThread,
				bool Verbose);
			// cache the genotype data
			CachingSNPData("EIGMIX", verbose);
			CalcEigMixGRM(IBD, nThread, verbose);
		} else if (strcmp(Method, "IndivBeta") == 0)
		{
			extern SEXP CalcIndivBetaGRM(CdMatTri<double> &grm, int NumThread,
				bool Verbose);
			// cache the genotype data
			CachingSNPData("Individual Beta", verbose);
			CalcIndivBetaGRM(IBD, nThread, verbose);
		} else
			throw ErrCoreArray("Invalid 'method'!");

		// Output
		PROTECT(rv_ans = Rf_allocMatrix(REALSXP, n, n));
		IBD.SaveTo(REAL(rv_ans));
		UNPROTECT(1);

	COREARRAY_CATCH
}

}
