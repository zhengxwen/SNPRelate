// ===========================================================
//
// genPCA.cpp: Principal component analysis on GWAS
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


#ifndef _HEADER_PCA_
#define _HEADER_PCA_

// CoreArray library header
#include <dGenGWAS.h>
#include <dVect.h>
#include <R_ext/Lapack.h>

// Standard library header
#include <cmath>
#include <memory>
#include <algorithm>


#ifdef COREARRAY_SIMD_SSE
#include <xmmintrin.h>
#endif
#ifdef COREARRAY_SIMD_SSE2
#include <emmintrin.h>
#endif


#ifndef COREARRAY_CALL_ALIGN16_ARG
#   define COREARRAY_CALL_ALIGN16_ARG    COREARRAY_CALL_ALIGN
#endif


namespace PCA
{
	// using namespace
	using namespace std;
	using namespace CoreArray;
	using namespace CoreArray::Vectorization;
	using namespace GWAS;


	// ---------------------------------------------------------------------
	// Vector Computing

	// Vectorization

	#define __ALIGNED16__ __attribute__((aligned(16)))

	template<typename tfloat,
		bool SSE = (FlagVectorization >= vtSSE),
		bool SSE2 = (FlagVectorization >= vtSSE2)
	> class CPCAMat
	{
	public:
		static const bool SSEFlag = SSE;
		static const bool SSE2Flag = SSE2;

		CPCAMat() { fN = fM = 0; }
		CPCAMat(size_t n, size_t m) { Reset(n, m); }

		void Reset(size_t n, size_t m)
		{
			fBuf.resize(n * m);
			fN = n; fM = m;
		}

		void COREARRAY_CALL_ALIGN16_ARG MulAdd(IdMatTri &Idx, size_t IdxCnt,
			size_t ArrayLength, tfloat *OutBuf)
		{
			for (; IdxCnt > 0; IdxCnt--, ++Idx)
			{
				tfloat *p1 = base() + Idx.Row() * fM;
				tfloat *p2 = base() + Idx.Column() * fM;
				tfloat rv = 0;
				// unroll loop
				size_t k = ArrayLength;
				while (k >= 4)
				{
					rv += p1[0] * p2[0]; rv += p1[1] * p2[1];
					rv += p1[2] * p2[2]; rv += p1[3] * p2[3];
					p1 += 4; p2 += 4; k -= 4;
				}
				while (k > 0)
				{
					rv += (*p1++) * (*p2++);
					k--;
				}
				(*OutBuf++) += rv;
			}
		}

		inline size_t N() const { return fN; }
		inline size_t M() const { return fM; }
		inline tfloat* base() { return &fBuf[0]; }

	private:
		vector<tfloat> fBuf;
		size_t fN, fM;
	};

	#ifdef COREARRAY_SIMD_SSE
	template<bool SSE2> class CPCAMat<float, true, SSE2>
	{
	public:
		static const bool SSEFlag = true;
		static const bool SSE2Flag = SSE2;

		CPCAMat() { fN = fM = 0; }
		CPCAMat(size_t n, size_t m) { Reset(n, m); }

		void Reset(size_t n, size_t m)
		{
			if (m & 0x03) m += 4 - (m & 0x03);
			fBuf.Reset(n * m);
			fN = n; fM = m;
		}

		void COREARRAY_CALL_ALIGN16_ARG MulAdd(IdMatTri &Idx, size_t IdxCnt,
			size_t ArrayLength, float *OutBuf)
		{
			for (; IdxCnt > 0; IdxCnt--, ++Idx)
			{
				float *p1 = base() + Idx.Row() * fM;
				float *p2 = base() + Idx.Column() * fM;
				// unroll loop
				float rv4[4] __ALIGNED16__ = {0, 0, 0, 0};
				size_t k = ArrayLength;
				__m128 rv128 = _mm_load_ps(rv4);

				while (k >= 16)
				{
					rv128 = _mm_add_ps(rv128, _mm_mul_ps(_mm_load_ps(&p1[0]),
						_mm_load_ps(&p2[0])));
					rv128 = _mm_add_ps(rv128, _mm_mul_ps(_mm_load_ps(&p1[4]),
						_mm_load_ps(&p2[4])));
					rv128 = _mm_add_ps(rv128, _mm_mul_ps(_mm_load_ps(&p1[8]),
						_mm_load_ps(&p2[8])));
					rv128 = _mm_add_ps(rv128, _mm_mul_ps(_mm_load_ps(&p1[12]),
						_mm_load_ps(&p2[12])));
					p1 += 16; p2 += 16; k -= 16;
				}
				while (k >= 4)
				{
					rv128 = _mm_add_ps(rv128, _mm_mul_ps(_mm_load_ps(p1),
						_mm_load_ps(p2)));
					p1 += 4; p2 += 4; k -= 4;
				}
				_mm_store_ps(rv4, rv128);
				rv4[0] += rv4[1] + rv4[2] + rv4[3];
				while (k > 0)
				{
					rv4[0] += (*p1++) * (*p2++);
					k--;
				}
				(*OutBuf++) += rv4[0];
			}
		}

		inline size_t N() const { return fN; }
		inline size_t M() const { return fM; }
		inline float* base() { return fBuf.get(); }
	private:
		TdAlignPtr<float> fBuf;
		size_t fN, fM;
	};
	#endif // COREARRAY_SIMD_SSE

	#ifdef COREARRAY_SIMD_SSE2
	template<> class CPCAMat<double, true, true>
	{
	public:
		static const bool SSEFlag = true;
		static const bool SSE2Flag = true;

		CPCAMat() { fN = fM = 0; }
		CPCAMat(size_t n, size_t m) { Reset(n, m); }

		void Reset(size_t n, size_t m)
		{
			if (m & 0x01) m += 2 - (m & 0x01);
			fBuf.Reset(n * m);
			fN = n; fM = m;
		}

		void COREARRAY_CALL_ALIGN16_ARG MulAdd(IdMatTri &Idx, size_t IdxCnt,
			size_t ArrayLength, double *OutBuf)
		{
			for (; IdxCnt > 0; IdxCnt--, ++Idx)
			{
				double *p1 = base() + Idx.Row() * fM;
				double *p2 = base() + Idx.Column() * fM;
				// unroll loop
				double rv2[2] __ALIGNED16__ = {0, 0};
				size_t k = ArrayLength;

				__m128d rv128 = _mm_load_pd(rv2);
				while (k >= 8)
				{
					rv128 = _mm_add_pd(rv128, _mm_mul_pd(_mm_load_pd(&p1[0]),
						_mm_load_pd(&p2[0])));
					rv128 = _mm_add_pd(rv128, _mm_mul_pd(_mm_load_pd(&p1[2]),
						_mm_load_pd(&p2[2])));
					rv128 = _mm_add_pd(rv128, _mm_mul_pd(_mm_load_pd(&p1[4]),
						_mm_load_pd(&p2[4])));
					rv128 = _mm_add_pd(rv128, _mm_mul_pd(_mm_load_pd(&p1[6]),
						_mm_load_pd(&p2[6])));
					p1 += 8; p2 += 8; k -= 8;
				}

				while (k >= 2)
				{
					rv128 = _mm_add_pd(rv128, _mm_mul_pd(_mm_load_pd(p1),
						_mm_load_pd(p2)));
					p1 += 2; p2 += 2; k -= 2;
				}
				_mm_store_pd(rv2, rv128);
				rv2[0] += rv2[1];
				if (k > 0) rv2[0] += (*p1) * (*p2);
				(*OutBuf++) += rv2[0];
			}
		}

		inline size_t N() const { return fN; }
		inline size_t M() const { return fM; }
		inline double* base() { return fBuf.get(); }
	private:
		TdAlignPtr<double> fBuf;
		size_t fN, fM;
	};
	#endif // COREARRAY_SIMD_SSE2



	// ---------------------------------------------------------------------
	// PCA parameters

	/// whether use Bayesian normalization
	bool BayesianNormal = false;
	/// The number of eigenvectors output
	long OutputEigenDim = 16;
	/// the pointer to the output buffer for SNP correlation and SNP/sample loadings
	double *Out_Buffer = NULL;
	/// the pointer to eigenvectors, snp loadings
	double *In_EigenVect = NULL;
	/// the pointer to the object storing the eigenvectors, snp loadings
	CPCAMat<double> *_EigenVectBuf = NULL;
	/// the pointer to allele frequency
	double *In_AveFreq = NULL;

	/// detect the effective value for BlockNumSNP
	void AutoDetectSNPBlockSize(int nSamp, bool Detect=true)
	{
		if (Detect)
		{
			C_UInt64 L2Cache = GDS_Mach_GetCPULevelCache(2);
			C_UInt64 L3Cache = GDS_Mach_GetCPULevelCache(3);
			C_UInt64 Cache = (L2Cache > L3Cache) ? L2Cache : L3Cache;
			if ((C_Int64)Cache <= 0) Cache = 1024*1024; // 1M
			BlockNumSNP = (Cache - 8*1024) / (sizeof(double)*nSamp);
		}
		BlockNumSNP = (BlockNumSNP / 4) * 4;
		if (BlockNumSNP < 16) BlockNumSNP = 16;
	}

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
	// ---------------------------------------------------------------------

	/// Thread variables
	const int N_MAX_THREAD = 256;
	IdMatTri PCA_Thread_MatIdx[N_MAX_THREAD];
	C_Int64 PCA_Thread_MatCnt[N_MAX_THREAD];

	// ================== PCA covariate matrix ==================

	/// The temparory variables used in computing
	static vector<int> PCA_gSum, PCA_gNum;
	/// The mean-adjusted genotype buffers
	static CPCAMat<double> PCA_Mat;
	/// The temporary variable
	static TdAlignPtr<double> tmpBuf;

	/// Convert the raw genotypes to the mean-adjusted genotypes
	static void _Do_PCA_ReadBlock(C_UInt8 *GenoBuf, long Start, long SNP_Cnt,
		void* Param)
	{
		// init ...
		const int n = MCWorkingGeno.Space.SampleNum();
		memset(&PCA_gSum[0], 0, SNP_Cnt*sizeof(int));
		memset(&PCA_gNum[0], 0, SNP_Cnt*sizeof(int));

		C_UInt8 *p;
		double *pf, *pp;

		// calculate the averages of genotypes for each SelSNP
		// store the values in PCA_gSum
		p = GenoBuf;
		pp = PCA_Mat.base();
		for (long iSample=0; iSample < n; iSample++)
		{
			pf = pp; pp += PCA_Mat.M();
			for (long iSNP=0; iSNP < SNP_Cnt; iSNP++)
			{
				if (*p < 3)
				{
					PCA_gSum[iSNP] += *p;
					PCA_gNum[iSNP]++;
				}
				*pf++ = *p++;
			}
		}

		// PCA_gSum = 2 * \bar{p}_l
		for (long iSNP=0; iSNP < SNP_Cnt; iSNP++)
		{
			double &fv = tmpBuf.get()[iSNP];
			long gN = PCA_gNum[iSNP];
			if (gN > 0)
				fv = (double)PCA_gSum[iSNP] / gN;
			else
				fv = 0;
		}

		// Averaging and set missing values to ZERO
		pp = PCA_Mat.base(); // G@ij - 2\bar{p}_l
		for (long iSample=0; iSample < n; iSample++)
		{
			vt<double, av16Align>::Sub(pp, pp, tmpBuf.get(), SNP_Cnt);
			pp += PCA_Mat.M();
		}
		pf = tmpBuf.get(); // 1/sqrt(p@j*(1-p@j))
		for (long iSNP=0; iSNP < SNP_Cnt; iSNP++, pf++)
		{
			double scale;
			if (BayesianNormal)
			{
				long gN = PCA_gNum[iSNP];
				scale = (gN > 0) ? ((PCA_gSum[iSNP] + 1.0) / (2*gN + 2)) : 0;
			} else
				scale = (*pf) * 0.5;
			if (0.0 < scale && scale < 1.0)
				*pf = 1.0 / sqrt(scale*(1.0-scale));
			else
				*pf = 0;
		}
		// (G@ij - 2p@j) / sqrt(p@j*(1-p@j))
		pp = PCA_Mat.base();
		for (long iSample=0; iSample < n; iSample++)
		{
			vt<double, av16Align>::Mul(pp, pp, tmpBuf.get(), SNP_Cnt);
			pp += PCA_Mat.M();
		}

		// missing values
		p = GenoBuf; pp = PCA_Mat.base();
		for (long iSample=0; iSample < n; iSample++)
		{
			pf = pp; pp += PCA_Mat.M();
			for (long iSNP=0; iSNP < SNP_Cnt; iSNP++)
			{
				if (*p > 2) *pf = 0;
				p++; pf++;
			}
		}
	}

	/// Compute the covariate matrix
	static void _Do_PCA_ComputeCov(int ThreadIndex, long Start, long SNP_Cnt,
		void* Param)
	{
		double *base = (double*)Param;
		IdMatTri I = PCA_Thread_MatIdx[ThreadIndex];
		PCA_Mat.MulAdd(I, PCA_Thread_MatCnt[ThreadIndex], SNP_Cnt,
			base + I.Offset());
	}

    /// Calculate the genetic covariace
	void DoCovCalculate(CdMatTri<double> &PublicCov, int NumThread,
		const char *Info, bool verbose)
	{
		// Initialize ...
		PCA_gSum.resize(BlockNumSNP);
		PCA_gNum.resize(BlockNumSNP);
		PCA_Mat.Reset(PublicCov.N(), BlockNumSNP);
		tmpBuf.Reset(PCA_Mat.M());
		memset(PublicCov.get(), 0, sizeof(double)*PublicCov.Size());

		MCWorkingGeno.Progress.Info = Info;
		MCWorkingGeno.Progress.Show() = verbose;
		MCWorkingGeno.InitParam(true, true, BlockNumSNP);

		MCWorkingGeno.SplitJobs(NumThread, PublicCov.N(), PCA_Thread_MatIdx,
			PCA_Thread_MatCnt);
		MCWorkingGeno.Run(NumThread, &_Do_PCA_ReadBlock, &_Do_PCA_ComputeCov,
			PublicCov.get());
	}


	// ================== SNP coefficients ==================

	// Correlation
	static double SNP_PC_Corr(double *pX, C_UInt8 *pY, long n)
	{
		long m = 0;
		double XY, X, XX, Y, YY;

		XY = X = XX = Y = YY = 0;
		while (n > 0)
		{
			if (*pY < 3)
			{
				XY += (*pX) * (*pY);
				X += *pX; XX += (*pX) * (*pX);
				Y += *pY; YY += (*pY) * (*pY);
				m++;
			}
			pX++; pY++; n--;
		}
		if (m > 1)
		{
			double c1 = XX-X*X/m, c2 = YY-Y*Y/m, val = c1*c2;
			if (val > 0)
				return (XY-X*Y/m) / sqrt(val);
			else
				return R_NaN;
		} else
			return R_NaN;
	}

	/// The thread entry for the calculation of SNP coefficients
	void Entry_SNPCorrCalc(PdThread Thread, int ThreadIndex, void *Param)
	{
		// The number of working samples
		const long n = MCWorkingGeno.Space.SampleNum();
		vector<C_UInt8> GenoBlock(n * BlockNumSNP);

		long _SNPstart, _SNPlen;
		while (RequireWork(&GenoBlock[0], _SNPstart, _SNPlen, false))
		{
			C_UInt8 *pGeno = &GenoBlock[0];
			double *Out = Out_Buffer + (_SNPstart*OutputEigenDim);
			for (long iSNP=0; iSNP < _SNPlen; iSNP++)
			{
				// SNP coefficients
				double *pEig = In_EigenVect;
				for (long iEig=0; iEig < OutputEigenDim; iEig++)
				{
					*Out++ = SNP_PC_Corr(pEig, pGeno, n);
					pEig += n;
				}
				pGeno += n;
			}
			// Update progress
			{
				TdAutoMutex _m(_Mutex);
				MCWorkingGeno.Progress.Forward(_SNPlen);
			}
		}
	}

	/// Calculate the SNP coefficients
	void DoSNPCoeffCalculate(int nEig, double *EigenVect, double *out_snpcorr,
		int NumThread, bool verbose, const char *Info)
	{
		// initialize mutex objects
		InitMutexObject();

		// Initialize progress information
		MCWorkingGeno.Progress.Info = Info;
		MCWorkingGeno.Progress.Show() = verbose;
		MCWorkingGeno.Progress.Init(MCWorkingGeno.Space.SNPNum());
		SNPStart = 0;
		OutputEigenDim = nEig;
		Out_Buffer = out_snpcorr; In_EigenVect = EigenVect;

		// Threads
		GDS_Parallel_RunThreads(Entry_SNPCorrCalc, NULL, NumThread);

		// destory the mutex objects
		DoneMutexObject();
	}


	// ================== SNP loadings ==================

	/// Normalize the genotype
	/** \tparam tfloat    template type, it should be floating number type
	 *  \param Geno   genotype: 0 -- BB, 1 -- AB, 2 -- AA, 3 -- missing
	 *  \param NormalizedGeno  normalized genotype
	 *  \param NumGeno    the number of total genotype
	 *  \param AdjNormal  whether or not use Bayerian correction factor
	**/
	inline static void NormalizeGeno(C_UInt8 *Geno, double *NormalizedGeno, long NumGeno)
	{
		// the sum of genotype, and the number of non-missing snps
		long n;
		long sum = GENO_Get_Sum_ValidNumSNP(Geno, NumGeno, &n);
		if (n > 0)
		{
			double ave = double(sum) / n;
			double scale = BayesianNormal ? ((sum+1.0) / (2*n+2)) : (ave/2);
			scale = ((0.0 < scale) && (scale < 1.0)) ?
				(1.0 / sqrt(scale*(1.0-scale))) : 0.0;
			C_UInt8 *pGeno = Geno;
			for (long i=NumGeno; i > 0; i--)
			{
				*NormalizedGeno++ = (*pGeno < 3) ? (*pGeno-ave)*scale : 0.0;
				pGeno++;
			}
		} else
			for (; NumGeno > 0; NumGeno--) *NormalizedGeno++ = 0.0;
	}

	/// The thread entry for the calculation of SNP loadings
	void Entry_SNPLoadingCalc(PdThread Thread, int ThreadIndex, void *Param)
	{
		// The number of working samples
		const long n = MCWorkingGeno.Space.SampleNum();
		vector<C_UInt8> GenoBlock(n * BlockNumSNP);
		TdAlignPtr<double> NormalGeno(n);

		long _SNPstart, _SNPlen;
		while (RequireWork(&GenoBlock[0], _SNPstart, _SNPlen, false))
		{
			for (long iSNP=0; iSNP < _SNPlen; iSNP++)
			{
				// Normalize genotypes
				NormalizeGeno((&GenoBlock[0]) + iSNP*n, NormalGeno.get(), n);
				// SNP loadings, X'W
				double *Out = Out_Buffer + ((_SNPstart+iSNP)*OutputEigenDim);
				double *pEig = _EigenVectBuf->base();
				for (long iEig=0; iEig < OutputEigenDim; iEig++)
				{
					*Out++ = vt<double, av16Align>::DotProd(pEig, NormalGeno.get(), n);
					pEig += _EigenVectBuf->M();
				}
			}
			// Update progress
			{
				TdAutoMutex _m(_Mutex);
				MCWorkingGeno.Progress.Forward(_SNPlen);
			}
		}
	}

	/// Calculate the SNP loadings
	void GetPCAFreqScale(double OutFreq[], double OutScale[])
	{
		if (MCWorkingGeno.Space.SNPOrder())
		{
			// initialize
			const int nsnp = MCWorkingGeno.Space.SNPNum();
			vector<C_UInt8> buf(nsnp);
			vector<int> n(nsnp);
			for (int i=0; i < nsnp; i++)
			{
				n[i] = 0;
				OutFreq[i] = 0;
			}

			// for-loop for each sample
			for (int iSamp=0; iSamp < MCWorkingGeno.Space.SampleNum(); iSamp++)
			{
				MCWorkingGeno.Space.sampleRead(iSamp, 1, &buf[0], true);
				for (int i=0; i < nsnp; i++)
				{
					C_UInt8 &v = buf[i];
					if (v <= 2)
					{
						OutFreq[i] += v;
						n[i] ++;
					}
				}
			}
			// average
			for (int i=0; i < MCWorkingGeno.Space.SNPNum(); i++)
			{
				const int nn = n[i];
				const double sum = OutFreq[i];
				double ave = sum / nn;
				double scale = BayesianNormal ? ((sum+1.0)/(2*nn+2)) : (0.5*ave);
				scale = ((0.0 < scale) && (scale < 1.0)) ?
					(1.0 / sqrt(scale*(1.0-scale))) : 0.0;
				OutFreq[i] = ave;
				OutScale[i] = scale;
			}
		} else {
			// initialize
			vector<C_UInt8> buf(MCWorkingGeno.Space.SampleNum());

			// for-loop for each snp
			for (int isnp=0; isnp < MCWorkingGeno.Space.SNPNum(); isnp++)
			{
				int n = 0;
				double sum = 0;
				MCWorkingGeno.Space.snpRead(isnp, 1, &buf[0], false);
				for (int i=0; i < MCWorkingGeno.Space.SampleNum(); i++)
				{
					C_UInt8 &v = buf[i];
					if (v <= 2) { sum += v; n ++; }
				}
				double ave = sum / n;
				double scale = BayesianNormal ? ((sum+1.0)/(2*n+2)) : (0.5*ave);
				scale = ((0.0 < scale) && (scale < 1.0)) ?
					(1.0 / sqrt(scale*(1.0-scale))) : 0.0;
				OutFreq[isnp] = ave;
				OutScale[isnp] = scale;
			}
		}
	}

	/// Calculate the SNP loadings
	void DoSNPLoadingCalculate(double *EigenVal, int nEig, double *EigenVect,
		double TraceXTX, double *out_snploading, int NumThread, bool verbose,
		const char *Info)
	{
		// initialize mutex objects
		InitMutexObject();

		// Initialize progress information
		MCWorkingGeno.Progress.Info = Info;
		MCWorkingGeno.Progress.Show() = verbose;
		MCWorkingGeno.Progress.Init(MCWorkingGeno.Space.SNPNum());
		SNPStart = 0;
		OutputEigenDim = nEig;
		Out_Buffer = out_snploading;

		// Array of eigenvectors
		const long n = MCWorkingGeno.Space.SampleNum();
		double Scale = double(n-1) / TraceXTX;
		_EigenVectBuf = new CPCAMat<double>(OutputEigenDim, n);
		for (long i=0; i < OutputEigenDim; i++)
		{
			vt<double>::Mul(_EigenVectBuf->base() + i*_EigenVectBuf->M(),
				EigenVect, sqrt(Scale / EigenVal[i]), n);
			EigenVect += n;
		}

		// Threads
		GDS_Parallel_RunThreads(Entry_SNPLoadingCalc, NULL, NumThread);

		// destory the mutex objects
		DoneMutexObject();
		delete _EigenVectBuf; _EigenVectBuf = NULL;
	}


	// ================== Sample loadings ==================

	/// The thread entry for the calculation of genetic covariace matrix
	void Entry_SampLoadingCalc(PdThread Thread, int ThreadIndex, void *Param)
	{
		// The number of working samples
		const long n = MCWorkingGeno.Space.SNPNum();
		vector<C_UInt8> GenoBlock(n * BlockSamp);
		vector<double> buf(n);

		long _SampStart, _SampLen;
		while (RequireWorkSamp(&GenoBlock[0], _SampStart, _SampLen, true))
		{
			// for-loop each sample
			for (long iS=0; iS < _SampLen; iS++)
			{
				// fill the buffer
				C_UInt8 *g = (&GenoBlock[0]) + iS*n;
				double *p = &buf[0];
				for (long i=0; i < n; i++, p++, g++)
					*p = (*g <= 2) ? (*g - In_AveFreq[i]) : 0.0;

				// for-loop each eigenvector
				double *out = Out_Buffer + _SampStart + iS;
				for (int i=0; i < OutputEigenDim; i++)
				{
					double *p1 = In_EigenVect + i;
					double *p2 = &buf[0];
					double sum = 0;
					for (int j=0; j < n; j++)
					{
						sum += (*p1) * (*p2);
						p1 += OutputEigenDim; p2 ++;
					}
					*out = sum;
					out += MCWorkingGeno.Space.SampleNum();
				}
			}
			// Update progress
			{
				TdAutoMutex _m(_Mutex);
				MCWorkingGeno.Progress.Forward(_SampLen);
            }
		}
	}

    /// Calculate the genetic covariace
	void DoSampLoadingCalculate(double *Ave_Freq, double *Scale, int EigenCnt,
		double *SNP_Loadings, double *EigenVal, int Num, double TraceXTX,
		double *out_SampLoadings, int NumThread, const char *Info,
		bool verbose)
	{
		// initialize mutex objects
		InitMutexObject();

		// Initialize progress information
		const int n = MCWorkingGeno.Space.SNPNum();
		MCWorkingGeno.Progress.Info = Info;
		MCWorkingGeno.Progress.Show() = verbose;
		MCWorkingGeno.Progress.Init(MCWorkingGeno.Space.SampleNum());

		double ss = double(Num-1) / TraceXTX;
		vector<double> eigen(EigenCnt);
		for (long i=0; i < EigenCnt; i++)
			eigen[i] = sqrt(ss / EigenVal[i]);

		OutputEigenDim = EigenCnt;
		In_AveFreq = Ave_Freq; In_EigenVect = SNP_Loadings;
		Out_Buffer = out_SampLoadings;
		double *p = SNP_Loadings;
		for (int i=0; i < n; i++)
		{
			for (int j=0; j < EigenCnt; j++)
				*p++ *= Scale[i] * eigen[j];
		}

		// threads
		SampStart = 0;
		GDS_Parallel_RunThreads(Entry_SampLoadingCalc, NULL, NumThread);

		// destory the mutex objects
		DoneMutexObject();
	}



	// ---------------------------------------------------------------------
	// Admixture Analyses
	// ---------------------------------------------------------------------

	// ================== Relative IBD matrix ==================

	/// Missing genotypic flags:
	//      0 -- missing value, other value -- valid genotype
	static vector<C_UInt8> Admix_Missing_Flag;

	/// Adjusted genotypic values for diagonal
	static vector<double> Admix_Adj_Geno;
	/// the product of allele frequencies: 4 * p_l * (1 - p_l)
	static vector<double> Admix_P_ONE_P;

	/// Convert the raw genotypes to the mean-adjusted genotypes
	static void _Do_Admix_RatioOfAvg_ReadBlock(C_UInt8 *GenoBuf,
		long Start, long SNP_Cnt, void* Param)
	{
		// init ...
		const int n = MCWorkingGeno.Space.SampleNum();
		memset(&PCA_gSum[0], 0, SNP_Cnt*sizeof(int));
		memset(&PCA_gNum[0], 0, SNP_Cnt*sizeof(int));

		C_UInt8 *p, *pMissing;
		double *pf, *pp;
		double *pAdjGeno = &Admix_Adj_Geno[0];

		// calculate the averages of genotypes for each SelSNP
		// store the values in PCA_gSum
		p = GenoBuf;
		pMissing = &Admix_Missing_Flag[0];
		pp = PCA_Mat.base();
		for (long iSample=0; iSample < n; iSample++)
		{
			pf = pp; pp += PCA_Mat.M();
			for (long iSNP=0; iSNP < SNP_Cnt; iSNP++)
			{
				if (*p < 3)
				{
					PCA_gSum[iSNP] += *p;
					PCA_gNum[iSNP]++;
					*pAdjGeno += (*p) * (2 - *p);
					*pMissing++ = true;
				} else {
					*pMissing++ = false;
				}
				*pf++ = *p++;
			}
			pAdjGeno ++;
		}

		// PCA_gSum = 2 * \bar{p}_l
		for (long iSNP=0; iSNP < SNP_Cnt; iSNP++)
		{
			double &fv = tmpBuf.get()[iSNP];
			long gN = PCA_gNum[iSNP];
			if (gN > 0)
				fv = (double)PCA_gSum[iSNP] / gN;
			else
				fv = 0;
		}

		// Averaging and set missing values to ZERO
		pp = PCA_Mat.base(); // G@ij - 2\bar{p}_l
		for (long iSample=0; iSample < n; iSample++)
		{
			vt<double, av16Align>::Sub(pp, pp, tmpBuf.get(), SNP_Cnt);
			pp += PCA_Mat.M();
		}

		// missing values
		p = GenoBuf; pp = PCA_Mat.base();
		for (long iSample=0; iSample < n; iSample++)
		{
			pf = pp; pp += PCA_Mat.M();
			for (long iSNP=0; iSNP < SNP_Cnt; iSNP++)
			{
				if (*p > 2) *pf = 0;
				p++; pf++;
			}
		}

		// 4 * p_l * (1 - p_l) -> tmpBuf
		for (long iSNP=0; iSNP < SNP_Cnt; iSNP++)
		{
			double &f = tmpBuf.get()[iSNP];
			f *= 0.5;
			f = 4 * f * (1 - f);
		}
	}

	/// Compute the relative IBD matrix
	static void _Do_Admix_RatioOfAvg_Compute(int ThreadIndex,
		long Start, long SNP_Cnt, void* Param)
	{
		double **Ptr = (double**)Param;

		// numerator
		IdMatTri I = PCA_Thread_MatIdx[ThreadIndex];
		PCA_Mat.MulAdd(I, PCA_Thread_MatCnt[ThreadIndex], SNP_Cnt,
			Ptr[0] + I.Offset());

		// denominator
		I = PCA_Thread_MatIdx[ThreadIndex];
		double *pAFreq = Ptr[1] + I.Offset();
		for (C_Int64 L = PCA_Thread_MatCnt[ThreadIndex]; L > 0; L--)
		{
			C_UInt8 *p1 = &(Admix_Missing_Flag[0]) + SNP_Cnt*I.Row();
			C_UInt8 *p2 = &(Admix_Missing_Flag[0]) + SNP_Cnt*I.Column();
			for (long i=0; i< SNP_Cnt; i++)
			{
				if (p1[i] && p2[i])
					*pAFreq += tmpBuf.get()[i];
			}
			pAFreq ++;
			++ I;
		}
	}

	/// Calculate the genetic covariace by ratio of averages (or ratio of sums)
	void DoAdmixCalc_RatioOfAvg(CdMatTri<double> &OutIBD, bool DiagAdj,
		int NumThread, bool verbose)
	{
		// initialize ...
		const long n = OutIBD.N();  // the number of individuals

		PCA_gSum.resize(BlockNumSNP);
		PCA_gNum.resize(BlockNumSNP);
		PCA_Mat.Reset(n, BlockNumSNP);
		tmpBuf.Reset(PCA_Mat.M());
		Admix_Missing_Flag.resize(BlockNumSNP * n);
		Admix_Adj_Geno.resize(n);

		memset(OutIBD.get(), 0, sizeof(double)*OutIBD.Size());
		memset(&Admix_Adj_Geno[0], 0, sizeof(double)*n);

		MCWorkingGeno.Progress.Info = "Eigen-analysis:";
		MCWorkingGeno.Progress.Show() = verbose;
		MCWorkingGeno.InitParam(true, true, BlockNumSNP);

		CdMatTri<double> AFreqProd(n);
		memset(AFreqProd.get(), 0, sizeof(double)*AFreqProd.Size());

		double *Ptr[2] = { OutIBD.get(), AFreqProd.get() };

		MCWorkingGeno.SplitJobs(NumThread, n, PCA_Thread_MatIdx,
			PCA_Thread_MatCnt);
		MCWorkingGeno.Run(NumThread, &_Do_Admix_RatioOfAvg_ReadBlock,
			&_Do_Admix_RatioOfAvg_Compute, (void*)Ptr);

		// output
		double *p = OutIBD.get();
		double *s = AFreqProd.get();
		for (long i=0; i < n; i++)
		{
			if (DiagAdj)
				(*p) -= Admix_Adj_Geno[i];
			for (long j=i; j < n; j++)
				{ *p /= *s; p++; s++; }
		}
	}


	// ==================================== //

	/// Convert the raw genotypes to the mean-adjusted genotypes
	static void _Do_Admix_AvgOfRatio_ReadBlock(C_UInt8 *GenoBuf,
		long Start, long SNP_Cnt, void* Param)
	{
		// init ...
		const int n = MCWorkingGeno.Space.SampleNum();
		memset(&PCA_gSum[0], 0, SNP_Cnt*sizeof(int));
		memset(&PCA_gNum[0], 0, SNP_Cnt*sizeof(int));

		C_UInt8 *p, *pMissing;
		double *pf, *pp;

		// calculate the averages of genotypes for each SelSNP
		// store the values in PCA_gSum
		p = GenoBuf;
		pMissing = &Admix_Missing_Flag[0];
		pp = PCA_Mat.base();
		for (long iSample=0; iSample < n; iSample++)
		{
			pf = pp; pp += PCA_Mat.M();
			for (long iSNP=0; iSNP < SNP_Cnt; iSNP++)
			{
				if (*p < 3)
				{
					PCA_gSum[iSNP] += *p;
					PCA_gNum[iSNP]++;
					*pMissing++ = true;
				} else {
					*pMissing++ = false;
				}
				*pf++ = *p++;
			}
		}

		// PCA_gSum = 2 * \bar{p}_l
		for (long iSNP=0; iSNP < SNP_Cnt; iSNP++)
		{
			double &fv = tmpBuf.get()[iSNP];
			long gN = PCA_gNum[iSNP];
			if (gN > 0)
				fv = (double)PCA_gSum[iSNP] / gN;
			else
				fv = 0;
		}

		// Averaging and set missing values to ZERO
		pp = PCA_Mat.base(); // G@ij - 2\bar{p}_l
		for (long iSample=0; iSample < n; iSample++)
		{
			vt<double, av16Align>::Sub(pp, pp, tmpBuf.get(), SNP_Cnt);
			pp += PCA_Mat.M();
		}

		// 1 / (2*sqrt(p@j*(1-p@j)))
		pf = tmpBuf.get();
		for (long iSNP=0; iSNP < SNP_Cnt; iSNP++, pf++)
		{
			double scale = (*pf) * 0.5;
			if ((0.0 < scale) && (scale < 1.0))
				*pf = 0.5 / sqrt(scale * (1.0-scale));
			else
				*pf = 0;
		}

		// G * (2 - G) / (4*p*(1-p))
		double *pAdjGeno = &Admix_Adj_Geno[0];
		p = GenoBuf;
		for (long iSample=0; iSample < n; iSample++)
		{
			pf = tmpBuf.get();
			for (long iSNP=0; iSNP < SNP_Cnt; iSNP++)
			{
				if (*p < 3)
					*pAdjGeno += (*p) * (2 - *p) * (*pf) * (*pf);
				p++; pf++;
			}
			pAdjGeno ++;
		}

		// (G@ij - 2p@j) / (2 * sqrt(p@j*(1-p@j)))
		pp = PCA_Mat.base();
		for (long iSample=0; iSample < n; iSample++)
		{
			vt<double, av16Align>::Mul(pp, pp, tmpBuf.get(), SNP_Cnt);
			pp += PCA_Mat.M();
		}

		// missing values
		p = GenoBuf; pp = PCA_Mat.base();
		for (long iSample=0; iSample < n; iSample++)
		{
			pf = pp; pp += PCA_Mat.M();
			for (long iSNP=0; iSNP < SNP_Cnt; iSNP++)
			{
				if (*p > 2) *pf = 0;
				p++; pf++;
			}
		}
	}

	/// Compute the relative IBD matrix
	static void _Do_Admix_AvgOfRatio_Compute(int ThreadIndex,
		long Start, long SNP_Cnt, void* Param)
	{
		void **Ptr = (void**)Param;
		double *Ptr1 = (double*)(Ptr[0]);

		// numerator
		IdMatTri I = PCA_Thread_MatIdx[ThreadIndex];
		PCA_Mat.MulAdd(I, PCA_Thread_MatCnt[ThreadIndex], SNP_Cnt,
			Ptr1 + I.Offset());

		// denominator
		I = PCA_Thread_MatIdx[ThreadIndex];
		int *Ptr2 = (int*)(Ptr[1]) + I.Offset();
		for (C_Int64 L = PCA_Thread_MatCnt[ThreadIndex]; L > 0; L--)
		{
			C_UInt8 *p1 = &Admix_Missing_Flag[0] + SNP_Cnt*I.Row();
			C_UInt8 *p2 = &Admix_Missing_Flag[0] + SNP_Cnt*I.Column();
			for (long i=0; i< SNP_Cnt; i++)
			{
				if ((*p1) && (*p2)) (*Ptr2) ++;
				p1 ++; p2 ++;
			}
			Ptr2 ++;
			++ I;
		}
	}

	/// Calculate the genetic covariace by averaging ratios
	void DoAdmixCalc_AvgOfRatios(CdMatTri<double> &OutIBD, bool DiagAdj,
		int NumThread, bool verbose)
	{
		// initialize ...
		const long n = OutIBD.N();  // the number of individuals

		PCA_gSum.resize(BlockNumSNP);
		PCA_gNum.resize(BlockNumSNP);
		PCA_Mat.Reset(n, BlockNumSNP);
		tmpBuf.Reset(PCA_Mat.M());
		Admix_Missing_Flag.resize(BlockNumSNP * n);
		Admix_Adj_Geno.resize(n);

		memset(OutIBD.get(), 0, sizeof(double)*OutIBD.Size());
		memset(&Admix_Adj_Geno[0], 0, sizeof(double)*n);

		MCWorkingGeno.Progress.Info = "Eigen-analysis:";
		MCWorkingGeno.Progress.Show() = verbose;
		MCWorkingGeno.InitParam(true, true, BlockNumSNP);

		CdMatTri<int> NumValid(n);
		memset(NumValid.get(), 0, sizeof(int)*NumValid.Size());

		void *Ptr[2] = { OutIBD.get(), NumValid.get() };

		MCWorkingGeno.SplitJobs(NumThread, n, PCA_Thread_MatIdx,
			PCA_Thread_MatCnt);
		MCWorkingGeno.Run(NumThread, &_Do_Admix_AvgOfRatio_ReadBlock,
			&_Do_Admix_AvgOfRatio_Compute, (void*)Ptr);

		// output
		double *p = OutIBD.get();
		int *s = NumValid.get();
		for (long i=0; i < n; i++)
		{
			if (DiagAdj)
				(*p) -= Admix_Adj_Geno[i];
			for (long j=i; j < n; j++)
				{ *p /= *s; p++; s++; }
		}
	}


	// ==================================== //

	/// Convert the raw genotypes to the mean-adjusted genotypes
	static void _Do_GRM_AvgOfRatio_ReadBlock(C_UInt8 *GenoBuf,
		long Start, long SNP_Cnt, void* Param)
	{
		const static double SCALE = 1.0 / sqrt(2.0);

		// init ...
		const int n = MCWorkingGeno.Space.SampleNum();
		memset(&PCA_gSum[0], 0, SNP_Cnt*sizeof(int));
		memset(&PCA_gNum[0], 0, SNP_Cnt*sizeof(int));

		C_UInt8 *p, *pMissing;
		double *pf, *pp;

		// calculate the averages of genotypes for each SelSNP
		// store the values in PCA_gSum
		p = GenoBuf;
		pMissing = &Admix_Missing_Flag[0];
		pp = PCA_Mat.base();
		for (long iSample=0; iSample < n; iSample++)
		{
			pf = pp; pp += PCA_Mat.M();
			for (long iSNP=0; iSNP < SNP_Cnt; iSNP++)
			{
				if (*p < 3)
				{
					PCA_gSum[iSNP] += *p;
					PCA_gNum[iSNP]++;
					*pMissing++ = true;
				} else {
					*pMissing++ = false;
				}
				*pf++ = *p++;
			}
		}

		// PCA_gSum = 2 * \bar{p}_l
		for (long iSNP=0; iSNP < SNP_Cnt; iSNP++)
		{
			double &fv = tmpBuf.get()[iSNP];
			long gN = PCA_gNum[iSNP];
			if (gN > 0)
				fv = (double)PCA_gSum[iSNP] / gN;
			else
				fv = 0;
		}

		// Averaging and set missing values to ZERO
		pp = PCA_Mat.base(); // G@ij - 2\bar{p}_l
		for (long iSample=0; iSample < n; iSample++)
		{
			vt<double, av16Align>::Sub(pp, pp, tmpBuf.get(), SNP_Cnt);
			pp += PCA_Mat.M();
		}

		// (1 - 2p) * (G - 2p) / (2*p*(1-p))
		double *pAdjGeno = &Admix_Adj_Geno[0];
		p = GenoBuf;
		for (long iSample=0; iSample < n; iSample++)
		{
			pf = tmpBuf.get();
			for (long iSNP=0; iSNP < SNP_Cnt; iSNP++)
			{
				double scale = (*pf) * 0.5;
				if ((0.0 < scale) && (scale < 1.0))
					scale = 1 / (2 * scale * (1 - scale));
				else
					scale = 0;
				if (*p < 3)
					*pAdjGeno += (1 - *pf) * (*p - *pf) * scale;
				p++; pf++;
			}
			pAdjGeno ++;
		}

		// 1 / (sqrt(2)*sqrt(p@j*(1-p@j)))
		pf = tmpBuf.get();
		for (long iSNP=0; iSNP < SNP_Cnt; iSNP++, pf++)
		{
			double scale = (*pf) * 0.5;
			if ((0.0 < scale) && (scale < 1.0))
				*pf = SCALE / sqrt(scale * (1 - scale));
			else
				*pf = 0;
		}

		// (G@ij - 2p@j) / (sqrt(2) * sqrt(p@j*(1-p@j)))
		pp = PCA_Mat.base();
		for (long iSample=0; iSample < n; iSample++)
		{
			vt<double, av16Align>::Mul(pp, pp, tmpBuf.get(), SNP_Cnt);
			pp += PCA_Mat.M();
		}

		// missing values
		p = GenoBuf; pp = PCA_Mat.base();
		for (long iSample=0; iSample < n; iSample++)
		{
			pf = pp; pp += PCA_Mat.M();
			for (long iSNP=0; iSNP < SNP_Cnt; iSNP++)
			{
				if (*p > 2) *pf = 0;
				p++; pf++;
			}
		}
	}

	/// Calculate the genetic relationship matrix (GRM)
	void DoGRMCalc(CdMatTri<double> &OutIBD, int NumThread,
		bool verbose)
	{
		// initialize ...
		const long n = OutIBD.N();  // the number of individuals

		PCA_gSum.resize(BlockNumSNP);
		PCA_gNum.resize(BlockNumSNP);
		PCA_Mat.Reset(n, BlockNumSNP);
		tmpBuf.Reset(PCA_Mat.M());
		Admix_Missing_Flag.resize(BlockNumSNP * n);
		Admix_Adj_Geno.resize(n);

		memset(OutIBD.get(), 0, sizeof(double)*OutIBD.Size());
		memset(&Admix_Adj_Geno[0], 0, sizeof(double)*n);

		MCWorkingGeno.Progress.Info = "GRM-analysis:";
		MCWorkingGeno.Progress.Show() = verbose;
		MCWorkingGeno.InitParam(true, true, BlockNumSNP);

		CdMatTri<int> NumValid(n);
		memset(NumValid.get(), 0, sizeof(int)*NumValid.Size());

		void *Ptr[2] = { OutIBD.get(), NumValid.get() };

		MCWorkingGeno.SplitJobs(NumThread, n, PCA_Thread_MatIdx,
			PCA_Thread_MatCnt);
		MCWorkingGeno.Run(NumThread, &_Do_GRM_AvgOfRatio_ReadBlock,
			&_Do_Admix_AvgOfRatio_Compute, (void*)Ptr);

		// output
		double *p = OutIBD.get();
		int *s = NumValid.get();
		for (long i=0; i < n; i++)
		{
			(*p) -= Admix_Adj_Geno[i];
			for (long j=i; j < n; j++)
				{ *p /= *s; p++; s++; }
		}
	}
}


using namespace PCA;

extern "C"
{

/// get the eigenvalues and eigenvectors, return 'nProtect'
static int GetEigen(double *pMat, int n, int nEig, const char *EigMethod,
	SEXP &EigVal, SEXP &EigVect)
{
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
			throw ErrCoreArray("LAPACK::DSPEV error (INFO: %d)!", info);

		// output eigenvalues
		vt<double>::Sub(REAL(EigVal), 0.0, REAL(EigVal), n);

		// output eigenvectors
		EigVect = PROTECT(allocMatrix(REALSXP, n, nEig));
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
		vector<double> tmp_Eigen(n);

		EigVal = PROTECT(NEW_NUMERIC(nEig));
		EigVect = PROTECT(allocMatrix(REALSXP, n, nEig));
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

		F77_NAME(dspevx)("V", "I", "L", &n, pMat,
			&VL, &VU, &IL, &IU, &ABSTOL,
			&M, &tmp_Eigen[0], REAL(EigVect), &LDZ,
			&tmp_Work[0], &tmp_IWork[0], &ifail[0], &info);
		if (info != 0)
			throw ErrCoreArray("LAPACK::DSPEVX error (INFO: %d)!", info);

		// output eigenvalues
		for (int i=0; i < nEig; i++)
			REAL(EigVal)[i] = -tmp_Eigen[i];
	} else
		throw ErrCoreArray("Unknown 'eigen.method'.");

	return nProtected;
}



// ========================================================================
// the functions for Principal Component Analysis (PCA)
// ========================================================================

/// to compute the eigenvalues and eigenvectors
COREARRAY_DLL_EXPORT SEXP gnrPCA(SEXP EigenCnt, SEXP NumThread,
	SEXP _BayesianNormal, SEXP NeedGenMat, SEXP GenMat_Only,
	SEXP EigenMethod, SEXP _Verbose)
{
	bool verbose = SEXP_Verbose(_Verbose);

	COREARRAY_TRY

		// ======== To cache the genotype data ========
		CachingSNPData("PCA", verbose);

		// ======== The calculation of genetic covariance matrix ========

		// the number of samples
		const R_xlen_t n = MCWorkingGeno.Space.SampleNum();
		// set parameters
		PCA::BayesianNormal = (LOGICAL(_BayesianNormal)[0] == TRUE);
		PCA::AutoDetectSNPBlockSize(n);
		// the upper-triangle genetic covariance matrix
		CdMatTri<double> Cov(n);

		// Calculate the genetic covariace
		PCA::DoCovCalculate(Cov, Rf_asInteger(NumThread), "PCA:", verbose);

		// Normalize
		double TraceXTX = Cov.Trace();
		double scale = double(n-1) / TraceXTX;
		vt<double, av16Align>::Mul(Cov.get(), Cov.get(), scale, Cov.Size());
		double TraceVal = Cov.Trace();

		// ======== output ========

		int nProtected = 0;
		PROTECT(rv_ans = NEW_LIST(5));
		nProtected ++;

		SET_ELEMENT(rv_ans, 0, ScalarReal(TraceXTX));
		SET_ELEMENT(rv_ans, 4, ScalarReal(TraceVal));

		if (LOGICAL(NeedGenMat)[0] == TRUE)
		{
			SEXP tmp;
			PROTECT(tmp = allocMatrix(REALSXP, n, n));
			nProtected ++;
			SET_ELEMENT(rv_ans, 1, tmp);
			Cov.SaveTo(REAL(tmp));
		}

		// ======== eigenvectors and eigenvalues ========

		if (LOGICAL(GenMat_Only)[0] != TRUE)
		{
			if (verbose)
			{
				Rprintf("PCA:\t%s\tBegin (eigenvalues and eigenvectors)\n",
					NowDateToStr().c_str());
			}

			vt<double>::Sub(Cov.get(), 0.0, Cov.get(), Cov.Size());

			int nEig = Rf_asInteger(EigenCnt);
			if (nEig <= 0)
				throw ErrCoreArray("Invalid 'eigen.cnt'.");
			if (nEig > n) nEig = n;

			SEXP EigVal  = R_NilValue;
			SEXP EigVect = R_NilValue;
			nProtected += GetEigen(Cov.get(), n, nEig,
				CHAR(STRING_ELT(EigenMethod, 0)), EigVal, EigVect);
			SET_ELEMENT(rv_ans, 2, EigVal);
			SET_ELEMENT(rv_ans, 3, EigVect);

			if (verbose)
			{
				Rprintf("PCA:\t%s\tEnd (eigenvalues and eigenvectors)\n",
					NowDateToStr().c_str());
			}
		}

		UNPROTECT(nProtected);

	COREARRAY_CATCH
}


/// to calculate the SNP correlations
COREARRAY_DLL_EXPORT SEXP gnrPCACorr(SEXP LenEig, SEXP EigenVect,
	SEXP NumThread, SEXP _Verbose)
{
	bool verbose = SEXP_Verbose(_Verbose);

	COREARRAY_TRY

		// ======== To cache the genotype data ========
		CachingSNPData("SNP Correlation", verbose);

		// ======== To compute the snp correlation ========
		PCA::AutoDetectSNPBlockSize(MCWorkingGeno.Space.SampleNum());

		PROTECT(rv_ans = allocMatrix(REALSXP, INTEGER(LenEig)[0],
			MCWorkingGeno.Space.SNPNum()));

		PCA::DoSNPCoeffCalculate(INTEGER(LenEig)[0], REAL(EigenVect),
			REAL(rv_ans), INTEGER(NumThread)[0], verbose, "SNP Correlation:");

		UNPROTECT(1);

	COREARRAY_CATCH
}


/// to calculate the SNP loadings
COREARRAY_DLL_EXPORT SEXP gnrPCASNPLoading(SEXP EigenVal, SEXP DimCnt,
	SEXP EigenVect, SEXP TraceXTX, SEXP NumThread, SEXP Bayesian,
	SEXP _Verbose)
{
	bool verbose = SEXP_Verbose(_Verbose);

	COREARRAY_TRY

		// ======== To cache the genotype data ========
		CachingSNPData("SNP Loading", verbose);

		// ======== To compute the snp correlation ========
		PCA::AutoDetectSNPBlockSize(MCWorkingGeno.Space.SampleNum());
		PCA::BayesianNormal = (LOGICAL(Bayesian)[0] == TRUE);

		PROTECT(rv_ans = NEW_LIST(3));

		SEXP loading, afreq, scale;
		PROTECT(loading = allocMatrix(REALSXP, INTEGER(DimCnt)[1],
			MCWorkingGeno.Space.SNPNum()));
		SET_ELEMENT(rv_ans, 0, loading);

		PROTECT(afreq = NEW_NUMERIC(MCWorkingGeno.Space.SNPNum()));
		SET_ELEMENT(rv_ans, 1, afreq);

		PROTECT(scale = NEW_NUMERIC(MCWorkingGeno.Space.SNPNum()));
		SET_ELEMENT(rv_ans, 2, scale);

		PCA::GetPCAFreqScale(REAL(afreq), REAL(scale));
		PCA::DoSNPLoadingCalculate(REAL(EigenVal), INTEGER(DimCnt)[1],
			REAL(EigenVect), REAL(TraceXTX)[0], REAL(loading),
			INTEGER(NumThread)[0], verbose,
			"SNP Loading:");

		UNPROTECT(4);

	COREARRAY_CATCH
}


/// to calculate the sample loadings from SNP loadings
COREARRAY_DLL_EXPORT SEXP gnrPCASampLoading(SEXP Num, SEXP EigenVal,
	SEXP EigenCnt, SEXP SNPLoadings, SEXP TraceXTX, SEXP AveFreq, SEXP Scale,
	SEXP NumThread, SEXP _Verbose)
{
	bool verbose = SEXP_Verbose(_Verbose);

	COREARRAY_TRY

		// ======== To cache the genotype data ========
		CachingSNPData("Sample Loading", verbose);

		PROTECT(rv_ans = allocMatrix(REALSXP,
			MCWorkingGeno.Space.SampleNum(), INTEGER(EigenCnt)[0]));

		// ======== To compute the snp correlation ========
		PCA::DoSampLoadingCalculate(REAL(AveFreq), REAL(Scale),
			INTEGER(EigenCnt)[0], REAL(SNPLoadings),
			REAL(EigenVal), INTEGER(Num)[0], REAL(TraceXTX)[0],
			REAL(rv_ans), INTEGER(NumThread)[0], "Sample Loading:",
			verbose);

		UNPROTECT(1);

	COREARRAY_CATCH
}


// ======================================================================*
// Genetic relationship matrix
// ======================================================================*

/// to compute the eigenvalues and eigenvectors
COREARRAY_DLL_EXPORT SEXP gnrGRM(SEXP _NumThread, SEXP _Verbose)
{
	bool verbose = SEXP_Verbose(_Verbose);

	COREARRAY_TRY

		// ======== To cache the genotype data ========
		CachingSNPData("GRM calculation", verbose);

		// ======== The calculation of genetic covariance matrix ========

		// the number of samples
		const R_xlen_t n = MCWorkingGeno.Space.SampleNum();

		// set internal parameters
		PCA::AutoDetectSNPBlockSize(n);

		// the upper-triangle IBD matrix
		CdMatTri<double> IBD(n);

		// Calculate the genetic covariace
		PCA::DoGRMCalc(IBD, INTEGER(_NumThread)[0], verbose);

		// Output
		PROTECT(rv_ans = allocMatrix(REALSXP, n, n));
		double *base = REAL(rv_ans);
		double *p = IBD.get();
		for (R_xlen_t i=0; i < n; i++)
		{
			for (R_xlen_t j=i; j < n; j++)
			{
				base[i*n + j] = base[j*n + i] = *p;
				p ++;
			}
		}
		UNPROTECT(1);

	COREARRAY_CATCH
}


// ======================================================================*
// the functions for eigen-analysis for admixtures
//

struct TEigPair
{
	double EigVal;
	int Index;
	TEigPair(double e, int i) { EigVal = e; Index = i; }
};

static bool _EigComp(const TEigPair &i, const TEigPair &j)
{
	return (fabs(i.EigVal) >= fabs(j.EigVal));
}

/// to compute the eigenvalues and eigenvectors
COREARRAY_DLL_EXPORT SEXP gnrEIGMIX(SEXP _EigenCnt, SEXP _NumThread,
	SEXP _NeedIBDMat, SEXP _IBDMatOnly, SEXP _Method, SEXP _DiagAdj,
	SEXP _Verbose)
{
	bool verbose = SEXP_Verbose(_Verbose);

	COREARRAY_TRY

		// ======== To cache the genotype data ========
		CachingSNPData("Eigen-analysis", verbose);

		// ======== The calculation of genetic covariance matrix ========

		// the number of samples
		const R_xlen_t n = MCWorkingGeno.Space.SampleNum();

		// set internal parameters
		PCA::AutoDetectSNPBlockSize(n);

		// the upper-triangle IBD matrix
		CdMatTri<double> IBD(n);

		// Calculate the genetic covariace
		switch (INTEGER(_Method)[0])
		{
			case 1:
				PCA::DoAdmixCalc_RatioOfAvg(IBD, LOGICAL(_DiagAdj)[0]==TRUE,
					INTEGER(_NumThread)[0], verbose);
				break;
			case 2:
				PCA::DoAdmixCalc_AvgOfRatios(IBD, LOGICAL(_DiagAdj)[0]==TRUE,
					INTEGER(_NumThread)[0], verbose);
				break;
			default:
				throw "Invalid method.";
		}

		// ======== The calculation of eigenvectors and eigenvalues ========

		int nProtected = 0;
		SEXP EigenVal=NULL, EigenVec=NULL, IBDMat=NULL;

		if (LOGICAL(_NeedIBDMat)[0] == TRUE)
		{
			PROTECT(IBDMat = allocMatrix(REALSXP, n, n));
			nProtected ++;

			double *base = REAL(IBDMat);
			double *p = IBD.get();
			for (R_xlen_t i=0; i < n; i++)
			{
				for (R_xlen_t j=i; j < n; j++)
				{
					base[i*n + j] = base[j*n + i] = *p;
					p ++;
				}
			}
		}

		if (LOGICAL(_IBDMatOnly)[0] != TRUE)
		{
			const size_t NN = n;
			vector<double> tmp_Work(NN*3);
			vector<double> tmp_EigenVec(NN*NN);

			vt<double>::Sub(IBD.get(), 0.0, IBD.get(), IBD.Size());
			if (verbose)
			{
				Rprintf("Eigen-analysis:\t%s\tBegin (eigenvalues and eigenvectors)\n",
					NowDateToStr().c_str());
			}

			int eigencnt = INTEGER(_EigenCnt)[0];
			if (eigencnt > n) eigencnt = n;

			PROTECT(EigenVal = NEW_NUMERIC(n));
			PROTECT(EigenVec = allocMatrix(REALSXP, n, eigencnt));
			nProtected += 2;

			{
				int info = 0;
				int _n = n;
				F77_NAME(dspev)("V", "L", &_n, IBD.get(), REAL(EigenVal),
					&tmp_EigenVec[0], &_n, &tmp_Work[0], &info);
				if (info != 0)
					throw "LAPACK::SPEV error!";
			}

			// sort in a decreasing order
			vector<TEigPair> lst;
			for (int i=0; i < n; i++)
				lst.push_back(TEigPair(REAL(EigenVal)[i], i));
			sort(lst.begin(), lst.end(), _EigComp);

			// output eigenvalues
			for (int i=0; i < n; i++)
				REAL(EigenVal)[i] = - lst[i].EigVal;

			// output eigenvectors
			for (int i=0; i < eigencnt; i++)
			{
				double *p = &tmp_EigenVec[0] + lst[i].Index * n;
				memmove(&(REAL(EigenVec)[i*n]), p, sizeof(double)*n);
			}

			if (verbose)
			{
				Rprintf("Eigen-analysis:\t%s\tEnd (eigenvalues and eigenvectors)\n",
					NowDateToStr().c_str());
			}
		}

		PROTECT(rv_ans = NEW_LIST(3));
		nProtected ++;

		if (EigenVal != NULL)
			SET_ELEMENT(rv_ans, 0, EigenVal);
		if (EigenVec != NULL)
			SET_ELEMENT(rv_ans, 1, EigenVec);
		if (IBDMat != NULL)
			SET_ELEMENT(rv_ans, 2, IBDMat);

		SEXP tmp;
		PROTECT(tmp = NEW_CHARACTER(3));
		nProtected ++;
			SET_STRING_ELT(tmp, 0, mkChar("eigenval"));
			SET_STRING_ELT(tmp, 1, mkChar("eigenvect"));
			SET_STRING_ELT(tmp, 2, mkChar("ibdmat"));
		SET_NAMES(rv_ans, tmp);

		UNPROTECT(nProtected);

	COREARRAY_CATCH
}

}

#endif  /* _HEADER_PCA_ */
