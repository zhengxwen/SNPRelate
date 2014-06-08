// ===========================================================
//     _/_/_/   _/_/_/  _/_/_/_/    _/_/_/_/  _/_/_/   _/_/_/
//      _/    _/       _/             _/    _/    _/   _/   _/
//     _/    _/       _/_/_/_/       _/    _/    _/   _/_/_/
//    _/    _/       _/             _/    _/    _/   _/
// _/_/_/   _/_/_/  _/_/_/_/_/     _/     _/_/_/   _/_/
// ===========================================================
//
// genPCA.cpp: Principal component analysis on genome-wide association studies
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


#ifndef _FuncPCA_H_
#define _FuncPCA_H_

// CoreArray library header
#include <dType.h>
#include <dVect.h>
#include <CoreGDSLink.h>
#include <dGenGWAS.h>

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
			fBuf.reset(new tfloat[n*m]);
			fN = n; fM = m;
		}

		void COREARRAY_CALL_ALIGN MulAdd(IdMatTri &Idx, size_t IdxCnt, size_t ArrayLength, tfloat *OutBuf)
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
		inline tfloat* base() { return fBuf.get(); }
	private:
		auto_ptr<tfloat> fBuf;
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

		void COREARRAY_CALL_ALIGN MulAdd(IdMatTri &Idx, size_t IdxCnt, size_t ArrayLength, float *OutBuf)
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
					rv128 = _mm_add_ps(rv128, _mm_mul_ps(_mm_load_ps(&p1[0]), _mm_load_ps(&p2[0])));
					rv128 = _mm_add_ps(rv128, _mm_mul_ps(_mm_load_ps(&p1[4]), _mm_load_ps(&p2[4])));
					rv128 = _mm_add_ps(rv128, _mm_mul_ps(_mm_load_ps(&p1[8]), _mm_load_ps(&p2[8])));
					rv128 = _mm_add_ps(rv128, _mm_mul_ps(_mm_load_ps(&p1[12]), _mm_load_ps(&p2[12])));
					p1 += 16; p2 += 16; k -= 16;
				}
				while (k >= 4)
				{
					rv128 = _mm_add_ps(rv128, _mm_mul_ps(_mm_load_ps(p1), _mm_load_ps(p2)));
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

		void COREARRAY_CALL_ALIGN MulAdd(IdMatTri &Idx, size_t IdxCnt, size_t ArrayLength, double *OutBuf)
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
					rv128 = _mm_add_pd(rv128, _mm_mul_pd(_mm_load_pd(&p1[0]), _mm_load_pd(&p2[0])));
					rv128 = _mm_add_pd(rv128, _mm_mul_pd(_mm_load_pd(&p1[2]), _mm_load_pd(&p2[2])));
					rv128 = _mm_add_pd(rv128, _mm_mul_pd(_mm_load_pd(&p1[4]), _mm_load_pd(&p2[4])));
					rv128 = _mm_add_pd(rv128, _mm_mul_pd(_mm_load_pd(&p1[6]), _mm_load_pd(&p2[6])));
					p1 += 8; p2 += 8; k -= 8;
				}

				while (k >= 2)
				{
					rv128 = _mm_add_pd(rv128, _mm_mul_pd(_mm_load_pd(p1), _mm_load_pd(p2)));
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

	/// detect the effective value for BlockSNP
	void AutoDetectSNPBlockSize(int nSamp, bool Detect=true)
	{
		if (Detect)
		{
			long L2Cache = conf_GetL2CacheMemory();
			if (L2Cache <= 0) L2Cache = 1024*1024; // 1M
			BlockSNP = (L2Cache - 8*1024) / (sizeof(double)*nSamp);
		}
		BlockSNP = (BlockSNP / 4) * 4;
		if (BlockSNP < 16) BlockSNP = 16;
	}

	/// init mutex objects
	void InitMutexObject()
	{
		_Mutex = plc_InitMutex();
	}
	/// destroy mutex objects
	void DoneMutexObject()
	{
		plc_DoneMutex(_Mutex); _Mutex = NULL;
	}


	// ****************** PCA covariate matrix ******************

	/// The temparory variables used in computing
	auto_ptr<int> PCA_gSum, PCA_gNum;
	/// The mean-adjusted genotype buffers
	CPCAMat<double> PCA_Mat;
	/// The temporary variable
	TdAlignPtr<double> tmpBuf;

	/// Thread variables
	const int N_MAX_THREAD = 256;
	IdMatTri PCA_Thread_MatIdx[N_MAX_THREAD];
	Int64 PCA_Thread_MatCnt[N_MAX_THREAD];

	/// Convert the raw genotypes to the mean-adjusted genotypes
	void _Do_PCA_ReadBlock(UInt8 *GenoBuf, long Start, long SNP_Cnt, void* Param)
	{
		// init ...
		const int n = MCWorkingGeno.Space.SampleNum();
		memset((void*)PCA_gSum.get(), 0, SNP_Cnt*sizeof(int));
		memset((void*)PCA_gNum.get(), 0, SNP_Cnt*sizeof(int));

		UInt8 *p;
		double *pf, *pp;

		// calculate the averages of genotypes for each SelSNP
		p = GenoBuf;
		pp = PCA_Mat.base();
		for (long iSample=0; iSample < n; iSample++)
		{
			pf = pp; pp += PCA_Mat.M();
			for (long iSNP=0; iSNP < SNP_Cnt; iSNP++)
			{
				if (*p < 3)
				{
					PCA_gSum.get()[iSNP] += *p;
					PCA_gNum.get()[iSNP]++;
				}
				*pf++ = *p++;
			}
		}

		for (long iSNP=0; iSNP < SNP_Cnt; iSNP++)
		{
			double &fv = tmpBuf.get()[iSNP];
			long gN = PCA_gNum.get()[iSNP];
			if (gN > 0)
				fv = (double)PCA_gSum.get()[iSNP] / gN;
			else
				fv = 0;
		}

		// Averaging and set missing values to ZERO
		pp = PCA_Mat.base(); // C@ij - 2p@j
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
				long gN = PCA_gNum.get()[iSNP];
				scale = (gN > 0) ? ((PCA_gSum.get()[iSNP] + 1.0) / (2*gN + 2)) : 0;
			} else
				scale = *pf / 2;
			if (0.0 < scale && scale < 1.0)
				*pf = 1.0 / sqrt(scale*(1.0-scale));
			else
				*pf = 0;
		}
		// (C@ij - 2p@j) / sqrt(p@j*(1-p@j))
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
	void _Do_PCA_ComputeCov(int ThreadIndex, long Start, long SNP_Cnt, void* Param)
	{
		double *base = (double*)Param;
		IdMatTri I = PCA_Thread_MatIdx[ThreadIndex];
		PCA_Mat.MulAdd(I, PCA_Thread_MatCnt[ThreadIndex], SNP_Cnt, base + I.Offset());
	}

    /// Calculate the genetic covariace
	void DoCovCalculate(CdMatTri<double> &PublicCov, int NumThread, const char *Info, bool verbose)
	{
		// Initialize ...
		PCA_gSum.reset(new int[BlockSNP]);
		PCA_gNum.reset(new int[BlockSNP]);
		PCA_Mat.Reset(PublicCov.N(), BlockSNP);
		tmpBuf.Reset(PCA_Mat.M());
		memset(PublicCov.get(), 0, sizeof(double)*PublicCov.Size());

		MCWorkingGeno.Progress.Info = Info;
		MCWorkingGeno.Progress.Show() = verbose;
		MCWorkingGeno.InitParam(true, true, BlockSNP);

		MCWorkingGeno.SplitJobs(NumThread, PublicCov.N(), PCA_Thread_MatIdx, PCA_Thread_MatCnt);
		MCWorkingGeno.Run(NumThread, &_Do_PCA_ReadBlock, &_Do_PCA_ComputeCov, PublicCov.get());
	}


	// ****************** SNP coefficients ******************

	// Correlation
	static double SNP_PC_Corr(double *pX, UInt8 *pY, long n)
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
				return conf_F64_NaN();
		} else
			return conf_F64_NaN();
	}

	/// The thread entry for the calculation of SNP coefficients
	void Entry_SNPCorrCalc(TdThread Thread, int ThreadIndex, void *Param)
	{
		// The number of working samples
		const long n = MCWorkingGeno.Space.SampleNum();
		auto_ptr<UInt8> GenoBlock(new UInt8[n*BlockSNP]);

		long _SNPstart, _SNPlen;
		while (RequireWork(GenoBlock.get(), _SNPstart, _SNPlen, false))
		{
			UInt8 *pGeno = GenoBlock.get();
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
	void DoSNPCoeffCalculate(int EigCnt, double *EigenVect, double *out_snpcorr,
		int NumThread, bool verbose, const char *Info)
	{
		// initialize mutex objects
		InitMutexObject();

		// Initialize progress information
		MCWorkingGeno.Progress.Info = Info;
		MCWorkingGeno.Progress.Show() = verbose;
		MCWorkingGeno.Progress.Init(MCWorkingGeno.Space.SNPNum());
		SNPStart = 0;
		OutputEigenDim = EigCnt;
		Out_Buffer = out_snpcorr; In_EigenVect = EigenVect;

		// Threads
		plc_DoBaseThread(Entry_SNPCorrCalc, NULL, NumThread);

		// destory the mutex objects
		DoneMutexObject();
	}


	// ****************** SNP loadings ******************

	/// Normalize the genotype
	/** \tparam tfloat    template type, it should be floating number type
	 *  \param Geno   genotype: 0 -- BB, 1 -- AB, 2 -- AA, 3 -- missing
	 *  \param NormalizedGeno  normalized genotype
	 *  \param NumGeno    the number of total genotype
	 *  \param AdjNormal  whether or not use Bayerian correction factor
	**/
	inline static void NormalizeGeno(UInt8 *Geno, double *NormalizedGeno, long NumGeno)
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
			UInt8 *pGeno = Geno;
			for (long i=NumGeno; i > 0; i--)
			{
				*NormalizedGeno++ = (*pGeno < 3) ? (*pGeno-ave)*scale : 0.0;
				pGeno++;
			}
		} else
			for (; NumGeno > 0; NumGeno--) *NormalizedGeno++ = 0.0;
	}

	/// The thread entry for the calculation of SNP loadings
	void Entry_SNPLoadingCalc(TdThread Thread, int ThreadIndex, void *Param)
	{
		// The number of working samples
		const long n = MCWorkingGeno.Space.SampleNum();
		auto_ptr<UInt8> GenoBlock(new UInt8[n*BlockSNP]);
		TdAlignPtr<double> NormalGeno(n);

		long _SNPstart, _SNPlen;
		while (RequireWork(GenoBlock.get(), _SNPstart, _SNPlen, false))
		{
			for (long iSNP=0; iSNP < _SNPlen; iSNP++)
			{
				// Normalize genotypes
				NormalizeGeno(GenoBlock.get() + iSNP*n, NormalGeno.get(), n);
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
			auto_ptr<UInt8> buf(new UInt8[nsnp]);
			auto_ptr<int> n(new int[nsnp]);
			for (int i=0; i < nsnp; i++)
				{ n.get()[i] = 0; OutFreq[i] = 0; }
			// for-loop for each sample
			for (int iSamp=0; iSamp < MCWorkingGeno.Space.SampleNum(); iSamp++)
			{
				MCWorkingGeno.Space.sampleRead(iSamp, 1, buf.get(), true);
				for (int i=0; i < nsnp; i++)
				{
					UInt8 &v = buf.get()[i];
					if (v <= 2) { OutFreq[i] += v; n.get()[i] ++; }
				}
			}
			// average
			for (int i=0; i < MCWorkingGeno.Space.SNPNum(); i++)
			{
				const int nn = n.get()[i];
				const double sum = OutFreq[i];
				double ave = sum / nn;
				double scale = BayesianNormal ? ((sum+1.0) / (2*nn+2)) : (0.5*ave);
				scale = ((0.0 < scale) && (scale < 1.0)) ? (1.0 / sqrt(scale*(1.0-scale))) : 0.0;
				OutFreq[i] = ave; OutScale[i] = scale;
			}
		} else {
			// initialize
			auto_ptr<UInt8> buf(new UInt8[MCWorkingGeno.Space.SampleNum()]);
			// for-loop for each snp
			for (int isnp=0; isnp < MCWorkingGeno.Space.SNPNum(); isnp++)
			{
				int n = 0;
				double sum = 0;
				MCWorkingGeno.Space.snpRead(isnp, 1, buf.get(), false);
				for (int i=0; i < MCWorkingGeno.Space.SampleNum(); i++)
				{
					UInt8 &v = buf.get()[i];
					if (v <= 2) { sum += v; n ++; }
				}
				double ave = sum / n;
				double scale = BayesianNormal ? ((sum+1.0) / (2*n+2)) : (0.5*ave);
				scale = ((0.0 < scale) && (scale < 1.0)) ? (1.0 / sqrt(scale*(1.0-scale))) : 0.0;
				OutFreq[isnp] = ave; OutScale[isnp] = scale;
			}
		}
	}

	/// Calculate the SNP loadings
	void DoSNPLoadingCalculate(double *EigenVal, int EigCnt, double *EigenVect,
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
		OutputEigenDim = EigCnt;
		Out_Buffer = out_snploading;

		// Array of eigenvectors
		const long n = MCWorkingGeno.Space.SampleNum();
		double Scale = double(n-1) / TraceXTX;
		_EigenVectBuf = new CPCAMat<double>(OutputEigenDim, n);
		for (long i=0; i < OutputEigenDim; i++)
		{
			vt<double>::Mul(_EigenVectBuf->base() + i*_EigenVectBuf->M(), EigenVect,
				sqrt(Scale / EigenVal[i]), n);
			EigenVect += n;
		}

		// Threads
		plc_DoBaseThread(Entry_SNPLoadingCalc, NULL, NumThread);

		// destory the mutex objects
		DoneMutexObject();
		delete _EigenVectBuf; _EigenVectBuf = NULL;
	}


	// ****************** Sample loadings ******************

	/// The thread entry for the calculation of genetic covariace matrix
	void Entry_SampLoadingCalc(TdThread Thread, int ThreadIndex, void *Param)
	{
		// The number of working samples
		const long n = MCWorkingGeno.Space.SNPNum();
		auto_ptr<UInt8> GenoBlock(new UInt8[n*BlockSamp]);
		auto_ptr<double> buf(new double[n]);

		long _SampStart, _SampLen;
		while (RequireWorkSamp(GenoBlock.get(), _SampStart, _SampLen, true))
		{
			// for-loop each sample
			for (long iS=0; iS < _SampLen; iS++)
			{
				// fill the buffer
				UInt8 *g = GenoBlock.get() + iS*n;
				double *p = buf.get();
				for (long i=0; i < n; i++, p++, g++)
					*p = (*g <= 2) ? (*g - In_AveFreq[i]) : 0.0;

				// for-loop each eigenvector
				double *out = Out_Buffer + _SampStart + iS;
				for (int i=0; i < OutputEigenDim; i++)
				{
					double *p1 = In_EigenVect + i, *p2 = buf.get();
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
		double *out_SampLoadings, int NumThread, const char *Info, bool verbose)
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

		// Threads
		SampStart = 0;
		plc_DoBaseThread(Entry_SampLoadingCalc, NULL, NumThread);

		// destory the mutex objects
		DoneMutexObject();
	}
}

#endif  /* _FuncPCA_H_ */
