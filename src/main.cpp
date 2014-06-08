// ===========================================================
//     _/_/_/   _/_/_/  _/_/_/_/    _/_/_/_/  _/_/_/   _/_/_/
//      _/    _/       _/             _/    _/    _/   _/   _/
//     _/    _/       _/_/_/_/       _/    _/    _/   _/_/_/
//    _/    _/       _/             _/    _/    _/   _/
// _/_/_/   _/_/_/  _/_/_/_/_/     _/     _/_/_/   _/_/
// ===========================================================
//
// main.cpp: Relatedness, Linkage Disequilibrium and Principal Component Analysis
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

#include <dGenGWAS.h>
#include <dVect.h>
#include <fstream>
#include <map>
#include <R.h>
#include <R_ext/Lapack.h>


using namespace std;
using namespace CoreArray;
using namespace CoreArray::Vectorization;
using namespace GWAS;

namespace PCA
{
	/// whether use Bayesian normalization
	extern bool BayesianNormal;

	void AutoDetectSNPBlockSize(int nSamp, bool Detect=true);

	void DoCovCalculate(CdMatTri<double> &PublicCov, int NumThread, const char *Info,
		bool verbose);

	void DoSNPCoeffCalculate(int EigCnt, double *EigenVect, double *out_snpcorr,
		int NumThread, bool verbose, const char *Info);

	void GetPCAFreqScale(double OutFreq[], double OutScale[]);
	void DoSNPLoadingCalculate(double *EigenVal, int EigCnt, double *EigenVect,
		double TraceXTX, double *out_snploading, int NumThread, bool verbose,
		const char *Info);

	void DoSampLoadingCalculate(double *Ave_Freq, double *Scale, int EigenCnt,
		double *SNP_Loadings, double *EigenVal, int Num, double TraceXTX,
		double *out_SampLoadings, int NumThread, const char *Info, bool verbose);
}

namespace IBS
{
	/// The structure of IBS states
	struct TIBSflag
	{
		UInt32 IBS0, IBS1, IBS2;
		TIBSflag() { IBS0 = IBS1 = IBS2 = 0; }
	};

	void AutoDetectSNPBlockSize(int nSamp, bool Detect=true);
	void DoIBSCalculate(CdMatTriDiag<TIBSflag> &PublicIBS, int NumThread,
		const char *Info, bool verbose);

	/// The structure of genetic distance
	struct TDistflag
	{
		Int64 SumGeno;
		double SumAFreq;
		TDistflag() { SumGeno = 0; SumAFreq = 0; }
	};

	/// Calculate the genetic distance matrix
	void DoDistCalculate(CdMatTriDiag<TDistflag> &PublicDist, int NumThread,
		const char *Info, bool verbose);
}

namespace IBD
{
	struct TIBD
	{
		double k0, k1;
		TIBD() { k0 = k1 = 0; }
	};

	/// The maximum number of iterations, 1000 by default
	extern long nIterMax;
	/// The reltol convergence tolerance, sqrt(machine.epsilon) by default
	extern double FuncRelTol;
	/// Whether constrict the estimates from PLINK
	extern bool KinshipConstraint;

	// PLINK method of moment
	void Init_EPrIBD_IBS(double *out_afreq);
	void Est_PLINK_Kinship(int IBS0, int IBS1, int IBS2, double &k0, double &k1,
		bool KinshipConstrict);

	// MLE
	extern long nIterMax;
	extern double FuncRelTol;
	extern int MethodMLE;
	extern bool Loglik_Adjust;

	void InitPackedGeno(void *buffer);
	void Do_MLE_IBD_Calculate(const double *AFreq, CdMatTriDiag<TIBD> &PublicIBD,
		CdMatTriDiag<int> *PublicNIter,
		double *out_AFreq, int NumThread, const char *Info, double *tmpAF, bool verbose);
	void Do_MLE_LogLik(const double *AFreq, const double *k0, const double *k1,
		double *tmp_AF, double *out_loglik);
	void Do_MLE_LogLik_k01(const double *AFreq, const double k0, const double k1,
		double *tmp_AF, double *out_loglik);
}

namespace LD
{
	extern int LD_Method;

	void InitPackedGeno();
	void DonePackedGeno();
	double calcLD(const int *snp1, const int *snp2, int cnt);
	void calcLD_mat(int nThread, double *out_LD);
	void calcLDpruning(int StartIdx, int *pos_bp, int slide_max_bp, int slide_max_n,
		const double LD_threshold, bool *out_SNP);
}

namespace INBREEDING
{
	// the individual inbreeding coefficients

	template<typename TYPE> static double _inb_mom(int n, TYPE snp[], double afreq[])
	{
		// get the initial values
		double F = 0;
		int nValid = 0;
		for (int i=0; i < n; i++)
		{
			int g = (int)(*snp++);
			if ((0 <= g) && (g <=2))
			{
				double p = afreq[i];
				double val = (g*g - (1+2*p)*g + 2*p*p) / (2*p*(1-p));
				if (conf_IsFinite64(val)) { F += val; nValid++; }
			}
		}
		if (nValid > 0) F /= nValid;
		return F;
	}

	template<typename TYPE> static double _inb_mom_ratio(int n, TYPE snp[], double afreq[])
	{
		// get the initial values
		double Den = 0, Num = 0;
		for (int i=0; i < n; i++)
		{
			int g = (int)(*snp++);
			if ((0 <= g) && (g <=2))
			{
				const double p = afreq[i];
				Num += g*g - (1+2*p)*g + 2*p*p;
				Den += 2*p*(1-p);
			}
		}
		return Num / Den;
	}

	template<typename TYPE> static double _inb_mle_loglik(double F, int n, TYPE snp[], double afreq[])
	{
		const double NaN = conf_F64_NaN();
		double rv = 0;
		for (int i=0; i < n; i++)
		{
			double p = afreq[i], val = NaN;
			switch (snp[i])
			{
				case 0:
					val = log((1-F)*(1-p)*(1-p) + F*(1-p)); break;
				case 1:
					val = log((1-F)*2*p*(1-p)); break;
				case 2:
					val = log((1-F)*p*p + F*p); break;
			}
			if (conf_IsFinite64(val)) rv += val;
		}
		return rv;
	}

	template<typename TYPE> static double _inb_mle(int n, TYPE snp[], double afreq[],
		const double reltol, int *out_iternum)
	{
		// initial value
		double F = _inb_mom_ratio(n, snp, afreq);
		if (conf_IsFinite64(F))
		{
			if (F < 0.001) F = 0.001;
			if (F > 1 - 0.001) F = 1 - 0.001;

			// MLE updating ...
			const int n_iter_max = 10000;

			double LogLik = _inb_mle_loglik(F, n, snp, afreq);
			double contol = fabs(LogLik) * reltol;
			int iter;
			for (iter=1; iter <= n_iter_max; iter++)
			{
				double OldLogLik = LogLik, sum = 0;
				int m = 0;
				TYPE *ps = snp;
				for (int i=0; i < n; i++)
				{
					double p = afreq[i], tmp;
					switch (*ps++)
					{
						case 0:
							tmp = F / (F + (1-p)*(1-F));
							if (conf_IsFinite64(tmp)) { sum += tmp; m++; }
							break;
						case 1:
							m++; break;
						case 2:
							tmp = F / (F + p*(1-F));
							if (conf_IsFinite64(tmp)) { sum += tmp; m++; }
							break;
					}
				}
				F = sum / m;
				LogLik = _inb_mle_loglik(F, n, snp, afreq);
				if (fabs(LogLik - OldLogLik) <= contol) break;
			}
			if (out_iternum) *out_iternum = iter;
		}
		return F;
	}
}



#define LongBool int


#ifdef COREARRAY_GNUG
#  ifdef COREARRAY_WINDOWS
#    define DLLEXPORT __attribute__((dllexport))
#  else
#    define DLLEXPORT
#  endif
#else
#  define DLLEXPORT __declspec(dllexport)
#endif

extern "C"
{
	#define CORETRY			try {
	#define CORECATCH(cmd)	} \
		catch (exception &E) { \
			gds_LastError() = E.what(); \
			cmd; \
		} \
		catch (const char *E) { \
			gds_LastError() = E; \
			cmd; \
		}


// ===========================================================
// the private functions
// ===========================================================

inline static void RStrAgn(const char *Text, char **rstr)
{
	*rstr = R_alloc(strlen(Text)+1, 1);
	if (*rstr == NULL)
		throw ErrCoreArray("R_alloc return NULL!");
	strcpy(*rstr, Text);
}



// ===========================================================
// the public functions
// ===========================================================

/// initialize the package
DLLEXPORT void gnrInit(char **lib_fn, char **rstr, LongBool *sse2)
{
	try {
		GDSInterface::InitGDSInterface(*lib_fn);
		RStrAgn("", rstr);

		#ifdef COREARRAY_SIMD_SSE2
		*sse2 = true;
		#else
		*sse2 = false;
		#endif
	}
	catch (exception &E)
	{
		RStrAgn(E.what(), rstr);
	}
}

/// finalize the package
DLLEXPORT void gnrDone()
{
	GDSInterface::DoneGDSInterface();
}



// the functions for the SNP genotype working space

/// set the genotype node
DLLEXPORT void gnrSetGenoSpace(TdSequenceX *Node,
	LongBool *SelSamp, LongBool *SelSampFlag, LongBool *SelSNP, LongBool *SelSNPFlag,
	int *out_nsnp, int *out_nsamp,
	LongBool *out_err)
{
	CORETRY
		MCWorkingGeno.Space.SetGeno(*Node, false);
		if (*SelSampFlag)
		{
			for (int i=0; i < MCWorkingGeno.Space.TotalSampleNum(); i++)
			{
				MCWorkingGeno.Space.SampleSelection()[i] = (SelSamp[i] != 0);
			}
		}
		if (*SelSNPFlag)
		{
			for (int i=0; i < MCWorkingGeno.Space.TotalSNPNum(); i++)
			{
				MCWorkingGeno.Space.SNPSelection()[i] = (SelSNP[i] != 0);
			}
		}
		MCWorkingGeno.Space.InitSelection();
		*out_nsnp = MCWorkingGeno.Space.SNPNum();
		*out_nsamp = MCWorkingGeno.Space.SampleNum();
		*out_err = 0;
	CORECATCH(*out_err = 1)
}

/// get the dimension of SNP genotypes
DLLEXPORT void gnrGetGenoDim(int *out_nsnp, int *out_nsamp)
{
	*out_nsnp = MCWorkingGeno.Space.SNPNum();
	*out_nsamp = MCWorkingGeno.Space.SampleNum();
}

/// set the SNP selection on the genotype node
DLLEXPORT void gnrSetGenoSelSNP(LongBool SNP_Flag[], LongBool *out_err)
{
	CORETRY
		auto_ptr<CBOOL> buf(new CBOOL[MCWorkingGeno.Space.SNPNum()]);
		for (int i=0; i < MCWorkingGeno.Space.SNPNum(); i++)
			buf.get()[i] = (SNP_Flag[i]!=0);
		MCWorkingGeno.Space.Set_SNPSelection(buf.get());
		*out_err = 0;
	CORECATCH(*out_err = 1)
}

/** to compute allele freq., minor allele freq. and missing rates
 *  \param remove_mono    whether remove monomorphic snps or not
 *  \param maf            the threshold of minor allele frequencies, keeping ">= maf"
 *  \param missrate       the threshold of missing rates, keeping "<= missing.rate"
**/
DLLEXPORT void gnrSelSNP_Base(LongBool *remove_mono, double *maf, double *missrate,
	int *out_num, LongBool *out_selection, LongBool *out_err)
{
	CORETRY
		const int n = MCWorkingGeno.Space.SNPNum();
		auto_ptr<CBOOL> sel(new CBOOL[n]);
		*out_num = MCWorkingGeno.Space.Select_SNP_Base(*remove_mono!=0,
			*maf, *missrate, sel.get());
		for (int i=0; i < n; i++) out_selection[i] = sel.get()[i];
		*out_err = 0;
	CORECATCH(*out_err = 1)
}

/// to compute allele freq., minor allele freq. and missing rates
DLLEXPORT void gnrSNPFreq(
	double out_allelefreq[], double out_minorfreq[], double out_missrate[],
	LongBool *out_err)
{
	CORETRY
		MCWorkingGeno.Space.GetAlleleFreqs(out_allelefreq);
		for (int i=0; i < MCWorkingGeno.Space.SNPNum(); i++)
			out_minorfreq[i] = min(out_allelefreq[i], 1-out_allelefreq[i]);
		MCWorkingGeno.Space.GetMissingRates(out_missrate);
		*out_err = 0;
	CORECATCH(*out_err = 1)
}

/// to compute sample missing rates
DLLEXPORT void gnrSampFreq(double out_missrate[], LongBool *out_err)
{
	CORETRY
		MCWorkingGeno.Space.GetSampMissingRates(out_missrate);
		*out_err = 0;
	CORECATCH(*out_err = 1)
}

/// to cache the genotype data
DLLEXPORT void gnrCacheGeno(double *out_GenoSum, LongBool *out_err)
{
	CORETRY
		if (out_GenoSum)
			*out_GenoSum = MCWorkingGeno.Space.GenoSum();
		if (out_err) *out_err = 0;
	CORECATCH(if (out_err) *out_err = 1)
}


/// add genotype buffer
DLLEXPORT void gnrInitGenoBuffer(LongBool *SNPorSamp, int *AF, int *out_obj, LongBool *out_err)
{
	CORETRY
		CdBufSpace *p = new CdBufSpace(MCWorkingGeno.Space,
			*SNPorSamp, CdBufSpace::TAccessFlag(*AF));
		memmove((void*)out_obj, (void*)&p, sizeof(CdBufSpace *));
		*out_err = 0;
	CORECATCH(*out_err = 1)
}

/// done genotype buffer
DLLEXPORT void gnrDoneGenoBuffer(int *buf_obj, LongBool *out_err)
{
	CORETRY
		CdBufSpace *p;
		memmove((void*)&p, (void*)buf_obj, sizeof(CdBufSpace *));
		delete p;
		*out_err = 0;
	CORECATCH(*out_err = 1)
}

/// get genotype
DLLEXPORT void gnrGetGenoBuffer(int *buf_obj, int *idx, int *out_buf, LongBool *out_err)
{
	CORETRY
		CdBufSpace *p;
		memmove((void*)&p, (void*)buf_obj, sizeof(CdBufSpace *));
		UInt8 *s = p->ReadGeno(*idx);
		for (long i=0; i < p->BufElmSize(); i++)
			out_buf[i] = s[i];
		*out_err = 0;
	CORECATCH(*out_err = 1)
}

/// copy genotypes
DLLEXPORT void gnrCopyGeno(TdSequenceX *Node, LongBool *snpfirstorder, LongBool *out_err)
{
	CORETRY
		*out_err = 1;
		if (*snpfirstorder)
		{
			CoreArray::Int32 cnt[2] = { 1, MCWorkingGeno.Space.SNPNum() };
			CdBufSpace buf(MCWorkingGeno.Space, false, CdBufSpace::acInc);
			for (long i=0; i < buf.IdxCnt(); i++)
			{
				UInt8 *p = buf.ReadGeno(i);
				CoreArray::Int32 st[2] = {i, 0};
				if (!gds_wData(*Node, st, cnt, p, svUInt8))
					return;
			}
		} else {
			CoreArray::Int32 cnt[2] = { 1, MCWorkingGeno.Space.SampleNum() };
			CdBufSpace buf(MCWorkingGeno.Space, true, CdBufSpace::acInc);
			for (long i=0; i < buf.IdxCnt(); i++)
			{
				UInt8 *p = buf.ReadGeno(i);
				CoreArray::Int32 st[2] = {i, 0};
				if (!gds_wData(*Node, st, cnt, p, svUInt8))
					return;
			}
		}
		*out_err = 0;
	CORECATCH(*out_err = 1)
}

/// copy genotypes to a memory buffer
DLLEXPORT void gnrCopyGenoMem(int *membuf, LongBool *snpfirstorder, LongBool *out_err)
{
	CORETRY
		*out_err = 1;
		if (*snpfirstorder)
		{
			int cnt = MCWorkingGeno.Space.SNPNum();
			CdBufSpace buf(MCWorkingGeno.Space, false, CdBufSpace::acInc);
			for (long i=0; i < buf.IdxCnt(); i++)
			{
				UInt8 *p = buf.ReadGeno(i);
				for (int k=0; k < cnt; k++) *membuf++ = *p++;
			}
		} else {
			int cnt = MCWorkingGeno.Space.SampleNum();
			CdBufSpace buf(MCWorkingGeno.Space, true, CdBufSpace::acInc);
			for (long i=0; i < buf.IdxCnt(); i++)
			{
				UInt8 *p = buf.ReadGeno(i);
				for (int k=0; k < cnt; k++) *membuf++ = *p++;
			}
		}
		*out_err = 0;
	CORECATCH(*out_err = 1)
}

/// Append genotypes
DLLEXPORT void gnrAppendGenoSpace(TdSequenceX *Node, LongBool *snpfirstorder, LongBool *out_err)
{
	CORETRY
		*out_err = 1;
		if (*snpfirstorder)
		{
			CdBufSpace buf(MCWorkingGeno.Space, false, CdBufSpace::acInc);
			for (long i=0; i < buf.IdxCnt(); i++)
			{
				UInt8 *p = buf.ReadGeno(i);
				if (!gds_AppendData(*Node, MCWorkingGeno.Space.SNPNum(), p, svUInt8))
					return;
			}
		} else {
			CdBufSpace buf(MCWorkingGeno.Space, true, CdBufSpace::acInc);
			for (long i=0; i < buf.IdxCnt(); i++)
			{
				UInt8 *p = buf.ReadGeno(i);
				if (!gds_AppendData(*Node, MCWorkingGeno.Space.SampleNum(), p, svUInt8))
					return;
			}
		}
		*out_err = 0;
	CORECATCH(*out_err = 1)
}

/// Append genotypes with switch strand
DLLEXPORT void gnrAppendGenoSpaceStrand(TdSequenceX *Node, LongBool *snpfirstorder,
	LongBool StrandFlag[], LongBool *out_err)
{
	CORETRY
		*out_err = 1;
		if (*snpfirstorder)
		{
			const int n = MCWorkingGeno.Space.SNPNum();
			CdBufSpace buf(MCWorkingGeno.Space, false, CdBufSpace::acInc);
			for (long i=0; i < buf.IdxCnt(); i++)
			{
				UInt8 *p = buf.ReadGeno(i);
				for (int k=0; k < n; k++)
					if (StrandFlag[k] && (p[k]<=2)) p[k] = 2 - p[k];
				if (!gds_AppendData(*Node, n, p, svUInt8))
					return;
			}
		} else {
			const int n = MCWorkingGeno.Space.SampleNum();
			CdBufSpace buf(MCWorkingGeno.Space, true, CdBufSpace::acInc);
			for (long i=0; i < buf.IdxCnt(); i++)
			{
				UInt8 *p = buf.ReadGeno(i);
				if (StrandFlag[i])
				{
					for (int k=0; k < n; k++)
						if (p[k] <= 2) p[k] = 2 - p[k];
				}
				if (!gds_AppendData(*Node, n, p, svUInt8))
					return;
			}
		}
		*out_err = 0;
	CORECATCH(*out_err = 1)
}




// the functions for Principal Component Analysis (PCA)

/// to compute the eigenvalues and eigenvectors
DLLEXPORT void gnrPCA(int *EigenCnt, int *NumThread, LongBool *_BayesianNormal,
	LongBool *NeedGenMat, LongBool *GenMat_Only, LongBool *Verbose, LongBool *DataCache,
	double *out_Eigenvalues, double *out_Eigenvectors,
	double *out_TraceXTX, double *out_GenMat, LongBool *out_err)
{
	CORETRY
		// ******** To cache the genotype data ********
		if (*DataCache)
		{
			double GenoSum=0;
			gnrCacheGeno(&GenoSum, NULL);
			if (*Verbose)
				Rprintf("PCA:\tthe sum of all working genotypes = %.0f\n", GenoSum);
		}

		// ******** The calculation of genetic covariance matrix ********

		// the number of samples
		const int n = MCWorkingGeno.Space.SampleNum();
		// set parameters
		PCA::BayesianNormal = (_BayesianNormal!=0);
		PCA::AutoDetectSNPBlockSize(n);
		// the upper-triangle genetic covariance matrix
		CdMatTri<double> Cov(n);

		// Calculate the genetic covariace
		PCA::DoCovCalculate(Cov, *NumThread, "PCA:", *Verbose);
		// Normalize
		double TraceXTX = Cov.Trace();
		double scale = double(n-1) / TraceXTX;
		vt<double, av16Align>::Mul(Cov.get(), Cov.get(), scale, Cov.Size());
		*out_TraceXTX = TraceXTX;
		if (*NeedGenMat) Cov.SaveTo(out_GenMat);

		// ******** The calculation of eigenvectors and eigenvalues ********

		if (!(*GenMat_Only))
		{
			const size_t NN = n;
			auto_ptr<double> tmp_Work(new double[NN*3]);
			auto_ptr<double> tmp_EigenVec(new double[NN*NN]);

			vt<double>::Sub(Cov.get(), 0.0, Cov.get(), Cov.Size());
			if (*Verbose)
				Rprintf("PCA:\t%s\tBegin (eigenvalues and eigenvectors)\n", NowDateToStr().c_str());
			{
				int info = 0;
				F77_NAME(dspev)("V", "L", &n, Cov.get(), out_Eigenvalues,
					tmp_EigenVec.get(), &n, tmp_Work.get(), &info);
				if (info != 0)
					Rprintf("LAPACK::SPEV error!");
			}

			// output eigenvalues
			vt<double>::Sub(out_Eigenvalues, 0.0, out_Eigenvalues, n);

			// output eigenvectors
			double *p = tmp_EigenVec.get();
			for (int i=0; i < *EigenCnt; i++)
			{
				memmove(out_Eigenvectors, p, sizeof(double)*n);
				out_Eigenvectors += n; p += n;
			}
			if (*Verbose)
				Rprintf("PCA:\t%s\tEnd (eigenvalues and eigenvectors)\n", NowDateToStr().c_str());
		}

		// output
		*out_err = 0;
	CORECATCH(*out_err = 1)
}

/// to calculate the SNP correlations
DLLEXPORT void gnrPCACorr(int *DimCnt, double *EigenVect, int *NumThread,
	LongBool *Verbose, LongBool *DataCache, double *out_snpcorr, LongBool *out_err)
{
	CORETRY
		// ******** To cache the genotype data ********
		if (*DataCache)
		{
			double GenoSum=0;
			gnrCacheGeno(&GenoSum, NULL);
			if (*Verbose)
				Rprintf("SNP Correlations:\tthe sum of all working genotypes = %.0f\n", GenoSum);
		}

		// ******** To compute the snp correlation ********
		PCA::AutoDetectSNPBlockSize(MCWorkingGeno.Space.SampleNum());
		PCA::DoSNPCoeffCalculate(DimCnt[1], EigenVect, out_snpcorr, *NumThread,
			*Verbose!=0, "SNP Correlations:");

		// output
		*out_err = 0;
	CORECATCH(*out_err = 1)
}

/// to calculate the SNP loadings
DLLEXPORT void gnrPCASNPLoading(double *EigenVal, int *DimCnt, double *EigenVect,
	double *TraceXTX, int *NumThread, LongBool *Bayesian, LongBool *Verbose,
	LongBool *DataCache, double *out_snploading, double *out_afreq, double *out_scale,
	LongBool *out_err)
{
	CORETRY
		// ******** To cache the genotype data ********
		if (*DataCache)
		{
			double GenoSum=0;
			gnrCacheGeno(&GenoSum, NULL);
			if (*Verbose)
				Rprintf("SNP Loadings:\tthe sum of all working genotypes = %.0f\n", GenoSum);
		}

		// ******** To compute the snp correlation ********
		PCA::AutoDetectSNPBlockSize(MCWorkingGeno.Space.SampleNum());
		PCA::BayesianNormal = (Bayesian!=0);
		PCA::GetPCAFreqScale(out_afreq, out_scale);
		PCA::DoSNPLoadingCalculate(EigenVal, DimCnt[1], EigenVect, *TraceXTX,
			out_snploading, *NumThread, *Verbose, "SNP Loadings:");

		// output
		*out_err = 0;
	CORECATCH(*out_err = 1)
}

/// to calculate the sample loadings from SNP loadings
DLLEXPORT void gnrPCASampLoading(int *Num, double *EigenVal, int *EigenCnt,
	double *SNPLoadings, double *TraceXTX, double *AveFreq, double *Scale,
	int *NumThread, LongBool *Verbose, LongBool *DataCache,
	double *out_samploading, LongBool *out_err)
{
	CORETRY
		// ******** To cache the genotype data ********
		if (*DataCache)
		{
			double GenoSum=0;
			gnrCacheGeno(&GenoSum, NULL);
			if (*Verbose)
				Rprintf("Sample Loadings:\tthe sum of all working genotypes = %.0f\n", GenoSum);
		}

		// ******** To compute the snp correlation ********
		PCA::DoSampLoadingCalculate(AveFreq, Scale, *EigenCnt, SNPLoadings,
			EigenVal, *Num, *TraceXTX, out_samploading,
			*NumThread, "Sample Loadings:", *Verbose);

		// output
		*out_err = 0;
	CORECATCH(*out_err = 1)
}



// the functions for Identity by state (IBS)

/// to compute the average IBS
DLLEXPORT void gnrIBSAve(LongBool *Verbose, LongBool *DataCache, int *NumThread,
	double *out_IBSMat, LongBool *out_err)
{
	CORETRY
		// ******** To cache the genotype data ********
		if (*DataCache)
		{
			double GenoSum=0;
			gnrCacheGeno(&GenoSum, NULL);
			if (*Verbose)
				Rprintf("IBS:\tthe sum of all working genotypes = %.0f\n", GenoSum);
		}

		// ******** The calculation of genetic covariance matrix ********

		// the number of samples
		const int n = MCWorkingGeno.Space.SampleNum();
		// to detect the block size
		IBS::AutoDetectSNPBlockSize(n);
		// the upper-triangle genetic covariance matrix
		CdMatTriDiag<IBS::TIBSflag> IBS(IBS::TIBSflag(), n);

		// Calculate the IBS matrix
		IBS::DoIBSCalculate(IBS, *NumThread, "IBS:", *Verbose);

		// output
		IBS::TIBSflag *p = IBS.get();
		for (int i=0; i < n; i++)
		{
			out_IBSMat[i*n + i] = 1.0;
			for (int j=i+1; j < n; j++, p++)
			{
				out_IBSMat[i*n + j] = out_IBSMat[j*n + i] =
					double(0.5*p->IBS1 + p->IBS2) / (p->IBS0 + p->IBS1 + p->IBS2);
			}
		}
		// output
		*out_err = 0;
	CORECATCH(*out_err = 1)
}

/// to compute the average IBS
DLLEXPORT void gnrIBSNum(LongBool *Verbose, LongBool *DataCache, int *NumThread,
	int *out_IBS0, int *out_IBS1, int *out_IBS2, LongBool *out_err)
{
	CORETRY
		// ******** To cache the genotype data ********
		if (*DataCache)
		{
			double GenoSum=0;
			gnrCacheGeno(&GenoSum, NULL);
			if (*Verbose)
				Rprintf("IBS:\tthe sum of all working genotypes = %.0f\n", GenoSum);
		}

		// ******** The calculation of genetic covariance matrix ********

		// the number of samples
		const int n = MCWorkingGeno.Space.SampleNum();
		// to detect the block size
		IBS::AutoDetectSNPBlockSize(n);
		// the upper-triangle genetic covariance matrix
		CdMatTriDiag<IBS::TIBSflag> IBS(IBS::TIBSflag(), n);

		// Calculate the IBS matrix
		auto_ptr<int> I(new int[n]);
		MCWorkingGeno.Space.GetSampValidNum(I.get());
		IBS::DoIBSCalculate(IBS, *NumThread, "IBS:", *Verbose);

		// output
		IBS::TIBSflag *p = IBS.get();
		for (int i=0; i < n; i++)
		{
			out_IBS0[i*n + i] = 0; out_IBS1[i*n + i] = 0;
			out_IBS2[i*n + i] = I.get()[i];
			for (int j=i+1; j < n; j++, p++)
			{
				out_IBS0[i*n + j] = out_IBS0[j*n + i] = p->IBS0;
				out_IBS1[i*n + j] = out_IBS1[j*n + i] = p->IBS1;
				out_IBS2[i*n + j] = out_IBS2[j*n + i] = p->IBS2;
			}
		}
		// output
		*out_err = 0;
	CORECATCH(*out_err = 1)
}


// the functions for genetic distances

/// to compute the average genetic distance
DLLEXPORT void gnrDist(LongBool *Verbose, LongBool *DataCache, int *NumThread,
	double *out_Dist, LongBool *out_err)
{
	CORETRY
		// ******** To cache the genotype data ********
		if (*DataCache)
		{
			double GenoSum=0;
			gnrCacheGeno(&GenoSum, NULL);
			if (*Verbose)
				Rprintf("Genetic Distance:\tthe sum of all working genotypes = %.0f\n", GenoSum);
		}

		// ******** The calculation of genetic covariance matrix ********

		// the number of samples
		const int n = MCWorkingGeno.Space.SampleNum();
		// to detect the block size
		IBS::AutoDetectSNPBlockSize(n);
		// the upper-triangle genetic covariance matrix
		CdMatTriDiag<IBS::TDistflag> Dist(IBS::TDistflag(), n);

		// Calculate the genetic distance matrix
		IBS::DoDistCalculate(Dist, *NumThread, "Genetic Distance:", *Verbose);

		// output
		IBS::TDistflag *p = Dist.get();
		for (int i=0; i < n; i++)
		{
			out_Dist[i*n + i] = 0;
			for (int j=i+1; j < n; j++, p++)
				out_Dist[i*n + j] = out_Dist[j*n + i] = (p->SumGeno / p->SumAFreq);
		}
		// output
		*out_err = 0;
	CORECATCH(*out_err = 1)
}


// the functions for identity-by-descent (IBD)

/// to compute the IBD coefficients by PLINK method of moment
DLLEXPORT void gnrIBD_PLINK(LongBool *Verbose, LongBool *DataCache, int *NumThread,
	LongBool *KinshipConstrict, double *out_k0, double *out_k1, double *out_afreq,
	LongBool *out_err)
{
	CORETRY
		// ******** to cache the genotype data ********
		if (*DataCache)
		{
			double GenoSum=0;
			gnrCacheGeno(&GenoSum, NULL);
			if (*Verbose)
				Rprintf("PLINK IBD:\tthe sum of all working genotypes = %.0f\n", GenoSum);
		}

		// the number of samples
		const int n = MCWorkingGeno.Space.SampleNum();
		// to detect the block size
		IBS::AutoDetectSNPBlockSize(n);

		// ******** PLINK method of moment ********

		// initialize the internal matrix
		IBD::Init_EPrIBD_IBS(out_afreq);
		// the upper-triangle genetic covariance matrix
		CdMatTriDiag<IBS::TIBSflag> IBS(IBS::TIBSflag(), n);
		// Calculate the IBS matrix
		IBS::DoIBSCalculate(IBS, *NumThread, "PLINK IBD:", *Verbose);

		// output
		IBS::TIBSflag *p = IBS.get();
		for (int i=0; i < n; i++)
		{
			out_k0[i*n + i] = out_k1[i*n + i] = 0;
			for (int j=i+1; j < n; j++, p++)
			{
				double k0, k1;
				IBD::Est_PLINK_Kinship(p->IBS0, p->IBS1, p->IBS2, k0, k1, *KinshipConstrict);
				out_k0[i*n + j] = out_k0[j*n + i] = k0;
				out_k1[i*n + j] = out_k1[j*n + i] = k1;
			}
		}

		// output
		*out_err = 0;
	CORECATCH(*out_err = 1)
}

/// get an error message
DLLEXPORT void gnrIBD_SizeInt(int *out_size, int *out_snp_size)
{
	long nSamp = MCWorkingGeno.Space.SampleNum();
	long nPackedSNP = (MCWorkingGeno.Space.SNPNum() % 4 > 0) ?
		(MCWorkingGeno.Space.SNPNum()/4 + 1) : (MCWorkingGeno.Space.SNPNum()/4);
	long nTotal = nSamp * nPackedSNP;
	*out_size = nTotal/sizeof(int) + ((nTotal % sizeof(int) > 0) ? 1 : 0);
	*out_snp_size = 4*nPackedSNP;
}

/// to compute the IBD coefficients by MLE
DLLEXPORT void gnrIBD_MLE(double *AlleleFreq, LongBool *UseSpecificAFreq,
	LongBool *KinshipConstraint, int *MaxIterCnt, double *RelTol,
	LongBool *CoeffCorrect, int *method, LongBool *Verbose, LongBool *DataCache,
	int *NumThread, LongBool *ifoutn,
	int *tmp_buffer, double *tmp_AF,
	double *out_k0, double *out_k1, double *out_afreq, int *out_niter,
	LongBool *out_err)
{
	CORETRY
		// ******** to cache the genotype data ********
		if (*DataCache)
		{
			double GenoSum=0;
			gnrCacheGeno(&GenoSum, NULL);
			if (*Verbose)
				Rprintf("MLE IBD:\tthe sum of all working genotypes = %.0f\n", GenoSum);
		}

		// ******** MLE IBD ********

		// initialize the internal matrix
		IBD::Init_EPrIBD_IBS(NULL);
		// initialize the packed genotypes
		IBD::InitPackedGeno(tmp_buffer);

		IBD::nIterMax = *MaxIterCnt; IBD::FuncRelTol = *RelTol;
		IBD::MethodMLE = *method; IBD::Loglik_Adjust = *CoeffCorrect;

		// the upper-triangle genetic covariance matrix
		const int n = MCWorkingGeno.Space.SampleNum();
		CdMatTriDiag<IBD::TIBD> IBD(IBD::TIBD(), n);
		CdMatTriDiag<int> niter;
		if (*ifoutn) niter.Reset(n);
		// Calculate the IBD matrix
		IBD::Do_MLE_IBD_Calculate(
			*UseSpecificAFreq ? AlleleFreq : NULL,
			IBD, (*ifoutn) ? &niter : NULL, out_afreq,
			*NumThread, "MLE IBD:", tmp_AF, *Verbose);

		// output
		IBD::TIBD *p = IBD.get();
		int *pn = niter.get();
		for (int i=0; i < n; i++)
		{
			out_k0[i*n + i] = out_k1[i*n + i] = 0;
			for (int j=i+1; j < n; j++, p++)
			{
				out_k0[i*n + j] = out_k0[j*n + i] = p->k0;
				out_k1[i*n + j] = out_k1[j*n + i] = p->k1;
				if (*ifoutn)
					out_niter[i*n + j] = out_niter[j*n + i] = *pn++;
			}
		}
		// output
		*out_err = 0;
	CORECATCH(*out_err = 1)
}

/// to compute log likelihood of MLE
DLLEXPORT void gnrIBD_LogLik(double *AlleleFreq, double *k0, double *k1,
	int *tmp_buffer, double *tmp_AF, double *out_loglik, LongBool *out_err)
{
	CORETRY
		// ******** MLE IBD ********
		// initialize the packed genotypes
		IBD::InitPackedGeno(tmp_buffer);
		// call
		IBD::Do_MLE_LogLik(AlleleFreq, k0, k1, tmp_AF, out_loglik);
		// output
		*out_err = 0;
	CORECATCH(*out_err = 1)
}

/// to compute log likelihood of MLE
DLLEXPORT void gnrIBD_LogLik_k01(double *AlleleFreq, double *k0, double *k1,
	int *tmp_buffer, double *tmp_AF, double *out_loglik, LongBool *out_err)
{
	CORETRY
		// ******** MLE IBD ********
		// initialize the packed genotypes
		IBD::InitPackedGeno(tmp_buffer);
		// call
		IBD::Do_MLE_LogLik_k01(AlleleFreq, *k0, *k1, tmp_AF, out_loglik);
		// output
		*out_err = 0;
	CORECATCH(*out_err = 1)
}


/// the functions for Linkage Disequilibrium (LD) analysis
DLLEXPORT void gnrLDpair(int *snp1, int *snp2, int *len, int *method,
	double *out_LD, LongBool *out_err)
{
	CORETRY
		LD::LD_Method = *method;
		*out_LD = LD::calcLD(snp1, snp2, *len);
		// output
		*out_err = 0;
	CORECATCH(*out_err = 1)
}

/// to compute the IBD coefficients by MLE
DLLEXPORT void gnrLDMat(int *method, LongBool *Verbose, LongBool *DataCache,
	int *NumThread, double *out_LD, LongBool *out_err)
{
	CORETRY
		// ******** to cache the genotype data ********
		if (*DataCache)
		{
			double GenoSum=0;
			gnrCacheGeno(&GenoSum, NULL);
			if (*Verbose)
				Rprintf("MLE IBD:\tthe sum of all working genotypes = %.0f\n", GenoSum);
		}

		// initialize the packed genotypes
		LD::InitPackedGeno();
		LD::LD_Method = *method;
		LD::calcLD_mat(*NumThread, out_LD);
		// output
		*out_err = 0;
	CORECATCH(*out_err = 1)
}

/// to compute the IBD coefficients by MLE
DLLEXPORT void gnrLDpruning(int *StartIdx, int *pos_bp, int *slide_max_bp, int *slide_max_n,
	double *LD_threshold, int *method, LongBool *out_SNPprune, LongBool *out_err)
{
	CORETRY
		auto_ptr<bool> flag(new bool[MCWorkingGeno.Space.SNPNum()]);
		LD::LD_Method = *method;
		LD::calcLDpruning(*StartIdx, pos_bp, *slide_max_bp, *slide_max_n, *LD_threshold,
			flag.get());
		for (int i=0; i < MCWorkingGeno.Space.SNPNum(); i++)
			out_SNPprune[i] = flag.get()[i];
		// output
		*out_err = 0;
	CORECATCH(*out_err = 1)
}


// the individual inbreeding coefficients

// to compute the inbreeding coefficient
DLLEXPORT void gnrIndInbCoef(int *n, int *snp, double *afreq, double *reltol, double *out_IC,
	LongBool *out_err)
{
	CORETRY
		*out_IC = INBREEDING::_inb_mle<int>(*n, snp, afreq, *reltol, NULL);
		*out_err = 0;
	CORECATCH(*out_err = 1)
}


// to compute the inbreeding coefficient
DLLEXPORT void gnrIndInb(double *allele_freq, int *method, double *reltol, double out_coeff[],
	int out_iternum[], LongBool *out_err)
{
	CORETRY
		const int n = MCWorkingGeno.Space.SNPNum();
		CdBufSpace buf(MCWorkingGeno.Space, false, CdBufSpace::acInc);
		for (long i=0; i < buf.IdxCnt(); i++)
		{
			UInt8 *p = buf.ReadGeno(i);
			switch (*method)
			{
				case 1:
					out_coeff[i] = INBREEDING::_inb_mom_ratio(n, p, allele_freq); break;
				case 2:
					out_coeff[i] = INBREEDING::_inb_mom(n, p, allele_freq); break;
				case 3:
					out_coeff[i] = INBREEDING::_inb_mle(n, p, allele_freq, *reltol, &out_iternum[i]); break;
			}
		}
		*out_err = 0;
	CORECATCH(*out_err = 1)
}






// conversion

/// to convert from GDS to PLINK PED
DLLEXPORT void gnrConvGDS2PED(char **pedfn, char **SampID, int *Sex,
	LongBool *verbose, LongBool *out_err)
{
	CORETRY
		MCWorkingGeno.Progress.Info = "\t\tOutput: ";
		MCWorkingGeno.Progress.Show() = *verbose;
		MCWorkingGeno.Progress.Init(MCWorkingGeno.Space.SampleNum());

		ofstream file(*pedfn);
		if (!file.good())
			throw ErrCoreArray("Fail to create the file %s.", *pedfn);
		CdBufSpace buf(MCWorkingGeno.Space, false, CdBufSpace::acInc);
		for (long i=0; i < buf.IdxCnt(); i++)
		{
			file << "0\t" << SampID[i] << "\t0\t0\t" << Sex[i] << "\t-9";
			UInt8 *g = buf.ReadGeno(i);
			for (long j=0; j < MCWorkingGeno.Space.SNPNum(); j++, g++)
			{
				const char *s = (*g==0) ? "B B" : ((*g==1)?"A B": ((*g==2)?"A A":"0 0"));
				file << "\t" << s;
			}
			file << endl;
			MCWorkingGeno.Progress.Forward(1);
		}
		*out_err = 0;
	CORECATCH(*out_err = 1)
}

/// to convert from GDS to PLINK BED
DLLEXPORT void gnrConvGDS2BED(char **bedfn, LongBool *SNPOrder,
	LongBool *verbose, LongBool *out_err)
{
	CORETRY
		MCWorkingGeno.Progress.Info = "\t\tOutput: ";
		MCWorkingGeno.Progress.Show() = *verbose;
		MCWorkingGeno.Progress.Init(MCWorkingGeno.Space.SampleNum());

		ofstream file(*bedfn, ios::binary);
		if (!file.good())
			throw ErrCoreArray("Fail to create the file %s.", *bedfn);
		// output prefix
		{
			char prefix[3];
			prefix[0] = 0x6C; prefix[1] = 0x1B;
			prefix[2] = (*SNPOrder) ? 0 : 1;
			file.write(prefix, 3);
		}

		CdBufSpace buf(MCWorkingGeno.Space, !(*SNPOrder), CdBufSpace::acInc);
		long nRe = buf.BufElmSize() % 4;
		long nPack = (nRe > 0) ? (buf.BufElmSize()/4 + 1) : (buf.BufElmSize()/4);
		vector<char> geno(nPack, 0);
		const UInt8 cvt[4] = { 3, 2, 0, 1 };

		for (long i=0; i < buf.IdxCnt(); i++)
		{
			UInt8 *s = buf.ReadGeno(i);
			char *p = &geno[0];
			for (long k=0; k < buf.BufElmSize()/4; k++, s+=4)
			{
				*p++ = cvt[s[0] & 0x03] | (cvt[s[1] & 0x03] << 2) |
					(cvt[s[2] & 0x03] << 4) | (cvt[s[3] & 0x03] << 6);
			}
			if (nRe > 0)
			{
				UInt8 b = 0;
				for (long k=0; k < nRe; k++, s++) b |= (cvt[*s & 0x03] << (2*k));
				*p++ = b;
			}
			file.write(&geno[0], nPack);
			MCWorkingGeno.Progress.Forward(1);
		}
		*out_err = 0;
	CORECATCH(*out_err = 1)
}

/// to convert from GDS to EIGENSOFT
DLLEXPORT void gnrConvGDS2EIGEN(char **pedfn, LongBool *verbose, LongBool *out_err)
{
	CORETRY
		MCWorkingGeno.Progress.Info = "\tOutput: ";
		MCWorkingGeno.Progress.Show() = *verbose;
		MCWorkingGeno.Progress.Init(MCWorkingGeno.Space.SNPNum());

		ofstream file(*pedfn);
		if (!file.good())
			throw ErrCoreArray("Fail to create the file %s.", *pedfn);
		CdBufSpace buf(MCWorkingGeno.Space, true, CdBufSpace::acInc);
		for (long i=0; i < buf.IdxCnt(); i++)
		{
			UInt8 *g = buf.ReadGeno(i);
			for (long j=0; j < MCWorkingGeno.Space.SampleNum(); j++, g++)
			{
				int geno = (*g <= 2) ? (*g) : 9;
				file << geno;
			}
			file << endl;
			MCWorkingGeno.Progress.Forward(1);
		}
		*out_err = 0;
	CORECATCH(*out_err = 1)
}

/// to detect PLINK BED
DLLEXPORT void gnrBEDFlag(char **bedfn, int *SNPOrder, LongBool *out_err)
{
	CORETRY
		ifstream file(*bedfn, ios::binary);
		if (!file.good())
			throw ErrCoreArray("Cannot open the file %s.", *bedfn);
		char prefix[3];
		file.read(prefix, 3);
		if ((prefix[0] != 0x6C) || (prefix[1] != 0x1B))
			throw ErrCoreArray("Invalid prefix in the bed file.");
		*SNPOrder = (UInt8)(prefix[2]);
		*out_err = 0;
	CORECATCH(*out_err = 1)
}

/// to convert from PLINK BED to GDS
DLLEXPORT void gnrConvBED2GDS(char **bedfn, TdSequenceX *Node, LongBool *verbose,
	LongBool *out_err)
{
	CORETRY
		TdSequenceX Seq = *Node;
		int DLen[2];
		gds_SeqGetDim(Seq, DLen);

		MCWorkingGeno.Progress.Info = " ";
		MCWorkingGeno.Progress.Show() = *verbose;
		MCWorkingGeno.Progress.Init(DLen[0]);

		ifstream file(*bedfn, ios::binary);
		if (!file.good())
			throw ErrCoreArray("Fail to open the file %s.", *bedfn);
		// read prefix
		char prefix[3];
		file.read(prefix, 3);

		long nRe = DLen[1] % 4;
		long nPack = (nRe > 0) ? (DLen[1]/4 + 1) : (DLen[1]/4);
		auto_ptr<char> srcgeno(new char[nPack]);
		auto_ptr<UInt8> dstgeno(new UInt8[DLen[1]]);
		CoreArray::Int32 st[2] = { 0, 0 }, cnt[2] = { 1, DLen[1] };
		const UInt8 cvt[4] = { 2, 3, 1, 0 };

		for (long i=0; i < DLen[0]; i++)
		{
			// read genotypes
			file.read(srcgeno.get(), nPack);
			// unpacked
			UInt8 *p = dstgeno.get();
			for (long k=0; k < DLen[1]/4; k++)
			{
				UInt8 g = srcgeno.get()[k];
				*p++ = cvt[g & 0x03]; g >>= 2;
				*p++ = cvt[g & 0x03]; g >>= 2;
				*p++ = cvt[g & 0x03]; g >>= 2;
				*p++ = cvt[g & 0x03];
			}
			if (nRe > 0)
			{
				UInt8 g = srcgeno.get()[DLen[1]/4];
				for (long k=0; k < nRe; k++)
				{
					*p++ = cvt[g & 0x03]; g >>= 2;
				}
			}
			// write
			st[0] = i;
			gds_wData(Seq, st, cnt, dstgeno.get(), svUInt8);
			MCWorkingGeno.Progress.Forward(1);
		}
		*out_err = 0;
	CORECATCH(*out_err = 1)
}


/// SNP functions

static inline bool ATGC(char ch)
	{ return (ch=='A') || (ch=='T') || (ch=='G') || (ch=='C'); }
static inline int ALLELE_MINOR(double freq)
	{ return (freq <= 0.5) ? 0 : 1; }

/// to detect strand problem
DLLEXPORT void gnrAlleleStrand(char *allele1[], double afreq1[], int I1[],
	char *allele2[], double afreq2[], int I2[], int *n,
	LongBool out_flag[], LongBool *out_err)
{
	CORETRY
		// init
		map<char, char> MAP;
		MAP['A'] = 'T'; MAP['C'] = 'G'; MAP['G'] = 'C'; MAP['T'] = 'A';

		// loop for each SNP
		for (int i=0; i < *n; i++)
		{
			bool switch_flag = false; // if true, need switch strand

			// ref / nonref alleles
			char s1=allele1[I1[i]-1][0], s2=allele1[I1[i]-1][2];
			char p1=allele2[I2[i]-1][0], p2=allele2[I2[i]-1][2];
			double F1=afreq1[I1[i]-1], F2=afreq2[I2[i]-1];

			if (ATGC(s1) && ATGC(s2) && ATGC(p1) && ATGC(p2))
			{
				// check
				if ( (s1 == p1) && (s2 == p2) )
				{
					// for example, + C/G <---> - C/G, strand ambi
					if (s1 == MAP[p2])
						switch_flag = (ALLELE_MINOR(F1) != ALLELE_MINOR(F2));
				} else if ( (s1 == p2) && (s2 == p1) )
				{
					// for example, + C/G <---> - G/C
					if (s1 == MAP[p1])
						switch_flag = (ALLELE_MINOR(F1) != ALLELE_MINOR(F2));
					else
						switch_flag = true;
				} else if ( (s1 == MAP[p1]) && (s2 == MAP[p2]) )
				{
					// for example, + C/G <---> - G/C
					if (s1 == p2)
						switch_flag = (ALLELE_MINOR(F1) != ALLELE_MINOR(F2));
				} else if ( (s1 == MAP[p2]) && (s2 == MAP[p1]) )
				{
					switch_flag = true;
				} else {
					// throw ErrHLA("Invalid strand in SNP %s", s.rsid.c_str());
				}
			}
			out_flag[i] = switch_flag;
		}
		*out_err = 0;
	CORECATCH(*out_err = 1)
}



/// get an error message
DLLEXPORT void gnrErrMsg(char **Msg)
{
	RStrAgn(gds_LastError().c_str(), Msg);
}


} // extern "C"
