// ===========================================================
//     _/_/_/   _/_/_/  _/_/_/_/    _/_/_/_/  _/_/_/   _/_/_/
//      _/    _/       _/             _/    _/    _/   _/   _/
//     _/    _/       _/_/_/_/       _/    _/    _/   _/_/_/
//    _/    _/       _/             _/    _/    _/   _/
// _/_/_/   _/_/_/  _/_/_/_/_/     _/     _/_/_/   _/_/
// ===========================================================
//
// SNPRelate.cpp: Relatedness, Linkage Disequilibrium and
//                Principal Component Analysis
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


#include <dType.h>
#include <fstream>
#include <map>
#include <dGenGWAS.h>
#include <dVect.h>
#include <R.h>
#include <Rdefines.h>
#include <R_ext/Lapack.h>


// R_XLEN_T_MAX is defined, R >= v3.0
#ifndef R_XLEN_T_MAX
#  define R_xlen_t	R_len_t
#  define XLENGTH	Rf_length
#endif


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
	struct TIBS_Flag
	{
		UInt32 IBS0, IBS1, IBS2;
		TIBS_Flag() { IBS0 = IBS1 = IBS2 = 0; }
	};

	void AutoDetectSNPBlockSize(int nSamp, bool Detect=true);
	void DoPLINKIBSCalculate(CdMatTriDiag<TIBS_Flag> &PublicIBS, int NumThread,
		const char *Info, bool verbose);
	void DoIBSCalculate(CdMatTri<TIBS_Flag> &PublicIBS, int NumThread,
		const char *Info, bool verbose);

	/// The structure of genetic distance
	struct TDissflag
	{
		Int64 SumGeno;
		double SumAFreq;
		TDissflag() { SumGeno = 0; SumAFreq = 0; }
	};

	/// Calculate the genetic distance matrix
	void DoDissCalculate(CdMatTri<TDissflag> &PublicDist, int NumThread,
		const char *Info, bool verbose);


	/// The structure of KING IBD estimator
	struct TKINGHomoFlag
	{
		UInt32 IBS0;       //< the number of loci sharing no allele
		UInt32 SumSq;      //< \sum_m (X_m^{(i)} - X_m^{(j)})^2
		double SumAFreq;   //< \sum_m p_m (1 - p_m)
		double SumAFreq2;  //< \sum_m p_m^2 (1 - p_m)^2
		TKINGHomoFlag() { IBS0 = SumSq = 0; SumAFreq = SumAFreq2 = 0; }
	};

	struct TKINGRobustFlag
	{
		UInt32 IBS0;       //< the number of loci sharing no allele
		UInt32 nLoci;      //< the total number of loci
		UInt32 SumSq;      //< \sum_m (X_m^{(i)} - X_m^{(j)})^2
		UInt32 N1_Aa;      //< the number of hetet loci for the first individual
		UInt32 N2_Aa;      //< the number of hetet loci for the second individual
		TKINGRobustFlag() { IBS0 = nLoci = SumSq = N1_Aa = N2_Aa = 0; }
	};

	/// Calculate KING IBD estimators
	void DoKINGCalculate(CdMatTri<TKINGHomoFlag> &PublicKING, int NumThread,
		const char *Info, bool verbose);
	/// Calculate KING IBD estimators
	void DoKINGCalculate(CdMatTri<TKINGRobustFlag> &PublicKING, int NumThread,
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
	void Init_EPrIBD_IBS(const double in_afreq[], double out_afreq[],
		bool CorrectFactor, long nSNP = -1);
	void Est_PLINK_Kinship(int IBS0, int IBS1, int IBS2, double &k0, double &k1,
		bool KinshipConstrict);

	// MLE
	extern double EPrIBS_IBD[3][3];
	extern long nIterMax;
	extern double FuncRelTol;
	extern int MethodMLE;
	extern bool Loglik_Adjust;
	extern double *MLEAlleleFreq;
	extern long nTotalSNP;


	void InitPackedGeno(void *buffer);
	void prIBDTable(int g1, int g2, double &t0, double &t1, double &t2, double p);

	void Do_MLE_IBD_Calculate(const double *AFreq, CdMatTriDiag<TIBD> &PublicIBD,
		CdMatTriDiag<int> *PublicNIter,
		double *out_AFreq, int NumThread, const char *Info, double *tmpAF, bool verbose);
	void Do_MLE_LogLik(const double *AFreq, const double *k0, const double *k1,
		double *tmp_AF, double *out_loglik);
	void Do_MLE_LogLik_k01(const double *AFreq, const double k0, const double k1,
		double *tmp_AF, double *out_loglik);

	void Do_MLE_IBD_Pair(int n, const int *geno1, const int *geno2, const double *AFreq,
		double &out_k0, double &out_k1, double &out_loglik, int &out_niter, double tmpprob[]);
}

namespace LD
{
	extern int LD_Method;

	void InitPackedGeno();
	void DonePackedGeno();
	double calcLD(const int *snp1, const int *snp2, int cnt,
		double &pA_A, double &pA_B, double &pB_A, double &pB_B);
	void calcLD_mat(int nThread, double *out_LD);
	void calcLD_slide_mat(int nThread, double *out_LD, int n_slide);
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
	// calling .C()
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

	// calling .Call()
	#define CORETRY_CALL	CORETRY
	#define CORECATCH_CALL	CORECATCH(error(gds_LastError().c_str()))


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
DLLEXPORT void gnrSetGenoSpace(PdSequenceX *Node,
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
		if (*out_nsnp <= 0)
			throw ErrCoreArray("There is no SNP!");
		*out_nsamp = MCWorkingGeno.Space.SampleNum();
		if (*out_nsamp <= 0)
			throw ErrCoreArray("There is no sample!");

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
		const R_xlen_t n = MCWorkingGeno.Space.SNPNum();
		auto_ptr<CBOOL> sel(new CBOOL[n]);
		*out_num = MCWorkingGeno.Space.Select_SNP_Base((*remove_mono) != 0,
			*maf, *missrate, sel.get());
		for (int i=0; i < n; i++) out_selection[i] = sel.get()[i];
		*out_err = 0;
	CORECATCH(*out_err = 1)
}

/** to compute missing rates, and specify allele frequency
 *  \param afreq          specifying the allele frequencies
 *  \param remove_mono    whether remove monomorphic snps or not
 *  \param maf            the threshold of minor allele frequencies, keeping ">= maf"
 *  \param missrate       the threshold of missing rates, keeping "<= missing.rate"
**/
DLLEXPORT void gnrSelSNP_Base_Ex(double *afreq, LongBool *remove_mono, double *maf,
	double *missrate, int *out_num, LongBool *out_selection, LongBool *out_err)
{
	CORETRY
		const R_xlen_t n = MCWorkingGeno.Space.SNPNum();
		auto_ptr<CBOOL> sel(new CBOOL[n]);
		*out_num = MCWorkingGeno.Space.Select_SNP_Base_Ex(afreq, (*remove_mono) != 0,
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
DLLEXPORT void gnrCopyGeno(PdSequenceX *Node, LongBool *snpfirstorder, LongBool *out_err)
{
	CORETRY
		*out_err = 1;
		if (*snpfirstorder)
		{
			CoreArray::Int32 cnt[2] = { 1, MCWorkingGeno.Space.SNPNum() };
			CdBufSpace buf(MCWorkingGeno.Space, false, CdBufSpace::acInc);
			for (int i=0; i < buf.IdxCnt(); i++)
			{
				UInt8 *p = buf.ReadGeno(i);
				CoreArray::Int32 st[2] = { (CoreArray::Int32)i, 0 };
				if (!gds_wData(*Node, st, cnt, p, svUInt8))
					return;
			}
		} else {
			CoreArray::Int32 cnt[2] = { 1, MCWorkingGeno.Space.SampleNum() };
			CdBufSpace buf(MCWorkingGeno.Space, true, CdBufSpace::acInc);
			for (int i=0; i < buf.IdxCnt(); i++)
			{
				UInt8 *p = buf.ReadGeno(i);
				CoreArray::Int32 st[2] = { (CoreArray::Int32)i, 0 };
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
DLLEXPORT void gnrAppendGenoSpace(PdSequenceX *Node, LongBool *snpfirstorder, LongBool *out_err)
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
DLLEXPORT void gnrAppendGenoSpaceStrand(PdSequenceX *Node, LongBool *snpfirstorder,
	LongBool StrandFlag[], LongBool *out_err)
{
	CORETRY
		*out_err = 1;
		if (*snpfirstorder)
		{
			const R_xlen_t n = MCWorkingGeno.Space.SNPNum();
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
			const R_xlen_t n = MCWorkingGeno.Space.SampleNum();
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
				Rprintf("PCA:\tthe sum of all working genotypes (0, 1 and 2) = %.0f\n", GenoSum);
		}

		// ******** The calculation of genetic covariance matrix ********

		// the number of samples
		const R_xlen_t n = MCWorkingGeno.Space.SampleNum();
		// set parameters
		PCA::BayesianNormal = ((*_BayesianNormal) == TRUE);
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
				int _n = n;
				F77_NAME(dspev)("V", "L", &_n, Cov.get(), out_Eigenvalues,
					tmp_EigenVec.get(), &_n, tmp_Work.get(), &info);
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
				Rprintf("SNP Correlations:\tthe sum of all working genotypes (0, 1 and 2) = %.0f\n", GenoSum);
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
				Rprintf("SNP Loadings:\tthe sum of all working genotypes (0, 1 and 2) = %.0f\n", GenoSum);
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
				Rprintf("Sample Loadings:\tthe sum of all working genotypes (0, 1 and 2) = %.0f\n", GenoSum);
		}

		// ******** To compute the snp correlation ********
		PCA::DoSampLoadingCalculate(AveFreq, Scale, *EigenCnt, SNPLoadings,
			EigenVal, *Num, *TraceXTX, out_samploading,
			*NumThread, "Sample Loadings:", *Verbose);

		// output
		*out_err = 0;
	CORECATCH(*out_err = 1)
}


// ***********************************************************************
// the functions for identity-by-state (IBS)
//

// the functions for Identity by state (IBS)

/// to compute the average IBS
DLLEXPORT SEXP gnrIBSAve(SEXP Verbose, SEXP DataCache, SEXP NumThread)
{
	CORETRY_CALL

		// ******** To cache the genotype data ********
		if (LOGICAL(DataCache)[0] == TRUE)
		{
			double GenoSum=0;
			gnrCacheGeno(&GenoSum, NULL);
			if (LOGICAL(Verbose)[0] == TRUE)
				Rprintf("IBS:\tthe sum of all working genotypes (0, 1 and 2) = %.0f\n", GenoSum);
		}

		// ******** The calculation of IBS matrix ********

		// the number of samples
		const R_xlen_t n = MCWorkingGeno.Space.SampleNum();
		// to detect the block size
		IBS::AutoDetectSNPBlockSize(n);
		// the upper-triangle genetic covariance matrix
		CdMatTri<IBS::TIBS_Flag> IBS(n);

		// initialize output variables
		SEXP dim;
		SEXP ibsmat = NEW_NUMERIC(n*n);  PROTECT(ibsmat);
		PROTECT(dim = NEW_INTEGER(2));
		INTEGER(dim)[0] = n; INTEGER(dim)[1] = n;
		setAttrib(ibsmat, R_DimSymbol, dim);
		double *out_IBSMat = REAL(ibsmat);

		// Calculate the IBS matrix
		IBS::DoIBSCalculate(IBS, INTEGER(NumThread)[0], "IBS:",
			LOGICAL(Verbose)[0] == TRUE);

		// output
		IBS::TIBS_Flag *p = IBS.get();
		for (int i=0; i < n; i++)
		{
			for (int j=i; j < n; j++, p++)
			{
				out_IBSMat[i*n + j] = out_IBSMat[j*n + i] =
					double(0.5*p->IBS1 + p->IBS2) / (p->IBS0 + p->IBS1 + p->IBS2);
			}
		}

		UNPROTECT(2);
		return ibsmat;

	CORECATCH_CALL
}

/// to compute the average IBS
DLLEXPORT SEXP gnrIBSNum(SEXP Verbose, SEXP DataCache, SEXP NumThread)
{
	CORETRY_CALL

		// ******** To cache the genotype data ********
		if (LOGICAL(DataCache)[0] == TRUE)
		{
			double GenoSum=0;
			gnrCacheGeno(&GenoSum, NULL);
			if (LOGICAL(Verbose)[0] == TRUE)
				Rprintf("IBS:\tthe sum of all working genotypes (0, 1 and 2) = %.0f\n", GenoSum);
		}

		// ******** The calculation of IBS matrix ********

		// the number of samples
		const R_xlen_t n = MCWorkingGeno.Space.SampleNum();
		// to detect the block size
		IBS::AutoDetectSNPBlockSize(n);

		// the upper-triangle IBS matrix
		CdMatTri<IBS::TIBS_Flag> IBS(n);

		// Calculate the IBS matrix
		IBS::DoIBSCalculate(IBS, INTEGER(NumThread)[0], "IBS:",
			LOGICAL(Verbose)[0] == TRUE);

		// initialize output variables
		SEXP dim;
		PROTECT(dim = NEW_INTEGER(2));
		INTEGER(dim)[0] = n; INTEGER(dim)[1] = n;

		SEXP IBS0    = NEW_NUMERIC(n*n);  PROTECT(IBS0);
		setAttrib(IBS0, R_DimSymbol, dim);
		SEXP IBS1    = NEW_NUMERIC(n*n);  PROTECT(IBS1);
		setAttrib(IBS1, R_DimSymbol, dim);
		SEXP IBS2    = NEW_NUMERIC(n*n);  PROTECT(IBS2);
		setAttrib(IBS2, R_DimSymbol, dim);

		double *out_IBS0    = REAL(IBS0);
		double *out_IBS1    = REAL(IBS1);
		double *out_IBS2    = REAL(IBS2);

		// output
		IBS::TIBS_Flag *p = IBS.get();
		for (int i=0; i < n; i++)
		{
			for (int j=i; j < n; j++, p++)
			{
				out_IBS0[i*n + j] = out_IBS0[j*n + i] = p->IBS0;
				out_IBS1[i*n + j] = out_IBS1[j*n + i] = p->IBS1;
				out_IBS2[i*n + j] = out_IBS2[j*n + i] = p->IBS2;
			}
		}

		// output
		SEXP ans;
		PROTECT(ans = NEW_LIST(3));
		SET_ELEMENT(ans, 0, IBS0);
		SET_ELEMENT(ans, 1, IBS1);
		SET_ELEMENT(ans, 2, IBS2);
		UNPROTECT(5);

		return ans;

	CORECATCH_CALL
}



// ***********************************************************************
// the functions for identity-by-descent (IBD)
//

/// to compute the IBD coefficients by PLINK method of moment
DLLEXPORT SEXP gnrIBD_PLINK(SEXP Verbose, SEXP DataCache, SEXP NumThread,
	SEXP AlleleFreq, SEXP UseSpecificAFreq, SEXP KinshipConstrict)
{
	CORETRY_CALL

		// ******** to cache the genotype data ********
		if (LOGICAL(DataCache)[0] == TRUE)
		{
			double GenoSum=0;
			gnrCacheGeno(&GenoSum, NULL);
			if (LOGICAL(Verbose)[0] == TRUE)
				Rprintf("PLINK IBD:\tthe sum of all working genotypes (0, 1 and 2) = %.0f\n", GenoSum);
		}

		// the number of individuals
		const R_xlen_t n = MCWorkingGeno.Space.SampleNum();
		// the number of SNPs
		const R_xlen_t m = MCWorkingGeno.Space.SNPNum();
		// to detect the block size
		IBS::AutoDetectSNPBlockSize(n);

		// ******** PLINK method of moment ********

		// initialize output variables
		SEXP dim;
		PROTECT(dim = NEW_INTEGER(2));
		INTEGER(dim)[0] = n; INTEGER(dim)[1] = n;

		SEXP k0    = NEW_NUMERIC(n*n);  PROTECT(k0);
		setAttrib(k0, R_DimSymbol, dim);

		SEXP k1    = NEW_NUMERIC(n*n);  PROTECT(k1);
		setAttrib(k1, R_DimSymbol, dim);

		SEXP afreq = NEW_NUMERIC(m);    PROTECT(afreq);
		double *out_k0    = REAL(k0);
		double *out_k1    = REAL(k1);
		double *out_afreq = REAL(afreq);

		// initialize the internal matrix
		IBD::Init_EPrIBD_IBS(LOGICAL(UseSpecificAFreq)[0] ? REAL(AlleleFreq) : NULL,
			out_afreq, !(LOGICAL(UseSpecificAFreq)[0]));
		// the upper-triangle genetic covariance matrix
		CdMatTriDiag<IBS::TIBS_Flag> IBS(IBS::TIBS_Flag(), n);
		// Calculate the IBS matrix
		IBS::DoPLINKIBSCalculate(IBS, INTEGER(NumThread)[0],
			"PLINK IBD:", LOGICAL(Verbose)[0]==TRUE);

		// output
		bool kc = LOGICAL(KinshipConstrict)[0] == TRUE;
		IBS::TIBS_Flag *p = IBS.get();
		for (int i=0; i < n; i++)
		{
			out_k0[i*n + i] = out_k1[i*n + i] = 0;
			for (int j=i+1; j < n; j++, p++)
			{
				double k0, k1;
				IBD::Est_PLINK_Kinship(p->IBS0, p->IBS1, p->IBS2, k0, k1, kc);
				out_k0[i*n + j] = out_k0[j*n + i] = k0;
				out_k1[i*n + j] = out_k1[j*n + i] = k1;
			}
		}

		// output
		SEXP ans;
		PROTECT(ans = NEW_LIST(3));
		SET_ELEMENT(ans, 0, k0);
		SET_ELEMENT(ans, 1, k1);
		SET_ELEMENT(ans, 2, afreq);
		UNPROTECT(5);

		return ans;

	CORECATCH_CALL
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
				Rprintf("MLE IBD:\tthe sum of all working genotypes (0, 1 and 2) = %.0f\n", GenoSum);
		}

		// ******** MLE IBD ********

		// initialize the packed genotypes
		IBD::InitPackedGeno(tmp_buffer);
		// initialize the internal matrix
		IBD::Init_EPrIBD_IBS(*UseSpecificAFreq ? AlleleFreq : NULL, NULL, false);

		IBD::nIterMax = *MaxIterCnt; IBD::FuncRelTol = *RelTol;
		IBD::MethodMLE = *method; IBD::Loglik_Adjust = *CoeffCorrect;
		IBD::KinshipConstraint = (*KinshipConstraint != 0);

		// the upper-triangle genetic covariance matrix
		const R_xlen_t n = MCWorkingGeno.Space.SampleNum();
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


/// to compute the IBD coefficients by MLE for a pair of individuals
DLLEXPORT void gnrPairIBD(int *n, int *geno1, int *geno2, double *AlleleFreq,
	LongBool *KinshipConstraint, int *MaxIterCnt, double *RelTol,
	LongBool *CoeffCorrect, int *method,
	double *out_k0, double *out_k1, double *out_loglik, int *out_niter,
	double *tmp_buffer, LongBool *out_err)
{
	CORETRY
		// initialize
		IBD::nIterMax = *MaxIterCnt; IBD::FuncRelTol = *RelTol;
		IBD::MethodMLE = *method; IBD::Loglik_Adjust = *CoeffCorrect;
		IBD::KinshipConstraint = (*KinshipConstraint != 0);
		IBD::Init_EPrIBD_IBS(AlleleFreq, NULL, false, *n);

		// get the initial values
		int IBS[3] = { 0, 0, 0 };
		for (int i=0; i < *n; i++)
		{
			if ((0<=geno1[i]) && (geno1[i]<=2) && (0<=geno2[i]) && (geno2[i]<=2))
			{
				IBS[2 - abs(geno1[i] - geno2[i])] ++;
			}
		}
		IBD::Est_PLINK_Kinship(IBS[0], IBS[1], IBS[2], *out_k0, *out_k1,
			IBD::KinshipConstraint);

		// compute
		if (*method >= 0)
		{
			IBD::Do_MLE_IBD_Pair(*n, geno1, geno2, AlleleFreq,
				*out_k0, *out_k1, *out_loglik, *out_niter, tmp_buffer);
		} else {
			*out_loglik = conf_F64_NaN();
			*out_niter = 0;
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


/// to compute the value of log likelihood for a pair of individuals
DLLEXPORT void gnrPairIBDLogLik(int *n, int *geno1, int *geno2, double *AlleleFreq,
	double *k0, double *k1, double *out_loglik, double *tmp_buffer, LongBool *out_err)
{
	CORETRY
		// initialize the probability table
		double *PrIBD = tmp_buffer;
		for (int i=0; i < *n; i++)
		{
			IBD::prIBDTable(geno1[i], geno2[i], PrIBD[0], PrIBD[1], PrIBD[2], AlleleFreq[i]);
			PrIBD += 3;
		}

		// calculate log likelihood value
		double k[3] = { *k0, *k1, 1 - *k0 - *k1 }, LogLik = 0;
		PrIBD = tmp_buffer;
		for (int i=0; i < *n; i++)
		{
			double sum = PrIBD[0]*k[0] + PrIBD[1]*k[1] + PrIBD[2]*k[2];
			if (sum > 0)
				LogLik += log(sum);
			PrIBD += 3;
		}

		// output
		*out_loglik = LogLik;
		*out_err = 0;
	CORECATCH(*out_err = 1)
}


/// to compute the IBD coefficients by KING method of moment (KING-homo)
DLLEXPORT SEXP gnrIBD_KING_Homo(SEXP Verbose, SEXP DataCache, SEXP NumThread)
{
	CORETRY_CALL

		// ******** To cache the genotype data ********
		if (LOGICAL(DataCache)[0] == TRUE)
		{
			double GenoSum=0;
			gnrCacheGeno(&GenoSum, NULL);
			if (LOGICAL(Verbose)[0] == TRUE)
				Rprintf("KING IBD:\tthe sum of all working genotypes (0, 1 and 2) = %.0f\n", GenoSum);
		}

		// ******** The calculation of genetic covariance matrix ********

		// the number of samples
		const R_xlen_t n = MCWorkingGeno.Space.SampleNum();
		// to detect the block size
		IBS::AutoDetectSNPBlockSize(n);

		// the upper-triangle IBD matrix
		CdMatTri<IBS::TKINGHomoFlag> IBD(n);
		// Calculate the IBD matrix
		IBS::DoKINGCalculate(IBD, INTEGER(NumThread)[0], "KING IBD:",
			LOGICAL(Verbose)[0] == TRUE);

		// initialize output variables
		SEXP dim;
		PROTECT(dim = NEW_INTEGER(2));
		INTEGER(dim)[0] = n; INTEGER(dim)[1] = n;

		SEXP k0    = NEW_NUMERIC(n*n);  PROTECT(k0);
		setAttrib(k0, R_DimSymbol, dim);
		double *out_k0 = REAL(k0);

		SEXP k1    = NEW_NUMERIC(n*n);  PROTECT(k1);
		setAttrib(k1, R_DimSymbol, dim);
		double *out_k1 = REAL(k1);

		// output
		IBS::TKINGHomoFlag *p = IBD.get();
		for (int i=0; i < n; i++)
		{
			out_k0[i*n + i] = out_k1[i*n + i] = 0;
			p ++;
			for (int j=i+1; j < n; j++, p++)
			{
				double theta = 0.5 - p->SumSq / (8 * p->SumAFreq);
				double k0 = p->IBS0 / (2 * p->SumAFreq2);
				double k1 = 2 - 2*k0 - 4*theta;
				out_k0[i*n + j] = out_k0[j*n + i] = k0;
				out_k1[i*n + j] = out_k1[j*n + i] = k1;
			}
		}

		// output
		SEXP ans;
		PROTECT(ans = NEW_LIST(2));
		SET_ELEMENT(ans, 0, k0);
		SET_ELEMENT(ans, 1, k1);
		UNPROTECT(4);

		return ans;

	CORECATCH_CALL
}

/// to compute the IBD coefficients by KING method of moment (KING-robust)
DLLEXPORT SEXP gnrIBD_KING_Robust(SEXP Verbose, SEXP DataCache,
	SEXP NumThread, SEXP FamilyID)
{
	CORETRY_CALL

		// ******** To cache the genotype data ********
		if (LOGICAL(DataCache)[0] == TRUE)
		{
			double GenoSum=0;
			gnrCacheGeno(&GenoSum, NULL);
			if (LOGICAL(Verbose)[0] == TRUE)
				Rprintf("KING IBD:\tthe sum of all working genotypes (0, 1 and 2) = %.0f\n", GenoSum);
		}

		// ******** The calculation of genetic covariance matrix ********

		// the number of samples
		const R_xlen_t n = MCWorkingGeno.Space.SampleNum();
		// to detect the block size
		IBS::AutoDetectSNPBlockSize(n);

		// IBD KING robust IBD estimator across family
		// the upper-triangle IBD matrix
		CdMatTri<IBS::TKINGRobustFlag> IBD(n);
		// Calculate the IBD matrix
		IBS::DoKINGCalculate(IBD, INTEGER(NumThread)[0], "KING IBD:",
			LOGICAL(Verbose)[0] == TRUE);


		// initialize output variables
		SEXP dim;
		PROTECT(dim = NEW_INTEGER(2));
		INTEGER(dim)[0] = n; INTEGER(dim)[1] = n;

		SEXP IBS0    = NEW_NUMERIC(n*n);  PROTECT(IBS0);
		setAttrib(IBS0, R_DimSymbol, dim);
		double *out_IBS0 = REAL(IBS0);

		SEXP Kinship = NEW_NUMERIC(n*n);  PROTECT(Kinship);
		setAttrib(Kinship, R_DimSymbol, dim);
		double *out_kinship = REAL(Kinship);

		int *FamID_Ptr = INTEGER(FamilyID);

		// output
		IBS::TKINGRobustFlag *p = IBD.get();
		for (int i=0; i < n; i++)
		{
			out_IBS0[i*n + i] = 0; out_kinship[i*n + i] = 0.5;
			p ++;
			for (int j=i+1; j < n; j++, p++)
			{
				out_IBS0[i*n + j] = out_IBS0[j*n + i] = double(p->IBS0) / p->nLoci;

				int f1 = FamID_Ptr[i];
				int f2 = FamID_Ptr[j];
				out_kinship[i*n + j] = out_kinship[j*n + i] =
						((f1 == f2) && (f1 != NA_INTEGER)) ?
					(0.5 - p->SumSq / (2.0 *(p->N1_Aa + p->N2_Aa))) :
					(0.5 - p->SumSq / (4.0 * min(p->N1_Aa, p->N2_Aa)));
			}
		}

		// output
		SEXP ans;
		PROTECT(ans = NEW_LIST(2));
		SET_ELEMENT(ans, 0, IBS0);
		SET_ELEMENT(ans, 1, Kinship);
		UNPROTECT(4);

		return ans;

	CORECATCH_CALL
}


/// to create a list of sample id
//  ID1 = matrix(sample.id, nrow=n, ncol=n, byrow=TRUE)
DLLEXPORT SEXP gnrIBDSelSampID_List1(SEXP SampID, SEXP Flag)
{
	// the number of samples
	const R_xlen_t n  = XLENGTH(SampID);
	const R_xlen_t nF = XLENGTH(Flag);
	R_xlen_t flag_cnt = 0;
	int *p_flag = LOGICAL(Flag);
	for (R_xlen_t i=0; i < nF; i++)
	{
		if (*p_flag == TRUE) flag_cnt ++;
		p_flag ++;
	}

	if (isFactor(SampID))
		SampID = Rf_asCharacterFactor(SampID);

	// output
	p_flag = LOGICAL(Flag);
	R_xlen_t idx = 0;
	SEXP ans;

	if (IS_CHARACTER(SampID))
	{
		PROTECT(ans = NEW_STRING(flag_cnt));
		for (R_xlen_t i=0; i < n; i++)
		{
			for (R_xlen_t j=0; j < n; j++, p_flag++)
			{
				if (*p_flag == TRUE)
					SET_STRING_ELT(ans, idx++, STRING_ELT(SampID, i));
			}
		}
	} else if (IS_NUMERIC(SampID))
	{
		PROTECT(ans = NEW_NUMERIC(flag_cnt));
		for (R_xlen_t i=0; i < n; i++)
		{
			for (R_xlen_t j=0; j < n; j++, p_flag++)
			{
				if (*p_flag == TRUE)
					REAL(ans)[idx++] = REAL(SampID)[i];
			}
		}
	} else if (IS_INTEGER(SampID))
	{
		PROTECT(ans = NEW_INTEGER(flag_cnt));
		for (R_xlen_t i=0; i < n; i++)
		{
			for (R_xlen_t j=0; j < n; j++, p_flag++)
			{
				if (*p_flag == TRUE)
					INTEGER(ans)[idx++] = INTEGER(SampID)[i];
			}
		}
	} else
		error("'sample.id' should be numeric- or character- type.");

	UNPROTECT(1);
	return ans;
}


/// to create a list of sample id
//  ID2 = matrix(sample.id, nrow=n, ncol=n)[flag]
DLLEXPORT SEXP gnrIBDSelSampID_List2(SEXP SampID, SEXP Flag)
{
	// the number of samples
	const R_xlen_t n  = XLENGTH(SampID);
	const R_xlen_t nF = XLENGTH(Flag);
	R_xlen_t flag_cnt = 0;
	int *p_flag = LOGICAL(Flag);
	for (R_xlen_t i=0; i < nF; i++)
	{
		if (*p_flag == TRUE) flag_cnt ++;
		p_flag ++;
	}

	if (isFactor(SampID))
		SampID = Rf_asCharacterFactor(SampID);

	// output
	p_flag = LOGICAL(Flag);
	R_xlen_t idx = 0;
	SEXP ans;

	if (IS_CHARACTER(SampID))
	{
		PROTECT(ans = NEW_STRING(flag_cnt));
		for (R_xlen_t i=0; i < n; i++)
		{
			for (R_xlen_t j=0; j < n; j++, p_flag++)
			{
				if (*p_flag == TRUE)
					SET_STRING_ELT(ans, idx++, STRING_ELT(SampID, j));
			}
		}
	} else if (IS_NUMERIC(SampID))
	{
		PROTECT(ans = NEW_NUMERIC(flag_cnt));
		for (R_xlen_t i=0; i < n; i++)
		{
			for (R_xlen_t j=0; j < n; j++, p_flag++)
			{
				if (*p_flag == TRUE)
					REAL(ans)[idx++] = REAL(SampID)[j];
			}
		}
	} else if (IS_INTEGER(SampID))
	{
		PROTECT(ans = NEW_INTEGER(flag_cnt));
		for (R_xlen_t i=0; i < n; i++)
		{
			for (R_xlen_t j=0; j < n; j++, p_flag++)
			{
				if (*p_flag == TRUE)
					INTEGER(ans)[idx++] = INTEGER(SampID)[j];
			}
		}
	} else
		error("'sample.id' should be numeric- or character- type.");

	UNPROTECT(1);
	return ans;
}



// ***********************************************************************
// the functions for linkage disequilibrium (LD)
//

/// the functions for Linkage Disequilibrium (LD) analysis
DLLEXPORT void gnrLDpair(int *snp1, int *snp2, int *len, int *method,
	double *out_LD, double &pA_A, double &pA_B, double &pB_A, double &pB_B,
	LongBool *out_err)
{
	CORETRY
		LD::LD_Method = *method;
		*out_LD = LD::calcLD(snp1, snp2, *len, pA_A, pA_B, pB_A, pB_B);
		// output
		*out_err = 0;
	CORECATCH(*out_err = 1)
}


/// to compute the IBD coefficients by MLE
DLLEXPORT void gnrLDMat(int *method, LongBool *Verbose, LongBool *DataCache,
	int *NumThread, int *n_slide, double *out_LD, LongBool *out_err)
{
	CORETRY
		// ******** to cache the genotype data ********
		if (*DataCache)
		{
			double GenoSum=0;
			gnrCacheGeno(&GenoSum, NULL);
			if (*Verbose)
				Rprintf("LD matrix:\tthe sum of all working genotypes (0, 1 and 2) = %.0f\n", GenoSum);
		}

		// initialize the packed genotypes
		LD::InitPackedGeno();
		LD::LD_Method = *method;
		if (*n_slide <= 0)
			LD::calcLD_mat(*NumThread, out_LD);
		else
			LD::calcLD_slide_mat(*NumThread, out_LD, *n_slide);
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
		const R_xlen_t n = MCWorkingGeno.Space.SNPNum();
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


// the functions for individual dissimilarity

/// to compute the individual dissimilarity
DLLEXPORT void gnrDiss(LongBool *Verbose, LongBool *DataCache, int *NumThread,
	double *out_Diss, LongBool *out_err)
{
	CORETRY
		// ******** To cache the genotype data ********
		if (*DataCache)
		{
			double GenoSum=0;
			gnrCacheGeno(&GenoSum, NULL);
			if (*Verbose)
				Rprintf("Dissimilarity:\tthe sum of all working genotypes (0, 1 and 2) = %.0f\n", GenoSum);
		}

		// ******** The calculation of genetic covariance matrix ********

		// the number of samples
		const R_xlen_t n = MCWorkingGeno.Space.SampleNum();
		// to detect the block size
		IBS::AutoDetectSNPBlockSize(n);
		// the upper-triangle genetic covariance matrix
		CdMatTri<IBS::TDissflag> Dist(n);

		// Calculate the genetic distance matrix
		IBS::DoDissCalculate(Dist, *NumThread, "", *Verbose);

		// output
		IBS::TDissflag *p = Dist.get();
		for (int i=0; i < n; i++)
		{
			out_Diss[i*n + i] = 2 * (p->SumGeno / p->SumAFreq);
			p ++;
			for (int j=i+1; j < n; j++, p++)
				out_Diss[i*n + j] = out_Diss[j*n + i] = (p->SumGeno / p->SumAFreq);
		}

		// output
		*out_err = 0;
	CORECATCH(*out_err = 1)
}



// calculate the dissimilarity between two individuals
static double _distance(double *dist, int n_dist, int I[], int N1, int N2)
{
	const int N = N1 + N2;
	double _mean = 0;
	for (int i=0; i < N1; i++)
	{
		for (int j=N1; j < N; j++)
		{
			_mean += dist[ I[i] * n_dist + I[j] ];
		}
	}
	return _mean / (N1*N2);
}

// return 0 .. (Range-1)
static inline int _RandomNum(int Range)
{
	int rv = (int)(unif_rand()*(Range-1) + 0.5);
	if (rv >= Range) rv = Range -1;
	return rv;
}

// to compute mean and standard deviation of individual dissimilarity by permutation
static void _RunDissZ(double *dist, int n_dist, int I[], int N1, int N2,
	double OutD[], int nOut)
{
	const int N  = N1 + N2;
	const int NSub1 = (N1 < N2) ? N1 : N2;
	const int NSub2 = N - NSub1;
	vector<int> I_Idx(&I[0], &I[N]);

	for (int cnt=0; cnt < nOut; cnt++)
	{
		// random subset
		for (int i=0; i < NSub1; i++)
		{
			int k = _RandomNum(N - i);
			swap(I_Idx[i], I_Idx[k+i]);
		}
		// calculate distance
		OutD[cnt] = _distance(dist, n_dist, &(I_Idx[0]), NSub1, NSub2);
	}
}


// to compute mean and standard deviation of individual dissimilarity by permutation
DLLEXPORT void gnrDistPerm(int *n_dist, double *dist, int *merge,
	int *n_perm, double *z_threshold,
	double Out_Merge_Z[], int Out_Merge_N1[], int Out_Merge_N2[],
	int Out_Ind_Grp[], int *out_err)
{
	GetRNGstate();

	CORETRY
		//
		vector< vector<int> > Array(*n_dist - 1);
		vector<double> d_buffer(*n_perm);

		for (int i_merge=0; i_merge < (*n_dist - 1); i_merge++)
		{
			// individual index
			int i1 = merge[i_merge];
			int i2 = merge[i_merge + (*n_dist - 1)];
			int n1, n2;
			vector<int> &A = Array[i_merge];

			if (i1 < 0)
			{
				A.push_back((-i1) - 1);
				n1 = 1;
			} else {
				A.insert(A.end(), Array[i1-1].begin(), Array[i1-1].end());
				n1 = Array[i1-1].size();
			}
			if (i2 < 0)
			{
				A.push_back((-i2) - 1);
				n2 = 1;
			} else {
				A.insert(A.end(), Array[i2-1].begin(), Array[i2-1].end());
				n2 = Array[i2-1].size();
			}

			Out_Merge_N1[i_merge] = n1;
			Out_Merge_N2[i_merge] = n2;

			if ((n1<=1) && (n2<=1))
			{
				Out_Merge_Z[i_merge] = 0;
			} else {

				double L = _distance(dist, *n_dist, &(A[0]), n1, n2);
				_RunDissZ(dist, *n_dist, &(A[0]), n1, n2, &(d_buffer[0]), *n_perm);

				// calculate mean
				double _mean = 0;
				for (int i=0; i < *n_perm; i++)
					_mean += d_buffer[i];
				_mean /= *n_perm;

				// standard deviation
				double _sd = 0;
				for (int i=0; i < *n_perm; i++)
					_sd += (d_buffer[i] - _mean)*(d_buffer[i] - _mean);
				_sd /= (*n_perm - 1);

				Out_Merge_Z[i_merge] = (_sd > 0) ? ((L - _mean)/sqrt(_sd)) : 0;
			}
		}

		// determine groups
		vector<int> grp_flag(*n_dist - 1, 0);
		for (int i=0; i < *n_dist; i++) Out_Ind_Grp[i] = 1;

		for (int i_merge=0; i_merge < (*n_dist - 1); i_merge++)
		{
			bool b = (Out_Merge_Z[i_merge] >= *z_threshold);
			if (!b)
			{
				// individual index
				int i1 = merge[i_merge];
				int i2 = merge[i_merge + (*n_dist - 1)];
				if ((i1 > 0) && (grp_flag[i1-1] != 0))
					b = true;
				if ((i2 > 0) && (grp_flag[i2-1] != 0))
					b = true;
			}
			if (b)
			{
				grp_flag[i_merge] = 1;

				const int N1 = Out_Merge_N1[i_merge];
				const int N2 = Out_Merge_N2[i_merge];
				const vector<int> &A = Array[i_merge];

				int max = 0;
				for (int i=0; i < N1; i++)
				{
					if (Out_Ind_Grp[A[i]] > max)
						max = Out_Ind_Grp[ A[i] ];
				}
				for (int i=N1; i < N1+N2; i++)
				{
					Out_Ind_Grp[A[i]] += max;
				}
			}
		}

		// output
		*out_err = 0;
	CORECATCH(*out_err = 1)

	PutRNGstate();
}




// conversion

/// to convert from GDS to PLINK PED
DLLEXPORT void gnrConvGDS2PED(char **pedfn, char *SampID[], char *Allele[],
	int *fmt_code, LongBool *verbose, LongBool *out_err)
{
	CORETRY
		MCWorkingGeno.Progress.Info = "\t\tOutput: ";
		MCWorkingGeno.Progress.Show() = *verbose;
		MCWorkingGeno.Progress.Init(MCWorkingGeno.Space.SampleNum());

		ofstream file(*pedfn);
		if (!file.good())
			throw ErrCoreArray("Fail to create the file '%s'.", *pedfn);

		const char *s;
		string s1, s2;
		const int fmt = fmt_code[0];
		CdBufSpace buf(MCWorkingGeno.Space, false, CdBufSpace::acInc);

		for (long i=0; i < buf.IdxCnt(); i++)
		{
			file << "0\t" << SampID[i] << "\t0\t0\t0\t-9";
			UInt8 *g = buf.ReadGeno(i);

			for (long j=0; j < MCWorkingGeno.Space.SNPNum(); j++, g++)
			{
				switch (fmt)
				{
				case 1:
					// SNP alleles: allelic codes
					s1.clear(); s2.clear();
					s = Allele[j];
					while ((*s != 0) && (*s != '/'))
					{
						s1.push_back(*s);
						s ++;
					}
					if (*s == '/') s++;
					while ((*s != 0) && (*s != '/'))
					{
						s2.push_back(*s);
						s ++;
					}
					if (s1.empty()) s1 = "0";
					if (s2.empty()) s2 = "0";

					switch (*g)
					{
						case 0:
							file << "\t" << s2 << " " << s2; break;
						case 1:
							file << "\t" << s1 << " " << s2; break;
						case 2:
							file << "\t" << s1 << " " << s1; break;
						default:
							file << "\t0 0";
					}
					break;

				case 2:
					// A/B codes
					s = (*g==0) ? "B B" : ((*g==1)?"A B": ((*g==2)?"A A":"0 0"));
					file << "\t" << s;
					break;

				case 3:
					// 1/2 codes
					s = (*g==0) ? "2 2" : ((*g==1)?"1 2": ((*g==2)?"1 1":"0 0"));
					file << "\t" << s;
					break;
				}
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

		ofstream file(*bedfn, ios::binary);
		if (!file.good())
			throw ErrCoreArray("Fail to create the file '%s'.", *bedfn);
		// output prefix
		{
			char prefix[3];
			prefix[0] = 0x6C; prefix[1] = 0x1B;
			prefix[2] = (*SNPOrder) ? 0 : 1;
			file.write(prefix, 3);
		}

		CdBufSpace buf(MCWorkingGeno.Space, !(*SNPOrder), CdBufSpace::acInc);
		MCWorkingGeno.Progress.Init(buf.IdxCnt());

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
			throw ErrCoreArray("Fail to create the file '%s'.", *pedfn);
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
			throw ErrCoreArray("Cannot open the file '%s'.", *bedfn);
		char prefix[3];
		file.read(prefix, 3);
		if ((prefix[0] != 0x6C) || (prefix[1] != 0x1B))
			throw ErrCoreArray("Invalid prefix in the bed file.");
		*SNPOrder = (UInt8)(prefix[2]);
		*out_err = 0;
	CORECATCH(*out_err = 1)
}

/// to convert from PLINK BED to GDS
DLLEXPORT void gnrConvBED2GDS(char **bedfn, PdSequenceX *Node, LongBool *verbose,
	LongBool *out_err)
{
	CORETRY
		PdSequenceX Seq = *Node;
		int DLen[2];
		gds_SeqGetDim(Seq, DLen);

		MCWorkingGeno.Progress.Info = " ";
		MCWorkingGeno.Progress.Show() = *verbose;
		MCWorkingGeno.Progress.Init(DLen[0]);

		ifstream file(*bedfn, ios::binary);
		if (!file.good())
			throw ErrCoreArray("Fail to open the file '%s'.", *bedfn);
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

/// to detect and correct strand problem

static inline bool ATGC(const string &s)
{
	return (s=="A") || (s=="T") || (s=="G") || (s=="C");
}
static inline int ALLELE_MINOR(double freq)
{
	return (freq <= 0.5) ? 0 : 1;
}
static inline void split_allele(const char *txt, string &allele1, string &allele2)
{
	const char *p = strchr(txt, '/');
	if (p != NULL)
	{
		// the first allele
		allele1.assign(txt, p);
		for (unsigned int i=0; i < allele1.size(); i++)
			allele1[i] = toupper(allele1[i]);

		// the second allele	
		allele2 = p + 1;
		for (unsigned int i=0; i < allele2.size(); i++)
			allele2[i] = toupper(allele2[i]);
	} else {
		// no a second allele
		allele1 = txt;
		for (unsigned int i=0; i < allele1.size(); i++)
			allele1[i] = toupper(allele1[i]);
		allele2.clear();
	}
}

DLLEXPORT void gnrAlleleStrand(char *allele1[], double afreq1[], int I1[],
	char *allele2[], double afreq2[], int I2[],
	int *if_same_strand, int *n,
	LongBool out_flag[], int *out_n_stand_amb, int *out_n_mismatch,
	LongBool *out_err)
{
	CORETRY
		// initialize: A-T pair, C-G pair
		map<string, string> MAP;
		MAP["A"] = "T"; MAP["C"] = "G"; MAP["G"] = "C"; MAP["T"] = "A";

		*out_n_stand_amb = 0;
		*out_n_mismatch = 0;

		const bool check_strand = (if_same_strand[0] == 0);

		// loop for each SNP
		for (int i=0; i < *n; i++)
		{
			// if true, need switch strand
			bool switch_flag = false;

			// if true, need to compare the allele frequencies
			//   0 -- no switch detect
			//   1 -- detect whether switch or not for stand ambiguity
			//   2 -- detect whether switch or not for mismatching alleles
			int switch_freq_detect = 0;

			// ``ref / nonref alleles''
			string s1, s2;
			string p1, p2;
			split_allele(allele1[I1[i]-1], s1, s2);
			split_allele(allele2[I2[i]-1], p1, p2);

			// allele frequency
			const double F1 = afreq1[I1[i]-1];
			const double F2 = afreq2[I2[i]-1];

			if (ATGC(s1) && ATGC(s2) && ATGC(p1) && ATGC(p2))
			{
				// check
				if ( (s1 == p1) && (s2 == p2) )
				{
					if (check_strand)
					{
						// for example, + C/G <---> - C/G, strand ambi
						if (s1 == MAP[p2])
							switch_freq_detect = 1;
					}
				} else if ( (s1 == p2) && (s2 == p1) )
				{
					if (check_strand)
					{
						// for example, + C/G <---> - G/C, strand ambi
						if (s1 == MAP[p1])
							switch_freq_detect = 1;
						else
							switch_flag = true;
					} else
						switch_flag = true;
				} else {
					if (check_strand)
					{
						if ( (s1 == MAP[p1]) && (s2 == MAP[p2]) )
						{
							// for example, + C/G <---> - G/C, strand ambi
							if (s1 == p2)
								switch_freq_detect = 1;
						} else if ( (s1 == MAP[p2]) && (s2 == MAP[p1]) )
							switch_flag = true;
						else
							switch_freq_detect = 2;
					} else
						switch_freq_detect = 2;
				}
			} else {
				if ((s1 == p1) && (s2 == p2))
				{
					if (s1 == s2)
						switch_freq_detect = 1;  // ambiguous
				} else if ((s1 == p2) && (s2 == p1))
				{
					if (s1 == s2)
						switch_freq_detect = 1;  // ambiguous
					else
						switch_flag = true;
				} else
					switch_freq_detect = 2;
			}

			if (switch_freq_detect != 0)
			{
				switch_flag = (ALLELE_MINOR(F1) != ALLELE_MINOR(F2));
				if (switch_freq_detect == 1)
					(*out_n_stand_amb) ++;
				else
					(*out_n_mismatch) ++;
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
