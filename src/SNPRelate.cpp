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
// Copyright (C) 2011 - 2014	Xiuwen Zheng [zhengxwen@gmail.com]
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


#include <dGenGWAS.h>
#include <dVect.h>

#include <Rinternals.h>
#include <R_ext/Lapack.h>
#include <R_ext/Rdynload.h>

#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <utility>
#include <algorithm>


using namespace std;
using namespace CoreArray;
using namespace CoreArray::Vectorization;
using namespace GWAS;


#define LongBool int

extern "C"
{

// ===========================================================
// the public functions
// ===========================================================

/// initialize the package, return the flag of SSE2
COREARRAY_DLL_EXPORT SEXP gnrInit()
{
	Init_GDS_Routines();

	SEXP rv_ans = NEW_LOGICAL(1);
#ifdef COREARRAY_SIMD_SSE2
	LOGICAL(rv_ans)[0] = TRUE;
#else
	LOGICAL(rv_ans)[0] = FALSE;
#endif
	return rv_ans;
}



// the functions for the SNP genotype working space

/// set the genotype node
COREARRAY_DLL_EXPORT void gnrSetGenoSpace(PdSequenceX *Node,
	LongBool *SelSamp, LongBool *SelSampFlag, LongBool *SelSNP,
	LongBool *SelSNPFlag, int *out_nsnp, int *out_nsamp,
	LongBool *out_err)
{
	CORE_TRY
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
	CORE_CATCH(*out_err = 1)
}

/// get the dimension of SNP genotypes
COREARRAY_DLL_EXPORT void gnrGetGenoDim(int *out_nsnp, int *out_nsamp)
{
	*out_nsnp = MCWorkingGeno.Space.SNPNum();
	*out_nsamp = MCWorkingGeno.Space.SampleNum();
}

/// set the SNP selection on the genotype node
COREARRAY_DLL_EXPORT void gnrSetGenoSelSNP(LongBool SNP_Flag[], LongBool *out_err)
{
	CORE_TRY
		vector<C_BOOL> buf(MCWorkingGeno.Space.SNPNum());
		for (int i=0; i < MCWorkingGeno.Space.SNPNum(); i++)
			buf[i] = (SNP_Flag[i]!=0);
		MCWorkingGeno.Space.Set_SNPSelection(&buf[0]);
		*out_err = 0;
	CORE_CATCH(*out_err = 1)
}

/** to compute allele freq., minor allele freq. and missing rates
 *  \param remove_mono    whether remove monomorphic snps or not
 *  \param maf            the threshold of minor allele frequencies, keeping ">= maf"
 *  \param missrate       the threshold of missing rates, keeping "<= missing.rate"
**/
COREARRAY_DLL_EXPORT void gnrSelSNP_Base(LongBool *remove_mono, double *maf, double *missrate,
	int *out_num, LongBool *out_selection, LongBool *out_err)
{
	CORE_TRY
		const R_xlen_t n = MCWorkingGeno.Space.SNPNum();
		vector<C_BOOL> sel(n);
		*out_num = MCWorkingGeno.Space.Select_SNP_Base((*remove_mono) != 0,
			*maf, *missrate, &sel[0]);
		for (int i=0; i < n; i++) out_selection[i] = sel[i];
		*out_err = 0;
	CORE_CATCH(*out_err = 1)
}

/** to compute missing rates, and specify allele frequency
 *  \param afreq          specifying the allele frequencies
 *  \param remove_mono    whether remove monomorphic snps or not
 *  \param maf            the threshold of minor allele frequencies, keeping ">= maf"
 *  \param missrate       the threshold of missing rates, keeping "<= missing.rate"
**/
COREARRAY_DLL_EXPORT void gnrSelSNP_Base_Ex(double *afreq, LongBool *remove_mono, double *maf,
	double *missrate, int *out_num, LongBool *out_selection, LongBool *out_err)
{
	CORE_TRY
		const R_xlen_t n = MCWorkingGeno.Space.SNPNum();
		vector<C_BOOL> sel(n);
		*out_num = MCWorkingGeno.Space.Select_SNP_Base_Ex(afreq, (*remove_mono) != 0,
			*maf, *missrate, &sel[0]);
		for (int i=0; i < n; i++) out_selection[i] = sel[i];
		*out_err = 0;
	CORE_CATCH(*out_err = 1)
}


/// to compute allele freq., minor allele freq. and missing rates
COREARRAY_DLL_EXPORT SEXP gnrSNPRateFreq()
{
	COREARRAY_TRY

		R_xlen_t L = MCWorkingGeno.Space.SNPNum();
		SEXP AF, MF, MR;
		PROTECT(rv_ans = NEW_LIST(3));

			PROTECT(AF = NEW_NUMERIC(L));
			SET_ELEMENT(rv_ans, 0, AF);
			PROTECT(MF = NEW_NUMERIC(L));
			SET_ELEMENT(rv_ans, 1, MF);
			PROTECT(MR = NEW_NUMERIC(L));
			SET_ELEMENT(rv_ans, 2, MR);

			MCWorkingGeno.Space.GetAlleleFreqs(REAL(AF));
			double *pAF = REAL(AF), *pMF = REAL(MF);
			for (R_xlen_t i=0; i < L; i++)
				pMF[i] = min(pAF[i], 1 - pAF[i]);
			MCWorkingGeno.Space.GetMissingRates(REAL(MR));

		UNPROTECT(4);

	COREARRAY_CATCH
}


/// to compute allele frequency
COREARRAY_DLL_EXPORT SEXP gnrSNPFreq()
{
	COREARRAY_TRY
		PROTECT(rv_ans = NEW_NUMERIC(MCWorkingGeno.Space.SNPNum()));
		MCWorkingGeno.Space.GetAlleleFreqs(REAL(rv_ans));
		UNPROTECT(1);
	COREARRAY_CATCH
}


/// to compute sample missing rates
COREARRAY_DLL_EXPORT SEXP gnrSampFreq()
{
	COREARRAY_TRY
		PROTECT(rv_ans = NEW_NUMERIC(MCWorkingGeno.Space.SampleNum()));
		MCWorkingGeno.Space.GetSampMissingRates(REAL(rv_ans));
		UNPROTECT(1);
	COREARRAY_CATCH
}


/// to cache the genotype data
COREARRAY_DLL_EXPORT void gnrCacheGeno(double *out_GenoSum, LongBool *out_err)
{
	CORE_TRY
		if (out_GenoSum)
			*out_GenoSum = MCWorkingGeno.Space.GenoSum();
		if (out_err) *out_err = 0;
	CORE_CATCH(if (out_err) *out_err = 1)
}


/// add genotype buffer
COREARRAY_DLL_EXPORT void gnrInitGenoBuffer(LongBool *SNPorSamp, int *AF, int *out_obj, LongBool *out_err)
{
	CORE_TRY
		CdBufSpace *p = new CdBufSpace(MCWorkingGeno.Space,
			*SNPorSamp, CdBufSpace::TAccessFlag(*AF));
		memmove((void*)out_obj, (void*)&p, sizeof(CdBufSpace *));
		*out_err = 0;
	CORE_CATCH(*out_err = 1)
}

/// done genotype buffer
COREARRAY_DLL_EXPORT void gnrDoneGenoBuffer(int *buf_obj, LongBool *out_err)
{
	CORE_TRY
		CdBufSpace *p;
		memmove((void*)&p, (void*)buf_obj, sizeof(CdBufSpace *));
		delete p;
		*out_err = 0;
	CORE_CATCH(*out_err = 1)
}


/// get genotype
COREARRAY_DLL_EXPORT void gnrGetGenoBuffer(int *buf_obj, int *idx, int *out_buf, LongBool *out_err)
{
	CORE_TRY
		CdBufSpace *p;
		memmove((void*)&p, (void*)buf_obj, sizeof(CdBufSpace *));
		C_UInt8 *s = p->ReadGeno(*idx);
		for (long i=0; i < p->BufElmSize(); i++)
			out_buf[i] = s[i];
		*out_err = 0;
	CORE_CATCH(*out_err = 1)
}


/// copy genotypes
COREARRAY_DLL_EXPORT SEXP gnrCopyGeno(SEXP Node, SEXP snpfirstorder)
{
	int snpdim = asLogical(snpfirstorder);
	if (snpdim == NA_LOGICAL)
		error("'snpfirstdim' must be TRUE or FALSE.");

	COREARRAY_TRY

		PdSequenceX obj = GDS_R_SEXP2Obj(Node);
		GDS_R_NodeValid(obj, TRUE);

		if (snpdim)
		{
			C_Int32 cnt[2] = { 1, MCWorkingGeno.Space.SNPNum() };
			CdBufSpace buf(MCWorkingGeno.Space, false, CdBufSpace::acInc);
			for (int i=0; i < buf.IdxCnt(); i++)
			{
				C_UInt8 *p = buf.ReadGeno(i);
				C_Int32 st[2] = { (C_Int32)i, 0 };
				GDS_Seq_wData(obj, st, cnt, p, svUInt8);
			}
		} else {
			C_Int32 cnt[2] = { 1, MCWorkingGeno.Space.SampleNum() };
			CdBufSpace buf(MCWorkingGeno.Space, true, CdBufSpace::acInc);
			for (int i=0; i < buf.IdxCnt(); i++)
			{
				C_UInt8 *p = buf.ReadGeno(i);
				C_Int32 st[2] = { (C_Int32)i, 0 };
				GDS_Seq_wData(obj, st, cnt, p, svUInt8);
			}
		}

	COREARRAY_CATCH
}


/// copy genotypes to a memory buffer
COREARRAY_DLL_EXPORT SEXP gnrCopyGenoMem(SEXP snpfirstdim)
{
	int if_snp = asLogical(snpfirstdim);
	if (if_snp == NA_LOGICAL)
		error("'snpfirstdim' must be TRUE, FALSE or NULL.");

	COREARRAY_TRY

		if (if_snp)
		{
			PROTECT(rv_ans = allocMatrix(INTSXP,
				MCWorkingGeno.Space.SNPNum(), MCWorkingGeno.Space.SampleNum()));
			int *pMem = INTEGER(rv_ans);

			const long n = MCWorkingGeno.Space.SNPNum();
			CdBufSpace buf(MCWorkingGeno.Space, false, CdBufSpace::acInc);
			for (long i=0; i < buf.IdxCnt(); i++)
			{
				C_UInt8 *p = buf.ReadGeno(i);
				for (long k=n; k > 0; k--, p++)
					*pMem++ = (*p <= 2) ? (*p) : NA_INTEGER;
			}
		} else {
			PROTECT(rv_ans = allocMatrix(INTSXP,
				MCWorkingGeno.Space.SampleNum(), MCWorkingGeno.Space.SNPNum()));
			int *pMem = INTEGER(rv_ans);

			const long n = MCWorkingGeno.Space.SampleNum();
			CdBufSpace buf(MCWorkingGeno.Space, true, CdBufSpace::acInc);
			for (long i=0; i < buf.IdxCnt(); i++)
			{
				C_UInt8 *p = buf.ReadGeno(i);
				for (long k=n; k > 0; k--, p++)
					*pMem++ = (*p <= 2) ? (*p) : NA_INTEGER;
			}
		}
		UNPROTECT(1);

	COREARRAY_CATCH
}


/// Append genotypes with switch strand
COREARRAY_DLL_EXPORT void gnrAppendGenoSpaceStrand(PdSequenceX *Node,
	LongBool *snpfirstorder, LongBool StrandFlag[], LongBool *out_err)
{
	CORE_TRY
		*out_err = 1;
		if (*snpfirstorder)
		{
			const R_xlen_t n = MCWorkingGeno.Space.SNPNum();
			CdBufSpace buf(MCWorkingGeno.Space, false, CdBufSpace::acInc);
			for (long i=0; i < buf.IdxCnt(); i++)
			{
				C_UInt8 *p = buf.ReadGeno(i);
				for (int k=0; k < n; k++)
					if (StrandFlag[k] && (p[k]<=2)) p[k] = 2 - p[k];
				GDS_Seq_AppendData(*Node, n, p, svUInt8);
			}
		} else {
			const R_xlen_t n = MCWorkingGeno.Space.SampleNum();
			CdBufSpace buf(MCWorkingGeno.Space, true, CdBufSpace::acInc);
			for (long i=0; i < buf.IdxCnt(); i++)
			{
				C_UInt8 *p = buf.ReadGeno(i);
				if (StrandFlag[i])
				{
					for (int k=0; k < n; k++)
						if (p[k] <= 2) p[k] = 2 - p[k];
				}
				GDS_Seq_AppendData(*Node, n, p, svUInt8);
			}
		}
		*out_err = 0;
	CORE_CATCH(*out_err = 1)
}


/// parse two alleles from a text with format '%s/%s'
static void _ParseAlleleString(const char *str, string &A1, string &A2)
{
	const char *p = str;
	while ((*p != 0) && (*p != '/'))
		p ++;
	A1.assign(str, p);

	if (*p == '/')
		A2.assign(p + 1);
	else
		A2.clear();
}

/// Strand-switching according to A alleles
/** \param Node             the GDS node 'genotype'
 *  \param AlleleNode       the GDS node 'snp.allele'
 *  \param NewAlleleNode    the GDS node '!snp.allele'
 *  \param SNP_First_Dim    TRUE/FALSE
 *  \param A_Allele         characters of A allele
 *  \return a logical vector with 'TRUE' for switching alleles and NA for being unable to determine
**/
COREARRAY_DLL_EXPORT SEXP gnrStrandSwitch(SEXP Node, SEXP AlleleNode,
	SEXP NewAlleleNode, SEXP SNP_First_Dim, SEXP A_Allele)
{
	int snpfirstdim = asLogical(SNP_First_Dim);

	COREARRAY_TRY

		PdGDSObj GenoObj   = GDS_R_SEXP2Obj(Node);
		GDS_R_NodeValid(GenoObj, TRUE);
		PdGDSObj AlleleObj = GDS_R_SEXP2Obj(AlleleNode);
		GDS_R_NodeValid(AlleleObj, TRUE);
		PdGDSObj NewAObj = GDS_R_SEXP2Obj(NewAlleleNode);
		GDS_R_NodeValid(NewAObj, FALSE);

		int NDim[2];
		GDS_Seq_GetDim(GenoObj, NDim, 2);
		long n_snp  = (snpfirstdim ? NDim[1] : NDim[0]);
		long n_samp = (snpfirstdim ? NDim[0] : NDim[1]);
		string val, A, B;

		PROTECT(rv_ans = NEW_LOGICAL(XLENGTH(A_Allele)));
		int *ll_ans = LOGICAL(rv_ans);

		// determine the values of rv_ans		
		for (long i=0; i < n_snp; i++)
		{
			C_Int32 st=i, cnt=1;
			GDS_Seq_rData(AlleleObj, &st, &cnt, &val, svStrUTF8);

			SEXP str = STRING_ELT(A_Allele, i);
			if (str != NA_STRING)
			{
				_ParseAlleleString(val.c_str(), A, B);

				const char *s = CHAR(str);
				if (A.compare(s) == 0)
					ll_ans[i] = FALSE;
				else if (B.compare(s) == 0)
					ll_ans[i] = TRUE;
				else
					ll_ans[i] = NA_LOGICAL;
			} else
				ll_ans[i] = NA_LOGICAL;

			if (ll_ans[i] == TRUE)
				GDS_Seq_AppendString(NewAObj, (B + "/" + A).c_str());
			else
				GDS_Seq_AppendString(NewAObj, val.c_str());
		}

		// strand-switching
		if (snpfirstdim)
		{
			vector<C_UInt8> buffer(n_snp);
			C_Int32 st[2] = { 0, 0 };
			C_Int32 cnt[2] = { 1, (int)n_snp };

			for (long i=0; i < n_samp; i++)
			{
				st[0] = i;
				GDS_Seq_rData(GenoObj, st, cnt, &buffer[0], svUInt8);
				for (long j=0; j < n_snp; j++)
				{
					if (ll_ans[j] == TRUE)
					{
						C_UInt8 &g = buffer[j];
						if (g <= 2) g = 2 - g;
					}
				}
				GDS_Seq_wData(GenoObj, st, cnt, &buffer[0], svUInt8);
			}
		} else {
			vector<C_UInt8> buffer(n_samp);
			C_Int32 st[2] = { 0, 0 };
			C_Int32 cnt[2] = { 1, (int)n_samp };

			for (long i=0; i < n_snp; i++)
			{
				if (ll_ans[i] == TRUE)
				{
					st[0] = i;
					GDS_Seq_rData(GenoObj, st, cnt, &buffer[0], svUInt8);
					for (long j=0; j < n_samp; j++)
					{
						C_UInt8 &g = buffer[j];
						if (g <= 2) g = 2 - g;
					}
					GDS_Seq_wData(GenoObj, st, cnt, &buffer[0], svUInt8);
				}
			}
		}

		UNPROTECT(1);

	COREARRAY_CATCH
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
COREARRAY_DLL_EXPORT void gnrDistPerm(int *n_dist, double *dist, int *merge,
	int *n_perm, double *z_threshold,
	double Out_Merge_Z[], int Out_Merge_N1[], int Out_Merge_N2[],
	int Out_Ind_Grp[], int *out_err)
{
	GetRNGstate();

	CORE_TRY
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
	CORE_CATCH(*out_err = 1)

	PutRNGstate();
}




// conversion

/// to convert from GDS to PLINK PED
COREARRAY_DLL_EXPORT void gnrConvGDS2PED(char **pedfn, char *SampID[], char *Allele[],
	int *fmt_code, LongBool *verbose, LongBool *out_err)
{
	CORE_TRY
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
			C_UInt8 *g = buf.ReadGeno(i);

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
	CORE_CATCH(*out_err = 1)
}


/// to convert from GDS to PLINK BED
COREARRAY_DLL_EXPORT SEXP gnrConvGDS2BED(SEXP bedfn, SEXP SNPOrder, SEXP Verbose)
{
	const char *fn = CHAR(STRING_ELT(bedfn, 0));
	int if_snp = (asLogical(SNPOrder) == TRUE);
	int if_verbose = asLogical(Verbose);
	if (if_verbose == NA_LOGICAL)
		error("'verbose' must be TRUE or FALSE.");

	COREARRAY_TRY

		MCWorkingGeno.Progress.Info = "\t";
		MCWorkingGeno.Progress.Show() = (if_verbose == TRUE);

		ofstream file(fn, ios::binary);
		if (!file.good())
			throw ErrCoreArray("Fail to create the file '%s'.", fn);
		// output prefix
		{
			char prefix[3];
			prefix[0] = 0x6C; prefix[1] = 0x1B;
			prefix[2] = (if_snp) ? 0 : 1;
			file.write(prefix, 3);
		}

		CdBufSpace buf(MCWorkingGeno.Space, !if_snp, CdBufSpace::acInc);
		MCWorkingGeno.Progress.Init(buf.IdxCnt());

		long nRe = buf.BufElmSize() % 4;
		long nPack = (nRe > 0) ? (buf.BufElmSize()/4 + 1) : (buf.BufElmSize()/4);
		vector<char> geno(nPack, 0);
		const C_UInt8 cvt[4] = { 3, 2, 0, 1 };

		for (long i=0; i < buf.IdxCnt(); i++)
		{
			C_UInt8 *s = buf.ReadGeno(i);
			char *p = &geno[0];
			for (long k=0; k < buf.BufElmSize()/4; k++, s+=4)
			{
				*p++ = cvt[s[0] & 0x03] | (cvt[s[1] & 0x03] << 2) |
					(cvt[s[2] & 0x03] << 4) | (cvt[s[3] & 0x03] << 6);
			}
			if (nRe > 0)
			{
				C_UInt8 b = 0;
				for (long k=0; k < nRe; k++, s++)
					b |= (cvt[*s & 0x03] << (2*k));
				*p++ = b;
			}
			file.write(&geno[0], nPack);
			MCWorkingGeno.Progress.Forward(1);
		}

	COREARRAY_CATCH
}


/// to convert from GDS to EIGENSOFT
COREARRAY_DLL_EXPORT void gnrConvGDS2EIGEN(char **pedfn, LongBool *verbose, LongBool *out_err)
{
	CORE_TRY
		MCWorkingGeno.Progress.Info = "\tOutput: ";
		MCWorkingGeno.Progress.Show() = *verbose;
		MCWorkingGeno.Progress.Init(MCWorkingGeno.Space.SNPNum());

		ofstream file(*pedfn);
		if (!file.good())
			throw ErrCoreArray("Fail to create the file '%s'.", *pedfn);
		CdBufSpace buf(MCWorkingGeno.Space, true, CdBufSpace::acInc);
		for (long i=0; i < buf.IdxCnt(); i++)
		{
			C_UInt8 *g = buf.ReadGeno(i);
			for (long j=0; j < MCWorkingGeno.Space.SampleNum(); j++, g++)
			{
				int geno = (*g <= 2) ? (*g) : 9;
				file << geno;
			}
			file << endl;
			MCWorkingGeno.Progress.Forward(1);
		}
		*out_err = 0;
	CORE_CATCH(*out_err = 1)
}


/// to detect PLINK BED
COREARRAY_DLL_EXPORT SEXP gnrBEDFlag(SEXP bedfn)
{
	const char *fn = CHAR(STRING_ELT(bedfn, 0));
	COREARRAY_TRY

		ifstream file(fn, ios::binary);
		if (!file.good())
			throw ErrCoreArray("Cannot open the file '%s'.", fn);

		char prefix[3];
		file.read(prefix, 3);
		if ((prefix[0] != 0x6C) || (prefix[1] != 0x1B))
			throw ErrCoreArray("Invalid prefix in the bed file.");
		rv_ans = NEW_INTEGER(1);
		INTEGER(rv_ans)[0] = (C_UInt8)(prefix[2]);

	COREARRAY_CATCH
}


/// to convert from PLINK BED to GDS
COREARRAY_DLL_EXPORT SEXP gnrConvBED2GDS(SEXP bedfn, SEXP Node, SEXP verbose)
{
	const char *fn = CHAR(STRING_ELT(bedfn, 0));
	COREARRAY_TRY

		PdSequenceX Seq = GDS_R_SEXP2Obj(Node);

		int DLen[2];
		GDS_Seq_GetDim(Seq, DLen, 2);

		MCWorkingGeno.Progress.Info = " ";
		MCWorkingGeno.Progress.Show() = (LOGICAL(verbose)[0] == TRUE);
		MCWorkingGeno.Progress.Init(DLen[0]);

		ifstream file(fn, ios::binary);
		if (!file.good())
			throw ErrCoreArray("Fail to open the file '%s'.", fn);
		// read prefix
		char prefix[3];
		file.read(prefix, 3);

		long nRe = DLen[1] % 4;
		long nPack = (nRe > 0) ? (DLen[1]/4 + 1) : (DLen[1]/4);
		vector<char> srcgeno(nPack);
		vector<C_UInt8> dstgeno(DLen[1]);
		C_Int32 st[2] = { 0, 0 }, cnt[2] = { 1, DLen[1] };
		const C_UInt8 cvt[4] = { 2, 3, 1, 0 };

		for (long i=0; i < DLen[0]; i++)
		{
			// read genotypes
			file.read(&srcgeno[0], nPack);
			// unpacked
			C_UInt8 *p = &dstgeno[0];
			for (long k=0; k < DLen[1]/4; k++)
			{
				C_UInt8 g = srcgeno[k];
				*p++ = cvt[g & 0x03]; g >>= 2;
				*p++ = cvt[g & 0x03]; g >>= 2;
				*p++ = cvt[g & 0x03]; g >>= 2;
				*p++ = cvt[g & 0x03];
			}
			if (nRe > 0)
			{
				C_UInt8 g = srcgeno[DLen[1]/4];
				for (long k=0; k < nRe; k++)
				{
					*p++ = cvt[g & 0x03]; g >>= 2;
				}
			}
			// write
			st[0] = i;
			GDS_Seq_wData(Seq, st, cnt, &dstgeno[0], svUInt8);
			MCWorkingGeno.Progress.Forward(1);
		}

	COREARRAY_CATCH
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
static inline void split_allele(const char *txt, string &a1, string &a2)
{
	const char *p = strchr(txt, '/');
	if (p != NULL)
	{
		// the first allele
		a1.assign(txt, p);
		for (unsigned int i=0; i < a1.size(); i++)
			a1[i] = toupper(a1[i]);

		// the second allele	
		a2 = p + 1;
		for (unsigned int i=0; i < a2.size(); i++)
			a2[i] = toupper(a2[i]);
	} else {
		// no a second allele
		a1 = txt;
		for (unsigned int i=0; i < a1.size(); i++)
			a1[i] = toupper(a1[i]);
		a2.clear();
	}
}

COREARRAY_DLL_EXPORT void gnrAlleleStrand(
	char *allele1[], double afreq1[], int I1[],
	char *allele2[], double afreq2[], int I2[],
	int *if_same_strand, int *n,
	LongBool out_flag[], int *out_n_stand_amb, int *out_n_mismatch,
	LongBool *out_err)
{
	CORE_TRY
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

	CORE_CATCH(*out_err = 1)
}


/// parse the chromosome information if it is character-type
COREARRAY_DLL_EXPORT SEXP gnrChromParse(SEXP gdsobj)
{
	COREARRAY_TRY

		PdSequenceX Obj = GDS_R_SEXP2Obj(gdsobj);
		GDS_R_NodeValid(Obj, TRUE);

		C_Int32 st, cnt, TotalCnt;
		GDS_Seq_GetDim(Obj, &TotalCnt, 1);

		set<string> set_val;
		string val;
		int chr_min=INT_MAX, chr_max=-INT_MAX;

		for (int i=0; i < TotalCnt; i++)
		{
			st = i; cnt = 1;
			GDS_Seq_rData(Obj, &st, &cnt, &val, svStrUTF8);

			// whether val is numeric
			char *endptr = (char*)(val.c_str());
			int v = strtol(val.c_str(), &endptr, 10);
			if (endptr != val.c_str())
			{
				// numeric
				if (v < chr_min) chr_min = v;
				if (v > chr_max) chr_max = v;
			} else {
				// non-numeric
				if (!val.empty())
					set_val.insert(val);
			}
		}
		if (chr_min == INT_MAX) chr_min = NA_INTEGER;
		if (chr_max == -INT_MAX) chr_max = NA_INTEGER;

		PROTECT(rv_ans = NEW_LIST(3));
		SET_ELEMENT(rv_ans, 0, ScalarInteger(chr_min));
		SET_ELEMENT(rv_ans, 1, ScalarInteger(chr_max));

		SEXP tmp = PROTECT(NEW_CHARACTER(set_val.size()));
		SET_ELEMENT(rv_ans, 2, tmp);
		int id = 0;
		for (set<string>::iterator i=set_val.begin(); i != set_val.end(); i++)
		{
			SET_STRING_ELT(tmp, id, mkChar(i->c_str()));
			id ++;
		}

		UNPROTECT(2);

	COREARRAY_CATCH
}


/// return a logical vector with 'TRUE' for being in the range
COREARRAY_DLL_EXPORT SEXP gnrChromRangeNumeric(SEXP gdsobj, SEXP ChrMin,
	SEXP ChrMax)
{
	int CMin = INTEGER(ChrMin)[0];
	int CMax = INTEGER(ChrMax)[0];

	COREARRAY_TRY

		PdSequenceX Obj = GDS_R_SEXP2Obj(gdsobj);
		GDS_R_NodeValid(Obj, TRUE);

		C_Int32 st, cnt, TotalCnt;
		GDS_Seq_GetDim(Obj, &TotalCnt, 1);

		PROTECT(rv_ans = NEW_LOGICAL(TotalCnt));
		int *rv_ptr = LOGICAL(rv_ans);

		for (int i=0; i < TotalCnt; i++)
		{
			C_Int32 val;
			st = i; cnt = 1;
			GDS_Seq_rData(Obj, &st, &cnt, &val, svInt32);
			rv_ptr[i] = ((CMin<=val) && (val<=CMax)) ? TRUE : FALSE;
		}

		UNPROTECT(1);

	COREARRAY_CATCH
}


/// return a logical vector with 'TRUE' for numeric values
COREARRAY_DLL_EXPORT SEXP gnrChromParseNumeric(SEXP gdsobj)
{
	COREARRAY_TRY

		PdSequenceX Obj = GDS_R_SEXP2Obj(gdsobj);
		GDS_R_NodeValid(Obj, TRUE);

		C_Int32 st, cnt, TotalCnt;
		GDS_Seq_GetDim(Obj, &TotalCnt, 1);

		PROTECT(rv_ans = NEW_LOGICAL(TotalCnt));
		int *rv_ptr = LOGICAL(rv_ans);

		string val;
		for (int i=0; i < TotalCnt; i++)
		{
			st = i; cnt = 1;
			GDS_Seq_rData(Obj, &st, &cnt, &val, svStrUTF8);

			// whether val is numeric
			char *endptr = (char*)(val.c_str());
			strtol(val.c_str(), &endptr, 10);
			rv_ptr[i] = (endptr != val.c_str()) ? TRUE : FALSE;
		}

		UNPROTECT(1);

	COREARRAY_CATCH
}


/// get an error message
COREARRAY_DLL_EXPORT SEXP gnrErrMsg()
{
	SEXP ans_rv = mkString(GDS_GetError());
	GDS_SetError(NULL);
	return ans_rv;
}


/// initialize the package
COREARRAY_DLL_EXPORT void R_init_SNPRelate(DllInfo *info)
{
	static R_CallMethodDef callMethods[]  = {
		{"gnrInit", (DL_FUNC)&gnrInit, 0},
		{NULL, NULL, 0}
	};

	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

} // extern "C"
