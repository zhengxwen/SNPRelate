// ===========================================================
//
// genFst.cpp: Fixation index (Fst) Estimation
//
// Copyright (C) 2015-2024    Xiuwen Zheng
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

#ifndef _HEADER_FST_
#define _HEADER_FST_

// Standard library header
#include <vector>
#include <cmath>

// CoreArray library header
#include "dGenGWAS.h"
#include "dVect.h"

namespace Fst
{
	class COREARRAY_DLL_LOCAL ClassFst_WH15
	{
	public:
		std::vector<double> H;

		ClassFst_WH15() {}
		void SetNumPop(int nPop) {}
	};
}


// using namespace
using namespace std;
using namespace GWAS;
using namespace Vectorization;


extern "C"
{

// Fst method (Weir & Cockerham, 1984)
static bool WC84(C_UInt8 *pGeno, size_t nSamp, size_t NumPop, int PopIdx[],
	int ACnt[], int Cnt[], double P[], double &Numerator, double &Denominator)
{
	// initialize
	int ACntTol=0, CntTol=0;
	memset(ACnt, 0, sizeof(int)*NumPop);
	memset(Cnt,  0, sizeof(int)*NumPop);

	// calculate allele count for each population
	for (size_t i=0; i < nSamp; i++)
	{
		int pop = PopIdx[i] - 1;
		C_UInt8 g = *pGeno ++;
		if (g <= 2)
		{
			ACnt[pop] += g; Cnt[pop] += 2;
			ACntTol   += g; CntTol   += 2;
		}
	}

	// check no missing allele frequency
	for (size_t k=0; k < NumPop; k++)
	{
		if (Cnt[k] > 0)
			P[k] = (double)ACnt[k] / Cnt[k];
		else
			return false;
	}

	double P_All = (double)ACntTol / CntTol;
	double MSB=0, MSW=0, n_c=0;
	for (size_t k=0; k < NumPop; k++)
	{
		MSB += Cnt[k] * (P[k] - P_All) * (P[k] - P_All);
		MSW += Cnt[k] * P[k] * (1 - P[k]);
		n_c += Cnt[k] * Cnt[k];
	}
	MSB /= (NumPop - 1);
	MSW /= (CntTol - NumPop);
	n_c = (CntTol - n_c/CntTol) / (NumPop - 1);
	Numerator = (MSB - MSW);
	Denominator = (MSB + (n_c-1)*MSW);
	return true;
}


// Fst method (Weir & Hill, 2002)
static bool WH02(C_UInt8 *pGeno, size_t nSamp, size_t NumPop, int PopIdx[],
	int ACnt[], int Cnt[], double P[], double OutH[])
{
	// initialize
	memset(ACnt, 0, sizeof(int)*NumPop);
	memset(Cnt,  0, sizeof(int)*NumPop);

	// calculate allele count for each population
	for (size_t i=0; i < nSamp; i++)
	{
		int pop = PopIdx[i] - 1;
		C_UInt8 g = *pGeno ++;
		if (g <= 2)
		{
			ACnt[pop] += g;
			Cnt[pop]  += 2;
		}
	}

	// check no missing allele frequency
	for (size_t k=0; k < NumPop; k++)
	{
		if (Cnt[k] > 0)
			P[k] = (double)ACnt[k] / Cnt[k];
		else
			return false;
	}

	for (size_t k1=0; k1 < NumPop; k1++)
	{
		OutH[k1*NumPop + k1] =
			2.0 * Cnt[k1] / (Cnt[k1] - 1) * P[k1] * (1 - P[k1]);
		for (size_t k2=k1+1; k2 < NumPop; k2++)
		{
			OutH[k1*NumPop + k2] = P[k1] + P[k2] - 2*P[k1]*P[k2];
		}
	}
	return true;
}

static double WH02_beta(size_t nSamp, size_t NumPop, double H[], double *beta)
{
	double H_W=0, H_B=0;
	for (size_t k1=0; k1 < NumPop; k1++)
	{
		H_W += H[k1*NumPop + k1];
		for (size_t k2=k1+1; k2 < NumPop; k2++)
			H_B += H[k1*NumPop + k2];
	}
	H_W /= NumPop;
	H_B /= NumPop * (NumPop-1) / 2;
	if (beta)
	{
		for (size_t k1=0; k1 < NumPop; k1++)
		{
			for (size_t k2=k1; k2 < NumPop; k2++)
			{
				beta[k1*NumPop + k2] = beta[k2*NumPop + k1] =
					1 - H[k1*NumPop + k2] / H_B;
			}
		}
	}
	return 1 - H_W/H_B;
}


/// Calculate Fst estimtes
COREARRAY_DLL_EXPORT SEXP gnrFst(SEXP Pop, SEXP nPop, SEXP Method)
{
	int *PopIdx = INTEGER(Pop);
	int NumPop = Rf_asInteger(nPop);
	const char *MethodText = CHAR(STRING_ELT(Method, 0));

	COREARRAY_TRY

		const int nSamp = MCWorkingGeno.Space().SampleNum();
		const int nSNP  = MCWorkingGeno.Space().SNPNum();
		CdBufSpace BufSNP(MCWorkingGeno.Space(), true, CdBufSpace::acInc);

		vector<int> ACnt(NumPop), Cnt(NumPop);
		vector<double> P(NumPop);

		if (strcmp(MethodText, "W&C84") == 0)
		{
			// Weir & Cockerham, 1984
			double Numerator=0, Denominator=0, num, denom;
			SEXP ratio = PROTECT(NEW_NUMERIC(nSNP));

			// for-loop each SNP
			for (int i=0; i < nSNP; i++)
			{
				double r = R_NaN;
				if (WC84(BufSNP.ReadGeno(i), nSamp, NumPop, PopIdx,
					&ACnt[0], &Cnt[0], &P[0], num, denom))
				{
					Numerator += num;
					Denominator += denom;
					r = num / denom;
				}
				REAL(ratio)[i] = r;
			}

			// output
			PROTECT(rv_ans = NEW_LIST(2));
			SET_ELEMENT(rv_ans, 0, ScalarReal(Numerator/Denominator));
			SET_ELEMENT(rv_ans, 1, ratio);
			UNPROTECT(2);

		} else if (strcmp(MethodText, "W&H02") == 0)
		{
			// Weir & Hill, 2002
			vector<double> H(NumPop*NumPop, 0);
			vector<double> SumH(NumPop*NumPop, 0);
			SEXP ratio = PROTECT(NEW_NUMERIC(nSNP));

			// for-loop each SNP
			for (int i=0; i < nSNP; i++)
			{
				double r = R_NaN;
				if (WH02(BufSNP.ReadGeno(i), nSamp, NumPop, PopIdx,
					&ACnt[0], &Cnt[0], &P[0], &H[0]))
				{
					r = WH02_beta(nSamp, NumPop, &H[0], NULL);
					vec_f64_add(&SumH[0], &H[0], NumPop*NumPop);
				}
				REAL(ratio)[i] = r;
			}

			// output
			PROTECT(rv_ans = NEW_LIST(3));
			SEXP beta = PROTECT(Rf_allocMatrix(REALSXP, NumPop, NumPop));
			double r = WH02_beta(nSamp, NumPop, &SumH[0], REAL(beta));
			SET_ELEMENT(rv_ans, 0, ScalarReal(r));
			SET_ELEMENT(rv_ans, 1, ratio);
			SET_ELEMENT(rv_ans, 2, beta);
			UNPROTECT(3);
		}

	COREARRAY_CATCH
}

}

#endif  /* _HEADER_FST_ */
