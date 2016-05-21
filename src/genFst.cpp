// ===========================================================
//
// genFst.cpp: Fixation index (Fst) Estimation
//
// Copyright (C) 2015-2016    Xiuwen Zheng
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

// CoreArray library header
#include "dGenGWAS.h"

// Standard library header
#include <vector>
#include <cmath>

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


extern "C"
{
/// to compute Fst
COREARRAY_DLL_EXPORT SEXP gnrFst(SEXP Pop, SEXP nPop, SEXP Method)
{
	int *PopIdx = INTEGER(Pop);
	int NumPop = asInteger(nPop);
	const char *MetText = CHAR(STRING_ELT(Method, 0));

	COREARRAY_TRY

		const int nSamp = MCWorkingGeno.Space().SampleNum();
		CdBufSpace BufSNP(MCWorkingGeno.Space(), true, CdBufSpace::acInc);

		if (strcmp(MetText, "W&H02") == 0)
		{
			// Weir & Hill, 2002
			vector<double> H(NumPop*NumPop, 0);
			vector<int> ACnt(NumPop), Cnt(NumPop);
			vector<double> P(NumPop);

			// for-loop each SNP
			for (int i=0; i < MCWorkingGeno.Space().SNPNum(); i++)
			{
				// read genotypes
				C_UInt8 *pg = BufSNP.ReadGeno(i);

				// calculate allele count
				memset(&(ACnt[0]), 0, sizeof(int)*NumPop);
				memset(&(Cnt[0]),  0, sizeof(int)*NumPop);
				for (int j=0; j < nSamp; j++)
				{
					int po = PopIdx[j];
					if (po == NA_INTEGER)
						throw "'pop' should not have any missing value.";
					po --;
					C_UInt8 g = *pg ++;
					if (g <= 2)
					{
						ACnt[po] += g;
						Cnt[po] += 2;
					}
				}

				// check no missing allele frequency
				bool valid=true;
				for (int k=0; k < NumPop; k++)
				{
					if (Cnt[k] > 0)
					{
						P[k] = (double)ACnt[k] / Cnt[k];
					} else {
						valid = false;
						break;
					}	
				}
				if (valid)
				{
					for (int k1=0; k1 < NumPop; k1++)
					{
						H[k1*NumPop + k1] +=
							2.0 * Cnt[k1] / (Cnt[k1] - 1) *
							P[k1] * (1 - P[k1]);
						for (int k2=k1+1; k2 < NumPop; k2++)
						{
							H[k1*NumPop + k2] +=
								P[k1] + P[k2] - 2*P[k1]*P[k2];
						}
					}
				}
			}

			// compute beta
			double H_W=0, H_B=0;
			for (int k1=0; k1 < NumPop; k1++)
			{
				H_W += H[k1*NumPop + k1];
				for (int k2=k1+1; k2 < NumPop; k2++)
					H_B += H[k1*NumPop + k2];
			}
			H_W /= NumPop;
			H_B /= NumPop * (NumPop-1) / 2;

			// output
			PROTECT(rv_ans = NEW_LIST(2));
			SET_ELEMENT(rv_ans, 0, ScalarReal(1 - H_W/H_B));
			SEXP beta_mat = PROTECT(Rf_allocMatrix(REALSXP, NumPop, NumPop));
			double *pmat = REAL(beta_mat);
			SET_ELEMENT(rv_ans, 1, beta_mat);
			for (int k1=0; k1 < NumPop; k1++)
			{
				for (int k2=k1; k2 < NumPop; k2++)
				{
					pmat[k1*NumPop + k2] = pmat[k2*NumPop + k1] =
						1 - H[k1*NumPop + k2] / H_B;
				}
			}
			UNPROTECT(2);

		} else if (strcmp(MetText, "W&C84") == 0)
		{
			// Weir & Cockerham, 1984
			vector<int> ACnt(NumPop), Cnt(NumPop);
			vector<double> P(NumPop);
			double Numerator=0, Denominator=0;

			// for-loop each SNP
			for (int i=0; i < MCWorkingGeno.Space().SNPNum(); i++)
			{
				// read genotypes
				C_UInt8 *pg = BufSNP.ReadGeno(i);

				// calculate allele count
				int ACntTol=0, CntTol=0;
				memset(&(ACnt[0]), 0, sizeof(int)*NumPop);
				memset(&(Cnt[0]),  0, sizeof(int)*NumPop);
				for (int j=0; j < nSamp; j++)
				{
					int po = PopIdx[j] - 1;
					C_UInt8 g = *pg ++;
					if (g <= 2)
					{
						ACnt[po] += g; Cnt[po] += 2;
						ACntTol += g; CntTol += 2;
					}
				}

				// check no missing allele frequency
				bool valid=true;
				for (int k=0; k < NumPop; k++)
				{
					if (Cnt[k] > 0)
					{
						P[k] = (double)ACnt[k] / Cnt[k];
					} else {
						valid = false;
						break;
					}	
				}
				if (valid)
				{
					double P_All = (double)ACntTol / CntTol;
					double MSB=0, MSW=0, n_c=0;
					for (int k=0; k < NumPop; k++)
					{
						MSB += Cnt[k] * (P[k]-P_All) * (P[k]-P_All);
						MSW += Cnt[k] * P[k] * (1-P[k]);
						n_c += Cnt[k] * Cnt[k];
					}
					MSB /= (NumPop - 1);
					MSW /= (CntTol - NumPop);
					n_c = (CntTol - n_c/CntTol) / (NumPop - 1);

					Numerator += (MSB - MSW);
					Denominator += (MSB + (n_c-1)*MSW);
				}
			}

			// output
			PROTECT(rv_ans = NEW_LIST(1));
			SET_ELEMENT(rv_ans, 0, ScalarReal(Numerator/Denominator));
			UNPROTECT(1);
		}

	COREARRAY_CATCH
}

}

#endif  /* _HEADER_FST_ */
