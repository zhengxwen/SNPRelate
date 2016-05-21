// ===========================================================
//
// genRec.cpp:
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


#ifndef _HEADER_REC_
#define _HEADER_REC_

// CoreArray library header
#include <dGenGWAS.h>

// Standard library header
#include <cmath>
#include <memory>
#include <algorithm>


namespace Recombination
{
	// using namespace
	using namespace std;
	using namespace CoreArray;
	using namespace Vectorization;
	using namespace GWAS;
}


using namespace GWAS;
using namespace Vectorization;

extern "C"
{
// =======================================================================
// Genetic IBD matrix
// =======================================================================

// Missing value
#define M          0x7F


// =======================================================================

// Genetic IBD matrix: (1 - \theta_i) / (1 - \theta_B)
//                     1 - M_i = (1 - (1 - g_i)^2)
static C_Int8 OneLocusDiag[4] = { 0, 2, 0, M };  // 2*(1 - M_i)

// Genetic IBD matrix: (1 - \theta_{ij}) / (1 - \theta_B)
//                     1 - M_{ij} = 1/2 * (1 - (1 - g_i)(1 - g_j))
static C_Int8 OneLocusOffDiag[4][4] = {  // 2*(1 - M_{ij})
	{ 0, 1, 2, M }, { 1, 1, 1, M }, { 2, 1, 0, M }, { M, M, M, M }
};

static void OneLocusRun(CdMatTri<double> &IBD, size_t nSamp, size_t nSNP,
	bool verbose)
{
	IBD.Clear(0);
	CdMatTri<double> Denom(nSamp, 0);
	CdMatTri<C_Int8> M_ij(nSamp);

	CdBufSpace BufSNP(MCWorkingGeno.Space(), true, CdBufSpace::acInc);
	CdProgression Progress(0, verbose);
	Progress.Init(nSNP, true);

	for (size_t iSNP=0; iSNP < nSNP; iSNP++)
	{
		const C_UInt8 *pGeno = BufSNP.ReadGeno(iSNP);
		C_Int8 *pM = M_ij.Get();
		C_Int64 Sum = 0;
		int nSum = 0;
		for (size_t i=0; i < nSamp; i++)
		{
			const C_UInt8 g_i = pGeno[i];
			*pM++ = OneLocusDiag[g_i];
			for (size_t j=i+1; j < nSamp; j++)
			{
				const C_UInt8 g_j = pGeno[j];
				C_Int8 v = (*pM++) = OneLocusOffDiag[g_i][g_j];
				if (v != M)
					{ Sum += v; nSum ++; }
			}
		}

		if (nSum > 0)
		{
			double OneMb = (double)Sum / nSum;
			double *pI = IBD.Get(), *pD = Denom.Get();
			pM = M_ij.Get();
			for (size_t i=IBD.Size(); i > 0; i--)
			{
				if (*pM != M)
				{
					*pI += *pM;
					*pD += OneMb;
				}
				pI ++; pD ++; pM ++;
			}
		}

		Progress.Forward(1);
	}

	// the ratio
	vt<double>::Div(IBD.Get(), IBD.Get(), Denom.Get(), IBD.Size());
}


// =======================================================================

// Genetic IBD matrix: (1 - \theta_i) / (1 - \theta_B)
//                     M_{ij} = 1/(n_i * n_j) * \sum_A n_{A_i} n_{A_j}
static C_Int8 TwoLociDiag[4][4] = {  // 4 * (1 - M_i)
		{ 0, 4, 0, M }, { 4, 4, 4, M }, { 0, 4, 0, M }, { M, M, M, M }
	};

static C_Int8 TwoLociOffDiag[4][4][4][4] = {  // 4*(1 - M_{ij})
	{
		// g_i = (0, 0)
		{ { 0, 2, 4, M }, { 2, 3, 4, M }, { 4, 4, 4, M }, { M, M, M, M } },
		// g_i = (0, 1)
		{ { 2, 2, 2, M }, { 3, 3, 3, M }, { 4, 4, 4, M }, { M, M, M, M } },
		// g_i = (0, 2)
		{ { 4, 2, 0, M }, { 4, 3, 2, M }, { 4, 4, 4, M }, { M, M, M, M } },
		// g_i = (0, NA)
		{ { M, M, M, M }, { M, M, M, M }, { M, M, M, M }, { M, M, M, M } }
	},
	{
		// g_i = (1, 0)
		{ { 2, 3, 4, M }, { 2, 3, 4, M }, { 2, 3, 4, M }, { M, M, M, M } },
		// g_i = (1, 1)
		{ { 3, 3, 3, M }, { 3, 3, 3, M }, { 3, 3, 3, M }, { M, M, M, M } },
		// g_i = (1, 2)
		{ { 4, 3, 2, M }, { 4, 3, 2, M }, { 4, 3, 2, M }, { M, M, M, M } },
		// g_i = (1, NA)
		{ { M, M, M, M }, { M, M, M, M }, { M, M, M, M }, { M, M, M, M } }
	},
	{
		// g_i = (2, 0)
		{ { 4, 4, 4, M }, { 2, 3, 4, M }, { 0, 2, 4, M }, { M, M, M, M } },
		// g_i = (2, 1)
		{ { 4, 4, 4, M }, { 3, 3, 3, M }, { 2, 2, 2, M }, { M, M, M, M } },
		// g_i = (2, 2)
		{ { 4, 4, 4, M }, { 4, 3, 2, M }, { 4, 2, 0, M }, { M, M, M, M } },
		// g_i = (2, NA)
		{ { M, M, M, M }, { M, M, M, M }, { M, M, M, M }, { M, M, M, M } }
	},
	{
		// g_i = (NA, 0)
		{ { M, M, M, M }, { M, M, M, M }, { M, M, M, M }, { M, M, M, M } },
		// g_i = (NA, 1)
		{ { M, M, M, M }, { M, M, M, M }, { M, M, M, M }, { M, M, M, M } },
		// g_i = (NA, 2)
		{ { M, M, M, M }, { M, M, M, M }, { M, M, M, M }, { M, M, M, M } },
		// g_i = (NA, NA)
		{ { M, M, M, M }, { M, M, M, M }, { M, M, M, M }, { M, M, M, M } }
	}
};

static void TwoLociRun(CdMatTri<double> &IBD, size_t nSamp, size_t nSNP,
	bool verbose)
{
	IBD.Clear(0);
	CdMatTri<double> Denom(nSamp, 0);
	CdMatTri<C_Int8> M_ij(nSamp);

	CdBufSpace BufSNP(MCWorkingGeno.Space(), true, CdBufSpace::acInc);
	CdProgression Progress(0, verbose);
	Progress.Init(nSNP-1, true);

	vector<C_UInt8> GenoBuffer(nSamp);
	size_t iSNP = 0;
	const C_UInt8 *pGeno1 = &GenoBuffer[0];
	const C_UInt8 *pGeno2 = BufSNP.ReadGeno(iSNP);
	iSNP ++;

	for (; iSNP < nSNP; iSNP++)
	{
		memcpy((void*)pGeno1, (void*)pGeno2, nSamp);
		pGeno2 = BufSNP.ReadGeno(iSNP);

		C_Int8 *pM = M_ij.Get();
		C_Int64 Sum = 0;
		int nSum = 0;
		for (size_t i=0; i < nSamp; i++)
		{
			const C_UInt8 g1_i = pGeno1[i];
			const C_UInt8 g2_i = pGeno2[i];
			*pM++ = TwoLociDiag[g1_i][g2_i];
			for (size_t j=i+1; j < nSamp; j++)
			{
				const C_UInt8 g1_j = pGeno1[j];
				const C_UInt8 g2_j = pGeno2[j];
				C_Int8 v = (*pM++) = TwoLociOffDiag[g1_i][g2_i][g1_j][g2_j];
				if (v != M)
					{ Sum += v; nSum ++; }
			}
		}

		if (nSum > 0)
		{
			double OneMb = (double)Sum / nSum;
			double *pI = IBD.Get(), *pD = Denom.Get();
			pM = M_ij.Get();
			for (size_t i=IBD.Size(); i > 0; i--)
			{
				if (*pM != M)
				{
					*pI += *pM;
					*pD += OneMb;
				}
				pI ++; pD ++; pM ++;
			}
		}

		Progress.Forward(1);
	}

	// the ratio
	vt<double>::Div(IBD.Get(), IBD.Get(), Denom.Get(), IBD.Size());
}



/// to compute the eigenvalues and eigenvectors
COREARRAY_DLL_EXPORT SEXP gnrRecombination(SEXP _NumThread, SEXP _Method,
	SEXP _Verbose)
{
	// const int nThread  = Rf_asInteger(_NumThread);
	const char *Method = CHAR(STRING_ELT(_Method, 0));
	const bool verbose = SEXP_Verbose(_Verbose);

	COREARRAY_TRY

		// ======== To cache the genotype data ========
		CachingSNPData("IBD Calculation", verbose);

		// ======== The calculation of IBD matrix ========

		// the number of samples
		const R_xlen_t nSamp = MCWorkingGeno.Space().SampleNum();
		const R_xlen_t nSNP  = MCWorkingGeno.Space().SNPNum();

		// the upper-triangle IBD matrix
		CdMatTri<double> IBD(nSamp);

		if (strcmp(Method, "OneLocus") == 0)
		{
			OneLocusRun(IBD, nSamp, nSNP, verbose);
		} else if (strcmp(Method, "TwoLoci") == 0)
		{
			TwoLociRun(IBD, nSamp, nSNP, verbose);
		} else
			throw ErrCoreArray("Invalid 'method'!");

		// Output
		PROTECT(rv_ans = Rf_allocMatrix(REALSXP, nSamp, nSamp));
		double *base = REAL(rv_ans);
		double *p = IBD.Get();
		for (R_xlen_t i=0; i < nSamp; i++)
		{
			for (R_xlen_t j=i; j < nSamp; j++)
			{
				base[i*nSamp + j] = base[j*nSamp + i] = *p;
				p ++;
			}
		}
		UNPROTECT(1);

	COREARRAY_CATCH
}

}

#endif  /* _HEADER_REC_ */
