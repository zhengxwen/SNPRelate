// ===========================================================
//
// genSlideWin.cpp: The Method of Sliding Windows
//
// Copyright (C) 2015    Xiuwen Zheng
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

// using namespace
using namespace std;
using namespace GWAS;


extern "C"
{
/// to compute Fst
extern SEXP gnrFst(SEXP Pop, SEXP nPop, SEXP Method);
extern SEXP gnrSNPRateFreq();

/// Get the list element named str, or return NULL
static SEXP GetListElement(SEXP list, const char *str)
{
	SEXP elmt = R_NilValue;
	SEXP names = getAttrib(list, R_NamesSymbol);
	R_xlen_t n = XLENGTH(list);
	for (R_xlen_t i=0; i < n; i++)
	{
		if (strcmp(CHAR(STRING_ELT(names, i)), str) == 0)
		{
			elmt = VECTOR_ELT(list, i);
			break;
		}
	}
	return elmt;
}

/// Get the average value
static double GetMean(SEXP list)
{
	double *p = REAL(list), sum = 0;
	R_xlen_t n = XLENGTH(list), m = 0;
	for (R_xlen_t i=0; i < n; i++)
	{
		if (R_FINITE(*p))
		{
			sum += *p;
			m ++;
		}
		p ++;
	}
	return sum / m;
}

/// Get the number of sliding windows, always > 0
static int SlidingNumWin(int start, int end, int winsize, int shift)
{
	int cnt = 0;
	end -= winsize;
	for (; start <= end; start += shift)
		cnt ++;
	return cnt + 1;
}


/// get the number of sliding windows
COREARRAY_DLL_EXPORT SEXP gnrSlidingNumWin(SEXP Start, SEXP End,
	SEXP WinSize, SEXP Shift)
{
	return ScalarInteger(
		SlidingNumWin(
			Rf_asInteger(Start), Rf_asInteger(End),
			Rf_asInteger(WinSize), Rf_asInteger(Shift))
	);
}


/// get an error message
COREARRAY_DLL_EXPORT SEXP gnrSlidingWindow(SEXP FUNIdx, SEXP WinSize,
	SEXP Shift, SEXP Unit, SEXP WinStart, SEXP AsIs, SEXP chflag,
	SEXP ChrPos, SEXP Param, SEXP Verbose)
{
	const int FunIdx   = Rf_asInteger(FUNIdx);
	const int winsize  = Rf_asInteger(WinSize);
	const int shift	   = Rf_asInteger(Shift);
	const char *c_Unit = CHAR(STRING_ELT(Unit, 0));
	const char *c_AsIs = CHAR(STRING_ELT(AsIs, 0));

	const size_t nChr  = XLENGTH(ChrPos);
	const int *pPos	   = INTEGER(ChrPos);
	const int *pFlag   = LOGICAL(chflag);
	const bool is_basepair = (strcmp(c_Unit, "basepair") == 0);

	if (MCWorkingGeno.Space.TotalSNPNum() != XLENGTH(chflag))
		error("Internal error in 'gnrSlidingWindow': invalid chflag.");

	int nProtected = 0;
	SEXP ans_rv = PROTECT(NEW_LIST(4));
	nProtected ++;

	// for-loop parameters
	int PosMin=INT_MAX, PosMax=-INT_MAX;
	// the number of sliding windows
	int nWin=0;

	// get the range of chpos
	const int *p = pPos;
	for (size_t n=nChr; n > 0; n--)
	{
		if (*p < PosMin) PosMin = *p;
		if (*p > PosMax) PosMax = *p;
		p ++;
	}
	if ((PosMin == NA_INTEGER) || (PosMax == NA_INTEGER))
		error("Internal error in 'gnrSlidingWindow': invalid position.");

	// posrange
	SEXP posrange = PROTECT(NEW_INTEGER(2));
	nProtected ++;
	INTEGER(posrange)[0] = PosMin;
	INTEGER(posrange)[1] = PosMax;
	SET_ELEMENT(ans_rv, 3, posrange);

	// calculate the number of windows
	if (is_basepair)
	{
		if (Rf_isInteger(WinStart))
			PosMin = Rf_asInteger(WinStart);
	} else {
		PosMax = (long)nChr - 1;
		PosMin = 0;
		if (Rf_isInteger(WinStart))
			PosMin = Rf_asInteger(WinStart) - 1;
	}
	nWin = SlidingNumWin(PosMin, PosMax, winsize, shift);

	bool verbose = (Rf_asLogical(Verbose) == TRUE);
	if (verbose)
		Rprintf(", %d windows\n", nWin);
	CdProgression Progress(1, verbose);
	Progress.Init(nWin, true);

	// get the parameter
	SEXP PL[64];
	for (int i=0; i < 64; i++) PL[i] = R_NilValue;
	switch (FunIdx)
	{
	case 1:	 // snpgdsFst
		PL[0] = GetListElement(Param, "population");
		PL[1] = GetListElement(Param, "npop");
		PL[2] = GetListElement(Param, "method");
		break;
	case 2:	 // snpgdsSNPRateFreq
		break;
	}

	// rvlist
	SEXP rvlist = R_NilValue;
	int as_i = -1;
	if (strcmp(c_AsIs, "list") == 0)
	{
		as_i = 0;
		rvlist = PROTECT(NEW_LIST(nWin));
	} else if (strcmp(c_AsIs, "numeric") == 0)
	{
		as_i = 1;
		rvlist = PROTECT(NEW_NUMERIC(nWin));
	} else if (strcmp(c_AsIs, "array") == 0)
	{
		switch (FunIdx)
		{
		case 1:	 // snpgdsFst
			rvlist = PROTECT(allocMatrix(REALSXP, Rf_asInteger(PL[1])+1, nWin));
			break;
		case 2:	 // snpgdsSNPRateFreq
			rvlist = PROTECT(allocMatrix(REALSXP, 3, nWin));
			break;
		default:
			rvlist = PROTECT(NEW_NUMERIC(nWin));
		}
		as_i = 2;
	}
	if (Rf_isNumeric(rvlist))
	{
		double *p = REAL(rvlist);
		size_t n = XLENGTH(rvlist);
		for (size_t i=0; i < n; i++) p[i] = R_NaN;
	}
	SET_ELEMENT(ans_rv, 0, rvlist);
	nProtected ++;

	// nlist
	SEXP nlist = PROTECT(NEW_INTEGER(nWin));
	SET_ELEMENT(ans_rv, 1, nlist);
	nProtected ++;

	// poslist
	SEXP poslist = PROTECT(NEW_NUMERIC(nWin));
	SET_ELEMENT(ans_rv, 2, poslist);
	nProtected ++;

	// for-loop
	int x = PosMin, iWin = 0;
	while (iWin < nWin)
	{
		C_BOOL *pb = MCWorkingGeno.Space.SNPSelection();
		size_t n=XLENGTH(chflag), ip = 0;
		int num = 0;
		double pos_sum = 0;

		if (is_basepair)
		{
			for (size_t i=0; i < n; i++)
			{
				C_BOOL val = 0;
				if (pFlag[i])
				{
					int P = pPos[ip];
					if ((x <= P) && (P < x+winsize))
					{
						pos_sum += P;
						val = 1; num ++;
					}
					ip ++;
				}
				*pb++ = val;
			}
		} else {
			for (size_t i=0; i < n; i++)
			{
				C_BOOL val = 0;
				if (pFlag[i])
				{
					if ((x <= (int)ip) && ((int)ip < x+winsize))
					{
						pos_sum += pPos[ip];
						val = 1; num ++;
					}
					ip ++;
				}
				*pb++ = val;
			}
		}

		MCWorkingGeno.Space.InitSelectionSNPOnly();

		INTEGER(nlist)[iWin] = num;
		if (num > 0)
		{
			SEXP rv;
			switch (FunIdx)
			{
			case 1:	 // snpgdsFst
				rv = gnrFst(PL[0], PL[1], PL[2]);
				switch (as_i)
				{
					case 0:  // list
						SET_ELEMENT(rvlist, iWin, rv);
						break;
					case 1:  // numeric
						REAL(rvlist)[iWin] = Rf_asReal(VECTOR_ELT(rv, 0));
						break;
					case 2:  // array
						int m  = Rf_asInteger(PL[1]);
						int mk = m + 1;
						double *p = REAL(rvlist);
						double *s = REAL(VECTOR_ELT(rv, 1));
						p[iWin*mk + 0] = Rf_asReal(VECTOR_ELT(rv, 0));
						for (int k=0; k < m; k++)
							p[iWin*mk + k + 1] = s[k*m + k];
						break;
				}
				break;

			case 2:	 // snpgdsSNPRateFreq
				rv = gnrSNPRateFreq();
				switch (as_i)
				{
					case 0:  // list
						SET_ELEMENT(rvlist, iWin, rv);
						break;
					case 1:  // numeric
						REAL(rvlist)[iWin] = GetMean(VECTOR_ELT(rv, 0));
						break;
					case 2:  // array
						double *p = REAL(rvlist);
						p[iWin*3 + 0] = GetMean(VECTOR_ELT(rv, 0));
						p[iWin*3 + 1] = GetMean(VECTOR_ELT(rv, 1));
						p[iWin*3 + 2] = GetMean(VECTOR_ELT(rv, 2));
						break;
				}
				break;
			}

			REAL(poslist)[iWin] = pos_sum / num;
		} else {
			REAL(poslist)[iWin] = R_NaN;
		}

		// poslist[i] <- median(ppos)
		x += shift; iWin ++;

		Progress.Forward(1);
	}

	UNPROTECT(nProtected);

	return ans_rv;
}

}

#endif  /* _HEADER_FST_ */
