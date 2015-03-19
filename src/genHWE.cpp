// ===========================================================
//
// genHWE.cpp: Hardy-Weinberg Equilibrium Test on SNPs
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

#ifndef _HEADER_HWE_
#define _HEADER_HWE_

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
//
// This function is a modification of code written by Janis E. Wigginton.
//
// This code implements an exact SNP test of Hardy-Weinberg Equilibrium as
// described in Wigginton, JE, Cutler, DJ, and Abecasis, GR (2005) A Note on
// Exact Tests of Hardy-Weinberg Equilibrium. American Journal of Human
// Genetics. 76, 887-93 (2005).
//
static double SNPHWE_pValue(int obs_hets, int obs_hom1, int obs_hom2,
	double *het_probs)
{
	int obs_homc = (obs_hom1 < obs_hom2) ? obs_hom2 : obs_hom1;
	int obs_homr = (obs_hom1 < obs_hom2) ? obs_hom1 : obs_hom2;

	int rare_copies = 2 * obs_homr + obs_hets;
	int genotypes   = obs_hets + obs_homc + obs_homr;

	if (genotypes <= 0) return R_NaN;

	// double * het_probs = (double *) malloc((size_t) (rare_copies + 1) * sizeof(double));
	memset(het_probs, 0, sizeof(double) * (rare_copies+1));

	// start at midpoint
	int mid = rare_copies * (2 * genotypes - rare_copies) / (2 * genotypes);

	// check to ensure that midpoint and rare alleles have same parity
	if ((rare_copies & 1) ^ (mid & 1))
		mid ++;

	int curr_hets = mid;
	int curr_homr = (rare_copies - mid) / 2;
	int curr_homc = genotypes - curr_hets - curr_homr;

	het_probs[mid] = 1.0;
	double sum = het_probs[mid];
	for (curr_hets = mid; curr_hets > 1; curr_hets -= 2)
	{
		het_probs[curr_hets - 2] = het_probs[curr_hets] * curr_hets *
			(curr_hets - 1.0) / (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0));
		sum += het_probs[curr_hets - 2];

		// 2 fewer heterozygotes for next iteration -> add one rare,
		// one common homozygote
		curr_homr++;
		curr_homc++;
	}

	curr_hets = mid;
	curr_homr = (rare_copies - mid) / 2;
	curr_homc = genotypes - curr_hets - curr_homr;
	for (curr_hets = mid; curr_hets <= rare_copies - 2; curr_hets += 2)
	{
		het_probs[curr_hets + 2] = het_probs[curr_hets] * 4.0 * curr_homr *
			curr_homc / ((curr_hets + 2.0) * (curr_hets + 1.0));
		sum += het_probs[curr_hets + 2];

		// add 2 heterozygotes for next iteration -> subtract one rare,
		// one common homozygote
		curr_homr--;
		curr_homc--;
	}

	for (int i=0; i <= rare_copies; i++)
		het_probs[i] /= sum;

	double p_hwe = 0.0;
	// p-value calculation for p_hwe
	for (int i=0; i <= rare_copies; i++)
	{
		if (het_probs[i] > het_probs[obs_hets])
			continue;
		p_hwe += het_probs[i];
	}

	return (p_hwe > 1.0) ? 1.0 : p_hwe;
}


/// to compute the p-value for the exact SNP test of Hardy-Weinberg Equilibrium
COREARRAY_DLL_EXPORT SEXP gnrHWE()
{
	COREARRAY_TRY

		// the number of samples
		size_t n = MCWorkingGeno.Space.SNPNum();

		vector<int> AA(n), AB(n), BB(n);
		MCWorkingGeno.Space.GetABNumPerSNP(&AA[0], &AB[0], &BB[0]);

		vector<double> het_probs(2*MCWorkingGeno.Space.SampleNum());
		PROTECT(rv_ans = NEW_NUMERIC(n));
		double *p = REAL(rv_ans);

		for (size_t i=0; i < n; i++)
			p[i] = SNPHWE_pValue(AB[i], AA[i], BB[i], &het_probs[0]);

		UNPROTECT(1);

	COREARRAY_CATCH
}

}

#endif  /* _HEADER_HWE_ */
