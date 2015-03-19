// ===========================================================
//
// genIBD.cpp: Identity by descent (IBD) Analysis on GWAS
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


#ifndef _HEADER_IBD_
#define _HEADER_IBD_

// CoreArray library header
#include "dGenGWAS.h"
#include "dVect.h"

// Standard library header
#include <limits>
#include <vector>
#include <cmath>
#include <cfloat>
#include <memory>
#include <algorithm>


#ifdef COREARRAY_SIMD_SSE
#include <xmmintrin.h>
#endif
#ifdef COREARRAY_SIMD_SSE2
#include <emmintrin.h>
#endif


namespace GWAS_Math
{
	using namespace std;

	// Nelder-Mead Simplex algorithm
	// J.A. Nelder and R. Mead, A simplex method for function minimization,
	//   Computer Journal vol. 7 (1965), 308®C315.
	//
	// Press, W. H., S. A. Teukolsky, W. T. Vetterling and B. P. Flannery, 2002.
	//   Numerical Recipes in C++: The Art of Scientific Computing,
	//   Ed. 2. Cambridge University Press, Cambridge, UK.

	/// Nelder-Mead Simplex algorithm
	/** Extrapolates by a factor fac through the face of the simplex across from
	 *  the high point, tries it, and replaces the high point if the new point
	 *  is better. **/
	template<typename tfloat, int ndim, typename Functor>
		inline tfloat Simplex_Point_Try(tfloat p[][ndim], tfloat y[],
			tfloat psum[], int ihi, tfloat fac, Functor funk, void *funkdata)
	{
		tfloat fac1 = (1.0-fac)/ndim, fac2 = fac1-fac;
		tfloat ytry, ptry[ndim];
		for (int j=0; j < ndim; j++)
			ptry[j] = psum[j]*fac1 - p[ihi][j]*fac2;
		// Evaluate the function at the trial point.
		ytry = funk(ptry, funkdata);
		if (ytry < y[ihi])
		{
			// If it's better than the highest, then replace the highest.
			y[ihi] = ytry;
			for (int j=0; j < ndim; j++)
			{
				psum[j] += ptry[j] - p[ihi][j];
				p[ihi][j] = ptry[j];
			}
		}
		return ytry;
	}

	/// Nelder-Mead Simplex algorithm
	/** Multidimensional minimization of the function funk(x) where x[1..ndim]
	 *  is a vector in ndim dimensions, by the downhill simplex method of Nelder
	 *  and Mead. The matrix p[1..ndim+1][1..ndim] is input. Its ndim+1 rows are
	 *  ndim-dimensional vectors which are the vertices of the starting simplex.
	 *  On output, p will have been reset to ndim+1 new points all within reltol
	 *  of a minimum function value, outx will be the point corresponding to the
	 *  minimum value outy, and nfunk gives the number of function evaluations
	 *  taken.
	**/
	template<typename tfloat, int ndim, typename Functor>
		void SimplexMin(tfloat p[][ndim], tfloat outx[], tfloat &outy,
		int &nfunk, Functor funk, void *funkdata, tfloat reltol, int nfunkmax)
	{
		int i, ihi, ilo, inhi, j;
		tfloat convtol, sum, ysave, ytry;
		tfloat y[ndim+1], psum[ndim];

		// The components of y are initialized to the values of funk evaluated
		// at the ndim+1 vertices (rows) of p.
		for (i=0; i <= ndim; i++)
			y[i] = funk(p[i], funkdata);
		nfunk = ndim;
		// Set converage tolenance with respect to reltol
		convtol = reltol * (fabs(y[0]) + fabs(reltol));
		if (convtol < numeric_limits<tfloat>::epsilon())
			convtol = numeric_limits<tfloat>::epsilon();

		// Get psum
		for (j=0; j<ndim; j++)
		{
			for (sum=0, i=0; i<=ndim; i++) sum += p[i][j];
			psum[j] = sum;
		}

		// do loop
		while (true)
		{
			ilo = 0;
			// First we must determine which point is the highest (worst),
			// next-highest, and lowest (best), by looping over the points
			// in the simplex.
			ihi = y[0]>y[1] ? (inhi=1,0) : (inhi=0,1);
			for (i=0; i<=ndim; i++)
			{
				if (y[i] <= y[ilo]) ilo=i;
				if (y[i] > y[ihi]) {
					inhi = ihi; ihi = i;
				} else
					if ((y[i] > y[inhi]) && (i!=ihi))
						inhi = i;
			}
			// Compute the fractional range from highest to lowest and return
			// if satisfactory. If returning, put best point and value in
			// outx and outy.
			if (((y[ihi]-y[ilo]) <= convtol) || (nfunk>=nfunkmax))
			{
				outy = y[ilo];
				for (i=0; i<ndim; i++) outx[i] = p[ilo][i];
				break;
			}

			nfunk += 2;
			// Begin a new iteration.
			// First extrapolate by a factor (-1.0)( through the face of the
			// simplex across from the high point, i.e., reflect the simplex from
			// the high point.
			ytry = Simplex_Point_Try(p, y, psum, ihi, (tfloat)-1.0, funk, funkdata);
			if (ytry <= y[ilo])
			{
				// Gives a result better than the best point, so try an additional
				// extrapolation by a factor (2.0).
				ytry = Simplex_Point_Try(p, y, psum, ihi, (tfloat)2.0, funk, funkdata);
			} else if (ytry >= y[inhi]) {
				// The reflected point is worse than the second-highest, so look for
				// an intermediate lower point, i.e., do a one-dimensional contraction.
				ysave = y[ihi];
				ytry = Simplex_Point_Try(p, y, psum, ihi, (tfloat)0.5, funk, funkdata);
				if (ytry >= ysave) {
					// Can't seem to get rid of that high point. Better contract
					// around the lowest (best) point.
					for (i=0; i<=ndim; i++)
					{
						if (i != ilo)
						{
							for (j=0; j<ndim; j++)
								p[i][j] = psum[j] = 0.5 * (p[i][j] + p[ilo][j]);
							y[i] = funk(psum, funkdata);
						}
					}
					// Keep track of function evaluations.
					nfunk += ndim;
					// Recompute psum
					for (j=0; j<ndim; j++)
					{
						for (sum=0, i=0; i<=ndim; i++) sum += p[i][j];
						psum[j] = sum;
					}
				}
			} else
				// Correct the evaluation count.
				-- nfunk;
		}
	}
}


namespace IBS
{
	using namespace CoreArray;

	/// The number of IBS 0 in the packed genotype
	extern C_UInt8 IBS0_Num_SNP[];
	/// The number of IBS 1 in the packed genotype
	extern C_UInt8 IBS1_Num_SNP[];
	/// The number of IBS 2 in the packed genotype
	extern C_UInt8 IBS2_Num_SNP[];
}


namespace IBD
{
	// using namespace
	using namespace std;
	using namespace CoreArray;
	using namespace CoreArray::Vectorization;
	using namespace GWAS;

	// ---------------------------------------------------------------------
	// public parameters

	/// The structure of 3 IBD coefficients
	struct TIBD
	{
		double k0, k1;
		TIBD() { k0 = k1 = 0; }
	};

	/// The structure of 9 IBD coefficients
	struct TIBD_Jacq
	{
		double D1, D2, D3, D4, D5, D6, D7, D8;
		TIBD_Jacq()
		{
			D1 = D2 = D3 = D4 = D5 = D6 = D7 = D8 = 0;
		}
		TIBD_Jacq(double f1, double f2, double f3, double f4, double f5,
			double f6, double f7, double f8)
		{
			D1 = f1; D2 = f2; D3 = f3; D4 = f4;
			D5 = f5; D6 = f6; D7 = f7; D8 = f8;
		}
	};


	/// Whether constrict the estimates from PLINK
	bool KinshipConstraint = false;


	// ---------------------------------------------------------------------
	// IBD: PLINK Method of Moment

	/// The expected probability of IBS i, given IBD j is EPrIBS_IBD[i][j]
	double EPrIBS_IBD[3][3];

	/// Initialzie the expected probability of IBS i, given by IBD
	// AlleleFreq should have enough buffer
	void Init_EPrIBD_IBS(const double in_afreq[], double out_afreq[],
		bool CorrectFactor, long nSNP = -1)
	{
		if (nSNP < 0)
			nSNP = MCWorkingGeno.Space.SNPNum();

		// clear EPrIBS_IBD
		memset((void*)EPrIBS_IBD, 0, sizeof(EPrIBS_IBD));
		vector<int> AA(nSNP), AB(nSNP), BB(nSNP);

		if (!in_afreq)
		{
			MCWorkingGeno.Space.GetABNumPerSNP(&AA[0], &AB[0], &BB[0]);
		}

		// for-loop each snp
		long nValid = 0;
		for (long i=0; i < nSNP; i++)
		{
			long n = 2 * (AA[i] + AB[i] + BB[i]);
			double p = (n > 0) ? (double(2*AA[i] + AB[i]) / n) : R_NaN;

			if (in_afreq)
			{
				p = in_afreq[i];
				if (R_FINITE(p))
				{
					if ((p < 0) || (p > 1))
						p = R_NaN;
				}
			}
			
			if (out_afreq) out_afreq[i] = p;

			// Second, the expected probability of IBS i, given by IBD
			double q = 1-p, Na = n;
			double x = 2*AA[i] + AB[i], y = 2*BB[i] + AB[i];
			double a00, a01, a02, a11, a12;

			if (CorrectFactor)
			{
				a00 =
					2*p*p*q*q * ( (x-1)/x * (y-1)/y *
					(Na/(Na-1)) *(Na/(Na-2)) * (Na/(Na-3)) );
				a01 =
					4*p*p*p*q * ( (x-1)/x * (x-2)/x * (Na/(Na-1)) *
					(Na/(Na-2)) * (Na/(Na-3)) ) + 4*p*q*q*q * ( (y-1)/y * (y-2)/y *
					(Na/(Na-1)) * (Na/(Na-2)) * (Na/(Na-3)) );
				a02 =
					q*q*q*q * ( (y-1)/y * (y-2)/y * (y-3)/y *
					(Na/(Na-1)) * (Na/(Na-2)) * (Na/(Na-3)) ) + p*p*p*p *
					( (x-1)/x * (x-2)/x * (x-3)/x * (Na/(Na-1)) * (Na/(Na-2)) *
					(Na/(Na-3)) ) + 4*p*p*q*q * ( (x-1)/x * (y-1)/y * (Na/(Na-1)) *
					(Na/(Na-2)) * (Na/(Na-3)) );
				a11 =
					2*p*p*q * ( (x-1)/x *  Na/(Na-1) * Na/(Na-2) ) +
					2*p*q*q * ( (y-1)/y *  Na/(Na-1) * Na/(Na-2) );
				a12 =
					p*p*p * ((x-1)/x * (x-2)/x *  Na/(Na-1) * Na/(Na-2)) +
					q*q*q * ( (y-1)/y * (y-2)/y *  Na/(Na-1) * Na/(Na-2)) +
					p*p*q * ( (x-1)/x * Na/(Na-1) * Na/(Na-2) ) +
					p*q*q * ((y-1)/y  * Na/(Na-1) * Na/(Na-2));
			} else {
				a00 = 2*p*p*q*q;
				a01 = 4*p*p*p*q + 4*p*q*q*q;
				a02 = q*q*q*q + p*p*p*p + 4*p*p*q*q;
				a11 = 2*p*p*q + 2*p*q*q;
				a12 = p*p*p + q*q*q + p*p*q + p*q*q;
			}

			if (R_FINITE(a00) && R_FINITE(a01) &&
				R_FINITE(a02) && R_FINITE(a11) && R_FINITE(a12))
			{
				EPrIBS_IBD[0][0] += a00;
				EPrIBS_IBD[0][1] += a01;
				EPrIBS_IBD[0][2] += a02;
				EPrIBS_IBD[1][1] += a11;
				EPrIBS_IBD[1][2] += a12;
				nValid++;
			}
		}

		EPrIBS_IBD[0][0] /= nValid; EPrIBS_IBD[1][0] = 0;       EPrIBS_IBD[2][0] = 0;
		EPrIBS_IBD[0][1] /= nValid; EPrIBS_IBD[1][1] /= nValid; EPrIBS_IBD[2][1] = 0;
		EPrIBS_IBD[0][2] /= nValid; EPrIBS_IBD[1][2] /= nValid; EPrIBS_IBD[2][2] = 1;
	}

	// estimate k0, k1 and k2
	void Est_PLINK_Kinship(int IBS0, int IBS1, int IBS2, double &k0, double &k1,
		bool constraint)
	{
		int nIBS012 = IBS0 + IBS1 + IBS2;
		double e00 = EPrIBS_IBD[0][0] * nIBS012;
		double e01 = EPrIBS_IBD[0][1] * nIBS012;
		double e11 = EPrIBS_IBD[1][1] * nIBS012;
		double e02 = EPrIBS_IBD[0][2] * nIBS012;
		double e12 = EPrIBS_IBD[1][2] * nIBS012;
		double e22 = EPrIBS_IBD[2][2] * nIBS012;

		k0 = IBS0 / e00;
		k1 = (IBS1 - k0 * e01) / e11;
		double k2 = (IBS2 - k0*e02 - k1*e12) / e22;

		// Bound IBD estimates to sum to 1, and fall within 0-1 range
		if (k0 > 1) { k0 = 1; k1 = k2 = 0; }
		if (k1 > 1) { k1 = 1; k0 = k2 = 0; }
		if (k2 > 1) { k2 = 1; k0 = k1 = 0; }
		if (k0 < 0) { double S = k1+k2; k1 /= S; k2 /= S; k0 = 0; }
		if (k1 < 0) { double S = k0+k2; k0 /= S; k2 /= S; k1 = 0; }
		if (k2 < 0) { double S = k0+k1; k0 /= S; k1 /= S; k2 = 0; }

		if (constraint)
		{
			// Possibly constrain IBD estimates to within possible triangle
			// i.e. 0.5 0.0 0.5 is invalid
			//
			// Constraint : z1^2 - 4 z0 z2 >= 0
			//            : x^2 - 2 pi x + z2  = 0
			//              where pi = (z1 + 2 z2) / 2
			//
			// So the constaint can also be written as
			//              pi^2 >=  z2
			k2 = 1 - k0 - k1;
			double pihat = k1 / 2 + k2 ;
			if (pihat*pihat < k2)
			{
				k0 = (1-pihat) * (1-pihat);
				k1 = 2 * pihat * (1-pihat);
			}
		}
	}



	// ---------------------------------------------------------------------
	// IBD: Maximum Likelihood Estimation (MLE)

	static const double IBDMLE_InitVal_Tol = 0.005;

	// public parameters
	/// The maximum number of iterations, 1000 by default
	long nIterMax = 1000;
	/// The reltol convergence tolerance, sqrt(machine.epsilon) by default
	double FuncRelTol = sqrt(DBL_EPSILON);
	/// the computational method for MLE, 0 -- EM, 1 -- NM
	int MethodMLE = 0;
	/// whether or not adjust loglik
	bool Loglik_Adjust = true;


	// internal parameters
	CdMatTriDiag<TIBD> *IBD = NULL;
	TIBD *pMatIBD = NULL;
	CdMatTriDiag<TIBD_Jacq> *IBD_Jacq = NULL;
	TIBD_Jacq *pMatIBD_Jacq = NULL;

	int *pNIter = NULL;
	IdMatTriD IBD_idx;


	/// Genotype, stored in a packed mode
	C_UInt8 *PackedGeno = NULL;
	/// the number of samples and snps
	long nSamp, nPackedSNP, nTotalSNP;
	///
	long nMatTriD, idxMatTriD;

	// AlleleFreq
	double *MLEAlleleFreq = NULL;

	/// initialize the variable "Geno" with packed genotypes
	void InitPackedGeno(void *buffer)
	{
		// set # of samples and snps
		nSamp = MCWorkingGeno.Space.SampleNum();
		nPackedSNP = (MCWorkingGeno.Space.SNPNum() % 4 > 0) ?
			(MCWorkingGeno.Space.SNPNum()/4 + 1) : (MCWorkingGeno.Space.SNPNum()/4);
		nTotalSNP = nPackedSNP * 4;
		PackedGeno = (C_UInt8*)buffer;

		// buffer
		CdBufSpace Buf(MCWorkingGeno.Space, false, CdBufSpace::acInc);
		C_UInt8 *p = PackedGeno;
		for (long i=0; i < MCWorkingGeno.Space.SampleNum(); i++)
		{
			p = Buf.ReadPackedGeno(i, p);
		}
	}



	// ================================================================ //

	/// Output the probability of IBS given by IBD
	/** \param g1  the first genotype (0 -- BB, 1 -- AB, 2 -- AA, other missing)
	 *  \param g2  the second genotype
	 *  \param t0  the probability of IBD given by IBD state = 0
	 *  \param t1  the probability of IBD given by IBD state = 1
	 *  \param t2  the probability of IBD given by IBD state = 2
	 *  \param p   the frequency of allele A
	**/
	void PrIBDTable(int g1, int g2,
		double &t0, double &t1, double &t2, const double p)
	{
		if ((0 < p) && (p < 1))
		{
			const double q = 1 - p;
			switch (g1) {
				case 0: /* mm */
					switch (g2) {
						case 0: /* mm,mm */
							t2 = q*q; t1 = t2*q; t0 = t1*q; break;
						case 1: /* mm, Mm */
							t1 = p*q*q; t0 = 2*t1*q; t2 = 0; break;
						case 2: /* mm, MM */
							t0 = p*p*q*q; t1 = t2 = 0; break;
						default:
							t0 = t1 = t2 = 0;
					}
					break;
				case 1: /* Mm */
					switch (g2) {
						case 0: /* Mm, mm */
							t1 = p*q*q; t0 = 2*t1*q; t2 = 0; break;
						case 1: /* Mm, Mm */
							t1 = p*q; t0 = 4*t1*t1; t2 = 2*t1; break;
						case 2: /* Mm, MM */
							t1 = p*p*q; t0 = 2*p*t1; t2 = 0; break;
						default:
							t0 = t1 = t2 = 0;
					}
					break;
				case 2: /* MM */
					switch (g2) {
						case 0: /* MM, mm */
							t0 = p*p*q*q; t1 = t2 = 0; break;
						case 1: /* MM, Mm */
							t1 = p*p*q; t0 = 2*p*t1; t2 = 0; break;
						case 2: /* MM, MM */
							t2 = p*p; t1 = t2*p; t0 = t1*p; break;
						default:
							t0 = t1 = t2 = 0;
					}
					break;
				default:
					t0 = t1 = t2 = 0;
			}
		} else
			t0 = t1 = t2 = 0;
	}


	#define LOGLIK_ADJUST(FUNC, _k0, _k1)	\
		_loglik = FUNC(PrIBD, _k0, _k1); \
		if (R_FINITE(_loglik)) { \
			if (out_loglik < _loglik) {	\
				out_loglik = _loglik; out_k0 = _k0; out_k1 = _k1; \
    		}	\
		}


	// ==================  MLE - EM algorithm  ==================

	static void EM_Prepare(double *PrIBD, C_UInt8 *p1, C_UInt8 *p2)
	{
		const double *Freq = &MLEAlleleFreq[0];
		for (long nPackSNP=nPackedSNP; nPackSNP > 0; nPackSNP--)
		{
			// first genotype
			C_UInt8 g1=*p1++, g2=*p2++;
			PrIBDTable(g1 & 0x03, g2 & 0x03, PrIBD[0], PrIBD[1], PrIBD[2], *Freq++);
			// second genotype
			g1 >>= 2; g2 >>= 2; PrIBD += 3;
			PrIBDTable(g1 & 0x03, g2 & 0x03, PrIBD[0], PrIBD[1], PrIBD[2], *Freq++);
			// third genotype
			g1 >>= 2; g2 >>= 2; PrIBD += 3;
			PrIBDTable(g1 & 0x03, g2 & 0x03, PrIBD[0], PrIBD[1], PrIBD[2], *Freq++);
			// fourth genotype
			g1 >>= 2; g2 >>= 2; PrIBD += 3;
			PrIBDTable(g1 & 0x03, g2 & 0x03, PrIBD[0], PrIBD[1], PrIBD[2], *Freq++);
			PrIBD += 3;
		}
	}

	/// return log likelihood value for the specified pair
	static double EM_LogLik(const double *PrIBD, const double k0, const double k1)
	{
		double k[3] = { k0, k1, 1-k0-k1 }, sum;
		const double *pr = PrIBD;
		double LogLik = 0;
		// unroll loop by 4
		for (long nPackSNP=nPackedSNP; nPackSNP > 0; nPackSNP--)
		{
			// first genotype
			sum = pr[0]*k[0] + pr[1]*k[1] + pr[2]*k[2];
			if (sum > 0)
				LogLik += log(sum);
			else if (pr[0] > 0)
				return R_NegInf;
			pr += 3;
			// second genotype
			sum = pr[0]*k[0] + pr[1]*k[1] + pr[2]*k[2];
			if (sum > 0)
				LogLik += log(sum);
			else if (pr[0] > 0)
				return R_NegInf;
			pr += 3;
			// third genotype
			sum = pr[0]*k[0] + pr[1]*k[1] + pr[2]*k[2];
			if (sum > 0)
				LogLik += log(sum);
			else if (pr[0] > 0)
				return R_NegInf;
			pr += 3;
			// fourth genotype
			sum = pr[0]*k[0] + pr[1]*k[1] + pr[2]*k[2];
			if (sum > 0)
				LogLik += log(sum);
			else if (pr[0] > 0)
				return R_NegInf;
			pr += 3;
		}
		return LogLik;
	}

	/// MLE -- EM algorithm
	/**   out_k0 and out_k1 are the initial parameter values, and they
	 *  should be within the region of { 0 < k0, 0 < k1, k0+k1 < 1 }.
	**/
	static void EMAlg(double *PrIBD, double &out_k0, double &out_k1, double &out_loglik,
		int *out_niter)
	{
		double k[3] = { out_k0, out_k1, 1-out_k0-out_k1 };
		// the old value of - log likelihood function
		double OldLogLik = 0;
		// the current value of - log likelihood function
		double LogLik = EM_LogLik(PrIBD, k[0], k[1]);
		// convergence tolerance
		double ConvTol;

		// check log likelihood
		if (R_FINITE(LogLik))
		{
			// the converage tolerance
			ConvTol = FuncRelTol * (fabs(LogLik) + fabs(FuncRelTol));
			if (ConvTol < 0) ConvTol = 0;
		} else {
			LogLik = 1e+8;
			ConvTol = FuncRelTol;
		}

		if (out_niter) *out_niter = nIterMax;

		// EM iteration
		for (long iIter=0; iIter <= nIterMax; iIter++)
		{
			double oldk[3] = { k[0], k[1], k[2] };
			double sum[2] = { 0, 0 };
			const double *pr = PrIBD;
			long nSNP = 0;

			LogLik = 0;
			for (long iSNP = nTotalSNP; iSNP > 0; iSNP--)
			{
				double mul[3] = { pr[0]*k[0], pr[1]*k[1], pr[2]*k[2] };
				double mulsum = mul[0] + mul[1] + mul[2];
				if (mulsum > 0)
				{
					sum[0] += mul[0] / mulsum; sum[1] += mul[1] / mulsum;
					nSNP ++;
					LogLik += log(mulsum);
				} else if (pr[0] > 0)
					// this will never happen
					throw "Invalid updated IBD coefficient parameters.";
				pr += 3;
			}
			// update k0, k1, k2 values
			k[0] = sum[0] / nSNP; k[1] = sum[1] / nSNP;
			k[2] = 1 - k[0] - k[1];

			// stopping rule
			if (fabs(LogLik - OldLogLik) <= ConvTol)
			{
				k[0] = oldk[0]; k[1] = oldk[1]; k[2] = oldk[2];
				if (out_niter) *out_niter = iIter;
				break;
			}

			OldLogLik = LogLik;
		}

		// fill the result
		out_k0 = k[0]; out_k1 = k[1]; out_loglik = LogLik;
		if (Loglik_Adjust)
		{
			double _loglik;
			LOGLIK_ADJUST(EM_LogLik, 0, 0);  // self
			LOGLIK_ADJUST(EM_LogLik, 0.25, 0.5);  // full sibs
			LOGLIK_ADJUST(EM_LogLik, 0, 1);  // offspring
			LOGLIK_ADJUST(EM_LogLik, 0.5, 0.5);  // half sibs
			LOGLIK_ADJUST(EM_LogLik, 0.75, 0.25);  // cousins
			LOGLIK_ADJUST(EM_LogLik, 1, 0);  // unrelated
		}
	}


	// ================== MLE - downhill simplex algorithm ==================

	static void NM_Prepare(double *pr, C_UInt8 *p1, C_UInt8 *p2)
	{
		const double *Freq = &MLEAlleleFreq[0];
		for (long nPackSNP=nPackedSNP; nPackSNP > 0; nPackSNP--)
		{
			// first genotype
			C_UInt8 g1=*p1++, g2=*p2++;
			PrIBDTable(g1 & 0x03, g2 & 0x03, pr[0], pr[1], pr[2], *Freq++);
			pr[0] -= pr[2]; pr[1] -= pr[2];
			// second genotype
			g1 >>= 2; g2 >>= 2; pr += 3;
			PrIBDTable(g1 & 0x03, g2 & 0x03, pr[0], pr[1], pr[2], *Freq++);
			pr[0] -= pr[2]; pr[1] -= pr[2];
			// third genotype
			g1 >>= 2; g2 >>= 2; pr += 3;
			PrIBDTable(g1 & 0x03, g2 & 0x03, pr[0], pr[1], pr[2], *Freq++);
			pr[0] -= pr[2]; pr[1] -= pr[2];
			// fourth genotype
			g1 >>= 2; g2 >>= 2; pr += 3;
			PrIBDTable(g1 & 0x03, g2 & 0x03, pr[0], pr[1], pr[2], *Freq++);
			pr[0] -= pr[2]; pr[1] -= pr[2];
			pr += 3;
		}
	}

	/// return log likelihood value for the specified pair
	static double NM_LogLik(const double *PrIBD, const double k0, const double k1)
	{
		// Check whether within the region
		if ((k0<0) || (k1<0) || (k0+k1>1)) return R_NegInf;

		const double *pr = PrIBD;
		double sum, LogLik = 0;
		// unroll loop by 4
		for (long nPackSNP=nPackedSNP; nPackSNP > 0; nPackSNP--)
		{
			// first genotype
			sum = pr[0]*k0 + pr[1]*k1 + pr[2];
			if (sum > 0)
				LogLik += log(sum);
			else if (pr[0] > 0)
				return R_NegInf;
			pr += 3;
			// second genotype
			sum = pr[0]*k0 + pr[1]*k1 + pr[2];
			if (sum > 0)
				LogLik += log(sum);
			else if (pr[0] > 0)
				return R_NegInf;
			pr += 3;
			// third genotype
			sum = pr[0]*k0 + pr[1]*k1 + pr[2];
			if (sum > 0)
				LogLik += log(sum);
			else if (pr[0] > 0)
				return R_NegInf;
			pr += 3;
			// fourth genotype
			sum = pr[0]*k0 + pr[1]*k1 + pr[2];
			if (sum > 0)
				LogLik += log(sum);
			else if (pr[0] > 0)
				return R_NegInf;
			pr += 3;
		}
		return LogLik;
	}

	/// the optimize function
	static double _optim(const double *x, void *ex)
	{
		double rv = -NM_LogLik((double*)ex, x[0], x[1]);
		if (!R_FINITE(rv)) rv = 1e+30;
		return rv;
	}

	/// MLE - downhill simplex algorithm
	/**   Out.k0 and Out.k1 are the initial parameter values, and they
	 *  should be within the region of { 0 < k0, 0 < k1, k0+k1 < 1 }.
	**/
	static void Simplex(double *PrIBD, double &out_k0, double &out_k1, double &out_loglik,
		int *out_niter)
	{
		// Determine three points in the simplex
		double p[3][2], f;
		// the first point
		p[0][0] = out_k0; p[0][1] = out_k1;
		// the second point
		p[1][0] = out_k0; f = (1 - out_k0) / 2;
		p[1][1] = (out_k1 <= f) ?
			(out_k1 + max(out_k1,   f-out_k1) / 2) :
			(out_k1 - max(out_k1-f, 1-out_k0-out_k1));
		// the third point
		p[2][1] = out_k1; f = (1 - out_k1) / 2;
		p[2][0] = (out_k0 <= f) ?
			(out_k0 + max(out_k0,   f-out_k0) / 2) :
			(out_k0 - max(out_k0-f, 1-out_k1-out_k0) / 2);

		double k[2], LogLik;
		int nIter;

		// Downhill simplex approach
		GWAS_Math::SimplexMin<double>(p, k, LogLik, nIter, _optim,
			(void*)PrIBD, FuncRelTol, nIterMax);
		if (out_niter) *out_niter = nIter;

		// fill the result
		out_k0 = k[0]; out_k1 = k[1]; out_loglik = -LogLik;
		if (Loglik_Adjust)
		{
			double _loglik;
			LOGLIK_ADJUST(NM_LogLik, 0, 0);  // self
			LOGLIK_ADJUST(NM_LogLik, 0.25, 0.5);  // full sibs
			LOGLIK_ADJUST(NM_LogLik, 0, 1);  // offspring
			LOGLIK_ADJUST(NM_LogLik, 0.5, 0.5);  // half sibs
			LOGLIK_ADJUST(NM_LogLik, 0.75, 0.25);  // cousins
			LOGLIK_ADJUST(NM_LogLik, 1, 0);  // unrelated
		}
	}


	/// The thread entry for the calculation of genetic covariace matrix
	static void Entry_MLEIBD(PdThread Thread, int ThreadIndex, void *Param)
	{
		// initialize buffer
		vector<double> PrIBD(3*nTotalSNP);

		// loop
		while (true)
		{
			// check
			bool WorkFlag;
			IdMatTriD idx(0);
			TIBD *pIBD = NULL;
			int *pniter = NULL;
			{
				TdAutoMutex _m(_Mutex);
				WorkFlag = idxMatTriD < nMatTriD;
				if (WorkFlag)
				{
					idx = IBD_idx; ++IBD_idx; idxMatTriD ++;
					pIBD = pMatIBD; pMatIBD ++;
					if (pNIter)
						{ pniter = pNIter; pNIter ++; }
					MCWorkingGeno.Progress.Forward(1, Thread==0);
				}
			}
			if (!WorkFlag) break;

			C_UInt8 *g1 = &PackedGeno[0] + nPackedSNP*idx.Row();
			C_UInt8 *g2 = &PackedGeno[0] + nPackedSNP*idx.Column();

			// calculate the initial values from PLINK
			C_UInt8 *p1 = g1, *p2 = g2;
			long IBS0=0, IBS1=0, IBS2=0;
			for (long i=0; i < nPackedSNP; i++, p1++, p2++)
			{
				size_t t = (size_t(*p1) << 8) | (*p2);
				IBS0 += IBS::IBS0_Num_SNP[t];
				IBS1 += IBS::IBS1_Num_SNP[t];
				IBS2 += IBS::IBS2_Num_SNP[t];
			}
			Est_PLINK_Kinship(IBS0, IBS1, IBS2, pIBD->k0, pIBD->k1, false);
			// adjust the initial values
			{
				double k0 = pIBD->k0, k1 = pIBD->k1, k2 = 1-k0-k1;
				if (k0 < IBDMLE_InitVal_Tol) k0 = IBDMLE_InitVal_Tol;
				if (k1 < IBDMLE_InitVal_Tol) k1 = IBDMLE_InitVal_Tol;
				if (k2 < IBDMLE_InitVal_Tol) k2 = IBDMLE_InitVal_Tol;
				double s = k0 + k1 + k2;
				pIBD->k0 = k0/s; pIBD->k1 = k1/s;
			}

			// MLE
			double _loglik;
			switch (MethodMLE)
			{
				case 0:
					// fill PrIBD, Pr(ibs state|ibd state j)
					EM_Prepare(&PrIBD[0], g1, g2);
					// MLE - EM algorithm
					EMAlg(&PrIBD[0], pIBD->k0, pIBD->k1, _loglik, pniter);
					break;

				case 1:
					// fill PrIBD, Pr(ibs state|ibd state j)
					NM_Prepare(&PrIBD[0], g1, g2);
					// MLE - downhill simplex algorithm
					Simplex(&PrIBD[0], pIBD->k0, pIBD->k1, _loglik, pniter);
					break;
			}
		}
	}


	// =============== MLE - EM algorithm for Jacquard's IBD ===============

	/// Output the probability of IBS given by IBD
	/** \param g1  the first genotype (0 -- BB, 1 -- AB, 2 -- AA, other missing)
	 *  \param g2  the second genotype
	 *  \param Pr  the probability of IBD given by IBD state
	 *  \param p   the frequency of allele A
	**/
	static void PrIBDTabJacq(int g1, int g2, double Pr[], const double p)
	{
		if ((0 < p) && (p < 1))
		{
			const double q = 1 - p;
			switch (g1) {
				case 0: /* mm */
					switch (g2) {
						case 0: /* mm,mm */
							Pr[0] = q;
							Pr[1] = Pr[2] = Pr[4] = Pr[6] = q*q;
							Pr[3] = Pr[5] = Pr[7] = q*q*q;
							Pr[8] = q*q*q*q;
							break;
						case 1: /* mm, Mm */
							Pr[0] = Pr[1] = Pr[4] = Pr[5] = Pr[6] = 0;
							Pr[2] = p*q; Pr[3] = 2*p*q*q;
							Pr[7] = p*q*q; Pr[8] = 2*p*q*q*q;
							break;
						case 2: /* mm, MM */
							Pr[0] = Pr[2] = Pr[4] = Pr[6] = Pr[7] = 0;
							Pr[1] = p*q; Pr[3] = p*p*q;
							Pr[5] = p*q*q; Pr[8] = p*p*q*q;
							break;
						default:
							Pr[0]=Pr[1]=Pr[2]=Pr[3]=Pr[4]=Pr[5]=Pr[6]=Pr[7]=Pr[8]=0;
					}
					break;
				case 1: /* Mm */
					switch (g2) {
						case 0: /* Mm, mm */
							Pr[0] = Pr[1] = Pr[2] = Pr[3] = Pr[6] = 0;
							Pr[4] = p*q; Pr[5] = 2*p*q*q;
							Pr[7] = p*q*q; Pr[8] = 2*p*q*q*q;
							break;
						case 1: /* Mm, Mm */
							Pr[0] = Pr[1] = Pr[2] = Pr[3] = Pr[4] = Pr[5] = 0;
							Pr[6] = 2*p*q; Pr[7] = p*q;
							Pr[8] = 4*p*p*q*q;
							break;
						case 2: /* Mm, MM */
							Pr[0] = Pr[1] = Pr[2] = Pr[3] = Pr[6] = 0;
							Pr[4] = p*q; Pr[5] = 2*p*p*q;
							Pr[7] = p*p*q; Pr[8] = 2*p*p*p*q;
							break;
						default:
							Pr[0]=Pr[1]=Pr[2]=Pr[3]=Pr[4]=Pr[5]=Pr[6]=Pr[7]=Pr[8]=0;
					}
					break;
				case 2: /* MM */
					switch (g2) {
						case 0: /* MM, mm */
							Pr[0] = Pr[2] = Pr[4] = Pr[6] = Pr[7] = 0;
							Pr[1] = p*q; Pr[3] = p*q*q;
							Pr[5] = p*p*q; Pr[8] = p*p*q*q;
							break;
						case 1: /* MM, Mm */
							Pr[0] = Pr[1] = Pr[4] = Pr[5] = Pr[6] = 0;
							Pr[2] = p*q; Pr[3] = 2*p*p*q;
							Pr[7] = p*p*q; Pr[8] = 2*p*p*p*q;
							break;
						case 2: /* MM, MM */
							Pr[0] = p;
							Pr[1] = Pr[2] = Pr[4] = Pr[6] = p*p;
							Pr[3] = Pr[5] = Pr[7] = p*p*p;
							Pr[8] = p*p*p*p;
						default:
							Pr[0]=Pr[1]=Pr[2]=Pr[3]=Pr[4]=Pr[5]=Pr[6]=Pr[7]=Pr[8]=0;
					}
					break;
				default:
					Pr[0]=Pr[1]=Pr[2]=Pr[3]=Pr[4]=Pr[5]=Pr[6]=Pr[7]=Pr[8]=0;
			}
		} else
			Pr[0]=Pr[1]=Pr[2]=Pr[3]=Pr[4]=Pr[5]=Pr[6]=Pr[7]=Pr[8]=0;
	}


	static void EM_Jacq_Prepare(double *PrIBD, C_UInt8 *p1, C_UInt8 *p2)
	{
		const double *Freq = &MLEAlleleFreq[0];
		for (long nPackSNP=nPackedSNP; nPackSNP > 0; nPackSNP--)
		{
			// first genotype
			C_UInt8 g1=*p1++, g2=*p2++;
			PrIBDTabJacq(g1 & 0x03, g2 & 0x03, PrIBD, *Freq++);
			// second genotype
			g1 >>= 2; g2 >>= 2; PrIBD += 9;
			PrIBDTabJacq(g1 & 0x03, g2 & 0x03, PrIBD, *Freq++);
			// third genotype
			g1 >>= 2; g2 >>= 2; PrIBD += 9;
			PrIBDTabJacq(g1 & 0x03, g2 & 0x03, PrIBD, *Freq++);
			// fourth genotype
			g1 >>= 2; g2 >>= 2; PrIBD += 9;
			PrIBDTabJacq(g1 & 0x03, g2 & 0x03, PrIBD, *Freq++);
			PrIBD += 9;
		}
	}

	/// return log likelihood value for the specified pair
	static double EM_Jacq_LogLik(const double *PrIBD,
		const TIBD_Jacq &par)
	{
		double D9 = 1 - par.D1 - par.D2 - par.D3 - par.D4 - par.D5 -
				par.D6 - par.D7 - par.D8;
		const double *p = PrIBD;
		double sum, LogLik=0;

		// unroll loop by 4
		for (long nPackSNP=nPackedSNP; nPackSNP > 0; nPackSNP--)
		{
			// first genotype
			sum = p[0]*par.D1 + p[1]*par.D2 + p[2]*par.D3 + p[3]*par.D4 +
				p[4]*par.D5 + p[5]*par.D6 + p[6]*par.D7 + p[7]*par.D8 + p[8]*D9;
			if (sum > 0)
				LogLik += log(sum);
			else if (p[8] > 0)
				return R_NegInf;
			p += 9;
			// second genotype
			sum = p[0]*par.D1 + p[1]*par.D2 + p[2]*par.D3 + p[3]*par.D4 +
				p[4]*par.D5 + p[5]*par.D6 + p[6]*par.D7 + p[7]*par.D8 + p[8]*D9;
			if (sum > 0)
				LogLik += log(sum);
			else if (p[8] > 0)
				return R_NegInf;
			p += 9;
			// third genotype
			sum = p[0]*par.D1 + p[1]*par.D2 + p[2]*par.D3 + p[3]*par.D4 +
				p[4]*par.D5 + p[5]*par.D6 + p[6]*par.D7 + p[7]*par.D8 + p[8]*D9;
			if (sum > 0)
				LogLik += log(sum);
			else if (p[8] > 0)
				return R_NegInf;
			p += 9;
			// fourth genotype
			sum = p[0]*par.D1 + p[1]*par.D2 + p[2]*par.D3 + p[3]*par.D4 +
				p[4]*par.D5 + p[5]*par.D6 + p[6]*par.D7 + p[7]*par.D8 + p[8]*D9;
			if (sum > 0)
				LogLik += log(sum);
			else if (p[8] > 0)
				return R_NegInf;
			p += 9;
		}
		return LogLik;
	}

	/// MLE -- EM algorithm for Jacquard's IBD
	/**   the initial parameter values, and they should be within
	 *  the region.
	**/
	static void EM_Jacq_Alg(double *PrIBD, TIBD_Jacq &par,
		double &out_loglik, int *out_niter)
	{
		double D[9] = { par.D1, par.D2, par.D3, par.D4,
			par.D5, par.D6, par.D7, par.D8,
			1 - par.D1 - par.D2 - par.D3 - par.D4 - par.D5 -
				par.D6 - par.D7 - par.D8
		};

		// the old value of - log likelihood function
		double OldLogLik = 0;
		// the current value of - log likelihood function
		double LogLik = EM_Jacq_LogLik(PrIBD, par);
		// convergence tolerance
		double ConvTol;

		// check log likelihood
		if (R_FINITE(LogLik))
		{
			// the converage tolerance
			ConvTol = FuncRelTol * (fabs(LogLik) + fabs(FuncRelTol));
			if (ConvTol < 0) ConvTol = 0;
		} else {
			LogLik = 1e+8;
			ConvTol = FuncRelTol;
		}

		if (out_niter) *out_niter = nIterMax;

		// EM iteration
		for (long iIter=0; iIter <= nIterMax; iIter++)
		{
			double OldD[9], sum[9];
			memcpy(OldD, D, sizeof(D));
			memset(sum, 0, sizeof(sum));

			const double *p = PrIBD;
			long nSNP = 0;

			LogLik = 0;
			for (long iSNP = nTotalSNP; iSNP > 0; iSNP--)
			{
				double m[9] = { p[0]*D[0], p[1]*D[1], p[2]*D[2],
					p[3]*D[3], p[4]*D[4], p[5]*D[5],
					p[6]*D[6], p[7]*D[7], p[8]*D[8] };
				double mulsum = m[0] + m[1] + m[2] + m[3] + m[4] +
					m[5] + m[6] + m[7] + m[8];
				if (mulsum > 0)
				{
					sum[0] += m[0] / mulsum; sum[1] += m[1] / mulsum;
					sum[2] += m[2] / mulsum; sum[3] += m[3] / mulsum;
					sum[4] += m[4] / mulsum; sum[5] += m[5] / mulsum;
					sum[6] += m[6] / mulsum; sum[7] += m[7] / mulsum;
					sum[8] += m[8] / mulsum;
					nSNP ++;
					LogLik += log(mulsum);
				} else if (p[8] > 0)
					// this will never happen
					throw "Invalid updated IBD coefficient parameters.";
				p += 9;
			}

			// update D values
			D[0] = sum[0] / nSNP; D[1] = sum[1] / nSNP;
			D[2] = sum[2] / nSNP; D[3] = sum[3] / nSNP;
			D[4] = sum[4] / nSNP; D[5] = sum[5] / nSNP;
			D[6] = sum[6] / nSNP; D[7] = sum[7] / nSNP;
			D[8] = sum[8] / nSNP;

			// stopping rule
			if (fabs(LogLik - OldLogLik) <= ConvTol)
			{
				memcpy(D, OldD, sizeof(D));
				if (out_niter) *out_niter = iIter;
				break;
			}

			OldLogLik = LogLik;
		}

		// fill the result
		par.D1 = D[0]; par.D2 = D[1]; par.D3 = D[2];
		par.D4 = D[3]; par.D5 = D[4]; par.D6 = D[5];
		par.D7 = D[6]; par.D8 = D[7];
		out_loglik = LogLik;
	}

	static TIBD_Jacq IBD_Jacq_InitVal(
		0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01);

	/// The thread entry for the calculation of genetic covariace matrix
	static void Entry_MLEIBD_Jacq(PdThread Thread, int ThreadIndex, void *Param)
	{
		// initialize buffer
		vector<double> PrIBD(9*nTotalSNP);

		// loop
		while (true)
		{
			// check
			bool WorkFlag;
			IdMatTriD idx(0);
			TIBD_Jacq *pIBD = NULL;
			int *pniter = NULL;
			{
				TdAutoMutex _m(_Mutex);
				WorkFlag = idxMatTriD < nMatTriD;
				if (WorkFlag)
				{
					idx = IBD_idx; ++IBD_idx; idxMatTriD ++;
					pIBD = pMatIBD_Jacq; pMatIBD_Jacq ++;
					if (pNIter)
						{ pniter = pNIter; pNIter ++; }
					MCWorkingGeno.Progress.Forward(1, Thread==0);
				}
			}
			if (!WorkFlag) break;

			C_UInt8 *g1 = &PackedGeno[0] + nPackedSNP*idx.Row();
			C_UInt8 *g2 = &PackedGeno[0] + nPackedSNP*idx.Column();

			// the initial values
			*pIBD = IBD_Jacq_InitVal;

			// MLE
			double _loglik;
			// fill PrIBD, Pr(ibs state|ibd state j)
			EM_Jacq_Prepare(&PrIBD[0], g1, g2);
			// MLE - EM algorithm
			EM_Jacq_Alg(&PrIBD[0], *pIBD, _loglik, pniter);
		}
	}


	/// initialize allele frequency
	static void InitAFreq(const double *AFreq, double *AF)
	{
		// initialize the allele frequency
		MLEAlleleFreq = AF;
		for (int i=0; i < nTotalSNP; i++)
			MLEAlleleFreq[i] = -1;
		if (AFreq)
		{
			for (int i=0; i < MCWorkingGeno.Space.SNPNum(); i++)
			{
				if (R_FINITE(AFreq[i]))
					MLEAlleleFreq[i] = AFreq[i];
			}
		} else {
			// init
			vector<int> n(nTotalSNP);
			for (int i=0; i < nTotalSNP; i++) n[i] = 0;
			for (int i=0; i < nTotalSNP; i++) MLEAlleleFreq[i] = 0;
			// for-loop for allele frequency
			C_UInt8 *p = &PackedGeno[0];
			for (int iSamp=0; iSamp < nSamp; iSamp++)
			{
				for (int i=0; i < nPackedSNP; i++)
				{
					C_UInt8 L = *p++;
					for (int k=0; k < 4; k++)
					{
						C_UInt8 B = L & 0x03; L >>= 2;
						if (B < 3)
						{
							n[i*4 + k] += 2;
							MLEAlleleFreq[i*4 + k] += B;
						}
					}
				}
			}
			// for-loop
			for (int i=0; i < nTotalSNP; i++)
			{
				double &v = MLEAlleleFreq[i];
				v = (n[i] > 0) ? (v / n[i]) : -1;
			}
		}
	}


	/// to conduct MLE IBD (3 coefficients)
	void Do_MLE_IBD_Calc(const double *AFreq,
		CdMatTriDiag<TIBD> &PublicIBD, CdMatTriDiag<int> *PublicNIter,
		double *out_AFreq, int NumThread, const char *Info, double *tmpAF,
		bool verbose)
	{
		InitAFreq(AFreq, tmpAF);
		for (int i=0; i < MCWorkingGeno.Space.SNPNum(); i++)
			out_AFreq[i] = MLEAlleleFreq[i];

		IBD = &PublicIBD; pMatIBD = PublicIBD.get();
		pNIter = (PublicNIter) ? PublicNIter->get() : NULL;
		IBD_idx.reset(nSamp);
		nMatTriD = PublicIBD.Size(); idxMatTriD = 0;

		// initialize the mutex object
		_Mutex = GDS_Parallel_InitMutex();

		// Initialize progress information
		MCWorkingGeno.Progress.Info = Info;
		MCWorkingGeno.Progress.Show() = verbose;
		MCWorkingGeno.Progress.Init(nMatTriD);

		// Threads
		GDS_Parallel_RunThreads(Entry_MLEIBD, NULL, NumThread);

		// destory the mutex objects
		GDS_Parallel_DoneMutex(_Mutex); _Mutex = NULL;
	}


	/// to conduct MLE IBD (9 coefficients)
	void Do_MLE_IBD_Jacq(const double *AFreq, CdMatTriDiag<TIBD_Jacq> &PublicIBD,
		CdMatTriDiag<int> *PublicNIter, double *out_AFreq, int NumThread,
		const char *Info, double *tmpAF, bool verbose)
	{
		InitAFreq(AFreq, tmpAF);
		for (int i=0; i < MCWorkingGeno.Space.SNPNum(); i++)
			out_AFreq[i] = MLEAlleleFreq[i];

		IBD_Jacq = &PublicIBD; pMatIBD_Jacq = PublicIBD.get();
		pNIter = (PublicNIter) ? PublicNIter->get() : NULL;
		IBD_idx.reset(nSamp);
		nMatTriD = PublicIBD.Size(); idxMatTriD = 0;

		// initialize the mutex object
		_Mutex = GDS_Parallel_InitMutex();

		// Initialize progress information
		MCWorkingGeno.Progress.Info = Info;
		MCWorkingGeno.Progress.Show() = verbose;
		MCWorkingGeno.Progress.Init(nMatTriD);

		// Threads
		GDS_Parallel_RunThreads(Entry_MLEIBD_Jacq, NULL, NumThread);

		// destory the mutex objects
		GDS_Parallel_DoneMutex(_Mutex); _Mutex = NULL;
	}


	/// to conduct MLE IBD for a pair of individuals, no missing genotypes (0, 1, 2)
	void Do_MLE_IBD_Pair(int n, const int *geno1, const int *geno2,
		const double *AFreq, double &out_k0, double &out_k1,
		double &out_loglik, int &out_niter, double tmpprob[])
	{
		// adjust the initial values
		{
			double k0 = out_k0, k1 = out_k1, k2 = 1 - out_k0 - out_k1;
			if (k0 < IBDMLE_InitVal_Tol) k0 = IBDMLE_InitVal_Tol;
			if (k1 < IBDMLE_InitVal_Tol) k1 = IBDMLE_InitVal_Tol;
			if (k2 < IBDMLE_InitVal_Tol) k2 = IBDMLE_InitVal_Tol;
			double s = k0 + k1 + k2;
			out_k0 = k0/s; out_k1 = k1/s;
		}

		// MLE
		nTotalSNP = n;
		nPackedSNP = n / 4;
		if ((n & 0x03) != 0) nPackedSNP ++;

		double *ptmp = tmpprob;
		switch (MethodMLE)
		{
			case 0:
				// fill PrIBD, Pr(ibs state | ibd state j)
				for (int i=0; i < n; i++)
				{
					PrIBDTable(geno1[i], geno2[i], ptmp[0], ptmp[1], ptmp[2], AFreq[i]);
					ptmp += 3;
				}
				for (int i=0; i < 4; i++)
				{
					ptmp[0] = ptmp[1] = ptmp[2] = 0;
					ptmp += 3;
				}
				// MLE - EM algorithm
				EMAlg(tmpprob, out_k0, out_k1, out_loglik, &out_niter);
				break;

			case 1:
				// fill PrIBD, Pr(ibs state | ibd state j)
				for (int i=0; i < n; i++)
				{
					PrIBDTable(geno1[i], geno2[i], ptmp[0], ptmp[1], ptmp[2], AFreq[i]);
					ptmp[0] -= ptmp[2]; ptmp[1] -= ptmp[2];
					ptmp += 3;
				}
				for (int i=0; i < 4; i++)
				{
					ptmp[0] = ptmp[1] = ptmp[2] = 0;
					ptmp += 3;
				}
				// MLE - downhill simplex algorithm
				Simplex(tmpprob, out_k0, out_k1, out_loglik, &out_niter);
				break;
		}
	}


	/// to compute the log likelihood
	void Do_MLE_LogLik(const double *AFreq, const double *k0, const double *k1,
		double *tmp_AF, double *out_loglik)
	{
		InitAFreq(AFreq, tmp_AF);
		// initialize buffer
		vector<double> PrIBD(3*nTotalSNP);

		for (int i=0; i < nSamp; i++)
		{
			for (int j=i; j < nSamp; j++)
			{
				EM_Prepare(&PrIBD[0],
					&PackedGeno[0] + nPackedSNP*i, &PackedGeno[0] + nPackedSNP*j);
				out_loglik[i*nSamp+j] = out_loglik[j*nSamp+i] =
					EM_LogLik(&PrIBD[0], k0[i*nSamp+j], k1[i*nSamp+j]);
			}
		}
	}

	/// to compute the log likelihood
	void Do_MLE_LogLik_k01(const double *AFreq, const double k0, const double k1,
		double *tmp_AF, double *out_loglik)
	{
		InitAFreq(AFreq, tmp_AF);
		// initialize buffer
		vector<double> PrIBD(3*nTotalSNP);

		for (int i=0; i < nSamp; i++)
		{
			for (int j=i; j < nSamp; j++)
			{
				EM_Prepare(&PrIBD[0],
					&PackedGeno[0] + nPackedSNP*i, &PackedGeno[0] + nPackedSNP*j);
				out_loglik[i*nSamp+j] = out_loglik[j*nSamp+i] =
					EM_LogLik(&PrIBD[0], k0, k1);
			}
		}
	}
}


namespace INBREEDING
{
	// the individual inbreeding coefficients

	template<typename TYPE> static double _inb_mom(int n, TYPE snp[],
		double afreq[])
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
				if (R_FINITE(val)) { F += val; nValid++; }
			}
		}
		if (nValid > 0) F /= nValid;
		return F;
	}

	template<typename TYPE> static double _inb_mom_ratio(int n, TYPE snp[],
		double afreq[])
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

	template<typename TYPE> static double _inb_mle_loglik(double F, int n,
		TYPE snp[], double afreq[])
	{
		double rv = 0;
		for (int i=0; i < n; i++)
		{
			double p = afreq[i], val = R_NaN;
			switch (snp[i])
			{
				case 0:
					val = log((1-F)*(1-p)*(1-p) + F*(1-p)); break;
				case 1:
					val = log((1-F)*2*p*(1-p)); break;
				case 2:
					val = log((1-F)*p*p + F*p); break;
			}
			if (R_FINITE(val)) rv += val;
		}
		return rv;
	}

	template<typename TYPE> static double _inb_mle(int n, TYPE snp[],
		double afreq[], const double reltol, int *out_iternum)
	{
		// initial value
		double F = _inb_mom_ratio(n, snp, afreq);
		if (R_FINITE(F))
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
							if (R_FINITE(tmp)) { sum += tmp; m++; }
							break;
						case 1:
							m++; break;
						case 2:
							tmp = F / (F + p*(1-F));
							if (R_FINITE(tmp)) { sum += tmp; m++; }
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


using namespace IBD;


extern "C"
{
/// internal IBD function
static void IBD_Init_Buffer(vector<int> &buf_geno, vector<double> &buf_afreq)
{
	size_t nSamp = MCWorkingGeno.Space.SampleNum();
	size_t nPackedSNP = (MCWorkingGeno.Space.SNPNum() % 4 > 0) ?
		(MCWorkingGeno.Space.SNPNum()/4 + 1) : (MCWorkingGeno.Space.SNPNum()/4);
	size_t nTotal = nSamp * nPackedSNP;

	size_t buf_size = nTotal/sizeof(int) +
		((nTotal % sizeof(int) > 0) ? 1 : 0);
	size_t buf_snp_size = 4*nPackedSNP;

	buf_geno.resize(buf_size);
	buf_afreq.resize(buf_snp_size);
}


/// to compute the IBD coefficients by MLE
COREARRAY_DLL_EXPORT SEXP gnrIBD_MLE(SEXP AlleleFreq, SEXP KinshipConstraint,
	SEXP MaxIterCnt, SEXP RelTol, SEXP CoeffCorrect, SEXP method,
	SEXP IfOutNum, SEXP NumThread, SEXP _Verbose)
{
	bool verbose = SEXP_Verbose(_Verbose);

	COREARRAY_TRY

		// ========  To cache the genotype data  ========
		CachingSNPData("MLE IBD", verbose);

		// ======== MLE IBD ========

		vector<int> tmp_buffer;
		vector<double> tmp_AF;
		IBD_Init_Buffer(tmp_buffer, tmp_AF);

		// initialize the packed genotypes
		IBD::InitPackedGeno(&(tmp_buffer[0]));
		// initialize the internal matrix
		IBD::Init_EPrIBD_IBS(isNull(AlleleFreq) ? NULL : REAL(AlleleFreq),
			NULL, false);

		IBD::nIterMax = INTEGER(MaxIterCnt)[0];
		IBD::FuncRelTol = REAL(RelTol)[0];
		IBD::MethodMLE = INTEGER(method)[0];
		IBD::Loglik_Adjust = (LOGICAL(CoeffCorrect)[0] == TRUE);
		IBD::KinshipConstraint = (LOGICAL(KinshipConstraint)[0] == TRUE);

		// the upper-triangle genetic covariance matrix
		const R_xlen_t n = MCWorkingGeno.Space.SampleNum();
		CdMatTriDiag<IBD::TIBD> IBD(IBD::TIBD(), n);
		CdMatTriDiag<int> niter;
		if (LOGICAL(IfOutNum)[0] == TRUE)
			niter.Reset(n);

		// R SEXP objects
		PROTECT(rv_ans = NEW_LIST(4));
		SEXP afreq = PROTECT(NEW_NUMERIC(MCWorkingGeno.Space.SNPNum()));
		SET_ELEMENT(rv_ans, 2, afreq);

		// Calculate the IBD matrix
		IBD::Do_MLE_IBD_Calc(
			isNull(AlleleFreq) ? NULL : REAL(AlleleFreq), IBD,
			(LOGICAL(IfOutNum)[0] == TRUE) ? &niter : NULL,
			REAL(afreq), INTEGER(NumThread)[0], "MLE IBD:",
			&(tmp_AF[0]), verbose);

		// output
		SEXP k0, k1, IterNum=NULL;
		PROTECT(k0 = allocMatrix(REALSXP, n, n));
		SET_ELEMENT(rv_ans, 0, k0);
		PROTECT(k1 = allocMatrix(REALSXP, n, n));
		SET_ELEMENT(rv_ans, 1, k1);
		if (LOGICAL(IfOutNum)[0] == TRUE)
		{
			PROTECT(IterNum = allocMatrix(INTSXP, n, n));
			SET_ELEMENT(rv_ans, 3, IterNum);
		}

		double *out_k0 = REAL(k0);
		double *out_k1 = REAL(k1);
		int *out_niter = (IterNum) ? INTEGER(IterNum) : NULL;
		IBD::TIBD *p = IBD.get();
		int *pn = niter.get();
		for (int i=0; i < n; i++)
		{
			size_t pp = i*n + i;
			out_k0[pp] = out_k1[pp] = 0;
			if (out_niter) out_niter[pp] = 0;

			for (int j=i+1; j < n; j++, p++)
			{
				out_k0[i*n + j] = out_k0[j*n + i] = p->k0;
				out_k1[i*n + j] = out_k1[j*n + i] = p->k1;
				if (out_niter)
					out_niter[i*n + j] = out_niter[j*n + i] = *pn++;
			}
		}

		UNPROTECT((IterNum) ? 5 : 4);

	COREARRAY_CATCH
}


/// to compute the IBD coefficients by MLE
COREARRAY_DLL_EXPORT SEXP gnrIBD_MLE_Jacquard(SEXP AlleleFreq, SEXP MaxIterCnt,
	SEXP RelTol, SEXP CoeffCorrect, SEXP method, SEXP IfOutNum, SEXP NumThread,
	SEXP _Verbose)
{
	bool verbose = SEXP_Verbose(_Verbose);

	COREARRAY_TRY

		// ======== To cache the genotype data ========
		CachingSNPData("MLE IBD", verbose);

		// ======== MLE IBD ========

		vector<int> tmp_buffer;
		vector<double> tmp_AF;
		IBD_Init_Buffer(tmp_buffer, tmp_AF);

		// initialize the packed genotypes
		IBD::InitPackedGeno(&(tmp_buffer[0]));

		IBD::nIterMax = INTEGER(MaxIterCnt)[0];
		IBD::FuncRelTol = REAL(RelTol)[0];
		IBD::MethodMLE = INTEGER(method)[0];
		IBD::Loglik_Adjust = (LOGICAL(CoeffCorrect)[0] == TRUE);

		// the upper-triangle genetic covariance matrix
		const R_xlen_t n = MCWorkingGeno.Space.SampleNum();
		CdMatTriDiag<IBD::TIBD_Jacq> IBD(IBD::TIBD_Jacq(), n);
		CdMatTriDiag<int> niter;
		if (LOGICAL(IfOutNum)[0] == TRUE)
			niter.Reset(n);

		// R SEXP objects
		PROTECT(rv_ans = NEW_LIST(10));
		SEXP afreq = PROTECT(NEW_NUMERIC(MCWorkingGeno.Space.SNPNum()));
		SET_ELEMENT(rv_ans, 8, afreq);

		// Calculate the IBD matrix
		IBD::Do_MLE_IBD_Jacq(
			isNull(AlleleFreq) ? NULL : REAL(AlleleFreq), IBD,
			(LOGICAL(IfOutNum)[0] == TRUE) ? &niter : NULL,
			REAL(afreq), INTEGER(NumThread)[0], "MLE IBD:",
			&(tmp_AF[0]), verbose);

		// output
		SEXP D[8], IterNum=NULL;
		for (int i=0; i < 8; i++)
		{
			PROTECT(D[i] = allocMatrix(REALSXP, n, n));
			SET_ELEMENT(rv_ans, i, D[i]);
		}
		if (LOGICAL(IfOutNum)[0] == TRUE)
		{
			PROTECT(IterNum = allocMatrix(INTSXP, n, n));
			SET_ELEMENT(rv_ans, 9, IterNum);
		}

		double *out_d[8] = { REAL(D[0]), REAL(D[1]), REAL(D[2]),
			REAL(D[3]), REAL(D[4]), REAL(D[5]), REAL(D[6]), REAL(D[7]) };

		int *out_niter = (IterNum) ? INTEGER(IterNum) : NULL;
		IBD::TIBD_Jacq *p = IBD.get();
		int *pn = niter.get();
		for (int i=0; i < n; i++)
		{
			size_t pp = i*n + i;
			out_d[0][pp] = 1;
			out_d[1][pp] = out_d[2][pp] = out_d[3][pp] = out_d[4][pp] =
				out_d[5][pp] = out_d[6][pp] = out_d[7][pp] = 0;
			if (out_niter) out_niter[pp] = 0;

			for (int j=i+1; j < n; j++, p++)
			{
				size_t p1 = i*n + j, p2 = j*n + i;
				out_d[0][p1] = out_d[0][p2] = p->D1;
				out_d[1][p1] = out_d[1][p2] = p->D2;
				out_d[2][p1] = out_d[2][p2] = p->D3;
				out_d[3][p1] = out_d[3][p2] = p->D4;
				out_d[4][p1] = out_d[4][p2] = p->D5;
				out_d[5][p1] = out_d[5][p2] = p->D6;
				out_d[6][p1] = out_d[6][p2] = p->D7;
				out_d[7][p1] = out_d[7][p2] = p->D8;
				if (out_niter)
					out_niter[p1] = out_niter[p2] = *pn++;
			}
		}

		UNPROTECT((IterNum) ? 11: 10);

	COREARRAY_CATCH
}


/// to compute the IBD coefficients by MLE for a pair of individuals
COREARRAY_DLL_EXPORT SEXP gnrPairIBD(SEXP geno1, SEXP geno2, SEXP AlleleFreq,
	SEXP KinshipConstraint, SEXP MaxIterCnt, SEXP RelTol,
	SEXP CoeffCorrect, SEXP method)
{
	COREARRAY_TRY

		// initialize
		const int n = XLENGTH(geno1);
		IBD::nIterMax = INTEGER(MaxIterCnt)[0];
		IBD::FuncRelTol = REAL(RelTol)[0];
		IBD::MethodMLE = INTEGER(method)[0];
		IBD::Loglik_Adjust = (LOGICAL(CoeffCorrect)[0] == TRUE);
		IBD::KinshipConstraint = (LOGICAL(KinshipConstraint)[0] == TRUE);
		IBD::Init_EPrIBD_IBS(REAL(AlleleFreq), NULL, false, n);

		// get the initial values
		int IBS[3] = { 0, 0, 0 };
		int *g1 = INTEGER(geno1), *g2 = INTEGER(geno2);
		for (int i=0; i < n; i++, g1++, g2++)
		{
			if ((0 <= *g1) && (*g1 <= 2) && (0 <= *g2) && (*g2 <= 2))
			{
				IBS[2 - abs(*g1 - *g2)] ++;
			}
		}

		double out_k0, out_k1, out_loglik;
		int out_niter;
		IBD::Est_PLINK_Kinship(IBS[0], IBS[1], IBS[2], out_k0, out_k1,
			IBD::KinshipConstraint);

		// compute
		if (INTEGER(method)[0] >= 0)
		{
			vector<double> tmp_buffer(3*n + 3*4);
			IBD::Do_MLE_IBD_Pair(n, INTEGER(geno1), INTEGER(geno2),
				REAL(AlleleFreq), out_k0, out_k1, out_loglik, out_niter,
				&tmp_buffer[0]);
		} else {
			out_loglik = R_NaN;
			out_niter = 0;
		}

		PROTECT(rv_ans = NEW_LIST(4));
		SET_ELEMENT(rv_ans, 0, ScalarReal(out_k0));
		SET_ELEMENT(rv_ans, 1, ScalarReal(out_k1));
		SET_ELEMENT(rv_ans, 2, ScalarReal(out_loglik));
		SET_ELEMENT(rv_ans, 3, ScalarInteger(out_niter));
		UNPROTECT(1);

	COREARRAY_CATCH
}


/// to compute log likelihood of MLE
COREARRAY_DLL_EXPORT SEXP gnrIBD_LogLik(SEXP AlleleFreq, SEXP k0, SEXP k1)
{
	COREARRAY_TRY

		vector<int> tmp_buffer;
		vector<double> tmp_AF;
		IBD_Init_Buffer(tmp_buffer, tmp_AF);

		// ======== MLE IBD ========
		// initialize the packed genotypes
		IBD::InitPackedGeno(&(tmp_buffer[0]));

		// call
		const int n = MCWorkingGeno.Space.SampleNum();
		PROTECT(rv_ans = allocMatrix(REALSXP, n, n));
		IBD::Do_MLE_LogLik(REAL(AlleleFreq), REAL(k0), REAL(k1),
			&(tmp_AF[0]), REAL(rv_ans));
		UNPROTECT(1);

	COREARRAY_CATCH
}


/// to compute log likelihood of MLE
COREARRAY_DLL_EXPORT SEXP gnrIBD_LogLik_k01(SEXP AlleleFreq, SEXP k0, SEXP k1)
{
	COREARRAY_TRY

		vector<int> tmp_buffer;
		vector<double> tmp_AF;
		IBD_Init_Buffer(tmp_buffer, tmp_AF);

		// ======== MLE IBD ========
		// initialize the packed genotypes
		IBD::InitPackedGeno(&(tmp_buffer[0]));

		// call
		const int n = MCWorkingGeno.Space.SampleNum();
		PROTECT(rv_ans = allocMatrix(REALSXP, n, n));
		IBD::Do_MLE_LogLik_k01(REAL(AlleleFreq), REAL(k0)[0], REAL(k1)[0],
			&(tmp_AF[0]), REAL(rv_ans));
		UNPROTECT(1);

	COREARRAY_CATCH
}


/// to compute the value of log likelihood for a pair of individuals
COREARRAY_DLL_EXPORT SEXP gnrPairIBDLogLik(SEXP geno1, SEXP geno2,
	SEXP AlleleFreq, SEXP k0, SEXP k1)
{
	COREARRAY_TRY

		const int n = XLENGTH(geno1);
		int *g1 = INTEGER(geno1);
		int *g2 = INTEGER(geno2);
		double *afreq = REAL(AlleleFreq);

		// initialize the probability table
		vector<double> tmp_buffer(3*n);
		double *PrIBD = &(tmp_buffer[0]);
		for (int i=0; i < n; i++)
		{
			IBD::PrIBDTable(g1[i], g2[i], PrIBD[0], PrIBD[1], PrIBD[2],
				afreq[i]);
			PrIBD += 3;
		}

		// calculate log likelihood value
		double k[3] = { REAL(k0)[0], REAL(k1)[0],
			1 - REAL(k0)[0] - REAL(k1)[0] };
		double LogLik = 0;
		PrIBD = &(tmp_buffer[0]);
		for (int i=0; i < n; i++)
		{
			double sum = PrIBD[0]*k[0] + PrIBD[1]*k[1] + PrIBD[2]*k[2];
			if (sum > 0)
				LogLik += log(sum);
			PrIBD += 3;
		}

		// output
		rv_ans = ScalarReal(LogLik);

	COREARRAY_CATCH
}

/// to create a list of sample id
//  ID1 = matrix(sample.id, nrow=n, ncol=n, byrow=TRUE)
COREARRAY_DLL_EXPORT SEXP gnrIBDSelSampID_List1(SEXP SampID, SEXP Flag)
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
COREARRAY_DLL_EXPORT SEXP gnrIBDSelSampID_List2(SEXP SampID, SEXP Flag)
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


// the individual inbreeding coefficients

/// to compute the inbreeding coefficient
COREARRAY_DLL_EXPORT SEXP gnrIndInbCoef(SEXP snp, SEXP afreq, SEXP reltol)
{
	int n       = XLENGTH(snp);
	int *G      = INTEGER(AS_INTEGER(snp));
	double *AF  = REAL(AS_NUMERIC(afreq));

	if (XLENGTH(reltol) != 1)
		error("`reltol' should a real number.");
	double rtol = REAL(AS_NUMERIC(reltol))[0];

	COREARRAY_TRY
		rv_ans = ScalarReal(INBREEDING::_inb_mle<int>(n, G, AF, rtol, NULL));
	COREARRAY_CATCH
}


// to compute the inbreeding coefficient
COREARRAY_DLL_EXPORT SEXP gnrIndInb(SEXP afreq, SEXP method, SEXP reltol,
	SEXP num_iter)
{
	double *AF = REAL(AS_NUMERIC(afreq));
	const char *met  = CHAR(STRING_ELT(method, 0));
	
	if (XLENGTH(reltol) != 1)
		error("`reltol' should a real number.");
	double rtol = REAL(AS_NUMERIC(reltol))[0];

	if (XLENGTH(num_iter) != 1)
		error("`out.num.iter' should a logical value.");
	bool if_iternum = LOGICAL(num_iter)[0];

	COREARRAY_TRY

		// the number of SNPs
		const R_xlen_t n = MCWorkingGeno.Space.SNPNum();
		// buffer object
		CdBufSpace buf(MCWorkingGeno.Space, false, CdBufSpace::acInc);
		// the number of samples
		const R_xlen_t m = buf.IdxCnt();

		SEXP OutCoeff;
		PROTECT(OutCoeff = NEW_NUMERIC(m));
		double *pCoeff = REAL(OutCoeff);

		PROTECT(rv_ans = NEW_LIST(2));
		SET_ELEMENT(rv_ans, 0, OutCoeff);

		if (strcmp(met, "mom.weir") == 0)
		{
			for (long i=0; i < buf.IdxCnt(); i++)
			{
				C_UInt8 *p = buf.ReadGeno(i);
				pCoeff[i] = INBREEDING::_inb_mom_ratio(n, p, AF);
			}
		} else if (strcmp(met, "mom.visscher") == 0)
		{
			for (long i=0; i < buf.IdxCnt(); i++)
			{
				C_UInt8 *p = buf.ReadGeno(i);
				pCoeff[i] = INBREEDING::_inb_mom(n, p, AF);
			}
		} else if (strcmp(met, "mle") == 0)
		{
			SEXP Num = NULL;
			if (if_iternum)
			{
				PROTECT(Num = NEW_INTEGER(m));
				SET_ELEMENT(rv_ans, 1, Num);
			}
			for (long i=0; i < buf.IdxCnt(); i++)
			{
				C_UInt8 *p = buf.ReadGeno(i);
				int iter_num = -1;
				pCoeff[i] = INBREEDING::_inb_mle(n, p, AF, rtol, &iter_num);
				if (Num)
					INTEGER(Num)[i] = iter_num;
			}
			if (Num) UNPROTECT(1);
		} else
			throw "Invalid 'method'.";

		UNPROTECT(2);

	COREARRAY_CATCH
}

}

#endif  /* _HEADER_IBD_ */
