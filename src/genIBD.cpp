// ===========================================================
//     _/_/_/   _/_/_/  _/_/_/_/    _/_/_/_/  _/_/_/   _/_/_/
//      _/    _/       _/             _/    _/    _/   _/   _/
//     _/    _/       _/_/_/_/       _/    _/    _/   _/_/_/
//    _/    _/       _/             _/    _/    _/   _/
// _/_/_/   _/_/_/  _/_/_/_/_/     _/     _/_/_/   _/_/
// ===========================================================
//
// genIBD.cpp: Identity by descent (IBD) analysis on genome-wide association studies
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


// CoreArray library header
#include <dType.h>
#include <dGWASMath.h>
#include <dVect.h>
#include <CoreGDSLink.h>
#include <dGenGWAS.h>

// Standard library header
#include <vector>
#include <cmath>
#include <cfloat>
#include <memory>
#include <algorithm>
#include <R.h>


#ifdef COREARRAY_SIMD_SSE
#include <xmmintrin.h>
#endif
#ifdef COREARRAY_SIMD_SSE2
#include <emmintrin.h>
#endif


#ifndef _FuncIBD_H_
#define _FuncIBD_H_

namespace IBS
{
	using namespace CoreArray;

	/// The number of IBS 0 in the packed genotype
	extern UInt8 IBS0_Num_SNP[];
	/// The number of IBS 1 in the packed genotype
	extern UInt8 IBS1_Num_SNP[];
	/// The number of IBS 2 in the packed genotype
	extern UInt8 IBS2_Num_SNP[];
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

	/// The structure of IBS states
	struct TIBD
	{
		double k0, k1;
		TIBD() { k0 = k1 = 0; }
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
		auto_ptr<int> AA(new int[nSNP]), AB(new int[nSNP]), BB(new int[nSNP]);

		if (!in_afreq)
		{
			MCWorkingGeno.Space.GetABNumPerSNP(AA.get(), AB.get(), BB.get());
		}

		// for-loop each snp
		long nValid = 0;
		for (long i=0; i < nSNP; i++)
		{
			long n = 2 * (AA.get()[i] + AB.get()[i] + BB.get()[i]);
			double p = (n > 0) ? (double(2*AA.get()[i] + AB.get()[i]) / n) : conf_F64_NaN();

			if (in_afreq)
			{
				p = in_afreq[i];
				if (conf_IsFinite64(p))
					if ((p < 0) || (p > 1)) p = conf_F64_NaN();
			}
			
			if (out_afreq) out_afreq[i] = p;

			// Second, the expected probability of IBS i, given by IBD
			double q = 1-p, Na = n, x = 2*AA.get()[i]+AB.get()[i], y = 2*BB.get()[i]+AB.get()[i];
			double a00, a01, a02, a11, a12;

			//
			// The following codes are from PLINK/genome.cpp
			//
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

			//
			// End: PLINK/genome.cpp
			//

			if (conf_IsFinite64(a00) && conf_IsFinite64(a01) &&
				conf_IsFinite64(a02) && conf_IsFinite64(a11) && conf_IsFinite64(a12))
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
	CdMatTriDiag<TIBD> *IBD;
	TIBD *pMatIBD = NULL;
	int *pNIter = NULL;
	IdMatTriD IBD_idx;

	/// Genotype, stored in a packed mode
	UInt8 *PackedGeno = NULL;
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
		PackedGeno = (UInt8*)buffer;

		// buffer
		CdBufSpace Buf(MCWorkingGeno.Space, false, CdBufSpace::acInc);
		UInt8 *p = PackedGeno;
		for (long i=0; i < MCWorkingGeno.Space.SampleNum(); i++)
		{
			p = Buf.ReadPackedGeno(i, p);
		}
	}



	// **************************************************************** //

	/// Output the probability of IBS given by IBD
	/** \param g1  the first genotype (0 -- BB, 1 -- AB, 2 -- AA, other missing)
	 *  \param g2  the second genotype
	 *  \param t0  the probability of IBD given by IBD state = 0
	 *  \param t1  the probability of IBD given by IBD state = 1
	 *  \param t2  the probability of IBD given by IBD state = 2
	 *  \param p   the frequency of allele A
	**/
	void prIBDTable(int g1, int g2, double &t0, double &t1, double &t2, double p)
	{
		if ((0 < p) && (p < 1))
		{
			double q = 1 - p;
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
		if (conf_IsFinite64(_loglik)) { \
			if (out_loglik < _loglik) {	\
				out_loglik = _loglik; out_k0 = _k0; out_k1 = _k1; \
    		}	\
		}


	// ****************** MLE - EM algorithm ******************

	static void EM_Prepare(double *PrIBD, UInt8 *p1, UInt8 *p2)
	{
		const double *Freq = &MLEAlleleFreq[0];
		for (long nPackSNP=nPackedSNP; nPackSNP > 0; nPackSNP--)
		{
			// first genotype
			UInt8 g1=*p1++, g2=*p2++;
			prIBDTable(g1 & 0x03, g2 & 0x03, PrIBD[0], PrIBD[1], PrIBD[2], *Freq++);
			// second genotype
			g1 >>= 2; g2 >>= 2; PrIBD += 3;
			prIBDTable(g1 & 0x03, g2 & 0x03, PrIBD[0], PrIBD[1], PrIBD[2], *Freq++);
			// third genotype
			g1 >>= 2; g2 >>= 2; PrIBD += 3;
			prIBDTable(g1 & 0x03, g2 & 0x03, PrIBD[0], PrIBD[1], PrIBD[2], *Freq++);
			// fourth genotype
			g1 >>= 2; g2 >>= 2; PrIBD += 3;
			prIBDTable(g1 & 0x03, g2 & 0x03, PrIBD[0], PrIBD[1], PrIBD[2], *Freq++);
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
				return conf_F64_NegInf();
			pr += 3;
			// second genotype
			sum = pr[0]*k[0] + pr[1]*k[1] + pr[2]*k[2];
			if (sum > 0)
				LogLik += log(sum);
			else if (pr[0] > 0)
				return conf_F64_NegInf();
			pr += 3;
			// third genotype
			sum = pr[0]*k[0] + pr[1]*k[1] + pr[2]*k[2];
			if (sum > 0)
				LogLik += log(sum);
			else if (pr[0] > 0)
				return conf_F64_NegInf();
			pr += 3;
			// fourth genotype
			sum = pr[0]*k[0] + pr[1]*k[1] + pr[2]*k[2];
			if (sum > 0)
				LogLik += log(sum);
			else if (pr[0] > 0)
				return conf_F64_NegInf();
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
		if (conf_IsFinite64(LogLik))
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


	// ****************** MLE - downhill simplex algorithm ******************

	static void NM_Prepare(double *pr, UInt8 *p1, UInt8 *p2)
	{
		const double *Freq = &MLEAlleleFreq[0];
		for (long nPackSNP=nPackedSNP; nPackSNP > 0; nPackSNP--)
		{
			// first genotype
			UInt8 g1=*p1++, g2=*p2++;
			prIBDTable(g1 & 0x03, g2 & 0x03, pr[0], pr[1], pr[2], *Freq++);
			pr[0] -= pr[2]; pr[1] -= pr[2];
			// second genotype
			g1 >>= 2; g2 >>= 2; pr += 3;
			prIBDTable(g1 & 0x03, g2 & 0x03, pr[0], pr[1], pr[2], *Freq++);
			pr[0] -= pr[2]; pr[1] -= pr[2];
			// third genotype
			g1 >>= 2; g2 >>= 2; pr += 3;
			prIBDTable(g1 & 0x03, g2 & 0x03, pr[0], pr[1], pr[2], *Freq++);
			pr[0] -= pr[2]; pr[1] -= pr[2];
			// fourth genotype
			g1 >>= 2; g2 >>= 2; pr += 3;
			prIBDTable(g1 & 0x03, g2 & 0x03, pr[0], pr[1], pr[2], *Freq++);
			pr[0] -= pr[2]; pr[1] -= pr[2];
			pr += 3;
		}
	}

	/// return log likelihood value for the specified pair
	static double NM_LogLik(const double *PrIBD, const double k0, const double k1)
	{
		// Check whether within the region
		if ((k0<0) || (k1<0) || (k0+k1>1)) return conf_F64_NegInf();

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
				return conf_F64_NegInf();
			pr += 3;
			// second genotype
			sum = pr[0]*k0 + pr[1]*k1 + pr[2];
			if (sum > 0)
				LogLik += log(sum);
			else if (pr[0] > 0)
				return conf_F64_NegInf();
			pr += 3;
			// third genotype
			sum = pr[0]*k0 + pr[1]*k1 + pr[2];
			if (sum > 0)
				LogLik += log(sum);
			else if (pr[0] > 0)
				return conf_F64_NegInf();
			pr += 3;
			// fourth genotype
			sum = pr[0]*k0 + pr[1]*k1 + pr[2];
			if (sum > 0)
				LogLik += log(sum);
			else if (pr[0] > 0)
				return conf_F64_NegInf();
			pr += 3;
		}
		return LogLik;
	}

	/// the optimize function
	static double _optim(const double *x, void *ex)
	{
		double rv = -NM_LogLik((double*)ex, x[0], x[1]);
		if (!conf_IsFinite64(rv)) rv = 1e+30;
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
	void Entry_MLEIBD(TdThread Thread, int ThreadIndex, void *Param)
	{
		// initialize buffer
		auto_ptr<double> PrIBD(new double[3*nTotalSNP]);

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

			UInt8 *g1 = &PackedGeno[0] + nPackedSNP*idx.Row();
			UInt8 *g2 = &PackedGeno[0] + nPackedSNP*idx.Column();

			// calculate the initial values from PLINK
			UInt8 *p1 = g1, *p2 = g2;
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
					EM_Prepare(PrIBD.get(), g1, g2);
					// MLE - EM algorithm
					EMAlg(PrIBD.get(), pIBD->k0, pIBD->k1, _loglik, pniter);
					break;
				case 1:
					// fill PrIBD, Pr(ibs state|ibd state j)
					NM_Prepare(PrIBD.get(), g1, g2);
					// MLE - downhill simplex algorithm
					Simplex(PrIBD.get(), pIBD->k0, pIBD->k1, _loglik, pniter);
					break;
			}
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
				if (conf_IsFinite64(AFreq[i]))
					MLEAlleleFreq[i] = AFreq[i];
			}
		} else {
			// init
			vector<int> n(nTotalSNP);
			for (int i=0; i < nTotalSNP; i++) n[i] = 0;
			for (int i=0; i < nTotalSNP; i++) MLEAlleleFreq[i] = 0;
			// for-loop for allele frequency
			UInt8 *p = &PackedGeno[0];
			for (int iSamp=0; iSamp < nSamp; iSamp++)
			{
				for (int i=0; i < nPackedSNP; i++)
				{
					UInt8 L = *p++;
					for (int k=0; k < 4; k++)
					{
						UInt8 B = L & 0x03; L >>= 2;
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


	/// to conduct MLE IBD
	void Do_MLE_IBD_Calculate(const double *AFreq, CdMatTriDiag<TIBD> &PublicIBD,
		CdMatTriDiag<int> *PublicNIter,
		double *out_AFreq, int NumThread, const char *Info, double *tmpAF, bool verbose)
	{
		InitAFreq(AFreq, tmpAF);
		for (int i=0; i < MCWorkingGeno.Space.SNPNum(); i++)
			out_AFreq[i] = MLEAlleleFreq[i];

		IBD = &PublicIBD; pMatIBD = PublicIBD.get();
		pNIter = (PublicNIter) ? PublicNIter->get() : NULL;
		IBD_idx.reset(nSamp);
		nMatTriD = PublicIBD.Size(); idxMatTriD = 0;

		// initialize the mutex object
		_Mutex = plc_InitMutex();

		// Initialize progress information
		MCWorkingGeno.Progress.Info = Info;
		MCWorkingGeno.Progress.Show() = verbose;
		MCWorkingGeno.Progress.Init(nMatTriD);

		// Threads
		plc_DoBaseThread(Entry_MLEIBD, NULL, NumThread);

		// destory the mutex objects
		plc_DoneMutex(_Mutex); _Mutex = NULL;
	}

	/// to conduct MLE IBD for a pair of individuals, no missing genotypes (0, 1, 2)
	void Do_MLE_IBD_Pair(int n, const int *geno1, const int *geno2, const double *AFreq,
		double &out_k0, double &out_k1, double &out_loglik, int &out_niter, double tmpprob[])
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
					prIBDTable(geno1[i], geno2[i], ptmp[0], ptmp[1], ptmp[2], AFreq[i]);
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
					prIBDTable(geno1[i], geno2[i], ptmp[0], ptmp[1], ptmp[2], AFreq[i]);
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
		auto_ptr<double> PrIBD(new double[3*nTotalSNP]);
		for (int i=0; i < nSamp; i++)
		{
			for (int j=i; j < nSamp; j++)
			{
				EM_Prepare(PrIBD.get(),
					&PackedGeno[0] + nPackedSNP*i, &PackedGeno[0] + nPackedSNP*j);
				out_loglik[i*nSamp+j] = out_loglik[j*nSamp+i] =
					EM_LogLik(PrIBD.get(), k0[i*nSamp+j], k1[i*nSamp+j]);
			}
		}
	}

	/// to compute the log likelihood
	void Do_MLE_LogLik_k01(const double *AFreq, const double k0, const double k1,
		double *tmp_AF, double *out_loglik)
	{
		InitAFreq(AFreq, tmp_AF);
		// initialize buffer
		auto_ptr<double> PrIBD(new double[3*nTotalSNP]);
		for (int i=0; i < nSamp; i++)
		{
			for (int j=i; j < nSamp; j++)
			{
				EM_Prepare(PrIBD.get(),
					&PackedGeno[0] + nPackedSNP*i, &PackedGeno[0] + nPackedSNP*j);
				out_loglik[i*nSamp+j] = out_loglik[j*nSamp+i] =
					EM_LogLik(PrIBD.get(), k0, k1);
			}
		}
	}
}

#endif  /* _FuncIBD_H_ */
