// ===========================================================
//     _/_/_/   _/_/_/  _/_/_/_/    _/_/_/_/  _/_/_/   _/_/_/
//      _/    _/       _/             _/    _/    _/   _/   _/
//     _/    _/       _/_/_/_/       _/    _/    _/   _/_/_/
//    _/    _/       _/             _/    _/    _/   _/
// _/_/_/   _/_/_/  _/_/_/_/_/     _/     _/_/_/   _/_/
// ===========================================================
//
// genLD.cpp: Linkage Disequilibrium (LD) analysis on genome-wide association studies
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
#include <dVect.h>
#include <CoreGDSLink.h>
#include <dGenGWAS.h>

// Standard library header
#include <cmath>
#include <cfloat>
#include <memory>
#include <list>
#include <algorithm>


#ifdef COREARRAY_SIMD_SSE
#include <xmmintrin.h>
#endif
#ifdef COREARRAY_SIMD_SSE2
#include <emmintrin.h>
#endif


#ifndef _FuncLD_H_
#define _FuncLD_H_

namespace LD
{
	// using namespace
	using namespace std;
	using namespace CoreArray;
	using namespace GDSInterface;
	using namespace GWAS;


	// ---------------------------------------------------------------------

	/// buffer size
	static const int _size = 256*256;

	/// The number of valid pair of SNPs in the packed genotypes
	UInt8 Valid_Num_SNP[_size];
	/// The number of aa in a pair of packed SNPs
	UInt8 Num_aa_SNP[_size];
	/// The number of aA in a pair of packed SNPs
	UInt8 Num_aA_SNP[_size];
	/// The number of AA in a pair of packed SNPs
	UInt8 Num_AA_SNP[_size];

	/// The number of AA BB pairs in the packed genotypes
	UInt8 Num_AA_BB_SNP[_size];
	/// The number of AA BB pairs in the packed genotypes
	UInt8 Num_aa_bb_SNP[_size];
	/// The number of aa BB pairs in the packed genotypes
	UInt8 Num_aa_BB_SNP[_size];
	/// The number of AA bb pairs in the packed genotypes
	UInt8 Num_AA_bb_SNP[_size];
	/// The number of double het pairs in the packed genotypes
	UInt8 Num_DH_SNP[_size];

	/// The sum of X for a pair of SNPs in the packed genotypes
	UInt8 Sum_X_SNP[_size];
	/// The sum of X^2 for a pair of SNPs in the packed genotypes
	UInt8 Sum_X_2_SNP[_size];
	/// The sum of X*Y for a pair of SNPs in the packed genotypes
	UInt8 Sum_XY_SNP[_size];

	/// The number of haplotype A / A in the packed genotypes
	UInt8 Num_A_A[_size];
	/// The number of haplotype A / B in the packed genotypes
	UInt8 Num_A_B[_size];
	/// The number of haplotype B / A in the packed genotypes
	UInt8 Num_B_A[_size];
	/// The number of haplotype B / B in the packed genotypes
	UInt8 Num_B_B[_size];
	/// The number of DH (double hets)
	UInt8 Num_DH2[_size];

	/// Genotype, stored in a packed mode
	auto_ptr<UInt8> PackedGeno;
	/// the number of samples and snps
	long nPackedSamp, nSNP;


	/// initial object
	class TInit
	{
	public:
		TInit()
		{
			#define PACKED_COND(cond, var, op)	\
				for (int s=0; s < _size; s++)	\
				{	\
					int g1 = s/256, g2 = s%256;	\
					int sum = 0;	\
					for (int i=0; i < 4; i++)	\
					{	\
						int b1 = g1 & 0x03, b2 = g2 & 0x03;	\
						if (cond) op;	\
						g1 >>= 2; g2 >>= 2;	\
					}	\
					var[s] = sum;	\
				}

			// Initialize the number of valid pair of SNPs in the packed genotypes
			PACKED_COND((b1 < 3) && (b2 < 3), Valid_Num_SNP, sum++);
			/// The number of aa in a pair of packed SNPs
			PACKED_COND((b1 == 0) && (b2 < 3), Num_aa_SNP, sum++);
			/// The number of aA in a pair of packed SNPs
			PACKED_COND((b1 == 1) && (b2 < 3), Num_aA_SNP, sum++);
			/// The number of AA in a pair of packed SNPs
			PACKED_COND((b1 == 2) && (b2 < 3), Num_AA_SNP, sum++);

			/// The number of AA BB pairs in the packed genotypes
			PACKED_COND((b1 == 2) && (b2 == 2), Num_AA_BB_SNP, sum++);
			/// The number of aa bb pairs in the packed genotypes
			PACKED_COND((b1 == 0) && (b2 == 0), Num_aa_bb_SNP, sum++);
			/// The number of aa BB pairs in the packed genotypes
			PACKED_COND((b1 == 0) && (b2 == 2), Num_aa_BB_SNP, sum++);
			/// The number of AA bb pairs in the packed genotypes
			PACKED_COND((b1 == 2) && (b2 == 0), Num_AA_bb_SNP, sum++);

			/// The sum of X for a pair of SNPs in the packed genotypes
			PACKED_COND((b1 < 3) && (b2 < 3), Sum_X_SNP, sum += b1);
			/// The sum of X^2 for a pair of SNPs in the packed genotypes
			PACKED_COND((b1 < 3) && (b2 < 3), Sum_X_2_SNP, sum += b1*b1);
			/// The sum of X*Y for a pair of SNPs in the packed genotypes
			PACKED_COND((b1 < 3) && (b2 < 3), Sum_XY_SNP, sum += b1*b2);

			static int IncArray[9][5] = {
				{ 0, 0, 0, 2, 0 }, // BB, BB
				{ 0, 0, 1, 1, 0 }, // BB, AB
				{ 0, 0, 2, 0, 0 }, // BB, AA
				{ 0, 1, 0, 1, 0 }, // AB, BB
				{ 0, 0, 0, 0, 2 }, // AB, AB (double hets)
				{ 1, 0, 1, 0, 0 }, // AB, AA
				{ 0, 2, 0, 0, 0 }, // AA, BB
				{ 1, 1, 0, 0, 0 }, // AA, AB
				{ 2, 0, 0, 0, 0 }  // AA, AA
			};
			/// The number of haplotype A / A in the packed genotypes
			PACKED_COND((b1 < 3) && (b2 < 3), Num_A_A, sum += IncArray[b1*3+b2][0]);
			/// The number of haplotype A / B in the packed genotypes
			PACKED_COND((b1 < 3) && (b2 < 3), Num_A_B, sum += IncArray[b1*3+b2][1]);
			/// The number of haplotype B / A in the packed genotypes
			PACKED_COND((b1 < 3) && (b2 < 3), Num_B_A, sum += IncArray[b1*3+b2][2]);
			/// The number of haplotype B / B in the packed genotypes
			PACKED_COND((b1 < 3) && (b2 < 3), Num_B_B, sum += IncArray[b1*3+b2][3]);
			/// The number of DH (double hets)
			PACKED_COND((b1 < 3) && (b2 < 3), Num_DH2, sum += IncArray[b1*3+b2][4]);
		}
	} InitObj;


	// ---------------------------------------------------------------------
	/// Composite LD estimation
	static double PairComposite(const int *snp1, const int *snp2, int cnt)
	{
		// Init data
		long n, naa, naA, nAA, nbb, nbB, nBB, nAABB, naabb, naaBB, nAAbb;
		n = naa = naA = nAA = nbb = nbB = nBB =
			nAABB = naabb = naaBB = nAAbb = 0;

		for (; cnt > 0; cnt--, snp1++, snp2++)
		{
			UInt8 g1 = (0<=*snp1 && *snp1<=2) ? (*snp1 | ~0x03) : 0xFF;
			UInt8 g2 = (0<=*snp2 && *snp2<=2) ? (*snp2 | ~0x03) : 0xFF;
			size_t _p = (size_t(g1) << 8) | (g2);
			size_t _q = (size_t(g2) << 8) | (g1);

			n += Valid_Num_SNP[_p];
			naa += Num_aa_SNP[_p]; nbb += Num_aa_SNP[_q];
			naA += Num_aA_SNP[_p]; nbB += Num_aA_SNP[_q];
			nAA += Num_AA_SNP[_p]; nBB += Num_AA_SNP[_q];
			nAABB += Num_AA_BB_SNP[_p]; naabb += Num_aa_bb_SNP[_p];
			naaBB += Num_aa_BB_SNP[_p]; nAAbb += Num_AA_bb_SNP[_p];
		}
		if (n > 0)
		{
			double delta =
				double(nAABB + naabb - naaBB - nAAbb) / (2*n) -
				double(naa-nAA)*double(nbb-nBB) / (2.0*n*n);
			double pa = double(2*naa + naA) / (2*n);
			double pA = 1 - pa, pAA = double(nAA)/n;
			double pb = double(2*nbb + nbB) / (2*n);
			double pB = 1 - pb, pBB = double(nBB)/n;
			double DA = pAA - pA*pA;
			double DB = pBB - pB*pB;
			double t = (pA*pa + DA) * (pB*pb + DB);
			if (t > 0)
				return delta / sqrt(t);
		}
		return conf_F64_NaN();
	}
	static double PairComposite(const UInt8 *snp1, const UInt8 *snp2)
	{
		// Init data
		long n, naa, naA, nAA, nbb, nbB, nBB, nAABB, naabb, naaBB, nAAbb;
		n = naa = naA = nAA = nbb = nbB = nBB =
			nAABB = naabb = naaBB = nAAbb = 0;

		for (long cnt=nPackedSamp; cnt > 0; cnt--, snp1++, snp2++)
		{
			size_t _p = (size_t(*snp1) << 8) | (*snp2);
			size_t _q = (size_t(*snp2) << 8) | (*snp1);
			n += Valid_Num_SNP[_p];
			naa += Num_aa_SNP[_p]; nbb += Num_aa_SNP[_q];
			naA += Num_aA_SNP[_p]; nbB += Num_aA_SNP[_q];
			nAA += Num_AA_SNP[_p]; nBB += Num_AA_SNP[_q];
			nAABB += Num_AA_BB_SNP[_p]; naabb += Num_aa_bb_SNP[_p];
			naaBB += Num_aa_BB_SNP[_p]; nAAbb += Num_AA_bb_SNP[_p];
		}
		if (n > 0)
		{
			double delta =
				double(nAABB + naabb - naaBB - nAAbb) / (2*n) -
				double(naa-nAA)*double(nbb-nBB) / (2.0*n*n);
			double pa = double(2*naa + naA) / (2*n);
			double pA = 1 - pa, pAA = double(nAA)/n;
			double pb = double(2*nbb + nbB) / (2*n);
			double pB = 1 - pb, pBB = double(nBB)/n;
			double DA = pAA - pA*pA;
			double DB = pBB - pB*pB;
			double t = (pA*pa + DA) * (pB*pb + DB);
			if (t > 0)
				return delta / sqrt(t);
		}
		return conf_F64_NaN();
	}

	// ---------------------------------------------------------------------
	/// return the value of log(val + espison)
	inline static double PLOG(const double val)
	{
		return log(val + DBL_EPSILON);
	}

	/// EM ---- estimate the proportions of the haplotype
	/** \param g1    the first packed genotypes
		\param g2    the second packed genotypes
		\param pA_A  output the proportion of A / A
		\param pA_B  output the proportion of A / B
		\param pB_A  output the proportion of B / A
		\param pB_B  output the proportion of B / B
		pA_A + pA_B + pB_A + pB_B = 1.
	*/
	inline static void ProportionHaplo(long nA_A, long nA_B, long nB_A, long nB_B, long nDH2,
		double &pA_A, double &pA_B, double &pB_A, double &pB_B)
	{
		// initial parameters
		const double EM_Init_Factor = 0.01;
		const long nIterMax = 1000;
		const double FuncRelTol = sqrt(DBL_EPSILON);

		// If we have missing data (i.e. un-phased double hets),
		double nTotal = nA_A + nA_B + nB_A + nB_B + nDH2;
		// then use EM algorithm.
		if ((nTotal > 0) && (nDH2 > 0))
		{
			// Set initial probs, the initial values will be always greater than 0
			double Div = nA_A + nA_B + nB_A + nB_B + 4.0*EM_Init_Factor;
			pA_A = (nA_A + EM_Init_Factor) / Div;
			pA_B = (nA_B + EM_Init_Factor) / Div;
			pB_A = (nB_A + EM_Init_Factor) / Div;
			pB_B = (nB_B + EM_Init_Factor) / Div;

			long nDH = nDH2 / 2;
			double OldLogLik = nA_A*PLOG(pA_A) + nA_B*PLOG(pA_B) +
				nB_A*PLOG(pB_A) + nB_B*PLOG(pB_B) +
				nDH*PLOG(pA_A*pB_B + pA_B*pB_A);
			double ConTol = fabs(FuncRelTol * OldLogLik);
			if (ConTol < DBL_EPSILON) ConTol = DBL_EPSILON;

			// do loop
			for (long IterNum=1; IterNum <= nIterMax; IterNum++)
			{
				// E-step
				double pAA_BB = pA_A * pB_B;
				double pAB_BA = pA_B * pB_A;
				// number of double hets which are AA + BB
				double nDH_AA_BB = pAA_BB/(pAA_BB + pAB_BA) * nDH;
				// number of double hets which are AB + BA
				double nDH_AB_BA = nDH - nDH_AA_BB;

				// M-step
				pA_A = (nA_A + nDH_AA_BB) / nTotal;
				pA_B = (nA_B + nDH_AB_BA) / nTotal;
				pB_A = (nB_A + nDH_AB_BA) / nTotal;
				pB_B = (nB_B + nDH_AA_BB) / nTotal;

				// check the likelihood function
				double LogLik = nA_A*PLOG(pA_A) + nA_B*PLOG(pA_B) +
					nB_A*PLOG(pB_A) + nB_B*PLOG(pB_B) +
					nDH*PLOG(pA_A*pB_B + pA_B*pB_A);

				// cout << i << ": " << FloatToStr(pA_A) << ", " << FloatToStr(pA_B)
				//	<< ", " << FloatToStr(pB_A) << ", " << FloatToStr(pB_B)
				//	<< endl << FloatToStr(LogLik) << ", " << FloatToStr(OldLogLik);
				// getchar();

				if (fabs(LogLik - OldLogLik) <= ConTol)
					break;
				OldLogLik = LogLik;
			}
		} else {
			pA_A = nA_A / nTotal; pA_B = nA_B / nTotal;
			pB_A = nB_A / nTotal; pB_B = nB_B / nTotal;
		}
	}

	/// r by EM algorithm assuming HWE
	static double PairR(const int *snp1, const int *snp2, int cnt,
		double &pA_A, double &pA_B, double &pB_A, double &pB_B)
	{
		// The number of haplotypes, nDH - double hets
		long nA_A, nA_B, nB_A, nB_B, nDH2;
		nA_A = nA_B = nB_A = nB_B = nDH2 = 0;

		// Calculate the numbers of haplotypes
		for (; cnt > 0; cnt--, snp1++, snp2++)
		{
			UInt8 g1 = (0<=*snp1 && *snp1<=2) ? (*snp1 | ~0x03) : 0xFF;
			UInt8 g2 = (0<=*snp2 && *snp2<=2) ? (*snp2 | ~0x03) : 0xFF;
			size_t _p = (size_t(g1) << 8) | (g2);
			nA_A += Num_A_A[_p]; nA_B += Num_A_B[_p];
			nB_A += Num_B_A[_p]; nB_B += Num_B_B[_p];
			nDH2 += Num_DH2[_p];
		}

		ProportionHaplo(nA_A, nA_B, nB_A, nB_B, nDH2, pA_A, pA_B, pB_A, pB_B);

		double pA = pA_A + pA_B;
		double p_A = pA_A + pB_A;
		double pB = pB_A + pB_B;
		double p_B = pA_B + pB_B;
		double D = pA_A - pA * p_A;
		return D / sqrt(pA * p_A * pB * p_B);
	}
	static double PairR(const UInt8 *snp1, const UInt8 *snp2)
	{
		// The number of haplotypes, nDH - double hets
		long nA_A, nA_B, nB_A, nB_B, nDH2;
		nA_A = nA_B = nB_A = nB_B = nDH2 = 0;

		// Calculate the numbers of haplotypes
		for (long cnt=nPackedSamp; cnt > 0; cnt--, snp1++, snp2++)
		{
			size_t _p = (size_t(*snp1) << 8) | (*snp2);
			nA_A += Num_A_A[_p]; nA_B += Num_A_B[_p];
			nB_A += Num_B_A[_p]; nB_B += Num_B_B[_p];
			nDH2 += Num_DH2[_p];
		}

		double pA_A, pA_B, pB_A, pB_B;
		ProportionHaplo(nA_A, nA_B, nB_A, nB_B, nDH2, pA_A, pA_B, pB_A, pB_B);

		double pA = pA_A + pA_B;
		double p_A = pA_A + pB_A;
		double pB = pB_A + pB_B;
		double p_B = pA_B + pB_B;
		double D = pA_A - pA * p_A;
		return D / sqrt(pA * p_A * pB * p_B);
	}

	// ---------------------------------------------------------------------
	/// D' by EM algorithm
	static double PairDPrime(const int *snp1, const int *snp2, int cnt,
		double &pA_A, double &pA_B, double &pB_A, double &pB_B)
	{
		// The number of haplotypes, nDH - double hets
		long nA_A, nA_B, nB_A, nB_B, nDH2;
		nA_A = nA_B = nB_A = nB_B = nDH2 = 0;

		// Calculate the numbers of haplotypes
		for (; cnt > 0; cnt--, snp1++, snp2++)
		{
			UInt8 g1 = (0<=*snp1 && *snp1<=2) ? (*snp1 | ~0x03) : 0xFF;
			UInt8 g2 = (0<=*snp2 && *snp2<=2) ? (*snp2 | ~0x03) : 0xFF;
			size_t _p = (size_t(g1) << 8) | (g2);
			nA_A += Num_A_A[_p]; nA_B += Num_A_B[_p];
			nB_A += Num_B_A[_p]; nB_B += Num_B_B[_p];
			nDH2 += Num_DH2[_p];
		}

		ProportionHaplo(nA_A, nA_B, nB_A, nB_B, nDH2, pA_A, pA_B, pB_A, pB_B);

		double pA = pA_A + pA_B;
		double pB = pB_A + pB_B;
		double p_A = pA_A + pB_A;
		double p_B = pA_B + pB_B;
		double D = pA_A - pA*p_A;
		// double r = D / sqrt(pA * p_A * pB * p_B);
		D = D / ((D>=0) ? min(pA*p_B, pB*p_A) : max(-pA*p_A, -pB*p_B));
		return D;
	}
	static double PairDPrime(const UInt8 *snp1, const UInt8 *snp2)
	{
		// The number of haplotypes, nDH - double hets
		long nA_A, nA_B, nB_A, nB_B, nDH2;
		nA_A = nA_B = nB_A = nB_B = nDH2 = 0;

		// Calculate the numbers of haplotypes
		for (long cnt=nPackedSamp; cnt > 0; cnt--, snp1++, snp2++)
		{
			size_t _p = (size_t(*snp1) << 8) | (*snp2);
			nA_A += Num_A_A[_p]; nA_B += Num_A_B[_p];
			nB_A += Num_B_A[_p]; nB_B += Num_B_B[_p];
			nDH2 += Num_DH2[_p];
		}

		double pA_A, pA_B, pB_A, pB_B;
		ProportionHaplo(nA_A, nA_B, nB_A, nB_B, nDH2, pA_A, pA_B, pB_A, pB_B);

		double pA = pA_A + pA_B;
		double pB = pB_A + pB_B;
		double p_A = pA_A + pB_A;
		double p_B = pA_B + pB_B;
		double D = pA_A - pA*p_A;
		// double r = D / sqrt(pA * p_A * pB * p_B);
		D = D / ((D>=0) ? min(pA*p_B, pB*p_A) : max(-pA*p_A, -pB*p_B));
		return D;
	}

	// ---------------------------------------------------------------------
	/// D' by EM algorithm
	static double PairCorr(const int *snp1, const int *snp2, int cnt)
	{
		// Init data
		long n, X, XX, Y, YY, XY;
		n = X = XX = Y = YY = XY = 0;

		// Calculate the numbers of haplotypes
		for (; cnt > 0; cnt--, snp1++, snp2++)
		{
			UInt8 g1 = (0<=*snp1 && *snp1<=2) ? (*snp1 | ~0x03) : 0xFF;
			UInt8 g2 = (0<=*snp2 && *snp2<=2) ? (*snp2 | ~0x03) : 0xFF;
			size_t _p = (size_t(g1) << 8) | (g2);
			size_t _q = (size_t(g2) << 8) | (g1);
			n += Valid_Num_SNP[_p];
			X += Sum_X_SNP[_p]; Y += Sum_X_SNP[_q];
			XX += Sum_X_2_SNP[_p]; YY += Sum_X_2_SNP[_q];
			XY += Sum_XY_SNP[_p];
		}
		if (n > 0)
		{
			double c1 = XX - double(X)*X/n;
			double c2 = YY - double(Y)*Y/n;
			double val = c1*c2;
			if (val > 0)
				return (XY - double(X)*Y/n) / sqrt(val);
		}
		return conf_F64_NaN();
	}
	static double PairCorr(const UInt8 *snp1, const UInt8 *snp2)
	{
		// Init data
		long n, X, XX, Y, YY, XY;
		n = X = XX = Y = YY = XY = 0;

		// Calculate the numbers of haplotypes
		for (long cnt=nPackedSamp; cnt > 0; cnt--, snp1++, snp2++)
		{
			size_t _p = (size_t(*snp1) << 8) | (*snp2);
			size_t _q = (size_t(*snp2) << 8) | (*snp1);
			n += Valid_Num_SNP[_p];
			X += Sum_X_SNP[_p]; Y += Sum_X_SNP[_q];
			XX += Sum_X_2_SNP[_p]; YY += Sum_X_2_SNP[_q];
			XY += Sum_XY_SNP[_p];
		}
		if (n > 0)
		{
			double c1 = XX - double(X)*X/n;
			double c2 = YY - double(Y)*Y/n;
			double val = c1*c2;
			if (val > 0)
				return (XY - double(X)*Y/n) / sqrt(val);
		}
		return conf_F64_NaN();
	}






	// ---------------------------------------------------------------------
	// public parameters

	// the method of Linkage Disequilibrium (LD)
	// 1 -- "composite"  Composite LD coefficients (by default)
	// 2 -- "r"          LD coefficient (by EM algorithm)
	// 3 -- "dprime"     D' coefficient
	// 4 -- "corr"       Correlation coefficient (BB, AB, AA are codes as 0, 1, 2)
	int LD_Method = 1;

	/// initialize the variable "Geno" with packed genotypes
	void InitPackedGeno()
	{
		// set # of samples and snps
		nSNP = MCWorkingGeno.Space.SNPNum();
		nPackedSamp = (MCWorkingGeno.Space.SampleNum() % 4 > 0) ?
			(MCWorkingGeno.Space.SampleNum()/4 + 1) : (MCWorkingGeno.Space.SampleNum()/4);
		PackedGeno.reset(new UInt8[nPackedSamp*nSNP]);

		// buffer
		CdBufSpace Buf(MCWorkingGeno.Space, true, CdBufSpace::acInc);
		UInt8 *p = PackedGeno.get();
		for (long i=0; i < MCWorkingGeno.Space.SNPNum(); i++)
		{
			p = Buf.ReadPackedGeno(i, p);
		}
	}

	/// destroy the variable "Geno" with packed genotypes
	void DonePackedGeno()
	{
		PackedGeno.reset(NULL);
	}


	// to compute pair LD
	double calcLD(const int *snp1, const int *snp2, int cnt,
		double &pA_A, double &pA_B, double &pB_A, double &pB_B)
	{
		switch (LD_Method)
		{
			case 1:
				return PairComposite(snp1, snp2, cnt);
			case 2:
				return PairR(snp1, snp2, cnt, pA_A, pA_B, pB_A, pB_B);
			case 3:
				return PairDPrime(snp1, snp2, cnt, pA_A, pA_B, pB_A, pB_B);
			case 4:
				return PairCorr(snp1, snp2, cnt);
		}
		return conf_F64_NaN();
	}

	// to compute LD matrix (n x n)
	void calcLD_mat(int nThread, double *out_LD)
	{
		for (int i=0; i < nSNP; i++)
		{
			out_LD[i*nSNP+i] = 1;
			for (int j=i+1; j < nSNP; j++)
			{
				double &p1 = out_LD[i*nSNP+j], &p2 = out_LD[j*nSNP+i];
				UInt8 *s1 = PackedGeno.get() + i*nPackedSamp;
				UInt8 *s2 = PackedGeno.get() + j*nPackedSamp;
				switch (LD_Method)
				{
					case 1:
						p1 = p2 = PairComposite(s1, s2); break;
					case 2:
						p1 = p2 = PairR(s1, s2); break;
					case 3:
						p1 = p2 = PairDPrime(s1, s2); break;
					case 4:
						p1 = p2 = PairCorr(s1, s2); break;
					default:
						p1 = p2 = conf_F64_NaN();
				}
			}
		}
	}

	// to compute LD matrix (n_slide x n)
	void calcLD_slide_mat(int nThread, double *out_LD, int n_slide)
	{
		for (int i=0; i < nSNP; i++)
		{
			for (int j=i+1; (j < nSNP) && ((j-i) <= n_slide); j++)
			{
				int d = j - i - 1;
				double &pVal = out_LD[i*n_slide + d];				
				UInt8 *s1 = PackedGeno.get() + i*nPackedSamp;
				UInt8 *s2 = PackedGeno.get() + j*nPackedSamp;
				switch (LD_Method)
				{
					case 1:
						pVal = PairComposite(s1, s2); break;
					case 2:
						pVal = PairR(s1, s2); break;
					case 3:
						pVal = PairDPrime(s1, s2); break;
					case 4:
						pVal = PairCorr(s1, s2); break;
					default:
						pVal = conf_F64_NaN();
				}
			}
		}
	}



	// ---------------------------------------------------------------------
	// to prune SNPs

	struct TSNP
	{
		int idx, pos_bp;
		vector<UInt8> genobuf;
		TSNP(int n=0): genobuf(n) {}
		TSNP(int _idx, int _pos, int n=0): genobuf(n) { idx = _idx; pos_bp = _pos; }
	};

	// to compute pair LD
	static double _CalcLD(const UInt8 *snp1, const UInt8 *snp2)
	{
		switch (LD_Method)
		{
			case 1:
				return PairComposite(snp1, snp2);
			case 2:
				return PairR(snp1, snp2);
			case 3:
				return PairDPrime(snp1, snp2);
			case 4:
				return PairCorr(snp1, snp2);
		}
		return conf_F64_NaN();
	}

	void calcLDpruning(int StartIdx, int *pos_bp, int slide_max_bp, int slide_max_n,
		const double LD_threshold, bool *out_SNP)
	{
		// initial variables
		nPackedSamp = (MCWorkingGeno.Space.SampleNum() % 4 > 0) ?
			(MCWorkingGeno.Space.SampleNum()/4 + 1) : (MCWorkingGeno.Space.SampleNum()/4);
		list<TSNP> ListGeno;
		out_SNP[StartIdx] = true;
		auto_ptr<UInt8> buf(new UInt8[nPackedSamp]);

		// -----------------------------------------------------
		// increasing searching

		CdBufSpace BufSNP(MCWorkingGeno.Space, true, CdBufSpace::acInc);
		ListGeno.push_back(TSNP(StartIdx, pos_bp[StartIdx], nPackedSamp));
		BufSNP.ReadPackedGeno(StartIdx, &(ListGeno.back().genobuf[0]));

		// for-loop, increasing
		for (int i=StartIdx+1; i < MCWorkingGeno.Space.SNPNum(); i++)
		{
			// load genotypes
			BufSNP.ReadPackedGeno(i, buf.get());
			// detect LD
			int TotalCnt = 0, ValidCnt = 0;
			list<TSNP>::iterator it;
			for (it=ListGeno.begin(); it != ListGeno.end(); TotalCnt++)
			{
				// check whether it is in the sliding window
				if ((abs(i - it->idx) <= slide_max_n) &&
					(abs(pos_bp[i] - it->pos_bp) <= slide_max_bp))
				{
					if (fabs(_CalcLD(&(it->genobuf[0]), buf.get())) <= LD_threshold)
						ValidCnt ++;
					it++;
				} else {
					ValidCnt ++;
					// delete it
					list<TSNP>::iterator tmp_it = it;
					it++;
					ListGeno.erase(tmp_it);
				}
			}
			// handle
			out_SNP[i] = (ValidCnt == TotalCnt);
			if (out_SNP[i])
			{
				ListGeno.push_back(TSNP(i, pos_bp[i], nPackedSamp));
				memmove(&(ListGeno.back().genobuf[0]), buf.get(), nPackedSamp);
			}
		}

		// -----------------------------------------------------
		// decreasing searching

		ListGeno.clear();
		for (int i=StartIdx; i < MCWorkingGeno.Space.SNPNum(); i++)
		{
			if (out_SNP[i])
			{
				// check whether it is in the sliding window
				if ((abs(i - StartIdx) <= slide_max_n) &&
					(abs(pos_bp[i] - pos_bp[StartIdx]) <= slide_max_bp))
				{
					ListGeno.push_back(TSNP(i, pos_bp[i], nPackedSamp));
					BufSNP.ReadPackedGeno(i, &(ListGeno.back().genobuf[0]));
				} else
					break;
			}
		}

		BufSNP.SetAccessFlag(CdBufSpace::acDec);
		// for-loop, descreasing
		for (int i=StartIdx-1; i >= 0; i--)
		{
			// load genotypes
			BufSNP.ReadPackedGeno(i, buf.get());
			// detect LD
			int TotalCnt = 0, ValidCnt = 0;
			list<TSNP>::iterator it;
			for (it=ListGeno.begin(); it != ListGeno.end(); TotalCnt++)
			{
				// check whether it is in the sliding window
				if ((abs(i - it->idx) <= slide_max_n) &&
					(abs(pos_bp[i] - it->pos_bp) <= slide_max_bp))
				{
					if (fabs(_CalcLD(&(it->genobuf[0]), buf.get())) <= LD_threshold)
						ValidCnt ++;
					it++;
				} else {
					ValidCnt ++;
					// delete it
					list<TSNP>::iterator tmp_it = it;
					it++;
					ListGeno.erase(tmp_it);
				}
			}
			// handle
			out_SNP[i] = (ValidCnt == TotalCnt);
			if (out_SNP[i])
			{
				ListGeno.push_front(TSNP(i, pos_bp[i], nPackedSamp));
				memmove(&(ListGeno.front().genobuf[0]), buf.get(), nPackedSamp);
			}
		}
	}
}


#endif  /* _FuncLD_H_ */
