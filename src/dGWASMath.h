// ===========================================================
//     _/_/_/   _/_/_/  _/_/_/_/    _/_/_/_/  _/_/_/   _/_/_/
//      _/    _/       _/             _/    _/    _/   _/   _/
//     _/    _/       _/_/_/_/       _/    _/    _/   _/_/_/
//    _/    _/       _/             _/    _/    _/   _/
// _/_/_/   _/_/_/  _/_/_/_/_/     _/     _/_/_/   _/_/
// ===========================================================
// Name        : dGWASMath
// Author      : Xiuwen Zheng
// Version     : 1.0.0.0
// Copyright   : Xiuwen Zheng (GPL v3.0)
// Created     : 04/22/2012
// Description :
// ===========================================================


#ifndef _dGWASMath_H_
#define _dGWASMath_H_

#include <limits>
#include <memory>
#include <iostream>
#include <cmath>


namespace GWAS_Math
{
	using namespace std;

	template<typename Tx> struct TWrapFloat {};

	template<> struct TWrapFloat<float>
		{ typedef double Type; };

	template<> struct TWrapFloat<double>
		{ typedef double Type; };

	template<> struct TWrapFloat<long double>
		{ typedef long double Type; };


	// Nelder-Mead Simplex algorithm
	// J.A. Nelder and R. Mead, A simplex method for function minimization,
	//   Computer Journal vol. 7 (1965), 308¨C315.
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
				--nfunk;
		}
	}



	/*  BFGS variable-metric method, based on Pascal code
	in J.C. Nash, `Compact Numerical Methods for Computers', 2nd edition,
	converted by p2c then re-crafted by B.D. Ripley */

	template<typename tfloat, int n,
		typename Functor1, typename Functor2>
	void BFGSMin(tfloat *b, Functor1 func, Functor2 fngr, void *FuncData,
		tfloat &OutFmin, int &Outfncount, int &Outgrcount, bool &fail,
		int maxit, tfloat abstol, tfloat reltol)
	{
		const tfloat stepredn = 0.2;
		const tfloat acctol = 0.0001;
		const tfloat reltest = 10.0;

		bool accpoint, enough;
		tfloat g[n], X[n], c[n];
		tfloat t[n], B[n][n];
		tfloat f, gradproj;
		tfloat s, steplength;
		tfloat D1, D2;
		int count, i, j, ilast, iter = 0;

		if (maxit <= 0)
		{
			fail = false;
			OutFmin = func(b, FuncData);
			Outfncount = Outgrcount = 0;
			return;
		}

		f = func(b, FuncData);
		OutFmin = f;
		Outfncount = Outgrcount = 1;
		fngr(b, g, FuncData);
		iter++;
		ilast = Outgrcount;

		do {
			if (ilast == Outgrcount)
				for (i = 0; i < n; i++)
				{
					for (j = 0; j < i; j++) B[i][j] = 0.0;
					B[i][i] = 1.0;
				}
			for (i = 0; i < n; i++)
				{ X[i] = b[i]; c[i] = g[i]; }
			gradproj = 0.0;
			for (i = 0; i < n; i++)
			{
				s = 0.0;
				for (j = 0; j <= i; j++) s -= B[i][j] * g[j];
				for (j = i + 1; j < n; j++) s -= B[j][i] * g[j];
				t[i] = s;
				gradproj += s * g[i];
			}

			if (gradproj < 0.0)
			{	// search direction is downhill
				steplength = 1.0;
				accpoint = true;
				do {
					count = 0;
					for (i = 0; i < n; i++)
					{
						b[i] = X[i] + steplength * t[i];
						if (reltest + X[i] == reltest + b[i]) /* no change */
							count++;
					}
					if (count < n)
					{
						f = func(b, FuncData);
						Outfncount++;
						accpoint = (f <= OutFmin + gradproj*steplength*acctol);
						if (!accpoint)
							steplength *= stepredn;
					}
				} while (!(count == n || accpoint));
				enough = (f > abstol) &&
					fabs(f - OutFmin) > reltol * (fabs(OutFmin) + reltol);
//					enough = fabs(f - OutFmin) > abstol;
				// stop if value if small or if relative change is low
				if (!enough)
					{ count = n; OutFmin = f; }
				if (count < n)
				{	// making progress
					OutFmin = f;
					fngr(b, g, FuncData);
					Outgrcount++;
					iter++;
					D1 = 0.0;
					for (i = 0; i < n; i++)
					{
						t[i] = steplength * t[i];
						c[i] = g[i] - c[i];
						D1 += t[i] * c[i];
					}
					if (D1 > 0)
					{
						D2 = 0.0;
						for (i = 0; i < n; i++)
						{
							s = 0.0;
							for (j = 0; j <= i; j++)
								s += B[i][j] * c[j];
							for (j = i + 1; j < n; j++)
								s += B[j][i] * c[j];
							X[i] = s;
							D2 += s * c[i];
						}
						D2 = 1.0 + D2 / D1;
						for (i = 0; i < n; i++)
						{
							for (j = 0; j <= i; j++)
								B[i][j] += (D2 * t[i] * t[j]
									- X[i] * t[j] - t[i] * X[j]) / D1;
						}
					} else {	/* D1 < 0 */
						ilast = Outgrcount;
					}
				} else {	/* no progress */
					if (ilast < Outgrcount)
						{ count = 0; ilast = Outgrcount; }
				}
			} else {		/* uphill search */
				count = 0;
				if (ilast == Outgrcount)
					count = n;
				else
					ilast = Outgrcount;
				/* Resets unless has just been reset */
			}
			if (iter >= maxit)
				break;
			if (Outgrcount - ilast > 2 * n)
				ilast = Outgrcount;	/* periodic restart */
		} while (count != n || ilast != Outgrcount);

		fail = (iter < maxit) ? false : true;
	}


	enum CGMethod { cgFletcherReeves, cgPolakRibiere, cgBealeSorenson };

	template<typename tfloat, int n,
		typename Functor1, typename Functor2>
	void cgmin(tfloat *StartX, tfloat *X,
		Functor1 fminfn, Functor2 fmingr, void *FuncData,
		CGMethod type, tfloat &OutFmin,
		int &fncount, int &grcount, bool &fail,
		int maxit, tfloat intol)
	{
		const tfloat stepredn = 0.2;
		const tfloat acctol = 0.0001;
		const tfloat reltest = 10.0;

		bool accpoint;
		tfloat g[n];
		tfloat c[n], t[n];
		int count, cycle, cyclimit, i;
		tfloat f;
		tfloat G1, G2, G3, gradproj;
		tfloat newstep, oldstep, setstep, steplength=1.0;
		tfloat tol;

		if (maxit <= 0)
		{
			OutFmin = fminfn(StartX, FuncData);
			fncount = grcount = 0;
			fail = false;
			return;
		}

		setstep = 1.7;
		fail = false;
		cyclimit = n;
		tol = intol * n * sqrt(intol);

		f = fminfn(StartX, FuncData);
		OutFmin = f;
		fncount = 1; grcount = 0;
		do {
			for (i = 0; i < n; i++)
				t[i] = c[i] = 0.0;
			cycle = 0;
			oldstep = 1.0;
			count = 0;
			do {
				cycle++;
				grcount++;
				if (grcount > maxit)
				{
					fail = 1;
					return;
				}
				fmingr(StartX, g, FuncData);
				G1 = G2 = 0.0;
				for (i = 0; i < n; i++)
				{
					X[i] = StartX[i];
					switch (type) {
						case cgFletcherReeves: // Fletcher-Reeves
							G1 += g[i] * g[i];
							G2 += c[i] * c[i];
							break;
						case cgPolakRibiere: // Polak-Ribiere
							G1 += g[i] * (g[i] - c[i]);
							G2 += c[i] * c[i];
							break;
						case cgBealeSorenson:  // Beale-Sorenson
							G1 += g[i] * (g[i] - c[i]);
							G2 += t[i] * (g[i] - c[i]);
							break;
					}
					c[i] = g[i];
				}
				if (G1 > tol)
				{
					G3 = (G2 > 0.0) ? (G1 / G2) : 1.0;
					gradproj = 0.0;
					for (i = 0; i < n; i++)
					{
						t[i] = t[i] * G3 - g[i];
						gradproj += t[i] * g[i];
					}
					steplength = oldstep;

					accpoint = false;
					do {
						count = 0;
						for (i = 0; i < n; i++)
						{
							StartX[i] = X[i] + steplength * t[i];
							if (reltest + X[i] == reltest + StartX[i])
								count++;
						}
						if (count < n)
						{	/* point is different */
							f = fminfn(StartX, FuncData);
							fncount++;
							accpoint =
								f <= OutFmin + gradproj * steplength * acctol;

							if (!accpoint)
								steplength *= stepredn;
							else
								OutFmin = f; /* we improved, so update value */
						}
					} while (!(count == n || accpoint));

					if (count < n)
					{
						newstep = 2 * (f - OutFmin - gradproj*steplength);
						if (newstep > 0)
						{
							newstep = -(gradproj * steplength * steplength / newstep);
							for (i = 0; i < n; i++)
								StartX[i] = X[i] + newstep * t[i];
							OutFmin = f;
							f = fminfn(StartX, FuncData);
							fncount++;
							if (f < OutFmin)
								OutFmin = f;
							else { /* reset StartX to match lowest point */
								for (i = 0; i < n; i++)
									StartX[i] = X[i] + steplength*t[i];
							}
						}
					}
				}
				oldstep = setstep * steplength;
				if (oldstep > 1.0)
					oldstep = 1.0;
			} while ((count != n) && (G1 > tol) && (cycle != cyclimit));

		} while ((cycle != 1) || ((count != n) && (G1 > tol)));
	}

}

#endif /* _dGWASMath_H_ */
