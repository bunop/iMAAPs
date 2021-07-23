/*
 * convolution.h
 *
 *  Created on: 2014.12.24
 *      Author: Yuan Kai
 *
 *  This head file contains two inline function. Both of them are used to calculate convolution.
 *
 */

#ifndef CONVOLUTION_H_
#define CONVOLUTION_H_

# include "fftw3.h"

inline void self_convolution(fftw_plan *plan, fftw_complex *out, fftw_complex *t, int n, double factor)
/*
 * Calculation for self convolution
 *
 * plan    FFTW plan pointer. Only one plan is needed for this function.
 * out     A pointer point to the output results of this function. At least n * sizeof ( fftw_complex ) space is needed.
 *         No boundary check for out.
 * t       A pointer point to a temporary space at least n * sizeof ( fftw_complex ).
 *         No boundary check for t.
 * n       The scale of FFTW plan. If n is equal to powers of 2, the calculation would be very fast.
 *         This n should be equal to the scale of FFTW plan.
 * factor  A coefficient for this convolution calculation.
 *
 */
{
	fftw_execute(plan[0]);

	for(int i=0;i<n;++i)
	{
		out[i][0]+=(t[i][0]*t[i][0]+t[i][1]*t[i][1])*factor;
	}
}

inline void convolution(fftw_plan *plan, fftw_complex *out, fftw_complex *t1, fftw_complex *t2, int n, double factor)
/*
 * Calculation for convolution of two plans
 *
 * plan    FFTW plan pointer. Two plans are needed for this function.
 * out     A pointer point to the output results of this function. At least n * sizeof ( fftw_complex ) space is needed.
 *         No boundary check for out.
 * t1,t2   Pointers point to temporary space at least n * sizeof ( fftw_complex ) for each.
 *         No boundary check for t1 and t2.
 * n       The scale of FFTW plan. If n is equal to powers of 2, the calculation would be very fast.
 *         This n should be equal to the scale of FFTW plan.
 * factor  A coefficient for this convolution calculation.
 *
 */
{
	fftw_execute(plan[0]);
	fftw_execute(plan[1]);

	for(int i=0;i<n;++i)
	{
		out[i][0]+=(t1[i][0]*t2[i][0]+t1[i][1]*t2[i][1])*factor;
		out[i][1]+=(t1[i][0]*t2[i][1]-t1[i][1]*t2[i][0])*factor;
	}
}

#endif /* CONVOLUTION_H_ */
