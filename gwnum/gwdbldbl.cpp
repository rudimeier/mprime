/**************************************************************
 *
 *  gwdbldbl.cpp
 *
 *  This file contains all the gwnum initialization routines that require
 *  extended precision floating point.  We want to initialize our sin/cos
 *  and FFT weight arrays with doubles that are as accurate as possible.
 *
 *  This is the only C++ routine in the gwnum library.  Since gwnum is
 *  a C based library, we declare all routines here as extern "C".
 * 
 *  Copyright 2005-2009 Mersenne Research, Inc.  All rights reserved.
 *
 **************************************************************/

/* Include files */

#include "gwdbldbl.h"

/* Pick which doubledouble package we will use. */

#define QD
//#define KEITH_BRIGGS

/* Turn on #define that will disable extended precision floating point */
/* registers, as they wreak havoc with double-double library routines. */

#ifndef X86_64
#define x86
#endif

/* Use Hida, Li & Bailey's QD doubledouble C++ package. */

#ifdef QD
#include "dd.cc"
#endif

/* Use Keith Briggs' doubledouble C++ package.  I find the QD package */
/* a better choice. */

#ifdef KEITH_BRIGGS
#define DD_INLINE
#include "doubledouble.cc"
#include "math.cc"
#define	dd_real doubledouble
#define	_2pi TwoPi
#define	_log2 Log2
#endif

/* Structure for "global" data that gwfft_weight_setup and passes */
/* to several gwdbldbl routines */
	
struct gwdbldbl_constants {
	dd_real	gw__b;
	dd_real	gw__logb;
	dd_real	gw__num_b_per_word;
	int	gw__c_is_one;
	dd_real gw__logb_abs_c_div_fftlen;
	dd_real gw__fftlen_inverse;
	dd_real gw__over_fftlen;
	double	gwdbl__b;
	double	gwdbl__b_inverse;
	double	gwdbl__num_b_per_word;
	double	gwdbl__logb_abs_c_div_fftlen;
#ifdef VERY_SLOPPY
	unsigned long last_sloppy_j;
	double	last_sloppy_result;
	double	fast_sloppy_multiplier;
#endif
	unsigned long last_inv_sloppy_j;
	double	last_inv_sloppy_result;
	double	fast_inv_sloppy_multiplier;
};

/* Macro routines below use to type the cast untyped data pointer */
/* input argument. */

#define dd_data		((struct gwdbldbl_constants *) dd_data_arg)

/* Now write all the routines that use the dd_real package. */


/* Utility routine to compute many of the constants used by the */
/* assembly language code */

extern "C"
void gwasm_constants (
	double	*asm_values)
{
	dd_real arg, sine1, cosine1, sine2, cosine2, sine3, cosine3;
#define P951			asm_values[0]
#define P618			asm_values[1]
#define P309			asm_values[2]
#define M262			asm_values[3]
#define P588			asm_values[4]
#define M162			asm_values[5]
#define M809			asm_values[6]
#define M382			asm_values[7]
#define P866			asm_values[8]
#define P433			asm_values[9]
#define P577			asm_values[10]
#define P975			asm_values[11]
#define P445			asm_values[12]
#define P180			asm_values[13]
#define P623			asm_values[14]
#define M358			asm_values[15]
#define P404			asm_values[16]
#define M223			asm_values[17]
#define M901			asm_values[18]
#define M691			asm_values[19]

/* Do some initial setup */
	
	x86_FIX

/* Initialize the five_reals sine-cosine data. */
/* NOTE: When computing cosine / sine, divide by the 64-bit sine not the */
/* extra precision sine since macros will multiply by the 64-bit sine. */

	arg = dd_real::_2pi / 5.0;		// 2*PI * 1 / 5
	sincos (arg, sine1, cosine1);

	arg = arg * 2.0;			// 2*PI * 2 / 5
	sincos (arg, sine2, cosine2);

	P951 = double (sine1);
	P618 = double (sine2 / P951);		// 0.588 / 0.951

	P309 = double (cosine1);
	M262 = double (cosine2 / P309);		// -0.809 / 0.309

	P588 = double (sine2);
	M162 = double (-sine1 / P588);		// -0.951 / 0.588

	M809 = double (cosine2);
	M382 = double (cosine1 / M809);		// 0.309 / -0.809

/* Initialize the six_reals sine-cosine data. */

	arg = dd_real::_2pi / 3.0;		// 2*PI / 3
	sine1 = sin (arg);			// Compute sine (0.866)

	P866 = double (sine1);
	P433 = double (sine1 * 0.5);		// 0.5 * P866
	P577 = double (dd_real (0.5) / sine1);	// 0.5 / 0.866

/* Initialize the seven_reals sine-cosine data. */

	arg = dd_real::_2pi / 7.0;		// 2*PI * 1 / 7
	sincos (arg, sine1, cosine1);		// cosine (0.623), sine (0.782)

	arg = arg * 2.0;			// 2*PI * 2 / 7
	sincos (arg, sine2, cosine2);		// cosine (-.223), sine (0.975)

	arg = arg * 1.5;			// 2*PI * 3 / 7
	sincos (arg, sine3, cosine3);		// cosine (-.901), sine (0.434)
		
	P975 = double (sine2);
	P445 = double (sine3 / P975);		// 0.434 / 0.975
	P180 = double (sine1 / P975 / P445);	// 0.782 / (0.975 * 0.445)

	P623 = double (cosine1);
	M358 = double (cosine2 / P623);		// -0.223 / 0.623
	P404 = double (cosine3 / P623 / M358);	// -.901 / (.623 * -.358)

	M223 = double (cosine2);
	M901 = double (cosine3);
	M691 = double (cosine1 / M901);		// 0.623 / -0.901

	END_x86_FIX
}

// Utility routine to compute a sin/cos premultiplier or a set of 3
// sine-cosine values.
// This is used during setup, formerly written in assembly language to take
// advantage of the extra precision in the FPU's 80-bit registers.
// NOTE: When computing cosine / sine, divide by the 64-bit sine
// not the 80-bit sine since macros will multiply by the 64-bit sine.

extern "C"
void gwsincos (
	unsigned long x,
	unsigned long N,
	double	*results)
{
	dd_real arg, sine, cosine;

	x86_FIX
	arg = dd_real::_2pi * (double) x / (double) N;
	sincos (arg, sine, cosine);
	sine += 1E-200;			/* Protect against divide by zero */
	results[0] = sine;
	results[1] = cosine / results[0];
	END_x86_FIX
}

extern "C"
void gwsincos3 (
	unsigned long x,
	unsigned long N,
	double	*results)
{
	dd_real arg1, arg2, arg3, sine, cosine;

	x86_FIX
	arg1 = dd_real::_2pi * (double) x / (double) N;
	sincos (arg1, sine, cosine);
	sine += 1E-200;			/* Protect against divide by zero */
	results[0] = sine;
	results[1] = cosine / results[0];
	arg2 = arg1 + arg1;
	sincos (arg2, sine, cosine);
	sine += 1E-200;			/* Protect against divide by zero */
	results[2] = sine;
	results[3] = cosine / results[2];
	arg3 = arg2 + arg1;
	sincos (arg3, sine, cosine);
	sine += 1E-200;			/* Protect against divide by zero */
	results[4] = sine;
	results[5] = cosine / results[4];
	END_x86_FIX
}


//
// Utility routines to compute fft weights
//
// The FFT weight for the j-th FFT word doing a b^n+c weighted transform is
//	b ^ (ceil (j*n/FFTLEN) - j*n/FFTLEN)   *    abs(c) ^ j/FFTLEN
//

extern "C"
void *gwdbldbl_data_alloc (void)
{
	return (malloc (sizeof (struct gwdbldbl_constants)));
}

extern "C"
void gwfft_weight_setup (
	void	*dd_data_arg,
	int	zero_pad,
	double	k,
	unsigned long b,
	unsigned long n,
	signed long c,
	unsigned long fftlen)
{
	x86_FIX
	dd_data->gw__b = dd_real ((double) b);
	dd_data->gw__logb = log (dd_real ((double) b));
	dd_data->gw__fftlen_inverse = dd_real (1.0) / dd_real ((double) fftlen);
	if (zero_pad) {
		dd_data->gw__num_b_per_word =
			dd_real ((double) (n + n)) * dd_data->gw__fftlen_inverse;
		dd_data->gw__c_is_one = TRUE;
		dd_data->gw__over_fftlen =
			dd_real (2.0) * dd_data->gw__fftlen_inverse;
	} else {
		dd_data->gw__num_b_per_word =
			(dd_real ((double) n) + log (dd_real (k)) / dd_data->gw__logb) *
			dd_data->gw__fftlen_inverse;
		dd_data->gw__c_is_one = (abs ((int) c) == 1);
		dd_data->gw__logb_abs_c_div_fftlen =
			log (dd_real (abs ((int) c))) / dd_data->gw__logb *
			dd_data->gw__fftlen_inverse;
		dd_data->gw__over_fftlen =
			dd_real (k * 2.0) * dd_data->gw__fftlen_inverse;
	}
	dd_data->gwdbl__b = (double) b;
	dd_data->gwdbl__b_inverse = 1.0 / (double) b;
	dd_data->gwdbl__num_b_per_word = (double) dd_data->gw__num_b_per_word;
	dd_data->gwdbl__logb_abs_c_div_fftlen = (double) dd_data->gw__logb_abs_c_div_fftlen;
#ifdef VERY_SLOPPY
	dd_data->last_sloppy_j = 0;
	dd_data->last_sloppy_result = 1.0;
	dd_data->fast_sloppy_multiplier = gwfft_weight (dd_data, 1);
#endif
	dd_data->last_inv_sloppy_j = 0;
	dd_data->last_inv_sloppy_result = 1.0;
	dd_data->fast_inv_sloppy_multiplier = gwfft_weight_inverse (dd_data, 1);
	END_x86_FIX
}

extern "C"
double gwfft_weight (
	void	*dd_data_arg,
	unsigned long j)
{
	dd_real temp, bpower, result;

	x86_FIX
	temp = dd_real ((double) j) * dd_data->gw__num_b_per_word;
	bpower = ceil (temp) - temp;
	if (! dd_data->gw__c_is_one)
		bpower += dd_data->gw__logb_abs_c_div_fftlen * dd_real ((double) j);
	result = exp (dd_data->gw__logb * bpower);
	END_x86_FIX
	return (double (result));
}

// Like the above, but faster and does not guarantee quite as much accuracy.

extern "C"
double gwfft_weight_sloppy (
	void	*dd_data_arg,
	unsigned long j)
{
	dd_real temp, bpower;

	x86_FIX
	temp = dd_real ((double) j) * dd_data->gw__num_b_per_word;
	bpower = ceil (temp) - temp;
	if (! dd_data->gw__c_is_one)
		bpower += dd_data->gw__logb_abs_c_div_fftlen * dd_real ((double) j);
	END_x86_FIX
	return (pow (dd_data->gwdbl__b, double (bpower)));

// We cannot be very sloppy as these weights are used to set FFT data values
// when reading save files.  If we are too sloppy we get roundoff errors > 0.5.

#ifdef VERY_SLOPPY
	double	temp, bpower, result;

// Our sequential sloppy optimizations won't work if abs(c) is not one.
// This is because the result is not in the range 1.0 to 2.0.

	if (! dd_data->gw__c_is_one) {
		temp = (double) j * dd_data->gwdbl__num_b_per_word;
		bpower = ceil (temp) - temp;
		if (bpower < 0.001 || bpower > 0.999)
			return (gwfft_weight (dd_data_arg, j));
		bpower += dd_data->gwdbl__logb_abs_c_div_fftlen * (double) j;
		return (pow (dd_data->gwdbl__b, bpower));
	}

// Compute weight from previous weight, but don't do too many of
// these in a row as floating point roundoff errors will accumulate

	if (j == dd_data->last_sloppy_j + 1 && (j & 0x7F)) {
		result = dd_data->last_sloppy_result * dd_data->fast_sloppy_multiplier;
		if (result >= dd_data->gwdbl__b) result = result * dd_data->gwdbl__b_inverse;
	}

// Use a slower sloppy technique

	else {
		temp = (double) j * dd_data->gwdbl__num_b_per_word;
		bpower = ceil (temp) - temp;
		result = pow (dd_data->gwdbl__b, double (bpower));
	}

// Just to be safe, if result is at all close to the boundaries return
// the carefully computed weight.

	if (result < 1.00001 || result > dd_data->gwdbl__b - 0.00001)
		result = gwfft_weight (dd_data_arg, j);
			 
// Save the result for faster sequential sloppy calls

	dd_data->last_sloppy_j  = j;
	dd_data->last_sloppy_result = result;
	return (result);
#endif
}

// Compute the inverse of the fft weight

extern "C"
double gwfft_weight_inverse (
	void	*dd_data_arg,
	unsigned long j)
{
	dd_real temp, bpower, result;

	x86_FIX
	temp = dd_real ((double) j) * dd_data->gw__num_b_per_word;
	bpower = ceil (temp) - temp;
	if (! dd_data->gw__c_is_one)
		bpower += dd_data->gw__logb_abs_c_div_fftlen * dd_real ((double) j);
	result = exp (dd_data->gw__logb * -bpower);
	END_x86_FIX
	return (double (result));
}

// Like the above, but faster and does not guarantee quite as much accuracy.
// We can be very sloppy as these weights are used to read FFT data values
// when writing save files.

extern "C"
double gwfft_weight_inverse_sloppy (
	void	*dd_data_arg,
	unsigned long j)
{
	double	temp, bpower, result;

// Our sequential sloppy optimizations won't work if abs(c) is not one.
// This is because the result is not in the range 0.5 to 1.0.

	if (! dd_data->gw__c_is_one) {
		temp = (double) j * dd_data->gwdbl__num_b_per_word;
		bpower = ceil (temp) - temp;
		if (bpower < 0.001 || bpower > 0.999)
			return (gwfft_weight_inverse (dd_data_arg, j));
		bpower += dd_data->gwdbl__logb_abs_c_div_fftlen * (double) j;
		return (pow (dd_data->gwdbl__b, - bpower));
	}

// Compute weight from previous weight, but don't do too many of
// these in a row as floating point roundoff errors will accumulate

	if (j == dd_data->last_inv_sloppy_j + 1 && (j & 0x7F)) {
		result = dd_data->last_inv_sloppy_result * dd_data->fast_inv_sloppy_multiplier;
		if (result <= dd_data->gwdbl__b_inverse) result = result * dd_data->gwdbl__b;
	}

// Use a slower sloppy technique

	else {
		temp = (double) j * dd_data->gwdbl__num_b_per_word;
		bpower = ceil (temp) - temp;
		result = pow (dd_data->gwdbl__b, - bpower);
	}

// Just to be safe, if result is at all close to the boundaries return
// the carefully computed weight.

	if (result < dd_data->gwdbl__b_inverse + 0.00001 || result > 0.99999)
		result = gwfft_weight_inverse (dd_data_arg, j);

// Save the result for faster sequential sloppy calls

	dd_data->last_inv_sloppy_j  = j;
	dd_data->last_inv_sloppy_result = result;
	return (result);
}

// This computes the inverse FFT weight multiplied by the appropriate constant
// to produce an integer during an FFT multiply's normalize stage.  This
// constant is 2/FFTLEN for a zero-padded FFT and k*2/FFTLEN for a
// non-zero-padded FFT.

extern "C"
double gwfft_weight_inverse_over_fftlen (
	void	*dd_data_arg,
	unsigned long j)
{
	dd_real temp, bpower, result;

	x86_FIX
	temp = dd_real ((double) j) * dd_data->gw__num_b_per_word;
	bpower = ceil (temp) - temp;
	if (! dd_data->gw__c_is_one)
		bpower += dd_data->gw__logb_abs_c_div_fftlen * dd_real ((double) j);
	result = exp (dd_data->gw__logb * -bpower) * dd_data->gw__over_fftlen;
	END_x86_FIX
	return (double (result));
}

// This computes the three FFT weights in one call.  It is faster than
// calling the above individually.

extern "C"
void gwfft_weights3 (
	void	*dd_data_arg,
	unsigned long j,
	double	*fft_weight,
	double	*fft_weight_inverse,
	double	*fft_weight_inverse_over_fftlen)
{
	dd_real temp, bpower, weight;

	x86_FIX
	temp = dd_real ((double) j) * dd_data->gw__num_b_per_word;
	bpower = ceil (temp) - temp;
	if (! dd_data->gw__c_is_one)
		bpower += dd_data->gw__logb_abs_c_div_fftlen * dd_real ((double) j);
	weight = exp (dd_data->gw__logb * bpower);
	*fft_weight = double (weight);
	weight = dd_real (1.0) / weight;
	if (fft_weight_inverse != NULL)
		*fft_weight_inverse = double (weight);
	if (fft_weight_inverse_over_fftlen != NULL)
		*fft_weight_inverse_over_fftlen = double (weight * dd_data->gw__over_fftlen);
	END_x86_FIX
}

// Returns logb(fft_weight).  This is used in determining the FFT weight
// fudge factor in two-pass FFTs.  This is much faster than computing the
// fft_weight because it eliminates a call to the double-double exp routine.

extern "C"
double gwfft_weight_exponent (
	void	*dd_data_arg,
	unsigned long j)
{
	double	tempdbl, bpowerdbl;
	dd_real temp, bpower;

// For speed, try this with plain old doubles first

	if (j == 0) return (0);
	tempdbl = (double) j * dd_data->gwdbl__num_b_per_word;
	bpowerdbl = ceil (tempdbl) - tempdbl;
	if (bpowerdbl > 0.001 && bpowerdbl < 0.999) return (bpowerdbl);

// If at all uncertain of the result, use doubledoubles to do the calculation

	x86_FIX
	temp = dd_real ((double) j) * dd_data->gw__num_b_per_word;
	bpower = ceil (temp) - temp;
	END_x86_FIX
	return (double (bpower));
}

//
// Utility routine to compute fft base for j-th fft word
//
// The FFT base for the j-th FFT word doing a b^n+c weighted transform is
//	ceil (j*n/FFTLEN)
// This routine returns ceil (j*n/FFTLEN) taking great care to return a
// value accurate to 53 bits.  This is important when j*n is really close to
// a multiple of FFTLEN (admittedly quite rare).  It would be really bad if
// rounding differences caused this routine to compute ceil (j*n/FFTLEN)
// differently than the weighting functions.
//

extern "C"
unsigned long gwfft_base (
	void	*dd_data_arg,
	unsigned long j)
{
	double	tempdbl, ceildbl, diffdbl;
	dd_real temp;
	unsigned long bpower;

// For speed, try this with plain old doubles first

	if (j == 0) return (0);
	tempdbl = (double) j * dd_data->gwdbl__num_b_per_word;
	ceildbl = ceil (tempdbl);
	diffdbl = ceildbl - tempdbl;
	if (diffdbl > 0.001 && diffdbl < 0.999)
		return ((unsigned long) ceildbl);

// If at all uncertain of the result, use doubledoubles to do the calculation

	x86_FIX
	temp = dd_real ((double) j) * dd_data->gw__num_b_per_word;
	bpower = (int) ceil (temp);
	END_x86_FIX
	return (bpower);
}
