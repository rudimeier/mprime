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
 *  Copyright 2005-2006 Mersenne Research, Inc.  All rights reserved.
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
	dd_real	gw__bits_per_word;
	int	gw__c_is_one;
	dd_real gw__log2_abs_c_div_fftlen;
	dd_real gw__fftlen_inverse;
	dd_real gw__over_fftlen;
	double	gwdbl__bits_per_word;
	double	gwdbl__log2_abs_c_div_fftlen;
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
// The FFT weight for the j-th FFT word doing a 2^q+c weighted transform is
//	2 ^ (ceil (j*q/FFTLEN) - j*q/FFTLEN)   *    abs(c) ^ j/FFTLEN
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
	unsigned long n,
	signed long c,
	unsigned long fftlen)
{
	x86_FIX
	dd_data->gw__fftlen_inverse = dd_real (1.0) / dd_real ((double) fftlen);
	if (zero_pad) {
		dd_data->gw__bits_per_word =
			dd_real ((double) (n+n)) * dd_data->gw__fftlen_inverse;
		dd_data->gw__c_is_one = TRUE;
		dd_data->gw__over_fftlen =
			dd_real (2.0) * dd_data->gw__fftlen_inverse;
	} else {
		dd_data->gw__bits_per_word =
			(dd_real ((double) n) +
			 log (dd_real (k)) / dd_real::_log2) *
			dd_data->gw__fftlen_inverse;
		dd_data->gw__c_is_one = (abs ((int) c) == 1);
		dd_data->gw__log2_abs_c_div_fftlen =
			log (dd_real (abs ((int) c))) / dd_real::_log2 *
			dd_data->gw__fftlen_inverse;
		dd_data->gw__over_fftlen =
			dd_real (k * 2.0) * dd_data->gw__fftlen_inverse;
	}
	dd_data->gwdbl__bits_per_word = (double) dd_data->gw__bits_per_word;
	dd_data->gwdbl__log2_abs_c_div_fftlen = (double) dd_data->gw__log2_abs_c_div_fftlen;
	END_x86_FIX
}

extern "C"
double gwfft_weight (
	void	*dd_data_arg,
	unsigned long j)
{
	dd_real temp, twopow, result;

	x86_FIX
	temp = dd_real ((double) j) * dd_data->gw__bits_per_word;
	twopow = ceil (temp) - temp;
	if (! dd_data->gw__c_is_one)
		twopow += dd_data->gw__log2_abs_c_div_fftlen * dd_real ((double) j);
	result = exp (dd_real::_log2 * twopow);
	END_x86_FIX
	return (double (result));
}

// Like the above, but faster and does not guarantee quite as much accuracy.

extern "C"
double gwfft_weight_sloppy (
	void	*dd_data_arg,
	unsigned long j)
{
	dd_real temp, twopow;

	x86_FIX
	temp = dd_real ((double) j) * dd_data->gw__bits_per_word;
	twopow = ceil (temp) - temp;
	if (! dd_data->gw__c_is_one)
		twopow += dd_data->gw__log2_abs_c_div_fftlen * dd_real ((double) j);
	END_x86_FIX
	return (pow (2.0, double (twopow)));
}

// Compute the inverse of the fft weight

extern "C"
double gwfft_weight_inverse (
	void	*dd_data_arg,
	unsigned long j)
{
	dd_real temp, twopow, result;

	x86_FIX
	temp = dd_real ((double) j) * dd_data->gw__bits_per_word;
	twopow = ceil (temp) - temp;
	if (! dd_data->gw__c_is_one)
		twopow += dd_data->gw__log2_abs_c_div_fftlen * dd_real ((double) j);
	result = exp (dd_real::_log2 * -twopow);
	END_x86_FIX
	return (double (result));
}

// Like the above, but faster and does not guarantee quite as much accuracy.

extern "C"
double gwfft_weight_inverse_sloppy (
	void	*dd_data_arg,
	unsigned long j)
{
	double	tempdbl, twopowdbl;

	tempdbl = (double) j * dd_data->gwdbl__bits_per_word;
	twopowdbl = ceil (tempdbl) - tempdbl;
	if (twopowdbl < 0.001 || twopowdbl > 0.999) {
		dd_real temp;

		x86_FIX
		temp = dd_real ((double) j) * dd_data->gw__bits_per_word;
		twopowdbl = (double) (ceil (temp) - temp);
		END_x86_FIX
	}

	if (! dd_data->gw__c_is_one)
		twopowdbl += dd_data->gwdbl__log2_abs_c_div_fftlen * (double) j;
	return (pow (2.0, - double (twopowdbl)));
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
	dd_real temp, twopow, result;

	x86_FIX
	temp = dd_real ((double) j) * dd_data->gw__bits_per_word;
	twopow = ceil (temp) - temp;
	if (! dd_data->gw__c_is_one)
		twopow += dd_data->gw__log2_abs_c_div_fftlen * dd_real ((double) j);
	result = exp (dd_real::_log2 * -twopow) * dd_data->gw__over_fftlen;
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
	dd_real temp, twopow, weight;

	x86_FIX
	temp = dd_real ((double) j) * dd_data->gw__bits_per_word;
	twopow = ceil (temp) - temp;
	if (! dd_data->gw__c_is_one)
		twopow += dd_data->gw__log2_abs_c_div_fftlen * dd_real ((double) j);
	weight = exp (dd_real::_log2 * twopow);
	*fft_weight = double (weight);
	weight = dd_real (1.0) / weight;
	if (fft_weight_inverse != NULL)
		*fft_weight_inverse = double (weight);
	if (fft_weight_inverse_over_fftlen != NULL)
		*fft_weight_inverse_over_fftlen = double (weight * dd_data->gw__over_fftlen);
	END_x86_FIX
}

// Returns log2(fft_weight).  This is used in determining the FFT weight
// fudge factor in two-pass FFTs.  This is much faster than computing the
// fft_weight because it eliminates a call to the double-double exp routine.

extern "C"
double gwfft_weight_exponent (
	void	*dd_data_arg,
	unsigned long j)
{
	double	tempdbl, twopowdbl;
	dd_real temp, twopow;

// For speed, try this with plain old doubles first

	if (j == 0) return (0);
	tempdbl = (double) j * dd_data->gwdbl__bits_per_word;
	twopowdbl = ceil (tempdbl) - tempdbl;
	if (twopowdbl > 0.001 && twopowdbl < 0.999) return (twopowdbl);

// If at all uncertain of the result, use doubledoubles to do the calculation

	x86_FIX
	temp = dd_real ((double) j) * dd_data->gw__bits_per_word;
	twopow = ceil (temp) - temp;
	END_x86_FIX
	return (double (twopow));
}

//
// Utility routine to compute fft base for j-th fft word
//
// The FFT base for the j-th FFT word doing a 2^q+c weighted transform is
//	ceil (j*q/FFTLEN)
// This routine returns ceil (j*q/FFTLEN) taking great care to return a
// value accurate to 53 bits.  This is important when j*q is really close to
// a multiple of FFTLEN (admittedly quite rare).  It would be really bad if
// rounding differences caused this routine to compute ceil (j*q/FFTLEN)
// differently than the weighting functions.
//

extern "C"
unsigned long gwfft_base (
	void	*dd_data_arg,
	unsigned long j)
{
	double	tempdbl, ceildbl, diffdbl;
	dd_real temp;
	unsigned long twopow;

// For speed, try this with plain old doubles first

	if (j == 0) return (0);
	tempdbl = (double) j * dd_data->gwdbl__bits_per_word;
	ceildbl = ceil (tempdbl);
	diffdbl = ceildbl - tempdbl;
	if (diffdbl > 0.001 && diffdbl < 0.999)
		return ((unsigned long) ceildbl);

// If at all uncertain of the result, use doubledoubles to do the calculation

	x86_FIX
	temp = dd_real ((double) j) * dd_data->gw__bits_per_word;
	twopow = (int) ceil (temp);
	END_x86_FIX
	return (twopow);
}
