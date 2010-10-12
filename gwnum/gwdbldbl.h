/*----------------------------------------------------------------------
| gwdbldbl.h
|
| This file contains the headers for the gwnum helper routines that use
| extended-precision floats.
| 
|  Copyright 2005-2010 Mersenne Research, Inc.  All rights reserved.
+---------------------------------------------------------------------*/

#ifndef _GWDBLDBL_H
#define _GWDBLDBL_H

/* This is a C library.  If used in a C++ program, don't let the C++ */
/* compiler mangle names. */

#ifdef __cplusplus
extern "C" {
#endif

/* Include common definitions */

#include "gwcommon.h"

/* Extended precision helper routines */

void gwasm_constants (double *);
void gwsincos (unsigned long, unsigned long, double *);
void gwsincos3 (unsigned long, unsigned long, double *);
void gwsincos5 (unsigned long, unsigned long, double *);
void gwsincos1by2 (unsigned long, unsigned long, double *);
void gwsincos12by2 (unsigned long, unsigned long, double *);
void gwsincos13by2 (unsigned long, unsigned long, double *);
void gwsincos15by2 (unsigned long, unsigned long, double *);
void gwsincos125by2 (unsigned long, unsigned long, double *);
void gwsincos1plus0123by2 (unsigned long, unsigned long, unsigned long, double *);
void gwsincos1plus01234567by2 (unsigned long, unsigned long, unsigned long, double *);
void *gwdbldbl_data_alloc (void);
void gwfft_weight_setup (void *, int, double, unsigned long, unsigned long, signed long, unsigned long);
double gwfft_weight (void *, unsigned long);
double gwfft_weight_sloppy (void *, unsigned long);
double gwfft_weight_inverse (void *, unsigned long);
double gwfft_weight_inverse_sloppy (void *, unsigned long);
double gwfft_weight_inverse_over_fftlen (void *, unsigned long);
void gwfft_weights3 (void *, unsigned long, double *, double *, double *);
double gwfft_weight_exponent (void *, unsigned long);
unsigned long gwfft_base (void *, unsigned long);
void gwsincos12by2_weighted (void *, unsigned long, unsigned long, unsigned long, unsigned long, double *);
void gwsincos15by2_weighted (void *, unsigned long, unsigned long, unsigned long, unsigned long, double *);
void gwfft_weights_fudged (void	*, unsigned long, unsigned long, double *, double *, double *, double *);
double gwfft_partial_weight (void *, unsigned long, unsigned long);
double gwfft_partial_weight_sloppy (void *, unsigned long, unsigned long);
double gwfft_partial_weight_inverse (void *, unsigned long, unsigned long);
double gwfft_partial_weight_inverse_sloppy (void *, unsigned long, unsigned long);

#ifdef __cplusplus
}
#endif

#endif
