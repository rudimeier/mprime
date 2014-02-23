/*----------------------------------------------------------------------
| gwdbldbl.h
|
| This file contains the headers for the gwnum helper routines that use
| extended-precision floats.
| 
|  Copyright 2005-2014 Mersenne Research, Inc.  All rights reserved.
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

#define gwsincos1by1_raw(a,b,c)	gwsincos1by_raw(a,b,c,1)
#define gwsincos1by2_raw(a,b,c)	gwsincos1by_raw(a,b,c,2)
#define gwsincos1by4_raw(a,b,c)	gwsincos1by_raw(a,b,c,4)
void gwsincos1by_raw (unsigned long, unsigned long, double *, int);

#define gwsincos1by1(a,b,c)	gwsincos1by(a,b,c,1)
#define gwsincos1by2(a,b,c)	gwsincos1by(a,b,c,2)
#define gwsincos1by4(a,b,c)	gwsincos1by(a,b,c,4)
void gwsincos1by (unsigned long, unsigned long, double *, int);

#define gwsincos12by1_raw(a,b,c) gwsincos12by_raw(a,b,c,1)
#define gwsincos12by2_raw(a,b,c) gwsincos12by_raw(a,b,c,2)
#define gwsincos12by4_raw(a,b,c) gwsincos12by_raw(a,b,c,4)
void gwsincos12by_raw (unsigned long, unsigned long, double *, int);

#define gwsincos12by1(a,b,c)	gwsincos12by(a,b,c,1)
#define gwsincos12by2(a,b,c)	gwsincos12by(a,b,c,2)
#define gwsincos12by4(a,b,c)	gwsincos12by(a,b,c,4)
void gwsincos12by (unsigned long, unsigned long, double *, int);

#define gwsincos13by2(a,b,c)	gwsincos13by(a,b,c,2)
#define gwsincos13by4(a,b,c)	gwsincos13by(a,b,c,4)
void gwsincos13by (unsigned long, unsigned long, double *, int);

#define gwsincos15by2(a,b,c)	gwsincos15by(a,b,c,2)
#define gwsincos15by4(a,b,c)	gwsincos15by(a,b,c,4)
void gwsincos15by (unsigned long, unsigned long, double *, int);

#define gwsincos15913by2(a,b,c)	gwsincos15913by(a,b,c,2)
#define gwsincos15913by4(a,b,c)	gwsincos15913by(a,b,c,4)
void gwsincos15913by (unsigned long, unsigned long, double *, int);

#define gwsincos125by2(a,b,c)	gwsincos125by(a,b,c,2)
#define gwsincos125by4(a,b,c)	gwsincos125by(a,b,c,4)
void gwsincos125by (unsigned long, unsigned long, double *, int);

#define gwsincos1234by2_raw(a,b,c)	gwsincos1234by_raw(a,b,c,2)
#define gwsincos1234by4_raw(a,b,c)	gwsincos1234by_raw(a,b,c,4)
void gwsincos1234by_raw (unsigned long, unsigned long, double *, int);

#define gwsincos1234by2(a,b,c)	gwsincos1234by(a,b,c,2)
#define gwsincos1234by4(a,b,c)	gwsincos1234by(a,b,c,4)
void gwsincos1234by (unsigned long, unsigned long, double *, int);

#define gwsincos1plus0123by1(a,b,c,d)	gwsincos1plus0123by(a,b,c,d,1)
#define gwsincos1plus0123by2(a,b,c,d)	gwsincos1plus0123by(a,b,c,d,2)
#define gwsincos1plus0123by4(a,b,c,d)	gwsincos1plus0123by(a,b,c,d,4)
void gwsincos1plus0123by (unsigned long, unsigned long, unsigned long, double *, int);

#define gwsincos1plus01234567by2(a,b,c,d)	gwsincos1plus01234567by(a,b,c,d,2)
#define gwsincos1plus01234567by4(a,b,c,d)	gwsincos1plus01234567by(a,b,c,d,4)
void gwsincos1plus01234567by (unsigned long, unsigned long, unsigned long, double *, int);

#define gwsincos012by2_weighted(a,b,c,d,e,f)	gwsincos012by_weighted(a,b,c,d,e,f,2)
#define gwsincos012by4_weighted(a,b,c,d,e,f)	gwsincos012by_weighted(a,b,c,d,e,f,4)
void gwsincos012by_weighted (void *, unsigned long, unsigned long, unsigned long, unsigned long, double *, int);

#define gwsincos15by2_weighted(a,b,c,d,e,f)	gwsincos15by_weighted(a,b,c,d,e,f,2)
#define gwsincos15by4_weighted(a,b,c,d,e,f)	gwsincos15by_weighted(a,b,c,d,e,f,4)
void gwsincos15by_weighted (void *, unsigned long, unsigned long, unsigned long, unsigned long, double *, int);

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
void gwfft_weights_fudged (void	*, unsigned long, unsigned long, double *, double *, double *, double *);
double gwfft_partial_weight (void *, unsigned long, unsigned long);
double gwfft_partial_weight_sloppy (void *, unsigned long, unsigned long);
double gwfft_partial_weight_inverse (void *, unsigned long, unsigned long);
double gwfft_partial_weight_inverse_sloppy (void *, unsigned long, unsigned long);

#ifdef __cplusplus
}
#endif

#endif
