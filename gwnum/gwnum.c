/*----------------------------------------------------------------------
| gwnum.c
|
| This file contains the C routines and global variables that are used
| in the multi-precision arithmetic routines.  That is, all routines
| that deal with the gwnum data type.
| 
|  Copyright 2002-2009 Mersenne Research, Inc.  All rights reserved.
+---------------------------------------------------------------------*/

/* Include files */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <memory.h>
#include <string.h>
#include "cpuid.h"
#include "gwnum.h"
#include "gwutil.h"
#include "gwdbldbl.h"

/* Include a random number generator.  For reasons we'll discuss later */
/* we do not want to use the C runtime library's random number generator */
/* to initialize GW_RANDOM */

#include "mt19937ar.c"

/* Handy macros to improve readability */

#define log2(n)		(log((double)(n)) / log (2.0))
#define logb(n)		(log((double)(n)) / log ((double)(b)))

/* MSVC6 has trouble with the pow function using integer arguments. */
/* For example, "(unsigned long) pow (5.0, 7.0)" returns 78124 instead */
/* of the correct 78125.  This macro, works around this trouble. */

#define intpow(b,n)	((long) floor (pow ((double)(b), (double)(n)) + 0.1))

/* global variables */

/* When debugging gwnum and giants, I sometimes write code that "cheats" */
/* by calling a routine that is part of prime95 rather than the gwnum and */
/* giants library.  Prime95 will set this routine pointer so that gwnum */
/* code can cheat while keeping the gwnum library interface clean. */

void (*OutputBothRoutine)(int, char *) = NULL;
#define OutputBoth(t,x)	(*OutputBothRoutine)(t,x)

/* Assembly helper routines */

struct gwinfo1_data {
	struct gwasm_jmptab *p4_cyclic_fft_info; /* P4 mersenne mod FFT info */
	struct gwasm_jmptab *p4_negacyclic_fft_info; /* P4 2^N+1 mod FFT */
	struct gwasm_jmptab *x86_cyclic_fft_info; /* x86 mersenne mod FFT */
	struct gwasm_jmptab *x86_negacyclic_fft_info; /* x86 2^N+1 mod FFT */
	uint32_t version;			/* Gwnum lib version number */
};

void gwinfo1 (struct gwinfo1_data *);

/* This structure declares all the global constants and temporaries needed */
/* by the assembly code.  It is allocated and initialized with C code. */
/* In 32-bit mode, the area before the structure also serves as stack space. */
/* This lets the assembly code which is starved for registers use the stack */
/* pointer to access the data in this structure. */

#ifdef X86_64
#define NEW_STACK_SIZE	0
#else
#define NEW_STACK_SIZE	(8192+256)
#endif
struct gwasm_data {
	void	*DESTARG;	/* Function argument */
	void	*DIST_TO_FFTSRCARG;/* SRCARG - DESTARG */
	void	*DIST_TO_MULSRCARG;/* SRC2ARG - DESTARG */
	void	*NORMRTN;	/* Post-multiply normalization routine */

	int32_t	ZERO_PADDED_FFT;/* True if doing a zero pad FFT */
	int32_t	POSTFFT;	/* True if assembly code can start the */
				/* FFT process on the result of a multiply */
	double	MAXERR;		/* Convolution error in a multiplication */

	void	*scratch_area;	/* Scratch area for pass 1 of SSE2 FFTs */
	void	*norm_grp_mults;/* Ptr to array #1 of normalize multipliers */
	void	*norm_col_mults;/* Ptr to array #2 of normalize multipliers */
	void	*norm_biglit_array;/* Ptr to byte array of big/lit flags */

	void	*carries;	/* Ptr to array of carries (2 pass FFT) */
	void	*plus1_premults;/* Address of 2^N+1 premultiplier data */
	intptr_t pass2gapsize;	/* Wasted bytes between pass2 data blks */
	void	*sincos1;

	void	*sincos2;
	void	*sincos3;
	void	*sincos4;
	void	*sincos5;

	intptr_t pass1blkdst;	/* Dist between blocks in pass 1 */
	void	*pass2_premults;/* Address of pass 2 premultiplier data */
	void	*xsincos_complex;/* Addr of pass2 complex sin/cos data */
	void	*sincos6;

	void	*sincos7;
	void	*sincos8;
	void	*sincos9;
	void	*sincos10;

	double	XMM_TWO[2];	/* 2.0 */
	double	XMM_HALF[2];	/* 0.5 */
	double	XMM_SQRTHALF[2];/* SQRT(0.5) */
	double	XMM_TMP1_8[16];	/* 8 XMM temporaries */

	int32_t	ffttype;	/* Type of fft (1, 2, 3, or 4) */
	int32_t	count1;		/* Counter used in common fft code */
	uint32_t ADDIN_ROW;	/* For adding a constant after multiply */
	uint32_t ADDIN_OFFSET;

	int32_t	addcount1;	/* Counter used in common fft code */
	int32_t	normval1;	/* Add / sub constants */
	int32_t	normval4;
	int32_t	ALL_COMPLEX_FFT;/* True if using all-complex FFTs */

	int32_t	cache_line_multiplier; /* Cache line multiplier from jmptab */
	int32_t	RATIONAL_FFT;	/* True if bits per FFT word is integer */
	int32_t	TOP_CARRY_NEEDS_ADJUSTING;/* True when carry out of top word */
				/* needs adjusting */
	int32_t	SPREAD_CARRY_OVER_4_WORDS;/* True when carry out of top word */
				/* must be spread over more than 2 words */

	int32_t	zero_fft;	/* TRUE if zero upper half in normalize */
	int32_t	const_fft;	/* TRUE if mul-by-const in normalize */
	double	ADDIN_VALUE;	/* Value to add in after a multiply */

	void	*norm_ptr1;
	void	*norm_ptr2;
	intptr_t normblkdst;	/* Dist between blocks in normalization */
	intptr_t normblkdst8;	/* Dist between 8 blocks in normalize */

	double	XMM_SUMOUT[2];	/* Used in normalization macros */
	double	XMM_MAXERR[2];	/* Used in normalization macros */
	double	XMM_LIMIT_BIGMAX[32]; /* Normalization constants */
	double	XMM_LIMIT_INVERSE[32];
	double	XMM_COL_MULTS[1024];
	double	XMM_TTP_FUDGE[32];
	double	XMM_TTMP_FUDGE[32];
	double	XMM_BIGVAL[2];	/* Used to round double to integer */
	double	XMM_MINUS_C[2];	/* -c stored as double */
	double	XMM_NORM012_FF[2]; /* Used in xnorm012 macros (FFTLEN/2) */
	double	XMM_LIMIT_BIGMAX_NEG[32];
	int32_t	XMM_ABSVAL[4];	/* Used to compute absolute values */
	double	XMM_P618[2];
	double	XMM_P309[2];
	double	XMM_M809[2];
	double	XMM_M262[2];
	double	XMM_M382[2];
	double	XMM_P951[2];
	double	XMM_P588[2];
	double	XMM_M162[2];
	double	DBLARG;
	double	CARRY_ADJUST7;
	double	XMM_P866[2];
	double	XMM_BIGVAL_NEG[2];
	double	CARRY_ADJUST1_HI;
	double	CARRY_ADJUST1_LO;
	double	XMM_UNUSED4[2];
	double	XMM_UNUSED5[2];
	double	XMM_P623[2];
	double	XMM_M358[2];
	double	XMM_P404[2];
	double	XMM_P975[2];
	double	XMM_P445[2];
	double	XMM_P180[2];

	uint32_t COPYZERO[8];	/* Offsets to help in gwcopyzero */

	uint32_t B_IS_2;
	uint32_t NUMARG;	/* Gwcopyzero assembly arg */
	double	ZPAD_INVERSE_K6;

	uint32_t *ASM_TIMERS;	/* Timers used for optimizing code */
	void	*sincos11;
	void	*sincos12;
	void	*zpad_addr;	/* Address of XMM_ZPAD0 */

	uint32_t FFTLEN;	/* The FFT size we are using */
	uint32_t ZPAD_TYPE;	/* 1,2,or 3 words in k (used by zero pad) */
	uint32_t BIGLIT_INCR2;	/* Offset to step in big/lit array */
	uint32_t BIGLIT_INCR4;	/* Offset to step in big/lit array */

	double	ZPAD_K6_HI;	 /* Zero padded FFT constants */
	double	ZPAD_K6_MID;
	double	ZPAD_K6_LO;
	double	ZPAD_SHIFT6;
	double	ZPAD_INVERSE_K5;
	double	ZPAD_K5_HI;
	double	ZPAD_K5_MID;
	double	ZPAD_K5_LO;
	double	ZPAD_SHIFT5;
	double	ZPAD_INVERSE_K4;
	double	ZPAD_K4_HI;
	double	ZPAD_K4_MID;
	double	ZPAD_K4_LO;
	double	ZPAD_SHIFT4;
	double	ZPAD_INVERSE_K3;
	double	ZPAD_K3_HI;
	double	ZPAD_K3_MID;
	double	ZPAD_K3_LO;
	double	ZPAD_SHIFT3;
	double	ZPAD_INVERSE_K2;
	double	ZPAD_K2_HI;
	double	ZPAD_K2_MID;
	double	ZPAD_K2_LO;
	double	ZPAD_SHIFT2;
	double	ZPAD_K1_HI;
	double	ZPAD_INVERSE_K1;
	double	ZPAD_K1_LO;
	double	ZPAD_SHIFT1;
	double	ZPAD0_7[8];

	double	XMM_BIGBIGVAL[2]; /* Used to round to big word multiple */
	double	XMM_MULCONST[2];
	double	XMM_K_TIMES_MULCONST_LO[2]; /* k*mulconst, low big word bits */
	double	XMM_K_TIMES_MULCONST_HI[2]; /* k*mulconst, top bits */
	double	XMM_MINUS_C_TIMES_MULCONST[2];
	double	XMM_K_TIMES_MULCONST_HI_1[2]; /* XMM_K_TIMES_MULCONST_HI, low bits */
	double	XMM_K_TIMES_MULCONST_HI_2[2]; /* XMM_K_TIMES_MULCONST_HI, top bits */
	double	XMM_K_LO[2];	/* k, bottom big word bits */
	double	XMM_K_HI[2];	/* k, top bits */
	double	XMM_K_HI_1[2];	/* XMM_K_HI, bottom big word bits */
	double	XMM_K_HI_2[2];	/* XMM_K_HI, top bits */

	double	INVERSE_K;	/* 1/K */
	double	K_HI;		/* Upper bits of K */
	double	K_LO;		/* Lower bits of K */
	double	CARRY_ADJUST1;	/* Adjustment constant #1 in wrapping carry */
	double	CARRY_ADJUST2;	/* Adjustment constant #2 in wrapping carry */
	double	CARRY_ADJUST3;	/* Adjustment constant #3 in wrapping carry */
	double	CARRY_ADJUST4;	/* Adjustment constant #4 in wrapping carry */
	double	CARRY_ADJUST5;	/* Adjustment constant #5 in wrapping carry */
	double	CARRY_ADJUST6;	/* Adjustment constant #6 in wrapping carry */

	uint32_t HIGH_WORD1_OFFSET;/* Offset of top FFT word */
	uint32_t HIGH_WORD2_OFFSET;/* Offset of second high FFT word */
	uint32_t HIGH_WORD3_OFFSET;/* Offset of third high FFT word */
	uint32_t HIGH_SCRATCH1_OFFSET;
				/* Offset of top FFT word in scratch area */
	uint32_t HIGH_SCRATCH2_OFFSET;
				/* Offset of second highest FFT word */
	uint32_t HIGH_SCRATCH3_OFFSET;
				/* Offset of third highest FFT word */

	void	*SAVED_RSP;	/* Saved stack pointer in 32-bit mode */
	void	*SRCARG;	/* Function argument */
	void	*SRC2ARG;	/* Function argument */
	void	*DEST2ARG;	/* Function argument */

	intptr_t normval2;	/* Add / sub temporaries */
	intptr_t normval3;
	intptr_t ZPAD_WORD5_OFFSET; /* Offset of the 5th word */
	intptr_t ZPAD_WORD5_RBP_OFFSET; /* Offset to 5th word's ttp data */

	int32_t	normcount1;
	int32_t	count2;
	int32_t	count3;
	int32_t	count4;

	int32_t	count5;
	float	HALF;		/* 0.5 */
	float	BIGVAL;		/* Rounding value */
	float	BIGBIGVAL;

	double	SQRTHALF;	/* SQRT(0.5) */
	double	MINUS_C;	/* -c */
	double	ttmp_ff_inv;	/* Inverse FFT adjust (2/FFTLEN) */
	double	P309;
	double	M809;
	double	M262;
	double	M382;
	double	P951;
	double	P588;
	double	M162;
	double	P618;
	double	P623;
	double	M358;
	double	P404;
	double	P975;
	double	P445;
	double	P180;
	double	M223;
	double	M901;
	double	M691;
	double	P866;
	double	P433;
	double	P577;
	float	P25;		/* 0.25 */
	float	P75;		/* 0.75 */
	float	P3;		/* 3.0 */

	int32_t	thread_num;	/* Thread num - so we can differentiate main */
				/* thread from auxillary threads */
	void	*data_addr;	/* FFT data to work on this pass */
	void	*data_prefetch;	/* FFT data to prefetch for next pass */
	void	*premult_addr;	/* Premult data to use this pass */
	void	*premult_prefetch;/* Premult data to prefetch for next pass */
	gwhandle *gwdata;	/* Allows callback routines to access gwdata */
	void	*pass1_wake_up_threads; /* Callback routine to wake up */
				/* auxillary threads */
	void	*pass1_wait_for_carries; /* Callback routine to ensure the */
				/* normalization routine accesses the carry */
				/* data in ascending block numbers */
	void	*pass1_done_with_carries; /* Callback routine to let the */
				/* next block access the carry data */
	void	*pass1_get_next_block; /* Callback routine to get next block */
				/* number for thread to process */
	void	*pass2_wake_up_threads; /* Callback routine to wake up */
				/* auxillary threads */
	void	*pass2_get_next_block; /* Callback routine to get next block */
				/* number for thread to process */
	void	(*thread_work_routine)(void*); /* Assembly routine to call */
				/* when auxillary thread wakes up */
	uint32_t this_block; /* Block currently being processed */
	uint32_t next_block; /* Next block to process */
	uint32_t last_pass1_block; /* Last block to process */
};

/* gwnum assembly routine pointers */

#define gw_fft(h,a)	(*(h)->GWPROCPTRS[0])(a)
#define gw_square(h,a)	(*(h)->GWPROCPTRS[1])(a)
#define gw_mul(h,a)	(*(h)->GWPROCPTRS[2])(a)
#define gw_mulf(h,a)	(*(h)->GWPROCPTRS[3])(a)
#define gw_add(h,a)	(*(h)->GWPROCPTRS[4])(a)
#define gw_addq(h,a)	(*(h)->GWPROCPTRS[5])(a)
#define gw_sub(h,a)	(*(h)->GWPROCPTRS[6])(a)
#define gw_subq(h,a)	(*(h)->GWPROCPTRS[7])(a)
#define gw_addsub(h,a)	(*(h)->GWPROCPTRS[8])(a)
#define gw_addsubq(h,a)	(*(h)->GWPROCPTRS[9])(a)
#define gw_copyzero(h,a) (*(h)->GWPROCPTRS[10])(a)
#define gw_addf(h,a)	(*(h)->GWPROCPTRS[11])(a)
#define gw_subf(h,a)	(*(h)->GWPROCPTRS[12])(a)
#define gw_addsubf(h,a)	(*(h)->GWPROCPTRS[13])(a)
#define gw_adds(h,a)	(*(h)->GWPROCPTRS[14])(a)
#define gw_muls(h,a)	(*(h)->GWPROCPTRS[15])(a)
#define norm_routines	16
#define copyzero_routines 20

/* Declare the hundreds of assembly normalization routines.  Naming of these */
/* routines is as follows: */
/* x or blank		(SSE2 or x87) */
/* r or i		(rational or irrational) */
/* 1 or 2		(1 or 2 pass FFTs) */
/* zp or blank		(zero padded FFT or not) */
/* e or blank		(roundoff error checking or not) */
/* c or blank		(mul by a small constant or not) */
/* b or blank		(b > 2 or not) */
/* s4 or blank		(SSE4 or not) */
/* k or blank		(XMM_K_HI is zero or not (zp only)) */
/* c1 or cm1 or blank	(c=1, c=-1, or abs(c)!=1 (zp only, not c only)) */
/* AMD or blank		(AMD prefetching or not (x only, 2 only)) */
#define extrn_explode()		extrn_explode1(r)	extrn_explode1(i)
#define extrn_explode1(name)	extrn_explode2(name)	extrn_explode2(x##name)
#define extrn_explode2(name)	extrn_explode3(name##1)	extrn_explode3(name##2)
#define extrn_explode3(name)	extrn_explode4(name)	extrn_explode4(name##zp) 
#define extrn_explode4(name)	extrn_explode5(name)	extrn_explode5(name##e)
#define extrn_explode5(name)	extrn_explode6(name)	extrn_explode6(name##c)
#define extrn_explode6(name)	extrn_explode7(name)	extrn_explode7(name##b)
#define extrn_explode7(name)	extrn_explode8(name)	extrn_explode8(name##s4)
#define extrn_explode8(name)	extrn_explode9(name)	extrn_explode9(name##k)
#define extrn_explode9(name)	extrn_explode10(name)	extrn_explode10(name##c1)	extrn_explode10(name##cm1)
#define extrn_explode10(name)	extrn_explode11(name)	extrn_explode11(name##AMD)
#define extrn_explode11(name)	extern void (*name)(void);
extrn_explode()

/* Now we put the normalization routines in an array so we can easily pick */
/* the normalization routines to use.  There is one table for SSE2 and one */
/* for x87.  The SSE2 table has these 480 combinations: */
/*	r or i */
/*	1 or 2 or 2AMD */
/*	e or blank */
/*	b or blank */
/*	s4 or blank */
/*	zpck or zpk or zpkc1 or zpkcm1 or zpc or zp or zpc1 or zpcm1 or c or blank */
/* We also define a macro that will pick the correct entry from the array. */

#define sse2_explode()			sse2_explode1(xr)		sse2_explode1(xi)		NULL
#define sse2_explode1(name)		sse2_explode2(name##1,notAMD)	sse2_explode2(name##2,notAMD)	sse2_explode2(name##2,AMD)
#define sse2_explode2(name,suff)	sse2_explode3(name,suff,notzp)	sse2_explode3(name##zp,suff,zp)
#define sse2_explode3(name,suff,zp)	sse2_explode4(name,suff,zp)	sse2_explode4(name##e,suff,zp)
#define sse2_explode4(name,suff,zp)	sse2_explode5(name,suff,zp,notc) sse2_explode5(name##c,suff,zp,c)
#define sse2_explode5(name,suff,zp,c)	sse2_explode6(name,suff,zp,c)	sse2_explode6(name##b,suff,zp,c)
#define sse2_explode6(name,suff,zp,c)	sse2_explode7##zp(name,suff,c)	sse2_explode7##zp(name##s4,suff,c)
#define sse2_explode7notzp(name,suff,c)	sse2_explode9##suff(name)
#define sse2_explode7zp(name,suff,c)	sse2_explode8##c(name,suff)	sse2_explode8##c(name##k,suff)
#define sse2_explode8c(name,suff)	sse2_explode9##suff(name)
#define sse2_explode8notc(name,suff)	sse2_explode9##suff(name)	sse2_explode9##suff(name##c1)	sse2_explode9##suff(name##cm1)
#define sse2_explode9notAMD(name)	&name,
#define sse2_explode9AMD(name)		&name##AMD,

void *sse2_prctab[] = { sse2_explode() };

int sse2_prctab_index (gwhandle *gwdata, int e, int c)
{
	int	index = 0;

	if (! gwdata->RATIONAL_FFT) index += 240;
	if (gwdata->PASS2_LEVELS) {
		if (gwdata->cpu_flags & CPU_3DNOW) index += 160;
		else index += 80;
	}
	if (! gwdata->ZERO_PADDED_FFT) {
		if (e) index += 8;
		if (c) index += 4;
		if (gwdata->b > 2) index += 2;
		if (gwdata->cpu_flags & CPU_SSE41) index += 1;
	} else {
		struct gwasm_data *asm_data = (struct gwasm_data *) gwdata->asm_data;
		index += 16;
		if (e) index += 32;
		if (!c) {
			if (gwdata->b > 2) index += 12;
			if (gwdata->cpu_flags & CPU_SSE41) index += 6;
			if (asm_data->XMM_K_HI[0] == 0.0) index += 3;
			if (gwdata->c == 1) index += 1;  if (gwdata->c == -1) index += 2;
		} else {
			index += 24;
			if (gwdata->b > 2) index += 4;
			if (gwdata->cpu_flags & CPU_SSE41) index += 2;
			if (asm_data->XMM_K_TIMES_MULCONST_HI[0] == 0.0) index += 1;
		}
	}
	 return (index);
}				    

/* The x87 normalization routines array has 32 combinations: */
/*	r or i */
/*	1 or 2 */
/*	zp or blank */
/*	e or blank */
/*	c or blank */
/* We also define a macro that will pick the correct entry from the array. */

#ifndef X86_64
#define x87_explode()			x87_explode1(r)		x87_explode1(i) NULL
#define x87_explode1(name)		x87_explode2(name##1)	x87_explode2(name##2)
#define x87_explode2(name)		x87_explode3(name)	x87_explode3(name##zp)
#define x87_explode3(name)		x87_explode4(name)	x87_explode4(name##e)
#define x87_explode4(name)		x87_explode5(name)	x87_explode5(name##c)
#define x87_explode5(name)		&name,

void *x87_prctab[] = { x87_explode() };

#define	x87_prctab_index(gwdata, e, c)  \
	    ((gwdata)->RATIONAL_FFT ? 0 : 16) + \
	    ((gwdata)->PASS2_LEVELS ? 8 : 0) + \
	    ((gwdata)->ZERO_PADDED_FFT ? 4 : 0) + \
	    (e ? 2 : 0) + \
	    (c ? 1 : 0)
#endif

/* Helper macros */

#ifndef isinf
#define isinf(x)		((x) != 0.0 && (x) == 2.0*(x))
#endif
#ifndef isnan
#define isnan(x)		((x) != (x))
#endif
#define is_valid_double(x)	(! isnan (x) && ! isinf (x))

/* Error codes returned */

#define GWERROR_BAD_FFT_DATA	-1

/* Special flag value set in the gwnum freeable field */

#define GWFREED_TEMPORARILY	0x80000000

/* Forward declarations */

int convert_giant_to_k2ncd (
	giant	g,		/* Giant to examine */
	double	*k,		/* K in (K*2^N+C)/D. */
	unsigned long *n,	/* N in (K*2^N+C)/D. */
	signed long *c,		/* C in (K*2^N+C)/D. */
	unsigned long *d);	/* D in (K*2^N+C)/D. */
int internal_gwsetup (
	gwhandle *gwdata,	/* Placeholder for gwnum global data */
	double	k,		/* K in K*B^N+C. Must be a positive integer. */
	unsigned long b,	/* B in K*B^N+C. */
	unsigned long n,	/* N in K*B^N+C. Exponent to test. */
	signed long c);		/* C in K*B^N+C. Must be rel. prime to K. */
long nonbase2_gianttogw (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	giant	a,
	gwnum	g,
	unsigned long limit,	/* How many FFT words to set */
	unsigned long offset,	/* Offset into FFT array of words to set */
	long	carry);		/* Carry to add into this section */
void raw_gwsetaddin (gwhandle *gwdata, unsigned long word, double val);
int multithread_init (gwhandle *gwdata);
void multithread_term (gwhandle *gwdata);
void pass1_aux_entry_point (void*);
void pass2_aux_entry_point (void*);

/* Find the power of two greater than or equal to N. */

unsigned long pow_two_above_or_equal (
	unsigned long n)
{
	unsigned long result;

	result = 1;
	for (n = n - 1; n; n = n >> 1) result = result << 1;
	return (result);
}


/* This routine builds a sin/cos table - used by gwsetup */

double *build_sin_cos_table (
	double	*table,		/* Pointer to the table to fill in */
	unsigned long N,	/* Number of DATA values processed by this */
				/* FFT level.  This explains the divide by 2 */
				/* for complex FFTs later in this routine */
	int	hermetian_skip,	/* True if some sin/cos values are skipped */
	int	type)		/* 0 = old style - a plain old array */
				/* 1 = SSE2 - data is duplicated */
				/* 2 = SSE2 - data is interleaved */
{
	unsigned long i;

/* Handle hermetian skip when interleaving.  First data slot is left */
/* undefined. */

	if (type == 2 && hermetian_skip) type = 3;

/* Special case the really small sin/cos tables.  If N is between 9 and 16 */
/* or between 33 and 64, then the assembly code is only doing one FFT level. */
/* In this case, the code just uses the middle sin/cos values of a 2N sized */
/* table.  We could optimize this inefficient memory usage at a later date. */

	if (N <= 8) return (table);
	if (N >= 9 && N <= 16) N = N * 2;
	if (N >= 33 && N <= 64 && type == 1 && hermetian_skip) N = N * 2;

/* In the all-complex case. build the same size table as the hermetian */
/* case which skips half the i values. */

	if (!hermetian_skip) N = N / 2;

/* Loop to build table. */

	for (i = hermetian_skip ? ((N & 4) ? 4 : 8) : 0; i < N; i += 4) {
		unsigned long shifted_i, shifted_N, flipped_i;
		double	sincos[6];

/* Flip the bits in i.  Our prime-factor-FFT makes this a little complex. */
/* The algorithm below works, but I've long since forgotten why. */

		shifted_i = i; shifted_N = N; flipped_i = 0;
		while ((shifted_N & 1) == 0) {
			flipped_i <<= 1;
			if (shifted_i & 1) flipped_i++;
			shifted_i >>= 1;
			shifted_N >>= 1;
		}
		flipped_i = (flipped_i * shifted_N) + shifted_i;

/* When the FFT is working on real data Hermetian symettry allows us to */
/* eliminate half of the FFT data and consequently half of the sin/cos data */
/* Case 1:  If shifted source is > shifted N/2, then we */
/* do not need these sin/cos values. */
/* Case 2:  If shifted source is zero, loop to find the top */
/* two bits.  Skip the number if the top two bits equal 3. */

		if (hermetian_skip) {
			if (shifted_i > shifted_N / 2) continue;
			if (shifted_i == 0) {
				unsigned long j;
				for (j = i; j > 3; j >>= 1);
				if (j == 3) continue;
			}
		}

/* Compute the 3 sin/cos values */

		gwsincos3 (flipped_i, N, (double *) &sincos);

/* Copy the sin/cos values in the appropriate way */

		if (type == 0) {
			memcpy (table, sincos, sizeof (sincos));
			table += 6;
		} else if (type == 1) {
			table[0] = table[1] = sincos[0];
			table[2] = table[3] = sincos[1];
			table[4] = table[5] = sincos[2];
			table[6] = table[7] = sincos[3];
			table[8] = table[9] = sincos[4];
			table[10] = table[11] = sincos[5];
			table += 12;
		} else if (type == 2) {
			table[0] = sincos[0];
			table[2] = sincos[1];
			table[4] = sincos[2];
			table[6] = sincos[3];
			table[8] = sincos[4];
			table[10] = sincos[5];
			type++;
		} else {
			table[1] = sincos[0];
			table[3] = sincos[1];
			table[5] = sincos[2];
			table[7] = sincos[3];
			table[9] = sincos[4];
			table[11] = sincos[5];
			type--;
			table += 12;
		}
	}
	return (table);
}

/* This routine builds a pass 2 premultiplier table - used by gwsetup */

double *build_premult_table (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	double	*table,		/* Pointer to the table to fill in */
	unsigned long pass2_size)
{
	unsigned long i, N, incr, type, m_fudge;

/* Build a premultiplier table for the second pass incrementing by */
/* the pre-calculated pass2_size. */

	N = gwdata->FFTLEN;
	incr = pass2_size;
	if (gwdata->ALL_COMPLEX_FFT) N = N / 2;

/* Mod 2^N+1 arithmetic starts at first data set, */
/* mod 2^N-1 skips some data sets */

 	if (gwdata->ALL_COMPLEX_FFT) i = 0;
	else i = incr * 4;

/* To add in the flipped_m component, we want the sin/cos of flipped_m */
/* over pass2_size.  This fudge factor will convert flipped_m into something */
/* that can be divided by N. */

	m_fudge = N / pass2_size;

/* Loop to build table. */

	type = 0;
	for ( ; i < N; i += incr) {
		unsigned long shifted_i, shifted_N, flipped_i, k, l, m;
		unsigned long grouping_size;
		double	*table_start;
		double	sincos[2];

/* Flip the bits in i.  Our prime-factor-FFT makes this a little complex. */
/* The algorithm below works, but I've long since forgotten why. */

		shifted_i = i; shifted_N = N; flipped_i = 0;
		while ((shifted_N & 1) == 0) {
			flipped_i <<= 1;
			if (shifted_i & 1) flipped_i++;
			shifted_i >>= 1;
			shifted_N >>= 1;
		}
		flipped_i = (flipped_i * shifted_N) + shifted_i;

/* When the FFT is working on real data Hermetian symettry allows us to */
/* eliminate half of the FFT data and consequently half of the sin/cos data */
/* Case 1:  If shifted source is > shifted N/2, then we */
/* do not need these sin/cos values. */
/* Case 2:  If shifted source is zero, loop to find the top */
/* two bits.  Skip the number if the top two bits equal 3. */

		if (!gwdata->ALL_COMPLEX_FFT) {
			if (shifted_i > shifted_N / 2) continue;
			if (shifted_i == 0) {
				unsigned long j;
				for (j = i; j > 3; j >>= 1);	
				if (j == 3) continue;
			}
		}

/* Generate the group multipliers.  We used to always create groups of 4, */
/* but to save memory we now group by different amounts based on pass 2 size */

		grouping_size = 4;
		if (pass2_size == 1024) grouping_size = 8;
		if (pass2_size == 2048) grouping_size = 8;
		if (pass2_size == 4096) grouping_size = 16;
		if (pass2_size == 8192) grouping_size = 16;
		table_start = table;
		for (k = 0; k < incr / 4; k += grouping_size) {

/* There are 4 multipliers in a XMM_PMD set */

			for (l = 0; l < 4; l++) {
				unsigned long real_k, pm;

/* Compute the sin/cos value (root of unity) */

				real_k = l * incr/4 + k;
				pm = real_k * flipped_i;
				if (!gwdata->ALL_COMPLEX_FFT) {
					gwsincos (pm % N, N, (double *) &sincos);
				}

/* If C > 0, then also multiply by the proper root of -1.  This is done */
/* by changing the value we are taking the sin/cos of */

				else {
					pm = pm * 4 + real_k;
					gwsincos (pm % (N*4), N*4, (double *) &sincos);
				}

/* Save the premultiplier value */

				table[l*4+type] = sincos[0];
				table[l*4+2+type] = sincos[1];
			}
			table += 16;
		}
	
/* Generate the 16 column multipliers * first 4 sin/cos values. */
/* Also multiply by the LAST 4 sin/cos values so that the xsincos_complex */
/* table can be 1/4 of it's usual size.  The extra room in the cache more */
/* than compensates for the 12 extra column multipliers. */

		for (m = 0; m < 4; m++) {
		    unsigned long flipped_m;
		    flipped_m = ((m & 1) << 1) + ((m & 2) >> 1);
		    for (k = 0; k < grouping_size; k++) {
			for (l = 0; l < 4; l++) {
				unsigned long pm;

/* Compute the sin/cos value (root of unity) */

				pm = k * flipped_i +
				     l * flipped_m * N/16 +
				     (k & 3) * flipped_m * m_fudge;
				if (!gwdata->ALL_COMPLEX_FFT) {
					gwsincos (pm % N, N, (double *) &sincos);
				}

/* If C > 0, then also multiply by the proper root of -1.  This is done */
/* by changing the value we are taking the sin/cos of */

				else {
					pm = pm * 4 + k;
					gwsincos (pm % (N*4), N*4, (double *) &sincos);
				}

/* Save the premultiplier value */

				table[l*4+type] = sincos[0];
				table[l*4+2+type] = sincos[1];
			}
			table += 16;
		    }
		}

		if (type == 0) table = table_start;
		type = 1 - type;
 	}

	return (table);
}

/* This routine builds a plus 1 premultiplier table - used by gwsetup */
/* when c is positive. */

double *build_plus1_table (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	double	*table,		/* Pointer to the table to fill in */
	unsigned long pass2_size)
{
	unsigned long i, j, k, l, N;
	int	pfa;

/* Set flag if this is a 3*2^n FFT */

	pfa = (gwdata->FFTLEN != pow_two_above_or_equal (gwdata->FFTLEN));

/* Adjust for two-pass FFTs */

	if (pass2_size == 1) N = gwdata->FFTLEN;
	else N = gwdata->FFTLEN / (pass2_size / 2);

/* Loop to build premultiplier table in the same order as the underlying */
/* assembly macro needs them.  The pfa macro operates on 3 cache lines */
/* while the power-of-two macro operates on 2 cache lines. */
/* A 64 length FFT needs 0,8,16,24 for the macro then 3 more iterations */
/* for the cache lines beginning with 2,4,6. */
/* A 48 length FFT needs 0,8,16 and 4,12,20 for the first macro then */
/* one more iteration for the cache lines beginning with 2. */

	for (i = 0; i < N / (pfa ? 24 : 32); i++) {
	for (l = 0; l < 2; l++) {
		double	sincos[2];

/* Generate the pre multipliers (roots of -1). */

		for (k = 0; k < (unsigned long) (pfa ? 3 : 4); k++) {
		for (j = 0; j < 2; j++) {
			long	temp;

/* Compute the sin/cos value */

			if (pfa)
				temp = (long) ((i * 2 + l * N/12 + j + k * N/6) % N);
			else
				temp = (long) ((i * 4 + l * 2 + j + k * N/8) % N);
			gwsincos (temp, N*2, (double *) &sincos);

/* Save the premultiplier value */

			table[0+j] = sincos[0];
			table[2+j] = sincos[1];

/* For two-pass FFTs we could apply the root of -1 for the upper SSE2 */
/* double here or in the pass 2 premultipliers.  We've arbitrarily chosen */
/* to do it in the pass 2 premults. */

			if (pass2_size > 1) {
				j = 1;
				table[0+j] = sincos[0];
				table[2+j] = sincos[1];
			}
		}
		table += 4;
		}
	}
 	}

	return (table);
}

/* This routine builds a normalization table - used by SSE2 normalizaion */
/* routines */

double *build_norm_table (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	double	*table,		/* Pointer to the table to fill in */
	unsigned long pass2_size, /* Size of pass2 */
	int	col)		/* TRUE if building column, not group, table */
{
	unsigned long i, k, num_cols;

/* Handle one-pass FFTs first, there are no group multipliers */

	if (pass2_size == 1) {
		if (!col) return (table);

/* Loop to build table */

		for (i = 0; i < gwdata->FFTLEN; i++) {
			unsigned long j, table_entry;
			double	ttp, ttmp;

/* Call asm routines to compute the two multipliers */

			gwfft_weights3 (gwdata->dd_data, i, &ttp, NULL, &ttmp);

/* Find where this data appears in the FFT array and in the table we are building. */

			j = addr_offset (gwdata, i) / sizeof (double);
			table_entry = j >> 1;

/* Now set the entry for the MSW or LSW in an SSE2 pair */

			table[table_entry*4+(j&1)] = ttmp;
			table[table_entry*4+2+(j&1)] = ttp;
		}
		return (table + gwdata->FFTLEN + gwdata->FFTLEN);
	}

/* Two pass FFTs are handled here */

	num_cols = pass2_size / 2;
	if (col) {

/* Loop to build table */

		for (i = 0; i < num_cols; i++) {
			double	ttp, ttmp;

/* Call asm routines to compute the two multipliers */

			gwfft_weights3 (gwdata->dd_data, i, &ttp, NULL, &ttmp);

/* Now set the entry for BOTH the MSW and LSW in an SSE2 pair */

			table[i*4] = ttmp;
			table[i*4+1] = ttmp;
			table[i*4+2] = ttp;
			table[i*4+3] = ttp;
		}
		return (table + num_cols * 4);
	}

/* Build the group multipliers table */

	else {
		unsigned long pfa, h, hlimit, haddin, m, mmult, u, umult;

/* Determine if this is a PFA 5, 6, 7, or 8 */

		for (pfa = gwdata->FFTLEN; pfa > 8; pfa >>= 1);

/* Loop to build table */

		umult = gwdata->FFTLEN / 2;
		hlimit = gwdata->FFTLEN / 4 / (2*num_cols);
		for (h = 0; h < hlimit; h++) {
			if (pfa == 5) {
				if (h < hlimit / 5) {
					haddin = h * 2 * num_cols;
					mmult = gwdata->FFTLEN / 20;
				} else {
					haddin = gwdata->FFTLEN/10 + (h - hlimit/5) * 2 * num_cols;
					mmult = gwdata->FFTLEN / 5;
				}
			} else if (pfa == 7) {
				if (h < hlimit / 7) {
					haddin = h * 2 * num_cols;
					mmult = gwdata->FFTLEN / 28;
				} else if (h < 3 * hlimit / 7) {
					haddin = gwdata->FFTLEN/14 + (h - hlimit/7) * 2 * num_cols;
					mmult = gwdata->FFTLEN / 14;
				} else {
					haddin = 3*gwdata->FFTLEN/14 + (h - 3*hlimit/7) * 2 * num_cols;
					mmult = gwdata->FFTLEN / 7;
				}
			} else {
				haddin = h * 2 * num_cols;
				mmult = gwdata->FFTLEN / 4;
			}
			for (u = 0; u < 2; u++) {
			for (m = 0; m < 2; m++) {
			for (k = 0; k < 2; k++) {
				double	ttp, ttmp;
				long	n;

/* Call asm routines to compute the two multipliers */

				n = haddin + u * umult + m * mmult + k * num_cols;
				gwfft_weights3 (gwdata->dd_data, n, &ttp, &ttmp, NULL);

/* Now set the entry for BOTH the MSW and LSW in an SSE2 pair */

				table[k] = ttmp;
				table[2+k] = ttp;
			}
			table += 4;
			}
			}
		}
		return (table);
	}
}

/* This routine builds a big/little flags table - used by SSE2 normalizaion */
/* routines */

double *build_biglit_table (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	double	*table,		/* Pointer to the table to fill in */
	unsigned long pass2_size)
{
	unsigned char *p;
	unsigned long h, i, j, k, m, u, gap;
	unsigned long pfa, hlimit, haddin, mmult, umult;

/* Handle one pass FFTs differently */

	if (pass2_size == 1) {

/* Loop to build table */

		p = (unsigned char *) table;
		for (i = 0; i < gwdata->FFTLEN; i++) {
			unsigned long table_entry;

/* Find where this data appears in the FFT array and in the table we are building. */

			j = addr_offset (gwdata, i) / sizeof (double);
			table_entry = j >> 1;

/* Now set the biglit table entry for a LSW in an SSE2 pair */

			if ((j & 1) == 0) {
				p[table_entry] = is_big_word (gwdata, i) * 16;
			}

/* Otherwise, set the biglit table entry for a MSW in an SSE2 pair */

			else {
				if (is_big_word (gwdata, i)) p[table_entry] += 32;
			}
		}
		return ((double *) (p + gwdata->FFTLEN / 2));
	}

/* Determine if this is a PFA 5, 6, 7, or 8 */

	for (pfa = gwdata->FFTLEN; pfa > 8; pfa >>= 1);

/* Determine the gap between XMM high and low words */

	gap = pass2_size / 2;

/* Loop to build table in exactly the same order that it will be */
/* used by the assembly code.  This is especially ugly in the PFA cases */

	p = (unsigned char *) table;
	umult = gwdata->FFTLEN / 2;
	hlimit = gwdata->FFTLEN / 4 / (2*gap);
	for (i = 0; i < gap; i += gwdata->PASS1_CACHE_LINES) {
	for (h = 0; h < hlimit; h++) {
		if (pfa == 5) {
			if (h < hlimit / 5) {
				haddin = h * 2 * gap;
				mmult = gwdata->FFTLEN / 20;
			} else {
				haddin = gwdata->FFTLEN/10 + (h - hlimit/5) * 2 * gap;
				mmult = gwdata->FFTLEN / 5;
			}
		} else if (pfa == 7) {
			if (h < hlimit / 7) {
				haddin = h * 2 * gap;
				mmult = gwdata->FFTLEN / 28;
			} else if (h < 3 * hlimit / 7) {
				haddin = gwdata->FFTLEN/14 + (h - hlimit/7) * 2 * gap;
				mmult = gwdata->FFTLEN / 14;
			} else {
				haddin = 3*gwdata->FFTLEN/14 + (h - 3*hlimit/7) * 2 * gap;
				mmult = gwdata->FFTLEN / 7;
			}
		} else {
			haddin = h * 2 * gap;
			mmult = gwdata->FFTLEN / 4;
		}
	for (j = 0; j < gwdata->PASS1_CACHE_LINES; j++) {
	for (u = 0; u < 2; u++) {
	for (m = 0; m < 2; m++) {
	for (k = 0; k < 2 * gap; k += gap) {
		unsigned long word;

/* Now set the big/little flag for a LSW in an SSE2 pair */
/* Otherwise, set the big/little flag for a MSW in an SSE2 pair */

		word = haddin + i + j + u * umult + m * mmult + k;
		if (k == 0) *p = is_big_word (gwdata, word) * 16;
		else if (is_big_word (gwdata, word)) *p += 32;

/* Set the ttp and ttmp fudge flags for two pass FFTs.  The fudge flag is */
/* set if the col mult * the grp mult is twice the correct fft_weight, */
/* meaning a mul by 0.5 is required to generate the correct multiplier. */
/* Since we can't do equality compares on floats, this test is a little bit */
/* cryptic. */

		if (gwfft_weight_exponent (gwdata->dd_data, word) + 0.5 <
		    gwfft_weight_exponent (gwdata->dd_data, word&(gap-1)) +
		    gwfft_weight_exponent (gwdata->dd_data, word&~(gap-1))) {
			if (k == 0) *p += 64;
			else *p += 128;
		}

/* Set some offsets that help the assembly code step through the big/lit */
/* array in a non-traditional order.  Two pass-FFTs step through the array */
/* in chunks of PASS1_CACHE_LINES, but the add, sub, and carry propagation */
/* code need to access the big/lit array linearly.  Set two global variables */
/* that tell the assembly code the big/lit array distance between words */
/* 0 and 2, and words 0 and 4. */

		if (word == 2)
			((struct gwasm_data *) gwdata->asm_data)->BIGLIT_INCR2 =
				(uint32_t) ((char *) p - (char *) table);
		if (word == 4)
			((struct gwasm_data *) gwdata->asm_data)->BIGLIT_INCR4 =
				(uint32_t) ((char *) p - (char *) table);
	}
	p++;
	}
	}
	}
	}
	}
	return ((double *) p);
}


/* This routine builds an x87 sin/cos table - used by gwsetup */

double *build_x87_sin_cos_table (
	double	*table,		/* Pointer to the table to fill in */
	unsigned long N,
	int	hermetian_skip)	/* True if some sin/cos values are skipped */
{
	unsigned long i;

/* Special case the really small sin/cos tables.  If N is between 9 and 16 */
/* then the assembly code is only doing one FFT level. */
/* In this case, the code just uses the middle sin/cos values of a 2N sized */
/* table.  We could optimize this inefficient memory usage at a later date. */

	if (N <= 8) return (table);
	if (N >= 9 && N <= 16) N = N * 2;

/* The N value passed in represents the number of real numbers that are */
/* processed in a section.  If heremetian_skip is not set, then we are */
/* instead dealing with complex numbers and there are half as many complex */
/* numbers in a section.  For example, when doing 8 levels in pass 2, this */
/* routine is called with N=512.  The first real section has 512 values, */
/* while the remaining pass 2 sections have 256 complex values. */

	if (!hermetian_skip) N = N / 2;

/* Loop to build table */

	for (i = hermetian_skip ? ((N & 4) ? 4 : 8) : 0; i < N; i += 4) {
		unsigned long shifted_i, shifted_N, flipped_i;
		double	sincos[6];

/* Flip the bits in i.  Our prime-factor-FFT makes this a little complex. */
/* The algorithm below works, but I've long since forgotten why. */

		shifted_i = i; shifted_N = N; flipped_i = 0;
		while ((shifted_N & 1) == 0) {
			flipped_i <<= 1;
			if (shifted_i & 1) flipped_i++;
			shifted_i >>= 1;
			shifted_N >>= 1;
		}
		flipped_i = (flipped_i * shifted_N) + shifted_i;

/* When the FFT is working on real data Hermetian symettry allows us to */
/* eliminate half of the FFT data and consequently half of the sin/cos data */
/* Case 1:  If shifted source is > shifted N/2, then we */
/* do not need these sin/cos values. */
/* Case 2:  If shifted source is zero, loop to find the top */
/* two bits.  Skip the number if the top two bits equal 3. */

		if (hermetian_skip) {
			if (shifted_i > shifted_N / 2) continue;
			if (shifted_i == 0) {
				unsigned long j;
				for (j = i; j > 3; j >>= 1);
				if (j == 3) continue;
			}
		}

/* Compute the 3 sin/cos values */

		gwsincos3 (flipped_i, N, (double *) &sincos);

/* Copy the sin/cos values to the table */

		memcpy (table, sincos, sizeof (sincos));
		table += 6;
	}
	return (table);
}

/* This routine builds a pass 2 premultiplier table - used by gwsetup */

double *build_x87_premult_table (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	double	*table,		/* Pointer to the table to fill in */
	unsigned long pass2_size)
{
	unsigned long i, N, incr;

/* Build a premultiplier table for the second pass incrementing by */
/* the pre-calculated pass2_size. */

	N = gwdata->FFTLEN;
	incr = pass2_size;
	if (gwdata->ALL_COMPLEX_FFT) N = N / 2;

/* Mod 2^N+1 arithmetic starts at first data set, */
/* mod 2^N-1 skips some data sets */

	if (gwdata->ALL_COMPLEX_FFT) i = 0;
	else i = incr * 2;

/* Loop to build table */

	for ( ; i < N; i += incr) {
		unsigned long shifted_i, shifted_N, flipped_i, k, l;
		double	sincos[2];

/* Flip the bits in i.  Our prime-factor-FFT makes this a little complex. */
/* The algorithm below works, but I've long since forgotten why. */

		shifted_i = i; shifted_N = N; flipped_i = 0;
		while ((shifted_N & 1) == 0) {
			flipped_i <<= 1;
			if (shifted_i & 1) flipped_i++;
			shifted_i >>= 1;
			shifted_N >>= 1;
		}
		flipped_i = (flipped_i * shifted_N) + shifted_i;

/* When the FFT is working on real data Hermetian symettry allows us to */
/* eliminate half of the FFT data and consequently half of the sin/cos data */
/* Case 1:  If shifted source is > shifted N/2, then we */
/* do not need these sin/cos values. */
/* Case 2:  If shifted source is zero, loop to find the top */
/* two bits.  Skip the number if the top two bits equal 3. */

		if (!gwdata->ALL_COMPLEX_FFT) {
			if (shifted_i > shifted_N / 2) continue;
			if (shifted_i == 0) {
				unsigned long j;
				for (j = i; j > 3; j >>= 1);	
				if (j == 3) continue;
			}
		}

/* Generate the group multipliers */

		for (k = 0; k < incr / 4; k += 4) {

/* There are 4 multipliers in a PMD set */

			for (l = 0; l < 4; l++) {

/* Compute the sin/cos value (root of unity) */

				if (!gwdata->ALL_COMPLEX_FFT) {
					gwsincos (((l * incr/4 + k) * flipped_i) % N, N, (double *) &sincos);
				}

/* If C > 0, then also multiply by the proper root of -1.  This is done */
/* by changing the value we are taking the sin/cos of */

				else {
					gwsincos (((l * incr/4 + k) * flipped_i * 4 + l*incr/4+k) % (N*4), N*4, (double *) &sincos);
				}

/* Save the premultiplier values */

				table[l*2] = sincos[0];
				table[l*2+1] = sincos[1];
			}
			table += 8;
		}
	
/* Generate the 4 column multipliers * 4 sin/cos values */

		for (k = 0; k < 4; k++) {
			for (l = 0; l < 4; l++) {

/* Compute the sin/cos value (root of unity) */

				if (!gwdata->ALL_COMPLEX_FFT) {
					gwsincos ((k * flipped_i + l * N/16) % N, N, (double *) &sincos);
				}

/* If C > 0, then also multiply by the proper root of -1.  This is done */
/* by changing the value we are taking the sin/cos of */

				else {
					gwsincos (((k * flipped_i * 2 + l * N/8) *2 + k) % (N*4), N*4, (double *) &sincos);
				}

/* Save the premultiplier value */

				table[l*2] = sincos[0];
				table[l*2+1] = sincos[1];
			}
			table += 8;
		}
 	}

	return (table);
}

/* This routine builds a plus 1 premultiplier table - used by gwsetup */
/* when c is positive. */

double *build_x87_plus1_table (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	double	*table,		/* Pointer to the table to fill in */
	unsigned long pass2_size)
{
	unsigned long i, k, N;
	int	pfa;

/* Set flag if this is a 3*2^n FFT */

	pfa = (gwdata->FFTLEN != pow_two_above_or_equal (gwdata->FFTLEN));

/* Adjust for two-pass FFTs */

	if (pass2_size == 1) N = gwdata->FFTLEN;
	else N = gwdata->FFTLEN / pass2_size;

/* Loop to build premultiplier table in the same order as the underlying */
/* assembly macro needs them. */

	for (i = 0; i < N / (pfa ? 6 : 8); i++) {
		double	sincos[2];

/* Generate the pre multipliers (roots of -1) used in one three_complex */
/* or four complex macro. */

		for (k = 0; k < (unsigned long) (pfa ? 3 : 4); k++) {
			long	temp;

/* Compute the sin/cos value */

			if (pfa)
				temp = (long) ((i + k * N/6) % N);
			else
				temp = (long) ((i + k * N/8) % N);
			gwsincos (temp, N*2, (double *) &sincos);

/* Save the premultiplier value */

			table[0] = sincos[0];
			table[1] = sincos[1];
			table += 2;
		}
	}

	return (table);
}

/* This routine builds a normalization table - used by x87 normalizaion */
/* routines */

double *build_x87_norm_table (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	double	*table,		/* Pointer to the table to fill in */
	unsigned long pass2_size, /* Size of pass2 */
	int	col)		/* TRUE if building column, not group, table */
{
	unsigned long i, k, num_cols;

/* Handle one-pass FFTs first, there are no group multipliers */

	if (pass2_size == 1) {
		if (!col) return (table);

/* Loop to build table */

		for (i = 0; i < gwdata->FFTLEN; i++) {
			unsigned long j;
			double	ttp, ttmp;

/* Call asm routines to compute the two multipliers */

			gwfft_weights3 (gwdata->dd_data, i, &ttp, NULL, &ttmp);

/* Find where this data appears in the FFT array and in the table we are building. */

			j = addr_offset (gwdata, i) / sizeof (double);

/* Now set the appropriate table entry.  These are put into the array */
/* in the same order that the normalization code needs them. */

			table[j*2] = ttmp;
			table[j*2+1] = ttp;
		}
		return (table + gwdata->FFTLEN + gwdata->FFTLEN);
	}

/* Two pass FFTs are handled here */

	num_cols = pass2_size;
	if (col) {

/* Loop to build columns table */

		for (i = 0; i < num_cols; i++) {
			double	ttp, ttmp;

/* Call asm routines to compute the two multipliers */

			gwfft_weights3 (gwdata->dd_data, i, &ttp, NULL, &ttmp);

/* Now set the appropriate table entry.  These are put into the array */
/* in the same order that the normalization code needs them. */

			table[i+i] = ttmp;
			table[i+i+1] = ttp;
		}
		return (table + num_cols * 2);
	}

/* Build the group multipliers table */

	else {
		unsigned long num_grps;
		
/* Loop to build group table */

		num_grps = gwdata->FFTLEN / num_cols;
		for (i = 0; i < num_grps; i++) {
			double	ttp, ttmp;

/* Call asm routines to compute the two multipliers */

			gwfft_weights3 (gwdata->dd_data, i * num_cols, &ttp, &ttmp, NULL);

/* Now set the appropriate table entry.  These are put into the array */
/* in the same order that the normalization code needs them. */

			if (i < num_grps / 2) k = i * 2;
			else k = (i - num_grps / 2) * 2 + 1;
			table[k+k] = ttmp;
			table[k+k+1] = ttp;
		}
		return (table + num_grps * 2);
	}
}

/* This routine builds a big/little flags table - used by x87 normalizaion */
/* routines */

double *build_x87_biglit_table (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	double	*table,		/* Pointer to the table to fill in */
	unsigned long pass2_size)
{
	unsigned char *p;
	unsigned long i, j, k, m;

/* Handle one pass FFTs differently */

	if (pass2_size == 1) {

/* Loop to build table */

		p = (unsigned char *) table;
		for (i = 0; i < gwdata->FFTLEN; i++) {
			unsigned long table_entry;

/* Find where this data appears in the FFT array and in the table we are building. */

			j = addr_offset (gwdata, i) / sizeof (double);
			table_entry = j >> 1;

/* Now set the biglit table entry for a LSW in a pair */

			if ((j & 1) == 0) {
				p[table_entry] = is_big_word (gwdata, i) * 16;
			}

/* Otherwise, set the biglit table entry for a MSW in a pair */

			else {
				if (is_big_word (gwdata, i)) p[table_entry] += 32;
			}
		}
		return ((double *) (p + gwdata->FFTLEN / 2));
	}

/* Loop to build table in exactly the same order that it will be */
/* used by the assembly code. */

	p = (unsigned char *) table;
	for (i = 0; i < pass2_size; i += gwdata->PASS1_CACHE_LINES * 2) {
	for (j = 0; j < gwdata->FFTLEN / 2; j += pass2_size) {
	for (k = 0; k < gwdata->PASS1_CACHE_LINES * 2; k++) {
	for (m = 0; m < gwdata->FFTLEN; m += gwdata->FFTLEN / 2) {
		unsigned long word;

/* Now set the big/little flag for a LSW in a pair */
/* Otherwise, set the big/little flag for a MSW in a pair */

		word = i + j + k + m;
		if (m == 0) *p = is_big_word (gwdata, word) * 16;
		else if (is_big_word (gwdata, word)) *p += 32;

/* Set the ttp and ttmp fudge flags for two pass FFTs */
/* The fudge flag is set if col mult * grp mult will be greater than 2 */

		if (gwfft_weight_exponent (gwdata->dd_data, word) + 0.5 <
		    gwfft_weight_exponent (gwdata->dd_data, word & (pass2_size-1)) +
		    gwfft_weight_exponent (gwdata->dd_data, word & ~(pass2_size-1))) {
			if (m == 0) *p += 64;
			else *p += 128;
		}

/* Set some offsets that help the assembly code step through the big/lit */
/* array in a non-traditional order.  Two pass-FFTs step through the array */
/* in chunks of PASS1_CACHE_LINES, but the add, sub, and carry propagation */
/* code need to access the big/lit array linearly.  Set two global variables */
/* that tell the assembly code the big/lit array distance between words */
/* 0 and 2, and words 0 and 4. */

		if (word == 2)
			((struct gwasm_data *) gwdata->asm_data)->BIGLIT_INCR2 =
				(uint32_t) ((char *) p - (char *) table);
		if (word == 4)
			((struct gwasm_data *) gwdata->asm_data)->BIGLIT_INCR4 =
				(uint32_t) ((char *) p - (char *) table);
	}
	p++;
	}
	}
	}
	return ((double *) p);
}

/* Return true if half of the words would have the same pattern of big */
/* and little words.  */

int match_pathological_pattern (
	unsigned long num_big_words,
	unsigned long total_length,
	double	pathological_fraction)
{
	double	half_length;
	unsigned long pathological_base, actual_base;

/* Compute the gwfft_base you would get of the word half way into the FFT if */
/* you had the pathological fraction of big and little words */

	half_length = (double) total_length * 0.5;
	pathological_base = (unsigned long) ceil (half_length * pathological_fraction);

/* Compute the base you would get given the actual fraction of big words */

	actual_base = (unsigned long) ceil (half_length * (double) num_big_words / (double) total_length);

/* Return pathological (true) if the actual_base is close to the pathological_base */

	return (actual_base >= pathological_base && actual_base <= pathological_base + 1);
}

/* Here is a particularly nasty routine.  It tries to detect whether the distribution */
/* of big and little words is "pathological".  We want the distribution to be random. */
/* If, for example, there are an equal number of big words and little words then the */
/* every other FFT word consists of big word * big word products, while the other half */
/* contains big word * small word products.  This greatly increases the round off error */
/* especially when b is large (big words are much larger than small words).  This */
/* ugliness was added to handle these cases that where the wrong FFT length was selected: */
/* 211*210^2047-1, 211*210^2687-1, 211*210^7679-1.  There are undoubtedly many others. */

int is_pathological_distribution (
	unsigned long num_big_words,
	unsigned long num_small_words)
{
	unsigned long total_length;

/* Handle cases that we really should never see (rational FFTs) */

	if (num_big_words == 0 || num_small_words == 0) return (FALSE);

/* While the remaining number of big words and small words is even, this */
/* represents a case of a big repeating pattern (the pattern in the upper half *
/* of the remaining words is the same as the pattern in the lower half). */

	total_length = num_big_words + num_small_words;
	while ((num_big_words & 1) == 0 && (total_length & 1) == 0) {
		num_big_words >>= 1;
		total_length >>= 1;
	}

/* The bad patterns occur when the number of big words divided by the FFT length */
/* is close to a small rational number like 1/2, 2/5, 3/4, etc.	 We'll define a */
/* pathological bit pattern as one where more than half of the FFT repeats the */
/* same cycle of big words and small words.  This definition may require some */
/* tweaking over time. */

	if (match_pathological_pattern (num_big_words, total_length, 1.0 / 2.0)) return (TRUE);
	if (match_pathological_pattern (num_big_words, total_length, 1.0 / 3.0)) return (TRUE);
	if (match_pathological_pattern (num_big_words, total_length, 2.0 / 3.0)) return (TRUE);
	if (match_pathological_pattern (num_big_words, total_length, 1.0 / 4.0)) return (TRUE);
	if (match_pathological_pattern (num_big_words, total_length, 3.0 / 4.0)) return (TRUE);
	if (match_pathological_pattern (num_big_words, total_length, 1.0 / 5.0)) return (TRUE);
	if (match_pathological_pattern (num_big_words, total_length, 2.0 / 5.0)) return (TRUE);
	if (match_pathological_pattern (num_big_words, total_length, 3.0 / 5.0)) return (TRUE);
	if (match_pathological_pattern (num_big_words, total_length, 4.0 / 5.0)) return (TRUE);
	if (match_pathological_pattern (num_big_words, total_length, 1.0 / 6.0)) return (TRUE);
	if (match_pathological_pattern (num_big_words, total_length, 5.0 / 6.0)) return (TRUE);
	if (match_pathological_pattern (num_big_words, total_length, 1.0 / 7.0)) return (TRUE);
	if (match_pathological_pattern (num_big_words, total_length, 2.0 / 7.0)) return (TRUE);
	if (match_pathological_pattern (num_big_words, total_length, 3.0 / 7.0)) return (TRUE);
	if (match_pathological_pattern (num_big_words, total_length, 4.0 / 7.0)) return (TRUE);
	if (match_pathological_pattern (num_big_words, total_length, 5.0 / 7.0)) return (TRUE);
	if (match_pathological_pattern (num_big_words, total_length, 6.0 / 7.0)) return (TRUE);
	if (match_pathological_pattern (num_big_words, total_length, 1.0 / 8.0)) return (TRUE);
	if (match_pathological_pattern (num_big_words, total_length, 3.0 / 8.0)) return (TRUE);
	if (match_pathological_pattern (num_big_words, total_length, 5.0 / 8.0)) return (TRUE);
	if (match_pathological_pattern (num_big_words, total_length, 7.0 / 8.0)) return (TRUE);

/* That's all the cases we test for now */

	return (FALSE);
}

/* This routine used to be in assembly language.  It scans the assembly */
/* code arrays looking for the best FFT size to implement our k*b^n+c FFT. */
/* Returns 0 for IBDWT FFTs, 1 for zero padded FFTs, or a gwsetup error */
/* code. */

int gwinfo (			/* Return zero-padded fft flag or error code */
	gwhandle *gwdata,	/* Gwnum global data */
	double	k,		/* K in K*B^N+C. Must be a positive integer. */
	unsigned long b,	/* N in K*B^N+C. */
	unsigned long n,	/* N in K*B^N+C. Exponent to test. */
	signed long c)		/* C in K*B^N+C. Must be rel. prime to K. */
{
	struct gwinfo1_data asm_info;
	struct gwasm_jmptab *jmptab, *zpad_jmptab, *generic_jmptab;
	double	log2k, logbk, log2b, log2c, log2maxmulbyconst;
	double	max_bits_per_input_word, max_bits_per_output_word;
	double	max_weighted_bits_per_output_word;
	int	num_b_in_big_word, num_small_words, num_big_words;
	double	b_per_input_word, bits_per_output_word;
	double	weighted_bits_per_output_word;
	unsigned long l2_cache_size, max_exp;
	int32_t *longp;
	char	buf[20];
	int	qa_nth_fft;
	unsigned long qa_fftlen;

/* If L2 cache size is unknown, assume it is 512KB.  We do this because */
/* it is probably a new processor, and most new processors will have this */
/* much cache. */

	if (CPU_L2_CACHE_SIZE >= 0)
		l2_cache_size = CPU_L2_CACHE_SIZE;
	else
		l2_cache_size = 512;

/* Nehalem has a small L2 cache and a large L3 cache.  We are better off */
/* choosing FFTs designed for larger L2 caches. */

	if (CPU_L3_CACHE_SIZE > 0)
		l2_cache_size += CPU_L3_CACHE_SIZE;

/* Since hyperthreaded CPUs have a shared L2 cache, we pretend we have */
/* half as much L2 cache when we are running two threads. */

	if (CPU_HYPERTHREADS > 1 && gwdata->num_threads > 1)
		l2_cache_size >>= 1;

/* The Celeron D (256K L2 cache) and Willamette P4 (256K L2 cache) have */
/* different optimal FFT implementations.  Adjust the cache size for the */
/* Celeron D so that we can distinguish between the two processors in our */
/* table lookup.  The Celeron D has a 4-way set associative cache, while */
/* the Willamette has an 8-way set-associative cache. */

	if (l2_cache_size == 256 && CPU_L2_SET_ASSOCIATIVE == 4)
		l2_cache_size = 257;

/* Get pointer to 4 assembly jmptables and the version number */

	gwinfo1 (&asm_info);

/* Make sure that the assembly code version number matches the C version */
/* number.  If they do not match, then the user linked in the wrong gwnum */
/* object files! */

	sprintf (buf, "%d.%d", asm_info.version / 100, asm_info.version % 100);
	if (strcmp (buf, GWNUM_VERSION)) return (GWERROR_VERSION);

/* Precalculate some needed values */

	log2k = log2 (k);
	logbk = logb (k);
	log2b = log2 (b);
	log2c = log2 (abs (c));
	log2maxmulbyconst = log2 (gwdata->maxmulbyconst);

/* First, see what FFT length we would get if we emulate the k*b^n+c modulo */
/* with a zero padded FFT.  If k is 1 and abs (c) is 1 then we can skip this */
/* loop as we're sure to find an IBDWT that will do the job. */

again:	zpad_jmptab = NULL;
	generic_jmptab = NULL;
	if (gwdata->specific_fftlen == 0 && (k > 1.0 || n < 500 || abs (c) > 1)) {

/* Use the proper 2^N-1 jmptable */

		if (gwdata->cpu_flags & CPU_SSE2)
			zpad_jmptab = asm_info.p4_cyclic_fft_info;
		else zpad_jmptab =
			zpad_jmptab = asm_info.x86_cyclic_fft_info;

/* Find the table entry for the FFT that can do a mod 2^2n FFT, handling */
/* k and c in the normalization routines.  We will compare this to the */
/* non-zero-padded FFT length later.  The zeroes in the upper half of FFT */
/* input data let us get about another 0.3 bits per input word. */

		qa_nth_fft = 10000;
		while ((max_exp = zpad_jmptab->max_exp) != 0) {

/* If QA'ing different FFT implementations, then skip the ones we've already */
/* tested and the ones that are a different FFT length than the ones already */
/* tested. */

			qa_nth_fft++;
			if (gwdata->qa_pick_nth_fft) {
				if (qa_nth_fft < gwdata->qa_pick_nth_fft)
					goto next1;
				if (qa_nth_fft == gwdata->qa_pick_nth_fft) {
					qa_fftlen = zpad_jmptab->fftlen;
					goto next1;
				}
				if (gwdata->qa_pick_nth_fft != 1 &&
				    qa_fftlen != zpad_jmptab->fftlen)
					goto next1;
				l2_cache_size = 9999999;
			}

/* Do a quick check on the suitability of this FFT */

			if ((double) n * log2b / (double) zpad_jmptab->fftlen > 26.0) goto next1;
			if (zpad_jmptab->fftlen < gwdata->minimum_fftlen) goto next1;

/* Check L2 cache size constraints */

			if (l2_cache_size < ((zpad_jmptab->flags_min_l2_cache_clm >> 16) & 0x3FFF))
				goto next1;

/* Check FFT requires prefetch capability */

			if (zpad_jmptab->flags_min_l2_cache_clm & 0x80000000 &&
			    ! (gwdata->cpu_flags & CPU_PREFETCH))
				goto next1;

/* Check FFT requires prefetchw (3DNow!) capability */

			if (zpad_jmptab->flags_min_l2_cache_clm & 0x40000000 &&
			    ! (gwdata->cpu_flags & CPU_3DNOW))
				goto next1;

/* Don't bother looking at this FFT length if the generic reduction would be faster */

			if (generic_jmptab != NULL &&
			    zpad_jmptab->timing > 3.0 * generic_jmptab->timing)
				goto next1;

/* See if this is the FFT length that would be used for a generic modulo reduction */

			if (generic_jmptab == NULL &&
			    2.0 * (log2k + n * log2b) + 128.0 < max_exp + 0.3 * zpad_jmptab->fftlen)
				generic_jmptab = zpad_jmptab;

/* Compare the maximum number of bits allowed in the FFT input word */
/* with the number of bits we would use.  Break when we find an acceptable */
/* FFT length. */
//  This is the old code which only supported b == 2
//			max_bits_per_word = (double) max_exp / zpad_jmptab->fftlen;
//			max_bits_per_word -= gwdata->safety_margin;
//			bits_per_word = (double) (n + n) * log2b / zpad_jmptab->fftlen;
//			if (bits_per_word < max_bits_per_word + 0.3) {
//				if (gwdata->qa_pick_nth_fft) goto use_zpad;
//				break;
//			}

/* In version 25.11, we need to handle b != 2.  See comments later on in this routine */
/* for a description of the concepts involved. */

/* Compute the maximum number of bits allowed in the FFT input word */

			max_bits_per_input_word = (double) max_exp / zpad_jmptab->fftlen;
			max_bits_per_input_word -= gwdata->safety_margin;

/* Apply our new formula (described later) to the maximum Mersenne exponent for this FFT length. */

			num_b_in_big_word = (int) ceil (max_bits_per_input_word);
			num_small_words = (int) ((num_b_in_big_word - max_bits_per_input_word) * zpad_jmptab->fftlen);
			num_big_words = zpad_jmptab->fftlen - num_small_words;
			max_bits_per_output_word =
				2 * (num_b_in_big_word - 1) +
				0.6 * log2 (num_big_words + num_small_words / 3.174802103936252);

/* Apply our new formula (described later) to the number we are planning to test.  */
/* This is different for the zero-pad case because only 4 words in the upper half */
/* of the FFT contain any data.  We can't use the FFT length if the k value will */
/* not fit in 4 words. */

			b_per_input_word = (double) (n + n) / zpad_jmptab->fftlen;
			if (logbk > 4.0 * b_per_input_word) goto next1;
			num_b_in_big_word = (int) ceil (b_per_input_word);
			num_small_words = (int) ((num_b_in_big_word - b_per_input_word) * (zpad_jmptab->fftlen / 2 + 4));
			num_big_words = (zpad_jmptab->fftlen / 2 + 4) - num_small_words;
			bits_per_output_word =
				2.0 * (num_b_in_big_word * log2b - 1.0) +
				0.6 * log2 (num_big_words + num_small_words / pow (2.0, log2b / 0.6));

/* And compute the weighted values as per the formulas described later */

			max_weighted_bits_per_output_word =
				2.0 * max_bits_per_input_word + 0.6 * log2 (zpad_jmptab->fftlen / 2 + 4);
			weighted_bits_per_output_word =
				2.0 * ((b_per_input_word + 1.0) * log2b - 1.0) +
				0.6 * log2 (zpad_jmptab->fftlen / 2 + 4);
			if ((n + n) % zpad_jmptab->fftlen == 0)
				weighted_bits_per_output_word -= ((log2b <= 4.0) ? log2b : 1.4 * log2b);
			else if (! is_pathological_distribution (num_big_words, num_small_words))
				weighted_bits_per_output_word -=
					((log2b <= 3.0) ? (log2b - 1.0) / 2.0 :
					 (log2b <= 6.0) ? 1.0 + (log2b - 3.0) / 3.0 :
							  2.0 + (log2b - 6.0) / 6.0);

/* See if this FFT length might work */

			if (weighted_bits_per_output_word <= max_weighted_bits_per_output_word &&

/* Result words are multiplied by k and the mul-by-const and any carry spread over 6 words. */
/* Thus, the multiplied FFT result word cannot be more than 7 times bits-per-input-word */
/* (bits-per-input-word are stored in the current word and the 6 words we propogate carries to). */

			    bits_per_output_word + log2k + log2maxmulbyconst <=
						    floor (7.0 * b_per_input_word) * log2b &&

/* The high part of upper result words are multiplied by c and the mul-by-const.  This must */
/* not exceed 51 bits. */

			    bits_per_output_word - floor (b_per_input_word) * log2b +
						    log2c + log2maxmulbyconst <= 51.0) {

/* We can use this FFT.  Use it if specifically requested, or */
/* look for a non-zero-padded FFT that might be even faster. */

				if (gwdata->qa_pick_nth_fft) goto use_zpad;
				break;
			}

/* Move to next jmptable entry */

next1:			for (longp = &zpad_jmptab->counts[4]; *longp; longp++);
			zpad_jmptab = (struct gwasm_jmptab *) (longp + 1);
		}
	}

/* Now see what FFT length we would use if a DWT does the k*b^n+c modulo. */

/* Use the proper 2^N+1 or 2^N-1 jmptable */

	if (c < 0) {
		if (gwdata->cpu_flags & CPU_SSE2)
			jmptab = asm_info.p4_cyclic_fft_info;
		else
			jmptab = asm_info.x86_cyclic_fft_info;
	} else {
		if (gwdata->cpu_flags & CPU_SSE2)
			jmptab = asm_info.p4_negacyclic_fft_info;
		else
			jmptab = asm_info.x86_negacyclic_fft_info;
	}

/* Find the table entry using either the specified fft length or */
/* the smallest FFT length that can handle the k,b,n,c being tested. */

	qa_nth_fft = 50000;
	while ((max_exp = jmptab->max_exp) != 0) {

/* If QA'ing different FFT implementations, then skip the ones we've */
/* already tested and the ones that are a different FFT length than the */
/* ones already tested. */

		qa_nth_fft++;
		if (gwdata->qa_pick_nth_fft) {
			if (qa_nth_fft < gwdata->qa_pick_nth_fft)
				goto next2;
			if (qa_nth_fft == gwdata->qa_pick_nth_fft) {
				qa_fftlen = jmptab->fftlen;
				goto next2;
			}
			if (gwdata->qa_pick_nth_fft >= 50000 &&
			    qa_fftlen != jmptab->fftlen)
				goto next2;
			l2_cache_size = 9999999;
		}

/* Do a quick check on the suitability of this FFT */

		if ((double) n * log2b / (double) jmptab->fftlen > 26.0) goto next2;
		if (jmptab->fftlen < gwdata->minimum_fftlen) goto next2;

/* Check if FFT requires prefetch capability */

		if (jmptab->flags_min_l2_cache_clm & 0x80000000 &&
		    ! (gwdata->cpu_flags & CPU_PREFETCH))
			goto next2;

/* Check if FFT requires prefetchw (3DNow!) capability */

		if (jmptab->flags_min_l2_cache_clm & 0x40000000 &&
		    ! (gwdata->cpu_flags & CPU_3DNOW))
			goto next2;

/* Handle benchmarking case that selects the nth FFT implementation */
/* regardless of cache size considerations. */

		if (gwdata->bench_pick_nth_fft) {
			l2_cache_size = 9999999;
			if (gwdata->specific_fftlen != jmptab->fftlen)
				goto next2;
			if (gwdata->cpu_flags & CPU_3DNOW &&
			    ! (jmptab->flags_min_l2_cache_clm & 0x40000000))
				goto next2;
			if (--gwdata->bench_pick_nth_fft) goto next2;
		}

/* Check L2 cache size constraints */

		if (l2_cache_size < ((jmptab->flags_min_l2_cache_clm >> 16) & 0x3FFF))
			goto next2;

/* Check if this table entry matches the specified FFT length. */

		if (gwdata->specific_fftlen) {
			if (gwdata->specific_fftlen == jmptab->fftlen) break;
		}

/* Or check that this FFT length will work with this k,n,c combo */

		else {
//  This is the old code which only supported b == 2
//			double max_bits_per_word;
//			double bits_per_word;
//
/* Compute the maximum number of bits allowed in the FFT input word */
//
//			max_bits_per_word = (double) max_exp / jmptab->fftlen;
//			max_bits_per_word -= gwdata->safety_margin;
//
/* For historical reasons, the jmptable computes maximum exponent based on */
/* a Mersenne-mod FFT (i.e k=1.0, c=-1).  Handle more complex cases here. */
/* A Mersenne-mod FFT produces 2 * bits_per_word in each FFT result word. */
/* The more general case yields 2 * bits_per_word + log2(k) + 1.7 * log2(c) */
/* in each FFT result word.  NOTE: From the data I've observed, doubling c */
/* about triples the roundoff error (that is log2(3) = 1.585 output bits). */
/* However, when I used 1.585 in the formula it was not hard to find cases */
/* where the roundoff error was too high, so we use 1.7 here for extra */
/* safety. */
//
//			bits_per_word = (log2k + n * log2b) / jmptab->fftlen;
//			if (2.0 * bits_per_word + log2k + 1.7 * log2c <=
//					2.0 * max_bits_per_word) {
/* Because carries are spread over 4 words, there is a minimum limit on */
/* the bits per word.  An FFT result word cannot be more than 5 times */
/* bits-per-word (bits-per-word are stored in the current word and the */
/* 4 words we propogate carries to).  How many bits are in an FFT result */
/* word?  Well, because of balanced representation the abs(input word) is */
/* (bits_per_word-1) bits long. An FFT result word contains multiplied data */
/* words, that's (bits_per_word-1)*2 bits.  Adding up many multiplied data */
/* words adds some bits proportional to the size of the FFT.  Experience */
/* has shown this to be 0.6 * log (FFTLEN).  This entire result is */
/* multiplied by k in the normalization code, so add another log2(k) bits. */
/* Finally, the mulbyconst does not affect our chance of getting a round off */
/* error, but does add to the size of the carry. */
//
//			loglen = log2 (jmptab->fftlen);
//			total_bits = (bits_per_word - 1.0) * 2.0 +
//				     1.7 * log2c + loglen * 0.6 +
//				     log2k + log2maxmulbyconst;
//			if (total_bits > 5.0 * bits_per_word) {

/* In version 25.11, we now need to handle b != 2.  Consider the case */
/* where b is ~4000.  If small words contain one b (~12 bits) and large words */
/* contain two b (~24 bits), then the number of bits in a result word is */
/* dominated by big words * big word products (~48 bits).  The old code above */
/* tested average bits per word (~18 bits) and underestimates a result word as */
/* containing ~36 bits.  So here's our upgraded model.  We calculate the number */
/* of big and little words.  A result word adds up many big word times big word */
/* products and big word times small word products.  Let base = b ^ num_b_in_big_word. */
/* Because of balanced representation, a big word times big word */
/* product is 2 * (log2(base) - 1) bits.  Summing them up adds about */
/* 0.6 * log2 (num_big_words) more bits.  Now for the big words times */
/* small words products that are also added in, the more bits in a small word the more */
/* it impacts the result word.  A big word times small word product has log2(b) fewer */
/* bits in it.  If we add two big word times small word products, the sum is */
/* about 0.6 bits bigger, add four products to get 1.2 bits bigger, etc. -- do this until you */
/* overcome the log2(b) bit difference.  That is, 2^(log2(b)/0.6) small */
/* products equals one big product.  Putting it all together, a result word contains */
/* 2 * (log2(base) - 1) + 0.6 * log2 (num_big_words + num_small_words / 2^(log2(b)/0.6)) */
/* bits plus the k and c adjustments noted above. */

/* Compute the maximum number of bits allowed in the FFT input word */

			max_bits_per_input_word = (double) max_exp / jmptab->fftlen;
			max_bits_per_input_word -= gwdata->safety_margin;

/* Apply our new formula above to the maximum Mersenne exponent for this FFT length. */

			num_b_in_big_word = (int) ceil (max_bits_per_input_word);
			num_small_words = (int) ((num_b_in_big_word - max_bits_per_input_word) * jmptab->fftlen);
			num_big_words = jmptab->fftlen - num_small_words;
			max_bits_per_output_word =
				2 * (num_b_in_big_word - 1) +
				0.6 * log2 (num_big_words + num_small_words / 3.174802103936252);

/* Apply our new formula to the number we are planning to test */

			b_per_input_word = (logbk + n) / jmptab->fftlen;
			num_b_in_big_word = (int) ceil (b_per_input_word);
			num_small_words = (int) ((num_b_in_big_word - b_per_input_word) * jmptab->fftlen);
			num_big_words = jmptab->fftlen - num_small_words;
			bits_per_output_word = 
				2.0 * (num_b_in_big_word * log2b - 1.0) +
				0.6 * log2 (num_big_words + num_small_words / pow (2.0, log2b / 0.6)) +
				log2k + 1.7 * log2c;

/* Unfortunately, the story does not end there.  The weights applied to each FFT word */
/* range from 1 to b.  These extra bits impact the round off error.  Thus, we calculate */
/* the weighted_bits_per_output_word for irrational FFTs as using another log2b bits. */

			max_weighted_bits_per_output_word = 2.0 * max_bits_per_input_word + 0.6 * log2 (jmptab->fftlen);
			weighted_bits_per_output_word =
				2.0 * ((b_per_input_word + 1.0) * log2b - 1.0) +
				0.6 * log2 (jmptab->fftlen) + log2k + 1.7 * log2c;

/* Also, testing shows that for small b an unweighted FFT saves about */
/* log2b output bits, and for larger b saves about 1.4 * log2b output bits. */

			if (k == 1.0 && n % jmptab->fftlen == 0)
				weighted_bits_per_output_word -= ((log2b <= 4.0) ? log2b : 1.4 * log2b);

/* A pathological case occurs when num_big_words is one and k is greater than one. */
/* The FFT weights for the small words will not range from 1 to b.  Depending on the */
/* fractional part of logb(k).  In the worst case scenario, the small word weights */
/* range from b - epsilon to b.  The example that raised this issue is 28*3^12285-1. */

			else if (num_big_words == 1 && k > 1.0)
				weighted_bits_per_output_word += log2b;
				
/* Furthermore, testing shows us that larger b values don't quite need the full log2b */
/* bits added (except for some pathological cases), probably because there are fewer */
/* extra bits generated by adding products because the smallest weighted words have */
/* fewer bits.  The correction is if log2b is 3 you can get 1 more output bit than */
/* expected, if log2b is 6 you get about 2 extra bits, if log2b is 12 you can get */
/* 3 extra bits. */

			else if (! is_pathological_distribution (num_big_words, num_small_words))
				weighted_bits_per_output_word -=
					((log2b <= 3.0) ? (log2b - 1.0) / 2.0 :
					 (log2b <= 6.0) ? 1.0 + (log2b - 3.0) / 3.0 :
							  2.0 + (log2b - 6.0) / 6.0);

/* If the bits in an output word is less than the maximum allowed, we can */
/* probably use this FFT length -- though we need to do a few more tests. */

			if (weighted_bits_per_output_word <= max_weighted_bits_per_output_word) {

/* Because carries are spread over 4 words, there is a minimum limit on */
/* the bits per word.  An FFT result word cannot be more than 5 times */
/* bits-per-input-word (bits-per-input-word are stored in the current */
/* word and the 4 words we propogate carries to).  The mul-by-const that */
/* during the normalization process adds to the size of the carry. */

				if (bits_per_output_word + log2maxmulbyconst >
							    floor (5.0 * b_per_input_word) * log2b) {
// This assert was designed to find any cases where using 5 or more carry words
// would use a shorter FFT than using a zero-padded FFT.  We did find a case:
// 15539*2^15095288+7 would use 1.5M FFT length if I implemented 5 carry words
// and requires a 1.75M FFT length for a zero-padded FFT.
//					ASSERTG (zpad_jmptab == NULL ||
//						 jmptab->fftlen >= zpad_jmptab->fftlen);
					goto next2;
				}

/* Because of limitations in the top_carry_adjust code, there is a limit */
/* on the size of k that can be handled.  This isn't a big deal since the */
/* zero-padded implementation should use the same FFT length.  Check to see */
/* if this k can be handled.  K must fit in the top three words for */
/* one-pass FFTs and within the top two words of two-pass FFTs. */

				if (jmptab->pass2_levels == 0 && logbk > floor (3.0 * b_per_input_word)) {
// This assert is designed to find any cases where using 3 or more carray adjust words
// would use a shorter FFT than using a zero-padded FFT.
					ASSERTG (zpad_jmptab != NULL && jmptab->fftlen >= zpad_jmptab->fftlen);
					goto next2;
				}
				if (jmptab->pass2_levels != 0 && logbk > floor (2.0 * b_per_input_word)) {
// This assert is designed to find any cases where using 3 or more carray adjust words
// would use a shorter FFT than using a zero-padded FFT.
					ASSERTG (zpad_jmptab != NULL && jmptab->fftlen >= zpad_jmptab->fftlen);
					goto next2;
				}

/* We've found an FFT length to use */

				break;
			}
		}

/* Move to next jmptable entry */

next2:		for (longp = &jmptab->counts[4]; *longp; longp++);
		jmptab = (struct gwasm_jmptab *) (longp + 1);
	}

/* If the zero pad FFT length is less than the DWT FFT length, then use */
/* the zero pad FFT length. */

	if (zpad_jmptab != NULL && zpad_jmptab->max_exp &&
	    (jmptab->max_exp == 0 || zpad_jmptab->fftlen < jmptab->fftlen)) {
use_zpad:	gwdata->ZERO_PADDED_FFT = TRUE;
		gwdata->ALL_COMPLEX_FFT = FALSE;
		jmptab = zpad_jmptab;
	}

/* If we found a DWT table entry then return the address in fft_info. */

	else if (jmptab->max_exp) {
		gwdata->ZERO_PADDED_FFT = FALSE;
		gwdata->ALL_COMPLEX_FFT = (c > 0);
	}

/* Error - neither method could handle this huge number */

	else
		return (GWERROR_TOO_LARGE);

/* See if the user requested a larger than normal FFT size */

	if (gwdata->larger_fftlen_count) {
		gwdata->larger_fftlen_count--;
		gwdata->minimum_fftlen = jmptab->fftlen + 1;
		goto again;
	}

/* Remember the FFT length and number of pass 2 levels and the gap */
/* between 2 blocks.  This is used by addr_offset called from */
/* gwmap_to_estimated_size without a call to gwsetup. */

	gwdata->FFTLEN = jmptab->fftlen;
	gwdata->PASS2_LEVELS = jmptab->pass2_levels; /* FFT levels in pass2 */
	gwdata->PASS2GAPSIZE = jmptab->counts[1] - 1;
	gwdata->GW_ALIGNMENT = 4096;	/* Guess an alignment so gwsize can */
					/* return a reasonable value for */
					/* large page support in gwsetup */

/* Remember the FFT returned so that we can return a different FFT to */
/* the QA code next time. */

	gwdata->qa_picked_nth_fft = qa_nth_fft;

/* Set pointer to assembler "jump" table */

	gwdata->jmptab = jmptab;
	return (0);
}


/* Initialize gwhandle for a future gwsetup call. */
/* The gwinit function has been superceeded by gwinit2.  By passing in the */
/* version number we can verify the caller used the same gwnum.h file as the */
/* one he eventually links with.  The sizeof (gwhandle) structure is used */
/* to verify he compiles with the same structure alignment options that */
/* were used when compiling gwnum.c.  For compatibility with existing code */
/* we delay reporting any compatibility problems until gwsetup is called. */

void gwinit2 (
	gwhandle *gwdata,	/* Placeholder for gwnum global data */
	int	struct_size,	/* Size of the gwdata structure */
	char	*version_string)
{

/* See if caller is using the same gwnum.h file that was used when */
/* this file was compiled.  Checking structure size also verifies he */
/* used the same compiler switches - especially regarding alignment. */
/* As a hack, we use GWERROR to record delayed error messages. */

	if (strcmp (version_string, GWNUM_VERSION)) {
		gwdata->GWERROR = GWERROR_VERSION_MISMATCH;
		return;
	}
	if (struct_size != sizeof (gwhandle)) {
		gwdata->GWERROR = GWERROR_STRUCT_SIZE_MISMATCH;
		return;
	}

/* Initialize gwhandle structure with the default values */

	memset (gwdata, 0, sizeof (gwhandle));
	gwdata->safety_margin = 0.0;
	gwdata->maxmulbyconst = 3;
	gwdata->specific_fftlen = 0;
	gwdata->minimum_fftlen = 0;
	gwdata->larger_fftlen_count = 0;
	gwdata->num_threads = 1;
	gwdata->force_general_mod = 0;
	gwdata->use_irrational_general_mod = 0;
	gwdata->use_large_pages = 0;

/* Init structure that allows giants and gwnum code to share */
/* allocated memory */
	
	init_ghandle (&gwdata->gdata);
	gwdata->gdata.allocate = &gwgiantalloc;
	gwdata->gdata.free = &gwgiantfree;
	gwdata->gdata.deallocate = &gwgiantdealloc;
	gwdata->gdata.handle = (void *) gwdata;

/* If CPU type and speed have not been initialized by the caller, do so now. */

	if (CPU_FLAGS == 0 && CPU_SPEED == 0.0) {
		guessCpuType ();
		guessCpuSpeed ();
	}
	gwdata->cpu_flags = CPU_FLAGS;
}

/* Allocate memory and initialize assembly code for arithmetic */
/* modulo k*b^n+c */

int gwsetup (
	gwhandle *gwdata,	/* Placeholder for gwnum global data */
	double	k,		/* K in K*B^N+C. Must be a positive integer. */
	unsigned long b,	/* B in K*B^N+C. */
	unsigned long n,	/* N in K*B^N+C. Exponent to test. */
	signed long c)		/* C in K*B^N+C. */
{
	int	gcd, error_code, setup_completed;
	double	orig_k;
	unsigned long orig_n;

/* Return delayed errors from gwinit2 */

	if (gwdata->GWERROR) return (gwdata->GWERROR);

/* Sanity check the k value */

	if (k < 1.0) return (GWERROR_K_TOO_SMALL);
	if (k > 9007199254740991.0) return (GWERROR_K_TOO_LARGE);

/* Init */

	setup_completed = FALSE;
	orig_k = k;
	orig_n = n;

/* Our code fails if k is a power of b.  For example, 3481*59^805-1 which */
/* equals 59^807-1.  I think this is because gwfft_base(FFTLEN) is off by one */
/* because even quad-precision floats won't calculate FFTLEN * num_b_per_word */
/* correctly.  There is an easy fix, if k is divisible by b we divide k by b */
/* and add one to n. */

	while (k > 1.0 && b > 1 && fmod (k, (double) b) == 0.0) {
		k = k / (double) b;
		n = n + 1;
	}

/* Our code fast code fails if k and c are not relatively prime.  This */
/* is because we cannot calculate 1/k.  Although the user shouldn't call */
/* us with this case, we handle it anyway by reverting to the slow general */
/* purpose multiply routines. */

	if (c == 0)
		gcd = 0;
	else if (k == 1.0 || abs (c) == 1)
		gcd = 1;
	else {
		stackgiant(kg,2);
		stackgiant(cg,2);
		dbltog (k, kg);
		itog (abs (c), cg);
		gcdg (kg, cg);
		gcd = cg->n[0];
	}

/* Call the internal setup routine when we can.  Gcd (k, c) must be 1, */
/* k * mulbyconst and c * mulbyconst cannot be too large.  Also, the FFT */
/* code has bugs when there are too few bits per FFT.  Rather than make */
/* difficult fixes we simply force these small numbers to use the generic */
/* reduction.  In truth, the caller should use a different math package for */
/* these small numbers. */

	if (gcd == 1 &&
	    k * gwdata->maxmulbyconst <= MAX_ZEROPAD_K &&
	    abs (c) * gwdata->maxmulbyconst <= MAX_ZEROPAD_C &&
	    log2(b) * (double) n >= 350.0 &&
	    (b == 2 || (gwdata->cpu_flags & CPU_SSE2)) &&
	    !gwdata->force_general_mod) {
		error_code = internal_gwsetup (gwdata, k, b, n, c);
		if (error_code == 0) setup_completed = TRUE;
		else if (b == 2) return (error_code);
		gwdata->GENERAL_MOD = FALSE;
	}

/* Emulate k not relatively prime to c, small n values, and */
/* large k or c values with a call to the general purpose modulo setup code. */

	if (!setup_completed) {
		double	bits;
		giant	g;

		bits = (double) n * log2 (b);
		g = allocgiant (((unsigned long) bits >> 5) + 4);
		if (g == NULL) return (GWERROR_MALLOC);
		ultog (b, g);
		power (g, n);
		dblmulg (k, g);
		iaddg (c, g);
		gwdata->force_general_mod = TRUE;
		error_code = gwsetup_general_mod_giant (gwdata, g);
		free (g);
		if (error_code) return (error_code);
	}

/* For future messages, format the input number as a string */

	gw_as_string (gwdata->GWSTRING_REP, orig_k, b, orig_n, c);

/* Return success */

	return (0);
}

/* These setup routines are for operations modulo an arbitrary binary number. */
/* This is three times slower than the special forms above. */

int gwsetup_general_mod (
	gwhandle *gwdata,	/* Placeholder for gwnum global data */
	const uint32_t *array,	/* The modulus as an array of 32-bit values */
	uint32_t arraylen)	/* Number of values in the array */
{
	giantstruct tmp;
	tmp.sign = arraylen;
	tmp.n = (uint32_t *) array;
	while (tmp.sign && tmp.n[tmp.sign-1] == 0) tmp.sign--;
	return (gwsetup_general_mod_giant (gwdata, &tmp));
}

int gwsetup_general_mod_64 (
	gwhandle *gwdata,	/* Placeholder for gwnum global data */
	const uint64_t *array,	/* The modulus as an array of 64-bit values */
	uint64_t arraylen)	/* Number of values in the array */
{
	giantstruct tmp;
	tmp.sign = (int) arraylen * 2;
	tmp.n = (uint32_t *) array;
	while (tmp.sign && tmp.n[tmp.sign-1] == 0) tmp.sign--;
	return (gwsetup_general_mod_giant (gwdata, &tmp));
}

/* Setup the FFT code for generic reduction */

int gwsetup_general_mod_giant (
	gwhandle *gwdata,	/* Placeholder for gwnum global data */
	giant	g)		/* The modulus */
{
	unsigned long bits;	/* Bit length of modulus */
	int	convertible;	/* Value can be converted to (k*2^n+c)/d */
	double	k;
	unsigned long n;
	signed long c;
	unsigned long d;
	unsigned long safety_bits;
	struct gwasm_jmptab *info;
	int	error_code;
	unsigned long fftlen, max_exponent, desired_n;
	giant	modified_modulus, tmp;

/* Return delayed errors from gwinit2 */

	if (gwdata->GWERROR) return (gwdata->GWERROR);

/* Init */

	bits = bitlen (g);

/* Examine the giant to see if it a (k*2^n+c)/d value that we can better optimize */
/* Also detect nasty bit patterns, like Phi (82730,2), where multiplying by a small d */
/* results in a less nasty bit pattern for the modulus. */

	d = 1;
	convertible = (!gwdata->force_general_mod &&
		       bits > 300 &&
		       convert_giant_to_k2ncd (g, &k, &n, &c, &d));

/* Copy the modulus except for k*2^n+c values */

	if (!convertible || d != 1) {
		gwdata->GW_MODULUS = allocgiant ((bits >> 5) + 1);
		if (gwdata->GW_MODULUS == NULL) {
			gwdone (gwdata);
			return (GWERROR_MALLOC);
		}
		gtog (g, gwdata->GW_MODULUS);
	}

/* Setup for values we are converting to use faster k*2^n+c FFTs. */

	if (convertible) {
		error_code = gwsetup (gwdata, k, 2, n, c);
		if (error_code) return (error_code);
		if (d != 1) {
			char	buf[60];
			strcpy (buf, gwdata->GWSTRING_REP);
			if (isdigit (buf[0]))
				sprintf (gwdata->GWSTRING_REP, "(%s)/%lu", buf, d);
			else
				sprintf (gwdata->GWSTRING_REP, "%s/%lu", buf, d);
		}
		return (0);
	}

/* If we need to multiply the modulus by a small d value, do so here */

	if (d != 1) {
		modified_modulus = allocgiant ((bits >> 5) + 2);
		if (modified_modulus == NULL) {
			gwdone (gwdata);
			return (GWERROR_MALLOC);
		}
		gtog (g, modified_modulus);
		ulmulg (d, modified_modulus);
		g = modified_modulus;
		bits = bitlen (g);
	} else
		modified_modulus = NULL;

/* We will need twice the number of input bits plus some padding */

	n = bits + bits + 128;

/* Setup the FFT code in much the same way that gwsetup_without_mod does. */
/* Unless the user insists, we try for an integral number of bits per word. */
/* There are pathological bit patterns that generate huge roundoff errors. */
/* For example, if we test (10^828809-1)/9 and put exactly 18 bits into */
/* each FFT word, then every FFT word in GW_MODULUS_FFT will contain the */
/* same value!  Not exactly, the random data our FFTs require for small */
/* roundoff errors.  Thus, the caller may need to insist we use an */
/* irrational FFT on occasion. */

/* Call gwinfo and have it figure out the FFT length we will use. */
/* Since we zero the upper half of FFT input data, the FFT */
/* outputs will be smaller.  This lets us get about another 0.3 bits */
/* per input word. */

	gwdata->safety_margin -= 0.3;
	error_code = gwinfo (gwdata, 1.0, 2, n, -1);
	gwdata->safety_margin += 0.3;
	if (error_code) return (error_code);
	info = gwdata->jmptab;
	fftlen = info->fftlen;
	max_exponent = info->max_exp;

/* Our FFTs don't handle cases where there are few bits per word because */
/* carries must be propagated over too many words.  Arbitrarily insist */
/* that n is at least 12 * fftlen.  */

	if (n < 12 * fftlen) n = 12 * fftlen;

/* Let the user request rational FFTs as they are a few percent faster */

	if (!gwdata->use_irrational_general_mod) {

/* If possible, increase n to the next multiple of FFT length.  This is */
/* because rational FFTs are faster than irrational FFTs (no FFT weights). */

		desired_n = ((n + fftlen - 1) / fftlen) * fftlen;
		if (desired_n < max_exponent) n = desired_n;
	}

/* If the user requested irrational FFTs, then make sure the bits */
/* per FFT word will distribute the big and little words of the modulus */
/* semi-randomly.  For example, in the (10^828809-1)/9 case above, if */
/* bits-per-word is 18.5 or 18.25 you will still get non-random patterns */
/* in the FFT words. */

	else {
		double	prime_number, bits_per_word;

/* Round bits_per_word up to the next half-multiple of 1/prime_number */

		prime_number = 53.0;
		bits_per_word = (double) n / (double) fftlen;
		bits_per_word = (ceil (bits_per_word * prime_number) + 0.5)/ prime_number;

/* If possible, use the n associated with the just-computed bits-per-word */

		desired_n = (unsigned long) ceil (bits_per_word * (double) fftlen);
		if (desired_n < max_exponent) n = desired_n;
	}

/* If possible, increase n to the next multiple of FFT length. */
/* The extra bits allow gwsmallmul to avoid emulate_mod calls more often. */
/* We hope the 0.3 safety_limit increase above will avoid getting too */
/* close to the FFT limit as many users of this library turn on error */
/* checking (slower) when near the FFT limit. */
/* If that doesn't work, try adding a half FFT length instead. */

	if (n + fftlen < max_exponent)
		n = n + fftlen;
	else if (gwdata->use_irrational_general_mod && n + fftlen / 2 < max_exponent)
		n = n + fftlen / 2;

/* Now setup the assembly code */

	gwdata->safety_margin -= 0.3;
	error_code = internal_gwsetup (gwdata, 1.0, 2, n, -1);
	gwdata->safety_margin += 0.3;
	if (error_code) return (error_code);

// BUG - setting the bit_length to the modulus size will break gwtogiant.
// we need a better/more-consistent way of dealing with the various 
// needed bit_lengths.  Also, PFGW should not be reading the bit_length
// value in integer.cpp.
//	gwdata->bit_length = bits;
	
/* Allocate memory for an FFTed copy of the modulus. */

	gwdata->GW_MODULUS_FFT = gwalloc (gwdata);
	if (gwdata->GW_MODULUS_FFT == NULL) {
		gwdone (gwdata);
		return (GWERROR_MALLOC);
	}
	gianttogw (gwdata, g, gwdata->GW_MODULUS_FFT);
	gwfft (gwdata, gwdata->GW_MODULUS_FFT, gwdata->GW_MODULUS_FFT);

/* Calculate number of words to zero during the copy prior to calculating */
/* the quotient. */

	gwdata->GW_ZEROWORDSLOW = (unsigned long)
			floor ((double) bits / gwdata->avg_num_b_per_word);

/* A quick emulate_mod refresher: */
/* 1) The maximum size of a quotient is gwdata->n/2 bits due to the */
/*    zeroing of high words in the normalization routine.  Obviously the */
/*    reciprocal needs to be accurate to at least gwdata->n/2 bits. */
/* 2) So that the quotient doesn't overflow, the maximum size of a value */
/*    entering emulate_mod is gwdata->n/2+bits bits */
/* 3) So that gwsquare and gwmul don't send emulate_mod a value that is */
/*    too large, the maximum input to these routines should be (allowing */
/*    for an 8 bit mulbyconstant) is (gwdata->n/2+bits-8)/2 bits.  This is */
/*    used by gwsmallmul to know when an emulate_mod is required */
/* 4) We cannot quite push to the limits calculated above because we have to */
/*    make sure the quotient calculation does not produce more than gwdata->n */
/*    bits of result -- otherwise the *high* order bits of the quotient will be */
/*    corrupted.  We allow ourselves a small safety margin by decreasing the */
/*    number of reciprocal bits calculated.  The safety margin must be larger */
/*    than the number of unzeroed bits caused by using the floor function in */
/*    calculating GW_ZEROWORDSLOW. */

	safety_bits = bits - (unsigned long) ((double) gwdata->GW_ZEROWORDSLOW * gwdata->avg_num_b_per_word) + 3;

/* Precompute the reciprocal of the modified modulus. */

	gwdata->GW_RECIP_FFT = gwalloc (gwdata);
	if (gwdata->GW_RECIP_FFT == NULL) {
		gwdone (gwdata);
		return (GWERROR_MALLOC);
	}
	tmp = allocgiant ((gwdata->n >> 5) + 1);
	if (tmp == NULL) {
		gwdone (gwdata);
		return (GWERROR_MALLOC);
	}
	itog (1, tmp);
	gshiftleft (bits + gwdata->n / 2 - safety_bits, tmp);
	divg (g, tmp);		/* computes gwdata->n/2-safety_margin+1 bits of reciprocal */
	gshiftleft (gwdata->n - (bits + gwdata->n / 2 - safety_bits), tmp);
				/* shift so gwmul routines wrap */
				/* quotient to lower end of fft */
	gianttogw (gwdata, tmp, gwdata->GW_RECIP_FFT);
	gwfft (gwdata, gwdata->GW_RECIP_FFT, gwdata->GW_RECIP_FFT);
	free (tmp);
	free (modified_modulus);

/* Calculate the maximum allowable size of a number used as input */
/* to gwmul.  We will make sure gwsmallmul does not generate any */
/* results bigger than this. */

	gwdata->GW_GEN_MOD_MAX = (unsigned long)
		 floor ((double)((gwdata->n/2-safety_bits+bits-8)/2) / gwdata->avg_num_b_per_word);
	gwdata->GW_GEN_MOD_MAX_OFFSET = addr_offset (gwdata, gwdata->GW_GEN_MOD_MAX-1);

/* Set flag indicating general-purpose modulo operations are in force */

	gwdata->GENERAL_MOD = TRUE;

/* It appears that when we multiply the modulus by a small d value, we are dealing */
/* with bit patterns that may generate a larger than usual SUM(INPUTS)/SUM(OUTPUTS) */
/* difference.  The test case is Phi(109965,2).  Thus we will increase MAXDIFF.	*/
	
	if (d != 1) gwdata->MAXDIFF *= 100.0;
		
/* Create dummy string representation. Calling gtoc to get the first */
/* several digits would be better, but it is too slow. */

	sprintf (gwdata->GWSTRING_REP, "A %ld-bit number", bits);

/* Return success */

	return (0);
}

/* This setup routine is for operations without a modulo. In essence, */
/* you are using gwnums as a general-purpose FFT multiply library. */

int gwsetup_without_mod (
	gwhandle *gwdata,	/* Placeholder for gwnum global data */
	unsigned long n)	/* Maximum number of bits in OUTPUT numbers. */
{
	struct gwasm_jmptab *info;
	unsigned long fftlen, max_exponent, desired_n;
	int	error_code;

/* Return delayed errors from gwinit2 */

	if (gwdata->GWERROR) return (gwdata->GWERROR);

/* Call gwinfo and have it figure out the FFT length we will use. */
/* Since the user must zero the upper half of FFT input data, the FFT */
/* outputs will be smaller.  This lets us get about another 0.3 bits */
/* per input word. */

	gwdata->safety_margin -= 0.3;
	error_code = gwinfo (gwdata, 1.0, 2, n, -1);
	gwdata->safety_margin += 0.3;
	info = gwdata->jmptab;
	if (error_code) return (error_code);

	max_exponent = info->max_exp;
	fftlen = info->fftlen;

/* If possible, increase n to the next multiple of FFT length.  This is */
/* because rational FFTs are faster than irrational FFTs (no FFT weights). */

	desired_n = ((n + fftlen - 1) / fftlen) * fftlen;
	if (desired_n < max_exponent) n = desired_n;

/* Our FFTs don't handle cases where there are few bits per word because */
/* carries must be propagated over too many words.  Arbitrarily insist */
/* that n is at least 12 * fftlen.  */

	if (n < 12 * fftlen) n = 12 * fftlen;

/* Now setup the assembly code */

	gwdata->safety_margin -= 0.3;
	error_code = internal_gwsetup (gwdata, 1.0, 2, n, -1);
	gwdata->safety_margin += 0.3;
	if (error_code) return (error_code);

/* Set flag indicating general-purpose modulo operations are not in force */

	gwdata->GENERAL_MOD = FALSE;

/* Create dummy string representation. */

	strcpy (gwdata->GWSTRING_REP, "No modulus");

/* Return success */

	return (0);
}


/* Examine a giant to see if it a (k*2^n+c)/d value. */
/* Returns TRUE if conversion was successful. */

int convert_giant_to_k2ncd (
	giant	g,		/* Giant to examine */
	double	*k,		/* K in (K*2^N+C)/D. */
	unsigned long *n,	/* N in (K*2^N+C)/D. */
	signed long *c,		/* C in (K*2^N+C)/D. */
	unsigned long *d)	/* D in (K*2^N+C)/D. */
{
	unsigned long less_nasty_d;
	int	i;
	uint32_t quick_test;
	giant	test_g, alloc_g;
	stackgiant(tmp,3);

/* Loop through a lot of small d values in hopes of finding a multiplier that */
/* will convert the input value into a k*2^n+c value.  We do this because */
/* small d values are the ones that generate repeating bit patterns in */
/* the modulus and unexpectedly large round off errors during operations */

	alloc_g = NULL;
	less_nasty_d = 1;
	for (*d = 1; *d <= 999; *d += 2) {

/* Do a quick test to see if this is a viable candidate */

		quick_test = (g->n[1] * *d) & 0xFFFFFC00;
		if (quick_test != 0 && quick_test != 0xFFFFFC00) continue;

/* Compute g * d to see if it has the proper k*2^n+c bit pattern */

		if (*d == 1) {
			test_g = g;
		} else {
			if (alloc_g == NULL) {
				alloc_g = allocgiant (((bitlen (g) + 10) >> 5) + 1);
				if (alloc_g == NULL) return (FALSE);
			}
			ultog (*d, alloc_g);
			mulg (g, alloc_g);
			test_g = alloc_g;
		}

/* See if this d value might result in a less nasty bit pattern for */
/* emulate_mod.  For example, Phi(82730,2) behaves much better if you */
/* multiply the modulus by 11. */

		if (*d >= 3 && less_nasty_d == 1) {
			int	count = 0;
			for (i = 0; i < test_g->sign; i++)
				if (test_g->n[i] == 0 || test_g->n[i] == -1) count++;
			if (count >= (test_g->sign >> 2)) less_nasty_d = *d;
		}

/* See if low order 2 words are viable for a k*2^n+c candidate */

		*c = (int32_t) test_g->n[0];

		if (test_g->n[1] == 0 && test_g->n[0] <= MAX_ZEROPAD_C);
		else if (test_g->n[1] == 0xFFFFFFFF && *c < 0 && *c >= -MAX_ZEROPAD_C);
		else continue;

/* Examine the middle words */
	
		for (i = 2; i < test_g->sign - 1; i++)
			if (test_g->n[i] != test_g->n[1]) break;

/* Now see if the high bits can fit in a 51-bit k */

		tmp->n[0] = test_g->n[i]; tmp->sign = 1;
		if (test_g->sign - i >= 2) { tmp->n[1] = test_g->n[i+1]; tmp->sign = 2; }
		if (test_g->sign - i >= 3) { tmp->n[2] = test_g->n[i+2]; tmp->sign = 3; }
		if (test_g->sign - i >= 4) continue;
		if (test_g->n[1] == 0xFFFFFFFF) iaddg (1, tmp);

		*n = i * 32;
		while ((tmp->n[0] & 0x1) == 0) {
			gshiftright(1, tmp);
			(*n)++;
		}
		if (bitlen (tmp) > 51) continue;

/* Set k and return success */

		*k = tmp->n[0];
		if (tmp->sign == 2) *k += (double) tmp->n[1] * 4294967296.0;
		free (alloc_g);
		return (TRUE);
	}

/* No luck in finding a (k*2^n+c)/d equivalent for the input value */

	*d = less_nasty_d;
	free (alloc_g);
	return (FALSE);
}

/* Common setup routine for the three different user-visible setup routines */
/* Allocate memory and initialize assembly code for arithmetic */
/* modulo k*b^n+c */

int internal_gwsetup (
	gwhandle *gwdata,	/* Placeholder for gwnum global data */
	double	k,		/* K in K*B^N+C. Must be a positive integer. */
	unsigned long b,	/* B in K*B^N+C. */
	unsigned long n,	/* N in K*B^N+C. Exponent to test. */
	signed long c)		/* C in K*B^N+C. Must be rel. prime to K. */
{
	struct gwasm_jmptab *info;
	void	*asm_data_alloc;
	struct gwasm_data *asm_data;
	int	error_code;
	unsigned long mem_needed;
	double	*tables;		/* Pointer tables we are building */
	unsigned long pass1_size, pass2_size;
	double	small_word, big_word, temp, asm_values[40];

/* Remember the arguments */

	gwdata->k = k;
	gwdata->b = b;
	gwdata->n = n;
	gwdata->c = c;

/* Init the FPU to assure we are in 64-bit precision mode */

	fpu_init ();

/* Select the proper FFT size for this k,b,n,c combination */

	error_code = gwinfo (gwdata, k, b, n, c);
	if (error_code) return (error_code);
	info = gwdata->jmptab;

/* Get pointer to fft info and allocate needed memory.  If we are */
/* trying to allocate large pages, then also allocate space for the */
/* one gwnum value that will be stored in large pages.  We allocate */
/* a little extra space to align the gwnum on a cache line boundary */
/* and because gwnum_size may not produce an accurate value prior */
/* to gwsetup completing. */

	tables = NULL;
	gwdata->large_pages_ptr = NULL;
	gwdata->large_pages_gwnum = NULL;
	mem_needed = info->mem_needed + info->scratch_size;
	if (gwdata->use_large_pages) {
		tables = (double *) large_pages_malloc (mem_needed + gwnum_size (gwdata) + 4096);
		if (tables != NULL) {
			/* Save large pages pointer for later freeing */
			gwdata->large_pages_ptr = tables;
			/* Save pointer to the gwnum we also allocated, so */
			/* that first gwalloc call can return it. */
			gwdata->large_pages_gwnum = align_ptr ((char *) tables + mem_needed, 128);
		}
	}
	if (tables == NULL) {
		tables = (double *) aligned_malloc (mem_needed, 4096);
		if (tables == NULL) return (GWERROR_MALLOC);
	}
	gwdata->gwnum_memory = tables;

/* Do a seemingly pointless memset! */
/* The memset will walk through the allocated memory sequentially, which */
/* increases the likelihood that contiguous virtual memory will map to */
/* contiguous physical memory. */

	memset (tables, 0, mem_needed);

/* Allocate space for the assembly code global data.  This area is preceded */
/* by a temporary stack.  This allows the assembly code to access the global */
/* data using offsets from the stack pointer. */

	asm_data_alloc = aligned_malloc (sizeof (struct gwasm_data) + NEW_STACK_SIZE, 4096);
	if (asm_data_alloc == NULL) {
		gwdone (gwdata);
		return (GWERROR_MALLOC);
	}
	gwdata->asm_data = (char *) asm_data_alloc + NEW_STACK_SIZE;
	asm_data = (struct gwasm_data *) gwdata->asm_data;

/* Copy values for asm code to use */

	asm_data->FFTLEN = gwdata->FFTLEN;
	asm_data->ZERO_PADDED_FFT = gwdata->ZERO_PADDED_FFT;
	asm_data->ALL_COMPLEX_FFT = gwdata->ALL_COMPLEX_FFT;
	asm_data->B_IS_2 = (b == 2);

/* Initialize the extended precision code that computes the FFT weights */

	gwdata->dd_data = gwdbldbl_data_alloc ();
	if (gwdata->dd_data == NULL) {
		gwdone (gwdata);
		return (GWERROR_MALLOC);
	}
	gwfft_weight_setup (gwdata->dd_data, gwdata->ZERO_PADDED_FFT,
			    k, b, n, c, gwdata->FFTLEN);

/* Calculate the number of bits in k*b^n.  This will be helpful in */
/* determining how much memory to allocate for giants. */

	gwdata->bit_length = log2 (k) + n * log2 (b);

/* Calculate the average number of base b's stored in each FFT word.  The total */
/* number of base b's the underlying FFT works with (i.e. the point at which data */
/* wraps around to the low FFT word) is 2*n for a zero pad FFT and logb(k) + n */
/* otherwise. */	

	gwdata->avg_num_b_per_word =
		(gwdata->ZERO_PADDED_FFT ? n * 2.0 : (logb (k) + n)) / gwdata->FFTLEN;

/* Calculate the number of base b's stored in each small FFT word. */

	gwdata->NUM_B_PER_SMALL_WORD = (unsigned long) gwdata->avg_num_b_per_word;

/* Set a flag if this is a rational FFT.  That is, an FFT where all the */
/* weighting factors are 1.0.  This happens when c is -1 and the every */
/* FFT word is has the same number of b's.  The assembly code can make some */
/* obvious optimizations when all the FFT weights are one. */

	gwdata->RATIONAL_FFT = asm_data->RATIONAL_FFT =
		((double) gwdata->NUM_B_PER_SMALL_WORD == gwdata->avg_num_b_per_word) && (c == -1);

/* Remember the maximum number of bits per word that this FFT length */
/* supports.  We this in gwnear_fft_limit.  Note that zero padded FFTs */
/* can support an extra 0.3 bits per word because of the all the zeroes. */

	gwdata->fft_max_bits_per_word = (double) info->max_exp / (double) gwdata->FFTLEN;
	if (gwdata->ZERO_PADDED_FFT) gwdata->fft_max_bits_per_word += 0.3;

/* Compute extra bits - the number of adds we can tolerate without */
/* a normalization operation. Under normal circumstances, max_bits */
/* will be greater than virtual bits, but playing with the safety margin */
/* or forcing use of a specific FFT length could change that. */

	gwdata->EXTRA_BITS = (unsigned long)
		pow (2.0, (gwdata->fft_max_bits_per_word - virtual_bits_per_word (gwdata)) / 2.0);

/* Remember the size of the scratch area */

	gwdata->SCRATCH_SIZE = info->scratch_size;

/* See how many cache lines are grouped in pass 1.  This will affect how */
/* we build the normalization tables.  Note that cache line sizes are */
/* different in the x87 (16 bytes) and SSE2 code (64 bytes). */

	gwdata->PASS1_CACHE_LINES = (info->flags_min_l2_cache_clm & 0xFFFF);

/* Determine the pass 1 & pass 2 sizes.  This affects how we build */
/* many of the sin/cos tables. */

	pass2_size = 1 << gwdata->PASS2_LEVELS;	/* Complex values in pass2 section */
	pass1_size = gwdata->FFTLEN / pass2_size; /* Real values in a pass1 section */

/* Initialize tables for the SSE2 assembly code. */

	if (gwdata->cpu_flags & CPU_SSE2) {

/* Build sin/cos and premultiplier tables used in pass 2 of two pass FFTs */
/* Remember that pass2_size is the number of complex values in a pass 2 */
/* section, but build_sin_cos_table wants the number of reals in a section. */
/* However, we build a 1/4-sized table by mixing some of the sin/cos */
/* data into the premultiplier table.  So, divide pass2_size by 2 instead of */
/* multiplying pass2_size by 2. */

/* For best prefetching, make sure tables remain on 128-byte boundaries */

		if (pass2_size > 1) {
			ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
			asm_data->pass2_premults = tables;
			tables = build_premult_table (gwdata, tables, pass2_size);

/* Build the rest of the tables */

			ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
			asm_data->xsincos_complex = tables;
			tables = build_sin_cos_table (tables, pass2_size/2, 0, 1);

			if (!gwdata->ALL_COMPLEX_FFT) {
				ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
				asm_data->sincos6 = tables;
				tables = build_sin_cos_table (tables, pass2_size * 4, 1, 2);
				asm_data->sincos7 = tables;
				tables = build_sin_cos_table (tables, pass2_size, 1, 1);
				tables += (16 - (tables - gwdata->gwnum_memory)) & 15;
			}

//			if (pass2_size == pow_two_above_or_equal (pass2_size)) {
				asm_data->sincos8 = asm_data->sincos7;
				asm_data->sincos9 = asm_data->sincos7;
				asm_data->sincos10 = asm_data->sincos7;
				asm_data->sincos11 = asm_data->sincos7;
				asm_data->sincos12 = asm_data->sincos7;
//			} else {
//				asm_data->sincos8 = tables;
//				tables = build_sin_cos_table (tables, pass2_size/4, !gwdata->ALL_COMPLEX_FFT, 1);
//				asm_data->sincos9 = tables;
//				tables = build_sin_cos_table (tables, pass2_size/16, !gwdata->ALL_COMPLEX_FFT, 1);
//				asm_data->sincos10 = tables;
//				tables = build_sin_cos_table (tables, pass2_size/64, !gwdata->ALL_COMPLEX_FFT, 1);
//				asm_data->sincos11 = tables;
//				tables = build_sin_cos_table (tables, pass2_size/256, !gwdata->ALL_COMPLEX_FFT, 1);
//				asm_data->sincos12 = tables;
//				tables = build_sin_cos_table (tables, pass2_size/1024, !gwdata->ALL_COMPLEX_FFT, 1);
//				tables += (16 - (tables - gwdata->gwnum_memory)) & 15;
//			}
		}

/* Allocate a table for carries.  Init with XMM_BIGVAL.  For best */
/* distribution of data in the L2 cache, make this table contiguous */
/* with all other data used in the first pass (scratch area, normalization */
/* tables, etc.)  Note that we put the tables that are only partly loaded */
/* (column multipliers and big/lit table after the tables that are */
/* loaded throughtout the first pass. */

		if (pass2_size > 1) {
			int	i, carry_table_size;
			double	xmm_bigval;
			ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
			asm_data->carries = tables;
			carry_table_size = gwdata->FFTLEN / (pass2_size / 2);
			xmm_bigval = 3.0 * 131072.0 * 131072.0 * 131072.0;
			for (i = 0; i < carry_table_size; i++)
				*tables++ = xmm_bigval;
			tables += (16 - (tables - gwdata->gwnum_memory)) & 15;
		}

/* Build the group muliplier normalization table.  Keep this table */
/* contiguous with other data used in pass 1. */

		ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
		asm_data->norm_grp_mults = tables;
		tables = build_norm_table (gwdata, tables, pass2_size, 0);

/* Build the plus1-pre-multiplier table (complex weights applied when c > 0 */
/* and we are doing a all-complex FFT rather than emulating it with a */
/* zero-padded FFT. */

		if (gwdata->ALL_COMPLEX_FFT) {
			ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
			asm_data->plus1_premults = tables;
			tables = build_plus1_table (gwdata, tables, pass2_size);
		}

/* Reserve room for the pass 1 scratch area. */

		if (gwdata->SCRATCH_SIZE) {
			ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
			asm_data->scratch_area = tables;
			tables = (double *) ((char *) tables + gwdata->SCRATCH_SIZE);
		}

/* Build sin/cos tables used in pass 1.  If FFTLEN is a power of two, */
/* many of the sin/cos tables can be shared. */

		ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
		asm_data->sincos1 = tables;
		tables = build_sin_cos_table (tables, pass1_size, !gwdata->ALL_COMPLEX_FFT, pass2_size == 1 ? 2 : 1);

		if (pass2_size > 1 && pass1_size == pow_two_above_or_equal (pass1_size))
			asm_data->sincos2 = asm_data->sincos1;
		else {
			asm_data->sincos2 = tables;
			tables = build_sin_cos_table (tables, pass1_size/4, !gwdata->ALL_COMPLEX_FFT, 1);
		}

		if (pass1_size == pow_two_above_or_equal (pass1_size)) {
			asm_data->sincos3 = asm_data->sincos2;
			asm_data->sincos4 = asm_data->sincos2;
			asm_data->sincos5 = asm_data->sincos2;
		} else {
			asm_data->sincos3 = tables;
			tables = build_sin_cos_table (tables, pass1_size/16, !gwdata->ALL_COMPLEX_FFT, 1);
			asm_data->sincos4 = tables;
			tables = build_sin_cos_table (tables, pass1_size/64, !gwdata->ALL_COMPLEX_FFT, 1);
			asm_data->sincos5 = tables;
			tables = build_sin_cos_table (tables, pass1_size/256, !gwdata->ALL_COMPLEX_FFT, 1);
		}
		tables += (16 - (tables - gwdata->gwnum_memory)) & 15;

/* Build the table of big vs. little flags.  This cannot be last table */
/* built as xnorm_2d macro reads 2 bytes past the end of the array. */

		ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
		asm_data->norm_biglit_array = tables;
		tables = build_biglit_table (gwdata, tables, pass2_size);
		tables += (16 - (tables - gwdata->gwnum_memory)) & 15;

/* Build the column normalization multiplier table. */

		ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
		asm_data->norm_col_mults = tables;
		tables = build_norm_table (gwdata, tables, pass2_size, 1);
	}

/* Initialze table for the x87 assembly code. */

	if (! (gwdata->cpu_flags & CPU_SSE2)) {

/* Allocate a table for carries.  Init with zero.  For best */
/* distribution of data in the L2 cache, make this table contiguous */
/* with the scratch area which is also used in the first pass. */

		if (pass2_size > 1) {
			int	i, carry_table_size;
			asm_data->carries = tables;
			carry_table_size = gwdata->FFTLEN / pass2_size;
			for (i = 0; i < carry_table_size; i++) *tables++ = 0.0;
		}

/* Reserve room for the pass 1 scratch area. */

		asm_data->scratch_area = tables;
		if (gwdata->SCRATCH_SIZE)
			tables = (double *) ((char *) tables + gwdata->SCRATCH_SIZE);

/* Build the group muliplier normalization table.  Keep this table */
/* contiguous with other data used in pass 1. */

		asm_data->norm_grp_mults = tables;
		tables = build_x87_norm_table (gwdata, tables, pass2_size, 0);

/* Build sin/cos tables used in pass 1.  If FFTLEN is a power of two, */
/* many of the sin/cos tables can be shared. */

		asm_data->sincos1 = tables;
		tables = build_x87_sin_cos_table (tables, pass1_size, !gwdata->ALL_COMPLEX_FFT);

		if (pass1_size == pow_two_above_or_equal (pass1_size)) {
			asm_data->sincos2 = asm_data->sincos1;
			asm_data->sincos3 = asm_data->sincos1;
			asm_data->sincos4 = asm_data->sincos1;
			asm_data->sincos5 = asm_data->sincos1;
		} else {
			asm_data->sincos2 = tables;
			tables = build_x87_sin_cos_table (tables, pass1_size/4, !gwdata->ALL_COMPLEX_FFT);
			asm_data->sincos3 = tables;
			tables = build_x87_sin_cos_table (tables, pass1_size/16, !gwdata->ALL_COMPLEX_FFT);
			asm_data->sincos4 = tables;
			tables = build_x87_sin_cos_table (tables, pass1_size/64, !gwdata->ALL_COMPLEX_FFT);
			asm_data->sincos5 = tables;
			tables = build_x87_sin_cos_table (tables, pass1_size/256, !gwdata->ALL_COMPLEX_FFT);
		}

/* Build sin/cos and premultiplier tables used in pass 2 of two pass FFTs */
/* Remember that pass2_size is the number of complex values in a pass 2 */
/* section, but build_x87_sin_cos_table wants the number of reals in */
/* a section. */

		if (pass2_size > 1) {
			asm_data->pass2_premults = tables;
			tables = build_x87_premult_table (gwdata, tables, pass2_size);
			asm_data->xsincos_complex = tables;
			tables = build_x87_sin_cos_table (tables, pass2_size*2, 0);

			if (!gwdata->ALL_COMPLEX_FFT) {
				asm_data->sincos6 = tables;
				tables = build_x87_sin_cos_table (tables, pass2_size*2, 1);
				asm_data->sincos7 = asm_data->sincos6;
				asm_data->sincos8 = asm_data->sincos6;
				asm_data->sincos9 = asm_data->sincos6;
				asm_data->sincos10 = asm_data->sincos6;
				asm_data->sincos11 = asm_data->sincos6;
				asm_data->sincos12 = asm_data->sincos6;
			}
		}

/* Build the plus1-pre-multiplier table (complex weights applied when c > 0 */
/* and we are doing a all-complex FFT rather than emulating it with a */
/* zero-padded FFT. */

		if (gwdata->ALL_COMPLEX_FFT) {
			asm_data->plus1_premults = tables;
			tables = build_x87_plus1_table (gwdata, tables, pass2_size);
		}

/* Build the column normalization multiplier table. */

		asm_data->norm_col_mults = tables;
		tables = build_x87_norm_table (gwdata, tables, pass2_size, 1);

/* Build the table of big vs. little flags. */

		asm_data->norm_biglit_array = tables;
		tables = build_x87_biglit_table (gwdata, tables, pass2_size);
	}

/* Finish verifying table size */

#ifdef GDEBUG
	{char buf[80];
	intptr_t mem = (intptr_t) tables - (intptr_t) gwdata->gwnum_memory;
	if (mem != mem_needed) {
		sprintf (buf, "FFTlen: %d, mem needed should be: %d\n",
			 gwdata->FFTLEN, (int) (mem - info->scratch_size));
		OutputBoth(0,buf);}}
#endif

/* Do more assembly initialization */

	asm_data->cache_line_multiplier = info->flags_min_l2_cache_clm & 0xFFFF;
	asm_data->addcount1 = info->counts[0];
	asm_data->normcount1 = info->counts[1];
	asm_data->count1 = info->counts[2];
	asm_data->count2 = info->counts[3];
	asm_data->count3 = info->counts[4];
	asm_data->count4 = info->counts[5];
	asm_data->count5 = info->counts[6];
	asm_data->zero_fft = 0;
	asm_data->const_fft = 0;
	asm_data->ZPAD0_7[7] = 0.0;
	asm_data->zpad_addr = &asm_data->ZPAD0_7[0];

/* Init constants.  I could have used true global variables but I did not */
/* so that these constants could be close to other variables used at the */
/* same time.  This might free up a cache line or two. */

	asm_data->XMM_TWO[0] = asm_data->XMM_TWO[1] = 2.0;
 	asm_data->HALF = 0.5;
 	asm_data->XMM_HALF[0] = asm_data->XMM_HALF[1] = 0.5;
	asm_data->SQRTHALF =
	asm_data->XMM_SQRTHALF[0] = asm_data->XMM_SQRTHALF[1] = sqrt (0.5);
	asm_data->P25 = 0.25;
	asm_data->P75 = 0.75;
	asm_data->P3 = 3.0;
	asm_data->XMM_ABSVAL[0] = asm_data->XMM_ABSVAL[2] = 0xFFFFFFFF;
	asm_data->XMM_ABSVAL[1] = asm_data->XMM_ABSVAL[3] = 0x7FFFFFFF;

	asm_data->XMM_TTP_FUDGE[0] = asm_data->XMM_TTP_FUDGE[1] =
	asm_data->XMM_TTP_FUDGE[2] = asm_data->XMM_TTP_FUDGE[3] =
	asm_data->XMM_TTP_FUDGE[4] = asm_data->XMM_TTP_FUDGE[5] =
	asm_data->XMM_TTP_FUDGE[6] = asm_data->XMM_TTP_FUDGE[7] =
	asm_data->XMM_TTP_FUDGE[9] = asm_data->XMM_TTP_FUDGE[11] =
	asm_data->XMM_TTP_FUDGE[13] = asm_data->XMM_TTP_FUDGE[15] =
	asm_data->XMM_TTP_FUDGE[16] = asm_data->XMM_TTP_FUDGE[18] =
	asm_data->XMM_TTP_FUDGE[20] = asm_data->XMM_TTP_FUDGE[22] = 1.0;
	asm_data->XMM_TTP_FUDGE[8] = asm_data->XMM_TTP_FUDGE[10] =
	asm_data->XMM_TTP_FUDGE[12] = asm_data->XMM_TTP_FUDGE[14] =
	asm_data->XMM_TTP_FUDGE[17] = asm_data->XMM_TTP_FUDGE[19] =
	asm_data->XMM_TTP_FUDGE[21] = asm_data->XMM_TTP_FUDGE[23] =
	asm_data->XMM_TTP_FUDGE[24] = asm_data->XMM_TTP_FUDGE[25] =
	asm_data->XMM_TTP_FUDGE[26] = asm_data->XMM_TTP_FUDGE[27] =
	asm_data->XMM_TTP_FUDGE[28] = asm_data->XMM_TTP_FUDGE[29] =
	asm_data->XMM_TTP_FUDGE[30] = asm_data->XMM_TTP_FUDGE[31] = 1.0 / (double) b;

	asm_data->XMM_TTMP_FUDGE[0] = asm_data->XMM_TTMP_FUDGE[1] =
	asm_data->XMM_TTMP_FUDGE[2] = asm_data->XMM_TTMP_FUDGE[3] =
	asm_data->XMM_TTMP_FUDGE[4] = asm_data->XMM_TTMP_FUDGE[5] =
	asm_data->XMM_TTMP_FUDGE[6] = asm_data->XMM_TTMP_FUDGE[7] =
	asm_data->XMM_TTMP_FUDGE[9] = asm_data->XMM_TTMP_FUDGE[11] =
	asm_data->XMM_TTMP_FUDGE[13] = asm_data->XMM_TTMP_FUDGE[15] =
	asm_data->XMM_TTMP_FUDGE[16] = asm_data->XMM_TTMP_FUDGE[18] =
	asm_data->XMM_TTMP_FUDGE[20] = asm_data->XMM_TTMP_FUDGE[22] = 1.0;
	asm_data->XMM_TTMP_FUDGE[8] = asm_data->XMM_TTMP_FUDGE[10] =
	asm_data->XMM_TTMP_FUDGE[12] = asm_data->XMM_TTMP_FUDGE[14] =
	asm_data->XMM_TTMP_FUDGE[17] = asm_data->XMM_TTMP_FUDGE[19] =
	asm_data->XMM_TTMP_FUDGE[21] = asm_data->XMM_TTMP_FUDGE[23] =
	asm_data->XMM_TTMP_FUDGE[24] = asm_data->XMM_TTMP_FUDGE[25] =
	asm_data->XMM_TTMP_FUDGE[26] = asm_data->XMM_TTMP_FUDGE[27] =
	asm_data->XMM_TTMP_FUDGE[28] = asm_data->XMM_TTMP_FUDGE[29] =
	asm_data->XMM_TTMP_FUDGE[30] = asm_data->XMM_TTMP_FUDGE[31] = (double) b;

/* Compute the x87 (64-bit) rounding constants */

	small_word = pow ((double) b, gwdata->NUM_B_PER_SMALL_WORD);
	big_word = (double) b * small_word;
	asm_data->BIGVAL = (float) (3.0 * pow (2.0, 62.0));
	asm_data->BIGBIGVAL = (float) (big_word * asm_data->BIGVAL);

/* Compute the SSE2 (53-bit) rounding constants */

	asm_data->XMM_BIGVAL[0] =
	asm_data->XMM_BIGVAL[1] = 3.0 * pow (2.0, 51.0);
	asm_data->XMM_BIGVAL_NEG[0] =
	asm_data->XMM_BIGVAL_NEG[1] = -asm_data->XMM_BIGVAL[0];
	asm_data->XMM_BIGBIGVAL[0] =
	asm_data->XMM_BIGBIGVAL[1] = big_word * asm_data->XMM_BIGVAL[0];

/* Negate c and store as a double */

	asm_data->MINUS_C =
	asm_data->XMM_MINUS_C[0] =
	asm_data->XMM_MINUS_C[1] = (double) -c;

/* The sumout value is FFTLEN/2 times larger than it should be.  Create an */
/* inverse to properly calculate the sumout when a multiplication ends. */

	asm_data->ttmp_ff_inv = 2.0 / (double) gwdata->FFTLEN;

/* Compute constant that converts fft_weight_over_fftlen found in the */
/* two-to-minus-phi tables into the true fft_weight value.  This is usually */
/* FFTLEN / 2, but when doing a non-zero-padded FFT this is FFTLEN / 2k. */

	asm_data->XMM_NORM012_FF[0] =
	asm_data->XMM_NORM012_FF[1] =
		(gwdata->ZERO_PADDED_FFT) ?
			(double) (gwdata->FFTLEN / 2) :
			(double) (gwdata->FFTLEN / 2) / k;

/* Split k for zero-padded FFTs emulating modulo k*b^n+c.  Note we only need */
/* to split k if k * big_word exceeds what a floating point register can hold. */
/* 2^49 should give us enough room to handle a carry and still round properly. */	

	if (gwdata->cpu_flags & CPU_SSE2 && k * big_word < 562949953421312.0) {
		asm_data->XMM_K_HI[0] = asm_data->XMM_K_HI[1] = 0.0;
		asm_data->XMM_K_LO[0] = asm_data->XMM_K_LO[1] = k;
	} else {
		asm_data->XMM_K_HI[0] =
		asm_data->XMM_K_HI[1] = floor (k / big_word) * big_word;
		asm_data->XMM_K_LO[0] =
		asm_data->XMM_K_LO[1] = k - asm_data->XMM_K_HI[0];
	}
	asm_data->XMM_K_HI_2[0] =
	asm_data->XMM_K_HI_2[1] = floor (k / big_word / big_word) * big_word * big_word;
	asm_data->XMM_K_HI_1[0] =
	asm_data->XMM_K_HI_1[1] = asm_data->XMM_K_HI[0] - asm_data->XMM_K_HI_2[0];
	gwsetmulbyconst (gwdata, gwdata->maxmulbyconst);

/* Compute the normalization constants indexed by biglit array entries */

	temp = 1.0 / small_word;	/* Compute lower limit inverse */
	asm_data->XMM_LIMIT_INVERSE[0] =
	asm_data->XMM_LIMIT_INVERSE[1] =
	asm_data->XMM_LIMIT_INVERSE[3] =
	asm_data->XMM_LIMIT_INVERSE[4] = temp;

					/* Compute lower limit bigmax */
	if (b > 2)
		temp = small_word;
	else if (gwdata->cpu_flags & CPU_SSE2)
		temp = small_word * asm_data->XMM_BIGVAL[0] - asm_data->XMM_BIGVAL[0];
	else
		temp = small_word * asm_data->BIGVAL - asm_data->BIGVAL;
	asm_data->XMM_LIMIT_BIGMAX[0] =
	asm_data->XMM_LIMIT_BIGMAX[1] =
	asm_data->XMM_LIMIT_BIGMAX[3] =
	asm_data->XMM_LIMIT_BIGMAX[4] = temp;

	temp = -temp;			/* Negative lower limit bigmax */
	asm_data->XMM_LIMIT_BIGMAX_NEG[0] =
	asm_data->XMM_LIMIT_BIGMAX_NEG[1] =
	asm_data->XMM_LIMIT_BIGMAX_NEG[3] =
	asm_data->XMM_LIMIT_BIGMAX_NEG[4] = temp;

	temp = 1.0 / big_word;		/* Compute upper limit inverse */
	asm_data->XMM_LIMIT_INVERSE[2] =
	asm_data->XMM_LIMIT_INVERSE[5] =
	asm_data->XMM_LIMIT_INVERSE[6] =
	asm_data->XMM_LIMIT_INVERSE[7] = temp;

					/* Compute upper limit bigmax */
	if (b > 2)
		temp = big_word;
	else if (gwdata->cpu_flags & CPU_SSE2)
		temp = big_word * asm_data->XMM_BIGVAL[0] - asm_data->XMM_BIGVAL[0];
	else
		temp = big_word * asm_data->BIGVAL - asm_data->BIGVAL;
	asm_data->XMM_LIMIT_BIGMAX[2] =
	asm_data->XMM_LIMIT_BIGMAX[5] =
	asm_data->XMM_LIMIT_BIGMAX[6] =
	asm_data->XMM_LIMIT_BIGMAX[7] = temp;

	temp = -temp;			/* Negative upper limit bigmax */
	asm_data->XMM_LIMIT_BIGMAX_NEG[2] =
	asm_data->XMM_LIMIT_BIGMAX_NEG[5] =
	asm_data->XMM_LIMIT_BIGMAX_NEG[6] =
	asm_data->XMM_LIMIT_BIGMAX_NEG[7] = temp;

	memcpy (asm_data->XMM_LIMIT_INVERSE+8, asm_data->XMM_LIMIT_INVERSE, 8 * sizeof (double));
	memcpy (asm_data->XMM_LIMIT_INVERSE+16, asm_data->XMM_LIMIT_INVERSE, 8 * sizeof (double));
	memcpy (asm_data->XMM_LIMIT_INVERSE+24, asm_data->XMM_LIMIT_INVERSE, 8 * sizeof (double));
	memcpy (asm_data->XMM_LIMIT_BIGMAX+8, asm_data->XMM_LIMIT_BIGMAX, 8 * sizeof (double));
	memcpy (asm_data->XMM_LIMIT_BIGMAX+16, asm_data->XMM_LIMIT_BIGMAX, 8 * sizeof (double));
	memcpy (asm_data->XMM_LIMIT_BIGMAX+24, asm_data->XMM_LIMIT_BIGMAX, 8 * sizeof (double));
	memcpy (asm_data->XMM_LIMIT_BIGMAX_NEG+8, asm_data->XMM_LIMIT_BIGMAX_NEG, 8 * sizeof (double));
	memcpy (asm_data->XMM_LIMIT_BIGMAX_NEG+16, asm_data->XMM_LIMIT_BIGMAX_NEG, 8 * sizeof (double));
	memcpy (asm_data->XMM_LIMIT_BIGMAX_NEG+24, asm_data->XMM_LIMIT_BIGMAX_NEG, 8 * sizeof (double));

/* Now compute a number of constants the assembly code needs.  These will */
/* be copied to properly aligned (SSE2 constants must be on 16-byte */
/* boundaries) and grouped (for better cache line locality) assembly */
/* global variables.  We compute these constants in C code because it is */
/* easier and because early Pentiums do not support SSE2 instructions and */
/* x86-64 does not support x87 instructions. */

	gwasm_constants ((double *) &asm_values);

	asm_data->P951 =
	asm_data->XMM_P951[0] = asm_data->XMM_P951[1] = asm_values[0];
	asm_data->P618 =
	asm_data->XMM_P618[0] = asm_data->XMM_P618[1] = asm_values[1];
	asm_data->P309 =
	asm_data->XMM_P309[0] = asm_data->XMM_P309[1] = asm_values[2];
	asm_data->M262 =
	asm_data->XMM_M262[0] = asm_data->XMM_M262[1] = asm_values[3];
	asm_data->P588 =
	asm_data->XMM_P588[0] = asm_data->XMM_P588[1] = asm_values[4];
	asm_data->M162 =
	asm_data->XMM_M162[0] = asm_data->XMM_M162[1] = asm_values[5];
	asm_data->M809 =
	asm_data->XMM_M809[0] = asm_data->XMM_M809[1] = asm_values[6];
	asm_data->M382 =
	asm_data->XMM_M382[0] = asm_data->XMM_M382[1] = asm_values[7];
	asm_data->P866 =
	asm_data->XMM_P866[0] = asm_data->XMM_P866[1] = asm_values[8];
	asm_data->P433 = asm_values[9];
	asm_data->P577 = asm_values[10];
	asm_data->P975 =
	asm_data->XMM_P975[0] = asm_data->XMM_P975[1] = asm_values[11];
	asm_data->P445 =
	asm_data->XMM_P445[0] = asm_data->XMM_P445[1] = asm_values[12];
	asm_data->P180 =
	asm_data->XMM_P180[0] = asm_data->XMM_P180[1] = asm_values[13];
	asm_data->P623 =
	asm_data->XMM_P623[0] = asm_data->XMM_P623[1] = asm_values[14];
	asm_data->M358 =
	asm_data->XMM_M358[0] = asm_data->XMM_M358[1] = asm_values[15];
	asm_data->P404 =
	asm_data->XMM_P404[0] = asm_data->XMM_P404[1] = asm_values[16];
	asm_data->M223 = asm_values[17];
	asm_data->M901 = asm_values[18];
	asm_data->M691 = asm_values[19];

/* Non-SSE2 initialization. Calculate constants used in two pass FFTs. */
/* Foremost is the pass 1 blkdst and normalize blkdst for auxillary mult */
/* routines.  The two values are the same except for larger FFTs which */
/* use a scratch area. */

	if (! (gwdata->cpu_flags & CPU_SSE2)) {
		if (info->pass2_levels) {
			/* Pass 2 data size: 2^pass2_levels complex numbers */
			/* Compute number of 4KB pages, used in normalized */
			/* add/sub */
			asm_data->normval4 = (16 << info->pass2_levels) >> 12;

			/* pad 64 bytes every 4KB + 64 bytes between blocks */
			asm_data->pass1blkdst =
			asm_data->normblkdst =
				asm_data->normval4 * (4096+64) + 64;
			asm_data->normblkdst8 = 0;

			if (info->scratch_size) {
				/* Compute scratch area normblkdst */
				asm_data->normblkdst =
					asm_data->cache_line_multiplier * 32;

				/* Only larger pass1's have padding */
				if (info->counts[0] >= 128)
					asm_data->normblkdst8 = 64;
			}

			/* Compute "cache lines in a page" / clm */
			asm_data->normval1 = (4096/32) / asm_data->cache_line_multiplier;

			/* Big/lit flags ptr fudge factor in add/sub */
			asm_data->normval2 =
				asm_data->cache_line_multiplier * 2 *
				(asm_data->addcount1 - 1);

			/* Big/lit flags ptr fudge factor #2 in add/sub */
			asm_data->normval3 =
				asm_data->cache_line_multiplier * 2 -
				((asm_data->cache_line_multiplier * 2 +
				  asm_data->normval2) *
				 asm_data->normval1 * asm_data->normval4);
		}
	}

/* SSE2 initialization code formerly done in assembly language */

	else {

		/* Data size = 2^pass2_levels complex numbers * 2 (for SSE2) */
		/* Calculate number of 8KB pages */
		asm_data->normval4 = (32 << info->pass2_levels) >> 13;

		/* Gap between blocks */
		asm_data->pass2gapsize = info->counts[1] - 1;

		/* Pass 1 block distance includes 128 pad bytes every 8KB */
		/* plus the gap between blocks */ 
		asm_data->pass1blkdst =
				asm_data->normval4 * (8192+128) +
				asm_data->pass2gapsize;

		if (info->scratch_size == 0) {
			/* Normblkdst = pass1blkdst - clm*64 */
			asm_data->normblkdst =
				asm_data->pass1blkdst -
				asm_data->cache_line_multiplier * 64;
			asm_data->normblkdst8 = 0;
		} else {
			/* Pad in clmblkdst is zero, clmblkdst8 is 128 */
   			asm_data->normblkdst = 0;
			asm_data->normblkdst8 = 128;
		}

		if (info->pass2_levels) {
			/* Calculate "cache lines in a chunk" / clm */
			asm_data->normval1 = (8192/64) / asm_data->cache_line_multiplier;

			/* Big/lit flags ptr fudge factor in add/sub */
			asm_data->normval2 =
				asm_data->cache_line_multiplier * 4 *
				(asm_data->addcount1 - 1);

			/* Big/lit flags ptr fudge factor #2 in add/sub */
			asm_data->normval3 =
				asm_data->cache_line_multiplier * 4 -
				(asm_data->cache_line_multiplier * 4 +
						asm_data->normval2) *
				    asm_data->normval1 * asm_data->normval4;
		}
	}

/* If the carry must be spread over more than 2 words, then set global */
/* so that assembly code knows this.  We check to see if three small FFT words */
/* can absorb the expected number of bits in a result word.  We are not */
/* aggressive in pushing this limit (we assume no big words will absorb */
/* any of the carry) as it is not a major performance penalty to do 4 word */
/* carries.  In fact, we might should do 4 words all the time. */

	if (gwdata->ZERO_PADDED_FFT ||
	    3.0 * gwdata->NUM_B_PER_SMALL_WORD * log2 (b) >
			2.0 * ((gwdata->NUM_B_PER_SMALL_WORD + 1) * log2 (b) - 1) +
			0.6 * log2 (gwdata->FFTLEN) + log2 (k) + 1.7 * log2 (c))
		asm_data->SPREAD_CARRY_OVER_4_WORDS = FALSE;
	else
		asm_data->SPREAD_CARRY_OVER_4_WORDS = TRUE;

/* Set some global variables that make life easier in the assembly code */
/* that wraps carry out of top FFT word into the bottom FFT word. */
/* This is needed when k > 1 and we are not doing a zero padded FFT. */

	asm_data->TOP_CARRY_NEEDS_ADJUSTING = (k > 1.0 && !gwdata->ZERO_PADDED_FFT);
	if (asm_data->TOP_CARRY_NEEDS_ADJUSTING) {
		unsigned long num_b_in_k, kbits, kbits_lo;
		unsigned long num_b_in_top_word, num_b_in_second_top_word, num_b_in_third_top_word;

/* Invert k and split k for computing top carry adjustment without */
/* precision problems. */

		asm_data->INVERSE_K = 1.0 / k;
		kbits = (unsigned long) ceil (log2 (k));
		kbits_lo = kbits / 2;
 		asm_data->K_HI = ((unsigned long) k) & ~((1 << kbits_lo) - 1);
		asm_data->K_LO = ((unsigned long) k) &  ((1 << kbits_lo) - 1);

/* Calculate top carry adjusting constants */

		num_b_in_top_word = gwdata->NUM_B_PER_SMALL_WORD;
		if (is_big_word (gwdata, gwdata->FFTLEN-1)) num_b_in_top_word++;
		num_b_in_second_top_word = gwdata->NUM_B_PER_SMALL_WORD;
		if (is_big_word (gwdata, gwdata->FFTLEN-2)) num_b_in_second_top_word++;
		num_b_in_third_top_word = gwdata->NUM_B_PER_SMALL_WORD;
		if (is_big_word (gwdata, gwdata->FFTLEN-3)) num_b_in_third_top_word++;

		num_b_in_k = (unsigned long) ceil (logb (k));
		asm_data->CARRY_ADJUST1 = pow ((double) b, num_b_in_k);
		asm_data->CARRY_ADJUST1_HI = ((unsigned long) asm_data->CARRY_ADJUST1) & ~((1 << kbits_lo) - 1);
		asm_data->CARRY_ADJUST1_LO = ((unsigned long) asm_data->CARRY_ADJUST1) &  ((1 << kbits_lo) - 1);
		asm_data->CARRY_ADJUST2 = pow ((double) b, num_b_in_top_word) / asm_data->CARRY_ADJUST1;
		asm_data->CARRY_ADJUST3 = gwfft_weight (gwdata->dd_data, gwdata->FFTLEN-1);
		asm_data->CARRY_ADJUST4 = pow ((double) b, num_b_in_second_top_word);
		asm_data->CARRY_ADJUST5 = gwfft_weight (gwdata->dd_data, gwdata->FFTLEN-2);
		asm_data->CARRY_ADJUST6 = pow ((double) b, num_b_in_third_top_word);
		asm_data->CARRY_ADJUST7 = gwfft_weight (gwdata->dd_data, gwdata->FFTLEN-3);
		// It's too hard to upgrade the old x87 code to match the SSE2 code
		// So, for x87 generate the same constants used prior to version 25.11
		if (! (gwdata->cpu_flags & CPU_SSE2)) {
			if (gwdata->PASS2_LEVELS)
				asm_data->CARRY_ADJUST4 *= asm_data->CARRY_ADJUST5;
			else
				asm_data->CARRY_ADJUST6 *= asm_data->CARRY_ADJUST7;
		}

/* In two-pass FFTs, we only support tweaking the top two words. In one-pass FFTs, */
/* we adjust the top three words.  Make sure this works.  A test case that fails: */
/* 489539*3^72778+1.  We should consider supporting tweaking the top 3 words. */

		ASSERTG ((gwdata->PASS2_LEVELS &&
			  num_b_in_k <= num_b_in_top_word + num_b_in_second_top_word) ||
			 (!gwdata->PASS2_LEVELS &&
			  num_b_in_k <= num_b_in_top_word + num_b_in_second_top_word + num_b_in_third_top_word));

/* Get the addr of the top three words.  This is funky because in two-pass */
/* FFTs we want the scratch area offset when normalizing after a multiply, */
/* but the FFT data when normalizing after an add/sub.  For one-pass FFTs, */
/* we always want the FFT data offset. */

		asm_data->HIGH_WORD1_OFFSET = addr_offset (gwdata, gwdata->FFTLEN-1);
		asm_data->HIGH_WORD2_OFFSET = addr_offset (gwdata, gwdata->FFTLEN-2);
		asm_data->HIGH_WORD3_OFFSET = addr_offset (gwdata, gwdata->FFTLEN-3);

		raw_gwsetaddin (gwdata, gwdata->FFTLEN-1, 0.0);
		asm_data->HIGH_SCRATCH1_OFFSET = asm_data->ADDIN_OFFSET;
		raw_gwsetaddin (gwdata, gwdata->FFTLEN-2, 0.0);
		asm_data->HIGH_SCRATCH2_OFFSET = asm_data->ADDIN_OFFSET;
		raw_gwsetaddin (gwdata, gwdata->FFTLEN-3, 0.0);
		asm_data->HIGH_SCRATCH3_OFFSET = asm_data->ADDIN_OFFSET;
	}

/* Set some global variables that make life easier in the assembly code */
/* that handles zero padded FFTs. */

	if (gwdata->ZERO_PADDED_FFT) {
		unsigned long num_b_in_k, num_b_in_word_0, num_b_in_word_1;
		unsigned long num_b_in_word_2, num_b_in_word_3, num_b_in_word_4, num_b_in_word_5;

		asm_data->ZPAD_WORD5_OFFSET = addr_offset (gwdata, 4);
		if (asm_data->ZPAD_WORD5_OFFSET == 8) {  /* FFTLEN = 80 and 112 */
			asm_data->ZPAD_WORD5_RBP_OFFSET = 8;
		} else {
			asm_data->ZPAD_WORD5_RBP_OFFSET = 256;
		}

		asm_data->HIGH_WORD1_OFFSET = addr_offset (gwdata, gwdata->FFTLEN/2-1);
		asm_data->HIGH_WORD2_OFFSET = addr_offset (gwdata, gwdata->FFTLEN/2-2);
		asm_data->HIGH_WORD3_OFFSET = addr_offset (gwdata, gwdata->FFTLEN/2-3);

		raw_gwsetaddin (gwdata, gwdata->FFTLEN/2-1, 0.0);
		asm_data->HIGH_SCRATCH1_OFFSET = asm_data->ADDIN_OFFSET;
		raw_gwsetaddin (gwdata, gwdata->FFTLEN/2-2, 0.0);
		asm_data->HIGH_SCRATCH2_OFFSET = asm_data->ADDIN_OFFSET;
		raw_gwsetaddin (gwdata, gwdata->FFTLEN/2-3, 0.0);
		asm_data->HIGH_SCRATCH3_OFFSET = asm_data->ADDIN_OFFSET;

		num_b_in_k = (unsigned long) ceil (logb (k));
		num_b_in_word_0 = gwdata->NUM_B_PER_SMALL_WORD; if (is_big_word (gwdata, 0)) num_b_in_word_0++;
		num_b_in_word_1 = gwdata->NUM_B_PER_SMALL_WORD; if (is_big_word (gwdata, 1)) num_b_in_word_1++;
		num_b_in_word_2 = gwdata->NUM_B_PER_SMALL_WORD; if (is_big_word (gwdata, 2)) num_b_in_word_2++;
		num_b_in_word_3 = gwdata->NUM_B_PER_SMALL_WORD; if (is_big_word (gwdata, 3)) num_b_in_word_3++;
		num_b_in_word_4 = gwdata->NUM_B_PER_SMALL_WORD; if (is_big_word (gwdata, 4)) num_b_in_word_4++;
		num_b_in_word_5 = gwdata->NUM_B_PER_SMALL_WORD; if (is_big_word (gwdata, 5)) num_b_in_word_5++;

		asm_data->ZPAD_SHIFT1 = pow ((double) b, (int) num_b_in_word_0);
		asm_data->ZPAD_SHIFT2 = pow ((double) b, (int) num_b_in_word_1);
		asm_data->ZPAD_SHIFT3 = pow ((double) b, (int) num_b_in_word_2);
		asm_data->ZPAD_SHIFT4 = pow ((double) b, (int) num_b_in_word_3);
		asm_data->ZPAD_SHIFT5 = pow ((double) b, (int) num_b_in_word_4);
		asm_data->ZPAD_SHIFT6 = pow ((double) b, (int) num_b_in_word_5);

		if (num_b_in_k <= num_b_in_word_0) asm_data->ZPAD_TYPE = 1;
		else if (num_b_in_k <= num_b_in_word_0 + num_b_in_word_1) asm_data->ZPAD_TYPE = 2;
		else asm_data->ZPAD_TYPE = 3;

		if (asm_data->ZPAD_TYPE == 1) {
			asm_data->ZPAD_K1_LO = k;
			asm_data->ZPAD_INVERSE_K1 = 1.0 / k;
		}

		if (asm_data->ZPAD_TYPE == 2) {
			asm_data->ZPAD_K1_HI = floor (k / asm_data->ZPAD_SHIFT1);
			asm_data->ZPAD_K1_LO = k - asm_data->ZPAD_K1_HI * asm_data->ZPAD_SHIFT1;
			asm_data->ZPAD_INVERSE_K1 = asm_data->ZPAD_SHIFT1 / k;
			asm_data->ZPAD_K2_HI = floor (k / asm_data->ZPAD_SHIFT2);
			asm_data->ZPAD_K2_LO = k - asm_data->ZPAD_K2_HI * asm_data->ZPAD_SHIFT2;
			asm_data->ZPAD_INVERSE_K2 = asm_data->ZPAD_SHIFT2 / k;
			asm_data->ZPAD_K3_HI = floor (k / asm_data->ZPAD_SHIFT3);
			asm_data->ZPAD_K3_LO = k - asm_data->ZPAD_K3_HI * asm_data->ZPAD_SHIFT3;
			asm_data->ZPAD_INVERSE_K3 = asm_data->ZPAD_SHIFT3 / k;
			asm_data->ZPAD_K4_HI = floor (k / asm_data->ZPAD_SHIFT4);
			asm_data->ZPAD_K4_LO = k - asm_data->ZPAD_K4_HI * asm_data->ZPAD_SHIFT4;
			asm_data->ZPAD_INVERSE_K4 = asm_data->ZPAD_SHIFT4 / k;
			asm_data->ZPAD_K5_HI = floor (k / asm_data->ZPAD_SHIFT5);
			asm_data->ZPAD_K5_LO = k - asm_data->ZPAD_K5_HI * asm_data->ZPAD_SHIFT5;
			asm_data->ZPAD_INVERSE_K5 = asm_data->ZPAD_SHIFT5 / k;
			asm_data->ZPAD_K6_HI = floor (k / asm_data->ZPAD_SHIFT6);
			asm_data->ZPAD_K6_LO = k - asm_data->ZPAD_K6_HI * asm_data->ZPAD_SHIFT6;
			asm_data->ZPAD_INVERSE_K6 = asm_data->ZPAD_SHIFT6 / k;
		}

		if (asm_data->ZPAD_TYPE == 3) {
			double	powb, bigpowb;
			powb = pow ((double) b, (int) num_b_in_word_0);
			bigpowb = pow ((double) b, (int) (num_b_in_word_0 + num_b_in_word_1));
			asm_data->ZPAD_K2_HI = floor (k / bigpowb);
			asm_data->ZPAD_K2_MID = floor ((k - asm_data->ZPAD_K2_HI*bigpowb) / powb);
			asm_data->ZPAD_K2_LO = k - asm_data->ZPAD_K2_HI*bigpowb - asm_data->ZPAD_K2_MID*powb;
			asm_data->ZPAD_INVERSE_K2 = powb / k;
			powb = pow ((double) b, (int) num_b_in_word_1);
			bigpowb = pow ((double) b, (int) (num_b_in_word_1 + num_b_in_word_2));
			asm_data->ZPAD_K3_HI = floor (k / bigpowb);
			asm_data->ZPAD_K3_MID = floor ((k - asm_data->ZPAD_K3_HI*bigpowb) / powb);
			asm_data->ZPAD_K3_LO = k - asm_data->ZPAD_K3_HI*bigpowb - asm_data->ZPAD_K3_MID*powb;
			asm_data->ZPAD_INVERSE_K3 = powb / k;
			powb = pow ((double) b, (int) num_b_in_word_2);
			bigpowb = pow ((double) b, (int) (num_b_in_word_2 + num_b_in_word_3));
			asm_data->ZPAD_K4_HI = floor (k / bigpowb);
			asm_data->ZPAD_K4_MID = floor ((k - asm_data->ZPAD_K4_HI*bigpowb) / powb);
			asm_data->ZPAD_K4_LO = k - asm_data->ZPAD_K4_HI*bigpowb - asm_data->ZPAD_K4_MID*powb;
			asm_data->ZPAD_INVERSE_K4 = powb / k;
			powb = pow ((double) b, (int) num_b_in_word_3);
			bigpowb = pow ((double) b, (int) (num_b_in_word_3 + num_b_in_word_4));
			asm_data->ZPAD_K5_HI = floor (k / bigpowb);
			asm_data->ZPAD_K5_MID = floor ((k - asm_data->ZPAD_K5_HI*bigpowb) / powb);
			asm_data->ZPAD_K5_LO = k - asm_data->ZPAD_K5_HI*bigpowb - asm_data->ZPAD_K5_MID*powb;
			asm_data->ZPAD_INVERSE_K5 = powb / k;
			powb = pow ((double) b, (int) num_b_in_word_4);
			bigpowb = pow ((double) b, (int) (num_b_in_word_4 + num_b_in_word_5));
			asm_data->ZPAD_K6_HI = floor (k / bigpowb);
			asm_data->ZPAD_K6_MID = floor ((k - asm_data->ZPAD_K6_HI*bigpowb) / powb);
			asm_data->ZPAD_K6_LO = k - asm_data->ZPAD_K6_HI*bigpowb - asm_data->ZPAD_K6_MID*powb;
			asm_data->ZPAD_INVERSE_K6 = bigpowb / k;
		}
	}

/* Set the procedure pointers from the asm structures and proc tabes */

	memcpy (gwdata->GWPROCPTRS, info->proc_ptrs, 4 * sizeof (void *));
	memcpy (gwdata->GWPROCPTRS+4, info->add_sub_norm_procs, 12 * sizeof (void *));
#ifndef X86_64
	if (! (gwdata->cpu_flags & CPU_SSE2)) {
		gwdata->GWPROCPTRS[norm_routines] = x87_prctab[x87_prctab_index (gwdata, 0, 0)];  // No error, no mulbyconst
		gwdata->GWPROCPTRS[norm_routines+1] = x87_prctab[x87_prctab_index (gwdata, 1, 0)];  // Error, no mulbyconst
		gwdata->GWPROCPTRS[norm_routines+2] = x87_prctab[x87_prctab_index (gwdata, 0, 1)];  // No error, mulbyconst
		gwdata->GWPROCPTRS[norm_routines+3] = x87_prctab[x87_prctab_index (gwdata, 1, 1)];  // Error, mulbyconst
	} else
#endif
	{
		gwdata->GWPROCPTRS[norm_routines] = sse2_prctab[sse2_prctab_index (gwdata, 0, 0)];  // No error, no mulbyconst
		gwdata->GWPROCPTRS[norm_routines+1] = sse2_prctab[sse2_prctab_index (gwdata, 1, 0)];  // Error, no mulbyconst
		gwdata->GWPROCPTRS[norm_routines+2] = sse2_prctab[sse2_prctab_index (gwdata, 0, 1)];  // No error, mulbyconst
		gwdata->GWPROCPTRS[norm_routines+3] = sse2_prctab[sse2_prctab_index (gwdata, 1, 1)];  // Error, mulbyconst

	}
	memcpy (gwdata->GWPROCPTRS+copyzero_routines,
		  info->add_sub_norm_procs + 12 +
		  (gwdata->RATIONAL_FFT ? 0 : 2),
		  2 * sizeof (void *));

/* Default normalization routines and behaviors */

	gwsetnormroutine (gwdata, 0, 0, 0);
	gwstartnextfft (gwdata, 0);
	raw_gwsetaddin (gwdata, 0, 0.0);
	if (gwdata->square_carefully_count == -1) gwset_square_carefully_count (gwdata, -1);

/* Clear globals */

	asm_data->MAXERR = 0.0;
	asm_data->COPYZERO[0] = 0;
	gwdata->GWERROR = 0;
	gwdata->GW_RANDOM = NULL;

/* Compute maximum allowable difference for error checking */
/* This error check is disabled for mod B^N+1 arithmetic */

	if (gwdata->ALL_COMPLEX_FFT)
		gwdata->MAXDIFF = 1.0E80;

/* We have observed that the difference seems to vary based on the size */
/* the FFT result word.  This is two times the number of bits per double. */
/* Subtract 1 from bits per double because one bit is the sign bit. */
/* Add log2(b) for the FFT weight that range from 1 to b. */
/* Add in a percentage of the log(FFTLEN) to account for carries. */
/* We use a different threshold for SSE2 which uses 64-bit instead of */
/* 80-bit doubles during the FFT */

	else {
		double bits_per_double, total_bits, loglen;
		bits_per_double = gwdata->avg_num_b_per_word * log2 (b) - 1.0;
		if (!gwdata->RATIONAL_FFT)
			bits_per_double += log2 (b);
		if (!gwdata->ZERO_PADDED_FFT)
			bits_per_double += log2 (-c);
		loglen = log2 (gwdata->FFTLEN);
		loglen *= 0.69;
		total_bits = bits_per_double * 2.0 + loglen * 2.0;
		gwdata->MAXDIFF = pow ((double) 2.0, total_bits -
				((gwdata->cpu_flags & CPU_SSE2) ? 49.08 : 49.65));
	}

/* Clear counters, init internal timers */

	gwdata->fft_count = 0.0;
	asm_data->ASM_TIMERS = (uint32_t *) &gwdata->ASM_TIMERS;

/* Default size of gwnum_alloc array is 50 */

	gwdata->gwnum_alloc = NULL;
	gwdata->gwnum_alloc_count = 0;
	gwdata->gwnum_alloc_array_size = 50;
	gwdata->gwnum_free = NULL;
	gwdata->gwnum_free_count = 0;

/* Activate giants / gwnum shared cached memory allocates */

	gwdata->gdata.blksize = gwnum_datasize (gwdata);

/* Compute alignment for allocated data.  Strangely enough assembly */
/* prefetching works best in pass 1 on a P4 if the data is allocated */
/* on an odd cache line.  An optimal 31 of the 32 cache lines on a 4KB */
/* page will be prefetchable.  Page aligned data would only prefetch */
/* 28 of the 32 cache lines. */

	if (gwdata->cpu_flags & CPU_SSE2) {
		if (gwdata->PASS2_LEVELS == 0) {	/* One pass */
			gwdata->GW_ALIGNMENT = 128;	/* P4 cache line alignment */
			gwdata->GW_ALIGNMENT_MOD = 0;
		} else if (gwdata->SCRATCH_SIZE == 0) {	/* Small two passes */
			gwdata->GW_ALIGNMENT = 4096;	/* Page alignment */
			gwdata->GW_ALIGNMENT_MOD = 0;
		} else {			/* Large two passes */
			gwdata->GW_ALIGNMENT = 1024;	/* Clmblkdst (up to 8) */
			gwdata->GW_ALIGNMENT_MOD = 128; /* + 1 cache line */
		}
	} else {
		if (gwdata->PASS2_LEVELS == 0)		/* One pass */
			gwdata->GW_ALIGNMENT = 128;	/* P4 cache line alignment */
		else				/* Two passes */
			gwdata->GW_ALIGNMENT = 4096;	/* Page alignment */
		gwdata->GW_ALIGNMENT_MOD = 0;
	}

/* If we are going to use multiple threads for multiplications, then do */
/* the required multi-thread initializations.  Someday, we might allow */
/* setting num_threads after gwsetup so we put all the multi-thread */
/* initialization in its own routine. */

	return (multithread_init (gwdata));
}


/* Utility routines to deal with the 7 words near the half-way point */
/* in a zero-padded SSE2 FFT.  This used to be all done in assembly code, */
/* but I moved it to C code when the multithread code was written. */

/* When doing zero-padded FFTs, the 7 words around the halfway point must */
/* be copied for later processing.  This macro does that before a typical */
/* forward FFT. */

void xcopy_7_words (
	struct gwasm_data *asm_data)
{
	double	*srcarg, *destarg;

	destarg = (double *) asm_data->DESTARG;
	srcarg = (double *)
		((char *) destarg + (size_t) asm_data->DIST_TO_FFTSRCARG);
	destarg[-5] = srcarg[4];	/* Copy 1st word above halfway point */
	destarg[-6] = srcarg[12];	/* Copy 2nd word */
	destarg[-7] = srcarg[20];	/* Copy 3rd word */
	destarg[-8] = srcarg[28];	/* Copy 4th word */
	destarg[-9] = * (double *)	/* Copy 1st word below halfway point */
		((char *) srcarg + asm_data->HIGH_WORD1_OFFSET);
	destarg[-10] = * (double *)	/* Copy 2nd word below */
		((char *) srcarg + asm_data->HIGH_WORD2_OFFSET);
	destarg[-11] = * (double *)	/* Copy 3rd word below */
		((char *) srcarg + asm_data->HIGH_WORD3_OFFSET);
}

/* When POSTFFT is set, we must copy the 7 words at two different spots. */
/* These two routines copy the four values above the half-way point after */
/* carries have been propagated and copy the three words just below the */
/* half-way point right after the last NORMRTN has been called. */

void xcopy_4_words (
	struct gwasm_data *asm_data)
{
	double	*srcarg;

	srcarg = (double *) asm_data->DESTARG;
	srcarg[-5] = srcarg[4];		/* Copy 1st word above halfway point */
	srcarg[-6] = srcarg[12];	/* Copy 2nd word */
	srcarg[-7] = srcarg[20];	/* Copy 3rd word */
	srcarg[-8] = srcarg[28];	/* Copy 4th word */
}

void xcopy_3_words (
	struct gwasm_data *asm_data)
{
	gwhandle *gwdata;
	double	*srcarg;

/* Get pointer to the gwdata structure */

	gwdata = asm_data->gwdata;

/* Copy 1st and 2nd words below the halfway point from the last block */

	srcarg = (double *) asm_data->DESTARG;
	if (asm_data->this_block == asm_data->last_pass1_block) {
		if (gwdata->SCRATCH_SIZE) {
			srcarg[-9] = * (double *)
					((char *) asm_data->scratch_area +
					 asm_data->HIGH_SCRATCH1_OFFSET);
			srcarg[-10] = * (double *)
					((char *) asm_data->scratch_area +
					 asm_data->HIGH_SCRATCH2_OFFSET);
		} else {
			srcarg[-9] = * (double *)
					((char *) srcarg +
					 asm_data->HIGH_WORD1_OFFSET);
			srcarg[-10] = * (double *)
					((char *) srcarg +
					 asm_data->HIGH_WORD2_OFFSET);
		}
	}

/* Copy 3rd word below the halfway point */

	if (asm_data->this_block <= gwdata->num_pass1_blocks - 3 &&
	    asm_data->this_block + asm_data->cache_line_multiplier >
						gwdata->num_pass1_blocks - 3) {
		if (gwdata->SCRATCH_SIZE) {
			srcarg[-11] = * (double *)
					((char *) asm_data->scratch_area +
					 asm_data->HIGH_SCRATCH3_OFFSET);
		} else {
			srcarg[-11] = * (double *)
					((char *) srcarg +
					 asm_data->HIGH_WORD3_OFFSET);
		}
	}
}


/* Inline routines for the processing blocks in multi-thread code below */

/* Calculate pass 1 block address.  Note that block addresses are */
/* computed pased on 128 pad bytes after every 128 cache lines.  That is: */
/* DESTARG + (block * 64) + (block >> 7 * 128) */

__inline void *pass1_data_addr (
	struct gwasm_data *asm_data,
	unsigned long block)
{
	return ((char *) asm_data->DESTARG +
		(block << 6) + ((block >> 7) << 7));
}

/* Calculate pass 2 block address */

__inline void *pass2_data_addr (
	gwhandle *gwdata,
	struct gwasm_data *asm_data,
	unsigned long block)
{
	return ((char *) asm_data->DESTARG +
		block * (unsigned long) asm_data->pass1blkdst);
}

/* Calculate pass 2 premultiplier address */

__inline void * pass2_premult_addr (
	gwhandle *gwdata,
	unsigned long block)
{
	return ((char *) gwdata->adjusted_pass2_premults +
		block * gwdata->pass2_premult_block_size);
}

/* Assign next available block in state 0 of pass 1.  These are assigned in */
/* sequential order. */

__inline void pass1_state0_assign_next_block (
	gwhandle *gwdata,
	struct gwasm_data *asm_data)
{
	asm_data->next_block = gwdata->next_block;
	asm_data->data_prefetch =
		pass1_data_addr (asm_data, asm_data->next_block);
	gwdata->next_block += asm_data->cache_line_multiplier;
}

/* Assign first block in pass 2. */

__inline void pass2_assign_first_block (
	gwhandle *gwdata,
	struct gwasm_data *asm_data)
{
	asm_data->this_block = 0;
	asm_data->data_addr = pass2_data_addr (gwdata, asm_data, 0);
	asm_data->premult_addr = pass2_premult_addr (gwdata, 0);
	asm_data->next_block = 1;
	asm_data->data_prefetch = pass2_data_addr (gwdata, asm_data, 1);
	asm_data->premult_prefetch = pass2_premult_addr (gwdata, 1);
	gwdata->next_block = 2;
}

/* Assign next available block in pass 2.  These are assigned in */
/* sequential order. */

__inline void pass2_assign_next_block (
	gwhandle *gwdata,
	struct gwasm_data *asm_data)
{
	asm_data->next_block = gwdata->next_block++;
	asm_data->data_prefetch =
		pass2_data_addr (gwdata, asm_data, asm_data->next_block);
	asm_data->premult_prefetch =
		pass2_premult_addr (gwdata, asm_data->next_block);
}


/* Routine for auxillary threads */

struct thread_data {
	gwhandle *gwdata;
	int	thread_num;
	void	*asm_data_alloc;
	void	*scratch_area;
};

void auxillary_thread (void *arg)
{
	gwhandle *gwdata;
	struct gwasm_data *asm_data, *main_thread_asm_data;
	struct thread_data *info;

/* Get pointers to various structures */

	info = (struct thread_data *) arg;
	gwdata = info->gwdata;
	main_thread_asm_data = (struct gwasm_data *) gwdata->asm_data;

/* Each thread needs its own copy of the asm_data.  Each thread needs */
/* its own stack too. */

	asm_data = (struct gwasm_data *)
		((char *) info->asm_data_alloc + NEW_STACK_SIZE);
	memcpy (asm_data, main_thread_asm_data, sizeof (struct gwasm_data));

/* Set the thread number so that the assembly code can differentiate between */
/* the main thread and an auxillary thread. */

	asm_data->thread_num = info->thread_num;

/* Each auxillary thread needs its own pass 1 scratch area */

	asm_data->scratch_area = info->scratch_area;

/* Call optional user provided callback routine so that the caller can set */
/* the thread's priority and affinity */

	if (gwdata->thread_callback != NULL)
		(*gwdata->thread_callback) (
			info->thread_num, 0,
			gwdata->thread_callback_data);

/* Loop waiting for work to do.  The main thread will signal the */
/* thread_work_to_do event whenever there is work for the auxillary */
/* thread to do. */

	for ( ; ; ) {
		gwevent_wait (&gwdata->thread_work_to_do, 0);
		gwmutex_lock (&gwdata->thread_lock);

/* If threads are to exit, break out of this work loop */

		if (gwdata->threads_must_exit) {
			gwmutex_unlock (&gwdata->thread_lock);
			break;
		}

/* Copy the main thread's asm_data's DESTARG for proper next_block */
/* address calculations.  We'll copy more asm_data later. */

		asm_data->DESTARG = main_thread_asm_data->DESTARG;

/* Get an available block for this thread to process (store it in the */
/* next_block field).  Note we set this_block to a dummy value so that */
/* get_next_block knows there is more work to do.  NOTE: There is no */
/* guarantee that there is an available block to process. */

		asm_data->this_block = 0;
		if (gwdata->pass1_state == 0) {
			if (gwdata->next_block >= gwdata->num_pass1_blocks)
				goto aux_done;
			pass1_state0_assign_next_block (gwdata, asm_data);
		} else if (gwdata->pass1_state == 1) {
			asm_data->next_block =
				asm_data->thread_num *
					asm_data->cache_line_multiplier;
			asm_data->data_prefetch =
				pass1_data_addr (asm_data, asm_data->next_block);
		} else {
			if (gwdata->next_block >= gwdata->num_pass2_blocks)
				goto aux_done;
			pass2_assign_next_block (gwdata, asm_data);
		}

/* Increment the number of active threads so that we can tell */
/* when all auxillary routines have finished. */

		gwdata->num_active_threads++;
		gwevent_reset (&gwdata->all_threads_done);

/* Copy some data from the main thread's asm_data to this thread.  We only */
/* need to copy data that can change for multiply call */

//bug - copy the necessary bytes using memcpy (and only once, not twice (or three times), per multiply)
asm_data->DIST_TO_FFTSRCARG = main_thread_asm_data->DIST_TO_FFTSRCARG;
asm_data->DIST_TO_MULSRCARG = main_thread_asm_data->DIST_TO_MULSRCARG;
asm_data->ffttype = main_thread_asm_data->ffttype;
///bug-and these are used in pass 1
asm_data->thread_work_routine = main_thread_asm_data->thread_work_routine;
asm_data->NORMRTN = main_thread_asm_data->NORMRTN;
asm_data->ADDIN_ROW = main_thread_asm_data->ADDIN_ROW;
asm_data->ADDIN_OFFSET = main_thread_asm_data->ADDIN_OFFSET;
asm_data->ADDIN_VALUE = main_thread_asm_data->ADDIN_VALUE;
asm_data->XMM_MULCONST[0] = main_thread_asm_data->XMM_MULCONST[0];
asm_data->XMM_MULCONST[1] = main_thread_asm_data->XMM_MULCONST[1];
//if zpad more mulconsts must be copied
asm_data->XMM_K_TIMES_MULCONST_HI[0] = main_thread_asm_data->XMM_K_TIMES_MULCONST_HI[0];
asm_data->XMM_K_TIMES_MULCONST_HI[1] = main_thread_asm_data->XMM_K_TIMES_MULCONST_HI[1];
asm_data->XMM_K_TIMES_MULCONST_LO[0] = main_thread_asm_data->XMM_K_TIMES_MULCONST_LO[0];
asm_data->XMM_K_TIMES_MULCONST_LO[1] = main_thread_asm_data->XMM_K_TIMES_MULCONST_LO[1];
asm_data->XMM_MINUS_C_TIMES_MULCONST[0] = main_thread_asm_data->XMM_MINUS_C_TIMES_MULCONST[0];
asm_data->XMM_MINUS_C_TIMES_MULCONST[1] = main_thread_asm_data->XMM_MINUS_C_TIMES_MULCONST[1];
//postfft is used in c callback routines
asm_data->POSTFFT = main_thread_asm_data->POSTFFT;

/* Now call the assembly code to do some work! */

		gwmutex_unlock (&gwdata->thread_lock);
		if (gwdata->pass1_state < 999)
			pass1_aux_entry_point (asm_data);
		else
			pass2_aux_entry_point (asm_data);
		gwmutex_lock (&gwdata->thread_lock);

/* The auxillary threads have run out of work.  Decrement the count of */
/* number of active auxillary threads.  Signal all threads done when */
/* last auxillary thread is done. */

		gwdata->num_active_threads--;
		if (gwdata->num_active_threads == 0)
			gwevent_signal (&gwdata->all_threads_done);

/* Reset thread_work_to_do event before looping to wait for more work. */

aux_done:	gwevent_reset (&gwdata->thread_work_to_do);
		gwmutex_unlock (&gwdata->thread_lock);
	}

/* Call optional user provided callback routine so that the caller can set */
/* do any necessary cleanup */

	if (gwdata->thread_callback != NULL)
		(*gwdata->thread_callback) (
			info->thread_num, 1,
			gwdata->thread_callback_data);

/* Free the allocated memory and exit the auxillary thread */

	aligned_free (info->scratch_area);
	aligned_free (info->asm_data_alloc);
	free (arg);
}



/* This routine is called by the main thread assembly code to */
/* fire up all the auxillary worker threads in pass 1. */

void pass1_wake_up_threads (
	struct gwasm_data *asm_data)
{
	gwhandle *gwdata;

/* Get pointer to the gwdata structure */

	gwdata = asm_data->gwdata;

/* Init pass1_state (kludge: passed in next_block parameter) */
/* State is either 0 (forward FFT only) or 1 (inverse FFT and */
/* optional next forward FFT) */

	gwdata->pass1_state = asm_data->next_block;

/* Prior to doing a forward zero-padded FFT copy 7 FFT elements */

	if (gwdata->pass1_state == 0 && gwdata->ZERO_PADDED_FFT)
		xcopy_7_words (asm_data);

/* Initialize the pointers that the normalization / carry propagation code */
/* use. */

	asm_data->norm_ptr1 = asm_data->norm_biglit_array;
	asm_data->norm_ptr2 = asm_data->norm_col_mults;
	asm_data->XMM_SUMOUT[0] = 0.0;
	asm_data->XMM_SUMOUT[1] = 0.0;
	asm_data->XMM_MAXERR[0] = asm_data->MAXERR;
	asm_data->XMM_MAXERR[1] = asm_data->MAXERR;

/* In the zero-padded FFT case prior to normalizing, the upper half of the */
/* carries are set to zero instead of XMM_BIGVAL.  This saves us a couple */
/* of clocks per FFT data element in xnorm_2d_zpad.  We zero the upper half */
/* of each cache line. */

	if (gwdata->pass1_state == 1 && gwdata->ZERO_PADDED_FFT) {
		int	i, carry_table_size;
		double	*table;
		carry_table_size = ((gwdata->FFTLEN >> gwdata->PASS2_LEVELS) << 1);
		table = (double *) asm_data->carries;
		for (i = 0; i < carry_table_size; i += 8, table += 8) {
			table[4] = 0.0;
			table[5] = 0.0;
			table[6] = 0.0;
			table[7] = 0.0;
		}
	}

/* Set up the this_block and next_block values for the main thread */

	asm_data->this_block = 0;
	asm_data->data_addr = pass1_data_addr (asm_data, 0);
	if (gwdata->pass1_state == 0) {
		gwdata->next_block = asm_data->cache_line_multiplier;
		pass1_state0_assign_next_block (gwdata, asm_data);
	} else {
		asm_data->next_block = asm_data->cache_line_multiplier;
		asm_data->data_prefetch =
			pass1_data_addr (asm_data, asm_data->next_block);
		gwdata->next_block = 0;
	}
}

void pass1_wake_up_threads_mt (
	struct gwasm_data *asm_data)
{
	gwhandle *gwdata;
	unsigned long i;

/* Get pointer to the gwdata structure and obtain lock */

	gwdata = asm_data->gwdata;
	gwmutex_lock (&gwdata->thread_lock);

/* Init pass1_state (kludge: passed in next_block parameter) */
/* State is either 0 (forward FFT only) or 1 (inverse FFT and */
/* optional next forward FFT) */

	gwdata->pass1_state = asm_data->next_block;

/* Prior to doing a forward zero-padded FFT copy 7 FFT elements */

	if (gwdata->pass1_state == 0 && gwdata->ZERO_PADDED_FFT)
		xcopy_7_words (asm_data);

/* Initialize the pointers that the normalization / carry propagation code */
/* use.  These are copied before and after each block is normalized so that */
/* each thread essentially shares these vallues. */

	asm_data->norm_ptr1 = asm_data->norm_biglit_array;
	asm_data->norm_ptr2 = asm_data->norm_col_mults;
	asm_data->XMM_SUMOUT[0] = 0.0;
	asm_data->XMM_SUMOUT[1] = 0.0;
	asm_data->XMM_MAXERR[0] = asm_data->MAXERR;
	asm_data->XMM_MAXERR[1] = asm_data->MAXERR;

/* In the zero-padded FFT case prior to normalizing, the upper half of the */
/* carries are set to zero instead of XMM_BIGVAL.  This saves us a couple */
/* of clocks per FFT data element in xnorm_2d_zpad.  We zero the upper half */
/* of each cache line. */

	if (gwdata->pass1_state == 1 && gwdata->ZERO_PADDED_FFT) {
		int	i, carry_table_size;
		double	*table;
		carry_table_size = ((gwdata->FFTLEN >> gwdata->PASS2_LEVELS) << 1);
		table = (double *) asm_data->carries;
		for (i = 0; i < carry_table_size; i += 8, table += 8) {
			table[4] = 0.0;
			table[5] = 0.0;
			table[6] = 0.0;
			table[7] = 0.0;
		}
	}

/* Set up the this_block and next_block values for the main thread */

	asm_data->this_block = 0;
	asm_data->data_addr = pass1_data_addr (asm_data, 0);
	if (gwdata->pass1_state == 0) {
		gwdata->next_block = asm_data->cache_line_multiplier;
		pass1_state0_assign_next_block (gwdata, asm_data);
	} else {
		asm_data->next_block = asm_data->cache_line_multiplier *
					gwdata->num_threads;
		asm_data->data_prefetch =
			pass1_data_addr (asm_data, asm_data->next_block);
		gwdata->next_block = 0;
	}

/* Set flags and events, set assembly routine ptr that auxilllary threads */
/* should call to do work, signal the auxillary threads to resume working */

	gwevent_signal (&gwdata->pass1_norm_events[0]);
	for (i = 1; i < gwdata->num_threads; i++)
		gwevent_reset (&gwdata->pass1_norm_events[i]);
	gwevent_reset (&gwdata->gwcarries_complete);
	gwevent_signal (&gwdata->thread_work_to_do);
	gwmutex_unlock (&gwdata->thread_lock);
}

/* This routine is called by the main thread assembly code to */
/* wait for all the previous blocks to finish with the carry data. */

void pass1_wait_for_carries (
	struct gwasm_data *asm_data)
{
}

void pass1_wait_for_carries_mt (
	struct gwasm_data *asm_data)
{
	gwhandle *gwdata;
	struct gwasm_data *main_thread_asm_data;

/* Get pointer to the gwdata structure */

	gwdata = asm_data->gwdata;

/* Wait until the next normalization block to process is the block */
/* that this thread is working on. */

	gwevent_wait (&gwdata->pass1_norm_events[asm_data->thread_num], 0);
	gwevent_reset (&gwdata->pass1_norm_events[asm_data->thread_num]);

/* Copy the required normalization values so that assembly code can */
/* access them */

	main_thread_asm_data = (struct gwasm_data *) gwdata->asm_data;
	asm_data->norm_ptr1 = main_thread_asm_data->norm_ptr1;
	asm_data->norm_ptr2 = main_thread_asm_data->norm_ptr2;
	asm_data->XMM_SUMOUT[0] = main_thread_asm_data->XMM_SUMOUT[0];
	asm_data->XMM_SUMOUT[1] = main_thread_asm_data->XMM_SUMOUT[1];
	asm_data->XMM_MAXERR[0] = main_thread_asm_data->XMM_MAXERR[0];
	asm_data->XMM_MAXERR[1] = main_thread_asm_data->XMM_MAXERR[1];
}

/* This routine is called by the pass 1 threads when they are done */
/* with the normalization code. */

#define PASS1_CARRIES_NO_FORWARD_FFT	0
#define PASS1_CARRIES_FORWARD_FFT	1

int pass1_done_with_carries (
	struct gwasm_data *asm_data)
{
	gwhandle *gwdata;

/* Get pointer to the gwdata structure */

	gwdata = asm_data->gwdata;

/* If postfft is not set, then do not forward FFT the result */

	if (!asm_data->POSTFFT) return (PASS1_CARRIES_NO_FORWARD_FFT);

/* If this is a zero padded FFT and we are doing the last block, then */
/* we need to copy a few of the data values before they are FFTed. */

	if (gwdata->ZERO_PADDED_FFT)
		xcopy_3_words (asm_data);

/* We must delay the forward FFT for the first few data blocks (until the */
/* carries are added back in). */

	if (asm_data->this_block < 8) return (PASS1_CARRIES_NO_FORWARD_FFT);
	return (PASS1_CARRIES_FORWARD_FFT);
}

int pass1_done_with_carries_mt (
	struct gwasm_data *asm_data)
{
	gwhandle *gwdata;
	struct gwasm_data *main_thread_asm_data;

/* Get pointer to the gwdata structure */

	gwdata = asm_data->gwdata;

/* Copy back the required normalization values so that the next thread */
/* can pick up where we left off.  This is a nop if this is the main thread. */

	main_thread_asm_data = (struct gwasm_data *) gwdata->asm_data;
	main_thread_asm_data->norm_ptr1 = asm_data->norm_ptr1;
	main_thread_asm_data->norm_ptr2 = asm_data->norm_ptr2;
	main_thread_asm_data->XMM_SUMOUT[0] = asm_data->XMM_SUMOUT[0];
	main_thread_asm_data->XMM_SUMOUT[1] = asm_data->XMM_SUMOUT[1];
	main_thread_asm_data->XMM_MAXERR[0] = asm_data->XMM_MAXERR[0];
	main_thread_asm_data->XMM_MAXERR[1] = asm_data->XMM_MAXERR[1];

/* Inform the next thread that it can now enter the normalization code */

	if (asm_data->this_block == asm_data->last_pass1_block)
		gwevent_signal (&gwdata->pass1_norm_events[0]);
	else if (asm_data->thread_num == gwdata->num_threads - 1)
		gwevent_signal (&gwdata->pass1_norm_events[0]);
	else
		gwevent_signal (&gwdata->pass1_norm_events[asm_data->thread_num + 1]);

/* If postfft is not set, then do not forward FFT the result */

	if (!asm_data->POSTFFT) return (PASS1_CARRIES_NO_FORWARD_FFT);
	
/* If this is a zero padded FFT and we are doing the last block, then */
/* we need to copy a few of the data values before they are FFTed. */

	if (gwdata->ZERO_PADDED_FFT)
		xcopy_3_words (asm_data);

/* We must delay the forward FFT for the first few data blocks (until the */
/* carries are added back in). */

	if (asm_data->this_block < gwdata->num_postfft_blocks)
		return (PASS1_CARRIES_NO_FORWARD_FFT);
	return (PASS1_CARRIES_FORWARD_FFT);
}

/* This routine is called by assembly code threads to get the next */
/* pass 1 block to process.  It returns a code indicating what part */
/* of the pass 1 to do next. */

/* Return codes: */

#define PASS1_DO_MORE_INVERSE_FFT	0
#define PASS1_DO_MORE_FORWARD_FFT	1
#define PASS1_COMPLETE			2
#define PASS1_EXIT_THREAD		3
#define PASS1_START_PASS2		4
#define PASS1_DO_GWCARRIES		5

/* Pass 1 states: */
/* 0 = forward fft */
/* 1 = inverse fft */
/* 2 = gwcarries has started */
/* 3 = gwcarries has completed  (may use event instead) */
/* 999 = not in pass 1, we're doing pass 2 */

int pass1_get_next_block (
	struct gwasm_data *asm_data)
{
	gwhandle *gwdata;

/* Get pointer to the gwdata structure */

	gwdata = asm_data->gwdata;

/* If state is zero, then we are doing the forward FFT.  In this case we */
/* can process the blocks in any order (much like we do in pass 2). */

	if (gwdata->pass1_state == 0) {

/* There are no more blocks to process when the next block is the same */
/* as the block just processed.  In that case, if this is the main thread */
/* wait for auxillary threads to complete.  If this is an auxillary thread */
/* then return code telling the assembly code to exit. */

		if (asm_data->this_block == asm_data->next_block) {
			asm_data->DIST_TO_FFTSRCARG = 0;
			return (PASS1_START_PASS2);
		}

/* There is more pass 1 work to do.  Get next available block (if any). */

		asm_data->this_block = asm_data->next_block;
		asm_data->data_addr = asm_data->data_prefetch;
		if (gwdata->next_block < gwdata->num_pass1_blocks) {
			pass1_state0_assign_next_block (gwdata, asm_data);
		}

		return (PASS1_DO_MORE_FORWARD_FFT);
	}

/* Otherwise, pass1_state is one and we are doing the inverse FFT (and */
/* if POSTFFT is set pass 1 of the forward FFT on the result).  In this */
/* case, the carry propgation code will force us to process blocks in order. */

/* Handle the uncommon cases.  These occur when there is no next block */
/* (next_block == this_block) or when we are processing the postfft data */
/* blocks (next_block < this_block) */

	if (asm_data->next_block <= asm_data->this_block) {

/* If gwcarries has not been called, then wait for all the carry */
/* propagations to complete and then have the assembly code do gwcarries */
/* For simplicity sake, we make the main thread do the gwcarries (I'd guess */
/* that 99.99% of the time the main thread will be idle by the time the */
/* last block has been normalized). */

		if (gwdata->pass1_state < 2) {
			gwdata->pass1_state = 2;
			return (PASS1_DO_GWCARRIES);
		}

/* Set state indicating gwcarries has completed. */ 
/* After the carries are added back in, save the lowest four words of a */
/* zero padded FFTs (before POSTFFT processing destroys the data). */

		if (gwdata->pass1_state < 3) {
			gwdata->pass1_state = 3;
			if (gwdata->ZERO_PADDED_FFT)
				xcopy_4_words (asm_data);
		}

/* There are no more blocks to process when the next block is the same */
/* as the block just processed.  Return code telling the assembly code */
/* to exit. */

		if (asm_data->this_block == asm_data->next_block)
			return (PASS1_COMPLETE);

/* We are processing a post FFT block.  Fall through to figure out next */
/* block to prefetch. */

	}

/* Calculate the next inverse FFT block to work on */

	asm_data->this_block = asm_data->next_block;
	asm_data->data_addr = asm_data->data_prefetch;
	if (gwdata->pass1_state == 1) {
		unsigned long next_block;
		next_block = asm_data->this_block +
				asm_data->cache_line_multiplier;
		if (next_block < gwdata->num_pass1_blocks) {
			asm_data->next_block = next_block;
			asm_data->data_prefetch =
				pass1_data_addr (asm_data, next_block);
			return (PASS1_DO_MORE_INVERSE_FFT);
		}
	}

/* Either there is no next block or it will be a POSTFFT block */

	if (asm_data->POSTFFT && gwdata->next_block < gwdata->num_postfft_blocks) {
		asm_data->next_block = gwdata->next_block;
		asm_data->data_prefetch =
			pass1_data_addr (asm_data, asm_data->next_block);
		gwdata->next_block += asm_data->cache_line_multiplier;
	}
	
/* Return indicating whether to jump to the inverse fft code or for */
/* POSTFFT processing jusmp to the forward fft code */

	if (gwdata->pass1_state == 1) {
		return (PASS1_DO_MORE_INVERSE_FFT);
	} else {
		return (PASS1_DO_MORE_FORWARD_FFT);
	}
}

int pass1_get_next_block_mt (
	struct gwasm_data *asm_data)
{
	gwhandle *gwdata;

/* Get pointer to the gwdata structure.  Grab lock before reading or */
/* writing any gwdata values. */

	gwdata = asm_data->gwdata;
	gwmutex_lock (&gwdata->thread_lock);

/* If state is zero, then we are doing the forward FFT.  In this case we */
/* can process the blocks in any order (much like we do in pass 2). */

	if (gwdata->pass1_state == 0) {

/* There are no more blocks to process when the next block is the same */
/* as the block just processed.  In that case, if this is the main thread */
/* wait for auxillary threads to complete.  If this is an auxillary thread */
/* then return code telling the assembly code to exit. */

		if (asm_data->this_block == asm_data->next_block) {
			if (asm_data->thread_num) {
				gwmutex_unlock (&gwdata->thread_lock);
				return (PASS1_EXIT_THREAD);
			}
			gwmutex_unlock (&gwdata->thread_lock);
			gwevent_wait (&gwdata->all_threads_done, 0);
			asm_data->DIST_TO_FFTSRCARG = 0;
			return (PASS1_START_PASS2);
		}

/* There is more pass 1 work to do.  Obtain lock and get next available */
/* block (if any). */

		asm_data->this_block = asm_data->next_block;
		asm_data->data_addr = asm_data->data_prefetch;
		if (gwdata->next_block < gwdata->num_pass1_blocks)
			pass1_state0_assign_next_block (gwdata, asm_data);
		gwmutex_unlock (&gwdata->thread_lock);

		return (PASS1_DO_MORE_FORWARD_FFT);
	}

/* Otherwise, pass1_state is one and we are doing the inverse FFT (and */
/* if POSTFFT is set pass 1 of the forward FFT on the result).  In this */
/* case, the carry propgation code will force us to process blocks in order. */

/* Handle the uncommon cases.  These occur when there is no next block */
/* (next_block == this_block) or when we are processing the postfft data */
/* blocks (next_block < this_block) */

	if (asm_data->next_block <= asm_data->this_block) {

/* If gwcarries has not been called, then wait for all the carry */
/* propagations to complete and then have the assembly code do gwcarries */
/* For simplicity sake, we make the main thread do the gwcarries (I'd guess */
/* that 99.99% of the time the main thread will be idle by the time the */
/* last block has been normalized). */

		if (gwdata->pass1_state < 2 && asm_data->thread_num == 0) {
			gwmutex_unlock (&gwdata->thread_lock);
			gwevent_wait (&gwdata->pass1_norm_events[0], 0);
			gwmutex_lock (&gwdata->thread_lock);
			gwdata->pass1_state = 2;
			gwmutex_unlock (&gwdata->thread_lock);
			return (PASS1_DO_GWCARRIES);
		}

/* Set state indicating gwcarries has completed. */ 
/* After the carries are added back in, save the lowest four words of a */
/* zero padded FFTs (before POSTFFT processing destroys the data). */

		if (gwdata->pass1_state < 3 && asm_data->thread_num == 0) {
			gwdata->pass1_state = 3;
			if (gwdata->ZERO_PADDED_FFT)
				xcopy_4_words (asm_data);
			gwevent_signal (&gwdata->gwcarries_complete);
		}

/* There are no more blocks to process when the next block is the same */
/* as the block just processed.  In that case, if this is the main thread */
/* wait for auxillary threads to complete.  If this is an auxillary thread */
/* then return code telling the assembly code to exit. */

		if (asm_data->this_block == asm_data->next_block) {
			if (asm_data->thread_num) {
				gwmutex_unlock (&gwdata->thread_lock);
				return (PASS1_EXIT_THREAD);
			}
			gwmutex_unlock (&gwdata->thread_lock);
			gwevent_wait (&gwdata->all_threads_done, 0);
			return (PASS1_COMPLETE);
		}

/* We are processing a post FFT block.  Wait for the gwcarries code to */
/* complete, then fall through to figure out next block to prefetch. */
		
		gwmutex_unlock (&gwdata->thread_lock);
		gwevent_wait (&gwdata->gwcarries_complete, 0);
		gwmutex_lock (&gwdata->thread_lock);
	}

/* Calculate the next inverse FFT block to work on */

	asm_data->this_block = asm_data->next_block;
	asm_data->data_addr = asm_data->data_prefetch;
	if (gwdata->pass1_state == 1) {
		unsigned long next_block;
		next_block = asm_data->this_block +
				asm_data->cache_line_multiplier *
				gwdata->num_threads;
		if (next_block < gwdata->num_pass1_blocks) {
			asm_data->next_block = next_block;
			asm_data->data_prefetch =
				pass1_data_addr (asm_data, next_block);
			gwmutex_unlock (&gwdata->thread_lock);
			return (PASS1_DO_MORE_INVERSE_FFT);
		}
	}

/* Either there is no next block or it will be a POSTFFT block */

	if (asm_data->POSTFFT && gwdata->next_block < gwdata->num_postfft_blocks) {
		asm_data->next_block = gwdata->next_block;
		asm_data->data_prefetch =
			pass1_data_addr (asm_data, asm_data->next_block);
		gwdata->next_block += asm_data->cache_line_multiplier;
	}
	
/* Return indicating whether to jump to the inverse fft code or for */
/* POSTFFT processing jump to the forward fft code */

	if (gwdata->pass1_state == 1) {
		gwmutex_unlock (&gwdata->thread_lock);
		return (PASS1_DO_MORE_INVERSE_FFT);
	} else {
		gwmutex_unlock (&gwdata->thread_lock);
		return (PASS1_DO_MORE_FORWARD_FFT);
	}
}


/* This callback routine is called by the main thread assembly code to */
/* fire up all the auxillary worker threads in pass 2. */

void pass2_wake_up_threads (
	struct gwasm_data *asm_data)
{
	gwhandle *gwdata;

/* Call assign_next_block twice to set up the main thread's this_block */
/* and next_block values. */

	gwdata = asm_data->gwdata;
	gwdata->pass1_state = 999;
	pass2_assign_first_block (gwdata, asm_data);
}

/* This is the multi-thread version of the routine above */

void pass2_wake_up_threads_mt (
	struct gwasm_data *asm_data)
{
	gwhandle *gwdata;

/* Get pointer to the gwdata structure and obtain lock */

	gwdata = asm_data->gwdata;
	gwmutex_lock (&gwdata->thread_lock);

/* Do all the same work the single thread version does to set up the */
/* main thread's data block and prefetch block. */

	gwdata->pass1_state = 999;
	pass2_assign_first_block (gwdata, asm_data);

/* Set flags and events, set assembly routine ptr that auxilllary threads */
/* should call to do work,  signal the auxillary threads to resume working */

	gwevent_signal (&gwdata->thread_work_to_do);
	gwmutex_unlock (&gwdata->thread_lock);
}

/* This routine is called by assembly code threads to get the next */
/* pass 2 block to process. */

int pass2_get_next_block (
	struct gwasm_data *asm_data)
{
	gwhandle *gwdata;

	if (asm_data->this_block == asm_data->next_block) {
		asm_data->DIST_TO_FFTSRCARG = 0;
		return (TRUE);
	}

	gwdata = asm_data->gwdata;

	asm_data->this_block = asm_data->next_block;
	asm_data->data_addr = asm_data->data_prefetch;
	asm_data->premult_addr = asm_data->premult_prefetch;
	if (gwdata->next_block < gwdata->num_pass2_blocks)
		pass2_assign_next_block (gwdata, asm_data);

	return (FALSE);
}

/* This is the multi-thread version of the routine above */

int pass2_get_next_block_mt (
	struct gwasm_data *asm_data)
{
	gwhandle *gwdata;

/* Get pointer to the gwdata structure.  Grab lock before reading or */
/* writing any gwdata values. */

	gwdata = asm_data->gwdata;
	gwmutex_lock (&gwdata->thread_lock);

/* There are no more blocks to process when the next block is the same */
/* as the block just processed.  In that case, if this is the main thread */
/* wait for auxillary threads to complete.  If this is an auxillary thread */
/* then return code telling the assembly code to exit. */

	if (asm_data->this_block == asm_data->next_block) {
		if (asm_data->thread_num == 0) {
			gwmutex_unlock (&gwdata->thread_lock);
			gwevent_wait (&gwdata->all_threads_done, 0);
			asm_data->DIST_TO_FFTSRCARG = 0;
		} else
			gwmutex_unlock (&gwdata->thread_lock);
		return (TRUE);
	}

/* Copy prefetched block and addresses to this block.  Get next available */
/* block to prefetch. */

	asm_data->this_block = asm_data->next_block;
	asm_data->data_addr = asm_data->data_prefetch;
	asm_data->premult_addr = asm_data->premult_prefetch;
	if (gwdata->next_block < gwdata->num_pass2_blocks) {
		pass2_assign_next_block (gwdata, asm_data);
	}
	gwmutex_unlock (&gwdata->thread_lock);

/* Return code indicating more work to do */

	return (FALSE);
}


/* Perform initializations required for multi-threaded operation */

int multithread_init (
	gwhandle *gwdata)
{
	struct gwasm_data *asm_data;
	unsigned long i;

/* Get pointer to assembly structure */

	asm_data = (struct gwasm_data *) gwdata->asm_data;

/* Save gwdata pointer in asm_data so that C callback routines can */
/* access gwdata.  Set flag indicating this is the main thread. */

	asm_data->gwdata = gwdata;
	asm_data->thread_num = 0;

/* Init other variables */

	gwdata->num_pass1_blocks = (1 << gwdata->PASS2_LEVELS) >> 1;
	gwdata->num_pass2_blocks = gwdata->FFTLEN >> gwdata->PASS2_LEVELS >> 2;
	asm_data->last_pass1_block =
		gwdata->num_pass1_blocks - asm_data->cache_line_multiplier;
//bug - minor optimization:  for small k/c we can probably set this to
//less than 8
	gwdata->num_postfft_blocks = 8;

/* Calculate the values used to compute pass 2 premultier pointers. */
/* We calculate an adjusted starting address of the premultiplier data */
/* so that both real and all-complex FFTs can use the same formula to */
/* calculate the proper address given a block number. */

	gwdata->pass2_premult_block_size =
		(gwdata->PASS2_LEVELS == 8) ? 32 * 128 :
		(gwdata->PASS2_LEVELS == 10) ? 64 * 128 :
		(gwdata->PASS2_LEVELS == 11) ? 96 * 128 :
		(gwdata->PASS2_LEVELS == 12) ? 128 * 128 :
		(gwdata->PASS2_LEVELS == 13) ? 192 * 128 : 0;
	if (gwdata->ALL_COMPLEX_FFT)
		gwdata->adjusted_pass2_premults = asm_data->pass2_premults;
	else
		gwdata->adjusted_pass2_premults =
			(char *) asm_data->pass2_premults -
			gwdata->pass2_premult_block_size;

/* Only two pass SSE2 FFTs support multi-threaded execution */

	if (gwdata->PASS2_LEVELS == 0 || !(gwdata->cpu_flags & CPU_SSE2))
		gwdata->num_threads = 1;

/* If we aren't multithreading, use the simpler version of routines */

	if (gwdata->num_threads <= 1) {
		asm_data->pass1_wake_up_threads = pass1_wake_up_threads;
		asm_data->pass1_wait_for_carries = pass1_wait_for_carries;
		asm_data->pass1_done_with_carries = pass1_done_with_carries;
		asm_data->pass1_get_next_block = pass1_get_next_block;
		asm_data->pass2_wake_up_threads = pass2_wake_up_threads;
		asm_data->pass2_get_next_block = pass2_get_next_block;
		return (0);
	}

/* Sanity check num_threads argument */

	if (gwdata->num_threads > MAX_AUXILLARY_THREADS+1)
		gwdata->num_threads = MAX_AUXILLARY_THREADS+1;

/* Init mutexes and events used to control auxillary threads */

	gwmutex_init (&gwdata->thread_lock);
	gwevent_init (&gwdata->thread_work_to_do);
	gwevent_init (&gwdata->all_threads_done);
	gwevent_init (&gwdata->gwcarries_complete);
	for (i = 0; i < gwdata->num_threads; i++)
		gwevent_init (&gwdata->pass1_norm_events[i]);
	gwdata->num_active_threads = 0;
	gwevent_signal (&gwdata->all_threads_done);

/* Set ptrs to call back routines in structure used by assembly code */

	asm_data->pass1_wake_up_threads = pass1_wake_up_threads_mt;
	asm_data->pass1_wait_for_carries = pass1_wait_for_carries_mt;
	asm_data->pass1_done_with_carries = pass1_done_with_carries_mt;
	asm_data->pass1_get_next_block = pass1_get_next_block_mt;
	asm_data->pass2_wake_up_threads = pass2_wake_up_threads_mt;
	asm_data->pass2_get_next_block = pass2_get_next_block_mt;

/* Pre-create each auxillary thread used in multiplication code. */
/* We allocate the memory here so that error recovery is easier. */

	gwdata->threads_must_exit = FALSE;
	memset (gwdata->thread_id, 0, sizeof (gwdata->thread_id));
	for (i = 0; i < gwdata->num_threads - 1; i++) {
		struct thread_data *info;

		info = (struct thread_data *)
			malloc (sizeof (struct thread_data));
		if (info == NULL) return (GWERROR_MALLOC);
		info->gwdata = gwdata;
		info->thread_num = i+1;
		info->asm_data_alloc =
			aligned_malloc (sizeof (struct gwasm_data) +
					NEW_STACK_SIZE, 4096);
		if (info->asm_data_alloc == NULL) {
			free (info);
			return (GWERROR_MALLOC);
		}
		info->scratch_area =
			aligned_malloc (gwdata->SCRATCH_SIZE, 128);
		if (info->scratch_area == NULL) {
			aligned_free (info->asm_data_alloc);
			free (info);
			return (GWERROR_MALLOC);
		}
		gwthread_create_waitable (&gwdata->thread_id[i],
					  &auxillary_thread, info);
	}

/* Return success */

	return (0);
}

/* Perform cleanup required by multi-threaded operation */

void multithread_term (
	gwhandle *gwdata)
{
	unsigned long i;

/* Return if we aren't multithreading or multithreading was not initialized */

	if (gwdata->num_threads <= 1) return;
	if (gwdata->thread_lock == NULL) return;

/* Set termination variable and fire up auxillary threads */

	gwmutex_lock (&gwdata->thread_lock);
	gwdata->threads_must_exit = TRUE;
	gwmutex_unlock (&gwdata->thread_lock);
	gwevent_signal (&gwdata->thread_work_to_do);

/* Wait for all the threads to exit.  We must do this so */
/* that this thread can safely delete the gwdata structure */

	for (i = 0; i < gwdata->num_threads - 1; i++)
		if (gwdata->thread_id[i])
			gwthread_wait_for_exit (&gwdata->thread_id[i]);

/* Now free up the multi-thread resources */

	gwmutex_destroy (&gwdata->thread_lock);
	gwevent_destroy (&gwdata->thread_work_to_do);
	gwevent_destroy (&gwdata->all_threads_done);
	gwevent_destroy (&gwdata->gwcarries_complete);
	for (i = 0; i < gwdata->num_threads; i++)
		gwevent_destroy (&gwdata->pass1_norm_events[i]);
	gwdata->thread_lock = NULL;
}

/* Cleanup any memory allocated for multi-precision math */

void gwdone (
	gwhandle *gwdata)	/* Handle returned by gwsetup */
{
	unsigned int i;

	multithread_term (gwdata);

	term_ghandle (&gwdata->gdata);
	if (gwdata->asm_data != NULL) {
		aligned_free ((char *) gwdata->asm_data - NEW_STACK_SIZE);
		gwdata->asm_data = NULL;
	}
	free (gwdata->dd_data);
	gwdata->dd_data = NULL;
	free (gwdata->gwnum_free);
	gwdata->gwnum_free = NULL;
	if (gwdata->gwnum_alloc != NULL) {
		for (i = 0; i < gwdata->gwnum_alloc_count; i++) {
			char	*p;
			int32_t	freeable;
			p = (char *) gwdata->gwnum_alloc[i];
			freeable = * (int32_t *) (p - 32) & ~GWFREED_TEMPORARILY;
			if (freeable) aligned_free ((char *) p - GW_HEADER_SIZE);
		}
		free (gwdata->gwnum_alloc);
		gwdata->gwnum_alloc = NULL;
	}
	free (gwdata->GW_MODULUS);
	gwdata->GW_MODULUS = NULL;
	if (gwdata->large_pages_ptr != NULL) {
		large_pages_free (gwdata->large_pages_ptr);
		gwdata->large_pages_ptr = NULL;
	} else
		aligned_free (gwdata->gwnum_memory);
	gwdata->gwnum_memory = NULL;
}

/* Routine to allocate aligned memory for our big numbers */
/* Memory is allocated on 128-byte boundaries, with an additional */
/* 32 bytes prior to the data for storing useful stuff */

gwnum gwalloc (
	gwhandle *gwdata)
{
	unsigned long size, aligned_size;
	char	*p, *q;
	int32_t	freeable;

/* Return cached gwnum if possible */

	if (gwdata->gwnum_free_count)
		return (gwdata->gwnum_free[--gwdata->gwnum_free_count]);

/* Allocate arrays if necessary */

	if (gwdata->gwnum_alloc == NULL) {
		gwdata->gwnum_free = (gwnum *)
			malloc (gwdata->gwnum_alloc_array_size * sizeof (gwnum));
		if (gwdata->gwnum_free == NULL) return (NULL);
		gwdata->gwnum_alloc = (gwnum *)
			malloc (gwdata->gwnum_alloc_array_size * sizeof (gwnum));
		if (gwdata->gwnum_alloc == NULL) return (NULL);
	} else if (gwdata->gwnum_alloc_count == gwdata->gwnum_alloc_array_size) {
		gwdata->gwnum_alloc_array_size += gwdata->gwnum_alloc_array_size >> 1;
		gwdata->gwnum_free = (gwnum *)
			realloc (gwdata->gwnum_free,
				 gwdata->gwnum_alloc_array_size * sizeof (gwnum));
		if (gwdata->gwnum_free == NULL) return (NULL);
		gwdata->gwnum_alloc = (gwnum *)
			realloc (gwdata->gwnum_alloc,
				 gwdata->gwnum_alloc_array_size * sizeof (gwnum));
		if (gwdata->gwnum_alloc == NULL) return (NULL);
	}

/* Use addr function on the last FFT value to compute the size. */
/* Allocate 96 extra bytes for header information and align the data */
/* appropriately.  When allocating memory out of the big buffer for */
/* the torture test, then only allocate data on 128 byte boundaries */
/* to maximize the number of gwnums allocated. */

	size = gwnum_datasize (gwdata);
	aligned_size = (size + GW_HEADER_SIZE + 127) & ~127;
	if (gwdata->large_pages_gwnum != NULL) {
		p = (char *) gwdata->large_pages_gwnum;
		gwdata->large_pages_gwnum = NULL;
		p += -GW_HEADER_SIZE & 127;
		freeable = 0;
	} else if (gwdata->GW_BIGBUF_SIZE >= size + aligned_size) {
		p = gwdata->GW_BIGBUF;
		gwdata->GW_BIGBUF += aligned_size;
		gwdata->GW_BIGBUF_SIZE -= aligned_size;
		p += -GW_HEADER_SIZE & 127;
		freeable = 0;
	} else {
		p = (char *) aligned_offset_malloc (
				size + GW_HEADER_SIZE, gwdata->GW_ALIGNMENT,
				(GW_HEADER_SIZE - gwdata->GW_ALIGNMENT_MOD) &
					(gwdata->GW_ALIGNMENT - 1));
		if (p == NULL) return (NULL);
		freeable = 1;
	}

/* Do a seemingly pointless memset!  This actually is very important. */
/* The memset will walk through the allocated memory sequentially, which */
/* increases the likelihood that contiguous virtual memory will map to */
/* contiguous physical memory.  The FFTs, especially the larger ones, */
/* optimizes L2 cache line collisions on the assumption that the FFT data */
/* is in contiguous physical memory.  Failure to do this results in as */
/* much as a 30% performance hit in an SSE2 2M FFT. */

	q = p + GW_HEADER_SIZE;
	memset (q, 0, size);

/* Initialize the header */

	* (uint32_t *) (q - 8) = size;	/* Size in bytes */
	* (uint32_t *) (q - 4) = 1;	/* Unnormalized adds count */
	* (uint32_t *) (q - 28) = 0;	/* Has-been-pre-ffted flag */
	* (int32_t *) (q - 32) = freeable; /* Mem should be freed flag */
	* (double *) (q - 16) = 0.0;
	* (double *) (q - 24) = 0.0;

/* Save pointer for easier cleanup */

	gwdata->gwnum_alloc[gwdata->gwnum_alloc_count++] = (gwnum) q;

/* Return the gwnum */

	return ((gwnum) q);
}

/* Free one of our special numbers */

void gwfree (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	q)		/* Number to free */
{
	if (gwdata->gwnum_free != NULL && q != NULL)
		gwdata->gwnum_free[gwdata->gwnum_free_count++] = q;
}

/* Typeless gwalloc and gwfree routines for giants code to call */

void *gwgiantalloc (
	void	*gwdata)
{
	return ((void *) gwalloc ((gwhandle *) gwdata));
}
void gwgiantfree (
	void	*gwdata,
	void	*q)
{
	gwfree ((gwhandle *) gwdata, (gwnum) q);
}

/* Specialized routine that allows giants code to deallocate one */
/* cached gwnum to free up space for allocating FFT giant's sincos data. */

void gwgiantdealloc (
	void	*gwdata_arg)
{
	gwhandle *gwdata;
	gwnum	p;
	int32_t	freeable;
	unsigned long i, j;

	gwdata = (gwhandle *) gwdata_arg;
	for (i = 0; i < gwdata->gwnum_free_count; i++) {
		p = gwdata->gwnum_free[i];
		freeable = * (int32_t *) ((char *) p - 32);
		if (freeable & GWFREED_TEMPORARILY) continue;
		if (!freeable) continue;

		for (j = 0; j < gwdata->gwnum_alloc_count; j++) {
			if (gwdata->gwnum_alloc[j] != p) continue;

			aligned_free ((char *) p - GW_HEADER_SIZE);
			gwdata->gwnum_free[i] = gwdata->gwnum_free[--gwdata->gwnum_free_count];
			gwdata->gwnum_alloc[j] = gwdata->gwnum_alloc[--gwdata->gwnum_alloc_count];
			return;
		}
	}
}


/* Specialized routines that let the giants code share the free */
/* memory pool used by gwnums. */

void gwfree_temporarily (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	q)
{
	* (int32_t *) ((char *) q - 32) |= GWFREED_TEMPORARILY;
	gwfree (gwdata, q);
}
void gwrealloc_temporarily (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	q)
{
	unsigned long i, j;

	ASSERTG (* (int32_t *) ((char *) q - 32) & GWFREED_TEMPORARILY);

	* (int32_t *) ((char *) q - 32) &= ~GWFREED_TEMPORARILY;

	for (i = j = 0; i < gwdata->gwnum_free_count; i++)
		if (gwdata->gwnum_free[i] != q)
			gwdata->gwnum_free[j++] = gwdata->gwnum_free[i];
	gwdata->gwnum_free_count = j;
}

/* Free all user allocated numbers */

void gwfreeall (
	gwhandle *gwdata)	/* Handle initialized by gwsetup */
{
	unsigned int i;
	char	*p;
	int32_t	freeable;

/* Go through the allocated list and free any user allocated gwnums that */
/* are freeable.  In other words, unless the user is using the BIGBUF */
/* kludge, free all possible memory. */

	gwdata->gwnum_free_count = 0;
	for (i = 0; i < gwdata->gwnum_alloc_count; i++) {
		if (gwdata->gwnum_alloc[i] == gwdata->GW_MODULUS_FFT) continue;
		if (gwdata->gwnum_alloc[i] == gwdata->GW_RECIP_FFT) continue;
		if (gwdata->gwnum_alloc[i] == gwdata->GW_RANDOM) continue;
		p = (char *) gwdata->gwnum_alloc[i];
		freeable = * (int32_t *) (p - 32) & ~GWFREED_TEMPORARILY;
		if (freeable) {
			aligned_free ((char *) p - GW_HEADER_SIZE);
			gwdata->gwnum_alloc[i--] = gwdata->gwnum_alloc[--gwdata->gwnum_alloc_count];
		}
		else
			gwdata->gwnum_free[gwdata->gwnum_free_count++] = gwdata->gwnum_alloc[i];
	}
}

/* Copy one gwnum to another gwnum */

void gwcopy (			/* Copy a gwnum */
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	s,		/* Source */
	gwnum	d)		/* Dest */
{
	uint32_t free_offset;

/* Load the one piece of information that should not be copied over */

	free_offset = ((uint32_t *) d)[-8];

/* Copy the data and 96-byte header */

	memcpy ((char *) d - GW_HEADER_SIZE,
		(char *) s - GW_HEADER_SIZE,
		((uint32_t *) s)[-2] + GW_HEADER_SIZE);

/* Restore the one piece of information that should not be copied over */

	((uint32_t *) d)[-8] = free_offset;
}

/* To optimize use of the L1 cache we scramble the FFT data. */
/* Consult the assembly language code for better descriptions of this */
/* shuffling process.  This C code must accurately reflect the shuffling */
/* the assembly language code is expecting. */

unsigned long addr_offset (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	unsigned long i)
{
	unsigned long fftlen;
	unsigned long addr, i1, i2, i3, i6;

	fftlen = gwdata->FFTLEN;

/* P4 uses a different memory layout - more suitable to SSE2 */

	if (gwdata->cpu_flags & CPU_SSE2) {
		unsigned long sets, pfa, temp;

/* Small FFTs use one pass, not very convoluted.  This the example for	*/
/* a length 2048 FFT:							*/
/*	0	512	1	513	1024	1536	1025	1537	*/
/*	2	...							*/
/*	...								*/
/*	510								*/
/* PFA-style FFTs are a little tricker.  See assembly code for example.	*/

		if (gwdata->PASS2_LEVELS == 0) {
			sets = fftlen >> 3;
			if (i >= (fftlen >> 1)) {
				i6 = 1;
				i -= (fftlen >> 1);
			} else
				i6 = 0;
			i1 = i & 1; i >>= 1;
			i3 = 0;
			for (pfa = sets; pfa > 8; pfa >>= 1);
			if (pfa == 5) {
				temp = sets / 5;
				if (i < temp * 2) {
					sets = temp;
				} else {
					i3 = temp; i -= temp * 2;
					sets = temp * 4;
				}
			} else if (pfa == 7) {
				temp = sets / 7;
				if (i < temp * 2) {
					sets = temp;
				} else if (i < temp * 6) {
					i3 = temp; i -= temp * 2;
					sets = temp * 2;
				} else {
					i3 = temp * 3; i -= temp * 6;
					sets = temp * 4;
				}
			}
			i3 += i % sets; i /= sets;
			addr = (((((i3 << 1) + i6) << 1) + i1) << 1) + i;
			addr = addr * sizeof (double);
		}

/* Larger FFTs use two passes.  This the example for a length 64K FFT:	*/
/*	0	1K	16K	17K	32K	33K	48K	49K	*/
/*	1	...							*/
/*	...								*/
/*	1023	...							*/
/*	2K	...							*/
/*	...								*/

		else {
			sets = fftlen >> (gwdata->PASS2_LEVELS + 2);
			if (i >= (fftlen >> 1)) {
				i6 = 1;
				i -= (fftlen >> 1);
			} else
				i6 = 0;
			i1 = i & ((1 << (gwdata->PASS2_LEVELS - 1)) - 1);
			i >>= (gwdata->PASS2_LEVELS - 1);
			i2 = i & 1; i >>= 1;
			i3 = 0;
			for (pfa = sets; pfa > 8; pfa >>= 1);
			if (pfa == 5) {
				temp = sets / 5;
				if (i < temp * 2) {
					sets = temp;
				} else {
					i3 = temp; i -= temp * 2;
					sets = temp * 4;
				}
			} else if (pfa == 7) {
				temp = sets / 7;
				if (i < temp * 2) {
					sets = temp;
				} else if (i < temp * 6) {
					i3 = temp; i -= temp * 2;
					sets = temp * 2;
				} else {
					i3 = temp * 3; i -= temp * 6;
					sets = temp * 4;
				}
			}
			i3 += i % sets; i /= sets;
			addr = i3 * (1 << (gwdata->PASS2_LEVELS - 1));
			addr = ((((((addr + i1) << 1) + i6) << 1) + i) << 1) + i2;
			addr = addr * sizeof (double);
			/* Now add 128 bytes every 8KB and one pass2gapsize */
			/* for every pass 2 block. */
			addr = addr + (addr >> 13) * 128 + i3 * gwdata->PASS2GAPSIZE;
		}
	}

/* One pass x87 FFTs use a near flat memory model. */

	else if (gwdata->PASS2_LEVELS == 0) {
		if (i >= (fftlen >> 1)) {
			i2 = 1;
			i -= (fftlen >> 1);
		} else
			i2 = 0;
		addr = i * 16 + i2 * 8;
	}

/* Two pass x87 FFTs use a near flat memory model.  Waste 64 bytes */
/* between 4KB.  Waste 64 bytes between every block (4KB, 16KB, or 64KB). */

	else {
		if (i >= (fftlen >> 1)) {
			i2 = 1;
			i -= (fftlen >> 1);
		} else
			i2 = 0;
		addr = i * 16 + i2 * 8 + (i >> 8) * 64 + (i >> gwdata->PASS2_LEVELS) * 64;
	}

/* Return the offset */

	return (addr);
}

/* Return the address of ith element in the FFT array */

double *addr (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	g,
	unsigned long i)
{
	return ((double *) ((char *) g + addr_offset (gwdata, i)));
}

/* Return the amount of data allocated by gwsetup */

unsigned long gwmemused (
	gwhandle *gwdata)	/* Handle initialized by gwsetup */
{
	return (gwdata->jmptab->mem_needed + gwdata->jmptab->scratch_size);
}

/* Get the amount of memory required for the gwnum's raw FFT data.  This */
/* does not include the GW_HEADER_SIZE bytes for the header or any pad */
/* bytes that might be allocated for alignment.  I see little need for */
/* programs to use this routine. */

unsigned long gwnum_datasize (
	gwhandle *gwdata)	/* Handle initialized by gwsetup */
{
	return (addr_offset (gwdata, gwdata->FFTLEN - 1) + sizeof (double));
}

/* Get the amount of memory likely to be allocated a gwnum.  This includes */
/* FFT data, headers, and pad bytes for alignment. */

unsigned long gwnum_size (
	gwhandle *gwdata)	/* Handle initialized by gwsetup */
{
	unsigned long mem;
	mem = gwnum_datasize (gwdata) + GW_HEADER_SIZE;
	mem += gwdata->GW_ALIGNMENT-1;
	return (mem - mem % (gwdata->GW_ALIGNMENT));
}

/* Each FFT word is multiplied by a two-to-phi value.  These */
/* routines set and get the FFT value without the two-to-phi */
/* multiplier. */

int get_fft_value (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	g,
	unsigned long i,
	long	*retval)
{
	double	val;

/* Get the FFT data and validate it */

	val = * addr (gwdata, g, i);
	if (! is_valid_double (val)) return (GWERROR_BAD_FFT_DATA);

/* Handle the rational FFT case quickly */

	if (gwdata->RATIONAL_FFT) {
		*retval = (long) val;
	}

/* Multiply by two-to-minus-phi to generate an integer. */

	else {
		val = val * gwfft_weight_inverse_sloppy (gwdata->dd_data, i);
		if (val < -0.5)
			*retval = (long) (val - 0.5);
		else
			*retval = (long) (val + 0.5);
	}

/* Return success */

	return (0);
}

void set_fft_value (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	g,
	unsigned long i,
	long	val)
{

/* Handle the rational FFT case quickly */

	if (gwdata->RATIONAL_FFT || val == 0) {
		* addr (gwdata, g, i) = val;
		return;
	}

/* Multiply by two-to-phi to generate the proper double. */

	* addr (gwdata, g, i) = val * gwfft_weight_sloppy (gwdata->dd_data, i);
}

/* Some words in the FFT data contain floor(p/N), some words contain */
/* floor(p/N)+1 bits.  This function returns TRUE in the latter case. */

int is_big_word (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	unsigned long i)
{
	unsigned long base, next_base;

/* Compute the number of b in this word.  It is a big word if */
/* the number of b is more than NUM_B_PER_SMALL_WORD. */

	base = gwfft_base (gwdata->dd_data, i);
	next_base = gwfft_base (gwdata->dd_data, i+1);
	return ((next_base - base) > gwdata->NUM_B_PER_SMALL_WORD);
}

/* Routine map a "bit" number into an FFT word and "bit" within that word */
/* If b != 2, this routine locates the nth b amongst the FFT words. */

void bitaddr (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	unsigned long bit,
	unsigned long *word,
	unsigned long *bit_in_word)
{

/* What word is the bit in? */

	*word = (unsigned long) ((double) bit / gwdata->avg_num_b_per_word);
	if (*word >= gwdata->FFTLEN) *word = gwdata->FFTLEN - 1;

/* Compute the bit within the word. */

	*bit_in_word = bit - gwfft_base (gwdata->dd_data, *word);
}

/* Return a description of the FFT type chosen */

void gwfft_description (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	char	*buf)		/* Buffer to return string in */
{
	sprintf (buf, "%s%sFFT length %lu%s",
		 gwdata->ALL_COMPLEX_FFT ? "all-complex " :
		 gwdata->ZERO_PADDED_FFT ? "zero-padded " :
		 gwdata->GENERAL_MOD ? "generic reduction " : "",
		 (gwdata->cpu_flags & CPU_SSE2) ? "" : "x87 ",
		 gwdata->FFTLEN > 4194304 ?	gwdata->FFTLEN / 1048576 :
		 gwdata->FFTLEN >= 4096 ?	gwdata->FFTLEN / 1024 :
						gwdata->FFTLEN,
		 gwdata->FFTLEN > 4194304 ?	"M" :
		 gwdata->FFTLEN >= 4096 ?	"K" :
						"");
	if (gwdata->num_threads > 1)
		sprintf (buf + strlen (buf), ", %d threads",
			 gwdata->num_threads);

	if (gw_using_large_pages (gwdata))
		strcat (buf, " using large pages");
}

/* Return a string representation of a k/b/n/c combination */

void gw_as_string (
	char	*buf,		/* Buffer to return string in */
	double	k,		/* K in K*B^N+C */
	unsigned long b,	/* B in K*B^N+C */
	unsigned long n,	/* N in K*B^N+C */
	signed long c)		/* C in K*B^N+C */
{
	if (k != 1.0)
		sprintf (buf, "%.0f*%lu^%lu%c%lu", k, b, n,
			 c < 0 ? '-' : '+', abs (c));
	else if (b == 2 && c == -1)
		sprintf (buf, "M%lu", n);
	else {
		unsigned long cnt, temp_n;
		for (cnt = 0, temp_n = n; !(temp_n & 1); temp_n >>= 1, cnt++);
		if (b == 2 && temp_n == 1 && c == 1)
			sprintf (buf, "F%lu", cnt);
		else
			sprintf (buf, "%lu^%lu%c%lu", b, n,
				 c < 0 ? '-' : '+', abs (c));
	}
}

/* Get or clear the roundoff error.  Remember that if the roundoff error */
/* exceeds 0.5 then the FFT results will be wrong.  It is prudent to watch */
/* the roundoff error to make sure the roundoff error does not get close */
/* to 0.5. */

double gw_get_maxerr (
	gwhandle *gwdata)
{
	return (((struct gwasm_data *) gwdata->asm_data)->MAXERR);
}
void gw_clear_maxerr (
	gwhandle *gwdata)
{
	((struct gwasm_data *) gwdata->asm_data)->MAXERR = 0.0;
}

/* Return TRUE if we are operating near the limit of this FFT length */
/* Input argument is the percentage to consider as near the limit. */
/* For example, if percent is 1.0 and the FFT can handle 20 bits per word, */
/* then if there are more than 19.98 bits per word this function will */
/* return TRUE. */

int gwnear_fft_limit (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	double	pct)
{

/* Return TRUE if the virtual bits per word is near the maximum bits */
/* per word. */

	return (virtual_bits_per_word (gwdata) >
			(100.0 - pct) / 100.0 * gwdata->fft_max_bits_per_word);
}

/* Compute the virtual bits per word.  That is, the mersenne-mod-equivalent */
/* bits that this k,b,c combination uses.  This code must carefully invert */
/* the calculations gwinfo uses in determining whether a k,b,n,c combination */
/* will work for a given FFT size. */

double virtual_bits_per_word (
	gwhandle *gwdata)	/* Handle initialized by gwsetup */
{
	double	log2b, b_per_input_word, weighted_bits_per_output_word;
	double	max_weighted_bits_per_output_word;
	int	num_b_in_big_word, num_small_words, num_big_words;

	log2b = log2 (gwdata->b);

/* Compute our bits per output word exactly like gwinfo does for a zero padded FFT. */

	if (gwdata->ZERO_PADDED_FFT) {
		b_per_input_word = (double) (gwdata->n + gwdata->n) / gwdata->FFTLEN;
		num_b_in_big_word = (int) ceil (b_per_input_word);
		num_small_words = (int) ((num_b_in_big_word - b_per_input_word) * (gwdata->FFTLEN / 2 + 4));
		num_big_words = (gwdata->FFTLEN / 2 + 4) - num_small_words;
		max_weighted_bits_per_output_word =
			2.0 * gwdata->fft_max_bits_per_word + 0.6 * log2 (gwdata->FFTLEN / 2 + 4);
		weighted_bits_per_output_word =
		       2.0 * ((b_per_input_word + 1.0) * log2b - 1.0) +
		       0.6 * log2 (gwdata->FFTLEN / 2 + 4);
		if ((gwdata->n + gwdata->n) % gwdata->FFTLEN == 0)
			weighted_bits_per_output_word -= ((log2b <= 4.0) ? log2b : 1.4 * log2b);
		else if (! is_pathological_distribution (num_big_words, num_small_words))
			weighted_bits_per_output_word -=
				((log2b <= 3.0) ? (log2b - 1.0) / 2.0 :
				 (log2b <= 6.0) ? 1.0 + (log2b - 3.0) / 3.0 :
						  2.0 + (log2b - 6.0) / 6.0);
	}

/* Compute our bits per output word exactly like gwinfo does for a non-zero-padded FFT. */

	else {
		b_per_input_word = gwdata->avg_num_b_per_word;
		num_b_in_big_word = (int) ceil (b_per_input_word);
		num_small_words = (int) ((num_b_in_big_word - b_per_input_word) * gwdata->FFTLEN);
		num_big_words = gwdata->FFTLEN - num_small_words;
		max_weighted_bits_per_output_word =
			2.0 * gwdata->fft_max_bits_per_word + 0.6 * log2 (gwdata->FFTLEN);
		weighted_bits_per_output_word =
			2.0 * ((b_per_input_word + 1.0) * log2b - 1.0) +
			0.6 * log2 (gwdata->FFTLEN) +
			log2 (gwdata->k) + 1.7 * log2 (abs (gwdata->c));
		if (gwdata->k == 1.0 && gwdata->n % gwdata->FFTLEN == 0)
			weighted_bits_per_output_word -= ((log2b <= 4.0) ? log2b : 1.4 * log2b);
		else if (num_big_words == 1 && gwdata->k > 1.0)
			weighted_bits_per_output_word += log2b;
		else if (! is_pathological_distribution (num_big_words, num_small_words))
			weighted_bits_per_output_word -=
				((log2b <= 3.0) ? (log2b - 1.0) / 2.0 :
				 (log2b <= 6.0) ? 1.0 + (log2b - 3.0) / 3.0 :
						  2.0 + (log2b - 6.0) / 6.0);
	}

/* Now generate a value that can compared to gwdata->fft_max_bits_per_word */

	return (weighted_bits_per_output_word / max_weighted_bits_per_output_word * gwdata->fft_max_bits_per_word);
}

/* Given k,b,n,c determine the fft length.  If k,b,n,c is not supported */
/* then return zero. */

unsigned long gwmap_to_fftlen (
	double	k,		/* K in K*B^N+C. Must be a positive integer. */
	unsigned long b,	/* B in K*B^N+C. */
	unsigned long n,	/* N in K*B^N+C. Exponent to test. */
	signed long c)		/* C in K*B^N+C. Must be rel. prime to K. */
{
	gwhandle gwdata;	/* Temporary gwnum handle */

/* Get pointer to fft info and return the FFT length */

	gwinit (&gwdata);
	if (gwinfo (&gwdata, k, b, n, c)) return (0);
	return (gwdata.jmptab->fftlen);
}

/* Given an fft length, determine the maximum allowable exponent.  If fftlen */
/* is not supported then return zero. */

unsigned long gwmap_fftlen_to_max_exponent (
	unsigned long fftlen)
{
	gwhandle gwdata;	/* Temporary gwnum handle */

/* Get pointer to fft info and return the FFT length */

	gwinit (&gwdata);
	gwset_specific_fftlen (&gwdata, fftlen);
	if (gwinfo (&gwdata, 1.0, 2, 0, -1)) return (0);
	return (gwdata.jmptab->max_exp);
}

/* Given an fft length, determine how much memory is used for normalization */
/* and sin/cos tables.  If k,b,n,c is not supported, then kludgily return */
/* 100 million bytes used. */

unsigned long gwmap_to_memused (
	double	k,		/* K in K*B^N+C. Must be a positive integer. */
	unsigned long b,	/* B in K*B^N+C. */
	unsigned long n,	/* N in K*B^N+C. Exponent to test. */
	signed long c)		/* C in K*B^N+C. Must be rel. prime to K. */
{
	gwhandle gwdata;	/* Temporary gwnum handle */

/* Get pointer to fft info and return the memory used */

	gwinit (&gwdata);
	if (gwinfo (&gwdata, k, b, n, c)) return (100000000L);
	return (gwdata.jmptab->mem_needed + gwdata.jmptab->scratch_size);
}

/* Return the estimated size of a gwnum */

unsigned long gwmap_to_estimated_size (
	double	k,		/* K in K*B^N+C. Must be a positive integer. */
	unsigned long b,	/* B in K*B^N+C. */
	unsigned long n,	/* N in K*B^N+C. Exponent to test. */
	signed long c)		/* C in K*B^N+C. Must be rel. prime to K. */
{
	gwhandle gwdata;	/* Temporary gwnum handle */

/* Get pointer to fft info and return the memory used */

	gwinit (&gwdata);
	if (gwinfo (&gwdata, k, b, n, c)) return (100000000L);
	return (addr_offset (&gwdata, gwdata.FFTLEN - 1) + sizeof (double));
}

/* Speed of other x87 processors compared to a Pentium II */

#define REL_486_SPEED	8.4	/* 486 is over 8 times slower than PII */
#define REL_K6_SPEED	3.0	/* K6 is 3 times slower than PII */
#define REL_P3_SPEED	0.8	/* Pentium III is 20% faster than PII */
#define REL_K7_SPEED	0.6	/* Athlons are much faster than a PII */

/* Speed of other SSE2 processors compared to a Pentium 4 */

#define REL_AMD64_SPEED	1.1	/* AMD64 is slightly slower than a P4 */
#define REL_PM_SPEED	1.4	/* Pentium M, Core are much slower than a P4 */
#define REL_ATOM_SPEED	5.0	/* Atoms are much, much slower than a P4 */
#define REL_CORE2_SPEED	0.625	/* Core 2 is much faster than a P4 */
#define REL_I7_SPEED	0.59	/* Core i7 is even faster than a Core 2 */
#define REL_PHENOM_SPEED 0.67	/* AMD Phenom is faster that a P4 */

/* Make a guess as to how long a squaring will take.  If the number cannot */
/* be handled, then kludgily return 100.0. */

double gwmap_to_timing (
	double	k,		/* K in K*B^N+C. Must be a positive integer. */
	unsigned long b,	/* B in K*B^N+C. */
	unsigned long n,	/* N in K*B^N+C. Exponent to test. */
	signed long c)		/* C in K*B^N+C. Must be rel. prime to K. */
{
	gwhandle gwdata;	/* Temporary gwnum handle */
	double	timing;

/* Get pointer to fft info */

	gwinit (&gwdata);
	if (gwinfo (&gwdata, k, b, n, c)) return (100.0);

/* Use my PII-400 or P4-1400 timings as a guide. */

	timing = gwdata.jmptab->timing;

/* Since the program is about 10% memory bound, the program will not */
/* speed up linearly with increase in chip speed.  Note, no attempt is */
/* made to differentiate between memory bus speed - we're */
/* just returning an educated guess here. */

/* Adjust timing for various CPU architectures. */
/* For Intel, 486s were very slow.  Pentium, Pentium Pro, Pentium II,
/* and old celerons were slow because they did not support prefetch. */
/* AMD64s and Pentium Ms are slower than P4s. */

	if (gwdata.cpu_flags & CPU_SSE2) {
		timing = 0.10 * timing + 0.90 * timing * 1400.0 / CPU_SPEED;
		if (strstr (CPU_BRAND, "Phenom")) timing *= REL_PHENOM_SPEED;
		else if (strstr (CPU_BRAND, "AMD")) timing *= REL_AMD64_SPEED;
		else if (strstr (CPU_BRAND, "Atom")) timing *= REL_ATOM_SPEED;
		else if (strstr (CPU_BRAND, "Core 2")) timing *= REL_CORE2_SPEED;
		else if (strstr (CPU_BRAND, "Core(TM)2")) timing *= REL_CORE2_SPEED;
		else if (strstr (CPU_BRAND, "Core(TM) i7")) timing *= REL_I7_SPEED;
		else if (strstr (CPU_BRAND, "Pentium(R) M")) timing *= REL_PM_SPEED;
		else if (strstr (CPU_BRAND, "Core")) timing *= REL_PM_SPEED;
	} else {
		timing = 0.10 * timing + 0.90 * timing * 400.0 / CPU_SPEED;
		if (strstr (CPU_BRAND, "486")) timing *= REL_486_SPEED;
		else if (strstr (CPU_BRAND, "Intel")) {
			if (gwdata.cpu_flags & CPU_PREFETCH) timing *= REL_P3_SPEED;
		} else if (strstr (CPU_BRAND, "AMD")) {
			if (strstr (CPU_BRAND, "Unknown")) timing *= REL_486_SPEED;
			else if (strstr (CPU_BRAND, "K5")) timing *= REL_486_SPEED;
			else if (strstr (CPU_BRAND, "K6")) timing *= REL_K6_SPEED;
			else timing *= REL_K7_SPEED;
		} else
			timing *= REL_486_SPEED;
	}
	return (timing);
}

/* Given k,b,n,c determine the fft length and zero-padding state to be */
/* used.  Caller peers into the gwdata structure to get this info. */

int gwmap_to_fft_info (
	gwhandle *gwdata,	/* Uninitialized gwnum global data */
	double	k,		/* K in K*B^N+C. Must be a positive integer. */
	unsigned long b,	/* B in K*B^N+C. */
	unsigned long n,	/* N in K*B^N+C. Exponent to test. */
	signed long c)		/* C in K*B^N+C. Must be rel. prime to K. */
{
	gwinit (gwdata);
	return (gwinfo (gwdata, k, b, n, c));
}


/* Internal routine to help gwcopyzero */

void calc8ptrs (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	unsigned long n,
	uint32_t *ptrs)
{
	unsigned long i, j, k;

/* This is a grossly inefficient way to do this.  However, it should */
/* be called rarely. */

	for (i = 0; i < 8; i++) ptrs[i] = 0;
	for (i = 0; i < n; i++) {
		j = addr_offset (gwdata, i);
		k = (j & 63) >> 3;
		if (j >= ptrs[k]) ptrs[k] = j - (k << 3) + 64;
	}
}

/* Routine that sets up and calls assembly code to copy a gwnum from */
/* source to dest while zeroing some lower FFT words */

void gwcopyzero (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	s,
	gwnum	d,
	unsigned long n)
{
	struct gwasm_data *asm_data;

	ASSERTG (((int32_t *) s)[-7] == 0);	// Number not partially FFTed?

/* Handle case where no words are zeroed.  Some of the assembly routines */
/* do not like a word count of zero. */

	if (n == 0) {
		gwcopy (gwdata, s, d);
		return;
	}

/* Call assembly language copy-and-zeroing routine */

	asm_data = (struct gwasm_data *) gwdata->asm_data;
	asm_data->SRCARG = s;
	asm_data->DESTARG = d;
	asm_data->NUMARG = n;
	if ((gwdata->cpu_flags & CPU_SSE2) &&
	    (asm_data->COPYZERO[0] == 0 || n != gwdata->saved_copyz_n)) {
		gwdata->saved_copyz_n = n;
		calc8ptrs (gwdata, n, (uint32_t *) asm_data->COPYZERO);
	}
	gw_copyzero (gwdata, asm_data);

/* Copy the unnormalized add counter and clear the */
/* has been partially FFTed flag. */

	((int32_t *) d)[-1] = ((int32_t *) s)[-1];
	((int32_t *) d)[-7] = 0;
}

/* Set the constant which the results of a multiplication should be */
/* multiplied by.  Use this macro in conjunction with the c argument of */
/* gwsetnormroutine. */

void gwsetmulbyconst (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	long	val)
{
	struct gwasm_data *asm_data;
	double	ktimesval, big_word;

/* Save mulbyconst and -c * mulbyconst as a double */

	asm_data = (struct gwasm_data *) gwdata->asm_data;
	asm_data->XMM_MULCONST[0] =
	asm_data->XMM_MULCONST[1] = (double) val;
	asm_data->XMM_MINUS_C_TIMES_MULCONST[0] =
	asm_data->XMM_MINUS_C_TIMES_MULCONST[1] = (double) -gwdata->c * (double) val;

/* Split k*mulconst for zero-padded FFTs emulating modulo k*b^n+c.  Note we only need */
/* to split k*mulconst if k*mulconst * big_word exceeds what a floating point register can hold. */
/* 2^49 should give us enough room to handle a carry and still round properly. */	

	ktimesval = gwdata->k * (double) val;
	big_word = pow ((double) gwdata->b, gwdata->NUM_B_PER_SMALL_WORD + 1);
	if (gwdata->cpu_flags & CPU_SSE2 && ktimesval * big_word < 562949953421312.0) {
		asm_data->XMM_K_TIMES_MULCONST_HI[0] = asm_data->XMM_K_TIMES_MULCONST_HI[1] = 0.0;
		asm_data->XMM_K_TIMES_MULCONST_LO[0] = asm_data->XMM_K_TIMES_MULCONST_LO[1] = ktimesval;
	} else {
		asm_data->XMM_K_TIMES_MULCONST_HI[0] =
		asm_data->XMM_K_TIMES_MULCONST_HI[1] = floor (ktimesval / big_word) * big_word;
		asm_data->XMM_K_TIMES_MULCONST_LO[0] =
		asm_data->XMM_K_TIMES_MULCONST_LO[1] = ktimesval - asm_data->XMM_K_TIMES_MULCONST_HI[0];
	}
	big_word = big_word * big_word;
	asm_data->XMM_K_TIMES_MULCONST_HI_2[0] =
	asm_data->XMM_K_TIMES_MULCONST_HI_2[1] =
		floor (asm_data->XMM_K_TIMES_MULCONST_HI[0] / big_word) * big_word;
	asm_data->XMM_K_TIMES_MULCONST_HI_1[0] =
	asm_data->XMM_K_TIMES_MULCONST_HI_1[1] =
		asm_data->XMM_K_TIMES_MULCONST_HI[0] - asm_data->XMM_K_TIMES_MULCONST_HI_2[0];
}

/* Set the flag which controls whether the multiply code should begin the */
/* forward FFT of the results of a multiply. This is a little faster than */
/* than doing a full forward FFT later on.  The downside is the caller */
/* cannot convert the results of the multiply to an integer. */

void gwstartnextfft (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	int	state)		/* New value for POSTFFT */
{
	struct gwasm_data *asm_data;

	if (!gwdata->GENERAL_MOD) {
		asm_data = (struct gwasm_data *) gwdata->asm_data;
		asm_data->POSTFFT = state;
	}
}

/* Add a small constant at the specified power of b after the */
/* next multiplication.  That is, value*b^power_of_b is added to */
/* the next multiplication result.  This only works if k=1. */

void gwsetaddinatpowerofb (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	long	value,
	unsigned long power_of_b)
{
	unsigned long word, b_in_word;

	ASSERTG (gwdata->k == 1.0);

/* Tell assembly code to add the shifted value to the multiplication result. */

	bitaddr (gwdata, power_of_b, &word, &b_in_word);
	raw_gwsetaddin (gwdata, word, value * pow ((double) gwdata->b, b_in_word));
}

/* Routine that tells the assembly code to add a small value to the */
/* results of each multiply. */

void gwsetaddin (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	long	value)
{
	unsigned long word, b_in_word;

	ASSERTG (gwdata->k == 1.0 || abs (gwdata->c) == 1);

/* In a zero-padded FFT, the value is added into ZPAD0 */

	if (gwdata->ZERO_PADDED_FFT) {
		((struct gwasm_data *) gwdata->asm_data)->ADDIN_VALUE = (double) value;
		return;
	}

/* If k is 1, add the value into the first FFT word */

	if (gwdata->k == 1.0) {
		raw_gwsetaddin (gwdata, 0, value);
		return;
	}

/* If value is a multiple of b, "shift" it right and increment b count.  This */
/* will ensure that we modify the proper FFT word. */

	for (b_in_word = 0; value && value % (long) gwdata->b == 0; value = value / (long) gwdata->b)
		b_in_word++;

/* Convert the input value to 1/k format.  Case 1 (k*b^n-1): Inverse of k */
/* is b^n.  Case 3 (k*b^n+1): Inverse of k is -b^n.  No other cases can */
/* be handled. */

	if (gwdata->c == -1) {
		bitaddr (gwdata, gwdata->n + b_in_word, &word, &b_in_word);
	}
	else if (gwdata->c == 1) {
		bitaddr (gwdata, gwdata->n + b_in_word, &word, &b_in_word);
		value = -value;
	}

/* Tell assembly code to add the shifted value to the multiplication result. */

	raw_gwsetaddin (gwdata, word, value * pow ((double) gwdata->b, b_in_word));
}

/* Routine that tells the assembly code to add a small value to the */
/* results of each multiply */

void raw_gwsetaddin (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	unsigned long word,
	double	val)
{
	struct gwasm_data *asm_data;
	unsigned long row;

/* Compute the offset to the FFT data value */

	asm_data = (struct gwasm_data *) gwdata->asm_data;
	asm_data->ADDIN_OFFSET = addr_offset (gwdata, word);

/* If this is a two-pass SSE2 FFT, then we need to tell the assembly code */
/* the affected "row", that is which set of pass 1 data the add-in will */
/* take place */

	if (gwdata->cpu_flags & CPU_SSE2) {
		if (gwdata->PASS2_LEVELS == 0) {
			row = asm_data->ADDIN_OFFSET & 31;
			if (row == 8) asm_data->ADDIN_OFFSET += 8;
			if (row == 16) asm_data->ADDIN_OFFSET -= 8;
		}

/* Factor in the blkdst value in xfft3.mac to compute the two pass */
/* SSE2 addin_offset. */

		else {
			unsigned long num_rows;

			num_rows = (1 << (gwdata->PASS2_LEVELS - 1));
			row = (word & (num_rows - 1));
			asm_data->ADDIN_ROW =
					row & ~(gwdata->PASS1_CACHE_LINES - 1);
			asm_data->ADDIN_OFFSET -= (row >> 7) * 128 +
					(row / gwdata->PASS1_CACHE_LINES) *
					gwdata->PASS1_CACHE_LINES * 64;

/* This case is particularly nasty as we have to convert the FFT data offset */
/* into a scratch area offset.  In assembly language terms, this means */
/* subtracting out multiples of blkdst and adding in multiples of clmblkdst */
/* and clmblkdst8. */

			if (gwdata->SCRATCH_SIZE) {
				unsigned long blkdst;

				blkdst = addr_offset (gwdata, 1<<gwdata->PASS2_LEVELS);
				row = asm_data->ADDIN_OFFSET / blkdst;
				asm_data->ADDIN_OFFSET -= row * blkdst;
				asm_data->ADDIN_OFFSET +=
					row * (gwdata->PASS1_CACHE_LINES * 64) +
					(row >> 3) * 128;
			}
		}
	}

/* And now x87 FFTs also can use a scratch area.  Like the SSE2 code */
/* we have to convert the FFT data offsets for two-pass FFTs. */

	if (! (gwdata->cpu_flags & CPU_SSE2) && gwdata->PASS2_LEVELS) {
		unsigned long num_cache_lines, cache_line;

		num_cache_lines = (1 << (gwdata->PASS2_LEVELS - 1));
		cache_line = ((word >> 1) & (num_cache_lines - 1));

		asm_data->ADDIN_ROW = ((num_cache_lines>>7) - (cache_line>>7)) * 65536 +
			    (128 / gwdata->PASS1_CACHE_LINES -
			     (cache_line & 127) / gwdata->PASS1_CACHE_LINES) * 256;
		asm_data->ADDIN_OFFSET -= (cache_line >> 7) * 64 +
				(cache_line / gwdata->PASS1_CACHE_LINES) *
				gwdata->PASS1_CACHE_LINES * 32;

/* This case is particularly nasty as we have to convert the FFT data offset */
/* into a scratch area offset.  In assembly language terms, this means */
/* subtracting out multiples of blkdst and adding in multiples of clmblkdst */
/* and clmblkdst32. */

		if (gwdata->SCRATCH_SIZE) {
			unsigned long blkdst;

			blkdst = addr_offset (gwdata, 1 << gwdata->PASS2_LEVELS);
			row = asm_data->ADDIN_OFFSET / blkdst;
			asm_data->ADDIN_OFFSET -= row * blkdst;
			asm_data->ADDIN_OFFSET += row * (gwdata->PASS1_CACHE_LINES * 32);

/* Handle the FFTs where clmblkdst32 is used */

			if ((gwdata->FFTLEN >> (gwdata->PASS2_LEVELS+1)) >= 128)
				asm_data->ADDIN_OFFSET += (row >> 5) * 64;
		}
	}

/* Set the addin value - multiply it by two-to-phi and FFTLEN/2/k. */

	asm_data->ADDIN_VALUE = val *
			gwfft_weight_sloppy (gwdata->dd_data, word) *
			gwdata->FFTLEN * 0.5 / gwdata->k;
}


/********************************************************/
/* Routines to convert between gwnums and other formats */
/********************************************************/

/* Convert a double to a gwnum */

void dbltogw (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	double	d,		/* Input number */
	gwnum	g)		/* Gwnum value to set */
{
	stackgiant(tmp, 2);

	dbltog (d, tmp);
	gianttogw (gwdata, tmp, g);
}

/* Convert a binary value to a gwnum */

void binarytogw (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	const uint32_t *array,	/* Array containing the binary value */
	uint32_t arraylen,	/* Length of the array */
	gwnum	n)		/* Destination gwnum */
{
	giantstruct tmp;
	tmp.sign = arraylen;
	tmp.n = (uint32_t *) array;
	while (tmp.sign && tmp.n[tmp.sign-1] == 0) tmp.sign--;
	gianttogw (gwdata, &tmp, n);
}

/* Convert a binary value to a gwnum */

void binary64togw (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	const uint64_t *array,	/* Array containing the binary value */
	uint64_t arraylen,	/* Length of the array */
	gwnum	n)		/* Destination gwnum */
{
	giantstruct tmp;
	tmp.sign = (int) arraylen * sizeof (uint64_t) / sizeof (uint32_t);
	tmp.n = (uint32_t *) array;
	while (tmp.sign && tmp.n[tmp.sign-1] == 0) tmp.sign--;
	gianttogw (gwdata, &tmp, n);
}

/* Convert a binary value to a gwnum */

void binarylongstogw (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	const unsigned long *array, /* Array containing the binary value */
	unsigned long arraylen,	/* Length of the array */
	gwnum	n)		/* Destination gwnum */
{
	giantstruct tmp;
	tmp.sign = arraylen * sizeof (unsigned long) / sizeof (uint32_t);
	tmp.n = (uint32_t *) array;
	while (tmp.sign && tmp.n[tmp.sign-1] == 0) tmp.sign--;
	gianttogw (gwdata, &tmp, n);
}

/* Convert a giant to gwnum FFT format.  Giant must be a positive number. */

void gianttogw (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	giant	a,
	gwnum	g)
{
	giant	newg = NULL;
	unsigned long i, limit, carry;

	ASSERTG (a->sign >= 0);		/* We only handle positive numbers */

/* To make the mod k*b^n+c step faster, gwnum's are pre-multiplied by 1/k */
/* If k is greater than 1, then we calculate the inverse of k, multiply */
/* the giant by the inverse of k, and do a mod k*b^n+c. */

	if (gwdata->k > 1.0) {
		newg = popg (&gwdata->gdata, (((unsigned long) gwdata->bit_length >> 5) + 1) * 2);

		/* Easy case 1 (k*b^n-1): Inverse of k is b^n */

		if (gwdata->c == -1) {
			if (gwdata->b == 2) {
				gtog (a, newg);
				gshiftleft (gwdata->n, newg);
			} else {
				itog (gwdata->b, newg);
				power (newg, gwdata->n);
				mulgi (&gwdata->gdata, a, newg);
			}
		}

		/* Easy case 2 (k*b^n+1): Inverse of k is -b^n */

		else if (gwdata->c == 1) {
			if (gwdata->b == 2) {
				gtog (a, newg);
				gshiftleft (gwdata->n, newg);
				negg (newg);
			} else {
				itog (gwdata->b, newg);
				power (newg, gwdata->n);
				negg (newg);
				mulgi (&gwdata->gdata, a, newg);
			}
		}

		else {				/* General inverse case */
			giant	n;
			n = popg (&gwdata->gdata, ((unsigned long) gwdata->bit_length >> 5) + 1);
			ultog (gwdata->b, n);	/* Compute k*b^n+c */
			power (n, gwdata->n);
			dblmulg (gwdata->k, n);
			iaddg (gwdata->c, n);
			dbltog (gwdata->k, newg);	/* Compute 1/k */
			invg (n, newg);
			ASSERTG (newg->sign > 0);  /* Assert inverse found */
			mulgi (&gwdata->gdata, a, newg);
						/* Multiply input num by 1/k */
			pushg (&gwdata->gdata, 1);
		}

		specialmodg (gwdata, newg);
		a = newg;
	}

/* Figure out how many FFT words we will need to set */

	limit = (unsigned long) ceil ((double) bitlen (a) / (gwdata->avg_num_b_per_word * log2 (gwdata->b)));
	if (limit > gwdata->FFTLEN) limit = gwdata->FFTLEN;

/* Now convert the giant to FFT format.  For base 2 we simply copy bits.  */

	if (gwdata->b == 2) {
		unsigned long mask1, mask2, e1len;
		int	bits1, bits2, bits_in_next_binval;
		unsigned long binval;
		uint32_t *e1;

		e1len = a->sign;
		e1 = a->n;

		bits1 = gwdata->NUM_B_PER_SMALL_WORD;
		bits2 = bits1 + 1;
		mask1 = (1L << bits1) - 1;
		mask2 = (1L << bits2) - 1;
		if (e1len) {binval = *e1++; e1len--; bits_in_next_binval = 32;}
		else binval = 0;
		carry = 0;
		for (i = 0; i < limit; i++) {
			int	big_word, bits;
			long	value, mask;
			big_word = is_big_word (gwdata, i);
			bits = big_word ? bits2 : bits1;
			mask = big_word ? mask2 : mask1;
			if (i == limit - 1) value = binval;
			else value = binval & mask;
			value = value + carry;
			if (value > (mask >> 1) && bits > 1 && i != gwdata->FFTLEN - 1) {
				value = value - (mask + 1);
				carry = 1;
			} else {
				carry = 0;
			}
			set_fft_value (gwdata, g, i, value);

			binval >>= bits;
			if (e1len == 0) continue;
			if (bits_in_next_binval < bits) {
				if (bits_in_next_binval)
					binval |= (*e1 >> (32 - bits_in_next_binval)) << (32 - bits);
				bits -= bits_in_next_binval;
				e1++; e1len--; bits_in_next_binval = 32;
				if (e1len == 0) continue;
			}
			if (bits) {
				binval |= (*e1 >> (32 - bits_in_next_binval)) << (32 - bits);
				bits_in_next_binval -= bits;
			}
		}
	}

/* Otherwise (non-base 2), we do a recursive divide and conquer radix conversion. */
/* The resursive routine writes on a, so make a copy before calling */

	else {
		if (a != newg) {
			newg = popg (&gwdata->gdata, a->sign * 2);
			gtog (a, newg);
			a = newg;
		}
		carry = nonbase2_gianttogw (gwdata, a, g, limit, 0, 0);
	}

/* Write carry, if any, to FFT data */

	if (carry) set_fft_value (gwdata, g, limit++, carry);

/* Clear the upper words */

	for (i = limit; i < gwdata->FFTLEN; i++)
		set_fft_value (gwdata, g, i, 0);

/* Clear various flags */

	((int32_t *) g)[-1] = 1; /* Set unnormalized add counter */
	((int32_t *) g)[-7] = 0; /* Clear has been partially FFTed flag */

/* Free allocated memory */

	if (a == newg) pushg (&gwdata->gdata, 1);
}

/* Internal recursive routine to convert a giant to gwnum FFT format. */

long nonbase2_gianttogw (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	giant	a,
	gwnum	g,
	unsigned long limit,	/* How many FFT words to set */
	unsigned long offset,	/* Offset into FFT array of words to set */
	long	carry)		/* Carry to add into this section */		 
{
	ASSERTG (a->sign >= 0);		/* We only handle positive numbers */

/* If we are converting a lot of words, divide and conquer. */

	if (limit >= 50) {
		giant	upper, tmp;
		int	num_b;
		unsigned long half_limit = limit >> 1;

		tmp = popg (&gwdata->gdata, limit * 2);
		upper = popg (&gwdata->gdata, a->sign * 2);
		num_b = gwfft_base (gwdata->dd_data, offset + half_limit) - gwfft_base (gwdata->dd_data, offset);
		itog (gwdata->b, tmp);
		power (tmp, num_b);
		gtog (a, upper);
		divg (tmp, upper);
		mulgi (&gwdata->gdata, upper, tmp);
		subg (tmp, a);
		carry = nonbase2_gianttogw (gwdata, a, g, half_limit, offset, carry);
		carry = nonbase2_gianttogw (gwdata, upper, g, limit - half_limit, offset + half_limit, carry);
		pushg (&gwdata->gdata, 2);
	}

/* Convert the giant to FFT format */

	else {
		giant	newg, tmp;
		unsigned long i, mask1, mask2;
		long	value;

		newg = popg (&gwdata->gdata, a->sign * 2);
		tmp = popg (&gwdata->gdata, a->sign * 2);

		mask1 = intpow (gwdata->b, gwdata->NUM_B_PER_SMALL_WORD);
		mask2 = gwdata->b * mask1;
		for (i = offset; i < offset + limit; i++) {
			unsigned long mask;

			mask = is_big_word (gwdata, i) ? mask2 : mask1;

			gtog (a, newg);
			ultog (mask, tmp);
			divg (tmp, a);
			mulgi (&gwdata->gdata, a, tmp);
			subg (tmp, newg);
			value = (newg->sign) ? newg->n[0] : 0;
			value += carry;

			if (value > (long) (mask >> 1) && i != gwdata->FFTLEN - 1) {
				value = value - mask;
				carry = 1;
			} else {
				carry = 0;
			}
			set_fft_value (gwdata, g, i, value);
		}
		pushg (&gwdata->gdata, 2);
	}

/* Return carry for next section */

	return (carry);
}


/* Convert a gwnum to a binary value.  Returns the number of 32-bit values */
/* written to the array.  The array is NOT zero-padded.  Returns a */
/* negative number if an error occurs during the conversion.  An error */
/* can happen if the FFT data contains a NaN or infinity value. */

long gwtobinary (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	n,		/* Source gwnum */
	uint32_t *array,	/* Array to contain the binary value */
	uint32_t arraylen)	/* Maximum size of the array */
{
	giant	tmp;
	int	err_code;

	tmp = popg (&gwdata->gdata, ((unsigned long) gwdata->bit_length >> 5) + 5);
	err_code = gwtogiant (gwdata, n, tmp);
	if (err_code < 0) return (err_code);
	ASSERTG ((unsigned long) tmp->sign <= arraylen);
	memcpy (array, tmp->n, tmp->sign * sizeof (uint32_t));
	pushg (&gwdata->gdata, 1);
	return (tmp->sign);
}

/* Convert a gwnum to a binary value.  Returns the number of 64-bit values */
/* written to the array.  The array is NOT zero-padded.  Returns a */
/* negative number if an error occurs during the conversion.  An error */
/* can happen if the FFT data contains a NaN or infinity value. */

long gwtobinary64 (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	n,		/* Source gwnum */
	uint64_t *array,	/* Array to contain the binary value */
	uint64_t arraylen)	/* Maximum size of the array */
{
	giant	tmp;
	int	err_code;

	tmp = popg (&gwdata->gdata, ((unsigned long) gwdata->bit_length >> 5) + 5);
	err_code = gwtogiant (gwdata, n, tmp);
	if (err_code < 0) return (err_code);
	tmp->n[tmp->sign] = 0;
	tmp->sign = (tmp->sign + 1) / 2;
	ASSERTG ((unsigned long) tmp->sign <= arraylen);
	memcpy (array, tmp->n, tmp->sign * sizeof (unsigned long));
	pushg (&gwdata->gdata, 1);
	return (tmp->sign);
}

/* Convert a gwnum to a binary value.  Returns the number of longs */
/* written to the array.  The array is NOT zero-padded.  Returns a */
/* negative number if an error occurs during the conversion.  An error */
/* can happen if the FFT data contains a NaN or infinity value. */

long gwtobinarylongs (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	n,		/* Source gwnum */
	unsigned long *array,	/* Array to contain the binary value */
	unsigned long arraylen)	/* Maximum size of the array */
{
	giant	tmp;
	int	err_code;

	tmp = popg (&gwdata->gdata, ((unsigned long) gwdata->bit_length >> 5) + 5);
	err_code = gwtogiant (gwdata, n, tmp);
	if (err_code < 0) return (err_code);
	if (sizeof (unsigned long) > sizeof (uint32_t)) {
		tmp->n[tmp->sign] = 0;
		tmp->sign = (tmp->sign + 1) / 2;
	}
	ASSERTG ((unsigned long) tmp->sign <= arraylen);
	memcpy (array, tmp->n, tmp->sign * sizeof (unsigned long));
	pushg (&gwdata->gdata, 1);
	return (tmp->sign);
}

/* Convert a gwnum to a giant.  WARNING: Caller must allocate an array that */
/* is several words larger than the maximum result that can be returned. */
/* This is a gross kludge that lets gwtogiant use the giant for intermediate */
/* calculations. */

int gwtogiant (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	gg,
	giant	v)
{
	int	err_code;
	unsigned long limit;

	ASSERTG (((uint32_t *) gg)[-7] == 0);	// Number not partially FFTed?

/* If this is a general-purpose mod, then only convert the needed words */
/* which will be less than half the FFT length.  If this is a zero padded */
/* FFT, then only convert a little more than half of the FFT data words. */
/* For a DWT, convert all the FFT data. */

	if (gwdata->GENERAL_MOD) limit = gwdata->GW_GEN_MOD_MAX + 4;
	else if (gwdata->ZERO_PADDED_FFT) limit = gwdata->FFTLEN / 2 + 4;
	else limit = gwdata->FFTLEN;

/* GENERAL_MOD has some strange cases we must handle.  In particular the */
/* last fft word translated can be 2^bits and the next word could be -1, */
/* this must be translated into zero, zero. */

	if (gwdata->GENERAL_MOD) {
		long	val, prev_val;
		err_code = get_fft_value (gwdata, gg, limit, &val);
		if (err_code) return (err_code);
		ASSERTG (val >= -1 && val <= 0);
		while (limit > 1) {		/* Find top word */
			err_code = get_fft_value (gwdata, gg, limit-1, &prev_val);
			if (err_code) return (err_code);
			if (val != prev_val || val < -1 || val > 0) break;
			limit--;
		}
		if (val > 0 || val < -1 || (val == -1 && prev_val >= 0)) {
			int	num_b;
			num_b = gwdata->NUM_B_PER_SMALL_WORD;
			if (is_big_word (gwdata, limit-1)) num_b++;
			prev_val += val * intpow (gwdata->b, num_b);
			set_fft_value (gwdata, gg, limit, (prev_val < 0) ? -1 : 0);
			set_fft_value (gwdata, gg, limit-1, prev_val);
		}
	}

/* If base is 2 we can simply copy the bits out of each FFT word */

	if (gwdata->b == 2) {
		long	val;
		int	j, bits, bitsout, carry;
		unsigned long i;
		uint32_t *outptr;

/* Collect bits until we have all of them */

		carry = 0;
		bitsout = 0;
		outptr = v->n;
		*outptr = 0;
		for (i = 0; i < limit; i++) {
			err_code = get_fft_value (gwdata, gg, i, &val);
			if (err_code) return (err_code);
			bits = gwdata->NUM_B_PER_SMALL_WORD;
			if (is_big_word (gwdata, i)) bits++;
			val += carry;
			for (j = 0; j < bits; j++) {
				*outptr >>= 1;
				if (val & 1) *outptr += 0x80000000;
				val >>= 1;
				bitsout++;
				if (bitsout == 32) {
					outptr++;
					bitsout = 0;
				}
			}
			carry = val;
		}

/* Finish outputting the last word and any carry data */

		while (bitsout || (carry != -1 && carry != 0)) {
			*outptr >>= 1;
			if (carry & 1) *outptr += 0x80000000;
			carry >>= 1;
			bitsout++;
			if (bitsout == 32) {
				outptr++;
				bitsout = 0;
			}
		}

/* Set the length */

		v->sign = (long) (outptr - v->n);
		while (v->sign && v->n[v->sign-1] == 0) v->sign--;

/* If carry is -1, the gwnum is negative.  Ugh.  Flip the bits and sign. */

		if (carry == -1) {
			for (j = 0; j < v->sign; j++) v->n[j] = ~v->n[j];
			while (v->sign && v->n[v->sign-1] == 0) v->sign--;
			iaddg (1, v);
			v->sign = -v->sign;
		}
	}

/* Otherwise (base is not 2) we must do a radix conversion */

	else {
		giantstruct *array = NULL;
		uint32_t *buf = NULL;
		giant	small_base = NULL;
		giant	large_base = NULL;
		unsigned long i, gap, small_size, last_small_size;

		array = (giantstruct *) malloc (limit * sizeof (giantstruct));
		if (array == NULL) { memerr: err_code =  GWERROR_MALLOC; goto err; }
		buf = (uint32_t *) malloc (limit * sizeof (uint32_t));
		if (buf == NULL) goto memerr;
		small_base = allocgiant (limit);
		if (small_base == NULL) goto memerr;
		large_base = allocgiant (limit);
		if (large_base == NULL) goto memerr;
		for (i = 0; i < limit; i++) {
			long	val;
			err_code = get_fft_value (gwdata, gg, i, &val);
			if (err_code) {
err:				free (array);
				free (buf);
				free (small_base);
				free (large_base);
				return (err_code);
			}
			array[i].n = &buf[i];
			setmaxsize(&array[i], limit);
			itog (val, &array[i]);
		}

/* Loop combining pairs into ever larger and larger numbers.  Do all but last combining pass. */

		gap = 1;
		while (gap + gap < limit) {
			small_size = gwfft_base (gwdata->dd_data, gap) - 1;
			if (gap == 1)
				itog (intpow (gwdata->b, small_size), small_base);
			else if (small_size == last_small_size * 2)
				squaregi (&gwdata->gdata, small_base);
			else
				mulgi (&gwdata->gdata, large_base, small_base);
			itog (gwdata->b, large_base);
			mulgi (&gwdata->gdata, small_base, large_base);
			for (i = 0; i + gap < limit; i += gap + gap) {
				gtog (&array[i+gap], v);
				if (gwfft_base (gwdata->dd_data, i+gap) - gwfft_base (gwdata->dd_data, i) == small_size)
					mulgi (&gwdata->gdata, small_base, v);
				else
					mulgi (&gwdata->gdata, large_base, v);
				addg (v, &array[i]);
			}
			gap = gap << 1;
			last_small_size = small_size;
		}

/* Do the last combining pass, outputting result directly to v. */

		if (gwfft_base (gwdata->dd_data, gap) == small_size * 2 + 1)
			mulgi (&gwdata->gdata, small_base, large_base);
		else
			squaregi (&gwdata->gdata, large_base);
		gtog (&array[gap], v);
		mulgi (&gwdata->gdata, large_base, v);
		addg (&array[0], v);

/* Clean up */

		free (array);
		free (buf);
		free (small_base);
		free (large_base);
	}

/* Since all gwnums are premultiplied by the inverse of k, we must now */
/* multiply by k to get the true result. */

	if (gwdata->k > 1.0) {
		stackgiant(k,2);
		dbltog (gwdata->k, k);
		mulgi (&gwdata->gdata, k, v);
	}

/* The gwnum is not guaranteed to be smaller than k*b^n+c.  Handle this */
/* possibility.  This also converts negative values to positive. */

	specialmodg (gwdata, v);

/* Return success */

	return (0);
}

/* Special modg.  This is a fast implementation of mod k*2^n+c using just */
/* shifts, adds, and divide and mul by small numbers.  All others moduli */
/* call the slow giants code. */

void specialmodg (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	giant	g)
{
	int	neg, count;
	giant	n;

/* If the modulus is a general-purpose number, then let the giants code */
/* do the work.  This is done for both GENERAL_MOD and (k*b^n+c)/d cases. */

	if (gwdata->GW_MODULUS != NULL) {
		modgi (&gwdata->gdata, gwdata->GW_MODULUS, g);
		return;
	}

/* Calculate the modulo number - k*b^n+c */

	n = popg (&gwdata->gdata, ((unsigned long) gwdata->bit_length >> 5) + 1);
	ultog (gwdata->b, n);
	power (n, gwdata->n);
	dblmulg (gwdata->k, n);
	iaddg (gwdata->c, n);

/* If b is not 2 let the giants code do the work. */

	if (gwdata->b != 2) {
		modgi (&gwdata->gdata, n, g);
		pushg (&gwdata->gdata, 1);
		return;
	}

/* Do the quick modulus code twice because in the case where */
/* abs(c) > k once won't get us close enough. */

	neg = FALSE;
	for (count = 0; count < 2; count++) {

/* Handle negative input values */

	    neg ^= (g->sign < 0);
	    g->sign = abs (g->sign);

/* If number is bigger than the modulus, do a mod using shifts and adds */
/* This will get us close to the right answer. */

	    if (gcompg (g, n) > 0) {
		giant	tmp1;

/* Allocate temporary */

		tmp1 = popg (&gwdata->gdata, ((unsigned long) gwdata->bit_length >> 5) + 5);

/* Calculate the modulo by dividing the upper bits of k, multiplying by */
/* c and subtracting that from the bottom bits. */

		gtogshiftright (gwdata->n, g, tmp1);	// Upper bits
		gmaskbits (gwdata->n, g);		// Lower bits

		if (gwdata->k == 1.0) {
			imulg (gwdata->c, tmp1);	// Upper bits times C
			subg (tmp1, g);
		} else {
			giant	tmp2, tmp3;

			tmp2 = popg (&gwdata->gdata, (((unsigned long) gwdata->bit_length >> 5) + 5) * 2);
			tmp3 = popg (&gwdata->gdata, ((unsigned long) gwdata->bit_length >> 5) + 5);

			gtog (tmp1, tmp2);
			dbltog (gwdata->k, tmp3);
			divgi (&gwdata->gdata, tmp3, tmp1);	// Upper bits over K
			mulgi (&gwdata->gdata, tmp1, tmp3);
			subg (tmp3, tmp2);	// Upper bits mod K

			gshiftleft (gwdata->n, tmp2);
			addg (tmp2, g);		// Upper bits mod K+lower bits

			imulg (gwdata->c, tmp1);// Upper bits over K times C
			subg (tmp1, g);
			pushg (&gwdata->gdata, 2);
		}

		pushg (&gwdata->gdata, 1);
	    }
	}

/* Add or subtract n until the g is between 0 and n-1 */

	while (g->sign < 0) addg (n, g);
	while (gcompg (g, n) >= 0) subg (n, g);

/* If input was negative, return k*b^n+c - g */

	if (neg && g->sign) {
		g->sign = -g->sign;
		addg (n, g);
	}

/* Free memory */

	pushg (&gwdata->gdata, 1);
}

/******************************************************************/
/* Wrapper routines for the multiplication assembly code routines */
/******************************************************************/

/* Internal wrapper routine to call fftmul assembly code. */
/* Caller must set NORMRTN prior to calling this routine.*/

void raw_gwfftmul (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	s,
	gwnum	d)
{
	struct gwasm_data *asm_data;
	uint32_t norm_count1, norm_count2;
	double	sumdiff;

	ASSERTG (((uint32_t *) s)[-1] >= 1);
	ASSERTG (((uint32_t *) d)[-1] >= 1);

/* Get the unnormalized add count for later use */

	norm_count1 = ((uint32_t *) s)[-1];
	norm_count2 = ((uint32_t *) d)[-1];

/* Call the assembly code */

	asm_data = (struct gwasm_data *) gwdata->asm_data;
	asm_data->SRCARG = s;
	asm_data->DESTARG = d;
	gw_mul (gwdata, asm_data);
	if (! is_valid_double (gwsumout (gwdata, d))) gwdata->GWERROR |= 1;
	gwdata->fft_count += 2.0;

/* Adjust if necessary the SUM(INPUTS) vs. SUM(OUTPUTS).  If norm_count */
/* is more than one, then the sums will be larger than normal.  This */
/* could trigger a spurious MAXDIFF warning.  Shrink the two SUMS to */
/* compensate. */

	if (norm_count1 != 1 || norm_count2 != 1) {
		double	adjustment;
		adjustment = 1.0 / ((double)norm_count1 * (double)norm_count2);
		gwsuminp (gwdata, d) *= adjustment;
		gwsumout (gwdata, d) *= adjustment;
	}

/* Test SUM(INPUTS) vs. SUM(OUTPUTS) */

	sumdiff = gwsuminp (gwdata, d) - gwsumout (gwdata, d);
	if (fabs (sumdiff) > gwdata->MAXDIFF) gwdata->GWERROR |= 2; 

/* Reset the unnormalized add count */

	((uint32_t *) d)[-1] = 1;
}

/* Common code to emulate the modulo with two multiplies in the */
/* general purpose case */

void emulate_mod (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	s)		/* Source and destination */
{
	struct gwasm_data *asm_data;
	gwnum	tmp;
	double	saved_addin_value;

	ASSERTG (* addr (gwdata, s, gwdata->FFTLEN-1) > -2.0 && * addr (gwdata, s, gwdata->FFTLEN-1) <= 0.0);

/* Save and clear the addin value */

	asm_data = (struct gwasm_data *) gwdata->asm_data;
	saved_addin_value = asm_data->ADDIN_VALUE;
	asm_data->ADDIN_VALUE = 0.0;

/* Copy the number and zero out the low words. */

	tmp = gwalloc (gwdata);
	gwcopyzero (gwdata, s, tmp, gwdata->GW_ZEROWORDSLOW);

/* Multiply by the reciprocal that has been carefully shifted so that the */
/* integer part of the result wraps to the lower FFT words.  Adjust the */
/* normalization routine so that the FFT code zeroes the high FFT words */
/* and we are left with just the quotient! */

	asm_data->NORMRTN = gwdata->GWPROCPTRS[copyzero_routines + (gwdata->NORMNUM & 1)];
	raw_gwfftmul (gwdata, gwdata->GW_RECIP_FFT, tmp);
	ASSERTG (* addr (gwdata, tmp, gwdata->FFTLEN/2-1) > -2.0 && * addr (gwdata, tmp, gwdata->FFTLEN/2-1) <= 0.0);

/* Muliply quotient and modulus.  Select normalization routine that does */
/* not zero the high FFT words. */

	asm_data->NORMRTN = gwdata->GWPROCPTRS[norm_routines + (gwdata->NORMNUM & 1)];
	raw_gwfftmul (gwdata, gwdata->GW_MODULUS_FFT, tmp);

/* Subtract from the original number to get the remainder */

	gwsub (gwdata, tmp, s);
	ASSERTG (* addr (gwdata, s, gwdata->FFTLEN-1) == 0.0);
	gwfree (gwdata, tmp);

/* Restore the addin value */

	asm_data->ADDIN_VALUE = saved_addin_value;
}

/* User-visible routines */

void gwfft (			/* Forward FFT */
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	s,		/* Source number */
	gwnum	d)		/* Destination (can overlap source) */
{
	struct gwasm_data *asm_data;

	ASSERTG (((uint32_t *) s)[-1] >= 1);

/* Copy the unnormalized add count */

	((uint32_t *) d)[-1] = ((uint32_t *) s)[-1];

/* If this is a zero-padded FFT and the source has been partially FFTed */
/* (by a prior POSTFFT setting) and the destination is different than the */
/* source, then we must also copy the 7 words around the halfway point from */
/* the source to the destination. */

	asm_data = (struct gwasm_data *) gwdata->asm_data;
	if (asm_data->ZERO_PADDED_FFT &&
	    ((uint32_t *) s)[-7] == 1 && s != d) {
		d[-5] = s[-5];
		d[-6] = s[-6];
		d[-7] = s[-7];
		d[-8] = s[-8];
		d[-9] = s[-9];
		d[-10] = s[-10];
		d[-11] = s[-11];
	}

/* Call the assembly code */

	asm_data->SRCARG = s;
	asm_data->DESTARG = d;
	gw_fft (gwdata, asm_data);
	gwdata->fft_count += 1.0;
}

void gwsquare (			/* Square a number */
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	s)		/* Source and destination */
{
	struct gwasm_data *asm_data;
	uint32_t norm_count;
	double	sumdiff;

	ASSERTG (((uint32_t *) s)[-1] >= 1);

/* If we are converting gwsquare calls into gwsquare_carefully calls */
/* do so now.  Turn off option to do a partial forward FFT on the result. */
/* NOTE: We must clear count since gwsquare_carefully calls back to this */
/* gwsquare routine. */

	if (gwdata->square_carefully_count) {
		int	n = gwdata->square_carefully_count;
		if (n > 1) gwstartnextfft (gwdata, 0);
		gwdata->square_carefully_count = 0;
		gwsquare_carefully (gwdata, s);
		gwdata->square_carefully_count = n - 1;
		return;
	}

/* Get the unnormalized add count for later use */

	norm_count = ((uint32_t *) s)[-1];

/* Call the assembly code */

	asm_data = (struct gwasm_data *) gwdata->asm_data;
	asm_data->NORMRTN = gwdata->GWPROCPTRS[norm_routines + gwdata->NORMNUM];
	asm_data->DESTARG = s;
	gw_square (gwdata, asm_data);
	if (! is_valid_double (gwsumout (gwdata, s))) gwdata->GWERROR |= 1;
	gwdata->fft_count += 2.0;

/* Adjust if necessary the SUM(INPUTS) vs. SUM(OUTPUTS).  If norm_count */
/* is more than one, then the sums will be larger than normal.  This */
/* could trigger a spurious MAXDIFF warning.  Shrink the two SUMS to */
/* compensate. */

	if (norm_count != 1) {
		double	adjustment;
		adjustment = 1.0 / ((double) norm_count * (double) norm_count);
		gwsuminp (gwdata, s) *= adjustment;
		gwsumout (gwdata, s) *= adjustment;
	}

/* Test SUM(INPUTS) vs. SUM(OUTPUTS) */

	sumdiff = gwsuminp (gwdata, s) - gwsumout (gwdata, s);
	if (fabs (sumdiff) > gwdata->MAXDIFF) gwdata->GWERROR |= 2; 

/* Reset the unnormalized add count */

	((uint32_t *) s)[-1] = 1;

/* Emulate mod with 2 multiplies case */

	if (gwdata->GENERAL_MOD) emulate_mod (gwdata, s);
}

void gwfftmul (			/* Multiply already FFTed source with dest */
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	s,		/* Already FFTed source number */
	gwnum	d)		/* Non-FFTed source. Also destination */
{
	struct gwasm_data *asm_data;

	ASSERTG (((uint32_t *) s)[-1] >= 1);
	ASSERTG (((uint32_t *) d)[-1] >= 1);

/* Call the assembly code */

	asm_data = (struct gwasm_data *) gwdata->asm_data;
	asm_data->NORMRTN = gwdata->GWPROCPTRS[norm_routines + gwdata->NORMNUM];
	raw_gwfftmul (gwdata, s, d);

/* Emulate mod with 2 multiplies case */

	if (gwdata->GENERAL_MOD) emulate_mod (gwdata, d);
}

void gwfftfftmul (		/* Multiply two already FFTed sources */
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	s,		/* Already FFTed source number */
	gwnum	s2,		/* Already FFTed source number */
	gwnum	d)		/* Destination (can overlap sources) */
{
	struct gwasm_data *asm_data;
	uint32_t norm_count1, norm_count2;
	double	sumdiff;

	ASSERTG (((uint32_t *) s)[-1] >= 1);
	ASSERTG (((uint32_t *) s2)[-1] >= 1);

/* Get the unnormalized add count for later use */

	norm_count1 = ((uint32_t *) s)[-1];
	norm_count2 = ((uint32_t *) s2)[-1];

/* Call the assembly code */

	asm_data = (struct gwasm_data *) gwdata->asm_data;
	asm_data->NORMRTN = gwdata->GWPROCPTRS[norm_routines + gwdata->NORMNUM];
	asm_data->SRCARG = s;
	asm_data->SRC2ARG = s2;
	asm_data->DESTARG = d;
	gw_mulf (gwdata, asm_data);
	if (! is_valid_double (gwsumout (gwdata, d))) gwdata->GWERROR |= 1;
	gwdata->fft_count += 1.0;

/* Adjust if necessary the SUM(INPUTS) vs. SUM(OUTPUTS).  If norm_count */
/* is more than one, then the sums will be larger than normal.  This */
/* could trigger a spurious MAXDIFF warning.  Shrink the two SUMS to */
/* compensate. */

	if (norm_count1 != 1 || norm_count2 != 1) {
		double	adjustment;
		adjustment = 1.0 / ((double)norm_count1 * (double)norm_count2);
		gwsuminp (gwdata, d) *= adjustment;
		gwsumout (gwdata, d) *= adjustment;
	}

/* Test SUM(INPUTS) vs. SUM(OUTPUTS) */

	sumdiff = gwsuminp (gwdata, d) - gwsumout (gwdata, d);
	if (fabs (sumdiff) > gwdata->MAXDIFF) gwdata->GWERROR |= 2; 

/* Reset the unnormalized add count */

	((uint32_t *) d)[-1] = 1;

/* Emulate mod with 2 multiplies case */

	if (gwdata->GENERAL_MOD) emulate_mod (gwdata, d);
}

void gwmul (			/* Multiply source with dest */
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	s,		/* Source number (changed to FFTed source!) */
	gwnum	d)		/* Source and destination */
{
	gwfft (gwdata, s, s);
	gwfftmul (gwdata, s, d);
}

void gwsafemul (		/* Multiply source with dest */
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	s,		/* Source number (not changed) */
	gwnum	d)		/* Source and destination */
{
	gwnum	qqq;

	qqq = gwalloc (gwdata);
	gwfft (gwdata, s, qqq);
	gwfftmul (gwdata, qqq, d);
	gwfree (gwdata, qqq);
}

/* Generate random FFT data.  We used to use the C runtime library. */
/* However, when a caller discovered a bug in gwsquare_carefully it */
/* very difficult to track down because the bug was no reproducible. */
/* We could make bugs reproducible by calling srand with a fixed value, */
/* but it is bad form for a library to do this.  Thus, we found a */
/* public domain random number generator to use. */

void gw_random_number (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	x)
{
	struct mt_state rand_info;
	giant	g;
	unsigned long i, len;

/* Init the random generator to a reproducible state */

	init_genrand (&rand_info, 5489);

/* Generate the random number */

	len = (((unsigned long) gwdata->bit_length) >> 5) + 1;
	g = popg (&gwdata->gdata, len);
	for (i = 0; i < len; i++) {
		g->n[i] = genrand_int32 (&rand_info);
	}
	g->sign = len;
	specialmodg (gwdata, g);
	gianttogw (gwdata, g, x);
	pushg (&gwdata->gdata, 1);
}

/* The FFT selection code assumes FFT data will essentially be random data */
/* yielding pretty well understood maximum round off errors.  When working */
/* with some numbers, especially at the start of a PRP exponentiation, the */
/* FFT data is decidedly not random, leading to much larger than expected */
/* roundoff errors.  In my own PRP code, I call gwsquare_carefully for the */
/* first 30 iterations.  To make this easier (and code more readable) you */
/* can call this routine and the next n gwsquare calls will be replaced by */
/* gwsquare_carefully calls.  If you pass an n of -1, the gwnum code will */
/* use a default value for n that should be suitable for getting a PRP */
/* exponentiation into a "random data state".  This routine can be called */
/* before gwsetup is called. */

void gwset_square_carefully_count (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	int	n)		/* Number of gwsquare calls to do carefully. */
				/* If n is -1, a default value is used */
{

/* Set the default count to enough iterations for a PRP exponentiation to */
/* generate a number larger than the modulus.  Then do an extra dozen */
/* iterations to hopefully scramble the data into a nice random looking */
/* pattern. */

	if (n == -1 && gwdata->FFTLEN) n = (int) ceil (log2 (gwdata->bit_length)) + 12;

/* Now remember the count */

	gwdata->square_carefully_count = n;
}

/* Square a number using a slower method that will have reduced */
/* round-off error on non-random input data.  Caller must make sure the */
/* input number has not been partially or fully FFTed. */

void gwsquare_carefully (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	s)		/* Source and destination */
{
	struct gwasm_data *asm_data;
	gwnum	tmp1, tmp2;
	double	saved_addin_value;
	unsigned long saved_extra_bits;

/* Generate a random number, if we have't already done so */

	if (gwdata->GW_RANDOM == NULL) {
		gwdata->GW_RANDOM = gwalloc (gwdata);
		gw_random_number (gwdata, gwdata->GW_RANDOM);
	}

/* Save and clear the addin value */

	asm_data = (struct gwasm_data *) gwdata->asm_data;
	saved_addin_value = asm_data->ADDIN_VALUE;
	asm_data->ADDIN_VALUE = 0.0;

/* Make sure we do not do addquick when computing s+random. */
/* If we do not do this, then the non-randomness of s can swamp */
/* the randomness of tmp1.  An example is the first PRP iterations */
/* of 2*3^599983-1 -- s is all positive values and gwadd3 thinks */
/* there are enough extra bits to do an add quick.  This generates */
/* a tmp1 with nearly all positive values -- very bad. */

	saved_extra_bits = gwdata->EXTRA_BITS;
	gwdata->EXTRA_BITS = 0;

/* Now do the squaring using three multiplies and adds */
/* Note that during the calculation of s*random we must relax the */
/* SUMINP != SUMOUT limit.  This is because s may be non-random. */
/* For example, in the first iterations of 2*3^500327-1, s is mostly */
/* large positive values.  This means we lose some of the lower bits */
/* of precision when we calculate SUMINP.  To combate this problem we */
/* increase MAXDIFF during the s*random calculation. */

	tmp1 = gwalloc (gwdata);
	tmp2 = gwalloc (gwdata);
	gwstartnextfft (gwdata, 0);		/* Disable POSTFFT */
	gwadd3 (gwdata, s, gwdata->GW_RANDOM, tmp1); /* Compute s+random */
	gwfft (gwdata, gwdata->GW_RANDOM, tmp2);
	gwdata->MAXDIFF *= 1024.0;
	gwfftmul (gwdata, tmp2, s);		/* Compute s*random */
	gwdata->MAXDIFF /= 1024.0;
	gwfftfftmul (gwdata, tmp2, tmp2, tmp2);	/* Compute random^2 */
	asm_data->ADDIN_VALUE = saved_addin_value;/* Restore the addin value */
	gwsquare (gwdata, tmp1);		/* Compute (s+random)^2 */
	gwsubquick (gwdata, tmp2, tmp1);	/* Calc s^2 from 3 results */
	gwaddquick (gwdata, s, s);
	gwsub3 (gwdata, tmp1, s, s);

/* Restore state, free memory and return */

	gwdata->EXTRA_BITS = saved_extra_bits;
	gwfree (gwdata, tmp1);
	gwfree (gwdata, tmp2);
}

/* Multiply numbers using a slower method that will have reduced */
/* round-off error on non-random input data.  Caller must make sure the */
/* input numbers have not been partially or fully FFTed. */

void gwmul_carefully (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	s,		/* Source */
	gwnum	t)		/* Source and destination */
{
	struct gwasm_data *asm_data;
	gwnum	tmp1, tmp2, tmp3, tmp4;
	double	saved_addin_value;
	unsigned long saved_extra_bits;

/* Generate a random number, if we have't already done so */

	if (gwdata->GW_RANDOM == NULL) {
		gwdata->GW_RANDOM = gwalloc (gwdata);
		gw_random_number (gwdata, gwdata->GW_RANDOM);
	}

/* Save and clear the addin value */

	asm_data = (struct gwasm_data *) gwdata->asm_data;
	saved_addin_value = asm_data->ADDIN_VALUE;
	asm_data->ADDIN_VALUE = 0.0;

/* Make sure we do not do addquick when computing s+random. */

	saved_extra_bits = gwdata->EXTRA_BITS;
	gwdata->EXTRA_BITS = 0;

/* Now do the multiply using four multiplies and adds */

	tmp1 = gwalloc (gwdata);
	tmp2 = gwalloc (gwdata);
	tmp3 = gwalloc (gwdata);
	tmp4 = gwalloc (gwdata);
	gwcopy (gwdata, s, tmp4);
	gwstartnextfft (gwdata, 0);		/* Disable POSTFFT */
	gwadd3 (gwdata, s, gwdata->GW_RANDOM, tmp1); /* Compute s+random */
	gwadd3 (gwdata, t, gwdata->GW_RANDOM, tmp3); /* Compute t+random */
	gwfft (gwdata, gwdata->GW_RANDOM, tmp2);
	gwdata->MAXDIFF *= 1024.0;
	gwfftmul (gwdata, tmp2, tmp4);		/* Compute s*random */
	gwfftmul (gwdata, tmp2, t);		/* Compute t*random */
	gwdata->MAXDIFF /= 1024.0;
	gwfftfftmul (gwdata, tmp2, tmp2, tmp2);	/* Compute random^2 */
	asm_data->ADDIN_VALUE = saved_addin_value; /* Restore addin value */
	gwmul (gwdata, tmp1, tmp3);	/* Compute (s+random)*(t+random) */
	gwsub (gwdata, tmp2, tmp3);		/* Subtract random^2 */
	gwsub (gwdata, t, tmp3);
	gwsub3 (gwdata, tmp3, tmp4, t);

/* Restore state, free memory and return */

	gwdata->EXTRA_BITS = saved_extra_bits;
	gwfree (gwdata, tmp1);
	gwfree (gwdata, tmp2);
	gwfree (gwdata, tmp3);
	gwfree (gwdata, tmp4);
}


/*********************************************************/
/* Wrapper routines for the add and sub assembly code    */
/*********************************************************/

void gwadd3quick (		/* Add two numbers without normalizing */
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	s1,		/* Source #1 */
	gwnum	s2,		/* Source #2 */
	gwnum	d)		/* Destination */
{
	struct gwasm_data *asm_data;

	ASSERTG (((uint32_t *) s1)[-1] >= 1);
	ASSERTG (((uint32_t *) s2)[-1] >= 1);

/* We could support partially FFTed inputs if we updated the 7 zero pad */
/* values.  Until then, assert inputs are not partially FFTed. */

	ASSERTG (((uint32_t *) s1)[-7] == 0);
	ASSERTG (((uint32_t *) s2)[-7] == 0);

/* Update the count of unnormalized adds and subtracts */

	((uint32_t *) d)[-1] = ((uint32_t *) s1)[-1] + ((uint32_t *) s2)[-1];

/* Copy the has-been-partially-FFTed flag */

	((uint32_t *) d)[-7] = ((uint32_t *) s1)[-7];

/* Now do the add */

	asm_data = (struct gwasm_data *) gwdata->asm_data;
	asm_data->SRCARG = s1;
	asm_data->SRC2ARG = s2;
	asm_data->DESTARG = d;
	gw_addq (gwdata, asm_data);
}

void gwsub3quick (		/* Compute s1 - s2 without normalizing */
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	s1,		/* Source #1 */
	gwnum	s2,		/* Source #2 */
	gwnum	d)		/* Destination */
{
	struct gwasm_data *asm_data;

	ASSERTG (((uint32_t *) s1)[-1] >= 1);
	ASSERTG (((uint32_t *) s2)[-1] >= 1);

/* We could support partially FFTed inputs if we updated the 7 zero pad */
/* values.  Until then, assert inputs are not partially FFTed. */

	ASSERTG (((uint32_t *) s1)[-7] == 0);
	ASSERTG (((uint32_t *) s2)[-7] == 0);

/* Update the count of unnormalized adds and subtracts */

	((uint32_t *) d)[-1] = ((uint32_t *) s1)[-1] + ((uint32_t *) s2)[-1];

/* Copy the has-been-partially-FFTed flag */

	((uint32_t *) d)[-7] = ((uint32_t *) s1)[-7];

/* Now do the subtract */

	asm_data = (struct gwasm_data *) gwdata->asm_data;
	asm_data->SRCARG = s2;
	asm_data->SRC2ARG = s1;
	asm_data->DESTARG = d;
	gw_subq (gwdata, asm_data);
}

void gwaddsub4quick (		/* Add & sub two numbers without normalizing */
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	s1,		/* Source #1 */
	gwnum	s2,		/* Source #2 */
	gwnum	d1,		/* Destination #1 */
	gwnum	d2)		/* Destination #2 */
{
	struct gwasm_data *asm_data;

	ASSERTG (((uint32_t *) s1)[-1] >= 1);
	ASSERTG (((uint32_t *) s2)[-1] >= 1);

/* We could support partially FFTed inputs if we updated the 7 zero pad */
/* values.  Until then, assert inputs are not partially FFTed. */

	ASSERTG (((uint32_t *) s1)[-7] == 0);
	ASSERTG (((uint32_t *) s2)[-7] == 0);

/* Update the counts of unnormalized adds and subtracts */

	((uint32_t *) d1)[-1] =
	((uint32_t *) d2)[-1] = ((uint32_t *) s1)[-1] + ((uint32_t *) s2)[-1];

/* Copy the has-been-partially-FFTed flag */

	((uint32_t *) d1)[-7] =
	((uint32_t *) d2)[-7] = ((uint32_t *) s1)[-7];

/* Now do the add & subtract */

	asm_data = (struct gwasm_data *) gwdata->asm_data;
	asm_data->SRCARG = s1;
	asm_data->SRC2ARG = s2;
	asm_data->DESTARG = d1;
	asm_data->DEST2ARG = d2;
	gw_addsubq (gwdata, asm_data);
}


void gwadd3 (			/* Add two numbers normalizing if needed */
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	s1,		/* Source #1 */
	gwnum	s2,		/* Source #2 */
	gwnum	d)		/* Destination */
{
	struct gwasm_data *asm_data;
	uint32_t normcnt1, normcnt2;

	ASSERTG (((uint32_t *) s1)[-1] >= 1);
	ASSERTG (((uint32_t *) s2)[-1] >= 1);
	ASSERTG (((uint32_t *) s1)[-7] == 0);
	ASSERTG (((uint32_t *) s2)[-7] == 0);

/* Get counts of unnormalized adds and subtracts */

	normcnt1 = ((uint32_t *) s1)[-1];
	normcnt2 = ((uint32_t *) s2)[-1];

/* Set the has-been-partially-FFTed flag */

	((uint32_t *) d)[-7] = 0;

/* Now do the add */

	asm_data = (struct gwasm_data *) gwdata->asm_data;
	asm_data->SRCARG = s1;
	asm_data->SRC2ARG = s2;
	asm_data->DESTARG = d;
	if (normcnt1 + normcnt2 <= gwdata->EXTRA_BITS) {
		gw_addq (gwdata, asm_data);
		((uint32_t *) d)[-1] = normcnt1 + normcnt2;
	} else {
		gw_add (gwdata, asm_data);
		((uint32_t *) d)[-1] = 1;
	}
}

void gwsub3 (			/* Compute s1 - s2 normalizing if needed */
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	s1,		/* Source #1 */
	gwnum	s2,		/* Source #2 */
	gwnum	d)		/* Destination */
{
	struct gwasm_data *asm_data;
	uint32_t normcnt1, normcnt2;

	ASSERTG (((uint32_t *) s1)[-1] >= 1);
	ASSERTG (((uint32_t *) s2)[-1] >= 1);
	ASSERTG (((uint32_t *) s1)[-7] == 0);
	ASSERTG (((uint32_t *) s2)[-7] == 0);

/* Get counts of unnormalized adds and subtracts */

	normcnt1 = ((uint32_t *) s1)[-1];
	normcnt2 = ((uint32_t *) s2)[-1];

/* Set the has-been-partially-FFTed flag */

	((uint32_t *) d)[-7] = 0;

/* Now do the subtract */

	asm_data = (struct gwasm_data *) gwdata->asm_data;
	asm_data->SRCARG = s2;
	asm_data->SRC2ARG = s1;
	asm_data->DESTARG = d;
	if (normcnt1 + normcnt2 <= gwdata->EXTRA_BITS) {
		gw_subq (gwdata, asm_data);
		((uint32_t *) d)[-1] = normcnt1 + normcnt2;
	} else {
		gw_sub (gwdata, asm_data);
		((uint32_t *) d)[-1] = 1;
	}
}

void gwaddsub4 (		/* Add & sub two nums normalizing if needed */
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	s1,		/* Source #1 */
	gwnum	s2,		/* Source #2 */
	gwnum	d1,		/* Destination #1 */
	gwnum	d2)		/* Destination #2 */
{
	struct gwasm_data *asm_data;
	uint32_t normcnt1, normcnt2;

	ASSERTG (((uint32_t *) s1)[-1] >= 1);
	ASSERTG (((uint32_t *) s2)[-1] >= 1);
	ASSERTG (((uint32_t *) s1)[-7] == 0);
	ASSERTG (((uint32_t *) s2)[-7] == 0);

/* Get counts of unnormalized adds and subtracts */

	normcnt1 = ((uint32_t *) s1)[-1];
	normcnt2 = ((uint32_t *) s2)[-1];

/* Set the has-been-partially-FFTed flag */

	((uint32_t *) d1)[-7] = ((uint32_t *) d2)[-7] = 0;

/* Now do the add & subtract */

	asm_data = (struct gwasm_data *) gwdata->asm_data;
	asm_data->SRCARG = s1;
	asm_data->SRC2ARG = s2;
	asm_data->DESTARG = d1;
	asm_data->DEST2ARG = d2;
	if (normcnt1 + normcnt2 <= gwdata->EXTRA_BITS) {
		gw_addsubq (gwdata, asm_data);
		((uint32_t *) d1)[-1] =
		((uint32_t *) d2)[-1] = normcnt1 + normcnt2;
	} else {
		gw_addsub (gwdata, asm_data);
		((uint32_t *) d1)[-1] =
		((uint32_t *) d2)[-1] = 1;
	}
}


void gwfftadd3 (		/* Add two FFTed numbers */
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	s1,		/* Source #1 */
	gwnum	s2,		/* Source #2 */
	gwnum	d)		/* Destination */
{
	struct gwasm_data *asm_data;

	ASSERTG (((uint32_t *) s1)[-1] >= 1);
	ASSERTG (((uint32_t *) s2)[-1] >= 1);
	ASSERTG (((uint32_t *) s1)[-7] == ((uint32_t *) s2)[-7]);

/* Update the count of unnormalized adds and subtracts */

	((uint32_t *) d)[-1] = ((uint32_t *) s1)[-1] + ((uint32_t *) s2)[-1];

/* Copy the has-been-partially-FFTed flag */

	((uint32_t *) d)[-7] = ((uint32_t *) s1)[-7];

/* If this is a zero-padded FFT, then also add the 7 copied doubles in */
/* the gwnum header */

	if (gwdata->ZERO_PADDED_FFT) {
		d[-5] = s1[-5] + s2[-5];
		d[-6] = s1[-6] + s2[-6];
		d[-7] = s1[-7] + s2[-7];
		d[-8] = s1[-8] + s2[-8];
		d[-9] = s1[-9] + s2[-9];
		d[-10] = s1[-10] + s2[-10];
		d[-11] = s1[-11] + s2[-11];
	}

/* Now do the add */

	asm_data = (struct gwasm_data *) gwdata->asm_data;
	asm_data->SRCARG = s1;
	asm_data->SRC2ARG = s2;
	asm_data->DESTARG = d;
	gw_addf (gwdata, asm_data);
}

void gwfftsub3 (		/* Compute FFTed s1 - FFTed s2 */
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	s1,		/* Source #1 */
	gwnum	s2,		/* Source #2 */
	gwnum	d)		/* Destination */
{
	struct gwasm_data *asm_data;

	ASSERTG (((uint32_t *) s1)[-1] >= 1);
	ASSERTG (((uint32_t *) s2)[-1] >= 1);
	ASSERTG (((uint32_t *) s1)[-7] == ((uint32_t *) s2)[-7]);

/* Update the count of unnormalized adds and subtracts */

	((uint32_t *) d)[-1] = ((uint32_t *) s1)[-1] + ((uint32_t *) s2)[-1];

/* Copy the has-been-partially-FFTed flag */

	((uint32_t *) d)[-7] = ((uint32_t *) s1)[-7];

/* If this is a zero-padded FFT, then also subtract the 7 copied doubles in */
/* the gwnum header */

	if (gwdata->ZERO_PADDED_FFT) {
		d[-5] = s1[-5] - s2[-5];
		d[-6] = s1[-6] - s2[-6];
		d[-7] = s1[-7] - s2[-7];
		d[-8] = s1[-8] - s2[-8];
		d[-9] = s1[-9] - s2[-9];
		d[-10] = s1[-10] - s2[-10];
		d[-11] = s1[-11] - s2[-11];
	}

/* Now do the subtract */

	asm_data = (struct gwasm_data *) gwdata->asm_data;
	asm_data->SRCARG = s2;
	asm_data->SRC2ARG = s1;
	asm_data->DESTARG = d;
	gw_subf (gwdata, asm_data);
}

void gwfftaddsub4 (		/* Add & sub two FFTed numbers */
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	s1,		/* Source #1 */
	gwnum	s2,		/* Source #2 */
	gwnum	d1,		/* Destination #1 */
	gwnum	d2)		/* Destination #2 */
{
	struct gwasm_data *asm_data;

	ASSERTG (((uint32_t *) s1)[-1] >= 1);
	ASSERTG (((uint32_t *) s2)[-1] >= 1);
	ASSERTG (((uint32_t *) s1)[-7] == ((uint32_t *) s2)[-7]);

/* Update the counts of unnormalized adds and subtracts */

	((uint32_t *) d1)[-1] =
	((uint32_t *) d2)[-1] = ((uint32_t *) s1)[-1] + ((uint32_t *) s2)[-1];

/* Copy the has-been-partially-FFTed flag */

	((uint32_t *) d1)[-7] =
	((uint32_t *) d2)[-7] = ((uint32_t *) s1)[-7];

/* If this is a zero-padded FFT, then also add & sub the 7 copied doubles in */
/* the gwnum header.  Copy data to temporaries first in case s1, s2 pointers */
/* are equal to the d1, d2 pointers! */

	if (gwdata->ZERO_PADDED_FFT) {
		double	v1, v2;
		v1 = s1[-5]; v2 = s2[-5]; d1[-5] = v1 + v2; d2[-5] = v1 - v2;
		v1 = s1[-6]; v2 = s2[-6]; d1[-6] = v1 + v2; d2[-6] = v1 - v2;
		v1 = s1[-7]; v2 = s2[-7]; d1[-7] = v1 + v2; d2[-7] = v1 - v2;
		v1 = s1[-8]; v2 = s2[-8]; d1[-8] = v1 + v2; d2[-8] = v1 - v2;
		v1 = s1[-9]; v2 = s2[-9]; d1[-9] = v1 + v2; d2[-9] = v1 - v2;
		v1 = s1[-10]; v2 = s2[-10]; d1[-10] = v1+v2; d2[-10] = v1-v2;
		v1 = s1[-11]; v2 = s2[-11]; d1[-11] = v1+v2; d2[-11] = v1-v2;
	}

/* Now do the add & subtract */

	asm_data = (struct gwasm_data *) gwdata->asm_data;
	asm_data->SRCARG = s1;
	asm_data->SRC2ARG = s2;
	asm_data->DESTARG = d1;
	asm_data->DEST2ARG = d2;
	gw_addsubf (gwdata, asm_data);
}

/* Routine to add a small number to a gwnum.  Some day, */
/* I might optimize this routine for the cases where just one or two */
/* doubles need to be modified in the gwnum */

void gwsmalladd (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	double	addin,		/* Small value to add to g */
	gwnum	g)		/* Gwnum to add a value into */
{

/* Assert unnormalized add count valid, input not partially FFTed. */

	ASSERTG (((uint32_t *) g)[-1] >= 1);
	ASSERTG (((uint32_t *) g)[-7] == 0);

/* A simple brute-force implementation.  We put k > 1 numbers through */
/* here because multiplying the addin value by 1/k complicates matters. */
/* If this routine is ever used much, we can try to optimize this. */
/* We also make sure there is at least one b value per FFT work, so */
/* that any carry can successfully spread over 3 FFT words. */

	if (gwdata->GWPROCPTRS[14] == NULL || gwdata->k > 1.0 ||
	    gwdata->NUM_B_PER_SMALL_WORD < 1) {
		gwnum	tmp;
		tmp = gwalloc (gwdata);
		if (addin >= 0.0) {
			dbltogw (gwdata, addin, tmp);
			gwaddquick (gwdata, tmp, g);
		} else {
			dbltogw (gwdata, -addin, tmp);
			gwsubquick (gwdata, tmp, g);
		}
		gwfree (gwdata, tmp);
	}

/* The assembler optimized version */

	else {
		struct gwasm_data *asm_data;

		asm_data = (struct gwasm_data *) gwdata->asm_data;
		asm_data->DESTARG = g;
		asm_data->DBLARG = addin;
		gw_adds (gwdata, asm_data);
	}

}

/* This routine multiplies a gwnum by a small positive value.  This lets us apply some */
/* optimizations that cannot be performed by a full FFT multiplication. */

void gwsmallmul (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	double	mult,		/* Small value to multiply g by */
	gwnum	g)		/* Gwnum to multiply */
{
	struct gwasm_data *asm_data = (struct gwasm_data *) gwdata->asm_data;

/* Assert unnormalized add count valid, input not partially FFTed. */

	ASSERTG (((uint32_t *) g)[-1] >= 1);
	ASSERTG (((uint32_t *) g)[-7] == 0);

/* The x87 assembly version won't spread carries over multiple words. */
/* Also, the x87 assembly version won't guard against carries out of the */
/* critical top words in a zero-padded number. */
/* A simple brute-force implementation.  Note we cannot call gwmul as this */
/* will use the caller's last gwsetnormroutine value which could incorrectly */
/* multiply by a constant. */
	
	if (! (gwdata->cpu_flags & CPU_SSE2) &&
	    (mult > 1024.0 || gwdata->ZERO_PADDED_FFT)) {
		gwnum	tmp;
		tmp = gwalloc (gwdata);
		if (mult == 1.0);
		else if (mult == 2.0) {
			gwadd (gwdata, g, g);
		}
		else if (mult == 3.0) {
			gwadd3quick (gwdata, g, g, tmp);
			gwadd (gwdata, tmp, g);
		}
		else if (mult == 4.0) {
			gwaddquick (gwdata, g, g);
			gwadd (gwdata, g, g);
		}
		else if (mult == 5.0) {
			gwadd3quick (gwdata, g, g, tmp);
			gwaddquick (gwdata, tmp, tmp);
			gwadd (gwdata, tmp, g);
		}
		else if (mult == 6.0) {
			gwadd3quick (gwdata, g, g, tmp);
			gwaddquick (gwdata, tmp, g);
			gwadd (gwdata, g, g);
		}
		else if (mult == 8.0) {
			gwaddquick (gwdata, g, g);
			gwaddquick (gwdata, g, g);
			gwadd (gwdata, g, g);
		}
		else if (mult == 9.0) {
			gwadd3quick (gwdata, g, g, tmp);
			gwaddquick (gwdata, tmp, tmp);
			gwaddquick (gwdata, tmp, tmp);
			gwadd (gwdata, tmp, g);
		}
		else {
			dbltogw (gwdata, mult, tmp);
			gwfft (gwdata, tmp, tmp);
			asm_data->NORMRTN = gwdata->GWPROCPTRS[norm_routines + (gwdata->NORMNUM & 1)];
			raw_gwfftmul (gwdata, tmp, g);
		}
		gwfree (gwdata, tmp);
	}

/* The assembler optimized version */

	else {
		asm_data->DESTARG = g;
		asm_data->DBLARG = mult;
		gw_muls (gwdata, asm_data);
		((uint32_t *) g)[-1] = 1;
	}

/* If the number has gotten too large (high words should all be */
/* weighted -1 or 0) then emulate general mod with 2 multiplies */

	if (gwdata->GENERAL_MOD &&
	    (* (double *) ((char *) g + gwdata->GW_GEN_MOD_MAX_OFFSET) <= 2.0 ||
	     * (double *) ((char *) g + gwdata->GW_GEN_MOD_MAX_OFFSET) > 0.0))
		emulate_mod (gwdata, g);
}
