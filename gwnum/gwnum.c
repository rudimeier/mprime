/*----------------------------------------------------------------------
| gwnum.c
|
| This file contains the C routines and global variables that are used
| in the multi-precision arithmetic routines.  That is, all routines
| that deal with the gwnum data type.
| 
|  Copyright 2002-2010 Mersenne Research, Inc.  All rights reserved.
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
	void	*DESTARG;	/* Function destination argument */
	intptr_t DIST_TO_FFTSRCARG;/* SRCARG - DESTARG */
	intptr_t DIST_TO_MULSRCARG;/* SRC2ARG - DESTARG */
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
	double	XMM_LIMIT_BIGMAX[96]; /* Normalization constants */
	double	XMM_LIMIT_INVERSE[96];
	double	XMM_COL_MULTS[1024];
	double	XMM_TTP_FUDGE[32];
	double	XMM_TTMP_FUDGE[32];
	double	XMM_BIGVAL[2];	/* Used to round double to integer */
	double	XMM_MINUS_C[2];	/* -c stored as double */
	double	XMM_NORM012_FF[2]; /* Used in xnorm012 macros (FFTLEN/2) */
	double	XMM_LIMIT_BIGMAX_NEG[96];
	int32_t	XMM_ABSVAL[4];	/* Used to compute absolute values */
	double	XMM_P309[2];	/* Used in radix-20 macros */
	double	XMM_P809[2];
	double	XMM_P951[2];
	double	XMM_P588[2];

	double	XMM_P618[2];	/* Used in old v25 home-grown PFA-5 macros */
	double	XMM_M809[2];
	double	XMM_M262[2];
	double	XMM_M382[2];
	double	XMM_M162[2];
	double	XMM_P866[2];	/* Used in PFA-6 */
	double	XMM_BIGVAL_NEG[2];
	double	CARRY_ADJUST1_HI;
	double	CARRY_ADJUST1_LO;
	double	XMM_P924[2];	/* Used in old v25 home-grown PFA-5 macros */
	double	XMM_P383[2];
	double	XMM_M358[2];
	double	XMM_P404[2];
	double	XMM_P445[2];
	double	XMM_P180[2];
	double	XMM_P975[2];	/* Used in radix-28 macros */
	double	XMM_P623[2];
	double	XMM_P901[2];
	double	XMM_P782[2];
	double	XMM_P434[2];
	double	XMM_P223[2];

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

	uint32_t UNUSED_UINT32;
	double	DBLARG;
	double	CARRY_ADJUST7;
};

/* gwnum assembly routine pointers */

#define gw_fft(h,a)	(*(h)->GWPROCPTRS[0])(a)
#define gw_add(h,a)	(*(h)->GWPROCPTRS[1])(a)
#define gw_addq(h,a)	(*(h)->GWPROCPTRS[2])(a)
#define gw_addf(h,a)	(*(h)->GWPROCPTRS[2])(a)
#define gw_sub(h,a)	(*(h)->GWPROCPTRS[3])(a)
#define gw_subq(h,a)	(*(h)->GWPROCPTRS[4])(a)
#define gw_subf(h,a)	(*(h)->GWPROCPTRS[4])(a)
#define gw_addsub(h,a)	(*(h)->GWPROCPTRS[5])(a)
#define gw_addsubq(h,a)	(*(h)->GWPROCPTRS[6])(a)
#define gw_addsubf(h,a)	(*(h)->GWPROCPTRS[6])(a)
#define gw_copyzero(h,a) (*(h)->GWPROCPTRS[7])(a)
#define gw_adds(h,a)	(*(h)->GWPROCPTRS[8])(a)
#define gw_muls(h,a)	(*(h)->GWPROCPTRS[9])(a)
#define norm_routines	10
#define zerohigh_routines 14

/* Macro to aid in build prctabs (a prctab is an array of pointers to assembly routines) */

#define extern_decl(name)		extern void (*name)(void);
#define array_entry(name)		&name,

/* Build an array of the assembly auxillary routines (addition, subtraction, etc.) */

#define aux_decl(name)			extern_decl(name##1)	extern_decl(name##2)
#define aux_decl3(name)			extern_decl(name##1)	extern_decl(name##2)	extern_decl(name##3)

#ifndef X86_64
aux_decl(gwadd)		aux_decl(gwaddq)
aux_decl(gwsub)		aux_decl(gwsubq)
aux_decl(gwaddsub)	aux_decl(gwaddsubq)
aux_decl(gwcopyzero)	aux_decl(gwmuls)

void *x87_aux_prctab[] = {
	&gwadd1, &gwaddq1, &gwsub1, &gwsubq1, &gwaddsub1, &gwaddsubq1, &gwcopyzero1, NULL, &gwmuls1,
	&gwadd2, &gwaddq2, &gwsub2, &gwsubq2, &gwaddsub2, &gwaddsubq2, &gwcopyzero2, NULL, &gwmuls2};
#endif

aux_decl3(gwxadd)	aux_decl(gwxaddq)
aux_decl3(gwxsub)	aux_decl(gwxsubq)
aux_decl3(gwxaddsub)	aux_decl(gwxaddsubq)
aux_decl(gwxcopyzero)	aux_decl3(gwxadds)	aux_decl3(gwxmuls)

void *sse2_aux_prctab[] = {
	&gwxadd1, &gwxaddq1, &gwxsub1, &gwxsubq1, &gwxaddsub1, &gwxaddsubq1, &gwxcopyzero1, &gwxadds1, &gwxmuls1,
	&gwxadd2, &gwxaddq2, &gwxsub2, &gwxsubq2, &gwxaddsub2, &gwxaddsubq2, &gwxcopyzero2, &gwxadds2, &gwxmuls2,
	&gwxadd3, &gwxaddq2, &gwxsub3, &gwxsubq2, &gwxaddsub3, &gwxaddsubq2, &gwxcopyzero2, &gwxadds3, &gwxmuls3};

/* Now we put the normalization routines in an array so we can easily pick */
/* the normalization routines to use.  There is one table for SSE2 and one */
/* for x87.  The SSE2 table has these 820 combinations: */
/*	r or i		(rational or irrational) */
/*	1 or 2 or 2AMD or 3 or 3AMD (1 pass or 2 pass or 2 pass with partial normalization - with optional AMD prefetching) */
/*	z or zp or blank (zero upper half of result or zero-padded FFT or normal FFT */
/*	e or blank	(roundoff error checking or not) */
/*	b or blank	(b > 2 or not) */
/*	s4 or blank	(SSE4 or not) */
/*	k or blank	(k for XMM_K_HI is zero or not) */
/*	c1 or cm1 or blank (c=1, c=-1, abs(c)!=1) */
/* We also define a macro that will pick the correct entry from the array. */

#define sse2_explode(macro)			sse2_explode1(macro,xr)			sse2_explode1(macro,xi)
#define sse2_explode1(macro,name)		sse2_explode2(macro,name##1,BLEND)	sse2_explode2(macro,name##2,CORE)	sse2_explode2(macro,name##2,K8)		sse2_explode2(macro,name##3,CORE)	sse2_explode2(macro,name##3,K8)
#define sse2_explode2(macro,name,suff)		sse2_zero_explode(macro,name##z,suff)	sse2_explode3(macro,name,suff,notzp)	sse2_explode3(macro,name##zp,suff,zp)
#define sse2_zero_explode(macro,name,suff)	sse2_explode9##suff(macro,name)		sse2_explode9##suff(macro,name##e)
#define sse2_explode3(macro,name,suff,zp)	sse2_explode4(macro,name,suff,zp)	sse2_explode4(macro,name##e,suff,zp)
#define sse2_explode4(macro,name,suff,zp)	sse2_explode5(macro,name,suff,zp,notc)	sse2_explode5(macro,name##c,suff,zp,c)
#define sse2_explode5(macro,name,suff,zp,c)	sse2_explode6(macro,name,suff,zp,c)	sse2_explode6(macro,name##b,suff,zp,c)
#define sse2_explode6(macro,name,suff,zp,c)	sse2_explode7##zp(macro,name,suff,c)	sse2_explode7##zp(macro,name##s4,suff,c)
#define sse2_explode7notzp(macro,name,suff,c)	sse2_explode9##suff(macro,name)
#define sse2_explode7zp(macro,name,suff,c)	sse2_explode8##c(macro,name,suff)	sse2_explode8##c(macro,name##k,suff)
#define sse2_explode8c(macro,name,suff)		sse2_explode9##suff(macro,name)
#define sse2_explode8notc(macro,name,suff)	sse2_explode9##suff(macro,name)		sse2_explode9##suff(macro,name##c1)	sse2_explode9##suff(macro,name##cm1)
#define sse2_explode9BLEND(macro,name)		macro(name##BLEND)
#define sse2_explode9CORE(macro,name)		macro(name##CORE)
#ifndef __APPLE__
#define sse2_explode9K8(macro,name)		macro(name##K8)
#else
#define sse2_explode9K8(macro,name)		macro(name##CORE)		// Macs do not have AMD CPUs
#endif

sse2_explode(extern_decl)
void *sse2_prctab[] = { sse2_explode(array_entry) NULL };

int sse2_prctab_index (gwhandle *gwdata, int z, int e, int c)
{
	int	index = 0;

	if (! gwdata->RATIONAL_FFT) index += 410;
	if (gwdata->PASS2_SIZE) {
		index += 82;
		if (gwdata->FFT_TYPE == FFT_TYPE_RADIX_4_DWPN) index += 164;
		if (gwdata->cpu_flags & CPU_3DNOW) index += 82;
	}
	if (z) {
		if (e) index += 1;
		return (index);
	}
	index += 2;
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

/* The x87 normalization routines array has 40 combinations: */
/*	r or i		(rational or irrational) */
/*	1 or 2		(1 or 2 pass FFTs) */
/*	z or zp or blank (zero upper half of result or zero-padded FFT or normal FFT */
/*	e or blank	(roundoff error checking or not) */
/*	c or blank	(mul by small const or not)  */
/* We also define a macro that will pick the correct entry from the array. */

#ifndef X86_64
#define x87_explode(macro)		x87_explode1(macro,r)		x87_explode1(macro,i)
#define x87_explode1(macro,name)	x87_explode2(macro,name##1)	x87_explode2(macro,name##2)
#define x87_explode2(macro,name)	x87_zero_explode(macro,name##z)	x87_explode3(macro,name)	x87_explode3(macro,name##zp)
#define x87_zero_explode(macro,name)	macro(name)			macro(name##e)
#define x87_explode3(macro,name)	x87_explode4(macro,name)	x87_explode4(macro,name##e)
#define x87_explode4(macro,name)	macro(name)			macro(name##c)

x87_explode(extern_decl)
void *x87_prctab[] = { x87_explode(array_entry) NULL };

#define	x87_prctab_index(gwdata, z, e, c)  \
	    ((gwdata)->RATIONAL_FFT ? 0 : 20) + \
	    ((gwdata)->PASS2_SIZE ? 10 : 0) + \
	    (z ? (e ? 1 : 0) : 2 + \
	     ((gwdata)->ZERO_PADDED_FFT ? 4 : 0) + \
	     (e ? 2 : 0) + \
	     (c ? 1 : 0))
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

/* Routine to split a r4dwpn FFT word into column and group multiplier indexes */

unsigned long dwpn_col (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	unsigned long word)
{
	unsigned long upper_sse2_word, high, low;

	upper_sse2_word = gwdata->PASS2_SIZE >> 1;
	high = word / upper_sse2_word;
	low = word - high * upper_sse2_word;

	high = high >> 1;
	high = high % gwdata->wpn_count;
	high = high * gwdata->PASS2_SIZE;

	return (high + low);
}

/* This routine builds the sin/cos table used in pass 1 by a traditional */
/* DJB radix-4 FFT  - called by gwsetup.   If this is an all-complex FFT, */
/* then the root-of-minus-1 premultipliers are also built. */

double *r4_build_pass1_table (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	double	*table)		/* Pointer to the table to fill in */
{
	unsigned long pass1_size, pass1_increment, first_levels_size;
	unsigned long group, i, j, k, N, temp, upper_sse2_word;
	int	pow2_count;
	double	sincos[6];

/* Initialize some needed constants */

	pass1_increment = gwdata->PASS2_SIZE;
	upper_sse2_word = pass1_increment / 2;
	pass1_size = gwdata->FFTLEN / gwdata->PASS2_SIZE; /* Real values in a pass1 section */

/* Determine size of the first levels. */

	if (pass1_size % 7 == 0)
		first_levels_size = 28;
	else if (pass1_size % 5 == 0 && !gwdata->ALL_COMPLEX_FFT)
		first_levels_size = 20;
	else
		first_levels_size = 8;

/* Count the power-of-two FFT levels after the initial FFT levels.  If odd, the */
/* the last levels will be a radix-8, if even the last levels wil be a radix-4. */

	pass1_size /= first_levels_size;
	for (pow2_count = 0; (pass1_size & 1) == 0; pass1_size /= 2) pow2_count++;

/* Set pointer to table of multipliers */

	gwdata->adjusted_pass1_premults = table;

/* Loop through all the pass 1 groups in the same order the assembly code will */
/* process the groups. */

	for (group = 0; group < upper_sse2_word; group += gwdata->PASS1_CACHE_LINES) {

		pass1_size = gwdata->FFTLEN / gwdata->PASS2_SIZE; /* Real values in a pass1 section */
		pass1_size /= first_levels_size; /* Complex values we're generating sin/cos data for */
		N = gwdata->PASS2_SIZE;

/* Output the sin/cos/premultiplier values for the radix-8 block that does the */
/* last 3 levels in pass 1.  NOTE:  We do not need the "j loop" (it would loop */
/* from zero to zero) when generating the sin/cos twiddle factors for the last */
/* levels of pass 1. */

		if (pow2_count & 1) {
			N = N * 8;

/* For the r4_g8cl_sixteen_reals_eight_complex_djbfft building block, output the extra */
/* sin/cos values needed for the sixteen_reals.  The sixteen_reals is done in three parts: */
/* 4 four_reals_fft, 1 eight_reals_fft, and 1 four-complex_fft.  The four_reals_fft needs */
/* four sin/cos twiddles using 2*N, the eight_reals_fft and four_complex_fft can use the same */
/* sin/cos twiddles that will be generated later for the r4r8_g8cl_eight_complex_fft8 macro. */

			if (!gwdata->ALL_COMPLEX_FFT) {
				for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
					temp = (group + i);
					gwsincos (temp, N*2, (double *) &sincos);
					table[0] = sincos[0];	/* temp for 2*N */
					table[2] = sincos[1];
					gwsincos (temp + pass1_increment, N*2, (double *) &sincos);
					table[4] = sincos[0];	/* temp + pass1_increment for 2*N */
					table[6] = sincos[1];
					gwsincos (temp + 2 * pass1_increment, N*2, (double *) &sincos);
					table[8] = sincos[0];	/* temp + 2 * pass1_increment for 2*N */
					table[10] = sincos[1];
					gwsincos (temp + 3 * pass1_increment, N*2, (double *) &sincos);
					table[12] = sincos[0];	/* temp + 3 * pass1_increment for 2*N */
					table[14] = sincos[1];

					temp += upper_sse2_word;
					gwsincos (temp, N*2, (double *) &sincos);
					table[1] = sincos[0];	/* temp for 2*N */
					table[3] = sincos[1];
					gwsincos (temp + pass1_increment, N*2, (double *) &sincos);
					table[5] = sincos[0];	/* temp + pass1_increment for 2*N */
					table[7] = sincos[1];
					gwsincos (temp + 2 * pass1_increment, N*2, (double *) &sincos);
					table[9] = sincos[0];	/* temp + 2 * pass1_increment for 2*N */
					table[11] = sincos[1];
					gwsincos (temp + 3 * pass1_increment, N*2, (double *) &sincos);
					table[13] = sincos[0];	/* temp + 3 * pass1_increment for 2*N */
					table[15] = sincos[1];

					table += 16;
				}
			}

/* Output the sin/cos values for the complex group -- specifically the r8_sg8cl_eight_complex_djbfft macro. */

			for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {

/* We now calculate the standard sin/cos twiddle factors (temp*0,1,2,3,4,5,6,7 roots of N) */

				temp = group + i;
				gwsincos (0, N, (double *) &sincos);
				table[0] = sincos[0];	/* temp*0 */
				table[2] = sincos[1];
				gwsincos (temp, N, (double *) &sincos);
				table[4] = sincos[0];	/* temp*1 */
				table[6] = sincos[1];
				gwsincos (temp * 2, N, (double *) &sincos);
				table[8] = sincos[0];	/* temp*2 */
				table[10] = sincos[1];
				gwsincos (temp * 3, N, (double *) &sincos);
				table[12] = sincos[0];	/* temp*3 */
				table[14] = sincos[1];
				gwsincos (temp * 4, N, (double *) &sincos);
				table[16] = sincos[0];	/* temp*4 */
				table[18] = sincos[1];
				gwsincos (temp * 5, N, (double *) &sincos);
				table[20] = sincos[0];	/* temp*5 */
				table[22] = sincos[1];
				gwsincos (temp * 6, N, (double *) &sincos);
				table[24] = sincos[0];	/* temp*6 */
				table[26] = sincos[1];
				gwsincos (temp * 7, N, (double *) &sincos);
				table[28] = sincos[0];	/* temp*7 */
				table[30] = sincos[1];

				temp = group + i + upper_sse2_word;
				gwsincos (0, N, (double *) &sincos);
				table[1] = sincos[0];	/* temp*0 */
				table[3] = sincos[1];
				gwsincos (temp, N, (double *) &sincos);
				table[5] = sincos[0];	/* temp*1 */
				table[7] = sincos[1];
				gwsincos (temp * 2, N, (double *) &sincos);
				table[9] = sincos[0];	/* temp*2 */
				table[11] = sincos[1];
				gwsincos (temp * 3, N, (double *) &sincos);
				table[13] = sincos[0];	/* temp*3 */
				table[15] = sincos[1];
				gwsincos (temp * 4, N, (double *) &sincos);
				table[17] = sincos[0];	/* temp*4 */
				table[19] = sincos[1];
				gwsincos (temp * 5, N, (double *) &sincos);
				table[21] = sincos[0];	/* temp*5 */
				table[23] = sincos[1];
				gwsincos (temp * 6, N, (double *) &sincos);
				table[25] = sincos[0];	/* temp*6 */
				table[27] = sincos[1];
				gwsincos (temp * 7, N, (double *) &sincos);
				table[29] = sincos[0];	/* temp*7 */
				table[31] = sincos[1];

				table += 32;
			}
			pass1_size /= 8;
		}

/* For the r4_four_complex_djbfft building block levels, output the sin/cos values. */

		while ((pass1_size & 3) == 0) {
			N = N * 4;

/* For the r4_eight_reals_four_complex_djbfft building block levels, output the extra */
/* sin/cos values needed for the eight_reals.  The eight_reals doubles N because */
/* the real part of the FFT is one level behind the complex part of the FFT. */

			if (!gwdata->ALL_COMPLEX_FFT) {
				for (j = 0; j < N*2 / 8; j += pass1_increment) {
					for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
						temp = (group + j + i);
						gwsincos5 (temp, N*2, (double *) &sincos);
						table[0] = sincos[0];	/* temp */
						table[2] = sincos[1];
						table[4] = sincos[2];	/* 5 * temp */
						table[6] = sincos[3];

						temp += upper_sse2_word;
						gwsincos5 (temp, N*2, (double *) &sincos);
						table[1] = sincos[0];	/* temp */
						table[3] = sincos[1];
						table[5] = sincos[2];	/* 5 * temp */
						table[7] = sincos[3];

						table += 8;
					}
				}
			}

/* Output the sin/cos value for the complex sections, used by the r4_four_complex_djbfft macro */

			for (j = 0; j < N / 4; j += pass1_increment) {
				for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
					temp = (group + j + i);
					gwsincos2 (temp, N, (double *) &sincos);
					table[0] = sincos[0];	/* temp */
					table[2] = sincos[1];
					table[4] = sincos[2];	/* 2 * temp */
					table[6] = sincos[3];

					gwsincos2 (temp + upper_sse2_word, N, (double *) &sincos);
					table[1] = sincos[0];	/* temp */
					table[3] = sincos[1];
					table[5] = sincos[2];	/* 2 * temp */
					table[7] = sincos[3];

					table += 8;
				}
			}
			pass1_size /= 4;
		}

/* For the r5_five_complex_djbfft building block levels, output the sin/cos values. */

		while (pass1_size % 5 == 0) {
			N = N * 5;

/* The r5_ten_reals building blocks require extra sin/cos */
/* values.  The ten_reals doubles N because the real part of the */
/* FFT is one level behind the complex part of the FFT. */

			if (!gwdata->ALL_COMPLEX_FFT) {
				for (j = 0; j < N / 5; j += pass1_increment) {
					for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
						temp = (group + j + i);
						gwsincos3 (temp, N*2, (double *) &sincos);
						table[0] = sincos[0];	/* temp */
						table[2] = sincos[1];
						table[4] = sincos[4];	/* 3 * temp */
						table[6] = sincos[5];

						temp += upper_sse2_word;
						gwsincos3 (temp, N*2, (double *) &sincos);
						table[1] = sincos[0];	/* temp */
						table[3] = sincos[1];
						table[5] = sincos[4];	/* 3 * temp */
						table[7] = sincos[5];

						table += 8;
					}
				}
			}

/* Output the sin/cos data for the complex sections, (the r5_five_complex_djbfft building block). */

			for (j = 0; j < N / 5; j += pass1_increment) {
				for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
					temp = (group + j + i);
					gwsincos2 (temp, N, (double *) &sincos);
					table[0] = sincos[0];	/* temp */
					table[2] = sincos[1];
					table[4] = sincos[2];	/* 2 * temp */
					table[6] = sincos[3];

					gwsincos2 (temp + upper_sse2_word, N, (double *) &sincos);
					table[1] = sincos[0];	/* temp */
					table[3] = sincos[1];
					table[5] = sincos[2];	/* 2 * temp */
					table[7] = sincos[3];

					table += 8;
				}
			}
			pass1_size /= 5;
		}

/* For the r3_three_complex_djbfft building block levels, output the sin/cos values. */

		while (pass1_size % 3 == 0) {
			N = N * 3;

/* The r3_six_reals building blocks require an extra sin/cos */
/* value.  The six_reals doubles N because the real part of the */
/* FFT is one level behind the complex part of the FFT. */

			if (!gwdata->ALL_COMPLEX_FFT) {
				for (j = 0; j < N / 3; j += pass1_increment) {
					for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
						temp = (group + j + i);
						gwsincos (temp, N*2, (double *) &sincos);
						table[0] = sincos[0];	/* temp */
						table[2] = sincos[1];

						temp += upper_sse2_word;
						gwsincos (temp, N*2, (double *) &sincos);
						table[1] = sincos[0];	/* temp */
						table[3] = sincos[1];

						table += 4;
					}
				}
			}

/* Output the sin/cos data for the complex sections (used by the r3_three_complex_djbfft building block). */

			for (j = 0; j < N / 3; j += pass1_increment) {
				for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
					temp = (group + j + i);
					gwsincos (temp, N, (double *) &sincos);
					table[0] = sincos[0];	/* temp */
					table[2] = sincos[1];

					gwsincos (temp + upper_sse2_word, N, (double *) &sincos);
					table[1] = sincos[0];	/* temp */
					table[3] = sincos[1];

					table += 4;
				}
			}
			pass1_size /= 3;
		}
		ASSERTG (pass1_size == 1);

/* Real FFTs output one last set of sin/cos values for the first 20-reals, 28-reals, or 8-reals FFT. */

		if (! gwdata->ALL_COMPLEX_FFT) {
			N = gwdata->FFTLEN;
			pass1_size = gwdata->FFTLEN / gwdata->PASS2_SIZE; /* Real values in a pass1 section */
			if (pass1_size % 20 == 0) {
				for (j = 0; j < N / 20; j += pass1_increment) {
					for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
						temp = (group + j + i);
						for (k = 1; k <= 9; k++) {	/* Create 9 twiddle factors */
							gwsincos (k * temp, N, (double *) &sincos);
							table[0] = sincos[0];
							table[2] = sincos[1];
							gwsincos (k * (temp + upper_sse2_word), N, (double *) &sincos);
							table[1] = sincos[0];
							table[3] = sincos[1];
							table += 4;
						}
					}
				}
			}
			else if (pass1_size % 28 == 0) {
				for (j = 0; j < N / 28; j += pass1_increment) {
					for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
						temp = (group + j + i);
						for (k = 1; k <= 13; k++) {	/* Create 13 twiddle factors */
							gwsincos (k * temp, N, (double *) &sincos);
							table[0] = sincos[0];
							table[2] = sincos[1];
							gwsincos (k * (temp + upper_sse2_word), N, (double *) &sincos);
							table[1] = sincos[0];
							table[3] = sincos[1];
							table += 4;
						}
					}
				}
			}
			else {
				for (j = 0; j < N / 8; j += pass1_increment) {
					for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
						temp = (group + j + i);
						gwsincos25 (temp, N, (double *) &sincos);
						table[0] = sincos[0];	/* temp */
						table[2] = sincos[1];
						table[4] = sincos[2];	/* 2 * temp */
						table[6] = sincos[3];
						table[8] = sincos[4];	/* 5 * temp */
						table[10] = sincos[5];

						temp += upper_sse2_word;
						gwsincos25 (temp, N, (double *) &sincos);
						table[1] = sincos[0];	/* temp */
						table[3] = sincos[1];
						table[5] = sincos[2];	/* 2 * temp */
						table[7] = sincos[3];
						table[9] = sincos[4];	/* 5 * temp */
						table[11] = sincos[5];

						table += 12;
					}
				}
			}
		}

/* For all-complex FFTs, the first FFT level converts from real values to all complex */
/* values by multiplying by a root of -1 weight and doing a butterfly.  This is */
/* simplified because weights in the bottom half are sqrt(-1) times the */
/* matching weights in the upper half.  Thus, we butterfly upper_word * weight */
/* with bottom_word * i * weight.  That equals (upper_word + i * lower_word) * weight */
/* That is just a complex multiply with the half of the butterfly output values */
/* unneeded thanks to Hermetian symmetry. */

/* The second and third levels do a standard radix-4 four-complex-FFT */
/* building block with a post-multiply by a sin/cos root of unity. */

		else {
			N = gwdata->FFTLEN / 2;
			for (j = 0; j < N / 4; j += pass1_increment) {
				for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
					temp = (group + j + i);

/* Here we compute the standard 0,1,2,3 * temp for the radix-4 sin/cos multipliers. */
/* Then we multiply in the roots-of-minus-one premultiplier.  The root-of-minus-one */
/* premultiplier was for 2N, and a root-of-minus-one-of-2N is the same as a root */
/* unity for 4N. */

					gwsincos (temp, N * 4, (double *) &sincos);
					table[0] = sincos[0];	/* premult  */
					table[2] = sincos[1];
					gwsincos (5 * temp, N * 4, (double *) &sincos);
					table[4] = sincos[0];	/* temp + premult */
					table[6] = sincos[1];
					gwsincos (9 * temp, N * 4, (double *) &sincos);
					table[8] = sincos[0];	/* 2 * temp + premult */
					table[10] = sincos[1];
					gwsincos (13 * temp, N * 4, (double *) &sincos);
					table[12] = sincos[0];	/* 3 * temp + premult */
					table[14] = sincos[1];

					temp += upper_sse2_word;
					gwsincos (temp, N * 4, (double *) &sincos);
					table[1] = sincos[0];	/* premult  */
					table[3] = sincos[1];
					gwsincos (5 * temp, N * 4, (double *) &sincos);
					table[5] = sincos[0];	/* temp + premult */
					table[7] = sincos[1];
					gwsincos (9 * temp, N * 4, (double *) &sincos);
					table[9] = sincos[0];	/* 2 * temp + premult */
					table[11] = sincos[1];
					gwsincos (13 * temp, N * 4, (double *) &sincos);
					table[13] = sincos[0];	/* 3 * temp + premult */
					table[15] = sincos[1];

					table += 16;
				}
			}
		}

/* Calculate the size of each group's sin/cos/premult data for pass1_get_next_block */

		if (group == 0)
			gwdata->pass1_premult_block_size = (unsigned long)
				((char *) table - (char *) gwdata->adjusted_pass1_premults) /
				gwdata->PASS1_CACHE_LINES;
	}

/* Return address of the end of the table */

	return (table);
}

/* This routine builds the sin/cos table used in all-complex pass 2 */
/* blocks in a traditional radix-4 FFT  - called by gwsetup. */

double *r4_build_pass2_complex_table (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	double	*table)		/* Pointer to the table to fill in */
{
	unsigned long i, N, limit, aux_table_size;
	double	sincos[6];

/* We also build a smaller auxillary table so the final several levels aren't using */
/* cache-unfriendly large strides to access sin/cos data.  If the last levels use */
/* eight_complex macros set auxillary table size to 512, otherwise the last levels */
/* use the four_complex macros and we'll build an auxillary table size of 256. */

	if (gwdata->PASS2_SIZE <= 640 && gwdata->PASS2_SIZE % 3 != 0)
		aux_table_size = 0;
	else {
		for (i = 0, N = gwdata->PASS2_SIZE; (N & 1) == 0; N >>= 1) i++;
		if (i < 8)
			aux_table_size = 1 << i;
		else
			aux_table_size = (i & 1) ? 512 : 256;
	}
	
/* If the pass 2 size is divisible by 3, then the first level does a */
/* radix-3 step which only requires one sin/cos value.  Alas, rather than */
/* computing the needed N/3 sin/cos values we must compute 1/2 or 2/5 times N */
/* sin/cos values for use in later r4_x4cl_2sc_four_complex_djbfft or */
/* r5_x5cl_2sc_five_complex_djbfft levels. */

	N = gwdata->PASS2_SIZE;
	if (N % 3 == 0) {
		if (!aux_table_size ||  N / aux_table_size % 2 == 0) limit = N / 2;
		else if (N % 5 == 0) limit = N * 2 / 5;
		else limit = N / 3;
		for (i = 0; i < limit; i++) {
			gwsincos (i, N, (double *) &sincos);
			table[1] = table[0] = sincos[0];	/* temp */
			table[3] = table[2] = sincos[1];
			table += 4;
		}
	}

/* For the radix-4 and radix-5 building blocks, output two sin/cos values. */
/* If the levels above the aux_table_size are all radix-5, then we can */
/* output a slightly smaller sin/cos table. */

	else {
		if (!aux_table_size ||  N / aux_table_size % 2 == 0) limit = N / 4;
		else limit = N / 5;
		for (i = 0; i < limit; i++) {
			gwsincos2 (i, N, (double *) &sincos);
			table[1] = table[0] = sincos[0];	/* temp */
			table[3] = table[2] = sincos[1];
			table[5] = table[4] = sincos[2];	/* 2 * temp */
			table[7] = table[6] = sincos[3];
			table += 8;
		}
	}

/* Build the smaller auxillary table */

	if (aux_table_size) {
		N = aux_table_size;
		for (i = 0; i < N / 4; i++) {
			gwsincos2 (i, N, (double *) &sincos);
			table[1] = table[0] = sincos[0];	/* temp */
			table[3] = table[2] = sincos[1];
			table[5] = table[4] = sincos[2];	/* 2 * temp */
			table[7] = table[6] = sincos[3];
			table += 8;
		}
	}

/* Return address of the end of the table */

	return (table);
}

/* This routine builds the sin/cos table used in pass 2 of real FFTs, */

double *r4_build_pass2_real_table (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	double	*table)		/* Pointer to the table to fill in */
{
	unsigned long i, N;
	double	sincos[6];

/* All complex FFTs, don't need these tables */

	if (gwdata->ALL_COMPLEX_FFT) return (table);

/* Real FFTs with an initial radix-3 step needs w^n for 2N. */

	N = gwdata->PASS2_SIZE;
	if (N % 3 == 0) {
		for (i = 0; i < N / 3; i++) {
			gwsincos (i, N * 2, (double *) &sincos);
			table[0] = sincos[0];	/* temp for 2*N */
			table[1] = sincos[1];
			table += 2;
		}
		while (N % 3 == 0) N = N / 3;
	}

/* Real FFTs with a radix-5 step needs w^n and w^3n for 2N. */

	if (N % 5 == 0) {
		for (i = 0; i < N / 5; i++) {
			gwsincos3 (i, N * 2, (double *) &sincos);
			table[0] = sincos[0];	/* temp for 2*N */
			table[1] = sincos[1];
			table[2] = sincos[4];	/* 3 * temp for 2*N */
			table[3] = sincos[5];
			table += 4;
		}
		while (N % 5 == 0) N = N / 5;
	}

/* Real FFTs with eight_reals macros need w^n and w^5n for 2N. */

	for (i = 0; i < N / 4; i++) {
		gwsincos5 (i, N * 2, (double *) &sincos);
		table[0] = sincos[0];	/* temp for 2*N */
		table[1] = sincos[1];
		table[2] = sincos[2];	/* 5 * temp for 2*N */
		table[3] = sincos[3];
		table += 4;
	}

/* Return address of the end of the table */

	return (table);
}


/* This routine builds the fixed sin/cos premultiplier table used in pass 1 by the radix-4 */
/* delayed multiplier FFT - called by gwsetup.  It is like the traditional radix-4 FFT */
/* but with the all-complex premultipliers split into two parts to reduce memory. */

double *r4delay_build_fixed_premult_table (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	double	*table)		/* Pointer to the table to fill in */
{
	unsigned long pass1_increment, j, N;
	double	sincos[6];

/* For all-complex FFTs, build the fixed roots-of-minus-one table. */
/* Output these roots in the same order they will be used the first two */
/* levels of pass 1. */

	if (gwdata->ALL_COMPLEX_FFT) {
		N = gwdata->FFTLEN / 2;
		pass1_increment = gwdata->PASS2_SIZE;
		for (j = 0; j < N / 4; j += pass1_increment) {

/* Here we compute the roots-of-minus-one premultiplier.  The root-of-minus-one */
/* premultiplier was for 2N, and a root-of-minus-one-of-2N is the same as a root */
/* unity for 4N. */					

			gwsincos (j, N * 4, (double *) &sincos);
			table[0] = table[1] = sincos[0];	/* premult  */
			table[2] = table[3] = sincos[1];
			gwsincos (j + N / 4, N * 4, (double *) &sincos);
			table[4] = table[5] = sincos[0];	/* premult * -1 ^ 1/8 */
			table[6] = table[7] = sincos[1];
			gwsincos (j + 2 * N / 4, N * 4, (double *) &sincos);
			table[8] = table[9] = sincos[0];	/* premult * -1 ^ 2/8 */
			table[10] = table[11] = sincos[1];
			gwsincos (j + 3 * N / 4, N * 4, (double *) &sincos);
			table[12] = table[13] = sincos[0];	/* premult * -1 ^ 3/8 */
			table[14] = table[15] = sincos[1];

			table += 16;
		}
	}

/* Return address of the end of the table */

	return (table);
}

/* This routine builds the sin/cos table used in pass 1 by the radix-4/8 DJB */
/* FFT with delayed sin/cos multiplies. */

double *r4delay_build_pass1_table (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	double	*table)		/* Pointer to the table to fill in */
{
	unsigned long pass1_size, pass1_increment, delay_count;
	unsigned long group, i, j, k, N, temp, upper_sse2_word;
	int	pow2_count;
	double	sincos[6];

/* Initialize some needed constants */

	pass1_increment = gwdata->PASS2_SIZE;
	upper_sse2_word = pass1_increment / 2;
	pass1_size = gwdata->FFTLEN / gwdata->PASS2_SIZE; /* Real values in a pass1 section */

/* Determine number of delay groups.  In a standard radix-4 FFT, there is only one sin/cos */
/* group in the last pass 1 level.  We reduce our memory usage by using just one fixed sin/cos */
/* table in the first FFT levels and having multiple groups of sin/cos data in the last pass 1 level. */
/* I call these groups of sin/cos data in the last pass 1 level "delay groups". */

	if (pass1_size % 7 == 0)
		delay_count = 14;
	else if (pass1_size % 5 == 0 && !gwdata->ALL_COMPLEX_FFT)
		delay_count = 10;
 	else if ((pass1_size == 512 && !gwdata->ALL_COMPLEX_FFT) ||
		 (pass1_size == 1024 && !gwdata->ALL_COMPLEX_FFT) ||
		 (pass1_size == 2560 && gwdata->ALL_COMPLEX_FFT) ||
		 (pass1_size == 5120 && gwdata->ALL_COMPLEX_FFT) ||
		 pass1_size == 1536 || pass1_size == 2048 || pass1_size == 3072 || pass1_size == 4096)
		delay_count = 16;
	else
		delay_count = 4;

/* Count the power-of-two FFT levels after the initial FFT levels.  If odd, the */
/* the last levels will be a radix-8, if even the last levels wil be a radix-4. */

	pass1_size /= (delay_count * 2);
	for (pow2_count = 0; (pass1_size & 1) == 0; pass1_size /= 2) pow2_count++;

/* Set pointer to table of multipliers */

	gwdata->adjusted_pass1_premults = table;

/* Loop through all the pass 1 groups in the same order the assembly code will */
/* process the groups. */

	for (group = 0; group < upper_sse2_word; group += gwdata->PASS1_CACHE_LINES) {

		pass1_size = gwdata->FFTLEN / gwdata->PASS2_SIZE; /* Real values in a pass1 section */
		pass1_size /= (delay_count * 2);	/* Complex values we're generating sin/cos data for */
		N = gwdata->PASS2_SIZE;

/* Output the sin/cos/premultiplier values for the radix-8 block that does the */
/* last 3 levels in pass 1.  NOTE:  We do not need the "j loop" (it would loop */
/* from zero to zero) when generating the sin/cos twiddle factors for the last */
/* levels of pass 1. */

		if (pow2_count & 1) {
			N = N * 8;

/* For the r4_g8cl_sixteen_reals_eight_complex_djbfft building block, output the extra */
/* sin/cos values needed for the sixteen_reals.  The sixteen_reals is done in three parts: */
/* 4 four_reals_fft, 1 eight_reals_fft, and 1 four-complex_fft.  The four_reals_fft needs */
/* four sin/cos twiddles using 2*N, the eight_reals_fft and four_complex_fft can use the same */
/* sin/cos twiddles that will be generated later for the r4r8_g8cl_eight_complex_fft8 macro. */

			if (!gwdata->ALL_COMPLEX_FFT) {
				for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
					temp = (group + i);
					gwsincos (temp, N*2, (double *) &sincos);
					table[0] = sincos[0];	/* temp for 2*N */
					table[2] = sincos[1];
					gwsincos (temp + pass1_increment, N*2, (double *) &sincos);
					table[4] = sincos[0];	/* temp + pass1_increment for 2*N */
					table[6] = sincos[1];
					gwsincos (temp + 2 * pass1_increment, N*2, (double *) &sincos);
					table[8] = sincos[0];	/* temp + 2 * pass1_increment for 2*N */
					table[10] = sincos[1];
					gwsincos (temp + 3 * pass1_increment, N*2, (double *) &sincos);
					table[12] = sincos[0];	/* temp + 3 * pass1_increment for 2*N */
					table[14] = sincos[1];

					temp += upper_sse2_word;
					gwsincos (temp, N*2, (double *) &sincos);
					table[1] = sincos[0];	/* temp for 2*N */
					table[3] = sincos[1];
					gwsincos (temp + pass1_increment, N*2, (double *) &sincos);
					table[5] = sincos[0];	/* temp + pass1_increment for 2*N */
					table[7] = sincos[1];
					gwsincos (temp + 2 * pass1_increment, N*2, (double *) &sincos);
					table[9] = sincos[0];	/* temp + 2 * pass1_increment for 2*N */
					table[11] = sincos[1];
					gwsincos (temp + 3 * pass1_increment, N*2, (double *) &sincos);
					table[13] = sincos[0];	/* temp + 3 * pass1_increment for 2*N */
					table[15] = sincos[1];

					table += 16;
				}
			}

/* Output the sin/cos values for the complex delay groups -- specifically the r8_sg8cl_eight_complex_fft8 macro. */

			for (k = 0; k < delay_count; k++) {
				for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
					unsigned long bigN, ktemp, actemp;

/* If this is an all-complex FFT, the roots of minus 1 (same as roots of FFTLEN*2) are */
/* split to reduce memory requirements.  We apply the part of the all-complex premultiplier here. */

					if (gwdata->ALL_COMPLEX_FFT) {
						bigN = gwdata->FFTLEN * 2;
						actemp = group + i;
					} else {
						bigN = gwdata->FFTLEN;
						actemp = 0;
					}

/* Factor in the delayed part of the sin/cos multiplies from the first 2 levels.  In the first 2 levels */
/* we use a fixed sin/cos table based only on j, leaving the group+i part to be applied here by */
/* creating delay_count table entries.  For an all-complex FFT we multiply by 0,2,1,-1*(group+i) with N = FFTLEN/2. */
/* For an all-real FFT we multiply by 0,2,1,5*(group+i) with N = FFTLEN. */

					if (gwdata->ALL_COMPLEX_FFT && delay_count == 4) {
						if (k == 0)
							ktemp = 0;
						else if (k == 1)
							ktemp = 2 * (group + i) * 4;
						else if (k == 2)
							ktemp = 1 * (group + i) * 4;
						else
							ktemp = bigN - 1 * (group + i) * 4;
					} else if (gwdata->ALL_COMPLEX_FFT) {
						/* 0,2,1,-1 combined with 0,2,1,-1 */
						int	kmap[16] = {0,8,4,-4,2,10,6,-2,1,9,5,-3,-1,7,3,-5};
						if (kmap[k] >= 0)
							ktemp = kmap[k] * (group + i) * 4;
						else
							ktemp = bigN + kmap[k] * (group + i) * 4;
					} else if (delay_count == 4) {
						if (k == 0)
							ktemp = 0;
						else if (k == 1)
							ktemp = 2 * (group + i);
						else if (k == 2)
							ktemp = 1 * (group + i);
						else
							ktemp = 5 * (group + i);
					} else if (delay_count == 16) {
						/* 0,2,1,5 combined with one 0,2,1,5 and three 0,2,1,-1 */
						int	kmap[16] = {0,8,4,20,2,18,10,-6,1,17,9,-7,5,21,13,-3};
						if (kmap[k] >= 0)
							ktemp = kmap[k] * (group + i);
						else
							ktemp = bigN + kmap[k] * (group + i);
					} else {			/* delay_count == 10 or 14 */
						/* Multipliers for the radix-20 or radix-28 step */
						ktemp = k * (group + i);
					}

/* We now calculate the standard sin/cos twiddle factors (temp*0,1,2,3,4,5,6,7 roots of N) */
/* combined with the group all-complex premultiplier roots of minus 1 (same as roots of FFTLEN*2) */
/* combined with the delayed group multipliers. */

					temp = (group + i) * (bigN / N);
					gwsincos (actemp + ktemp, bigN, (double *) &sincos);
					table[0] = sincos[0];	/* premult,delay and temp*0 */
					table[2] = sincos[1];
					gwsincos (actemp + ktemp + temp, bigN, (double *) &sincos);
					table[4] = sincos[0];	/* premult,delay and temp*1 */
					table[6] = sincos[1];
					gwsincos (actemp + ktemp + temp * 2, bigN, (double *) &sincos);
					table[8] = sincos[0];	/* premult,delay and temp*2 */
					table[10] = sincos[1];
					gwsincos (actemp + ktemp + temp * 3, bigN, (double *) &sincos);
					table[12] = sincos[0];	/* premult,delay and temp*3 */
					table[14] = sincos[1];
					gwsincos (actemp + ktemp + temp * 4, bigN, (double *) &sincos);
					table[16] = sincos[0];	/* premult,delay and temp*4 */
					table[18] = sincos[1];
					gwsincos (actemp + ktemp + temp * 5, bigN, (double *) &sincos);
					table[20] = sincos[0];	/* premult,delay and temp*5 */
					table[22] = sincos[1];
					gwsincos (actemp + ktemp + temp * 6, bigN, (double *) &sincos);
					table[24] = sincos[0];	/* premult,delay and temp*6 */
					table[26] = sincos[1];
					gwsincos (actemp + ktemp + temp * 7, bigN, (double *) &sincos);
					table[28] = sincos[0];	/* premult,delay and temp*7 */
					table[30] = sincos[1];

					temp = (group + i + upper_sse2_word) * (bigN / N);
					if (gwdata->ALL_COMPLEX_FFT) actemp += upper_sse2_word;
					gwsincos (actemp + ktemp, bigN, (double *) &sincos);
					table[1] = sincos[0];	/* premult,delay and temp*0 */
					table[3] = sincos[1];
					gwsincos (actemp + ktemp + temp, bigN, (double *) &sincos);
					table[5] = sincos[0];	/* premult,delay and temp*1 */
					table[7] = sincos[1];
					gwsincos (actemp + ktemp + temp * 2, bigN, (double *) &sincos);
					table[9] = sincos[0];	/* premult,delay and temp*2 */
					table[11] = sincos[1];
					gwsincos (actemp + ktemp + temp * 3, bigN, (double *) &sincos);
					table[13] = sincos[0];	/* premult,delay and temp*3 */
					table[15] = sincos[1];
					gwsincos (actemp + ktemp + temp * 4, bigN, (double *) &sincos);
					table[17] = sincos[0];	/* premult,delay and temp*4 */
					table[19] = sincos[1];
					gwsincos (actemp + ktemp + temp * 5, bigN, (double *) &sincos);
					table[21] = sincos[0];	/* premult,delay and temp*5 */
					table[23] = sincos[1];
					gwsincos (actemp + ktemp + temp * 6, bigN, (double *) &sincos);
					table[25] = sincos[0];	/* premult,delay and temp*6 */
					table[27] = sincos[1];
					gwsincos (actemp + ktemp + temp * 7, bigN, (double *) &sincos);
					table[29] = sincos[0];	/* premult,delay and temp*7 */
					table[31] = sincos[1];

					table += 32;
				}
			}
			pass1_size /= 8;
		}

/* Output the sin/cos/premultiplier values for the radix-4 block that does the */
/* last 2 levels in pass 1. */

		else {
			N = N * 4;

/* Output the extra sin/cos values needed for the eight_reals FFT work done */
/* on the last pass 1 level.  We double N because the real part of the FFT */
/* is one level behind the complex part of the FFT. */

			if (!gwdata->ALL_COMPLEX_FFT) {
				for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
					temp = (group + i);
					gwsincos5 (temp, N*2, (double *) &sincos);
					table[0] = sincos[0];	/* temp for 2*N */
					table[2] = sincos[1];
					table[4] = sincos[2];	/* temp * 5 for 2*N */
					table[6] = sincos[3];

					temp += upper_sse2_word;
					gwsincos5 (temp, N*2, (double *) &sincos);
					table[1] = sincos[0];	/* temp for 2*N */
					table[3] = sincos[1];
					table[5] = sincos[2];	/* temp * 5 for 2*N */
					table[7] = sincos[3];

					table += 8;
				}
			}

/* Output the sin/cos values for the complex delay groups -- specifically the r4_sg4cl_eight_reals_four_complex_fft4 macro. */

			for (k = 0; k < delay_count; k++) {
				for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
					unsigned long bigN, ktemp, actemp;

/* If this is an all-complex FFT, the roots of minus 1 (same as roots of FFTLEN*2) are */
/* split to reduce memory requirements.  We apply the part of the all-complex premultiplier here. */

					if (gwdata->ALL_COMPLEX_FFT) {
						bigN = gwdata->FFTLEN * 2;
						actemp = group + i;
					} else {
						bigN = gwdata->FFTLEN;
						actemp = 0;
					}

/* Factor in the delayed part of the sin/cos multiplies from the first 2 levels.  In the first 2 levels */
/* we use a fixed sin/cos table based only on j, leaving the group+i part to be applied here by */
/* creating delay_count table entries.  For an all-complex FFT we multiply by 0,2,1,-1*(group+i) with N = FFTLEN/2. */
/* For an all-real FFT we multiply by 0,2,1,5*(group+i) with N = FFTLEN. */

					if (gwdata->ALL_COMPLEX_FFT && delay_count == 4) {
						if (k == 0)
							ktemp = 0;
						else if (k == 1)
							ktemp = 2 * (group + i) * 4;
						else if (k == 2)
							ktemp = 1 * (group + i) * 4;
						else
							ktemp = bigN - 1 * (group + i) * 4;
					} else if (gwdata->ALL_COMPLEX_FFT) {
						/* 0,2,1,-1 combined with 0,2,1,-1 */
							int	kmap[16] = {0,8,4,-4,2,10,6,-2,1,9,5,-3,-1,7,3,-5};
							if (kmap[k] >= 0)
								ktemp = kmap[k] * (group + i) * 4;
							else
								ktemp = bigN + kmap[k] * (group + i) * 4;
					} else if (delay_count == 4) {
						if (k == 0)
							ktemp = 0;
						else if (k == 1)
							ktemp = 2 * (group + i);
						else if (k == 2)
							ktemp = 1 * (group + i);
						else
							ktemp = 5 * (group + i);
					} else if (delay_count == 16) {
						/* 0,2,1,5 combined with one 0,2,1,5 and three 0,2,1,-1 */
						int	kmap[16] = {0,8,4,20,2,18,10,-6,1,17,9,-7,5,21,13,-3};
						if (kmap[k] >= 0)
							ktemp = kmap[k] * (group + i);
						else
							ktemp = bigN + kmap[k] * (group + i);
					} else {			/* delay_count == 10 or 14 */
						/* Multipliers for the radix-20 or radix-28 step */
						ktemp = k * (group + i);
					}

/* We now calculate the standard sin/cos twiddle factors (temp*0,1,2,3 roots of N) */
/* combined with the group all-complex premultiplier roots of minus 1 (same as roots of FFTLEN*2) */
/* combined with the delayed group multipliers. */

					temp = (group + i) * (bigN / N);
					gwsincos (actemp + ktemp, bigN, (double *) &sincos);
					table[0] = sincos[0];	/* premult,delay and temp*0 */
					table[2] = sincos[1];
					gwsincos (actemp + ktemp + temp, bigN, (double *) &sincos);
					table[4] = sincos[0];	/* premult,delay and temp*1 */
					table[6] = sincos[1];
					gwsincos (actemp + ktemp + temp * 2, bigN, (double *) &sincos);
					table[8] = sincos[0];	/* premult,delay and temp*2 */
					table[10] = sincos[1];
					gwsincos (actemp + ktemp + temp * 3, bigN, (double *) &sincos);
					table[12] = sincos[0];	/* premult,delay and temp*3 */
					table[14] = sincos[1];

					temp = (group + i + upper_sse2_word) * (bigN / N);
					if (gwdata->ALL_COMPLEX_FFT) actemp += upper_sse2_word;
					gwsincos (actemp + ktemp, bigN, (double *) &sincos);
					table[1] = sincos[0];	/* premult,delay and temp*0 */
					table[3] = sincos[1];
					gwsincos (actemp + ktemp + temp, bigN, (double *) &sincos);
					table[5] = sincos[0];	/* premult,delay and temp*1 */
					table[7] = sincos[1];
					gwsincos (actemp + ktemp + temp * 2, bigN, (double *) &sincos);
					table[9] = sincos[0];	/* premult,delay and temp*2 */
					table[11] = sincos[1];
					gwsincos (actemp + ktemp + temp * 3, bigN, (double *) &sincos);
					table[13] = sincos[0];	/* premult,delay and temp*3 */
					table[15] = sincos[1];

					table += 16;
				}
			}
			pass1_size /= 4;
		}

/* For the r4_four_complex_djbfft building block levels, output the sin/cos values. */

		while ((pass1_size & 3) == 0) {
			N = N * 4;

/* For the r4_eight_reals_four_complex_djbfft building block levels, output the extra */
/* sin/cos values needed for the eight_reals.  The eight_reals doubles N because */
/* the real part of the FFT is one level behind the complex part of the FFT. */

			if (!gwdata->ALL_COMPLEX_FFT) {
				for (j = 0; j < N*2 / 8; j += pass1_increment) {
					for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
						temp = (group + j + i);
						gwsincos5 (temp, N*2, (double *) &sincos);
						table[0] = sincos[0];	/* temp */
						table[2] = sincos[1];
						table[4] = sincos[2];	/* 5 * temp */
						table[6] = sincos[3];

						temp += upper_sse2_word;
						gwsincos5 (temp, N*2, (double *) &sincos);
						table[1] = sincos[0];	/* temp */
						table[3] = sincos[1];
						table[5] = sincos[2];	/* 5 * temp */
						table[7] = sincos[3];

						table += 8;
					}
				}
			}

/* Output the sin/cos value for the complex sections, used by the r4_four_complex_djbfft macro */

			for (j = 0; j < N / 4; j += pass1_increment) {
				for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
					temp = (group + j + i);
					gwsincos2 (temp, N, (double *) &sincos);
					table[0] = sincos[0];	/* temp */
					table[2] = sincos[1];
					table[4] = sincos[2];	/* 2 * temp */
					table[6] = sincos[3];

					gwsincos2 (temp + upper_sse2_word, N, (double *) &sincos);
					table[1] = sincos[0];	/* temp */
					table[3] = sincos[1];
					table[5] = sincos[2];	/* 2 * temp */
					table[7] = sincos[3];

					table += 8;
				}
			}
			pass1_size /= 4;
		}

/* For the r5_five_complex_djbfft building block levels, output the sin/cos values. */

		while (pass1_size % 5 == 0) {
			N = N * 5;

/* The r5_ten_reals building blocks require extra sin/cos */
/* values.  The ten_reals doubles N because the real part of the */
/* FFT is one level behind the complex part of the FFT. */

			if (!gwdata->ALL_COMPLEX_FFT) {
				for (j = 0; j < N / 5; j += pass1_increment) {
					for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
						temp = (group + j + i);
						gwsincos3 (temp, N*2, (double *) &sincos);
						table[0] = sincos[0];	/* temp */
						table[2] = sincos[1];
						table[4] = sincos[4];	/* 3 * temp */
						table[6] = sincos[5];

						temp += upper_sse2_word;
						gwsincos3 (temp, N*2, (double *) &sincos);
						table[1] = sincos[0];	/* temp */
						table[3] = sincos[1];
						table[5] = sincos[4];	/* 3 * temp */
						table[7] = sincos[5];

						table += 8;
					}
				}
			}

/* Output the sin/cos data for the complex sections, (the r5_five_complex_djbfft building block). */

			for (j = 0; j < N / 5; j += pass1_increment) {
				for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
					temp = (group + j + i);
					gwsincos2 (temp, N, (double *) &sincos);
					table[0] = sincos[0];	/* temp */
					table[2] = sincos[1];
					table[4] = sincos[2];	/* 2 * temp */
					table[6] = sincos[3];

					gwsincos2 (temp + upper_sse2_word, N, (double *) &sincos);
					table[1] = sincos[0];	/* temp */
					table[3] = sincos[1];
					table[5] = sincos[2];	/* 2 * temp */
					table[7] = sincos[3];

					table += 8;
				}
			}
			pass1_size /= 5;
		}

/* For the r3_three_complex_djbfft building block levels, output the sin/cos values. */

		while (pass1_size % 3 == 0) {
			N = N * 3;

/* The r3_six_reals building blocks require an extra sin/cos */
/* value.  The six_reals doubles N because the real part of the */
/* FFT is one level behind the complex part of the FFT. */

			if (!gwdata->ALL_COMPLEX_FFT) {
				for (j = 0; j < N / 3; j += pass1_increment) {
					for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
						temp = (group + j + i);
						gwsincos (temp, N*2, (double *) &sincos);
						table[0] = sincos[0];	/* temp */
						table[2] = sincos[1];

						temp += upper_sse2_word;
						gwsincos (temp, N*2, (double *) &sincos);
						table[1] = sincos[0];	/* temp */
						table[3] = sincos[1];

						table += 4;
					}
				}
			}

/* Output the sin/cos data for the complex sections (used by the r3_three_complex_djbfft building block). */

			for (j = 0; j < N / 3; j += pass1_increment) {
				for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
					temp = (group + j + i);
					gwsincos (temp, N, (double *) &sincos);
					table[0] = sincos[0];	/* temp */
					table[2] = sincos[1];

					gwsincos (temp + upper_sse2_word, N, (double *) &sincos);
					table[1] = sincos[0];	/* temp */
					table[3] = sincos[1];

					table += 4;
				}
			}
			pass1_size /= 3;
		}
		ASSERTG (pass1_size == 1);

/* Calculate the size of each group's sin/cos/premult data for pass1_get_next_block */

		if (group == 0)
			gwdata->pass1_premult_block_size = (unsigned long)
				((char *) table - (char *) gwdata->adjusted_pass1_premults) /
				gwdata->PASS1_CACHE_LINES;
	}

/* Return address of the end of the table */

	return (table);
}

/* This routine builds the fixed postmultiplier table used in pass 1 of the */
/* radix-4/8 delayed DJB FFT - called by gwsetup. */

double *r4delay_build_fixed_pass1_table (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	double	*table)		/* Pointer to the table to fill in */
{
	unsigned long pass1_increment, upper_sse2_word, pass1_size, i, j, temp, N;
	double	sincos[6];

/* Initialize some needed constants */

	pass1_increment = gwdata->PASS2_SIZE;
	upper_sse2_word = pass1_increment / 2;
	pass1_size = gwdata->FFTLEN / gwdata->PASS2_SIZE; /* Real values in a pass1 section */

/* Real FFTs output one shared set of sin/cos values for the first 20-reals, 28-reals, or 8-reals FFT. */

	if (! gwdata->ALL_COMPLEX_FFT) {
		N = gwdata->FFTLEN;
		if (pass1_size % 20 == 0) {
			for (j = 0; j < N / 20; j += pass1_increment) {
				for (i = 1; i <= 9; i++) {	/* Create 9 twiddle factors */
					gwsincos (i * j, N, (double *) &sincos);
					table[0] = sincos[0];	/* i * j */
					table[2] = sincos[1];
					gwsincos (i * (j + upper_sse2_word), N, (double *) &sincos);
					table[1] = sincos[0];
					table[3] = sincos[1];
					table += 4;
				}
			}
			N = N / 20;
		}
		else if (pass1_size % 28 == 0) {
			for (j = 0; j < N / 28; j += pass1_increment) {
				for (i = 1; i <= 13; i++) {	/* Create 13 twiddle factors */
					gwsincos (i * j, N, (double *) &sincos);
					table[0] = sincos[0];	/* i * j */
					table[2] = sincos[1];
					gwsincos (i * (j + upper_sse2_word), N, (double *) &sincos);
					table[1] = sincos[0];
					table[3] = sincos[1];
					table += 4;
				}
			}
			N = N / 28;
		}
		else {
			for (j = 0; j < N / 8; j += pass1_increment) {
				gwsincos25 (j, N, (double *) &sincos);
				table[0] = sincos[0];	/* temp */
				table[2] = sincos[1];
				table[4] = sincos[2];	/* 2 * temp */
				table[6] = sincos[3];
				table[8] = sincos[4];	/* 5 * temp */
				table[10] = sincos[5];

				temp = j + upper_sse2_word;
				gwsincos25 (temp, N, (double *) &sincos);
				table[1] = sincos[0];	/* temp */
				table[3] = sincos[1];
				table[5] = sincos[2];	/* 2 * temp */
				table[7] = sincos[3];
				table[9] = sincos[4];	/* 5 * temp */
				table[11] = sincos[5];

				table += 12;
			}
			N = N / 8;
			/* Sometimes we also use a fixed sin/cos table for */
			/* the next FFT levels to further reduce memory usage. */
			if (pass1_size == 512 || pass1_size == 1024 || pass1_size == 1536 ||
			    pass1_size == 2048 || pass1_size == 3072 || pass1_size == 4096) {
				/* Output the sin/cos values for the real data */
				for (j = 0; j < N*2 / 8; j += pass1_increment) {
					gwsincos5 (j, N*2, (double *) &sincos);
					table[0] = sincos[0];	/* temp */
					table[2] = sincos[1];
					table[4] = sincos[2];	/* 5 * temp */
					table[6] = sincos[3];

					temp = j + upper_sse2_word;
					gwsincos5 (temp, N*2, (double *) &sincos);
					table[1] = sincos[0];	/* temp */
					table[3] = sincos[1];
					table[5] = sincos[2];	/* 5 * temp */
					table[7] = sincos[3];

					table += 8;
				}
				/* Output the sin/cos values for the complex data */
				for (j = 0; j < N / 4; j += pass1_increment) {
					gwsincos2 (j, N, (double *) &sincos);
					table[0] = sincos[0];	/* temp */
					table[2] = sincos[1];
					table[4] = sincos[2];	/* 2 * temp */
					table[6] = sincos[3];

					temp = j + upper_sse2_word;
					gwsincos2 (temp, N, (double *) &sincos);
					table[1] = sincos[0];	/* temp */
					table[3] = sincos[1];
					table[5] = sincos[2];	/* 2 * temp */
					table[7] = sincos[3];

					table += 8;
				}
				N = N / 4;
			}
		}
	}

/* For all-complex FFTs, also build the fixed roots-of-minus-one table. */
/* Output these roots in the same order they will be used in the first two */
/* levels of pass 1. */

	else {
		N = gwdata->FFTLEN / 2;
		for (j = 0; j < N / 4; j += pass1_increment) {
			gwsincos2 (j, N, (double *) &sincos);
			table[0] = sincos[0];	/* temp */
			table[2] = sincos[1];
			table[4] = sincos[2];	/* 2 * temp */
			table[6] = sincos[3];

			gwsincos2 (j + upper_sse2_word, N, (double *) &sincos);
			table[1] = sincos[0];	/* temp */
			table[3] = sincos[1];
			table[5] = sincos[2];	/* 2 * temp */
			table[7] = sincos[3];

			table += 8;
		}
		N = N / 4;
		/* Sometimes we also use a fixed sin/cos table for */
		/* the next FFT levels to further reduce memory usage. */
		if (pass1_size == 1536 || pass1_size == 2048 || pass1_size == 2560 ||
		    pass1_size == 3072 || pass1_size == 4096 || pass1_size == 5120) {
			/* Output the sin/cos values for the complex data */
			for (j = 0; j < N / 4; j += pass1_increment) {
				gwsincos2 (j, N, (double *) &sincos);
				table[0] = sincos[0];	/* temp */
				table[2] = sincos[1];
				table[4] = sincos[2];	/* 2 * temp */
				table[6] = sincos[3];

				temp = j + upper_sse2_word;
				gwsincos2 (temp, N, (double *) &sincos);
				table[1] = sincos[0];	/* temp */
				table[3] = sincos[1];
				table[5] = sincos[2];	/* 2 * temp */
				table[7] = sincos[3];

				table += 8;
			}
			N = N / 4;
		}
	}

/* Return address of the end of the table */

	return (table);
}

/* This routine builds the sin/cos table used in pass 1 by the radix-4/8 DJB */
/* FFT with delayed sin/cos multiplies and with partial normalization. */

double *r4dwpn_build_pass1_table (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	double	*table)		/* Pointer to the table to fill in */
{
	struct gwasm_data *asm_data = (struct gwasm_data *) gwdata->asm_data;
	unsigned long pass1_size, pass1_increment, delay_count;
	unsigned long group, i, j, k, N, temp, upper_sse2_word;
	int	pow2_count;
	double	sincos[6];

/* Initialize some needed constants */

	pass1_increment = gwdata->PASS2_SIZE;
	upper_sse2_word = pass1_increment / 2;
	pass1_size = gwdata->FFTLEN / gwdata->PASS2_SIZE; /* Real values in a pass1 section */

/* Determine number of delay groups.  In a standard radix-4 FFT, there is only one sin/cos */
/* group in the last pass 1 level.  We reduce our memory usage by using just one fixed sin/cos */
/* table in the first FFT levels and having multiple groups of sin/cos data in the last pass 1 level. */
/* I call these groups of sin/cos data in the last pass 1 level "delay groups". */

	if (pass1_size % 7 == 0)
		delay_count = 14;
	else if (pass1_size % 5 == 0 && !gwdata->ALL_COMPLEX_FFT)
		delay_count = 10;
 	else if ((pass1_size == 512 && !gwdata->ALL_COMPLEX_FFT) ||
		 (pass1_size == 1024 && !gwdata->ALL_COMPLEX_FFT) ||
		 (pass1_size == 2560 && gwdata->ALL_COMPLEX_FFT) ||
		 (pass1_size == 5120 && gwdata->ALL_COMPLEX_FFT) ||
		 pass1_size == 1536 || pass1_size == 2048 || pass1_size == 3072 || pass1_size == 4096)
		delay_count = 16;
	else
		delay_count = 4;

/* Count the power-of-two FFT levels after the initial FFT levels.  If odd, the */
/* the last levels will be a radix-8, if even the last levels wil be a radix-4. */

	pass1_size /= (delay_count * 2);
	for (pow2_count = 0; (pass1_size & 1) == 0; pass1_size /= 2) pow2_count++;

/* Set count of pass 1 blocks share one set of two-to-phi grp multipliers */

	if (pow2_count & 1) gwdata->wpn_count = 8;
	else if (gwdata->FFTLEN / gwdata->PASS2_SIZE == 1792 ||
		 gwdata->FFTLEN / gwdata->PASS2_SIZE == 2048) gwdata->wpn_count = 16;
	else gwdata->wpn_count = 4;
	asm_data->count2 = gwdata->wpn_count;
	asm_data->count3 = asm_data->addcount1 / gwdata->wpn_count;

/* Set pointer to table of multipliers */

	gwdata->adjusted_pass1_premults = table;

/* Loop through all the pass 1 groups in the same order the assembly code will */
/* process the groups. */

	for (group = 0; group < upper_sse2_word; group += gwdata->PASS1_CACHE_LINES) {

		pass1_size = gwdata->FFTLEN / gwdata->PASS2_SIZE; /* Real values in a pass1 section */
		pass1_size /= (delay_count * 2);	/* Complex values we're generating sin/cos data for */
		N = gwdata->PASS2_SIZE;

/* Output the sin/cos/premultiplier values for the radix-8 block that does the */
/* last 3 levels in pass 1.  NOTE:  We do not need the "j loop" (it would loop */
/* from zero to zero) when generating the sin/cos twiddle factors for the last */
/* levels of pass 1. */

		if (pow2_count & 1) {
			N = N * 8;

/* For the r4_g8cl_sixteen_reals_eight_complex_djbfft building block, output the extra */
/* sin/cos values needed for the sixteen_reals.  The sixteen_reals is done in three parts: */
/* 4 four_reals_fft, 1 eight_reals_fft, and 1 four-complex_fft.  The four_reals_fft needs */
/* four sin/cos twiddles using 2*N, the eight_reals_fft and four_complex_fft can use the same */
/* sin/cos twiddles that will be generated later for the r4r8_g8cl_eight_complex_fft8 macro. */

			if (!gwdata->ALL_COMPLEX_FFT) {
				for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
					temp = (group + i);
					gwsincos (temp, N*2, (double *) &sincos);
					table[0] = sincos[0];	/* temp for 2*N */
					table[2] = sincos[1];
					gwsincos (temp + pass1_increment, N*2, (double *) &sincos);
					table[4] = sincos[0];	/* temp + pass1_increment for 2*N */
					table[6] = sincos[1];
					gwsincos (temp + 2 * pass1_increment, N*2, (double *) &sincos);
					table[8] = sincos[0];	/* temp + 2 * pass1_increment for 2*N */
					table[10] = sincos[1];
					gwsincos (temp + 3 * pass1_increment, N*2, (double *) &sincos);
					table[12] = sincos[0];	/* temp + 3 * pass1_increment for 2*N */
					table[14] = sincos[1];

					temp += upper_sse2_word;
					gwsincos (temp, N*2, (double *) &sincos);
					table[1] = sincos[0];	/* temp for 2*N */
					table[3] = sincos[1];
					gwsincos (temp + pass1_increment, N*2, (double *) &sincos);
					table[5] = sincos[0];	/* temp + pass1_increment for 2*N */
					table[7] = sincos[1];
					gwsincos (temp + 2 * pass1_increment, N*2, (double *) &sincos);
					table[9] = sincos[0];	/* temp + 2 * pass1_increment for 2*N */
					table[11] = sincos[1];
					gwsincos (temp + 3 * pass1_increment, N*2, (double *) &sincos);
					table[13] = sincos[0];	/* temp + 3 * pass1_increment for 2*N */
					table[15] = sincos[1];

					table += 16;
				}
			}

/* Output the sin/cos values for the complex delay groups -- specifically the r8_sg8cl_eight_complex_fft8 macro. */

			for (k = 0; k < delay_count; k++) {
				for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
					unsigned long bigN, ktemp, actemp;

/* If this is an all-complex FFT, the roots of minus 1 (same as roots of FFTLEN*2) are */
/* split to reduce memory requirements.  We apply the part of the all-complex premultiplier here. */

					if (gwdata->ALL_COMPLEX_FFT) {
						bigN = gwdata->FFTLEN * 2;
						actemp = group + i;
					} else {
						bigN = gwdata->FFTLEN;
						actemp = 0;
					}

/* Factor in the delayed part of the sin/cos multiplies from the first 2 levels.  In the first 2 levels */
/* we use a fixed sin/cos table based only on j, leaving the group+i part to be applied here by */
/* creating delay_count table entries.  For an all-complex FFT we multiply by 0,2,1,-1*(group+i) with N = FFTLEN/2. */
/* For an all-real FFT we multiply by 0,2,1,5*(group+i) with N = FFTLEN. */

					if (gwdata->ALL_COMPLEX_FFT && delay_count == 4) {
						if (k == 0)
							ktemp = 0;
						else if (k == 1)
							ktemp = 2 * (group + i) * 4;
						else if (k == 2)
							ktemp = 1 * (group + i) * 4;
						else
							ktemp = bigN - 1 * (group + i) * 4;
					} else if (gwdata->ALL_COMPLEX_FFT) {
						/* 0,2,1,-1 combined with 0,2,1,-1 */
						int	kmap[16] = {0,8,4,-4,2,10,6,-2,1,9,5,-3,-1,7,3,-5};
						if (kmap[k] >= 0)
							ktemp = kmap[k] * (group + i) * 4;
						else
							ktemp = bigN + kmap[k] * (group + i) * 4;
					} else if (delay_count == 4) {
						if (k == 0)
							ktemp = 0;
						else if (k == 1)
							ktemp = 2 * (group + i);
						else if (k == 2)
							ktemp = 1 * (group + i);
						else
							ktemp = 5 * (group + i);
					} else if (delay_count == 16) {
						/* 0,2,1,5 combined with one 0,2,1,5 and three 0,2,1,-1 */
						int	kmap[16] = {0,8,4,20,2,18,10,-6,1,17,9,-7,5,21,13,-3};
						if (kmap[k] >= 0)
							ktemp = kmap[k] * (group + i);
						else
							ktemp = bigN + kmap[k] * (group + i);
					} else {			/* delay_count == 10 or 14 */
						/* Multipliers for the radix-20 or radix-28 step */
						ktemp = k * (group + i);
					}

/* We now calculate the standard sin/cos twiddle factors (temp*0,1,2,3,4,5,6,7 roots of N) */
/* combined with the group all-complex premultiplier roots of minus 1 (same as roots of FFTLEN*2) */
/* combined with the delayed group multipliers. */

					temp = (group + i) * (bigN / N);
					gwsincos (actemp + ktemp, bigN, (double *) &sincos);
					table[0] = sincos[0];	/* premult,delay and temp*0 */
					table[2] = sincos[1];
					gwsincos (actemp + ktemp + temp, bigN, (double *) &sincos);
					table[4] = sincos[0];	/* premult,delay and temp*1 */
					table[6] = sincos[1];
					gwsincos (actemp + ktemp + temp * 2, bigN, (double *) &sincos);
					table[8] = sincos[0];	/* premult,delay and temp*2 */
					table[10] = sincos[1];
					gwsincos (actemp + ktemp + temp * 3, bigN, (double *) &sincos);
					table[12] = sincos[0];	/* premult,delay and temp*3 */
					table[14] = sincos[1];
					gwsincos (actemp + ktemp + temp * 4, bigN, (double *) &sincos);
					table[16] = sincos[0];	/* premult,delay and temp*4 */
					table[18] = sincos[1];
					gwsincos (actemp + ktemp + temp * 5, bigN, (double *) &sincos);
					table[20] = sincos[0];	/* premult,delay and temp*5 */
					table[22] = sincos[1];
					gwsincos (actemp + ktemp + temp * 6, bigN, (double *) &sincos);
					table[24] = sincos[0];	/* premult,delay and temp*6 */
					table[26] = sincos[1];
					gwsincos (actemp + ktemp + temp * 7, bigN, (double *) &sincos);
					table[28] = sincos[0];	/* premult,delay and temp*7 */
					table[30] = sincos[1];

					temp = (group + i + upper_sse2_word) * (bigN / N);
					if (gwdata->ALL_COMPLEX_FFT) actemp += upper_sse2_word;
					gwsincos (actemp + ktemp, bigN, (double *) &sincos);
					table[1] = sincos[0];	/* premult,delay and temp*0 */
					table[3] = sincos[1];
					gwsincos (actemp + ktemp + temp, bigN, (double *) &sincos);
					table[5] = sincos[0];	/* premult,delay and temp*1 */
					table[7] = sincos[1];
					gwsincos (actemp + ktemp + temp * 2, bigN, (double *) &sincos);
					table[9] = sincos[0];	/* premult,delay and temp*2 */
					table[11] = sincos[1];
					gwsincos (actemp + ktemp + temp * 3, bigN, (double *) &sincos);
					table[13] = sincos[0];	/* premult,delay and temp*3 */
					table[15] = sincos[1];
					gwsincos (actemp + ktemp + temp * 4, bigN, (double *) &sincos);
					table[17] = sincos[0];	/* premult,delay and temp*4 */
					table[19] = sincos[1];
					gwsincos (actemp + ktemp + temp * 5, bigN, (double *) &sincos);
					table[21] = sincos[0];	/* premult,delay and temp*5 */
					table[23] = sincos[1];
					gwsincos (actemp + ktemp + temp * 6, bigN, (double *) &sincos);
					table[25] = sincos[0];	/* premult,delay and temp*6 */
					table[27] = sincos[1];
					gwsincos (actemp + ktemp + temp * 7, bigN, (double *) &sincos);
					table[29] = sincos[0];	/* premult,delay and temp*7 */
					table[31] = sincos[1];

					table += 32;
				}
			}
			pass1_size /= 8;
		}

/* Output the sin/cos/premultiplier values for the radix-4 block that does the */
/* last 2 levels in pass 1. */

		else {
			N = N * 4;

/* Output the extra sin/cos values needed for the eight_reals FFT work done */
/* on the last pass 1 level.  We double N because the real part of the FFT */
/* is one level behind the complex part of the FFT. */

			if (!gwdata->ALL_COMPLEX_FFT) {
				for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
					temp = (group + i);
					gwsincos5 (temp, N*2, (double *) &sincos);
					table[0] = sincos[0];	/* temp for 2*N */
					table[2] = sincos[1];
					table[4] = sincos[2];	/* temp * 5 for 2*N */
					table[6] = sincos[3];

					temp += upper_sse2_word;
					gwsincos5 (temp, N*2, (double *) &sincos);
					table[1] = sincos[0];	/* temp for 2*N */
					table[3] = sincos[1];
					table[5] = sincos[2];	/* temp * 5 for 2*N */
					table[7] = sincos[3];

					table += 8;
				}
			}

/* Output the sin/cos values for the complex delay groups -- specifically the r4_sg4cl_eight_reals_four_complex_fft4 macro. */

			for (k = 0; k < delay_count; k++) {
				for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
					unsigned long bigN, ktemp, actemp;

/* If this is an all-complex FFT, the roots of minus 1 (same as roots of FFTLEN*2) are */
/* split to reduce memory requirements.  We apply the part of the all-complex premultiplier here. */

					if (gwdata->ALL_COMPLEX_FFT) {
						bigN = gwdata->FFTLEN * 2;
						actemp = group + i;
					} else {
						bigN = gwdata->FFTLEN;
						actemp = 0;
					}

/* Factor in the delayed part of the sin/cos multiplies from the first 2 levels.  In the first 2 levels */
/* we use a fixed sin/cos table based only on j, leaving the group+i part to be applied here by */
/* creating delay_count table entries.  For an all-complex FFT we multiply by 0,2,1,-1*(group+i) with N = FFTLEN/2. */
/* For an all-real FFT we multiply by 0,2,1,5*(group+i) with N = FFTLEN. */

					if (gwdata->ALL_COMPLEX_FFT && delay_count == 4) {
						if (k == 0)
							ktemp = 0;
						else if (k == 1)
							ktemp = 2 * (group + i) * 4;
						else if (k == 2)
							ktemp = 1 * (group + i) * 4;
						else
							ktemp = bigN - 1 * (group + i) * 4;
					} else if (gwdata->ALL_COMPLEX_FFT) {
						/* 0,2,1,-1 combined with 0,2,1,-1 */
							int	kmap[16] = {0,8,4,-4,2,10,6,-2,1,9,5,-3,-1,7,3,-5};
							if (kmap[k] >= 0)
								ktemp = kmap[k] * (group + i) * 4;
							else
								ktemp = bigN + kmap[k] * (group + i) * 4;
					} else if (delay_count == 4) {
						if (k == 0)
							ktemp = 0;
						else if (k == 1)
							ktemp = 2 * (group + i);
						else if (k == 2)
							ktemp = 1 * (group + i);
						else
							ktemp = 5 * (group + i);
					} else if (delay_count == 16) {
						/* 0,2,1,5 combined with one 0,2,1,5 and three 0,2,1,-1 */
						int	kmap[16] = {0,8,4,20,2,18,10,-6,1,17,9,-7,5,21,13,-3};
						if (kmap[k] >= 0)
							ktemp = kmap[k] * (group + i);
						else
							ktemp = bigN + kmap[k] * (group + i);
					} else {			/* delay_count == 10 or 14 */
						/* Multipliers for the radix-20 or radix-28 step */
						ktemp = k * (group + i);
					}

/* We now calculate the standard sin/cos twiddle factors (temp*0,1,2,3 roots of N) */
/* combined with the group all-complex premultiplier roots of minus 1 (same as roots of FFTLEN*2) */
/* combined with the delayed group multipliers. */

					temp = (group + i) * (bigN / N);
					gwsincos (actemp + ktemp, bigN, (double *) &sincos);
					table[0] = sincos[0];	/* premult,delay and temp*0 */
					table[2] = sincos[1];
					gwsincos (actemp + ktemp + temp, bigN, (double *) &sincos);
					table[4] = sincos[0];	/* premult,delay and temp*1 */
					table[6] = sincos[1];
					gwsincos (actemp + ktemp + temp * 2, bigN, (double *) &sincos);
					table[8] = sincos[0];	/* premult,delay and temp*2 */
					table[10] = sincos[1];
					gwsincos (actemp + ktemp + temp * 3, bigN, (double *) &sincos);
					table[12] = sincos[0];	/* premult,delay and temp*3 */
					table[14] = sincos[1];

					temp = (group + i + upper_sse2_word) * (bigN / N);
					if (gwdata->ALL_COMPLEX_FFT) actemp += upper_sse2_word;
					gwsincos (actemp + ktemp, bigN, (double *) &sincos);
					table[1] = sincos[0];	/* premult,delay and temp*0 */
					table[3] = sincos[1];
					gwsincos (actemp + ktemp + temp, bigN, (double *) &sincos);
					table[5] = sincos[0];	/* premult,delay and temp*1 */
					table[7] = sincos[1];
					gwsincos (actemp + ktemp + temp * 2, bigN, (double *) &sincos);
					table[9] = sincos[0];	/* premult,delay and temp*2 */
					table[11] = sincos[1];
					gwsincos (actemp + ktemp + temp * 3, bigN, (double *) &sincos);
					table[13] = sincos[0];	/* premult,delay and temp*3 */
					table[15] = sincos[1];

					table += 16;
				}
			}
			pass1_size /= 4;
		}

/* Output multipliers for the four complex building blocks. */

		while ((pass1_size & 3) == 0) {

			N = N * 4;

/* For the non-wpn levels, output the sin/cos values. */

			if (N != gwdata->PASS2_SIZE * gwdata->wpn_count * 4) {

/* For the r4_eight_reals_four_complex_djbfft building block levels, output the extra */
/* sin/cos values needed for the eight_reals.  The eight_reals doubles N because */
/* the real part of the FFT is one level behind the complex part of the FFT. */

				if (!gwdata->ALL_COMPLEX_FFT) {
					for (j = 0; j < N*2 / 8; j += pass1_increment) {
					    for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
						temp = (group + j + i);
						gwsincos5 (temp, N*2, (double *) &sincos);
						table[0] = sincos[0];	/* temp */
						table[2] = sincos[1];
						table[4] = sincos[2];	/* 5 * temp */
						table[6] = sincos[3];

						temp += upper_sse2_word;
						gwsincos5 (temp, N*2, (double *) &sincos);
						table[1] = sincos[0];	/* temp */
						table[3] = sincos[1];
						table[5] = sincos[2];	/* 5 * temp */
						table[7] = sincos[3];

						table += 8;
					    }
					}
				}

/* Output the sin/cos value for the complex sections, used by the r4_four_complex_djbfft macro */

				for (j = 0; j < N / 4; j += pass1_increment) {
				    for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
					temp = (group + j + i);
					gwsincos2 (temp, N, (double *) &sincos);
					table[0] = sincos[0];	/* temp */
					table[2] = sincos[1];
					table[4] = sincos[2];	/* 2 * temp */
					table[6] = sincos[3];

					gwsincos2 (temp + upper_sse2_word, N, (double *) &sincos);
					table[1] = sincos[0];	/* temp */
					table[3] = sincos[1];
					table[5] = sincos[2];	/* 2 * temp */
					table[7] = sincos[3];

					table += 8;
				    }
				}
			}

/* For the wpn building block level, output the sin/cos values. */

			else {

/* For the r4_wpn_eight_reals_four_complex_djbfft building block levels, output the extra */
/* sin/cos values needed for the eight_reals.  The eight_reals doubles N because */
/* the real part of the FFT is one level behind the complex part of the FFT. */

				if (!gwdata->ALL_COMPLEX_FFT) {
					for (j = 0; j < N*2 / 8; j += pass1_increment) {
					    for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
						temp = (group + j + i);
						gwsincos5_weighted (gwdata->dd_data, temp, N*2, temp, (double *) &sincos);
						table[0] = sincos[0];	/* temp */
						table[2] = sincos[1];
						table[4] = sincos[2];
						table[6] = sincos[3];	/* 5 * temp */
						table[8] = sincos[4];
						table[10] = sincos[5];

						gwsincos5_weighted (gwdata->dd_data, temp + upper_sse2_word, N*2, temp, (double *) &sincos);
						table[1] = sincos[0];	/* temp */
						table[3] = sincos[1];
						table[5] = sincos[2];
						table[7] = sincos[3];	/* 5 * temp */
						table[9] = sincos[4];
						table[11] = sincos[5];

						table += 12;
					    }
					}
				}

/* Output the sin/cos value for the complex sections, used by the r4_wpn_four_complex_djbfft macro */
/* We apply the two-to-phi weight for the upper SSE2 word in the group multipliers.  There is a */
/* reason for doing it there rather than here (it reduces the number of valid fudge factor combinations */
/* for each SSE2 word from 4 to 3). */

				for (j = 0; j < N / 4; j += pass1_increment) {
				    for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
					double	ttp, ttmp;

					temp = (group + j + i);
					gwfft_weights3 (gwdata->dd_data, temp, &ttp, NULL, &ttmp);
					table[0] = ttp;
					table[2] = ttmp;
					gwsincos2_weighted (gwdata->dd_data, temp, N, temp, (double *) &sincos);
					table[4] = sincos[0];	/* temp */
					table[6] = sincos[1];
					table[8] = sincos[2];
					table[10] = sincos[3];	/* 2 * temp */
					table[12] = sincos[4];
					table[14] = sincos[5];

					table[1] = ttp;
					table[3] = ttmp;
					gwsincos2_weighted (gwdata->dd_data, temp + upper_sse2_word, N, temp, (double *) &sincos);
					table[5] = sincos[0];	/* temp */
					table[7] = sincos[1];
					table[9] = sincos[2];
					table[11] = sincos[3];	/* 2 * temp */
					table[13] = sincos[4];
					table[15] = sincos[5];

					table += 16;
				    }
				}
			}

			pass1_size /= 4;
		}

/* For the r5_five_complex_djbfft building block levels, output the sin/cos values. */

		while (pass1_size % 5 == 0) {
			N = N * 5;

/* The r5_ten_reals building blocks require extra sin/cos */
/* values.  The ten_reals doubles N because the real part of the */
/* FFT is one level behind the complex part of the FFT. */

			if (!gwdata->ALL_COMPLEX_FFT) {
				for (j = 0; j < N / 5; j += pass1_increment) {
					for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
						temp = (group + j + i);
						gwsincos3 (temp, N*2, (double *) &sincos);
						table[0] = sincos[0];	/* temp */
						table[2] = sincos[1];
						table[4] = sincos[4];	/* 3 * temp */
						table[6] = sincos[5];

						temp += upper_sse2_word;
						gwsincos3 (temp, N*2, (double *) &sincos);
						table[1] = sincos[0];	/* temp */
						table[3] = sincos[1];
						table[5] = sincos[4];	/* 3 * temp */
						table[7] = sincos[5];

						table += 8;
					}
				}
			}

/* Output the sin/cos data for the complex sections, (the r5_five_complex_djbfft building block). */

			for (j = 0; j < N / 5; j += pass1_increment) {
				for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
					temp = (group + j + i);
					gwsincos2 (temp, N, (double *) &sincos);
					table[0] = sincos[0];	/* temp */
					table[2] = sincos[1];
					table[4] = sincos[2];	/* 2 * temp */
					table[6] = sincos[3];

					gwsincos2 (temp + upper_sse2_word, N, (double *) &sincos);
					table[1] = sincos[0];	/* temp */
					table[3] = sincos[1];
					table[5] = sincos[2];	/* 2 * temp */
					table[7] = sincos[3];

					table += 8;
				}
			}
			pass1_size /= 5;
		}

/* For the r3_three_complex_djbfft building block levels, output the sin/cos values. */

		while (pass1_size % 3 == 0) {
			N = N * 3;

/* The r3_six_reals building blocks require an extra sin/cos */
/* value.  The six_reals doubles N because the real part of the */
/* FFT is one level behind the complex part of the FFT. */

			if (!gwdata->ALL_COMPLEX_FFT) {
				for (j = 0; j < N / 3; j += pass1_increment) {
					for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
						temp = (group + j + i);
						gwsincos (temp, N*2, (double *) &sincos);
						table[0] = sincos[0];	/* temp */
						table[2] = sincos[1];

						temp += upper_sse2_word;
						gwsincos (temp, N*2, (double *) &sincos);
						table[1] = sincos[0];	/* temp */
						table[3] = sincos[1];

						table += 4;
					}
				}
			}

/* Output the sin/cos data for the complex sections (used by the r3_three_complex_djbfft building block). */

			for (j = 0; j < N / 3; j += pass1_increment) {
				for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
					temp = (group + j + i);
					gwsincos (temp, N, (double *) &sincos);
					table[0] = sincos[0];	/* temp */
					table[2] = sincos[1];

					gwsincos (temp + upper_sse2_word, N, (double *) &sincos);
					table[1] = sincos[0];	/* temp */
					table[3] = sincos[1];

					table += 4;
				}
			}
			pass1_size /= 3;
		}

/* Calculate the size of each group's sin/cos/premult data for pass1_get_next_block */

		if (group == 0)
			gwdata->pass1_premult_block_size = (unsigned long)
				((char *) table - (char *) gwdata->adjusted_pass1_premults) /
				gwdata->PASS1_CACHE_LINES;
	}

/* Return address of the end of the table */

	return (table);
}

/* This routine builds a normalization table - used by radix-4 normalizaion */
/* routines.  It differs from the home-grown SSE2 routines in that the different */
/* memory layout for PFA makes this routine much simpler. */

double *r4_build_norm_table (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	double	*table,		/* Pointer to the table to fill in */
	int	col)		/* TRUE if building column, not group, table */
{
	unsigned long i, k, num_cols;

/* Handle one-pass FFTs first, there are no group multipliers */

	if (gwdata->PASS2_SIZE == 0) {
		if (!col) return (table);

/* Loop to build table */

		for (i = 0; i < gwdata->FFTLEN; i++) {
			unsigned long j, table_entry;
			double	ttp, ttmp;

/* Call double-precision routine to compute the two multipliers */

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

	num_cols = gwdata->PASS2_SIZE / 2;
	if (col) {

/* Loop to build table */

		for (i = 0; i < num_cols; i++) {
			double	ttp, ttmp;

/* Call double-precision routine to compute the two multipliers */

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
		unsigned long h, hlimit, haddin, m, mmult, u, umult;

/* Loop to build table */

		umult = gwdata->FFTLEN / 2;
		hlimit = gwdata->FFTLEN / 4 / (2*num_cols);
		for (h = 0; h < hlimit; h++) {
			haddin = h * 2 * num_cols;
			mmult = gwdata->FFTLEN / 4;
			for (u = 0; u < 2; u++) {
				for (m = 0; m < 2; m++) {
					for (k = 0; k < 2; k++) {
						double	ttp, ttmp;
						long	n;

/* Call double-precision routine to compute the two multipliers */

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

/* This routine builds a normalization table - used by radix-4 with partial */
/* normalization FFTs. */

double *r4dwpn_build_norm_table (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	double	*table,		/* Pointer to the table to fill in */
	int	col)		/* TRUE if building column, not group, table */
{
	unsigned long i, j, num_cols, upper_sse2_word, second_word_in_cache_line;

/* Build the group multipliers table */

	upper_sse2_word = gwdata->PASS2_SIZE / 2;
	second_word_in_cache_line = gwdata->FFTLEN / 4;
	num_cols = gwdata->PASS2_SIZE * gwdata->wpn_count;
	for (i = 0; i < gwdata->FFTLEN / 4; i += num_cols) {
		for (j = 0; j < gwdata->FFTLEN; j += gwdata->FFTLEN / 2) {
			double	word1_lower_weight, word1_upper_weight, word2_lower_weight, word2_upper_weight;
			int	word1_lower_sort, word1_upper_sort, word2_lower_sort, word2_upper_sort;
			double	ttp, ttmp, ttp_over_b, ttmp_times_b, *tab10;

/* The sort order of the weights determines which fudge factor combination never occur. */

			word1_lower_weight = gwfft_weight_exponent (gwdata->dd_data, i + j);
			word1_upper_weight = gwfft_weight_exponent (gwdata->dd_data, i + j + upper_sse2_word);
			word2_lower_weight = gwfft_weight_exponent (gwdata->dd_data, i + j + second_word_in_cache_line);
			word2_upper_weight = gwfft_weight_exponent (gwdata->dd_data, i + j + second_word_in_cache_line + upper_sse2_word);

/* Now sort the four weights */

			word1_lower_sort = (word1_lower_weight > word1_upper_weight) +
					   (word1_lower_weight > word2_lower_weight) +
					   (word1_lower_weight > word2_upper_weight);
			word1_upper_sort = (word1_upper_weight > word1_lower_weight) +
					   (word1_upper_weight > word2_lower_weight) +
					   (word1_upper_weight > word2_upper_weight);
			word2_lower_sort = (word2_lower_weight > word1_lower_weight) +
					   (word2_lower_weight > word1_upper_weight) +
					   (word2_lower_weight > word2_upper_weight);
			word2_upper_sort = (word2_upper_weight > word1_lower_weight) +
					   (word2_upper_weight > word1_upper_weight) +
					   (word2_upper_weight > word2_lower_weight);

/* Call double-precision routine to compute first set of multipliers. */
/* We compute two-to-phi multiplied by the fudge factor so that normalize won't have to. */

			gwfft_weights_fudged (gwdata->dd_data, i + j, gwdata->b, &ttp, &ttmp, &ttp_over_b, &ttmp_times_b);

/* Set the LSW entries in an SSE2 pair */

			table[0] = ttmp;
			table[2] = word1_lower_sort < 3 ? ttmp : ttmp_times_b;
			table[4] = word1_lower_sort < 2 ? ttmp : ttmp_times_b;
			table[6] = word1_lower_sort < 1 ? ttmp : ttmp_times_b;
			table[8] = ttmp_times_b;
			tab10 = table + 10;
			tab10[0] = ttp;
			tab10[2] = word1_lower_sort < 3 ? ttp : ttp_over_b;
			tab10[4] = word1_lower_sort < 2 ? ttp : ttp_over_b;
			tab10[6] = word1_lower_sort < 1 ? ttp : ttp_over_b;
			tab10[8] = ttp_over_b;

/* Set the MSW entries in an SSE2 pair */

			gwfft_weights_fudged (gwdata->dd_data, i + j + upper_sse2_word, gwdata->b, &ttp, &ttmp, &ttp_over_b, &ttmp_times_b);

			table[1] = ttmp;
			table[3] = word1_upper_sort < 3 ? ttmp : ttmp_times_b;
			table[5] = word1_upper_sort < 2 ? ttmp : ttmp_times_b;
			table[7] = word1_upper_sort < 1 ? ttmp : ttmp_times_b;
			table[9] = ttmp_times_b;
			tab10[1] = ttp;
			tab10[3] = word1_upper_sort < 3 ? ttp : ttp_over_b;
			tab10[5] = word1_upper_sort < 2 ? ttp : ttp_over_b;
			tab10[7] = word1_upper_sort < 1 ? ttp : ttp_over_b;
			tab10[9] = ttp_over_b;

			table += 20;

/* Repeat for the second word in the cache line */

			gwfft_weights_fudged (gwdata->dd_data, i + j + second_word_in_cache_line, gwdata->b, &ttp, &ttmp, &ttp_over_b, &ttmp_times_b);

			table[0] = ttmp;
			table[2] = word2_lower_sort < 3 ? ttmp : ttmp_times_b;
			table[4] = word2_lower_sort < 2 ? ttmp : ttmp_times_b;
			table[6] = word2_lower_sort < 1 ? ttmp : ttmp_times_b;
			table[8] = ttmp_times_b;
			tab10 = table + 10;
			tab10[0] = ttp;
			tab10[2] = word2_lower_sort < 3 ? ttp : ttp_over_b;
			tab10[4] = word2_lower_sort < 2 ? ttp : ttp_over_b;
			tab10[6] = word2_lower_sort < 1 ? ttp : ttp_over_b;
			tab10[8] = ttp_over_b;

			gwfft_weights_fudged (gwdata->dd_data, i + j + second_word_in_cache_line + upper_sse2_word, gwdata->b, &ttp, &ttmp, &ttp_over_b, &ttmp_times_b);

			table[1] = ttmp;
			table[3] = word2_upper_sort < 3 ? ttmp : ttmp_times_b;
			table[5] = word2_upper_sort < 2 ? ttmp : ttmp_times_b;
			table[7] = word2_upper_sort < 1 ? ttmp : ttmp_times_b;
			table[9] = ttmp_times_b;
			tab10[1] = ttp;
			tab10[3] = word2_upper_sort < 3 ? ttp : ttp_over_b;
			tab10[5] = word2_upper_sort < 2 ? ttp : ttp_over_b;
			tab10[7] = word2_upper_sort < 1 ? ttp : ttp_over_b;
			tab10[9] = ttp_over_b;

			table += 20;
		}
	}
	return (table);
}

/* This routine builds a big/little flags table - used by radix-4 normalizaion */
/* routines.  It differs from the home-grown SSE2 routines in that the different */
/* memory layout for PFA makes this routine much simpler. */

double *r4_build_biglit_table (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	double	*table)		/* Pointer to the table to fill in */
{
	unsigned char *p;
	unsigned long h, i, j, k, m, u, gap;
	unsigned long hlimit, haddin, mmult, umult;

/* Handle one pass FFTs differently */

	if (gwdata->PASS2_SIZE == 0) {

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

/* Determine the gap between XMM high and low words */

	gap = gwdata->PASS2_SIZE / 2;

/* Loop to build table in exactly the same order that it will be */
/* used by the assembly code. */

	p = (unsigned char *) table;
	umult = gwdata->FFTLEN / 2;
	hlimit = gwdata->FFTLEN / 4 / (2*gap);
	for (i = 0; i < gap; i += gwdata->PASS1_CACHE_LINES) {
		for (h = 0; h < hlimit; h++) {
			haddin = h * 2 * gap;
			mmult = gwdata->FFTLEN / 4;
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
							    gwfft_weight_exponent (gwdata->dd_data, word % gap) +
							    gwfft_weight_exponent (gwdata->dd_data, word - word % gap)) {
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

/* This routine builds the big/little flags table for a r4dwpn (radix-4 with partial normalization) FFT */

double *r4dwpn_build_biglit_table (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	double	*table)		/* Pointer to the table to fill in */
{
const	int	START_OF_CHAIN = 0x8000;
const	int	END_OF_CHAIN = 0x4000;
const	int	CHAIN_3_COMMON = 0x2000;
const	int	CHAIN_2_COMMON = 0x1000;
const	int	CHAIN_1_COMMON = 0x0800;
	int	combos[16], next_combo[16];
	unsigned char *p;
	unsigned long upper_sse2_word, num_cols;
	unsigned long group, i, j, k, m, n, num_combos;

/* Big/lit flags form a very regular pattern.  For example, if there are 18.3 b's */
/* per FFT word then you get either a big word followed by two little words or a */
/* big  followed by three little words.  Here we determine which patterns of big/lit */
/* are possible in an xnorm_2d_wpn macro which processes 4 SSE2 words.  There are */
/* 16 possible valid combinations which can be represented by indexing into an array */
/* of 32 SSE2 values. */

/* Generate all 16 possible valid combinations of big/lit flags */

	upper_sse2_word = gwdata->PASS2_SIZE / 2;
	num_combos = 0;
	for (i = 0; i < 16; i++) {
		int	combo;
		double	base[8], next_base[8], incr;

		/* Generate the bases for the first cache line */
		if (i == 0) {
			for (j = 0; j < 4; j++) {
				base[2*j+1] = (j * gwdata->FFTLEN / 4 + upper_sse2_word) * gwdata->avg_num_b_per_word;
				base[2*j] = (j * gwdata->FFTLEN / 4) * gwdata->avg_num_b_per_word;
				next_base[2*j+1] = base[2*j+1] + gwdata->avg_num_b_per_word;
				next_base[2*j] = base[2*j] + gwdata->avg_num_b_per_word;
			}
		}

		/* For i=0-7, adjust the weights so that one of bases is an integer */
		/* For i=8-15, adjust the weights so that one of the next_bases is an integer */
		if (i < 8) incr = ceil (base[i]) - base[i];
		else incr = ceil (next_base[i-8]) - next_base[i-8];
		for (j = 0; j < 8; j++) {
			base[j] += incr;
			next_base[j] += incr;
		}
		ASSERTG (ceil (i < 8 ? base[i] : next_base[i-8]) - (i < 8 ? base[i] : next_base[i-8]) < 0.01);

		/* Generate this combination */
		combo = 0;
		for (j = 0; j < 4; j++) {
			combo <<= 2;
			combo += ((int) ceil (next_base[2*j+1]) - (int) ceil (base[2*j+1]) > (int) gwdata->NUM_B_PER_SMALL_WORD) * 2 +
				 ((int) ceil (next_base[2*j]) - (int) ceil (base[2*j]) > (int) gwdata->NUM_B_PER_SMALL_WORD);
		}

		/* Ignore this combo if it is a duplicate.  Otherwise, add it to our combos collection. */
		for (j = 0; ; j++) {
			if (j == num_combos) {
				combos[num_combos++] = combo;
				break;
			}
			if (combo == combos[j]) break;
		}
	}

/* Now concatentate them to save space.  Let's hope they fit in 48 entries. */

	/* Init the next-in-chain array */
	for (i = 0; i < num_combos; i++)
		next_combo[i] = START_OF_CHAIN + END_OF_CHAIN;

	/* Look for 2 chains where the end of one chain has elements in common with */
	/* the start of the other chain,  First look for 3 elements in common, then */
	/* two, then one. */
	for (n = 3; n != 0; n--) {

		/* Examine all chain starts */
		for (i = 0; i < num_combos; i++) {
			int	chain_end;

			/* Skip if not the start of a chain */
			if (! (next_combo[i] & START_OF_CHAIN)) continue;

			/* Find end of chain */
			for (chain_end = i; ! (next_combo[chain_end] & END_OF_CHAIN); chain_end = next_combo[chain_end] & 0xFF);

			/* Now look at all chain ends */
			for (j = 0; j < num_combos; j++) {

				/* Skip if not a chain end */
				if (! (next_combo[j] & END_OF_CHAIN)) continue;
				/* Can't chain to ourselves! */
				if (j == chain_end) continue;

				/* See if chain end has the proper number of common elements with the chain start */
				if (n == 3 && (combos[i] >> 2) == (combos[j] & 0x3F)) {
					next_combo[j] = (next_combo[j] & START_OF_CHAIN) + CHAIN_3_COMMON + i;
					next_combo[i] &= ~START_OF_CHAIN;
					break;
				}
				if (n == 2 && (combos[i] >> 4) == (combos[j] & 0x0F)) {
					next_combo[j] = (next_combo[j] & START_OF_CHAIN) + CHAIN_2_COMMON + i;
					next_combo[i] &= ~START_OF_CHAIN;
					break;
				}
				if (n == 1 && (combos[i] >> 6) == (combos[j] & 0x03)) {
					next_combo[j] = (next_combo[j] & START_OF_CHAIN) + CHAIN_1_COMMON + i;
					next_combo[i] &= ~START_OF_CHAIN;
					break;
				}
			}
		}
	}

/* HACK: Remember the chains in ASM_TIMERS so that we can properly build LIMIT_INVERSE,
/* LIMIT_BIGMAX, and LIMIT_BIGMAX_NEG at a later time. */

	n = 0;
	for (i = 0; i < num_combos; i++) {
		if (! (next_combo[i] & START_OF_CHAIN)) continue;

		/* Output up to 4 bytes for each entry in the chain */
		for (j = i; ; j = next_combo[j] & 0xFF) {
			((char *)gwdata->ASM_TIMERS)[n++] = (combos[j] >> 6);
			if (! (next_combo[j] & CHAIN_3_COMMON))
				((char *)gwdata->ASM_TIMERS)[n++] = (combos[j] >> 4) & 3;
			if (! (next_combo[j] & (CHAIN_3_COMMON + CHAIN_2_COMMON)))
				((char *)gwdata->ASM_TIMERS)[n++] = (combos[j] >> 2) & 3;
			if (! (next_combo[j] & (CHAIN_3_COMMON + CHAIN_2_COMMON + CHAIN_1_COMMON)))
				((char *)gwdata->ASM_TIMERS)[n++] = combos[j] & 3;
			if (next_combo[j] & END_OF_CHAIN) break;
		}
	}
	GWASSERT (n <= 48);		// I've seen n get as high as 35

/* Determine the number of column two-to-phi multipliers */

	num_cols = gwdata->PASS2_SIZE * gwdata->wpn_count;

/* Loop to build table in exactly the same order that it will be */
/* used by the assembly code. */

	p = (unsigned char *) table;
	for (group = 0; group < upper_sse2_word; group += gwdata->PASS1_CACHE_LINES) {
	    for (n = 0; n < gwdata->FFTLEN / 4; n += num_cols) {
	        for (j = 0; j < num_cols; j += gwdata->PASS2_SIZE) {
		    for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
			for (k = 0; k < gwdata->FFTLEN; k += gwdata->FFTLEN / 4) {
			    for (m = 0; m < gwdata->PASS2_SIZE; m += upper_sse2_word) {
				unsigned long word;

/* Now set the big/little flag for a LSW in an SSE2 pair */
/* Otherwise, set the big/little flag for a MSW in an SSE2 pair */

				word = group + j + n + i + k + m;
				if (m == 0) *p = is_big_word (gwdata, word) * 16;
				else if (is_big_word (gwdata, word)) *p += 32;

/* Set the ttp and ttmp fudge flags for two pass FFTs.  The fudge flag is */
/* set if the col mult * the grp mult is b times the correct fft_weight, */
/* meaning a mul by 1/b is required to generate the correct multiplier. */
/* Since we can't do equality compares on floats, this test is a little bit */
/* cryptic. */

				if (gwfft_weight_exponent (gwdata->dd_data, word) + 0.5 <
				    gwfft_weight_exponent (gwdata->dd_data, n + k + m) +
				    gwfft_weight_exponent (gwdata->dd_data, group + j + i)) {
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

/* Apply our method for reducing fudge factor data by 25%.  We've observed that one */
/* of the 4 combinations never occurs. */

			    if (*p >= 128) *p -= 64;

			    p++;
			}

/* Combine last four big/lit 2-bit flag values into one 6-bit flags value. */

			p -= 4;
			for (m = 0; m <= 44; m++) {
				if (((char *)gwdata->ASM_TIMERS)[m]   == (p[0] & 48) >> 4 &&
				    ((char *)gwdata->ASM_TIMERS)[m+1] == (p[1] & 48) >> 4 &&
				    ((char *)gwdata->ASM_TIMERS)[m+2] == (p[2] & 48) >> 4 &&
				    ((char *)gwdata->ASM_TIMERS)[m+3] == (p[3] & 48) >> 4)
					break;
			}
			ASSERTG (m != 45);

/* Combine 1st and 2nd fudge factor flags.  Combine 3rd and 4th fudge factor flags. */
/* Output the 2 combined fudge factors and the 6-bit big/lit flags value using just 2 bytes. */
/* This is a more compact encoding than in previous versions. */

			p[0] = ((p[1] >> 6) + (p[0] >> 6)) * 16 + ((p[3] >> 6) + (p[2] >> 6)) * 2;
			p[1] = (unsigned char) (m << 2);
			p += 2;
		    }
		}
	    }
	}
	return ((double *) p);
}

/* This routine builds a sin/cos table - used by gwsetup */

double *hg_build_sin_cos_table (
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

double *hg_build_premult_table (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	double	*table)		/* Pointer to the table to fill in */
{
	unsigned long i, N, incr, type, m_fudge;

/* Build a premultiplier table for the second pass incrementing by */
/* the pre-calculated pass2_size. */

	N = gwdata->FFTLEN;
	incr = gwdata->PASS2_SIZE;
	if (gwdata->ALL_COMPLEX_FFT) N = N / 2;

/* Mod 2^N+1 arithmetic starts at first data set, */
/* mod 2^N-1 skips some data sets */

 	if (gwdata->ALL_COMPLEX_FFT) i = 0;
	else i = incr * 4;

/* To add in the flipped_m component, we want the sin/cos of flipped_m */
/* over pass2_size.  This fudge factor will convert flipped_m into something */
/* that can be divided by N. */

	m_fudge = N / gwdata->PASS2_SIZE;

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
		if (gwdata->PASS2_SIZE == 1024) grouping_size = 8;
		if (gwdata->PASS2_SIZE == 2048) grouping_size = 8;
		if (gwdata->PASS2_SIZE == 4096) grouping_size = 16;
		if (gwdata->PASS2_SIZE == 8192) grouping_size = 16;
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

double *hg_build_plus1_table (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	double	*table)		/* Pointer to the table to fill in */
{
	unsigned long i, j, k, l, N;
	int	pfa;

/* Set flag if this is a 3*2^n FFT */

	pfa = (gwdata->FFTLEN != pow_two_above_or_equal (gwdata->FFTLEN));

/* Adjust for two-pass FFTs */

	if (gwdata->PASS2_SIZE == 0) N = gwdata->FFTLEN;
	else N = gwdata->FFTLEN / (gwdata->PASS2_SIZE / 2);

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

			if (gwdata->PASS2_SIZE) {
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

double *hg_build_norm_table (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	double	*table,		/* Pointer to the table to fill in */
	int	col)		/* TRUE if building column, not group, table */
{
	unsigned long i, k, num_cols;

/* Handle one-pass FFTs first, there are no group multipliers */

	if (gwdata->PASS2_SIZE == 0) {
		if (!col) return (table);

/* Loop to build table */

		for (i = 0; i < gwdata->FFTLEN; i++) {
			unsigned long j, table_entry;
			double	ttp, ttmp;

/* Call double-precision routine to compute the two multipliers */

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

	num_cols = gwdata->PASS2_SIZE / 2;
	if (col) {

/* Loop to build table */

		for (i = 0; i < num_cols; i++) {
			double	ttp, ttmp;

/* Call double-precision routine to compute the two multipliers */

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

/* Call double-precision routine to compute the two multipliers */

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

double *hg_build_biglit_table (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	double	*table)		/* Pointer to the table to fill in */
{
	unsigned char *p;
	unsigned long h, i, j, k, m, u, gap;
	unsigned long pfa, hlimit, haddin, mmult, umult;

/* Handle one pass FFTs differently */

	if (gwdata->PASS2_SIZE == 0) {

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

	gap = gwdata->PASS2_SIZE / 2;

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
		    gwfft_weight_exponent (gwdata->dd_data, word % gap) +
		    gwfft_weight_exponent (gwdata->dd_data, word - word % gap)) {
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

/* Ancient x87 setup routines */

#ifndef X86_64

/* This routine builds an x87 sin/cos table - used by gwsetup */

double *x87_build_sin_cos_table (
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

double *x87_build_premult_table (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	double	*table)		/* Pointer to the table to fill in */
{
	unsigned long i, N, incr;

/* Build a premultiplier table for the second pass incrementing by */
/* the pre-calculated pass2_size. */

	N = gwdata->FFTLEN;
	incr = gwdata->PASS2_SIZE;
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

double *x87_build_plus1_table (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	double	*table)		/* Pointer to the table to fill in */
{
	unsigned long i, k, N;
	int	pfa;

/* Set flag if this is a 3*2^n FFT */

	pfa = (gwdata->FFTLEN != pow_two_above_or_equal (gwdata->FFTLEN));

/* Adjust for two-pass FFTs */

	if (gwdata->PASS2_SIZE == 0) N = gwdata->FFTLEN;
	else N = gwdata->FFTLEN / gwdata->PASS2_SIZE;

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

double *x87_build_norm_table (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	double	*table,		/* Pointer to the table to fill in */
	int	col)		/* TRUE if building column, not group, table */
{
	unsigned long i, k, num_cols;

/* Handle one-pass FFTs first, there are no group multipliers */

	if (gwdata->PASS2_SIZE == 0) {
		if (!col) return (table);

/* Loop to build table */

		for (i = 0; i < gwdata->FFTLEN; i++) {
			unsigned long j;
			double	ttp, ttmp;

/* Call double-precision routine to compute the two multipliers */

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

	num_cols = gwdata->PASS2_SIZE;
	if (col) {

/* Loop to build columns table */

		for (i = 0; i < num_cols; i++) {
			double	ttp, ttmp;

/* Call double-precision routine to compute the two multipliers */

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

/* Call double-precision routine to compute the two multipliers */

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

double *x87_build_biglit_table (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	double	*table)		/* Pointer to the table to fill in */
{
	unsigned char *p;
	unsigned long i, j, k, m;

/* Handle one pass FFTs differently */

	if (gwdata->PASS2_SIZE == 0) {

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
	for (i = 0; i < gwdata->PASS2_SIZE; i += gwdata->PASS1_CACHE_LINES * 2) {
	for (j = 0; j < gwdata->FFTLEN / 2; j += gwdata->PASS2_SIZE) {
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
		    gwfft_weight_exponent (gwdata->dd_data, word % gwdata->PASS2_SIZE) +
		    gwfft_weight_exponent (gwdata->dd_data, word - word % gwdata->PASS2_SIZE)) {
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

/* End ancient x87 setup routines */

#endif

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

/* Determine the "bif" value we will look for.  Often this is a straight-forward mapping from */
/* the CPU_ARCHITECTURE.  However, for some CPU architectures, like Pentium M and Core Solo, we */
/* don't have jmptable entries detailing the fastest FFT implementations for those architectures. */

int calculate_bif (
	gwhandle *gwdata,	/* Gwnum global data */
	unsigned long fftlen)
{
	int	retval;

/* Map the CPU architecture as determined by CPUID to one of the CPU architectures */
/* that the FFT assembly code is optimized for. */

	switch (CPU_ARCHITECTURE) {
	case CPU_ARCHITECTURE_PENTIUM_M:	/* Not sure what is best for these three architectures */
	case CPU_ARCHITECTURE_CORE:
	case CPU_ARCHITECTURE_ATOM:
		retval = 2;			/* Look for FFTs optimized for large cache P4s */
		break;
	case CPU_ARCHITECTURE_PENTIUM_4:
		if (CPU_L2_CACHE_SIZE <= 128)	/* We haven't optimized for these yet */
			retval = 4;		/* Look for FFTs optimized for 256K cache P4s */
		else if (CPU_L2_CACHE_SIZE <= 256)
			retval = 4;		/* Look for FFTs optimized for 256K cache P4s */
		else if (CPU_L2_CACHE_SIZE <= 512)
			retval = 3;		/* Look for FFTs optimized for 512K cache P4s */
		else if (gwdata->cpu_flags & CPU_TLB_PRIMING)
			retval = 3;		/* Look for FFTs optimized for 512K cache P4s */
		else
			retval = 2;		/* Look for FFTs optimized for large cache P4s */
		break;
	case CPU_ARCHITECTURE_CORE_2:
		retval = 0;			/* Look for FFTs optimized for Core 2 */
		break;
	case CPU_ARCHITECTURE_CORE_I7:
		retval = 1;			/* Look for FFTs optimized for Core i3/i5/i7/i9 */
		break;
	case CPU_ARCHITECTURE_INTEL_OTHER:	/* This is probably one of Intel's next generation CPUs */ 
		retval = 1;			/* Look for FFTs optimized for Core i3/i5/i7/i9 */
		break;
	case CPU_ARCHITECTURE_AMD_K8:
		retval = 6;			/* Look for FFTs optimized for K8 */
		break;
	case CPU_ARCHITECTURE_AMD_K10:
		retval = 7;			/* Look for FFTs optimized for K10 */
		break;
	case CPU_ARCHITECTURE_PRE_SSE2:		/* Cannot happen, gwinfo should have selected x87 FFTs */
	case CPU_ARCHITECTURE_AMD_OTHER:
	case CPU_ARCHITECTURE_OTHER:
	default:
		retval = 0;			/* For no particularly good reason, look for FFTs optimized for Core 2 */
		break;
	}

/* For slower CPU architectures we didn't bother to find the best FFT implementation */
/* for the larger FFTs.  This was done to reduce the size of the executable.  If we */
/* are asked to run one of these large FFTs, select an FFT optimized for a different */
/* CPU architecture. */

	if (fftlen > 4194304 && (retval == 3 || retval == 4))
		retval = 2;		/* Small cache P4s have best FFT implementations up to 4M */
	if (fftlen > 6291456 && retval == 2)
		retval = 0;		/* P4s have best FFT implementations up to 6M */
	if (fftlen > 6291456 && retval == 6)
		retval = 7;		/* K8s have best FFT implementations up to 6M */

/* Return the result */

	return (retval);
}

/* Ugly little macros to bump jmptable pointer to next procedure entry or to next count */
#define INC_JMPTAB_1(x)	x = (struct gwasm_jmptab *) ((char *)(x) + sizeof(uint32_t) + sizeof(void*) + sizeof(uint32_t))
#define DEC_JMPTAB_1(x)	x = (struct gwasm_jmptab *) ((char *)(x) - sizeof(uint32_t) - sizeof(void*) - sizeof(uint32_t))
#define INC_JMPTAB_2(x)	x = (struct gwasm_jmptab *) ((char *)(x) + sizeof(int32_t))

/* This routine checks to see if there is an FFT implementation for this FFT length and */
/* CPU architecture.  For example, when the FFT length is just less than a power of two, on */
/* some CPUs it may be better to use the larger power-of-two FFT length and thus there */
/* will not be an FFT implementation for this slightly smaller FFT length. */ 

int is_fft_implemented (
	gwhandle *gwdata,	/* Gwnum global data */
	struct gwasm_jmptab *jmptab)
{
	int	desired_bif;		/* The "best implementation for" value we will look for. */
					/* See mult.asm for defined BIF_ values. */

/* For FFT lengths below 7K there is only one FFT implementation and it is always available */

	if (jmptab->fftlen < 7168) return (TRUE);

/* If we are benchmarking all FFT implementations or we are doing QA, then we want to test */
/* this FFT length even if it isn't optimal */

	if (gwdata->bench_pick_nth_fft || gwdata->qa_pick_nth_fft) return (TRUE);

/* Determine the "bif" value we will look for.  Often this is a straight-forward mapping from */
/* the CPU_ARCHITECTURE.  However, for some CPU architectures, like Pentium M and Core Solo, we */
/* don't have jmptable entries detailing the fastest FFT implementations for those architectures. */

	desired_bif = calculate_bif (gwdata, jmptab->fftlen);

/* Loop through the FFT implementations to see if we find an implementation */
/* that matches our desired "bif" value. */

	while (jmptab->flags & 0x80000000) {
		if (((jmptab->flags >> 13) & 0xF) == desired_bif) return (TRUE);
		INC_JMPTAB_1 (jmptab);
	}

/* FFT implementation not found.  A larger FFT length should be faster. */

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
	unsigned long max_exp;
	char	buf[20];
	int	qa_nth_fft, desired_bif;
	void	*prev_proc_ptrs[5];
	uint32_t flags, no_prefetch, in_place;
	struct gwasm_data *asm_data;

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
	if (gwdata->specific_fftlen == 0 &&
	    (k > 1.0 || n < 500 || abs (c) > 1) &&
	    gwdata->qa_pick_nth_fft < 1000) {

/* Use the proper 2^N-1 jmptable */

		if (gwdata->cpu_flags & CPU_SSE2)
			zpad_jmptab = asm_info.p4_cyclic_fft_info;
		else zpad_jmptab =
			zpad_jmptab = asm_info.x86_cyclic_fft_info;

/* Find the table entry for the FFT that can do a mod 2^2n FFT, handling */
/* k and c in the normalization routines.  We will compare this to the */
/* non-zero-padded FFT length later.  The zeroes in the upper half of FFT */
/* input data let us get about another 0.3 bits per input word. */

		while ((max_exp = zpad_jmptab->max_exp) != 0) {

/* Do a quick check on the suitability of this FFT */

			if ((double) n * log2b / (double) zpad_jmptab->fftlen > 26.0) goto next1;
			if (zpad_jmptab->fftlen < gwdata->minimum_fftlen) goto next1;
			if (! is_fft_implemented (gwdata, zpad_jmptab)) goto next1;

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

/* We can use this FFT. */
/* Look for a non-zero-padded FFT that might be even faster. */

				break;
			}

/* Move past procedure entries and counts to next jmptable entry */

next1:			while (zpad_jmptab->flags & 0x80000000) INC_JMPTAB_1 (zpad_jmptab);
			DEC_JMPTAB_1 (zpad_jmptab);
			while (zpad_jmptab->counts[0]) INC_JMPTAB_2 (zpad_jmptab);
			zpad_jmptab = (struct gwasm_jmptab *) &zpad_jmptab->counts[1];
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

	while ((max_exp = jmptab->max_exp) != 0) {

/* Check if this table entry matches the specified FFT length. */

		if (gwdata->specific_fftlen) {
			if (gwdata->specific_fftlen == jmptab->fftlen) break;
			goto next2;
		}

/* Do a quick check on the suitability of this FFT */

		if ((double) n * log2b / (double) jmptab->fftlen > 26.0) goto next2;
		if (jmptab->fftlen < gwdata->minimum_fftlen) goto next2;
		if (! is_fft_implemented (gwdata, jmptab)) goto next2;

/* Check if this FFT length will work with this k,n,c combo */

//  This is the old code which only supported b == 2
//		double max_bits_per_word;
//		double bits_per_word;
//
/* Compute the maximum number of bits allowed in the FFT input word */
//
//		max_bits_per_word = (double) max_exp / jmptab->fftlen;
//		max_bits_per_word -= gwdata->safety_margin;
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
//		bits_per_word = (log2k + n * log2b) / jmptab->fftlen;
//		if (2.0 * bits_per_word + log2k + 1.7 * log2c <=
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
//		loglen = log2 (jmptab->fftlen);
//		total_bits = (bits_per_word - 1.0) * 2.0 +
//				     1.7 * log2c + loglen * 0.6 +
//				     log2k + log2maxmulbyconst;
//		if (total_bits > 5.0 * bits_per_word) {

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
/* Also, some examples such as 19464*19^31895+1 still raise round off errors. */
/* For added safety we assume an extra 0.25 bits of output are needed when */
/* base is not 2. */

		else if (! is_pathological_distribution (num_big_words, num_small_words)) {
			weighted_bits_per_output_word -=
					((log2b <= 3.0) ? (log2b - 1.0) / 2.0 :
					 (log2b <= 6.0) ? 1.0 + (log2b - 3.0) / 3.0 :
							  2.0 + (log2b - 6.0) / 6.0);
			if (b != 2) weighted_bits_per_output_word += 0.25;
		}

/* If the bits in an output word is less than the maximum allowed, we can */
/* probably use this FFT length -- though we need to do a few more tests. */

		if (weighted_bits_per_output_word <= max_weighted_bits_per_output_word) {

/* Unfortunately, users uncovered some cases that the spreading carry over 4 words */
/* test in the next paragraph did not push us to a higher FFT length.  The problem */
/* is that if the carry doesn't fit, the next FFT multiply will have some words that */
/* are 1 bit larger than expected.  This leads to an even greater chance that the next */
/* carry won't spread.  This can spiral out of control.  Since this situation is much */
/* more serious than the reason above where we choose a higher FFT, we will be even more */
/* conservative in our estimate of bits_per_output_word.  We suspect our estimate is */
/* more likely to be wrong when b != 2, simply because we have a much longer history */
/* of dealing with the b = 2 case.  Someday we should consider implementing spreading */
/* the carry over 5 or 6 words.  The test cases that necessitated this fix were */
/* 9238*619^619+1, 9276*626^626+1, 9876*626^626+1, 7504*627^627+1, 8004*627^627+1, */
/* 9056*627^627+1, 9256*627^627+1, 9386*635^635+1, 8619*650^650+1, 9732*650^650+1. */
/* No one has yet reported a failure for the b = 2 case but we'll add some safety anyway. */

			if (b == 2) bits_per_output_word += 0.3;
			else bits_per_output_word += 0.9;
			
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
//				ASSERTG (zpad_jmptab == NULL ||
//					 jmptab->fftlen >= zpad_jmptab->fftlen);
				goto next2;
			}

/* Because of limitations in the top_carry_adjust code, there is a limit */
/* on the size of k that can be handled.  This isn't a big deal since the */
/* zero-padded implementation should use the same FFT length.  Check to see */
/* if this k can be handled.  K must fit in the top three words for */
/* one-pass FFTs and within the top two words of two-pass FFTs. */

			if ((jmptab->flags & 0x3F) == 0 && logbk > floor (3.0 * b_per_input_word)) {
// This assert is designed to find any cases where using 3 or more carray adjust words
// would use a shorter FFT than using a zero-padded FFT.
				ASSERTG (zpad_jmptab != NULL && jmptab->fftlen >= zpad_jmptab->fftlen);
				goto next2;
			}
			if ((jmptab->flags & 0x3F) != 0 && logbk > floor (2.0 * b_per_input_word)) {
// This assert is designed to find any cases where using 3 or more carray adjust words
// would use a shorter FFT than using a zero-padded FFT.
				ASSERTG (zpad_jmptab != NULL && jmptab->fftlen >= zpad_jmptab->fftlen);
				goto next2;
			}

/* We've found an FFT length to use */

			break;
		}

/* Move past procedure entries and counts to next jmptable entry */

next2:		while (jmptab->flags & 0x80000000) INC_JMPTAB_1 (jmptab);
		DEC_JMPTAB_1 (jmptab);
		while (jmptab->counts[0]) INC_JMPTAB_2 (jmptab);
		jmptab = (struct gwasm_jmptab *) &jmptab->counts[1];
	}

/* If the zero pad FFT length is less than the DWT FFT length OR we */
/* are QA'ing every FFT implementation, then use the zero pad FFT length. */

	if (zpad_jmptab != NULL && zpad_jmptab->max_exp &&
	    (jmptab->max_exp == 0 || zpad_jmptab->fftlen < jmptab->fftlen || gwdata->qa_pick_nth_fft)) {
		gwdata->ZERO_PADDED_FFT = TRUE;
		gwdata->ALL_COMPLEX_FFT = FALSE;
		jmptab = zpad_jmptab;
	}

/* If we found a DWT table entry then use it. */

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

/* We've found the right "jump" table entry, save the pointer and FFT length */

	gwdata->jmptab = jmptab;
	gwdata->FFTLEN = jmptab->fftlen;

/************************************************************************/
/* Decide which implementation of this FFT length is best for this CPU. */
/************************************************************************/

/* Loop through all the implementations for this FFT length until we find */
/* the one best suited to this CPU. */

	qa_nth_fft = gwdata->ZERO_PADDED_FFT ? 100 : 1000;
	desired_bif = calculate_bif (gwdata, gwdata->FFTLEN);
	prev_proc_ptrs[0] = NULL;
	prev_proc_ptrs[1] = NULL;
	prev_proc_ptrs[2] = NULL;
	prev_proc_ptrs[3] = NULL;
	prev_proc_ptrs[4] = NULL;
	for ( ; ; ) {
		int	arch;		/* (ppro=0,p3=1,p4=2,core=3,k8=4,k10=5) */
		int	best_impl_for;	/* (CORE2=0,I7=1,etc.  See BIF_ definitions in mult.asm) */
		int	fft_type;	/* (home-grown=0, radix-4=1, r4delay=2, r4dwpn=3) */

/* Handle an FFT implementation not found condition.  Should only happen */
/* if we're benchmarking or QA'ing and we've tested every implementation. */

		if (! (jmptab->flags & 0x80000000)) {

/* If we are QA'ing every FFT implementation and we did not find another */
/* zero-padded FFT implementation, then go find a non-zero-padded one. */

			if (gwdata->qa_pick_nth_fft && gwdata->qa_pick_nth_fft < 1000) {
				gwdata->qa_pick_nth_fft = 1000;
				goto again;
			}

/* Else return an error */

			return (GWERROR_TOO_LARGE);
		}

/* If this CPU will crash running this FFT then skip this entry. */
/* Our K8 and K10 optimized FFTs requires prefetchw (3DNow!) capability. */
/* If this is an Intel CPU, skip over these implementations. */

		arch = (jmptab->flags >> 17) & 0xF;
		if ((arch == 4 || arch == 5) && ! (gwdata->cpu_flags & CPU_3DNOW))
			goto next3;

/* Handle benchmarking case that selects the nth FFT implementation */
/* without regard to any other consideration.  NOTE: Due to the extreme */
/* penalty a K8 pays for using the movaps instruction that the Core and P4 */
/* implementations use, we will not benchmark these on a K8. */

		if (gwdata->bench_pick_nth_fft) {
			if (jmptab->proc_ptr == prev_proc_ptrs[0]) goto next3;
			if (jmptab->proc_ptr == prev_proc_ptrs[1]) goto next3;
			if (jmptab->proc_ptr == prev_proc_ptrs[2]) goto next3;
			if (jmptab->proc_ptr == prev_proc_ptrs[3]) goto next3;
			if (jmptab->proc_ptr == prev_proc_ptrs[4]) goto next3;
			if (CPU_ARCHITECTURE == CPU_ARCHITECTURE_AMD_K8 &&
			    (jmptab->flags & 0x3F) != 0 &&
			    (arch == 2 || arch == 3))
				goto next3;
			gwdata->bench_pick_nth_fft--;
			if (gwdata->bench_pick_nth_fft) goto next3;
			break;
		}

/* Handle the QA case that tests every possible FFT implementation */
/* Remember the FFT returned so that we can return a different FFT to */
/* the QA code next time. */

		if (gwdata->qa_pick_nth_fft) {
			if (jmptab->proc_ptr == prev_proc_ptrs[0]) goto next3;
			if (jmptab->proc_ptr == prev_proc_ptrs[1]) goto next3;
			if (jmptab->proc_ptr == prev_proc_ptrs[2]) goto next3;
			if (jmptab->proc_ptr == prev_proc_ptrs[3]) goto next3;
			if (jmptab->proc_ptr == prev_proc_ptrs[4]) goto next3;
			if (CPU_ARCHITECTURE == CPU_ARCHITECTURE_AMD_K8 &&
			    (arch == 2 || arch == 3))
				goto next3;
			qa_nth_fft++;
			if (qa_nth_fft <= gwdata->qa_pick_nth_fft) goto next3;
			gwdata->qa_picked_nth_fft = qa_nth_fft;
			break;
		}

/* If this is a small FFT or an x87 FFT, then there is only one implementation */

		if (gwdata->FFTLEN < 7168) break;
		if (! (gwdata->cpu_flags & CPU_SSE2)) break;

/* The Radix-4/8 DJB FFT with partial normalization saves a few multiplies by doing part */
/* of the normalization during the forward and inverse FFT.  Unfortunately, this optimization */
/* makes the SUM(INPUTS) != SUM(OUTPUTS) error check impossible.  If the user prefers */
/* SUM(INPUTS) != SUM(OUTPUTS) error checking, then skip r4dwpn FFT type (except all-complex */
/* FFTs which never could support SUM(INPUTS) != SUM(OUTPUTS) error checking. */

		fft_type = (jmptab->flags >> 21) & 0xF;
		if (gwdata->sum_inputs_checking &&
		    ! gwdata->ALL_COMPLEX_FFT &&
		    fft_type == FFT_TYPE_RADIX_4_DWPN)
			goto next3;

/* See if this is the best implementation for this CPU architecture */

		best_impl_for = (jmptab->flags >> 13) & 0xF;
		if (best_impl_for == desired_bif) break;

/* Move onto the next FFT implementation */

next3:		prev_proc_ptrs[4] = prev_proc_ptrs[3];
		prev_proc_ptrs[3] = prev_proc_ptrs[2];
		prev_proc_ptrs[2] = prev_proc_ptrs[1];
		prev_proc_ptrs[1] = prev_proc_ptrs[0];
		prev_proc_ptrs[0] = jmptab->proc_ptr;
		INC_JMPTAB_1 (jmptab);
	}

/* Remember the information from the chosen FFT implementation */

	flags = jmptab->flags;
	gwdata->GWPROCPTRS[0] = jmptab->proc_ptr;
	gwdata->mem_needed = jmptab->mem_needed;

/* Break the flags word from the jmptable entry into its constituent parts */
/* The 32-bit flags word is as follows (copied from mult.asm): */
/*	80000000h		always on */
/*	2 SHL 26		(no prefetching - not used by gwnum) */
/*	1 SHL 26		(in_place) */
/*	fft_type SHL 21		(hg=0, r4=1, r4delay=2) */
/*	arch SHL 17		(ppro=0,p3=1,p4=2,core=3,k8=4,k10=5) */
/*	best_impl_for SHL 13	(CORE2=0,I7=1,P4_1024=2,etc.) */
/*	clm SHL 9		(1,2,4,8) */
/*	pass2size_over_64	(many valid values) */

	no_prefetch = (flags >> 27) & 0x00000001;
	in_place = (flags >> 26) & 0x00000001;
	gwdata->FFT_TYPE = (flags >> 21) & 0x0000000F;
	gwdata->ARCH = (flags >> 17) & 0x0000000F;
	gwdata->PASS1_CACHE_LINES = ((flags >> 9) & 0x0000000F) * 2;
	gwdata->PASS2_SIZE = (flags & 0x0000001FF) << 6;

/* Set more info so that addr_offset called from */
/* gwmap_to_estimated_size can work without a call to gwsetup. */

	gwdata->GW_ALIGNMENT = 4096;	/* Guess an alignment so gwsize can */
					/* return a reasonable value for */
					/* large page support in gwsetup */

/* Calculate the scratch area size -- needed by gwmemused without calling gwsetup */

	if (gwdata->PASS2_SIZE && !in_place) {
		if (! (gwdata->cpu_flags & CPU_SSE2)) {	// x87 scratch area size
			int	pass1_chunks, gaps;
			// x87 scratch area size, pads 64 bytes every 32 chunks
			pass1_chunks = (gwdata->FFTLEN / gwdata->PASS2_SIZE) >> 1;
			gaps = pass1_chunks / 32 - 1;
			if (gaps < 3) gaps = 0;
			// x87 pads 64 bytes every 32 chunks
			gwdata->SCRATCH_SIZE = pass1_chunks * gwdata->PASS1_CACHE_LINES * 32 + gaps * 64;
		} else {				// SSE2 scratch area size
			int	pass1_chunks, gaps;
			// SSE2 pads 128 bytes every 8 chunks
			pass1_chunks = (gwdata->FFTLEN / gwdata->PASS2_SIZE) >> 2;
			gaps = pass1_chunks / 8 - 1;
			gwdata->SCRATCH_SIZE = pass1_chunks * gwdata->PASS1_CACHE_LINES * 64 + gaps * 128;
		}
	} else
		gwdata->SCRATCH_SIZE = 0;

/* Calculate the gap between data blocks.  This is used by addr_offset called from */
/* gwmap_to_estimated_size without a call to gwsetup.  This must match the calculation */
/* done in the set_FFT_constants macro in xmult.mac */
	
	if (gwdata->PASS2_SIZE * 2 * 2 * 8 / 8192 * 128 % 256 == 0)
		gwdata->PASS2GAPSIZE = -128;
	else
		gwdata->PASS2GAPSIZE = 0;

/* Copy the counts */

	asm_data = (struct gwasm_data *) gwdata->asm_data;
	if (asm_data != NULL) {
		if (! (gwdata->cpu_flags & CPU_SSE2)) {
			if (gwdata->PASS2_SIZE == 0) {
				/* 5 counts for one pass x87 FFTs */
				asm_data->count1 = jmptab->counts[0];
				asm_data->count2 = jmptab->counts[1];
				asm_data->count3 = jmptab->counts[2];
				asm_data->count4 = jmptab->counts[3];
				asm_data->count5 = jmptab->counts[4];
			} else {
				/* Count of pass 2 sections */
				asm_data->addcount1 = gwdata->FFTLEN / gwdata->PASS2_SIZE / 2;
				/* Count of all-complex pass 2 sections */
				asm_data->count1 = asm_data->addcount1;
				if (gwdata->ZERO_PADDED_FFT || c < 0) asm_data->count1--;
			}
		} else {
			if (gwdata->PASS2_SIZE == 0) {
				/* 7 counts for one pass SSE2 FFTs */
				asm_data->addcount1 = jmptab->counts[0];
				asm_data->normcount1 = jmptab->counts[1];
				asm_data->count1 = jmptab->counts[2];
				asm_data->count2 = jmptab->counts[3];
				asm_data->count3 = jmptab->counts[4];
				asm_data->count4 = jmptab->counts[5];
				asm_data->count5 = jmptab->counts[6];
			} else if (gwdata->FFT_TYPE == FFT_TYPE_HOME_GROWN) {
				int	pfa, pfa_shift;
				/* Count of pass 2 sections */
				asm_data->addcount1 = gwdata->FFTLEN / gwdata->PASS2_SIZE / 4;
				/* Count of all-complex pass 2 sections */
				asm_data->count1 = asm_data->addcount1;
				if (gwdata->ZERO_PADDED_FFT || c < 0) asm_data->count1--;
				/* Up to three section counts for gw_carries.  Examples: */
				/* 320 pass 2 sections: (256 << 11) + 64 */
				/* 384 pass 2 sections: 384 */
				/* 448 pass 2 sections: (256 << 22) + (128 << 11) + 64 */
				/* 512 pass 2 sections: 512 */
				for (pfa = asm_data->addcount1, pfa_shift = 0;
				     pfa > 8;
				     pfa >>= 1, pfa_shift++);
				if (pfa == 5)
					asm_data->count3 = ((4 << 11) + 1) << pfa_shift;
				else if (pfa == 7)
					asm_data->count3 = ((4 << 22) + (2 << 11) + 1) << pfa_shift;
				else
					asm_data->count3 = asm_data->addcount1;
			} else {
				/* Count of pass 2 sections */
				asm_data->addcount1 = gwdata->FFTLEN / gwdata->PASS2_SIZE / 4;
				/* Count of all-complex pass 2 sections */
				asm_data->count1 = asm_data->addcount1;
				if (gwdata->ZERO_PADDED_FFT || c < 0) asm_data->count1--;
				/* Only one section counts for gw_carries */
				asm_data->count3 = asm_data->addcount1;
			}
		}
	}

/* All done */

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

/* Sanity check the k,b,n values */

	if (k < 1.0) return (GWERROR_K_TOO_SMALL);
	if (k > 9007199254740991.0) return (GWERROR_K_TOO_LARGE);
	if (log2(b) * (double) n > MAX_PRIME_SSE2) return (GWERROR_TOO_LARGE);
	if (! (gwdata->cpu_flags & CPU_SSE2) && log2(b) * (double) n > MAX_PRIME) return (GWERROR_TOO_LARGE);

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
		/* If we've already copied the modulus, use it.  For example, */
		/* gwsetup_general_mod_giant on (2^313+1)/3 will call this routine */
		/* to try an IBDWT on 2^313+1.  This number is too small and */
		/* we need to revert back to a general mod on (2^313+1)/3. */
		if (gwdata->GW_MODULUS != NULL) {
			gwdata->force_general_mod = TRUE;
			error_code = gwsetup_general_mod_giant (gwdata, gwdata->GW_MODULUS);
			if (error_code) return (error_code);
		} else {
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
			/* If GCD(k,c) was not 1, then the reciprocal in generic modular */
			/* reduction may well have a nasty bit pattern.  The test case that */
			/* brought this to light is 3*2^77574+3.  We don't have a lot of */
			/* experience with these kinds of composite numbers.  For now, we */
			/* just increase the MAXDIFF setting.  We may also need to increment */
			/* safety margin to select a larger FFT size in some cases. */
			if (gcd > 1) gwdata->MAXDIFF *= 100.0;
		}
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

/* Examine the giant to see if it a (k*2^n+c)/d value that we can better optimize. */
/* Also detect nasty bit patterns, like Phi (82730,2), where multiplying by a small d */
/* results in a less nasty bit pattern for the modulus. */

	d = 1;
	convertible = (!gwdata->force_general_mod &&
		       bits > 300 &&
		       convert_giant_to_k2ncd (g, &k, &n, &c, &d));

/* Copy the modulus except for convertible k*2^n+c values */

	if (!convertible || d != 1) {
		/* Yuk, we may already have saved the modulus.  For example, (2^313+1)/3 */
		/* will come through here and save the modulus.  But gwsetup of 2^313+1 */
		/* is too small for IBDWT, so this routine is recursively called. */
		/* We cannot reallocate because that will cause a memory leak. */
		if (gwdata->GW_MODULUS == NULL) {
			gwdata->GW_MODULUS = allocgiant ((bits >> 5) + 1);
			if (gwdata->GW_MODULUS == NULL) {
				gwdone (gwdata);
				return (GWERROR_MALLOC);
			}
			gtog (g, gwdata->GW_MODULUS);
		}
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
/* multiply the modulus by 11.  We'll assume that d produces a better */
/* modulus candidate if more than a quarter of the words are zero or */
/* minus one -- at least until someone improves on this scheme. */

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
	unsigned long pass1_size;
	double	small_word, big_word, temp, asm_values[40];

/* Remember the arguments */

	gwdata->k = k;
	gwdata->b = b;
	gwdata->n = n;
	gwdata->c = c;

/* Init the FPU to assure we are in 64-bit precision mode */

	fpu_init ();

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

/* Select the proper FFT size for this k,b,n,c combination */

	error_code = gwinfo (gwdata, k, b, n, c);
	if (error_code) {
		gwdone (gwdata);
		return (error_code);
	}
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
	mem_needed = gwdata->mem_needed + gwdata->SCRATCH_SIZE;
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
/* weighting factors are 1.0.  This happens when abs(c) is 1 and every */
/* FFT word has the same number of b's.  The assembly code can make some */
/* obvious optimizations when all the FFT weights are one. */

	gwdata->RATIONAL_FFT = asm_data->RATIONAL_FFT =
		((double) gwdata->NUM_B_PER_SMALL_WORD == gwdata->avg_num_b_per_word) && (abs (c) == 1);

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

/* Determine the pass 1 size.  This affects how we build */
/* many of the sin/cos tables. */

	if (gwdata->PASS2_SIZE)
		pass1_size = gwdata->FFTLEN / gwdata->PASS2_SIZE; /* Real values in a pass1 section */
	else
		pass1_size = gwdata->FFTLEN;

/* Initialize tables for the radix-4 FFT assembly code.  These are guaranteed to be */
/* two-pass FFTs as we use the older FFT code for one pass FFTs. */

	if (gwdata->FFT_TYPE == FFT_TYPE_RADIX_4) {

/* Build sin/cos and premultiplier tables used in pass 1 of two pass FFTs. */
/* For best prefetching, make sure tables remain on 128-byte boundaries */

		ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
		tables = r4_build_pass1_table (gwdata, tables);

/* Build the sin/cos table used in complex pass 2 blocks */

		ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
		asm_data->xsincos_complex = tables;
		tables = r4_build_pass2_complex_table (gwdata, tables);
		ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
		asm_data->sincos3 = tables;
		tables = r4_build_pass2_real_table (gwdata, tables);

/* Allocate a table for carries.  Init with XMM_BIGVAL.  For best */
/* distribution of data in the L2 cache, make this table contiguous */
/* with all other data used in the first pass (scratch area, normalization */
/* tables, etc.)  Note that we put the tables that are only partly loaded */
/* (column multipliers and big/lit table after the tables that are */
/* loaded throughout the first pass. */

		if (gwdata->PASS2_SIZE) {
			int	i, carry_table_size;
			double	xmm_bigval;
			ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
			asm_data->carries = tables;
			carry_table_size = gwdata->FFTLEN / (gwdata->PASS2_SIZE / 2);
			xmm_bigval = 3.0 * 131072.0 * 131072.0 * 131072.0;
			for (i = 0; i < carry_table_size; i++) *tables++ = xmm_bigval;
			tables += (16 - (tables - gwdata->gwnum_memory)) & 15;
		}

/* Build the group muliplier normalization table.  Keep this table */
/* contiguous with other data used in pass 1. */

		ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
		asm_data->norm_grp_mults = tables;
		tables = r4_build_norm_table (gwdata, tables, 0);

/* Reserve room for the pass 1 scratch area. */

		if (gwdata->SCRATCH_SIZE) {
			ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
			asm_data->scratch_area = tables;
			tables = (double *) ((char *) tables + gwdata->SCRATCH_SIZE);
		}

/* Build the table of big vs. little flags.  This cannot be last table */
/* built as xnorm_2d macro reads 2 bytes past the end of the array. */

		ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
		asm_data->norm_biglit_array = tables;
		tables = r4_build_biglit_table (gwdata, tables);
		tables += (16 - (tables - gwdata->gwnum_memory)) & 15;

/* Build the column normalization multiplier table. */

		ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
		asm_data->norm_col_mults = tables;
		tables = r4_build_norm_table (gwdata, tables, 1);
	}

/* Initialize tables for a modified radix-4 FFT.  This particular FFT */
/* uses a radix-8 building block when there are an odd number of FFT */
/* levels in a pass. */

	else if (gwdata->FFT_TYPE == FFT_TYPE_RADIX_4_DELAYED) {

/* Build sin/cos and premultiplier tables used in pass 1 of two pass FFTs. */
/* For best prefetching, make sure tables remain on 128-byte boundaries */

		ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
		tables = r4delay_build_pass1_table (gwdata, tables);
		ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
		asm_data->sincos1 = tables;
		tables = r4delay_build_fixed_premult_table (gwdata, tables);
		ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
		asm_data->sincos2 = tables;
		tables = r4delay_build_fixed_pass1_table (gwdata, tables);

/* Build the sin/cos table used in complex pass 2 blocks */
/* The pass 2 tables are the same as for a traditional radix-4 FFT */		

		ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
		asm_data->xsincos_complex = tables;
		tables = r4_build_pass2_complex_table (gwdata, tables);
		ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
		asm_data->sincos3 = tables;
		tables = r4_build_pass2_real_table (gwdata, tables);

/* Allocate a table for carries.  Init with XMM_BIGVAL.  For best */
/* distribution of data in the L2 cache, make this table contiguous */
/* with all other data used in the first pass (scratch area, normalization */
/* tables, etc.)  Note that we put the tables that are only partly loaded */
/* (column multipliers and big/lit table after the tables that are */
/* loaded throughout the first pass. */

		if (gwdata->PASS2_SIZE) {
			int	i, carry_table_size;
			double	xmm_bigval;
			ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
			asm_data->carries = tables;
			carry_table_size = gwdata->FFTLEN / (gwdata->PASS2_SIZE / 2);
			xmm_bigval = 3.0 * 131072.0 * 131072.0 * 131072.0;
			for (i = 0; i < carry_table_size; i++) *tables++ = xmm_bigval;
			tables += (16 - (tables - gwdata->gwnum_memory)) & 15;
		}

/* Build the group muliplier normalization table.  Keep this table */
/* contiguous with other data used in pass 1. */

		ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
		asm_data->norm_grp_mults = tables;
		tables = r4_build_norm_table (gwdata, tables, 0);

/* Reserve room for the pass 1 scratch area. */

		if (gwdata->SCRATCH_SIZE) {
			ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
			asm_data->scratch_area = tables;
			tables = (double *) ((char *) tables + gwdata->SCRATCH_SIZE);
		}

/* Build the table of big vs. little flags.  This cannot be last table */
/* built as xnorm_2d macro reads 2 bytes past the end of the array. */

		ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
		asm_data->norm_biglit_array = tables;
		tables = r4_build_biglit_table (gwdata, tables);
		tables += (16 - (tables - gwdata->gwnum_memory)) & 15;

/* Build the column normalization multiplier table. */

		ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
		asm_data->norm_col_mults = tables;
		tables = r4_build_norm_table (gwdata, tables, 1);
	}

/* Initialize tables for an r4delay FFT with partial normalization (r4dwpn). */

	else if (gwdata->FFT_TYPE == FFT_TYPE_RADIX_4_DWPN) {

/* Build sin/cos and premultiplier tables used in pass 1 of two pass FFTs. */
/* For best prefetching, make sure tables remain on 128-byte boundaries */

		ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
		tables = r4dwpn_build_pass1_table (gwdata, tables);
		ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
		asm_data->sincos1 = tables;
		tables = r4delay_build_fixed_premult_table (gwdata, tables);
		ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
		asm_data->sincos2 = tables;
		tables = r4delay_build_fixed_pass1_table (gwdata, tables);

/* Build the sin/cos table used in complex pass 2 blocks */
/* The pass 2 tables are the same as for a traditional radix-4 FFT */		

		ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
		asm_data->xsincos_complex = tables;
		tables = r4_build_pass2_complex_table (gwdata, tables);
		ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
		asm_data->sincos3 = tables;
		tables = r4_build_pass2_real_table (gwdata, tables);

/* Allocate a table for carries.  Init with XMM_BIGVAL.  For best */
/* distribution of data in the L2 cache, make this table contiguous */
/* with all other data used in the first pass (scratch area, normalization */
/* tables, etc.)  Note that we put the tables that are only partly loaded */
/* (column multipliers and big/lit table after the tables that are */
/* loaded throughout the first pass. */

		if (gwdata->PASS2_SIZE) {
			int	i, carry_table_size;
			double	xmm_bigval;
			ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
			asm_data->carries = tables;
			carry_table_size = gwdata->FFTLEN / (gwdata->PASS2_SIZE / 2);
			xmm_bigval = 3.0 * 131072.0 * 131072.0 * 131072.0;
			for (i = 0; i < carry_table_size; i++) *tables++ = xmm_bigval;
			tables += (16 - (tables - gwdata->gwnum_memory)) & 15;
		}

/* Build the group muliplier normalization table.  Keep this table */
/* contiguous with other data used in pass 1. */

		ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
		asm_data->norm_grp_mults = tables;
		tables = r4dwpn_build_norm_table (gwdata, tables, 0);

/* Reserve room for the pass 1 scratch area. */

		if (gwdata->SCRATCH_SIZE) {
			ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
			asm_data->scratch_area = tables;
			tables = (double *) ((char *) tables + gwdata->SCRATCH_SIZE);
		}

/* Build the table of big vs. little flags.  This cannot be last table */
/* built as xnorm_2d macro reads 2 bytes past the end of the array. */

		ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
		asm_data->norm_biglit_array = tables;
		tables = r4dwpn_build_biglit_table (gwdata, tables);
		tables += (16 - (tables - gwdata->gwnum_memory)) & 15;

/* Build the column normalization multiplier table. */

		asm_data->norm_col_mults = NULL;
	}

/* Initialize tables for the home-grown SSE2 FFT code. */

	else if (gwdata->cpu_flags & CPU_SSE2) {

/* Build sin/cos and premultiplier tables used in pass 2 of two pass FFTs */
/* Remember that pass2_size is the number of complex values in a pass 2 */
/* section, but build_sin_cos_table wants the number of reals in a section. */
/* However, we build a 1/4-sized table by mixing some of the sin/cos */
/* data into the premultiplier table.  So, divide pass2_size by 2 instead of */
/* multiplying pass2_size by 2. */

/* For best prefetching, make sure tables remain on 128-byte boundaries */

		if (gwdata->PASS2_SIZE) {
			ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
			asm_data->pass2_premults = tables;
			tables = hg_build_premult_table (gwdata, tables);

/* Build the rest of the tables */

			ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
			asm_data->xsincos_complex = tables;
			tables = hg_build_sin_cos_table (tables, gwdata->PASS2_SIZE/2, 0, 1);

			if (!gwdata->ALL_COMPLEX_FFT) {
				ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
				asm_data->sincos6 = tables;
				tables = hg_build_sin_cos_table (tables, gwdata->PASS2_SIZE * 4, 1, 2);
				asm_data->sincos7 = tables;
				tables = hg_build_sin_cos_table (tables, gwdata->PASS2_SIZE, 1, 1);
				tables += (16 - (tables - gwdata->gwnum_memory)) & 15;
			}

			asm_data->sincos8 = asm_data->sincos7;
			asm_data->sincos9 = asm_data->sincos7;
			asm_data->sincos10 = asm_data->sincos7;
			asm_data->sincos11 = asm_data->sincos7;
			asm_data->sincos12 = asm_data->sincos7;
		}

/* Allocate a table for carries.  Init with XMM_BIGVAL.  For best */
/* distribution of data in the L2 cache, make this table contiguous */
/* with all other data used in the first pass (scratch area, normalization */
/* tables, etc.)  Note that we put the tables that are only partly loaded */
/* (column multipliers and big/lit table after the tables that are */
/* loaded throughout the first pass. */

		if (gwdata->PASS2_SIZE) {
			int	i, carry_table_size;
			double	xmm_bigval;
			ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
			asm_data->carries = tables;
			carry_table_size = gwdata->FFTLEN / (gwdata->PASS2_SIZE / 2);
			xmm_bigval = 3.0 * 131072.0 * 131072.0 * 131072.0;
			for (i = 0; i < carry_table_size; i++)
				*tables++ = xmm_bigval;
			tables += (16 - (tables - gwdata->gwnum_memory)) & 15;
		}

/* Build the group muliplier normalization table.  Keep this table */
/* contiguous with other data used in pass 1. */

		ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
		asm_data->norm_grp_mults = tables;
		tables = hg_build_norm_table (gwdata, tables, 0);

/* Build the plus1-pre-multiplier table (complex weights applied when c > 0 */
/* and we are doing a all-complex FFT rather than emulating it with a */
/* zero-padded FFT. */

		if (gwdata->ALL_COMPLEX_FFT) {
			ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
			asm_data->plus1_premults = tables;
			tables = hg_build_plus1_table (gwdata, tables);
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
		tables = hg_build_sin_cos_table (tables, pass1_size, !gwdata->ALL_COMPLEX_FFT, gwdata->PASS2_SIZE == 0 ? 2 : 1);

		if (gwdata->PASS2_SIZE && pass1_size == pow_two_above_or_equal (pass1_size))
			asm_data->sincos2 = asm_data->sincos1;
		else {
			asm_data->sincos2 = tables;
			tables = hg_build_sin_cos_table (tables, pass1_size/4, !gwdata->ALL_COMPLEX_FFT, 1);
		}

		if (pass1_size == pow_two_above_or_equal (pass1_size)) {
			asm_data->sincos3 = asm_data->sincos2;
			asm_data->sincos4 = asm_data->sincos2;
			asm_data->sincos5 = asm_data->sincos2;
		} else {
			asm_data->sincos3 = tables;
			tables = hg_build_sin_cos_table (tables, pass1_size/16, !gwdata->ALL_COMPLEX_FFT, 1);
			asm_data->sincos4 = tables;
			tables = hg_build_sin_cos_table (tables, pass1_size/64, !gwdata->ALL_COMPLEX_FFT, 1);
			asm_data->sincos5 = tables;
			tables = hg_build_sin_cos_table (tables, pass1_size/256, !gwdata->ALL_COMPLEX_FFT, 1);
		}
		tables += (16 - (tables - gwdata->gwnum_memory)) & 15;

/* Build the table of big vs. little flags.  This cannot be last table */
/* built as xnorm_2d macro reads 2 bytes past the end of the array. */

		ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
		asm_data->norm_biglit_array = tables;
		tables = hg_build_biglit_table (gwdata, tables);
		tables += (16 - (tables - gwdata->gwnum_memory)) & 15;

/* Build the column normalization multiplier table. */

		ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
		asm_data->norm_col_mults = tables;
		tables = hg_build_norm_table (gwdata, tables, 1);
	}

/* Initialze table for the x87 assembly code. */

#ifndef X86_64

	else {

/* Allocate a table for carries.  Init with zero.  For best */
/* distribution of data in the L2 cache, make this table contiguous */
/* with the scratch area which is also used in the first pass. */

		if (gwdata->PASS2_SIZE) {
			int	i, carry_table_size;
			asm_data->carries = tables;
			carry_table_size = gwdata->FFTLEN / gwdata->PASS2_SIZE;
			for (i = 0; i < carry_table_size; i++) *tables++ = 0.0;
		}

/* Reserve room for the pass 1 scratch area. */

		asm_data->scratch_area = tables;
		if (gwdata->SCRATCH_SIZE)
			tables = (double *) ((char *) tables + gwdata->SCRATCH_SIZE);

/* Build the group muliplier normalization table.  Keep this table */
/* contiguous with other data used in pass 1. */

		asm_data->norm_grp_mults = tables;
		tables = x87_build_norm_table (gwdata, tables, 0);

/* Build sin/cos tables used in pass 1.  If FFTLEN is a power of two, */
/* many of the sin/cos tables can be shared. */

		asm_data->sincos1 = tables;
		tables = x87_build_sin_cos_table (tables, pass1_size, !gwdata->ALL_COMPLEX_FFT);

		if (pass1_size == pow_two_above_or_equal (pass1_size)) {
			asm_data->sincos2 = asm_data->sincos1;
			asm_data->sincos3 = asm_data->sincos1;
			asm_data->sincos4 = asm_data->sincos1;
			asm_data->sincos5 = asm_data->sincos1;
		} else {
			asm_data->sincos2 = tables;
			tables = x87_build_sin_cos_table (tables, pass1_size/4, !gwdata->ALL_COMPLEX_FFT);
			asm_data->sincos3 = tables;
			tables = x87_build_sin_cos_table (tables, pass1_size/16, !gwdata->ALL_COMPLEX_FFT);
			asm_data->sincos4 = tables;
			tables = x87_build_sin_cos_table (tables, pass1_size/64, !gwdata->ALL_COMPLEX_FFT);
			asm_data->sincos5 = tables;
			tables = x87_build_sin_cos_table (tables, pass1_size/256, !gwdata->ALL_COMPLEX_FFT);
		}

/* Build sin/cos and premultiplier tables used in pass 2 of two pass FFTs */
/* Remember that pass2_size is the number of complex values in a pass 2 */
/* section, but x87_build_sin_cos_table wants the number of reals in */
/* a section. */

		if (gwdata->PASS2_SIZE) {
			asm_data->pass2_premults = tables;
			tables = x87_build_premult_table (gwdata, tables);
			asm_data->xsincos_complex = tables;
			tables = x87_build_sin_cos_table (tables, gwdata->PASS2_SIZE*2, 0);

			if (!gwdata->ALL_COMPLEX_FFT) {
				asm_data->sincos6 = tables;
				tables = x87_build_sin_cos_table (tables, gwdata->PASS2_SIZE*2, 1);
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
			tables = x87_build_plus1_table (gwdata, tables);
		}

/* Build the column normalization multiplier table. */

		asm_data->norm_col_mults = tables;
		tables = x87_build_norm_table (gwdata, tables, 1);

/* Build the table of big vs. little flags. */

		asm_data->norm_biglit_array = tables;
		tables = x87_build_biglit_table (gwdata, tables);
	}

#endif

/* Finish verifying table size */

#ifdef GDEBUG
	{char buf[80];
	intptr_t mem = (intptr_t) tables - (intptr_t) gwdata->gwnum_memory;
	if (mem != mem_needed) {
		sprintf (buf, "FFTlen: %d, mem needed should be: %d\n",
			 (int) gwdata->FFTLEN, (int) (mem - gwdata->SCRATCH_SIZE));
		OutputBoth(0,buf);}}
#endif

/* Do more assembly initialization */

	asm_data->zero_fft = 0;
	asm_data->const_fft = 0;
	asm_data->ZPAD0_7[7] = 0.0;
	asm_data->zpad_addr = &asm_data->ZPAD0_7[0];

/* Copy the count of cache lines grouped in pass 1.  This affects how */
/* we build the normalization tables.  Note that cache line sizes are */
/* different in the x87 (16 bytes) and SSE2 code (64 bytes).  In the assembly */
/* code this is called clm or cache_line_multiplier. */
	
	asm_data->cache_line_multiplier = gwdata->PASS1_CACHE_LINES;

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

/* Newer FFTs use a clever mechanism to reduce big/lit flags data. */
/* Output the valid combinations that were precomputed in r4dwpn_build_biglit_table. */

	if (gwdata->FFT_TYPE == FFT_TYPE_RADIX_4_DWPN) {
		int	i;
		double	inv1, inv2, bmax1, bmax2;

		inv1 = asm_data->XMM_LIMIT_INVERSE[0];
		inv2 = asm_data->XMM_LIMIT_INVERSE[7];
		bmax1 = asm_data->XMM_LIMIT_BIGMAX[0];
		bmax2 = asm_data->XMM_LIMIT_BIGMAX[7];
		for (i = 0; i < 48; i++) {
			asm_data->XMM_LIMIT_INVERSE[i+i] = (((char *)gwdata->ASM_TIMERS)[i] & 1) ? inv2 : inv1;
			asm_data->XMM_LIMIT_INVERSE[i+i+1] = (((char *)gwdata->ASM_TIMERS)[i] & 2) ? inv2 : inv1;
			asm_data->XMM_LIMIT_BIGMAX[i+i] = (((char *)gwdata->ASM_TIMERS)[i] & 1) ? bmax2 : bmax1;
			asm_data->XMM_LIMIT_BIGMAX[i+i+1] = (((char *)gwdata->ASM_TIMERS)[i] & 2) ? bmax2 : bmax1;
			asm_data->XMM_LIMIT_BIGMAX_NEG[i+i] = (((char *)gwdata->ASM_TIMERS)[i] & 1) ? -bmax2 : -bmax1;
			asm_data->XMM_LIMIT_BIGMAX_NEG[i+i+1] = (((char *)gwdata->ASM_TIMERS)[i] & 2) ? -bmax2 : -bmax1;
		}
	}

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
	asm_data->XMM_P809[0] = asm_data->XMM_P809[1] = -asm_values[6];
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
	asm_data->XMM_P223[0] = asm_data->XMM_P223[1] = -asm_values[17];
	asm_data->M901 = asm_values[18];
	asm_data->XMM_P901[0] = asm_data->XMM_P901[1] = -asm_values[18];
	asm_data->M691 = asm_values[19];
	asm_data->XMM_P924[0] = asm_data->XMM_P924[1] = asm_values[20];
	asm_data->XMM_P383[0] = asm_data->XMM_P383[1] = asm_values[21];
	asm_data->XMM_P782[0] = asm_data->XMM_P782[1] = asm_values[22];
	asm_data->XMM_P434[0] = asm_data->XMM_P434[1] = asm_values[23];

/* Non-SSE2 initialization. Calculate constants used in two pass FFTs. */
/* Foremost is the pass 1 blkdst and normalize blkdst for auxillary mult */
/* routines.  The two values are the same except for larger FFTs which */
/* use a scratch area. */

	if (! (gwdata->cpu_flags & CPU_SSE2)) {
		if (gwdata->PASS2_SIZE) {
			/* Pass 2 data size: pass2_size complex numbers */
			/* Compute number of 4KB pages, used in normalized */
			/* add/sub */
			asm_data->normval4 = (gwdata->PASS2_SIZE * 16) >> 12;

			/* pad 64 bytes every 4KB + 64 bytes between blocks */
			asm_data->pass1blkdst =
			asm_data->normblkdst =
				asm_data->normval4 * (4096+64) + 64;
			asm_data->normblkdst8 = 0;

			if (gwdata->SCRATCH_SIZE) {
				/* Compute scratch area normblkdst */
				asm_data->normblkdst =
					asm_data->cache_line_multiplier * 32;

				/* Only larger pass1's have padding */
				if (asm_data->addcount1 >= 128)
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

		/* Data size = pass2_size complex numbers * 2 (for SSE2) */
		/* Calculate number of 8KB pages */
		asm_data->normval4 = (gwdata->PASS2_SIZE * 16 * 2) >> 13;

		/* Pass 1 block distance includes 128 pad bytes every 8KB */
		asm_data->pass2gapsize = gwdata->PASS2GAPSIZE;
		asm_data->pass1blkdst = gwdata->PASS2_SIZE * 16 * 2;
		asm_data->pass1blkdst += (asm_data->pass1blkdst >> 13) * 128;
		asm_data->pass1blkdst += asm_data->pass2gapsize;

		/* Calculate normblk distances */
		if (gwdata->SCRATCH_SIZE == 0) {
			/* Normblkdst = pass1blkdst - clm*64 */
			asm_data->normblkdst =
				asm_data->pass1blkdst -
				asm_data->cache_line_multiplier * 64;
			asm_data->normblkdst8 = 0;
		} else if ((gwdata->FFT_TYPE == FFT_TYPE_RADIX_4 ||
			    gwdata->FFT_TYPE == FFT_TYPE_RADIX_4_DELAYED ||
			    gwdata->FFT_TYPE == FFT_TYPE_RADIX_4_DWPN) &&
			   (gwdata->FFTLEN / gwdata->PASS2_SIZE == 80 ||
			    gwdata->FFTLEN / gwdata->PASS2_SIZE == 96 ||
			    gwdata->FFTLEN / gwdata->PASS2_SIZE == 112 ||
			    gwdata->FFTLEN / gwdata->PASS2_SIZE == 224)) {
			/* Small pass 1's with PFA have no padding every 8 clmblkdsts */
   			asm_data->normblkdst = 0;
			asm_data->normblkdst8 = 0;
		} else {
			/* Pad in clmblkdst is zero, clmblkdst8 is 128 */
			asm_data->normblkdst = 0;
			asm_data->normblkdst8 = 128;
		}

		if (gwdata->PASS2_SIZE) {
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

			/* This FFT type uses only half the big/lit data */
			if (gwdata->FFT_TYPE == FFT_TYPE_RADIX_4_DWPN) {
				asm_data->normval2 /= 2;
				asm_data->normval3 /= 2;
			}
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
			0.6 * log2 (gwdata->FFTLEN) + log2 (k) + 1.7 * log2 (abs (c)))
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
		asm_data->CARRY_ADJUST4 = pow ((double) b, num_b_in_second_top_word);
		asm_data->CARRY_ADJUST6 = pow ((double) b, num_b_in_third_top_word);
		if (gwdata->FFT_TYPE == FFT_TYPE_RADIX_4_DWPN) {
			asm_data->CARRY_ADJUST3 = gwfft_weight (gwdata->dd_data, gwdata->FFTLEN-1) /
						  gwfft_weight (gwdata->dd_data, dwpn_col (gwdata, gwdata->FFTLEN-1));
			asm_data->CARRY_ADJUST5 = gwfft_weight (gwdata->dd_data, gwdata->FFTLEN-2) /
						  gwfft_weight (gwdata->dd_data, dwpn_col (gwdata, gwdata->FFTLEN-2));
			asm_data->CARRY_ADJUST7 = gwfft_weight (gwdata->dd_data, gwdata->FFTLEN-3) /
						  gwfft_weight (gwdata->dd_data, dwpn_col (gwdata, gwdata->FFTLEN-3));
		} else {
			asm_data->CARRY_ADJUST3 = gwfft_weight (gwdata->dd_data, gwdata->FFTLEN-1);
			asm_data->CARRY_ADJUST5 = gwfft_weight (gwdata->dd_data, gwdata->FFTLEN-2);
			asm_data->CARRY_ADJUST7 = gwfft_weight (gwdata->dd_data, gwdata->FFTLEN-3);
		}

		// It's too hard to upgrade the old x87 code to match the SSE2 code
		// So, for x87 generate the same constants used prior to version 25.11
		if (! (gwdata->cpu_flags & CPU_SSE2)) {
			if (gwdata->PASS2_SIZE)
				asm_data->CARRY_ADJUST4 *= asm_data->CARRY_ADJUST5;
			else
				asm_data->CARRY_ADJUST6 *= asm_data->CARRY_ADJUST7;
		}

/* In two-pass FFTs, we only support tweaking the top two words. In one-pass FFTs, */
/* we adjust the top three words.  Make sure this works.  A test case that fails: */
/* 489539*3^72778+1.  We should consider supporting tweaking the top 3 words. */

		ASSERTG ((gwdata->PASS2_SIZE &&
			  num_b_in_k <= num_b_in_top_word + num_b_in_second_top_word) ||
			 (!gwdata->PASS2_SIZE &&
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

/* Pre-compute the adjustments to copying the 7 words around the halfway point */
/* In a radix-4 delay with partial normalization, we need to apply an adjustment */
/* so that the copied words are fully normalized. */

		if (gwdata->FFT_TYPE == FFT_TYPE_RADIX_4_DWPN) {
			int	i;
			for (i = 0; i < 7; i++) {
				gwdata->ZPAD_COPY7_ADJUST[i] = gwfft_weight (gwdata->dd_data, dwpn_col (gwdata, gwdata->FFTLEN / 2 - 3 + i));
			}
			/* Adjustments for ZPAD0 - ZPAD6 */
			for (i = 0; i < 7; i++) {
				gwdata->ZPAD_0_7_ADJUST[i] = gwfft_weight_inverse (gwdata->dd_data, dwpn_col (gwdata, i));
			}
		} else {
			int	i;
			for (i = 0; i < 7; i++) gwdata->ZPAD_COPY7_ADJUST[i] = 1.0;
		}
	}

/* Set the procedure pointers from the proc tables */

#ifndef X86_64
	if (! (gwdata->cpu_flags & CPU_SSE2)) {
		memcpy (gwdata->GWPROCPTRS+1, &x87_aux_prctab[gwdata->PASS2_SIZE ? 9 : 0], 9 * sizeof (void *));
		gwdata->GWPROCPTRS[norm_routines] = x87_prctab[x87_prctab_index (gwdata, 0, 0, 0)];  // No error, no mulbyconst
		gwdata->GWPROCPTRS[norm_routines+1] = x87_prctab[x87_prctab_index (gwdata, 0, 1, 0)];  // Error, no mulbyconst
		gwdata->GWPROCPTRS[norm_routines+2] = x87_prctab[x87_prctab_index (gwdata, 0, 0, 1)];  // No error, mulbyconst
		gwdata->GWPROCPTRS[norm_routines+3] = x87_prctab[x87_prctab_index (gwdata, 0, 1, 1)];  // Error, mulbyconst
		gwdata->GWPROCPTRS[zerohigh_routines] = x87_prctab[x87_prctab_index (gwdata, 1, 0, 0)];  // Zero, no error
		gwdata->GWPROCPTRS[zerohigh_routines+1] = x87_prctab[x87_prctab_index (gwdata, 1, 1, 0)];  // Zero, error
	} else
#endif
	{
		memcpy (gwdata->GWPROCPTRS+1, &sse2_aux_prctab[gwdata->FFT_TYPE == FFT_TYPE_RADIX_4_DWPN ? 18 : gwdata->PASS2_SIZE ? 9 : 0], 9 * sizeof (void *));
		gwdata->GWPROCPTRS[norm_routines] = sse2_prctab[sse2_prctab_index (gwdata, 0, 0, 0)];  // No error, no mulbyconst
		gwdata->GWPROCPTRS[norm_routines+1] = sse2_prctab[sse2_prctab_index (gwdata, 0, 1, 0)];  // Error, no mulbyconst
		gwdata->GWPROCPTRS[norm_routines+2] = sse2_prctab[sse2_prctab_index (gwdata, 0, 0, 1)];  // No error, mulbyconst
		gwdata->GWPROCPTRS[norm_routines+3] = sse2_prctab[sse2_prctab_index (gwdata, 0, 1, 1)];  // Error, mulbyconst
		gwdata->GWPROCPTRS[zerohigh_routines] = sse2_prctab[sse2_prctab_index (gwdata, 1, 0, 0)];  // Zero, no error
		gwdata->GWPROCPTRS[zerohigh_routines+1] = sse2_prctab[sse2_prctab_index (gwdata, 1, 1, 0)];  // Zero, error
	}

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
/* and for radix-4 delay with partial normalization FFTs. */

	if (gwdata->ALL_COMPLEX_FFT || gwdata->FFT_TYPE == FFT_TYPE_RADIX_4_DWPN)
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
		if (gwdata->PASS2_SIZE == 0) {		/* One pass */
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
		if (gwdata->PASS2_SIZE == 0)		/* One pass */
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
	gwhandle *gwdata;
	double	*srcarg, *destarg;

	gwdata = asm_data->gwdata;
	destarg = (double *) asm_data->DESTARG;
	srcarg = (double *) ((char *) destarg + asm_data->DIST_TO_FFTSRCARG);
	destarg[-5] = srcarg[4] * gwdata->ZPAD_COPY7_ADJUST[3];	/* Copy 1st word above halfway point */
	destarg[-6] = srcarg[12] * gwdata->ZPAD_COPY7_ADJUST[4]; /* Copy 2nd word */
	destarg[-7] = srcarg[20] * gwdata->ZPAD_COPY7_ADJUST[5]; /* Copy 3rd word */
	destarg[-8] = srcarg[28] * gwdata->ZPAD_COPY7_ADJUST[6]; /* Copy 4th word */
	destarg[-9] = * (double *)	/* Copy 1st word below halfway point */
		((char *) srcarg + asm_data->HIGH_WORD1_OFFSET) * gwdata->ZPAD_COPY7_ADJUST[2];
	destarg[-10] = * (double *)	/* Copy 2nd word below */
		((char *) srcarg + asm_data->HIGH_WORD2_OFFSET) * gwdata->ZPAD_COPY7_ADJUST[1];
	destarg[-11] = * (double *)	/* Copy 3rd word below */
		((char *) srcarg + asm_data->HIGH_WORD3_OFFSET) * gwdata->ZPAD_COPY7_ADJUST[0];
}

/* When POSTFFT is set, we must copy the 7 words at two different spots. */
/* These two routines copy the four values above the half-way point after */
/* carries have been propagated and copy the three words just below the */
/* half-way point right after the last NORMRTN has been called. */

void xcopy_4_words (
	struct gwasm_data *asm_data)
{
	gwhandle *gwdata;
	double	*srcarg;

	gwdata = asm_data->gwdata;
	srcarg = (double *) asm_data->DESTARG;
	srcarg[-5] = srcarg[4] * gwdata->ZPAD_COPY7_ADJUST[3]; /* Copy 1st word above halfway point */
	srcarg[-6] = srcarg[12] * gwdata->ZPAD_COPY7_ADJUST[4]; /* Copy 2nd word */
	srcarg[-7] = srcarg[20] * gwdata->ZPAD_COPY7_ADJUST[5]; /* Copy 3rd word */
	srcarg[-8] = srcarg[28] * gwdata->ZPAD_COPY7_ADJUST[6]; /* Copy 4th word */
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
					 asm_data->HIGH_SCRATCH1_OFFSET) * gwdata->ZPAD_COPY7_ADJUST[2];
			srcarg[-10] = * (double *)
					((char *) asm_data->scratch_area +
					 asm_data->HIGH_SCRATCH2_OFFSET) * gwdata->ZPAD_COPY7_ADJUST[1];
		} else {
			srcarg[-9] = * (double *)
					((char *) srcarg +
					 asm_data->HIGH_WORD1_OFFSET) * gwdata->ZPAD_COPY7_ADJUST[2];
			srcarg[-10] = * (double *)
					((char *) srcarg +
					 asm_data->HIGH_WORD2_OFFSET) * gwdata->ZPAD_COPY7_ADJUST[1];
		}
	}

/* Copy 3rd word below the halfway point */

	if (asm_data->this_block <= gwdata->num_pass1_blocks - 3 &&
	    asm_data->this_block + asm_data->cache_line_multiplier >
						gwdata->num_pass1_blocks - 3) {
		if (gwdata->SCRATCH_SIZE) {
			srcarg[-11] = * (double *)
					((char *) asm_data->scratch_area +
					 asm_data->HIGH_SCRATCH3_OFFSET) * gwdata->ZPAD_COPY7_ADJUST[0];
		} else {
			srcarg[-11] = * (double *)
					((char *) srcarg +
					 asm_data->HIGH_WORD3_OFFSET) * gwdata->ZPAD_COPY7_ADJUST[0];
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

/* Calculate pass 1 sin/cos/premult address (for those FFTs that do not use */
/* the same sin/cos table for every pass 1 group. */

__inline void *pass1_premult_addr (
	gwhandle *gwdata,
	unsigned long block)
{
	if (gwdata->FFT_TYPE == FFT_TYPE_RADIX_4 ||
	    gwdata->FFT_TYPE == FFT_TYPE_RADIX_4_DELAYED ||
	    gwdata->FFT_TYPE == FFT_TYPE_RADIX_4_DWPN)
		return ((char *) gwdata->adjusted_pass1_premults + block * gwdata->pass1_premult_block_size);
	return (NULL);
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
	asm_data->premult_prefetch =
		pass1_premult_addr (gwdata, asm_data->next_block);
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
			asm_data->premult_prefetch =
				pass1_premult_addr (gwdata, asm_data->next_block);
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
		carry_table_size = (gwdata->FFTLEN / gwdata->PASS2_SIZE) << 1;
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
	asm_data->premult_addr = pass1_premult_addr (gwdata, 0);
	if (gwdata->pass1_state == 0) {
		gwdata->next_block = asm_data->cache_line_multiplier;
		pass1_state0_assign_next_block (gwdata, asm_data);
	} else {
		asm_data->next_block = asm_data->cache_line_multiplier;
		asm_data->data_prefetch =
			pass1_data_addr (asm_data, asm_data->next_block);
		asm_data->premult_prefetch =
			pass1_premult_addr (gwdata, asm_data->next_block);
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
/* each thread essentially shares these values. */

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
		carry_table_size = (gwdata->FFTLEN / gwdata->PASS2_SIZE) << 1;
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
	asm_data->premult_addr = pass1_premult_addr (gwdata, 0);
	if (gwdata->pass1_state == 0) {
		gwdata->next_block = asm_data->cache_line_multiplier;
		pass1_state0_assign_next_block (gwdata, asm_data);
	} else {
		asm_data->next_block = asm_data->cache_line_multiplier *
					gwdata->num_threads;
		asm_data->data_prefetch =
			pass1_data_addr (asm_data, asm_data->next_block);
		asm_data->premult_prefetch =
			pass1_premult_addr (gwdata, asm_data->next_block);
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
	gwhandle *gwdata;

/* Get pointer to the gwdata structure */

	gwdata = asm_data->gwdata;

/* In the radix-4 delay with partial normalization zero-padded FFT case prior to */
/* normalizing the first data block, apply the partial normalization multipliers */
/* to ZPAD0 - ZPAD6. */

	if (asm_data->this_block == 0 && gwdata->ZERO_PADDED_FFT && gwdata->FFT_TYPE == FFT_TYPE_RADIX_4_DWPN) {
		int	i;
		for (i = 0; i < 7; i++) asm_data->ZPAD0_7[i] *= gwdata->ZPAD_0_7_ADJUST[i];
	}
}

void pass1_wait_for_carries_mt (
	struct gwasm_data *asm_data)
{
	gwhandle *gwdata;
	struct gwasm_data *main_thread_asm_data;

/* Get pointer to the gwdata structure */

	gwdata = asm_data->gwdata;

/* In the radix-4 delay with partial normalization zero-padded FFT case prior to */
/* normalizing the first data block, apply the partial normalization multipliers */
/* to ZPAD0 - ZPAD6. */

	if (asm_data->this_block == 0 && gwdata->ZERO_PADDED_FFT && gwdata->FFT_TYPE == FFT_TYPE_RADIX_4_DWPN) {
		int	i;
		for (i = 0; i < 7; i++) asm_data->ZPAD0_7[i] *= gwdata->ZPAD_0_7_ADJUST[i];
	}

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
		asm_data->premult_addr = asm_data->premult_prefetch;
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
	asm_data->premult_addr = asm_data->premult_prefetch;
	if (gwdata->pass1_state == 1) {
		unsigned long next_block;
		next_block = asm_data->this_block +
				asm_data->cache_line_multiplier;
		if (next_block < gwdata->num_pass1_blocks) {
			asm_data->next_block = next_block;
			asm_data->data_prefetch =
				pass1_data_addr (asm_data, next_block);
			asm_data->premult_prefetch =
				pass1_premult_addr (gwdata, next_block);
			return (PASS1_DO_MORE_INVERSE_FFT);
		}
	}

/* Either there is no next block or it will be a POSTFFT block */

	if (asm_data->POSTFFT && gwdata->next_block < gwdata->num_postfft_blocks) {
		asm_data->next_block = gwdata->next_block;
		asm_data->data_prefetch =
			pass1_data_addr (asm_data, asm_data->next_block);
		asm_data->premult_prefetch =
			pass1_premult_addr (gwdata, asm_data->next_block);
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
		asm_data->premult_addr = asm_data->premult_prefetch;
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
	asm_data->premult_addr = asm_data->premult_prefetch;
	if (gwdata->pass1_state == 1) {
		unsigned long next_block;
		next_block = asm_data->this_block +
				asm_data->cache_line_multiplier *
				gwdata->num_threads;
		if (next_block < gwdata->num_pass1_blocks) {
			asm_data->next_block = next_block;
			asm_data->data_prefetch =
				pass1_data_addr (asm_data, next_block);
			asm_data->premult_prefetch =
				pass1_premult_addr (gwdata, next_block);
			gwmutex_unlock (&gwdata->thread_lock);
			return (PASS1_DO_MORE_INVERSE_FFT);
		}
	}

/* Either there is no next block or it will be a POSTFFT block */

	if (asm_data->POSTFFT && gwdata->next_block < gwdata->num_postfft_blocks) {
		asm_data->next_block = gwdata->next_block;
		asm_data->data_prefetch =
			pass1_data_addr (asm_data, asm_data->next_block);
		asm_data->premult_prefetch =
			pass1_premult_addr (gwdata, asm_data->next_block);
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

/* Only two pass SSE2 FFTs support multi-threaded execution */

	if (gwdata->PASS2_SIZE == 0 || !(gwdata->cpu_flags & CPU_SSE2)) {
		gwdata->num_threads = 1;
		return (0);
	}

/* Get pointer to assembly structure */

	asm_data = (struct gwasm_data *) gwdata->asm_data;

/* Save gwdata pointer in asm_data so that C callback routines can */
/* access gwdata.  Set flag indicating this is the main thread. */

	asm_data->gwdata = gwdata;
	asm_data->thread_num = 0;

/* Init other variables */

	gwdata->num_pass1_blocks = gwdata->PASS2_SIZE >> 1;
	gwdata->num_pass2_blocks = (gwdata->FFTLEN / gwdata->PASS2_SIZE) >> 2;
	asm_data->last_pass1_block =
	    gwdata->num_pass1_blocks - asm_data->cache_line_multiplier;

//bug - minor optimization:  for small k/c we can probably set this to
//less than 8
	gwdata->num_postfft_blocks = 8;

/* Calculate the values used to compute pass 2 premultier pointers. */
/* This only happens for our home-grown FFTs.  Non-power-of-2 pass 2 */
/* sizes are not supported. */
/* We calculate an adjusted starting address of the premultiplier data */
/* so that both real and all-complex FFTs can use the same formula to */
/* calculate the proper address given a block number. */

	gwdata->pass2_premult_block_size =
		(gwdata->PASS2_SIZE == 256) ? 32 * 128 :
		(gwdata->PASS2_SIZE == 1024) ? 64 * 128 :
		(gwdata->PASS2_SIZE == 2048) ? 96 * 128 :
		(gwdata->PASS2_SIZE == 4096) ? 128 * 128 :
		(gwdata->PASS2_SIZE == 8192) ? 192 * 128 : 0;
	if (gwdata->ALL_COMPLEX_FFT)
		gwdata->adjusted_pass2_premults = asm_data->pass2_premults;
	else
		gwdata->adjusted_pass2_premults =
			(char *) asm_data->pass2_premults -
			gwdata->pass2_premult_block_size;

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

		if (gwdata->PASS2_SIZE == 0) {
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

/* Larger FFTs use two passes.  This the example for a length 64K FFT (pass2_size = 2K): */
/*	0	1K	16K	17K	32K	33K	48K	49K	*/
/*	1	...							*/
/*	...								*/
/*	1023	...							*/
/*	2K	...							*/
/*	...								*/
/* and PFA layouts are even funkier.					*/		

		else if (gwdata->FFT_TYPE == FFT_TYPE_HOME_GROWN) {
			sets = (fftlen / gwdata->PASS2_SIZE) >> 2;
			if (i >= (fftlen >> 1)) {
				i6 = 1;
				i -= (fftlen >> 1);
			} else
				i6 = 0;
			i1 = i % (gwdata->PASS2_SIZE >> 1);
			i = i / (gwdata->PASS2_SIZE >> 1);
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
			addr = i3 * (gwdata->PASS2_SIZE >> 1);
			addr = ((((((addr + i1) << 1) + i6) << 1) + i) << 1) + i2;
			addr = addr * sizeof (double);
			/* Now add 128 bytes every 8KB and one pass2gapsize */
			/* for every pass 2 block. */
			addr = addr + (addr >> 13) * 128 + i3 * gwdata->PASS2GAPSIZE;
		}

/* Newer traditional radix-4 large FFTs use don't have a special layout for PFA. */

		else {
			unsigned long top2, row;
			top2 = (i << 2) / fftlen; i -= (top2 * fftlen) >> 2;
			i1 = i % (gwdata->PASS2_SIZE >> 1); i = i / (gwdata->PASS2_SIZE >> 1);
			i2 = i & 1; i >>= 1;
			row = i * (gwdata->PASS2_SIZE >> 1) + i1;
			addr = (((row << 2) + top2) << 1) + i2;
			addr = addr * sizeof (double);
			/* Now add 128 bytes every 8KB and one pass2gapsize */
			/* for every pass 2 block. */
			addr = addr + (addr >> 13) * 128 + i * gwdata->PASS2GAPSIZE;
		}
	}

/* One pass x87 FFTs use a near flat memory model. */

	else if (gwdata->PASS2_SIZE == 0) {
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
		addr = i * 16 + i2 * 8 + (i >> 8) * 64 + (i / gwdata->PASS2_SIZE) * 64;
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
	return (gwdata->mem_needed + gwdata->SCRATCH_SIZE);
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

/* Handle r4dwpn FFTs which are only partially normalized */

	else if (gwdata->FFT_TYPE == FFT_TYPE_RADIX_4_DWPN) {
		val = val * gwfft_weight_sloppy (gwdata->dd_data, dwpn_col (gwdata, i))
		          * gwfft_weight_inverse_sloppy (gwdata->dd_data, i);
		if (val < -0.5)
			*retval = (long) (val - 0.5);
		else
			*retval = (long) (val + 0.5);
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

/* Handle r4dwpn FFTs which are only partially normalized */

	if (gwdata->FFT_TYPE == FFT_TYPE_RADIX_4_DWPN) {
		* addr (gwdata, g, i) = val * gwfft_weight_sloppy (gwdata->dd_data, i) / 
					      gwfft_weight_sloppy (gwdata->dd_data, dwpn_col (gwdata, i));
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

/* Map a gwerror code into human readable text */

void gwerror_text (
	gwhandle *gwdata,	/* Handle used in all gwnum calls */
	int	error_code,	/* Error code to turn ino text */
	char	*buf,		/* Buffer to write text to */
	int	buflen)		/* Sizeof the text buffer */
{
	char	localbuf[512];

/* Map error code to text */

	switch (error_code) {
	case GWERROR_VERSION:
		strcpy (localbuf, "Improperly compiled and linked.  Gwnum.h and FFT assembly code version numbers do not match.");
		break;
	case GWERROR_TOO_LARGE:
		strcpy (localbuf, "Number sent to gwsetup is too large for the FFTs to handle.");
		break;
	case GWERROR_K_TOO_SMALL:
		strcpy (localbuf, "Value of k in k*b^n+c is too small.  Values less than one are not supported.");
		break;
	case GWERROR_K_TOO_LARGE:
		strcpy (localbuf, "Value of k in k*b^n+c is too large.  Values greater than 2251799813685247 are not supported.");
		break;
	case GWERROR_MALLOC:
		strcpy (localbuf, "Unable to allocate memory.  One possible cause is the operating system's swap area is too small.");
		break;
	case GWERROR_VERSION_MISMATCH:
		strcpy (localbuf, "GWNUM_VERSION from gwinit call doesn't match GWNUM_VERSION when gwnum.c was compiled.  Recompile and relink.");
		break;
	case GWERROR_STRUCT_SIZE_MISMATCH:
		strcpy (localbuf, "Gwhandle structure size from gwinit call doesn't match size when gwnum.c was compiled.  Check compiler alignment switches, recompile and relink.");
		break;
	default:
		sprintf (localbuf, "Unknown gwnum error code: %d", error_code);
		break;
	}

/* Copy error message to caller's buffer */

	if ((int) strlen (localbuf) >= buflen) localbuf[buflen-1] = 0;
	strcpy (buf, localbuf);
}

/* Return a description of the FFT type chosen */

void gwfft_description (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	char	*buf)		/* Buffer to return string in */
{
	sprintf (buf, "%s%s%sFFT length %lu%s",
		 gwdata->ALL_COMPLEX_FFT ? "all-complex " :
		 gwdata->ZERO_PADDED_FFT ? "zero-padded " :
		 gwdata->GENERAL_MOD ? "generic reduction " : "",
		 (! (gwdata->cpu_flags & CPU_SSE2)) ? "x87 " :
		 (gwdata->PASS2_SIZE == 0) ? "" :
		 (gwdata->ARCH == 2) ? "Pentium4 " :
		 (gwdata->ARCH == 3) ? "Core2 " :
		 (gwdata->ARCH == 4) ? "AMD K8 " :
		 (gwdata->ARCH == 5) ? "AMD K10 " : " ",
		 (gwdata->PASS2_SIZE == 0) ? "" :
		 (gwdata->FFT_TYPE == 0) ? "type-0 " :
		 (gwdata->FFT_TYPE == 1) ? "type-1 " :
		 (gwdata->FFT_TYPE == 2) ? "type-2 " : "type-3 ",
		 gwdata->FFTLEN >= 1048576 && (gwdata->FFTLEN & 0xFFFFF) == 0 ?
			 gwdata->FFTLEN / 1048576 :
		 gwdata->FFTLEN >= 1024 && (gwdata->FFTLEN & 0x3FF) == 0 ?
			 gwdata->FFTLEN / 1024 : gwdata->FFTLEN,
		 gwdata->FFTLEN >= 1048576 && (gwdata->FFTLEN & 0xFFFFF) == 0 ? "M" :
		 gwdata->FFTLEN >= 1024 && (gwdata->FFTLEN & 0x3FF) == 0 ? "K" : "");

	if (gwdata->PASS2_SIZE) {
		int	p1size;
		char	p1buf[20], p2buf[20];

		p1size = gwdata->FFTLEN / gwdata->PASS2_SIZE;
		if (p1size % 1024) sprintf (p1buf, "%d", p1size);
		else sprintf (p1buf, "%dK", p1size / 1024);
		if (gwdata->PASS2_SIZE % 1024) sprintf (p2buf, "%d", (int) gwdata->PASS2_SIZE);
		else sprintf (p2buf, "%dK", (int) (gwdata->PASS2_SIZE / 1024));
		sprintf (buf + strlen (buf), ", Pass1=%s, Pass2=%s", p1buf, p2buf);
	}

	if (gwdata->num_threads > 1)
		sprintf (buf + strlen (buf), ", %d threads", (int) gwdata->num_threads);

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
			 c < 0 ? '-' : '+', (unsigned long) abs (c));
	else if (b == 2 && c == -1)
		sprintf (buf, "M%lu", n);
	else {
		unsigned long cnt, temp_n;
		for (cnt = 0, temp_n = n; !(temp_n & 1); temp_n >>= 1, cnt++);
		if (b == 2 && temp_n == 1 && c == 1)
			sprintf (buf, "F%lu", cnt);
		else
			sprintf (buf, "%lu^%lu%c%lu", b, n,
				 c < 0 ? '-' : '+', (unsigned long) abs (c));
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

/* Get pointer to fft info and return the maximum exponent for the FFT length */

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
	return (gwdata.mem_needed + gwdata.SCRATCH_SIZE);
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
		if (gwdata->PASS2_SIZE == 0) {
			row = asm_data->ADDIN_OFFSET & 31;
			if (row == 8) asm_data->ADDIN_OFFSET += 8;
			if (row == 16) asm_data->ADDIN_OFFSET -= 8;
		}

/* Factor in the blkdst value in xfft3.mac to compute the two pass */
/* SSE2 addin_offset. */

		else {
			unsigned long num_rows;

			num_rows = gwdata->PASS2_SIZE >> 1;
			row = word % num_rows;
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

				blkdst = addr_offset (gwdata, gwdata->PASS2_SIZE);
				row = asm_data->ADDIN_OFFSET / blkdst;
				asm_data->ADDIN_OFFSET -= row * blkdst;
				asm_data->ADDIN_OFFSET += row * (gwdata->PASS1_CACHE_LINES * 64);
				asm_data->ADDIN_OFFSET += (row >> 3) * (uint32_t) asm_data->normblkdst8;
			}
		}
	}

/* And now x87 FFTs also can use a scratch area.  Like the SSE2 code */
/* we have to convert the FFT data offsets for two-pass FFTs. */

	if (! (gwdata->cpu_flags & CPU_SSE2) && gwdata->PASS2_SIZE) {
		unsigned long num_cache_lines, cache_line;

		num_cache_lines = gwdata->PASS2_SIZE >> 1;
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

			blkdst = addr_offset (gwdata, gwdata->PASS2_SIZE);
			row = asm_data->ADDIN_OFFSET / blkdst;
			asm_data->ADDIN_OFFSET -= row * blkdst;
			asm_data->ADDIN_OFFSET += row * (gwdata->PASS1_CACHE_LINES * 32);

/* Handle the FFTs where clmblkdst32 is used */

			if (((gwdata->FFTLEN / gwdata->PASS2_SIZE) >> 1) >= 128)
				asm_data->ADDIN_OFFSET += (row >> 5) * 64;
		}
	}

/* Handle r4dwpn FFTs which are only partially normalized */

	if (gwdata->FFT_TYPE == FFT_TYPE_RADIX_4_DWPN) {
		asm_data->ADDIN_VALUE = val *
					gwfft_weight_sloppy (gwdata->dd_data, word) /
					gwfft_weight_sloppy (gwdata->dd_data, dwpn_col (gwdata, word));
	}

/* Set the addin value - multiply it by two-to-phi and FFTLEN/2/k. */

	else
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

	if (gwdata->GENERAL_MOD) limit = gwdata->GW_GEN_MOD_MAX + 3;
	else if (gwdata->ZERO_PADDED_FFT) limit = gwdata->FFTLEN / 2 + 4;
	else limit = gwdata->FFTLEN;

/* GENERAL_MOD has some strange cases we must handle.  In particular the */
/* last fft word translated can be 2^bits and the next word could be -1, */
/* this must be translated into zero, zero. */

	if (gwdata->GENERAL_MOD) {
		long	val, prev_val;
		while (limit < gwdata->FFTLEN) {
			err_code = get_fft_value (gwdata, gg, limit, &val);
			if (err_code) return (err_code);
			if (val == -1 || val == 0) break;
			limit++;
			ASSERTG (limit <= gwdata->FFTLEN / 2 + 2);
		}
		while (limit > 1) {		/* Find top word */
			err_code = get_fft_value (gwdata, gg, limit-1, &prev_val);
			if (err_code) return (err_code);
			if (val != prev_val || val < -1 || val > 0) break;
			limit--;
		}
		limit++;
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
		if (array == NULL) { memerr: err_code = GWERROR_MALLOC; goto err; }
		buf = (uint32_t *) malloc (limit * sizeof (uint32_t));
		if (buf == NULL) goto memerr;
		small_base = popg (&gwdata->gdata, limit);
		if (small_base == NULL) goto memerr;
		large_base = popg (&gwdata->gdata, limit);
		if (large_base == NULL) goto memerr;
		for (i = 0; i < limit; i++) {
			long	val;
			err_code = get_fft_value (gwdata, gg, i, &val);
			if (err_code) {
err:				free (array);
				free (buf);
				if (small_base != NULL) pushg (&gwdata->gdata, 1);
				if (large_base != NULL) pushg (&gwdata->gdata, 1);
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
		pushg (&gwdata->gdata, 2);
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
	asm_data->DESTARG = d;
	asm_data->DIST_TO_FFTSRCARG = 0;
	asm_data->DIST_TO_MULSRCARG = (intptr_t) s - (intptr_t) d;
	asm_data->ffttype = 3;
	gw_fft (gwdata, asm_data);
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

	asm_data->NORMRTN = gwdata->GWPROCPTRS[zerohigh_routines + (gwdata->NORMNUM & 1)];
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

	asm_data->DESTARG = d;
	asm_data->DIST_TO_FFTSRCARG = (intptr_t) s - (intptr_t) d;
	asm_data->DIST_TO_MULSRCARG = 0;
	asm_data->ffttype = 1;
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
	asm_data->DIST_TO_FFTSRCARG = 0;
	asm_data->DIST_TO_MULSRCARG = 0;
	asm_data->ffttype = 2;
	gw_fft (gwdata, asm_data);
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
	asm_data->DESTARG = d;
	asm_data->DIST_TO_MULSRCARG = (intptr_t) s - (intptr_t) d;
	asm_data->DIST_TO_FFTSRCARG = (intptr_t) s2 - (intptr_t) d;
	asm_data->ffttype = 4;
	gw_fft (gwdata, asm_data);
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

	if (gwdata->GWPROCPTRS[8] == NULL || gwdata->k > 1.0 ||
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
	    (* (double *) ((char *) g + gwdata->GW_GEN_MOD_MAX_OFFSET) <= -2.0 ||
	     * (double *) ((char *) g + gwdata->GW_GEN_MOD_MAX_OFFSET) > 0.0))
		emulate_mod (gwdata, g);
}
