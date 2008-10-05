/**************************************************************
 *
 *	giants.h
 *
 *	Header file for large-integer arithmetic library giants.c.
 *
 *	Updates:
 *          30 Apr 98  JF   USE_ASSEMBLER_MUL removed
 *          29 Apr 98  JF   Function prototypes cleaned up
 *	    20 Apr 97  RDW
 *
 *	c. 1997 Perfectly Scientific, Inc.
 *	c. 1998-2007 Just For Fun Software, Inc.
 *	All Rights Reserved.
 *
 **************************************************************/

#ifndef _GIANTS_H_
#define _GIANTS_H_

/* This is a C library.  If used in a C++ program, don't let the C++ */
/* compiler mangle names. */

#ifdef __cplusplus
extern "C" {
#endif

/* Include common definitions */

#include "gwcommon.h"

/**************************************************************
 *
 * Mul options
 *
 **************************************************************/

#define AUTO_MUL 0
#define GRAMMAR_MUL 1
#define FFT_MUL 2
#define KARAT_MUL 3

/**************************************************************
 *
 * Structure definitions
 *
 **************************************************************/

typedef struct
{
	 int	sign;
	 uint32_t *n;		/* ptr to array of longs */
#ifdef GDEBUG
	int	maxsize;
#endif
} giantstruct;

typedef giantstruct *giant;

#ifdef GDEBUG
#define setmaxsize(g,s)	(g)->maxsize = s
#else
#define setmaxsize(g,s)
#endif

/**************************************************************
 *
 * Function Prototypes
 *
 **************************************************************/

/* Creates a new giant */
giant 	newgiant(int numshorts);

/* Returns the bit-length n; e.g. n=7 returns 3. */
int 	bitlen(giant n);

/* Returns the value of the pos bit of g. */
#define bitval(g,pos)	((g)->n[(pos)>>5] & (1 << ((pos)&31)))

/* Returns whether g is zero. */
#define isZero(g)	((g)->sign == 0)

/* Returns whether g is one. */
#define isone(g)	(((g)->sign == 1) && ((g)->n[0] == 1))

/* Returns whether g is one. */
#define istwo(g)	(((g)->sign == 1) && ((g)->n[0] == 2))

/* Copies one giant to another. */
void 	gtog(giant src, giant dest);

/* Integer -> giant. */
#define setzero(g)	(g)->sign = 0
#define setone(g)	(g)->sign = (g)->n[0] = 1
void 	itog(int n, giant g);
void 	ultog(uint32_t n, giant g);

/* Long long -> giant. */
void 	ulltog(uint64_t n, giant g);

/* Double -> giant. */
void 	dbltog(double n, giant g);

/* Character string <-> giant. */
void 	ctog (char *s, giant g);
void 	gtoc (giant g, char *s, int sbufsize);

/* Returns the sign of g: -1, 0, 1. */
int	gsign(giant g);

/* Returns 1, 0, -1 as a>b, a=b, a<b. */
int	gcompg(giant a, giant b);

/* Set AUTO_MUL for automatic FFT crossover (this is the
 * default), set FFT_MUL for forced FFT multiply, set
 * GRAMMAR_MUL for forced grammar school multiply. */
void	setmulmode(int mode);

/**************************************************************
 *
 * Math Functions
 *
 **************************************************************/

/* g := -g. */
#define negg(g)		((g)->sign = -(g)->sign)

/* g := |g|. */
#define absg(g)		if ((g)->sign < 0) (g)->sign = -(g)->sign

/* g += i, with i an int (for compatability with old giants.h */
#define	iaddg(i,g)	sladdg ((int32_t)(i), g)

/* g += i, where i is an unsigned long */
void 	uladdg (uint32_t i, giant g);

/* g += i, where i is a signed long */
void 	sladdg (int32_t i, giant g);

/* b += a. */
void 	addg(giant a, giant b);

/* g -= i, where i is an unsigned long */
void 	ulsubg (uint32_t i, giant g);

/* b -= a. */
void 	subg(giant a, giant b);

/* g *= g. */
void	squareg(giant g);

/* b *= a. */
void	mulg(giant a, giant b);

/* g *= i, where i is an unsigned long */
void 	ulmulg (uint32_t i, giant g);
void 	imulg (int32_t i, giant g);

/* g *= i, where i is an unsigned long long */
void 	ullmulg (uint64_t i, giant g);

/* g *= n, where n is a double */
void 	dblmulg (double n, giant g);

/* num := num % den, any positive den. */
void 	modg(giant den, giant num);

/* num := [num/den], any positive den. */
void 	divg(giant den, giant num);
void 	dbldivg(double den, giant num);

/* Mask rightmost bits in g, same as g = g % 2^bits. */
void 	gmaskbits(int bits, giant g);

/* Shift g left by bits, introducing zeros on the right. */
void 	gshiftleft(int bits, giant g);

/* Shift g right by bits, losing bits on the right. */
#define	gshiftright(n,g)	{if (n) gtogshiftright (n, g, g);}
void 	gtogshiftright (int bits, giant src, giant dest);

/* If 1/x exists (mod n), then x := 1/x.  If
 * inverse does not exist, then x := - GCD(n, x). */
void 	invg(giant n, giant x);

/* General GCD, x:= GCD(n, x). */
void 	gcdg(giant n, giant x);

/* x becomes x^n, NO mod performed. */
void power (giant x, int n);
void powerg (giant x, giant n);

/* num := [num/den], any positive den. */
void 	powermod(giant x, int n, giant z);

/* x := x^n (mod z). */
void 	powermodg(giant x, giant n, giant z);


/* Alternate memory allocator functions and other global data. */
/* Prime95 uses these alternate routines to cut down on the amount of */
/* memory allocated by having giant temporaries and gwnum temporaries */
/* share the same cached memory pool. */

typedef struct {
	 void	*memblk;		/* Available memory block */
	 long	offset;			/* Offset into the block */
	 giant	prev;			/* Last giant allocated by popg */
} gstacknode;

typedef struct {
	 unsigned long blksize;		/* Gwnum size */
	 void*	(*allocate)(void *);	/* Gwnum allocator */
	 void	(*free)(void *,void *);	/* Gwnum free routine */
	 void*	handle;			/* Gwdata handle */
	 int	num_popgs;		/* Number of popg allocations */
	 gstacknode stack;		/* Linked list of popg allocations */
	 double	*ooura_sincos;
	 int	*ooura_ip;
	 int	ooura_fft_size;
	 giant	cur_recip;
	 giant	cur_den;
} ghandle;

void	init_ghandle (ghandle *);
void	term_ghandle (ghandle *);

/* Allocate and free data in a stack-like manner.  These routines have */
/* funky names for historical reasons. */

giant	popg(ghandle *,int);	/* Number of longs in data area */
void	pushg(ghandle *,int);/* Number of items to return to stack */

/* Interruptable and alternate memory allocator versions of some routines */
void	squaregi(ghandle *, giant g);
void	mulgi(ghandle *, giant a, giant b);
void 	modgi(ghandle *, giant den, giant num);
void 	divgi(ghandle *, giant den, giant num);
int 	invgi(ghandle *, int, giant n, giant x);
int 	gcdgi(ghandle *, int, giant n, giant x);


/* Low-level math routines the caller can use for multi-precision */
/* arithmetic */

extern void addhlp (uint32_t *res, uint32_t *carry, uint32_t val);
extern void subhlp (uint32_t *res, uint32_t *carry, uint32_t val);
extern void muladdhlp (uint32_t *res, uint32_t *carryl,
		       uint32_t *carryh, uint32_t val1, uint32_t val2);
extern void muladd2hlp (uint32_t *res, uint32_t *carryl,
			uint32_t *carryh, uint32_t val1, uint32_t val2);
extern void mulsubhlp (uint32_t *res, uint32_t *carryl,
		       uint32_t *carryh, uint32_t val1, uint32_t val2);

/* External routine pointers. */

extern int (*StopCheckRoutine)(int);

#ifdef __cplusplus
}
#endif

#endif