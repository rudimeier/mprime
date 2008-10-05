/*----------------------------------------------------------------------
| Copyright 1995-2008 Mersenne Research, Inc.  All rights reserved
| Author:  George Woltman
| Email: woltman@alum.mit.edu
|
| This file contains routines to determine the CPU type and speed.
| Plus, as a bonus, you get 3 routines to portably access the high
| resolution timer.
+---------------------------------------------------------------------*/

#ifndef _CPUID_H
#define _CPUID_H

/* This is a C library.  If used in a C++ program, don't let the C++ */
/* compiler mangle names. */

#ifdef __cplusplus
extern "C" {
#endif

/* Include common definitions */

#include "gwcommon.h"

/* Routines that the user should call */

void guessCpuType (void);
void guessCpuSpeed (void);

/* Routines to access the high resolution timer */

int isHighResTimerAvailable (void);
double getHighResTimer (void);
double getHighResTimerFrequency (void);

/* Global variables describing the CPU we are running on */

extern char CPU_BRAND[49];		/* Text description of CPU */
extern double CPU_SPEED;		/* Actual CPU Speed in MHz */
#define CPU_RDTSC	0x0001
#define CPU_CMOV	0x0002
#define CPU_PREFETCH	0x0004
#define CPU_SSE		0x0008
#define CPU_SSE2	0x0010
#define CPU_MMX		0x0020
#define CPU_3DNOW	0x0040
#define CPU_SSE3	0x0080
#define CPU_SSSE3	0x0100		/* Supplemental SSE3 */
#define CPU_SSE41	0x0200
#define CPU_SSE42	0x0400
extern unsigned int CPU_FLAGS;		/* Cpu capabilities */
extern unsigned int CPU_HYPERTHREADS;	/* Number of virtual processors */
					/* that each CPU core supports */
extern int CPU_L1_CACHE_SIZE;		/* In KB */
extern int CPU_L2_CACHE_SIZE;		/* In KB */
extern int CPU_L3_CACHE_SIZE;		/* In KB */
extern int CPU_L1_CACHE_LINE_SIZE;
extern int CPU_L2_CACHE_LINE_SIZE;
extern int CPU_L3_CACHE_LINE_SIZE;
extern int CPU_L1_DATA_TLBS;
extern int CPU_L2_DATA_TLBS;
extern int CPU_L3_DATA_TLBS;
extern int CPU_L1_SET_ASSOCIATIVE;
extern int CPU_L2_SET_ASSOCIATIVE;
extern int CPU_L3_SET_ASSOCIATIVE;

extern unsigned int CPU_SIGNATURE;	/* Vendor-specific family number, */
					/* model number, stepping ID, etc. */

/* Assembly language structures and routines */

struct cpuid_data {
	uint32_t EAX;		/* For communicating with asm routines */
	uint32_t EBX;
	uint32_t ECX;
	uint32_t EDX;
};

unsigned long ecpuidsupport (void);
void ecpuid (struct cpuid_data *);

/* Cleaner access to cpuid assembly code */

#define isCpuidSupported()	ecpuidsupport ()
#define Cpuid(a,s)		{ (s)->EAX=a; ecpuid (s); }

/* Routine used to time code chunks */
void erdtsc (uint32_t *hi, uint32_t *lo);
#define rdtsc(hi,lo)		erdtsc(hi,lo)

/* Init the x87 FPU, probably not needed any longer */
void fpu_init (void);

#ifdef __cplusplus
}
#endif

#endif
