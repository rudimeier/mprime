/*----------------------------------------------------------------------
| Copyright 1995-2008 Mersenne Research, Inc.  All rights reserved
| Author:  George Woltman
| Email: woltman@alum.mit.edu
|
| This file contains routines to determine the CPU type and speed.
+---------------------------------------------------------------------*/

/* Include files */

#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include <string.h>
#if defined (__linux__) || defined (__FreeBSD__) || defined (__EMX__) || defined (__APPLE__)
#include <sys/time.h>
#define _timeb	timeb
#define _ftime	ftime
#endif
#include <sys/timeb.h>
#include "cpuid.h"

/* Global variables describing the CPU we are running on */

char	CPU_BRAND[49] = "";
double	CPU_SPEED = 0.0;
unsigned int CPU_FLAGS = 0;
unsigned int CPU_HYPERTHREADS = 1;	/* Number of virtual processors */
					/* that each CPU core supports */
int	CPU_L1_CACHE_SIZE = -1;
int	CPU_L2_CACHE_SIZE = -1;
int	CPU_L3_CACHE_SIZE = -1;
int	CPU_L1_CACHE_LINE_SIZE = -1;
int	CPU_L2_CACHE_LINE_SIZE = -1;
int	CPU_L3_CACHE_LINE_SIZE = -1;
int	CPU_L1_DATA_TLBS = -1;
int	CPU_L2_DATA_TLBS = -1;
int	CPU_L3_DATA_TLBS = -1;
int	CPU_L1_SET_ASSOCIATIVE = -1;
int	CPU_L2_SET_ASSOCIATIVE = -1;
int	CPU_L3_SET_ASSOCIATIVE = -1;

unsigned int CPU_SIGNATURE = 0;		/* Vendor-specific family number, */
					/* model number, stepping ID, etc. */

/* Internal routines to see if CPU-specific instructions (RDTSC, CMOV */
/* SSE, SSE2) are supported.  CPUID could report them as supported yet */
/* the OS might not support them.   Use signals in non-MSVC compiles. */

#if !defined(_MSC_VER) && !defined(__WATCOMC__)

#include <setjmp.h>
#include <signal.h>
int	boom;
jmp_buf	env;
void sigboom_handler (int i)
{
	boom = TRUE;
	longjmp (env, 1);
}
int canExecInstruction (
	unsigned long cpu_flag)
{
	boom = FALSE;
	(void) signal (SIGILL, sigboom_handler);
	if (setjmp (env) == 0) {
		switch (cpu_flag) {
		case CPU_RDTSC:		/* RDTSC */
			__asm__ __volatile__ (".byte 0x0F\n .byte 0x31\n");
			break;
		case CPU_CMOV:		/* CMOV */
			__asm__ __volatile__ (".byte 0x0F\n .byte 0x42\n .byte 0xC0\n");
			break;
		case CPU_MMX:		/* PADDB */
			__asm__ __volatile__ (".byte 0x0F\n .byte 0xFC\n .byte 0xC0\n");
			break;
		case CPU_SSE:		/* ORPS */
			__asm__ __volatile__ (".byte 0x0F\n .byte 0x56\n .byte 0xC0\n");
			break;
		case CPU_SSE2:		/* ADDPD */
			__asm__ __volatile__ (".byte 0x66\n .byte 0x0F\n .byte 0x58\n .byte 0xC0\n");
			break;
		case CPU_PREFETCH:	/* PREFETCHT1 */
			__asm__ __volatile__ (".byte 0x0F\n .byte 0x18\n .byte 0x16\n");
			break;
		case CPU_3DNOW:		/* PREFETCHW */
			__asm__ __volatile__ (".byte 0x0F\n .byte 0x0D\n .byte 0x0E\n");
			break;
		}
	}
	(void) signal (SIGILL, SIG_DFL);
	fpu_init ();
	return (!boom);
}

#else

/* The MS 64-bit compiler does not allow inline assembly.  Fortunately, any */
/* CPU capable of running x86-64 bit code can execute these instructions. */

#ifdef X86_64

int canExecInstruction (
	unsigned long cpu_flag)
{
	return (TRUE);
}

#else

#include <excpt.h>
#ifdef __WATCOMC__
#define __asm	_asm
#define __emit	db
#endif

int canExecInstruction (
	unsigned long cpu_flag)
{
	int	succeeded;
	__try {
		switch (cpu_flag) {
		case CPU_RDTSC:		/* RDTSC */
			__asm __emit 0x0F
			__asm __emit 0x31
			break;
		case CPU_CMOV:		/* CMOV */
			__asm __emit 0x0F
			__asm __emit 0x42
			__asm __emit 0xC0
			break;
		case CPU_MMX:		/* PADDB */
			__asm __emit 0x0F
			__asm __emit 0xFC
			__asm __emit 0xC0
			break;
		case CPU_SSE:		/* ORPS */
			__asm __emit 0x0F
			__asm __emit 0x56
			__asm __emit 0xC0
			break;
		case CPU_SSE2:		/* ADDPD */
			__asm __emit 0x66
			__asm __emit 0x0F
			__asm __emit 0x58
			__asm __emit 0xC0
			break;
		case CPU_PREFETCH:	/* PREFETCHT1 */
			__asm __emit 0x0F
			__asm __emit 0x18
			__asm __emit 0x16
			break;
		case CPU_3DNOW:		/* PREFETCHW */
			__asm __emit 0x0F
			__asm __emit 0x0D
			__asm __emit 0x0E
			break;
		}
		succeeded = TRUE;
	}
	__except (EXCEPTION_EXECUTE_HANDLER) {
		succeeded = FALSE;
	}
	fpu_init ();
	return (succeeded);
}
#endif
#endif


/* Work with CPUID instruction to guess the cpu type and features */
/* See Intel's document AP-485 for using CPUID on Intel processors */
/* AMD and VIA have similar documents */

void guessCpuType (void)
{
	struct cpuid_data reg;
	unsigned long max_cpuid_value;
	unsigned long max_extended_cpuid_value;
	unsigned long extended_family, extended_model, type, family_code;
	unsigned long model_number, stepping_id, brand_index;
	char	vendor_id[13];
static	char *	BRAND_NAMES[] = {	/* From Intel Ap-485 */
			"",	/* brand_index = 0 */
			"Intel(R) Celeron(R) processor",
			"Intel(R) Pentium(R) III processor",
			"Intel(R) Pentium(R) III Xeon(TM) processor",
			"Intel(R) Pentium(R) III processor",
			"",	/* brand_index = 5 */
			"Mobile Intel(R) Pentium(R) III Processor - M",
			"Mobile Intel(R) Celeron(R) processor",
			"Intel(R) Pentium(R) 4 processor",
			"Intel(R) Pentium(R) 4 processor",
			"Intel(R) Celeron(R) processor",
			"Intel(R) Xeon(TM) processor",
			"Intel(R) Xeon(TM) Processor MP",
			"",	/* brand_index = 13 */
			"Mobile Intel(R) Pentium(R) 4 Processor - M",
			"Mobile Intel(R) Celeron(R) processor",
			"",	/* brand_index = 16 */
			"Mobile Intel(R) processor",
			"Mobile Intel(R) Celeron(R) M processor",
			"Mobile Intel(R) Celeron(R) processor",
			"Intel(R) Celeron(R) processor",
			"Mobile Intel(R) processor",
			"Intel(R) Pentium(R) M processor"
			"Mobile Intel(R) Celeron(R) processor"
	};
#define NUM_BRAND_NAMES	(sizeof (BRAND_NAMES) / sizeof (char *))

/* Set up default values for features we cannot determine with CPUID */

	CPU_BRAND[0] = 0;
	CPU_SPEED = 100.0;
	CPU_FLAGS = 0;
	CPU_HYPERTHREADS = 1;
	CPU_L1_CACHE_SIZE = -1;
	CPU_L2_CACHE_SIZE = -1;
	CPU_L3_CACHE_SIZE = -1;
	CPU_L1_CACHE_LINE_SIZE = -1;
	CPU_L2_CACHE_LINE_SIZE = -1;
	CPU_L3_CACHE_LINE_SIZE = -1;
	CPU_L1_DATA_TLBS = -1;
	CPU_L2_DATA_TLBS = -1;
	CPU_L3_DATA_TLBS = -1;
	CPU_SIGNATURE = 0;

/* If CPUID instruction is not supported, assume we have a 486 (not all */
/* 486 chips supported CPUID.  The CPU might be something else, but that */
/* isn't particularly important. */

	if (! isCpuidSupported ()) {
		strcpy (CPU_BRAND, "CPUID not supported - 486 CPU assumed");
		return;
	}

/* Call CPUID with 0 argument.  It returns how highest argument CPUID */
/* can accept as well as the vendor string */

	Cpuid (0, &reg);
	max_cpuid_value = reg.EAX;
	memcpy (vendor_id, &reg.EBX, 4);
	memcpy (vendor_id+4, &reg.EDX, 4);
	memcpy (vendor_id+8, &reg.ECX, 4);
	vendor_id[12] = 0;

/* So far all vendors have adopted Intel's definition of CPUID with 1 as an */
/* argument.  Let's assume future vendors will do the same.  CPUID returns */
/* the processor family, stepping, etc.  It also returns the feature flags. */

	if (max_cpuid_value >= 1) {
		Cpuid (1, &reg);
		CPU_SIGNATURE = reg.EAX & 0x0FFF3FFF;
		extended_family = (reg.EAX >> 20) & 0xFF;
		extended_model = (reg.EAX >> 16) & 0xF;
		type = (reg.EAX >> 12) & 0x3;
		family_code = (reg.EAX >> 8) & 0xF;
		model_number = (reg.EAX >> 4) & 0xF;
		stepping_id = reg.EAX & 0xF;
		brand_index = reg.EBX & 0xFF;
		if ((reg.EDX >> 4) & 0x1 && canExecInstruction (CPU_RDTSC))
			CPU_FLAGS |= CPU_RDTSC;
		if ((reg.EDX >> 15) & 0x1 && canExecInstruction (CPU_CMOV))
			CPU_FLAGS |= CPU_CMOV;
		if ((reg.EDX >> 23) & 0x1 && canExecInstruction (CPU_MMX))
			CPU_FLAGS |= CPU_MMX;
		if ((reg.EDX >> 25) & 0x1 && canExecInstruction (CPU_PREFETCH))
			CPU_FLAGS |= CPU_PREFETCH;
		if ((reg.EDX >> 25) & 0x1 && canExecInstruction (CPU_SSE))
			CPU_FLAGS |= CPU_SSE;
		if ((reg.EDX >> 26) & 0x1 && canExecInstruction (CPU_SSE2))
			CPU_FLAGS |= CPU_SSE2;
		if ((reg.ECX >> 0) & 0x1)
			CPU_FLAGS |= CPU_SSE3;
		if ((reg.ECX >> 9) & 0x1)
			CPU_FLAGS |= CPU_SSSE3;
		if ((reg.ECX >> 19) & 0x1)
			CPU_FLAGS |= CPU_SSE41;
		if ((reg.ECX >> 20) & 0x1)
			CPU_FLAGS |= CPU_SSE42;
	}

/* Call CPUID with 0x80000000 argument.  It tells us how many extended CPU */
/* functions are supported. */

	Cpuid (0x80000000, &reg);
	max_extended_cpuid_value = reg.EAX;

/* Although not guaranteed, all vendors have standardized on putting the */
/* brand string (if supported) at cpuid calls 0x8000002, 0x80000003, and */
/* 0x80000004.  We'll assume future vendors will do the same. */

	if (max_extended_cpuid_value >= 0x80000004) {
		Cpuid (0x80000002, &reg);
		memcpy (CPU_BRAND, &reg.EAX, 4);
		memcpy (CPU_BRAND+4, &reg.EBX, 4);
		memcpy (CPU_BRAND+8, &reg.ECX, 4);
		memcpy (CPU_BRAND+12, &reg.EDX, 4);
		Cpuid (0x80000003, &reg);
		memcpy (CPU_BRAND+16, &reg.EAX, 4);
		memcpy (CPU_BRAND+20, &reg.EBX, 4);
		memcpy (CPU_BRAND+24, &reg.ECX, 4);
		memcpy (CPU_BRAND+28, &reg.EDX, 4);
		Cpuid (0x80000004, &reg);
		memcpy (CPU_BRAND+32, &reg.EAX, 4);
		memcpy (CPU_BRAND+36, &reg.EBX, 4);
		memcpy (CPU_BRAND+40, &reg.ECX, 4);
		memcpy (CPU_BRAND+44, &reg.EDX, 4);
		CPU_BRAND[48] = 0;
		while (CPU_BRAND[0] == ' ') strcpy (CPU_BRAND, CPU_BRAND+1);
	}

/*-------------------------------------------------------------------+
| Check for INTEL vendor string.  Perform INTEL-specific operations. |
+-------------------------------------------------------------------*/

	if (strcmp ((const char *) vendor_id, "GenuineIntel") == 0) {

/* Try to determine if hyperthreading is supported.  I think this code */
/* only tells us if the hardware supports hyperthreading.  If the feature */
/* is turned off in the BIOS, we don't detect this. */

		if (max_cpuid_value >= 1) {
			Cpuid (1, &reg);
			if ((reg.EDX >> 28) & 0x1) {
				CPU_HYPERTHREADS = (reg.EBX >> 16) & 0xFF;
				if (max_cpuid_value >= 4) {
					reg.ECX = 0;
					Cpuid (4, &reg);
					CPU_HYPERTHREADS /= (reg.EAX >> 26) + 1;
				}
				if (CPU_HYPERTHREADS < 1) CPU_HYPERTHREADS = 1;
			}
		}

/* Call CPUID with 2 argument.  It returns the cache size and structure */
/* in a series of 8-bit descriptors */

		if (max_cpuid_value >= 2) {
			Cpuid (2, &reg);
			if ((reg.EAX & 0xFF) > 0) {
				unsigned int descriptors[15];
				int i, count;
				count = 0;
				if (! (reg.EAX & 0x80000000)) {
					descriptors[count++] = (reg.EAX >> 24) & 0xFF;
					descriptors[count++] = (reg.EAX >> 16) & 0xFF;
					descriptors[count++] = (reg.EAX >> 8) & 0xFF;
				}
				if (! (reg.EBX & 0x80000000)) {
					descriptors[count++] = (reg.EBX >> 24) & 0xFF;
					descriptors[count++] = (reg.EBX >> 16) & 0xFF;
					descriptors[count++] = (reg.EBX >> 8) & 0xFF;
					descriptors[count++] = reg.EBX & 0xFF;
				}
				if (! (reg.ECX & 0x80000000)) {
					descriptors[count++] = (reg.ECX >> 24) & 0xFF;
					descriptors[count++] = (reg.ECX >> 16) & 0xFF;
					descriptors[count++] = (reg.ECX >> 8) & 0xFF;
					descriptors[count++] = reg.ECX & 0xFF;
				}
				if (! (reg.EDX & 0x80000000)) {
					descriptors[count++] = (reg.EDX >> 24) & 0xFF;
					descriptors[count++] = (reg.EDX >> 16) & 0xFF;
					descriptors[count++] = (reg.EDX >> 8) & 0xFF;
					descriptors[count++] = reg.EDX & 0xFF;
				}
				for (i = 0; i < count; i++) {
					switch (descriptors[i]) {
					case 0x03:
						CPU_L2_DATA_TLBS = 64;
						break;
					case 0x0A:
						CPU_L1_CACHE_SIZE = 8;
						CPU_L1_CACHE_LINE_SIZE = 32;
						CPU_L1_SET_ASSOCIATIVE = 2;
						break;
					case 0x0C:
						CPU_L1_CACHE_SIZE = 16;
						CPU_L1_CACHE_LINE_SIZE = 32;
						CPU_L1_SET_ASSOCIATIVE = 4;
						break;
					case 0x22:
						CPU_L3_CACHE_SIZE = 512;
						CPU_L3_CACHE_LINE_SIZE = 64;
						CPU_L3_SET_ASSOCIATIVE = 4;
						break;
					case 0x23:
						CPU_L3_CACHE_SIZE = 1024;
						CPU_L3_CACHE_LINE_SIZE = 64;
						CPU_L3_SET_ASSOCIATIVE = 8;
						break;
					case 0x25:
						CPU_L3_CACHE_SIZE = 2048;
						CPU_L3_CACHE_LINE_SIZE = 64;
						CPU_L3_SET_ASSOCIATIVE = 8;
						break;
					case 0x29:
						CPU_L3_CACHE_SIZE = 4096;
						CPU_L3_CACHE_LINE_SIZE = 64;
						CPU_L3_SET_ASSOCIATIVE = 8;
						break;
					case 0x2C:
						CPU_L1_CACHE_SIZE = 32;
						CPU_L1_CACHE_LINE_SIZE = 64;
						CPU_L1_SET_ASSOCIATIVE = 8;
						break;
					case 0x39:
						CPU_L2_CACHE_SIZE = 128;
						CPU_L2_CACHE_LINE_SIZE = 128;
						CPU_L2_SET_ASSOCIATIVE = 4;
						break;
					case 0x3A:
						CPU_L2_CACHE_SIZE = 192;
						CPU_L2_CACHE_LINE_SIZE = 128;
						CPU_L2_SET_ASSOCIATIVE = 6;
						break;
					case 0x3B:
						CPU_L2_CACHE_SIZE = 128;
						CPU_L2_CACHE_LINE_SIZE = 128;
						CPU_L2_SET_ASSOCIATIVE = 2;
						break;
					case 0x3C:
						CPU_L2_CACHE_SIZE = 256;
						CPU_L2_CACHE_LINE_SIZE = 128;
						CPU_L2_SET_ASSOCIATIVE = 4;
						break;
					case 0x3D:
						CPU_L2_CACHE_SIZE = 384;
						CPU_L2_CACHE_LINE_SIZE = 128;
						CPU_L2_SET_ASSOCIATIVE = 6;
						break;
					case 0x3E:
						CPU_L2_CACHE_SIZE = 512;
						CPU_L2_CACHE_LINE_SIZE = 128;
						CPU_L2_SET_ASSOCIATIVE = 6;
						break;
					case 0x40:
						if (family_code == 15) {
							/* no L3 cache */
						} else {
							CPU_L2_CACHE_SIZE = 0;
						}
						break;
					case 0x41:
						CPU_L2_CACHE_SIZE = 128;
						CPU_L2_CACHE_LINE_SIZE = 32;
						CPU_L2_SET_ASSOCIATIVE = 4;
						break;
					case 0x42:
						CPU_L2_CACHE_SIZE = 256;
						CPU_L2_CACHE_LINE_SIZE = 32;
						CPU_L2_SET_ASSOCIATIVE = 4;
						break;
					case 0x43:
						CPU_L2_CACHE_SIZE = 512;
						CPU_L2_CACHE_LINE_SIZE = 32;
						CPU_L2_SET_ASSOCIATIVE = 4;
						break;
					case 0x44:
						CPU_L2_CACHE_SIZE = 1024;
						CPU_L2_CACHE_LINE_SIZE = 32;
						CPU_L2_SET_ASSOCIATIVE = 4;
						break;
					case 0x45:
						CPU_L2_CACHE_SIZE = 2048;
						CPU_L2_CACHE_LINE_SIZE = 32;
						CPU_L2_SET_ASSOCIATIVE = 4;
						break;
					case 0x46:
						CPU_L3_CACHE_SIZE = 4096;
						CPU_L3_CACHE_LINE_SIZE = 64;
						CPU_L3_SET_ASSOCIATIVE = 4;
						break;
					case 0x47:
						CPU_L3_CACHE_SIZE = 8192;
						CPU_L3_CACHE_LINE_SIZE = 64;
						CPU_L3_SET_ASSOCIATIVE = 8;
						break;
					case 0x49:
						if (family_code == 0x0F && model_number == 0x06) {
							CPU_L3_CACHE_SIZE = 4096;
							CPU_L3_CACHE_LINE_SIZE = 64;
							CPU_L3_SET_ASSOCIATIVE = 16;
						} else {
							CPU_L2_CACHE_SIZE = 4096;
							CPU_L2_CACHE_LINE_SIZE = 64;
							CPU_L2_SET_ASSOCIATIVE = 16;
						}
						break;
					case 0x4A:
						CPU_L3_CACHE_SIZE = 6144;
						CPU_L3_CACHE_LINE_SIZE = 64;
						CPU_L3_SET_ASSOCIATIVE = 12;
						break;
					case 0x4B:
						CPU_L3_CACHE_SIZE = 8192;
						CPU_L3_CACHE_LINE_SIZE = 64;
						CPU_L3_SET_ASSOCIATIVE = 16;
						break;
					case 0x4C:
						CPU_L3_CACHE_SIZE = 12288;
						CPU_L3_CACHE_LINE_SIZE = 64;
						CPU_L3_SET_ASSOCIATIVE = 12;
						break;
					case 0x4D:
						CPU_L3_CACHE_SIZE = 16384;
						CPU_L3_CACHE_LINE_SIZE = 64;
						CPU_L3_SET_ASSOCIATIVE = 16;
						break;
					case 0x4E:
						CPU_L2_CACHE_SIZE = 6144;
						CPU_L2_CACHE_LINE_SIZE = 64;
						CPU_L2_SET_ASSOCIATIVE = 24;
						break;
					case 0x5B:
						CPU_L2_DATA_TLBS = 64;
						break;
					case 0x5C:
					case 0xB3:
						CPU_L2_DATA_TLBS = 128;
						break;
					case 0x5D:
					case 0xB4:
						CPU_L2_DATA_TLBS = 256;
						break;
					case 0x60:
						CPU_L1_CACHE_SIZE = 16;
						CPU_L1_CACHE_LINE_SIZE = 64;
						CPU_L1_SET_ASSOCIATIVE = 8;
						break;
					case 0x66:
						CPU_L1_CACHE_SIZE = 8;
						CPU_L1_CACHE_LINE_SIZE = 64;
						CPU_L1_SET_ASSOCIATIVE = 4;
						break;
					case 0x67:
						CPU_L1_CACHE_SIZE = 16;
						CPU_L1_CACHE_LINE_SIZE = 64;
						CPU_L1_SET_ASSOCIATIVE = 4;
						break;
					case 0x68:
						CPU_L1_CACHE_SIZE = 32;
						CPU_L1_CACHE_LINE_SIZE = 64;
						CPU_L1_SET_ASSOCIATIVE = 4;
						break;
					case 0x78:
						CPU_L2_CACHE_SIZE = 1024;
						CPU_L2_CACHE_LINE_SIZE = 64;
						CPU_L2_SET_ASSOCIATIVE = 4;
						break;
					case 0x79:
						CPU_L2_CACHE_SIZE = 128;
						CPU_L2_CACHE_LINE_SIZE = 128;
						CPU_L2_SET_ASSOCIATIVE = 8;
						break;
					case 0x7A:
						CPU_L2_CACHE_SIZE = 256;
						CPU_L2_CACHE_LINE_SIZE = 128;
						CPU_L2_SET_ASSOCIATIVE = 8;
						break;
					case 0x7B:
						CPU_L2_CACHE_SIZE = 512;
						CPU_L2_CACHE_LINE_SIZE = 128;
						CPU_L2_SET_ASSOCIATIVE = 8;
						break;
					case 0x7C:
						CPU_L2_CACHE_SIZE = 1024;
						CPU_L2_CACHE_LINE_SIZE = 128;
						CPU_L2_SET_ASSOCIATIVE = 8;
						break;
					case 0x7D:
						CPU_L2_CACHE_SIZE = 2048;
						CPU_L2_CACHE_LINE_SIZE = 64;
						CPU_L2_SET_ASSOCIATIVE = 8;
						break;
					case 0x7F:
						CPU_L2_CACHE_SIZE = 512;
						CPU_L2_CACHE_LINE_SIZE = 64;
						CPU_L2_SET_ASSOCIATIVE = 2;
						break;
					case 0x82:
						CPU_L2_CACHE_SIZE = 256;
						CPU_L2_CACHE_LINE_SIZE = 32;
						CPU_L2_SET_ASSOCIATIVE = 8;
						break;
					case 0x83:
						CPU_L2_CACHE_SIZE = 512;
						CPU_L2_CACHE_LINE_SIZE = 32;
						CPU_L2_SET_ASSOCIATIVE = 8;
						break;
					case 0x84:
						CPU_L2_CACHE_SIZE = 1024;
						CPU_L2_CACHE_LINE_SIZE = 32;
						CPU_L2_SET_ASSOCIATIVE = 8;
						break;
					case 0x85:
						CPU_L2_CACHE_SIZE = 2048;
						CPU_L2_CACHE_LINE_SIZE = 32;
						CPU_L2_SET_ASSOCIATIVE = 8;
						break;
					case 0x86:
						CPU_L2_CACHE_SIZE = 512;
						CPU_L2_CACHE_LINE_SIZE = 64;
						CPU_L2_SET_ASSOCIATIVE = 4;
						break;
					case 0x87:
						CPU_L2_CACHE_SIZE = 1024;
						CPU_L2_CACHE_LINE_SIZE = 64;
						CPU_L2_SET_ASSOCIATIVE = 8;
						break;
					}
				}
			}
		}

/* If we haven't figured out the L2 cache size, use 0x80000006 to deduce */
/* the cache size. */

		if (CPU_L2_CACHE_SIZE == -1 &&
		    max_extended_cpuid_value >= 0x80000006) {
			Cpuid (0x80000006, &reg);
			CPU_L2_CACHE_LINE_SIZE = reg.ECX & 0xFF;
			if ((reg.ECX >> 12 & 0xF) == 2)
				CPU_L2_SET_ASSOCIATIVE = 2;
			if ((reg.ECX >> 12 & 0xF) == 4)
				CPU_L2_SET_ASSOCIATIVE = 4;
			if ((reg.ECX >> 12 & 0xF) == 6)
				CPU_L2_SET_ASSOCIATIVE = 8;
			if ((reg.ECX >> 12 & 0xF) == 8)
				CPU_L2_SET_ASSOCIATIVE = 16;
			CPU_L2_CACHE_SIZE = (reg.ECX >> 16);
		}

/* If we still haven't figured out the L2 or L3 cache size, use 0x00000004 to deduce */
/* the cache size. */

		if ((CPU_L2_CACHE_SIZE == -1 || CPU_L3_CACHE_SIZE == -1) &&
		    max_cpuid_value >= 4) {
			int	i;
			for (i = 0; i < 10; i++) {
				int	cache_type, cache_level, cores, threads;
				int	associativity, partitions, line_size, sets;
				int	prefetch_stride;

				reg.ECX = i;
				Cpuid (4, &reg);
				cache_type = reg.EAX & 0x1F;
				if (cache_type == 0) break;
				if (cache_type == 2) continue;  /* Ignore instruction caches */
				cores = (reg.EAX >> 26) + 1;
				threads = ((reg.EAX >> 14) & 0xFFF) + 1;
				cache_level = (reg.EAX >> 5) & 0x7;
				associativity = (reg.EBX >> 22) + 1;
				partitions = ((reg.EBX >> 12) & 0x3FF) + 1;
				line_size = (reg.EBX & 0xFFF);
				sets = reg.ECX + 1;
				prefetch_stride = (reg.EDX & 0x3FF);
				if (prefetch_stride == 0) prefetch_stride = 64;
				if (CPU_L2_CACHE_SIZE == -1 && cache_level == 2) {
					CPU_L2_CACHE_SIZE = (associativity * partitions * line_size * sets) >> 10;
					CPU_L2_SET_ASSOCIATIVE = associativity;
					CPU_L2_CACHE_LINE_SIZE = prefetch_stride;
				}
				if (CPU_L3_CACHE_SIZE == -1 && cache_level == 3) {
					CPU_L3_CACHE_SIZE = (associativity * partitions * line_size * sets) >> 10;
					CPU_L3_SET_ASSOCIATIVE = associativity;
					CPU_L3_CACHE_LINE_SIZE = prefetch_stride;
				}
			}
		}

/* If we haven't figured out the brand string, create one based on the */
/* sample code in Intel's AP-485 document. */

		if (CPU_BRAND[0] == 0) {

			if (family_code == 4)
				strcpy (CPU_BRAND, "Intel 486 processor");

			else if (family_code == 5) {
				if (type == 0 && model_number <= 2)
					strcpy (CPU_BRAND, "Intel Pentium processor");
				else if (type == 1 && model_number <= 3)
					strcpy (CPU_BRAND, "Intel Pentium OverDrive processor");
				else if (type == 0 && model_number >= 4)
					strcpy (CPU_BRAND, "Intel Pentium MMX processor");
				else if (type == 1 && model_number >= 4)
					strcpy (CPU_BRAND, "Intel Pentium MMX OverDrive processor");
				else
					strcpy (CPU_BRAND, "Intel Pentium processor");
			}

			else if (family_code == 6 && model_number == 1)
				strcpy (CPU_BRAND, "Intel Pentium Pro processor");

			else if (family_code == 6 && model_number == 3) {
				if (type == 0)
					strcpy (CPU_BRAND, "Intel Pentium II processor");
				else
					strcpy (CPU_BRAND, "Intel Pentium II OverDrive processor");
			}

			else if (family_code == 6 && model_number == 5) {
				if (CPU_L2_CACHE_SIZE == 0)
					strcpy (CPU_BRAND, "Intel Celeron processor");
				else if (CPU_L2_CACHE_SIZE >= 1024)
					strcpy (CPU_BRAND, "Intel Pentium II Xeon processor");
				else
					strcpy (CPU_BRAND, "Intel Pentium II or Pentium II Xeon processor");
			}

			else if (family_code == 6 && model_number == 6)
				strcpy (CPU_BRAND, "Intel Celeron processor");

			else if (family_code == 6 && model_number == 7) {
				if (CPU_L2_CACHE_SIZE >= 1024)
					strcpy (CPU_BRAND, "Intel Pentium III Xeon processor");
				else
					strcpy (CPU_BRAND, "Intel Pentium III or Pentium III Xeon processor");
			}

			else if (family_code == 6 && model_number == 0xB &&
				 stepping_id == 1 && brand_index == 0x3)
				strcpy (CPU_BRAND, "Intel(R) Celeron(R) processor");

			else if (family_code == 15 && model_number == 1 &&
				 stepping_id == 3 && brand_index == 0xB)
				strcpy (CPU_BRAND, "Intel(R) Xeon(TM) processor MP");

			else if (family_code == 15 && model_number == 1 &&
				 stepping_id == 3 && brand_index == 0xE)
				strcpy (CPU_BRAND, "Intel(R) Xeon(TM) processor");

			else if (brand_index != 0 &&
				 brand_index < NUM_BRAND_NAMES)
				strcpy (CPU_BRAND, BRAND_NAMES[brand_index]);

/* If we've failed to figure out the brand string, create a default. */

			else if (CPU_BRAND[0] == 0)
				strcpy (CPU_BRAND, "Unknown Intel CPU");
		}
	}

/*---------------------------------------------------------------+
| Check for AMD vendor string.  Perform AMD-specific operations. |
+---------------------------------------------------------------*/

	else if (strcmp ((const char *) vendor_id, "AuthenticAMD") == 0) {
		unsigned long extended_feature_bits = 0;
		unsigned long advanced_power_mgmt = 0;

		if (max_extended_cpuid_value >= 0x80000001) {
			Cpuid (0x80000001, &reg);
			extended_feature_bits = reg.EDX;
		}

		if (max_extended_cpuid_value >= 0x80000007) {
			Cpuid (0x80000007, &reg);
			advanced_power_mgmt = reg.EDX;
		}

/* Deduce the cpu type given the family, model, stepping, etc.  If we */
/* haven't figured out the brand string, create one based on the cpu type. */

/* Information from AMD Processor Recognition Application Note, */
/* Publication # 20734 Revision: 3.07 February 2004 */
/* Chap.3 Table 4: Summary of Processor Signatures for AMD Processors */

		if (family_code == 4) {
			strcpy (CPU_BRAND, "AMD Am486 or Am5x86 processor");
		}
		if ((CPU_BRAND[0] == 0) || (strstr (CPU_BRAND, "Unknown"))) {
			if (family_code == 5) {
				if (model_number <= 3) {
					strcpy (CPU_BRAND, "AMD K5 processor");
				}
				else if ((model_number >= 6) && (model_number <= 7)) {
					strcpy (CPU_BRAND, "AMD K6 processor");
				}
				else if (model_number == 8) {
					strcpy (CPU_BRAND, "AMD K6-2 processor");
				}
				else if (model_number == 9) {
					strcpy (CPU_BRAND, "AMD K6-III processor");
				}
			}
			if (family_code == 6) {
				int mobileCPU = (advanced_power_mgmt & 7) == 7;
				int mp = 0 != (extended_feature_bits & (1 << 19));
				char *szCPU = "";

				switch (model_number)
				{
					case 1:
					case 2:
					case 4:
						szCPU = "AMD Athlon processor";
						break;
					case 3:
						szCPU = "AMD Duron processor";
						break;
					case 6:
						if (mobileCPU) {
							szCPU = "Mobile AMD Athlon 4 processor";
						}
						else {
							szCPU = (mp) ? "AMD Athlon MP processor" : "AMD Athlon XP processor";
						}
						break;
					case 7:
						szCPU = (mobileCPU) ? "Mobile AMD Duron processor" : "AMD Duron processor";
						break;
					case 8:
						szCPU = (mp) ? "AMD Athlon MP processor" : "AMD Athlon XP processor";
						break;
					case 10:
						if (mobileCPU) {
							szCPU = "Mobile AMD Athlon XP-M processor";
						}
						else {
							szCPU = (mp) ? "AMD Athlon MP processor" : "AMD Athlon XP processor";
						}
				}
				strcpy (CPU_BRAND, szCPU);
			}
			if (family_code == 15) {
				if (model_number == 4) {
					strcpy (CPU_BRAND, "AMD Athlon 64 processor");
				}
				else if (model_number == 5) {
					strcpy (CPU_BRAND, "AMD Opteron processor");
				}
			}
		}

/* Early Athlon CPUs support the SSE prefetch instructions even though */
/* they do not support the full SSE instruction set.  I think testing for */
/* the AMD MMX extensions capability will detect this case. */

		if (max_extended_cpuid_value >= 0x80000001 &&
		    ! (CPU_FLAGS & CPU_PREFETCH)) {
			Cpuid (0x80000001, &reg);
			if ((reg.EDX >> 22) & 0x1 &&
			    canExecInstruction (CPU_PREFETCH))
				CPU_FLAGS |= CPU_PREFETCH;
		}

/* Check for support of 3DNow! instructions.  The prefetchw instruction */
/* from the 3DNow! instruction set is used by the assembly code. */

		if (max_extended_cpuid_value >= 0x80000001 &&
		    ! (CPU_FLAGS & CPU_3DNOW)) {
			Cpuid (0x80000001, &reg);
			if ((reg.EDX >> 31) & 0x1 &&
			    canExecInstruction (CPU_3DNOW))
				CPU_FLAGS |= CPU_3DNOW;
		}

/* Get the L1 cache size and number of data TLBs */

		if (max_extended_cpuid_value >= 0x80000005) {
			Cpuid (0x80000005, &reg);
			CPU_L1_DATA_TLBS = (reg.EBX >> 16) & 0xFF;
			CPU_L1_CACHE_SIZE = (reg.ECX >> 24) & 0xFF;
			CPU_L1_CACHE_LINE_SIZE = reg.ECX & 0xFF;
		}

/* Get the L2 cache size */

		if (max_extended_cpuid_value >= 0x80000006) {
			Cpuid (0x80000006, &reg);
			CPU_L2_DATA_TLBS = (reg.EBX >> 16) & 0xFFF;
			CPU_L2_CACHE_SIZE = (reg.ECX >> 16) & 0xFFFF;
			if (CPU_L2_CACHE_SIZE == 1) /* Workaround Duron bug */
				CPU_L2_CACHE_SIZE = 64;
			CPU_L2_CACHE_LINE_SIZE = reg.ECX & 0xFF;
		}

/* If we haven't figured out the brand string, create a default one */

		if (CPU_BRAND[0] == 0)
			strcpy (CPU_BRAND, "Unknown AMD CPU");
	}

/*---------------------------------------------------------------+
| Check for VIA vendor string.  Perform VIA-specific operations. |
+---------------------------------------------------------------*/

	else if (strcmp ((const char *) vendor_id, "CentaurHauls") == 0) {

/* Get the L1 cache size and number of data TLBs */

		if (max_extended_cpuid_value >= 0x80000005) {
			Cpuid (0x80000005, &reg);
			CPU_L2_DATA_TLBS = (reg.EBX >> 16) & 0xFF;
			CPU_L1_CACHE_SIZE = (reg.ECX >> 24) & 0xFF;
			CPU_L1_CACHE_LINE_SIZE = reg.ECX & 0xFF;
		}

/* Get the L2 cache size */

		if (max_extended_cpuid_value >= 0x80000006) {
			Cpuid (0x80000006, &reg);
			CPU_L2_CACHE_SIZE = (reg.ECX >> 24) & 0xFF;
			CPU_L2_CACHE_LINE_SIZE = reg.ECX & 0xFF;
		}

/* If we haven't figured out the brand string, create a default one */

		if (CPU_BRAND[0] == 0)
			strcpy (CPU_BRAND, "Unknown VIA/CYRIX CPU");
	}

/*--------------------------------------------------------+
| An unknown CPU vendor.  Fill in defaults as best we can |
+--------------------------------------------------------*/

	else {
		if (CPU_BRAND[0] == 0) {
			strcpy (CPU_BRAND, "Unrecognized CPU vendor: ");
			strcat (CPU_BRAND, vendor_id);
		}
	}
}

void guessCpuSpeed (void)
{

/* If RDTSC is not supported, then measuring the CPU speed is real hard */
/* so we just assume the CPU speed is 100 MHz.  This isn't a big deal */
/* since CPUs not supporting RDTSC aren't powerful enough to run prime95. */

	if (! (CPU_FLAGS & CPU_RDTSC)) {
		CPU_SPEED = 100.0;
		return;
	}

/* If this machine supports a high resolution counter, use that for timing */

	if (isHighResTimerAvailable ()) {
		uint32_t start_hi, start_lo, end_hi, end_lo;
		double	frequency, temp, start_time, end_time;
		unsigned long iterations;
		double	speed1, speed2, speed3, avg_speed;
		int	tries;

/* Compute the number of high resolution ticks in one millisecond */
/* This should give us good accuracy while hopefully avoiding time slices */

		frequency = getHighResTimerFrequency ();
		iterations = (unsigned long) (frequency / 1000.0);

/* Do up to 20 iterations (idea lifted from Intel's code) until 3 straight */
/* speed calculations are within 1 MHz of each other.  This is good since */
/* outside forces can interfere with this calculation. */

		tries = 0; speed1 = 0.0; speed2 = 0.0;
		do {

/* Shuffle the last calculations, bump counter */

			speed3 = speed2;
			speed2 = speed1;
			tries++;

/* Loop waiting for high resolution timer to change */

			temp = getHighResTimer ();
			while ((start_time = getHighResTimer ()) == temp);
			rdtsc (&start_hi, &start_lo);

/* Now loop waiting for timer to tick off about 1 millisecond */

			temp = start_time + (double) iterations;
			while ((end_time = getHighResTimer ()) < temp);
			rdtsc (&end_hi, &end_lo);

/* Compute speed based on number of clocks in the time interval */

			speed1 = (end_hi * 4294967296.0 + end_lo -
				  start_hi * 4294967296.0 - start_lo) *
				 frequency /
				 (end_time - start_time) / 1000000.0;

/* Calculate average of last 3 speeds.  Loop if this average isn't */
/* very close to all of the last three speed calculations. */

			avg_speed = (speed1 + speed2 + speed3) / 3.0;
		} while (tries < 3 ||
		         (tries < 20 &&
		          (fabs (speed1 - avg_speed) > 1.0 ||
		           fabs (speed2 - avg_speed) > 1.0 ||
		           fabs (speed3 - avg_speed) > 1.0)));

/* Final result is average speed of last three calculations */

		CPU_SPEED = avg_speed;
	}

/* Otherwise use the low resolution timer to measure CPU speed */

	else {
		struct _timeb temp, start_time, end_time;
		uint32_t start_hi, start_lo, end_hi, end_lo;
		double	speed1, speed2, speed3, avg_speed, elapsed_time;
		int	tries;

/* Do up to 10 iterations until 3 straight speed calculations are within */
/* 1 MHz of each other. */

		tries = 0; speed1 = 0.0; speed2 = 0.0;
		do {

/* Shuffle the last calculations, bump counter */

			speed3 = speed2;
			speed2 = speed1;
			tries++;

/* Loop waiting for low resolution timer to change */

			_ftime (&temp);
			do
				_ftime (&start_time);
			while (temp.millitm == start_time.millitm);
			rdtsc (&start_hi, &start_lo);

/* Now loop waiting for timer to change again */

			do
				_ftime (&end_time);
			while (start_time.millitm == end_time.millitm);
			rdtsc (&end_hi, &end_lo);

/* Compute elapsed time.  Since most PCs have a low resolution clock */
/* that ticks every 18.20648193 seconds, then if elapsed time is close */
/* to 1 / 18.2 seconds, then assume the elapsed time is */
/* 1 / 18.206... = 0.054925493 seconds. */

			elapsed_time = (end_time.time - start_time.time) +
				  ((int) end_time.millitm -
				   (int) start_time.millitm) / 1000.0;
			if (elapsed_time >= 0.049 && elapsed_time <= 0.061)
				elapsed_time = 0.054925493;

/* Compute speed based on number of clocks in the time interval */

			speed1 = (end_hi * 4294967296.0 + end_lo -
				  start_hi * 4294967296.0 - start_lo) /
				 elapsed_time / 1000000.0;

/* Caclulate average of last 3 speeds.  Loop if this average isn't */
/* very close to all of the last three speed calculations. */

			avg_speed = (speed1 + speed2 + speed3) / 3.0;
		} while (tries < 3 ||
		         (tries < 10 &&
		          (fabs (speed1 - avg_speed) > 1.0 ||
		           fabs (speed2 - avg_speed) > 1.0 ||
		           fabs (speed3 - avg_speed) > 1.0)));

/* Final result is average speed of last three calculations */

		CPU_SPEED = avg_speed;
	}
}


/*--------------------------------------------------------------------------
| And now, the routines that access the high resolution performance counter.
+-------------------------------------------------------------------------*/

#ifdef _WIN32
#include <windows.h>
#endif

#ifdef __OS2__
#define INCL_DOSPROFILE
#include <os2.h>
#endif

int isHighResTimerAvailable (void)
{
#ifdef _WIN32
	LARGE_INTEGER large;
	return (QueryPerformanceCounter (&large));
#endif
#ifdef __OS2__
/* DosTmrQueryTime/DosTmrQueryFreq functions use the 8253/4 timer chip to */
/* obtain a timestamp.  DosTmrQueryTime returns the current count, and */
/* DosTmrQueryFreq returns the frequency at which the counter increments. */
/* This frequency is about 1.1MHz, which gives you a timer that's accurate */
/* to the microsecond. */
        return (TRUE);
#endif
#if defined (__linux__) || defined (__FreeBSD__) || defined (__APPLE__)
	struct timeval start, end;
	struct timezone tz;
	int	i;

/* Return true if gettimeofday is more accurate than 1/10 millisecond. */
/* Try 10 times to see if gettimeofday returns two values less than */
/* 100 microseconds apart. */

	for (i = 0; i < 10; i++) {
		gettimeofday (&start, &tz);
		for ( ; ; ) {
			gettimeofday (&end, &tz);
			if (start.tv_sec != end.tv_sec) break;
			if (start.tv_usec == end.tv_usec) continue;
			if (end.tv_usec - start.tv_usec < 100) return (TRUE);
		}
	}
	return (FALSE);
#endif
}

double getHighResTimer (void)
{
#ifdef _WIN32
	LARGE_INTEGER large;

	QueryPerformanceCounter (&large);
	return ((double) large.HighPart * 4294967296.0 +
		(double) /*(unsigned long)*/ large.LowPart);
#endif
#ifdef __OS2__
        unsigned long long qwTmrTime;
        DosTmrQueryTime((PQWORD)&qwTmrTime);
        return (qwTmrTime);
#endif
#if defined (__linux__) || defined (__FreeBSD__) || defined (__APPLE__)
	struct timeval x;
	struct timezone tz;

	gettimeofday (&x, &tz);
	return ((double) x.tv_sec * 1000000.0 + (double) x.tv_usec);
#endif
}

double getHighResTimerFrequency (void)
{
#ifdef _WIN32
	LARGE_INTEGER large;

	QueryPerformanceFrequency (&large);
	return ((double) large.HighPart * 4294967296.0 +
		(double) /*(unsigned long)*/ large.LowPart);
#endif
#ifdef __OS2__
	ULONG ulTmrFreq;
	DosTmrQueryFreq(&ulTmrFreq);
	return (ulTmrFreq);
#endif
#if defined (__linux__) || defined (__FreeBSD__) || defined (__APPLE__)
	return (1000000.0);
#endif
}
