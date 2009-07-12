/*----------------------------------------------------------------------
| This file contains various utility routines that may be used by gwnum
| routines, prime95, or PRP.
| 
|  Copyright 2004-2009 Mersenne Research, Inc.  All rights reserved.
+---------------------------------------------------------------------*/

#ifndef _GWUTIL_H
#define _GWUTIL_H

/* This is a C library.  If used in a C++ program, don't let the C++ */
/* compiler mangle names. */

#ifdef __cplusplus
extern "C" {
#endif

/* Align a pointer to the given boundary (boundary must be a power of 2) */

#ifdef _WIN64
#define align_ptr(p,n)	(void *) (((uint64_t)(p) + (n)-1) & ~((n)-1))
#else
#define align_ptr(p,n)	(void *) (((long)(p) + (n)-1) & ~((n)-1))
#endif

/* Aligned malloc routines.  MSVC 8 supports these in the C runtime library. */
/* Emulate these routines for other ports. */

void * aligned_offset_malloc (size_t size, size_t alignment, size_t mod);
void * aligned_malloc (size_t size, size_t alignment);
void  aligned_free (void *ptr);

/* Large page allocation routines */

void * large_pages_malloc (size_t size);
void large_pages_free (void *ptr);

#ifdef __cplusplus
}
#endif

#endif
