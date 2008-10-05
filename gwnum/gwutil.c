/*----------------------------------------------------------------------
| gwutil.c
|
| This file contains various utility routines that may be used by gwnum
| routines, prime95, or PRP.
| 
|  Copyright 2004-2007 Mersenne Research, Inc.  All rights reserved.
+---------------------------------------------------------------------*/

/* Include files */

#include <stdlib.h>
#ifndef __APPLE__
#include <malloc.h>
#endif
#include "gwutil.h"


/* Aligned malloc routines.  MSVC 8 supports these in the C runtime library. */
/* Emulate these routines for other ports. */

void * aligned_offset_malloc (
	size_t	size,
	size_t	alignment,
	size_t	mod)
{
#ifdef _WIN64
	return (_aligned_offset_malloc (size, alignment, mod));
#else
	char	*p, *q;
	p = (char *) malloc (sizeof (void *) + size + alignment);
	if (p == NULL) return (NULL);
	q = (char *) (((long) p + sizeof (void *) + mod + alignment - 1) & ~(alignment - 1)) - mod;
	* (void **) ((char *) q - sizeof (void *)) = p;
	return (q);
#endif
}

void * aligned_malloc (
	size_t	size,
	size_t	alignment)
{
#ifdef _WIN64
	return (_aligned_malloc (size, alignment));
#else
	return (aligned_offset_malloc (size, alignment, 0));
#endif
}

void aligned_free (
	void	*ptr)
{
#ifdef _WIN64
	_aligned_free (ptr);
#else
	if (ptr == NULL) return;
	free (* (void **) ((char *) ptr - sizeof (void *)));
#endif
}
