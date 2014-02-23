/*----------------------------------------------------------------------
| common.h
|
| This file contains handy #defines that I use in all my projects
| 
|  Copyright 2005-2010 Mersenne Research, Inc.
|  All Rights Reserved.
+---------------------------------------------------------------------*/

#ifndef _COMMON_H
#define _COMMON_H

#ifdef __cplusplus
#define EXTERNC extern "C"
#else
#define EXTERNC
#endif

#define TRUE	1
#define FALSE	0

/* Define the world/group/owner read/write attributes for creating files */
/* I've always used 0666 in Unix (everyone gets R/W access), but MSVC 8 */
/* now refuses to work with that setting -- insisting on 0600 instead. */

#ifdef _WIN32
#define	CREATE_FILE_ACCESS	0600
#else
#define	CREATE_FILE_ACCESS	0666
#endif

/* Define the ASSERT macro I use while debugging */

#ifdef GDEBUG
#include <assert.h>
#define ASSERTG		assert
#else
#define ASSERTG(a)
#endif

/* Define a "safe" strcpy.  The official C runtime library says that overlapping */
/* buffers produce undefined results.  This safe strcpy allows overlapping */
/* buffers by using memmove instead. */

#define safe_strcpy(d,s)	memmove (d, s, strlen (s) + 1)
#ifdef GDEBUG
#undef strcpy
#define strcpy(d,s)	assert((d) >= ((s)+strlen(s)+1) || (s) >= (d)+strlen(s)+1), safe_strcpy(d,s)
#endif

#endif
