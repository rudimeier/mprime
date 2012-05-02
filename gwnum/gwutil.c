/*----------------------------------------------------------------------
| gwutil.c
|
| This file contains various utility routines that may be used by gwnum
| routines, prime95, or PRP.
| 
|  Copyright 2004-2009 Mersenne Research, Inc.  All rights reserved.
+---------------------------------------------------------------------*/

/* Include files */

#ifdef _WIN32
#include "windows.h"
#endif
#include <stdlib.h>
#if defined (__linux__) || defined (__HAIKU__)
#include <malloc.h>
#endif
#include "gwcommon.h"
#include "gwutil.h"


/* Aligned malloc routines.  MSVC 8 supports these in the C runtime library. */
/* Emulate these routines for other ports. */

void * aligned_offset_malloc (
	size_t	size,
	size_t	alignment,
	size_t	mod)
{
#ifdef _WIN64
	char	*p, *q;
	p = (char *) malloc (sizeof (void *) + size + alignment);
	if (p == NULL) return (NULL);
	q = (char *) (((uint64_t) p + sizeof (void *) + mod + alignment - 1) & ~(alignment - 1)) - mod;
	* (void **) ((char *) q - sizeof (void *)) = p;
	return (q);
// For compatibility with 32-bit versions of the program, I've elected to use
// my own implementation rather than the Microsoft routine below.
//	return (_aligned_offset_malloc (size, alignment, mod));
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
	return (aligned_offset_malloc (size, alignment, 0));
// For compatibility with 32-bit versions of the program, I've elected to use
// my own implementation rather than the Microsoft routine below.
//	return (_aligned_malloc (size, alignment));
#else
	return (aligned_offset_malloc (size, alignment, 0));
#endif
}

void aligned_free (
	void	*ptr)
{
#ifdef _WIN64
	if (ptr == NULL) return;
	free (* (void **) ((char *) ptr - sizeof (void *)));
// For compatibility with 32-bit versions of the program, I've elected to use
// my own implementation rather than the Microsoft routine below.
//	_aligned_free (ptr);
#else
	if (ptr == NULL) return;
	free (* (void **) ((char *) ptr - sizeof (void *)));
#endif
}

void * large_pages_malloc (
	size_t	size)
{

// Jean Penne reports that MSVC6 does not define MEM_LARGE_PAGES.
// Simple fix - we don't support large pages for older compilers.
// This is fine since we barely support large pages right now - they were
// added primarily for me to test if they might provide a performance benefit.

#if defined (_WIN32) && defined (MEM_LARGE_PAGES)
	static int first_call = 1;
	static size_t large_page_size = 0;
	static DWORD lasterr = 0;
	LPVOID p;

	if (first_call) {
		HANDLE hToken;
		LUID luid;
		TOKEN_PRIVILEGES tp;
		HINSTANCE  hDll;      
		int (*pGetLargePageMinimum)(void);

		first_call = 0;

		// Grant large page access
		OpenProcessToken (GetCurrentProcess(),
				  TOKEN_ADJUST_PRIVILEGES, &hToken);
		lasterr = GetLastError ();
		LookupPrivilegeValue (NULL, "SeLockMemoryPrivilege", &luid);
		lasterr = GetLastError ();

		tp.PrivilegeCount = 1;
		tp.Privileges[0].Luid = luid;
		tp.Privileges[0].Attributes = SE_PRIVILEGE_ENABLED;
		AdjustTokenPrivileges (hToken, FALSE, &tp,
				       sizeof (TOKEN_PRIVILEGES),
				       (PTOKEN_PRIVILEGES) NULL,
				       (PDWORD) NULL);
		lasterr = GetLastError ();

		// Dynamic link to get large page size
		// Call succeeds only on Windows Server 2003 SP1 or later
		hDll = LoadLibrary (TEXT ("kernel32.dll"));
		pGetLargePageMinimum = (int (*)(void)) GetProcAddress (hDll, "GetLargePageMinimum");
		if (pGetLargePageMinimum != NULL) 
			large_page_size = (*pGetLargePageMinimum)();
		FreeLibrary(hDll);
	}

	// If this system does not support large pages, return NULL
	if (large_page_size == 0) return (NULL);

	// Now allocate the memory
	p = VirtualAlloc (NULL, (size + large_page_size - 1) & ~(large_page_size - 1),
			  MEM_RESERVE | MEM_COMMIT | MEM_LARGE_PAGES, PAGE_READWRITE);
	lasterr = GetLastError ();
	return (p);
#else
	return (NULL);
#endif
}

void large_pages_free (
	void	*ptr)
{
#ifdef _WIN32
	VirtualFree (ptr, 0, MEM_RELEASE);
#endif
}
