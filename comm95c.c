/*
 * Common routines and variables used by Prime95 and NTPrime
 *
 * Comm95a contains information used only during setup
 * Comm95b contains information used only during execution
 * Comm95c contains information used during setup and execution
 */ 

#include <ras.h>
#include <winsock.h>
#include <wininet.h>
#include <process.h>

int	SOCKETS_INITIALIZED = 0;

/* Common routines */

/* Load the PrimeNet DLL, make sure an internet connection is active */

int LoadPrimeNet (void)
{
static	int	RAS_NOT_AVAILABLE = 0;
static	HMODULE	HRAS = 0;
static	DWORD (APIENTRY *RAS_ENUM)(LPRASCONNA, LPDWORD, LPDWORD);
static	DWORD (APIENTRY *RAS_STAT)(HRASCONN, LPRASCONNSTATUSA);
	RASCONN connections[10];
	DWORD	bufsize;
	DWORD	i, num_connections;
	DWORD	ret;

/* Special handling prior to first primenet call. */
/* Init Winsock, requesting version 1.1 */

	if (! SOCKETS_INITIALIZED) {
		static WSADATA zz;
		int	res;
		res = WSAStartup (MAKEWORD (1, 1), &zz);
		if (res != 0) {
			char buf[80];
			sprintf (buf, "ERROR: Winsock initialization returned %d.\n", res);
			OutputStr (COMM_THREAD_NUM, buf);
			return (FALSE);
		}
		SOCKETS_INITIALIZED = 1;
	}

/* If we're not using a dial-up connection, let primenet try */
/* to contact the server. */

	if (!DIAL_UP) return (TRUE);

/* Since Windows 95 can bring up a "Connect To" dialog box */
/* on any call to primenet, we try to make sure we are */
/* already connected before we call primenet.  Otherwise, if */
/* no one is at the computer to respond to the "Connect To" */
/* dialog, the thread hangs until some one does respond. */

/* RAS calls, see below, is no longer the MS-prefered method of detecting */
/* an Internet connection.  Starting in version 22.10 we offer a way for */
/* for users to use the prefered wininet.dll method. */
/* InternetGetConnectedState should return FALSE if the modem is not */
/* connected to the Internet. */
/* Starting in version 25.1, this became the default detection method. */

	if (IniGetInt (INI_FILE, "AlternateModemDetection", 1)) {
		DWORD	flags;
		if (InternetGetConnectedState (&flags, 0)) return (TRUE);
		goto no_modem_connection;
	}

// Unfortunately, the RASAPI32.DLL is not installed on every
// system.  We must load it dynamically.  If the RAS library
// is not found, let primenet.dll try to contact the server.

	if (RAS_NOT_AVAILABLE) return (TRUE);
	if (HRAS == 0) {
		RAS_NOT_AVAILABLE = 1;
		HRAS = LoadLibrary ("rasapi32.dll");
		if (HRAS == 0) return (TRUE);
		RAS_ENUM = (DWORD (APIENTRY *)(LPRASCONNA, LPDWORD, LPDWORD))
			GetProcAddress (HRAS, "RasEnumConnectionsA");
		if (RAS_ENUM == NULL) return (TRUE);
		RAS_STAT = (DWORD (APIENTRY *)(HRASCONN, LPRASCONNSTATUSA))
			GetProcAddress (HRAS, "RasGetConnectStatusA");
		if (RAS_STAT == NULL) return (TRUE);
		RAS_NOT_AVAILABLE = 0;
	}

// Call RAS to see if there are any active connections to the Internet

	connections[0].dwSize = sizeof (RASCONN);
	bufsize = sizeof (connections);
        ret = (*RAS_ENUM) ((RASCONN *) &connections, &bufsize, &num_connections);

// If RAS returns an error who knows what went wrong. 
// Let primenet try to connect anyway.

	if (ret) return (TRUE);

// See if any of these connections are really connected

	for (i = 0; i < num_connections; i++) {
		RASCONNSTATUS status;
		status.dwSize = sizeof (RASCONNSTATUS);
		ret = (*RAS_STAT) (connections[i].hrasconn, &status);
		if (ret) continue;
		if (status.rasconnstate == RASCS_Connected) return (TRUE);
	}

// Print error message if no there are no connections

no_modem_connection:
	OutputStr (COMM_THREAD_NUM, "Dial-up connection not active.\n");
	return (FALSE);
}

/* Unload the PrimeNet DLL */

void UnloadPrimeNet (void)
{

/* Tell winsock we are done. */

	if (SOCKETS_INITIALIZED) {
		// Should we call WSACancelBlockingCall first??
		// Should we check error code from WSACleanup?
		// Should we call WSACleanup after each communication session
		// with the server?  That is, are we tying up any resources?
		WSACleanup ();
		SOCKETS_INITIALIZED = 0;
	}
}

/* Get Windows Serial Number */

void getWindowsSerialNumber (
	char	*output)
{
	HKEY	hkey = 0;
	char	buf[256];
	DWORD	type, disposition;
	DWORD	bufsize = sizeof (buf);

	*output = 0;
	if (RegCreateKeyEx (
			HKEY_LOCAL_MACHINE,
			"Software\\Microsoft\\Windows\\CurrentVersion",
			0,
			NULL,
			REG_OPTION_NON_VOLATILE,
			KEY_ALL_ACCESS,
			NULL,
			&hkey,
			&disposition) == ERROR_SUCCESS &&
	    RegQueryValueEx (hkey, "ProductId", NULL, &type,
			(BYTE *) buf, &bufsize) == ERROR_SUCCESS &&
	    type == REG_SZ)
		strcpy (output, buf);
	if (hkey) RegCloseKey (hkey);
}

/* Get Window's SID.  Hopefully this will combined with the Window's */
/* serial # will generate a unique computer ID.  The SID code */
/* came courtesy of www.sysinternals.com with this copyright */
/* and restriction. */
// Copyright (c) 1997-2002 Mark Russinovich and Bryce Cogswell
//
// Changes the computer SID. 
//
// This code is protected under copyright law. You do not have 
// permission to use this code in a commercial SID-changing product.

PSECURITY_DESCRIPTOR GetRegSecDesc (HKEY Root, TCHAR *Path, 
				    SECURITY_INFORMATION Information)
{
	HKEY					hKey;
	LONG					Status;
	DWORD					nb = 0;
	PSECURITY_DESCRIPTOR	SecDesc;

	//
	// Open the key with no access requests, since we don't need
	// any.
	// SECURITY_DESCRIPTOR
	if (RegOpenKeyEx (Root, Path, 0, KEY_READ, &hKey) != ERROR_SUCCESS)
		return NULL;

	//
	// Grab a copy of the security for key
	//
	if (RegGetKeySecurity (hKey, Information, NULL, &nb) 
					!= ERROR_INSUFFICIENT_BUFFER)
		return NULL;

	SecDesc = malloc (nb);
	Status = RegGetKeySecurity (hKey, Information, SecDesc, &nb);

	//
	// Close the key anyway
	//
	RegCloseKey (hKey);
	if (Status != ERROR_SUCCESS) {
		free (SecDesc);
		return NULL;
	}
	return SecDesc;
}
PSECURITY_DESCRIPTOR GetRegAccess (HKEY hKey)
{
	DWORD		nb = 0;
	PSECURITY_DESCRIPTOR SecDesc;
	//
	// Get access
	//
	if (RegGetKeySecurity (hKey, DACL_SECURITY_INFORMATION, NULL, &nb) != ERROR_INSUFFICIENT_BUFFER)
		return NULL;
	SecDesc = (PSECURITY_DESCRIPTOR) malloc (nb);
	if (RegGetKeySecurity (hKey, DACL_SECURITY_INFORMATION, SecDesc, &nb) != ERROR_SUCCESS) {
		free (SecDesc);
		return (NULL);
	}
	return (SecDesc);
}
LONG SetRegAccess (HKEY hKey, LPCTSTR lpSubKey,
		   PSECURITY_DESCRIPTOR SecDesc, PHKEY phKey)
{
	//
	// Grant requested access
	//
	if (RegSetKeySecurity (*phKey, DACL_SECURITY_INFORMATION, SecDesc) 
					!= ERROR_SUCCESS)
		return FALSE;

	//
	// Re-open the key if requested
	//
	if (! hKey) return TRUE;

	RegCloseKey (*phKey);
	return (RegOpenKey (hKey, lpSubKey, phKey) == ERROR_SUCCESS);
}
PBYTE IsSubAuthValid( PBYTE SidData, DWORD SidLength )
{
	PBYTE	sidPtr;

	sidPtr = NULL;
	if ( SidLength % sizeof(DWORD) == 0 )  {
		for ( sidPtr = SidData + SidLength - 5*sizeof(DWORD); sidPtr >= SidData; sidPtr -= sizeof(DWORD) )
			if ( ((PDWORD)sidPtr)[1] == 0x05000000  &&  ((PDWORD)sidPtr)[2] == 0x00000015 )
				break;
		if ( sidPtr < SidData )
			sidPtr = NULL;
	}
	return sidPtr;
}
void GetTextualSid(
    PSID pSid,            // binary Sid
    PTCHAR TextualSid)    // buffer for Textual representation of Sid
{
    PSID_IDENTIFIER_AUTHORITY psia;
    DWORD dwSubAuthorities;
    DWORD dwSidRev=SID_REVISION;
    DWORD dwCounter;
    DWORD dwSidSize;

    // Validate the binary SID.

    if(!IsValidSid(pSid)) return;

    // Get the identifier authority value from the SID.

    psia = GetSidIdentifierAuthority(pSid);

    // Get the number of subauthorities in the SID.

    dwSubAuthorities = *GetSidSubAuthorityCount(pSid);

    // Compute the buffer length.
    // S-SID_REVISION- + IdentifierAuthority- + subauthorities- + NULL

    dwSidSize=(15 + 12 + (12 * dwSubAuthorities) + 1) * sizeof(TCHAR);

    // Add 'S' prefix and revision number to the string.

    dwSidSize= _stprintf(TextualSid, _T("S-%lu-"), dwSidRev );

    // Add SID identifier authority to the string.

    if ( (psia->Value[0] != 0) || (psia->Value[1] != 0) ) {

        dwSidSize += _stprintf(TextualSid + lstrlen(TextualSid),
                    _T("0x%02hx%02hx%02hx%02hx%02hx%02hx"),
                    (USHORT)psia->Value[0],
                    (USHORT)psia->Value[1],
                    (USHORT)psia->Value[2],
                    (USHORT)psia->Value[3],
                    (USHORT)psia->Value[4],
                    (USHORT)psia->Value[5]);

    } else {

        dwSidSize += _stprintf(TextualSid + lstrlen(TextualSid),
                     _T("%lu"),
                    (ULONG)(psia->Value[5]      )   +
                    (ULONG)(psia->Value[4] <<  8)   +
                    (ULONG)(psia->Value[3] << 16)   +
                    (ULONG)(psia->Value[2] << 24)   );
    }

    // Add SID subauthorities to the string.
    //
    for (dwCounter=0 ; dwCounter < dwSubAuthorities ; dwCounter++) {
        dwSidSize+= _stprintf(TextualSid + dwSidSize, _T("-%lu"),
                    *GetSidSubAuthority(pSid, dwCounter) );
    }
}
void getWindowsSID (
	char	*output)
{
	PSECURITY_DESCRIPTOR	newSecDesc, oldSecDesc;
	DWORD					valType;
	PBYTE					vData;
	HKEY					hKey;
	DWORD					nb;
	PBYTE					sidPtr;
	DWORD					Status;

	*output = 0;
	//
	// Now, get the descriptor of HKLM\SOFTWARE and apply this to SECURITY
	//
	newSecDesc = GetRegSecDesc( HKEY_LOCAL_MACHINE, "SOFTWARE",
					DACL_SECURITY_INFORMATION );

	//
	// Read the last subauthority of the current computer SID
	//
	if( RegOpenKey( HKEY_LOCAL_MACHINE, "SECURITY\\SAM\\Domains\\Account", 
			&hKey) != ERROR_SUCCESS ) {
		free (newSecDesc);
		return;
	}
	oldSecDesc = GetRegAccess( hKey );
	SetRegAccess( HKEY_LOCAL_MACHINE, "SECURITY\\SAM\\Domains\\Account", 
		newSecDesc, &hKey );
	nb = 0;
	vData = NULL;
	RegQueryValueEx( hKey, "V", NULL, &valType, vData, &nb );
	vData = (PBYTE) malloc( nb );
	Status = RegQueryValueEx( hKey, "V", NULL, &valType, vData, &nb );
	if( Status != ERROR_SUCCESS ) {
		SetRegAccess( HKEY_LOCAL_MACHINE, "SECURITY\\SAM\\Domains\\Account",
				oldSecDesc, &hKey );
		free (vData);
		free (oldSecDesc);
		free (newSecDesc);
		return;
	}
	SetRegAccess( NULL, NULL, oldSecDesc, &hKey );
	RegCloseKey( hKey );
	free (oldSecDesc);
	free (newSecDesc);

	//
	// Make sure that we're dealing with a SID we understand
	//
	if( !(sidPtr = IsSubAuthValid( vData, nb ))) {
		free (vData);
		return;
	}

	GetTextualSid( sidPtr, output );
	free (vData);
}



/* Return the number of MB of physical memory.  Windows doesn't return the full amount */
/* of physical memory - probably shadow ram or some such.  Compensate by rounding up. */

unsigned long physical_memory (void)
{
	MEMORYSTATUS mem;

	GlobalMemoryStatus (&mem);
	return ((unsigned long) ((mem.dwTotalPhys + 1000000) >> 20));
}

/* Return a better guess for amount of memory to use in a torture test. */
/* Caller passes in its guess for amount of memory to use, but this routine */
/* can reduce that guess based on OS-specific code that looks at amount */
/* of available physical memory. */
/* This code was written by an anonymous GIMPS user. */

unsigned long GetSuggestedMemory (unsigned long nDesiredMemory)
{
	MEMORYSTATUS ms = {0};

	// In-use Physical RAM in bytes
	SIZE_T dwUsedMem	= ms.dwTotalPhys - ms.dwAvailPhys;
	// Desired memory in bytes
	SIZE_T dwDesiredMem	= nDesiredMemory << 20;
	SIZE_T dwDesiredMemNew	= dwDesiredMem;

	GlobalMemoryStatus (&ms);

	// if very small/no page-file (pagefile <= total RAM) and
	// in-use memory + desired memory > total RAM, then
	// we have to set desired memory to free RAM,
	// because the OS can't page out other apps to
	// reclaim enough free memory since
	// there's not enough space in the pagefile
	// to store the paged-out apps
	if ((ms.dwTotalPageFile <= ms.dwTotalPhys) &&
		  ((dwUsedMem + dwDesiredMem) >= ms.dwTotalPhys))
	{
		dwDesiredMemNew = ms.dwAvailPhys;
	}

	return ((unsigned long) (dwDesiredMemNew >> 20));
}

/* Return the number of CPUs in the system */

unsigned long num_cpus (void)
{
	SYSTEM_INFO sys;

	GetSystemInfo (&sys);
	return (sys.dwNumberOfProcessors);
}

/* Return 1 to print time in AM/PM format.  Return 2 to print */
/* times using a 24-hour clock. */

int getDefaultTimeFormat (void)
{
	char	buf[10];

	GetLocaleInfo (
		LOCALE_USER_DEFAULT, LOCALE_ITIME, (LPTSTR) buf, sizeof (buf));
	return (buf[0] == '0' ? 1 : 2);
}

