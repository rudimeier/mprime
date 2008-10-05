/*
 * Primenet communication routines for all operating systems
 * Uses sockets and HTTP
 */ 

/*
// THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF
// ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A
// PARTICULAR PURPOSE.
//
// Copyright (c) 1997-2008 Mersenne Research, Inc. All Rights Reserved.
//
//  MODULE:   primenet.c
//
//  PURPOSE:  Implements PrimeNet Version 4 and 5 API as HTTP network client
//
//  AUTHOR:   Peter Hunter, on the basis of work by Scott Kurowski (v3 API)
//            Michiel van Loon, OS/2 adaptations 
//            Kurowski 5/1998, 4.0 API support for MPrime 16.x
//            Kurowski 9/1999, 4.0 API changes for MPrime 19.x
//	      Woltman 1/2002, Windows support and bug fixes
//	      Woltman 10/2005, Version 5 API support, CURL library
//
//  ASSUMPTIONS: 1. less than 4k of data is sent or received per call
//               2. HTTP/1.1
//               3. PrimeNet Version 5 or later API on server and client
*/

/* Linux defines, adapted for OS/2, FreeBSD, and Windows */

#define CURL_STATICLIB

#ifdef __WATCOMC__
#include <types.h>
#endif
#include <stdlib.h>
#include <sys/types.h>
#include <string.h>
#include <ctype.h>
#include <stdio.h>

#include "secure5.c"

#ifdef _WINDOWS_
#include <curl.h>
#include <winsock.h>
#else
#include <curl/curl.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <errno.h>
#include <arpa/inet.h>
#ifndef __IBMC__
#include <unistd.h>
#endif
typedef int SOCKET;
#endif

#if defined (__EMX__) && ! defined (AOUT)
#define htonl(x) (lswap(x))
#define ntohl(x) (lswap(x))
#define htons(x) (bswap(x))
#define ntohs(x) (bswap(x))
unsigned short bswap(unsigned short);
#endif

int v4_PRIMENET (short operation, void *pkt);
int v4_format_args (char* args, short operation, void* pkt);
int v4_parse_page (char *response_buf, short operation, void *pkt);

char iniSection[] = "PrimeNet";
char hx[] = "0123456789ABCDEF";
char szSITE[] = "v5.mersenne.org";	/* PrimeNet Server's home domain */
#define nHostPort 80			/* Internet PrimeNet port */
char szFILE[] = "/v5server/?";		/* HTTP GET string */

#define PROXY_HOST_BUFSIZE	120
#define PROXY_USER_BUFSIZE	50
#define PROXY_PASSWORD_BUFSIZE	50

/* implement the missing inet_aton call */

#if defined (_WINDOWS_) || defined (__EMX__)
int inet_aton (char *cp, struct in_addr *inp)
{
	u_long temp;
 
	temp = inet_addr(cp);
	if (temp == -1) return (0);
	inp->s_addr = temp;
	return (1);
}
#endif

/* Routine to get the error number after a failed sockets call */

int getLastSocketError (void)
{
#ifdef _WINDOWS_
	return (WSAGetLastError ());
#else
	return (errno);
#endif
}

/* simple password de/scrambler */

char SCRAMBLE_STRING[] = "/cgi-bin/pnHttp.exe";

void scramble (char *s)
{
	char	out[100];
	char	*p = s, *z = out;
	unsigned int i, c = (unsigned int) strlen (SCRAMBLE_STRING);

	for (i = 0; i < strlen (s); i++) {
		int b = (unsigned char) *p++ ^ SCRAMBLE_STRING[i % c];
		*z++ = hx[b >> 4];
		*z++ = hx[b % 16];
	}
	*z = 0;
	strcpy (s, out);
}

void unscramble (char *s)
{
	char	out[50];
	char	*q = s, *z = out;
	unsigned int i, c = (unsigned int) strlen (SCRAMBLE_STRING);

	for (i = 0; i < strlen (s) >> 1; i++) {
		*z = (char) (strchr (hx, *q++) - hx) * 16;
		*z += (char) (strchr (hx, *q++) - hx);
		*z++ ^= SCRAMBLE_STRING[i % c];
	}
	*z = 0;
	strcpy (s, out);
}

/* base64 encode for basic proxy passwords */

static int encode[] = {
  'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H',
  'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P',
  'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X',
  'Y', 'Z', 'a', 'b', 'c', 'd', 'e', 'f',
  'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n',
  'o', 'p', 'q', 'r', 's', 't', 'u', 'v',
  'w', 'x', 'y', 'z', '0', '1', '2', '3',
  '4', '5', '6', '7', '8', '9', '+', '/'
};

/* Base64-encode a null-terminated string. */

void encode64 (
	char	*data)
{
	char	*s, *end, *buf;
	unsigned int x, length;
	int	i, j;
	char	temp[400];

	length = (unsigned int) strlen (data);
	if (length == 0) return;

	end = data + length - 3;

	buf = temp;
	for (s = data; s < end; ) {
		x = *s++ << 24;
		x |= *s++ << 16;
		x |= *s++ << 8;

		*buf++ = encode[x >> 26];
		x <<= 6;
		*buf++ = encode[x >> 26];
		x <<= 6;
		*buf++ = encode[x >> 26];
		x <<= 6;
		*buf++ = encode[x >> 26];
	}
	end += 3;

	x = 0;
	for (i = 0; s < end; i++)
		x |= *s++ << (24 - 8 * i);

	for (j = 0; j < 4; j++) {
		if (8 * i >= 6 * j) {
			*buf++ = encode [x >> 26];
			x <<= 6;
		} else {
			*buf++ = '=';
		}
	}

	*buf = 0;

	strcpy (data, temp);
}


/* Get proxy information from INI file */

void getProxyInfo (
	char	*szProxyHost,
	unsigned short *nProxyPort,
	char	*szProxyUser,
	char	*szProxyPassword)
{
	char	*colon;

/* Get the host name of the optional proxy server.  If using a proxy */
/* server strip the optional http:// prefix. */

	IniSectionGetString (INI_FILE, iniSection, "ProxyHost",
			     szProxyHost, PROXY_HOST_BUFSIZE, NULL);
	if (szProxyHost[0] == 0) return;

	if ((szProxyHost[0] == 'H' || szProxyHost[0] == 'h') &&
	    (szProxyHost[1] == 'T' || szProxyHost[1] == 't') &&
	    (szProxyHost[2] == 'T' || szProxyHost[2] == 't') &&
	    (szProxyHost[3] == 'P' || szProxyHost[3] == 'p') &&
	    szProxyHost[4] == ':' && szProxyHost[5] == '/' &&
	    szProxyHost[6] == '/')
		strcpy (szProxyHost, szProxyHost + 7);

/* Get optional port number */

	if ((colon = strchr (szProxyHost, ':'))) {
		*nProxyPort = (unsigned short) atoi (colon + 1);
		*colon = 0;
	} else
		*nProxyPort = 8080;

/* Secure proxy - get username and password to negotiate access */

	IniSectionGetString (INI_FILE, iniSection, "ProxyUser",
			     szProxyUser, PROXY_USER_BUFSIZE, NULL);
	IniSectionGetString (INI_FILE, iniSection, "ProxyPass",
			     szProxyPassword, PROXY_PASSWORD_BUFSIZE, NULL);

/* Scramble or unscramble the password as necessary */

	if (!IniSectionGetInt (INI_FILE, iniSection, "ProxyMask", 0)) {
		scramble (szProxyPassword);
		IniSectionWriteString (INI_FILE, iniSection,
					"ProxyPass", szProxyPassword);
		IniSectionWriteInt (INI_FILE, iniSection, "ProxyMask", 1);
	}
	unscramble (szProxyPassword);
}

/*///////////////////////////////////////////////////////////////////////////
//
// HTTP GET procedure (Sockets Implementation)
//
///////////////////////////////////////////////////////////////////////////*/

/*
// pnHttpServer: Uses GET to send a formatted HTTP argument string
//               and downloads the server result page
*/

int pnHttpServer (char *pbuf, unsigned cbuf, char* postargs)
{
	char	szBuffer[1000];
	char	szURL[4096];			/* URL assembly buffer */
	char	buf[4096];
	struct in_addr defaddr;
	struct hostent *hp, def;
	struct sockaddr_in sn;
	SOCKET	s;
	int	res, debug, url_format;
	int	timeout;
	char	*alist[1];
	unsigned int count;
	char	szProxyHost[PROXY_HOST_BUFSIZE];
	char	szProxyUser[PROXY_USER_BUFSIZE];
	char	szProxyPassword[PROXY_PASSWORD_BUFSIZE];
	char	*con_host;
	char	szOtherGetInfo[256];
	unsigned short nProxyPort, con_port;

/* Get debug logging and URL format flags */

	debug = IniSectionGetInt (INI_FILE, iniSection, "Debug", 0);
	url_format = IniSectionGetInt (INI_FILE, iniSection, "UseFullURL", 2);
 
/* Get information about the optional proxy server */

	getProxyInfo (szProxyHost, &nProxyPort, szProxyUser, szProxyPassword);
	if (szProxyHost[0]) {
		con_host = szProxyHost;
		con_port = nProxyPort;
	}

/* No proxy server - use default site (mersenne.org) */

	else {
		con_host = szSITE;
		con_port = nHostPort;
	}

/* Output debug info */

redirect:
	if (debug) {
		sprintf (buf, "host = %s, port = %d\n", con_host, con_port);
		LogMsg (buf);
	}

/* Convert host name into an IP address */

	hp = gethostbyname (con_host);
	if (!hp) {
		char	szAltSiteAddr[30];
		if (con_host == szSITE) {
			IniSectionGetString (INI_FILE, iniSection,
					"MersenneIP", szAltSiteAddr, 29, NULL);
			con_host = szAltSiteAddr;
		}
		if (!inet_aton (con_host, &defaddr)) {
			sprintf (buf, "Unknown host: %s\n", con_host);
			if (debug) LogMsg (buf);
			else OutputStr (COMM_THREAD_NUM, buf);
			return (PRIMENET_ERROR_CONNECT_FAILED);
		}
		alist[0] = (char *) &defaddr;
		def.h_name = szSITE;
		def.h_addr_list = alist;
		def.h_length = sizeof (struct in_addr);
		def.h_addrtype = AF_INET;
		def.h_aliases = 0;
		hp = &def;
	}

	if (debug) {
		sprintf (buf, "IP-addr = %s\n",
			 inet_ntoa (*(struct in_addr *)hp->h_addr));
		LogMsg (buf);
	}

	memset (&sn, 0, sizeof (sn));
	sn.sin_family = hp->h_addrtype;
	if (hp->h_length > (int) sizeof (sn.sin_addr)) {
		hp->h_length = sizeof (sn.sin_addr);
	}
	memcpy (&sn.sin_addr, hp->h_addr, hp->h_length);
	sn.sin_port = htons (con_port);

/* Create a socket and connect to server */

rel_url:
	if ((s = socket (hp->h_addrtype, SOCK_STREAM, 0)) < 0) {
		sprintf (buf, "Error in socket call: %d\n",
			 getLastSocketError ());
		if (debug) LogMsg (buf);
		else OutputStr (COMM_THREAD_NUM, buf);
		return (PRIMENET_ERROR_CONNECT_FAILED);
	}

	if (connect (s, (struct sockaddr *) &sn, sizeof (sn)) < 0) {
		sprintf (buf, "Error in connect call: %d\n",
			 getLastSocketError ());
		if (debug) LogMsg (buf);
		else OutputStr (COMM_THREAD_NUM, buf);
		closesocket (s);
		return (PRIMENET_ERROR_CONNECT_FAILED);
	}

/* Prevent SIGPIPE signals in Linux (and other) environments */

#ifdef SO_NOSIGPIPE
	{
		int	i = 1;
		res = setsockopt (s, SOL_SOCKET, SO_NOSIGPIPE, &i, sizeof(i));
		if (res < 0) {
			if (debug) {
				sprintf (buf, "Error in NOSIGPIPE call: %d\n",
					 getLastSocketError ());
				LogMsg (buf);
			}
		}
	}
#endif

/* GET method, data follows ? in URL */

	strcpy (szURL, "GET ");
	if (*szProxyHost || url_format == 1) {
		strcat (szURL, "http://");
		strcat (szURL, szSITE);
		if (IniSectionGetInt (INI_FILE, iniSection, "SendPortNumber", 1))
			sprintf (szURL + strlen (szURL), ":%d", nHostPort);
	}
	strcat (szURL, szFILE);
	strcat (szURL, postargs);
	strcat (szURL, " HTTP/1.1\r\n");

/* Append host: header */

	if (!*szProxyHost)
		sprintf (szURL+strlen(szURL), "Host: %s\r\n", szSITE);

/* Persistent connections are not supported */

	strcat (szURL, "Connection: close\r\n");

/* Append other GET info */

	IniSectionGetString (INI_FILE, iniSection, "OtherGetInfo", szOtherGetInfo,
		sizeof (szOtherGetInfo), NULL);
	if (*szOtherGetInfo) {
		strcat (szURL, szOtherGetInfo);
		strcat (szURL, "\r\n");
	}

/* Append proxy authorization here */

	if (*szProxyHost && *szProxyUser) {
		char	buf[200];
		strcat (szURL, "Proxy-Authorization: Basic ");
		sprintf (buf, "%s:%s", szProxyUser, szProxyPassword);
		encode64 (buf);
		strcat (szURL, buf);
		strcat (szURL, "\r\n");
	}

	strcat (szURL, "\r\n");
	if (debug) LogMsg (szURL);

/* Send the URL request */

	timeout = 90000;		/* 90 seconds */
	res = setsockopt (s, SOL_SOCKET, SO_SNDTIMEO, (char *) &timeout,
			  sizeof (timeout));
	if (res < 0) {
		if (debug) {
			sprintf (buf, "Error in send timeout call: %d\n",
				 getLastSocketError ());
			LogMsg (buf);
		}
	}
	res = send (s, szURL, (int) strlen (szURL), 0);
	if (res < 0) {
		sprintf (buf, "Error in send call: %d\n",
			 getLastSocketError ());
		if (debug) LogMsg (buf);
		else OutputStr (COMM_THREAD_NUM, buf);
		closesocket (s);
		if (url_format == 2 && *szProxyHost == 0) {
			if (debug) LogMsg ("Trying full URL\n");
			url_format = 1;
			goto rel_url;
		}
		return (PRIMENET_ERROR_SEND_FAILED);
	}

/* Now accumulate the response */

	timeout = 90000;		/* 90 seconds */
	res = setsockopt (s, SOL_SOCKET, SO_RCVTIMEO, (char *) &timeout,
			  sizeof (timeout));
	if (res < 0) {
		if (debug) {
			sprintf (buf, "Error in receive timeout call: %d\n",
				 getLastSocketError ());
			LogMsg (buf);
		}
	}
	*pbuf = 0; count = 1;
	while (count < cbuf) {
		res = recv (s, szBuffer, 999, 0);
		if (res < 0) {
			sprintf (buf, "Error in recv call: %d\n",
				 getLastSocketError ());
			if (debug) LogMsg (buf);
			else OutputStr (COMM_THREAD_NUM, buf);
			closesocket (s);
			if (url_format == 2 && *szProxyHost == 0) {
				if (debug) LogMsg ("Trying full URL\n");
				url_format = 1;
				goto rel_url;
			}
			return (PRIMENET_ERROR_RECV_FAILED);
		}
		if (res == 0) break;
		szBuffer[res] = 0;
		if (debug) {
			sprintf (buf, "RECV: %s\n", szBuffer);
			LogMsg (buf);
		}
		if ((count + strlen(szBuffer)) >= cbuf)
			szBuffer[cbuf - count] = 0;
		strcat (pbuf, szBuffer);
		count += (unsigned int) strlen (szBuffer);
	}

	closesocket (s);

/* pbuf + 9 is where the message code following HTTP/1.1 starts */

	if (count <= 10) res = -1;
	else res = atoi (pbuf + 9);

/* Some proxy servers can redirect us to another host.  These are */
/* the 300 series of error codes.  This code probably isn't right. */
/* We can improve it as people find problems. */

	if (res >= 300 && res <= 399) {
		char	*location, *colon;
		location = strstr (pbuf, "Location:");
		if (location != NULL) {
			if (debug) LogMsg ("Attempting a redirect.\n");
			location += 9;
			while (isspace (*location)) location++;
			strcpy (szProxyHost, location);

/* Parse the redirection address */

			if ((location[0] == 'H' || location[0] == 'h') &&
			    (location[1] == 'T' || location[1] == 't') &&
			    (location[2] == 'T' || location[2] == 't') &&
			    (location[3] == 'P' || location[3] == 'p') &&
			    location[4] == ':' && location[5] == '/' &&
			    location[6] == '/')
				strcpy (location, location + 7);
			con_host = location;

/* Get optional port number */

			if ((colon = strchr (location, ':')) != NULL) {
				con_port = (unsigned short) atoi (colon + 1);
				*colon = 0;
			} else
				con_port = 80;
			goto redirect;
		}
	}

/* Any return code other than 200 is an error.  We've had problems using */
/* both full URLs and relative URLs.  Thus, our default behavior is to try */
/* both before giving up. */

	if (res != 200) {
		if (url_format == 2 && *szProxyHost == 0) {
			if (debug) LogMsg ("Trying full URL\n");
			url_format = 1;
			goto rel_url;
		}
		strcpy (buf, "HTTP return code is not 200\n");
		if (debug) LogMsg (buf);
		else OutputStr (COMM_THREAD_NUM, buf);
		return (PRIMENET_ERROR_SERVER_UNSPEC);
	}

	return (PRIMENET_NO_ERROR);
}


/*///////////////////////////////////////////////////////////////////////////
//
// HTTP GET procedure (cURL Implementation)
//
///////////////////////////////////////////////////////////////////////////*/

/* This callback routine assembles the server's response to our request */

size_t WriteMemoryCallback (
	void	*ptr,
	size_t	size,
	size_t	nmemb,
	void	*data)
{
	size_t realsize = size * nmemb;
	char	*buf = (char *) data;
	int	buflen = (int) strlen (buf);

/* Truncate response to fit in a 4096 byte buffer */

	if (buflen + realsize <= 4095) {
		memcpy (buf + buflen, ptr, realsize);
		buf[buflen + realsize] = 0;
	} else {
		memcpy (buf + buflen, ptr, 4095 - buflen);
		buf[4095] = 0;
	}
	return realsize;
}

/* This callback is for dumping out cURL debug information */

int curl_trace (
	CURL	*handle,
	curl_infotype type,
	unsigned char *data,
	size_t	size,
	void	*userp)
{
	char	*text;
	char	buf[4096];
	int	len;

	switch (type) {
	case CURLINFO_TEXT:
		text = "== Info";
		break;
	case CURLINFO_HEADER_OUT:
		text = "=> Send header";
		break;
	case CURLINFO_DATA_OUT:
		text = "=> Send data";
		break;
	case CURLINFO_HEADER_IN:
		text = "<= Recv header";
		break;
	case CURLINFO_DATA_IN:
		text = "<= Recv data";
		break;
	case CURLINFO_SSL_DATA_IN:
		text = "<= Recv SSL data";
		break;
	case CURLINFO_SSL_DATA_OUT:
		text = "<= Send SSL data";
		break;
	default: /* in case a new one is introduced to shock us */
		return 0;
	}

/* Output the data. */

	strcpy (buf, text);
	strcat (buf, ": ");
	len = (int) strlen (buf);
	memcpy (buf + len, data, size);
	if (data[size-1] != '\n') {
		buf[len + size] = '\n';
		buf[len + size + 1] = 0;
	} else
		buf[len + size] = 0;
	LogMsg (buf);

	return 0;
}

/*
// pnHttpServerCURL: Uses GET to send a formatted HTTP argument string
//               and downloads the server result page
*/

int pnHttpServerCURL (char *pbuf, unsigned cbuf, char* postargs)
{
	CURL	*curl;
	CURLcode res;
	int	debug;
	char	szAltSiteAddr[120];
	char	url[4096], buf[4150], errbuf[CURL_ERROR_SIZE];
	char	szProxyHost[PROXY_HOST_BUFSIZE];
	char	szProxyUser[PROXY_USER_BUFSIZE];
	char	szProxyPassword[PROXY_PASSWORD_BUFSIZE];
	unsigned short nProxyPort;

/* Init the cURL structures */

	curl = curl_easy_init ();
	if (curl == NULL) return (PRIMENET_ERROR_CURL_INIT);
	curl_easy_setopt (curl, CURLOPT_NOPROGRESS, 1);

/* Get debug logging flag */

	debug = IniSectionGetInt (INI_FILE, iniSection, "Debug", 0);
 
/* Give curl library the HTTP string to send */

	strcpy (url, "http://");
	IniSectionGetString (INI_FILE, iniSection, "MersenneIP",
			     szAltSiteAddr, sizeof (szAltSiteAddr), NULL);
	if (szAltSiteAddr[0])
		strcat (url, szAltSiteAddr);
	else if (USE_V4)
		strcat (url, "mersenne.org");
	else
		strcat (url, szSITE);
	if (IniSectionGetInt (INI_FILE, iniSection, "SendPortNumber", 0))
		sprintf (url + strlen (url), ":%d", nHostPort);
	if (USE_V4)
		strcat (url, "/cgi-bin/pnHttp.exe?");
	else
		strcat (url, szFILE);
	strcat (url, postargs);
	curl_easy_setopt (curl, CURLOPT_URL, url);
	if (debug) {
		sprintf (buf, "URL: %s\n", url);
		LogMsg (buf);
	}

/* Get information about the optional proxy server */

	getProxyInfo (szProxyHost, &nProxyPort, szProxyUser, szProxyPassword);
	if (szProxyHost[0]) {
		curl_easy_setopt (curl, CURLOPT_PROXY, szProxyHost);
		curl_easy_setopt (curl, CURLOPT_PROXYPORT, (long) nProxyPort);
//bug?		curl_easy_setopt (curl, CURLOPT_PROXYTYPE, ???);
		if (szProxyUser[0]) {
			sprintf (buf, "%s:%s", szProxyUser, szProxyPassword);
			curl_easy_setopt (curl, CURLOPT_PROXYUSERPWD, buf);
			curl_easy_setopt (curl, CURLOPT_PROXYAUTH, CURLAUTH_ANY);
		}
	}

/* Setup function to receive the response */

	curl_easy_setopt (curl, CURLOPT_WRITEFUNCTION, WriteMemoryCallback);
	curl_easy_setopt (curl, CURLOPT_WRITEDATA, (void *) pbuf);
	pbuf[0] = 0;

/* Output verbose debug information */

	if (debug >= 2) {
		curl_easy_setopt (curl, CURLOPT_DEBUGFUNCTION, curl_trace);
		curl_easy_setopt (curl, CURLOPT_DEBUGDATA, NULL);
		/* the DEBUGFUNCTION has no effect until we enable VERBOSE */
		curl_easy_setopt (curl, CURLOPT_VERBOSE, 1);
	}

/* Send the URL request */

	curl_easy_setopt (curl, CURLOPT_FOLLOWLOCATION, 1);
	curl_easy_setopt (curl, CURLOPT_NOSIGNAL, 1);
	curl_easy_setopt (curl, CURLOPT_CONNECTTIMEOUT, 180);
	curl_easy_setopt (curl, CURLOPT_TIMEOUT, 180);
	curl_easy_setopt (curl, CURLOPT_ERRORBUFFER, errbuf);
	res = curl_easy_perform (curl);
	if (res != CURLE_OK) {
		sprintf (buf, "Unexpected CURL library error: %s\n", errbuf);
		LogMsg (buf);
		OutputStr (COMM_THREAD_NUM, buf);
		curl_easy_cleanup (curl);
		return (PRIMENET_ERROR_CURL_PERFORM);
	}

/* Cleanup */

	curl_easy_cleanup (curl);

/* Log a debug message */

	if (debug) {
		LogMsg ("RESPONSE:\n");
		LogMsg (pbuf);
	}

/* Return success */

	return (PRIMENET_NO_ERROR);
}



/*///////////////////////////////////////////////////////////////////////////
//
// HTTP GET argument formatting procedures
//
////////////////////////////////////////////////////////////////////////////*/

/* armor parameter control chars as hex codes for transport */

#define ARMOR_CHARS		"&+%\r\n"

char *armor (char *d, char *s)
{

/* & is token delimiter, '+' is space char */

	while (*s) {
		if (strchr (ARMOR_CHARS, *s)) {	
			*d++ = '%';	/* convert chars to %nn hex codes */
			*d++ = hx[(*s) / 16];
			*d++ = hx[(*s) % 16];
		} else if (*s == ' ')	/* convert spaces to '+' */
			*d++ = '+';
		else *d++ = *s;		/* copy normal character */
		s++;
	}
	*d = 0;
	return (d);
}

/*
// format_args: format a HTTP argument string from a PrimeNet v5 packet
*/

int format_args (char* args, short operation, void* pkt)
{
	char	*p;

/* Format the common header */

	sprintf (args, "v=%.2f&px=GIMPS", PRIMENET_TRANSACTION_API_VERSION);
	p = args + strlen (args);

/* Format the message dependent args */

	switch (operation) {
	case PRIMENET_UPDATE_COMPUTER_INFO:	/* update computer info */
		{
		struct primenetUpdateComputerInfo *z;

		z = (struct primenetUpdateComputerInfo *) pkt;

//		if (!_stricmp (z->user_id, "ANONYMOUS"))
//			strcpy (z->user_id, "admin_user_anon");

		strcpy (p, "&t=uc&g=");
		p = armor (p + strlen (p), z->computer_guid);
		strcpy (p, "&hg=");
		p = armor (p + strlen (p), z->hardware_guid);
		strcpy (p, "&wg=");
		p = armor (p + strlen (p), z->windows_guid);
		strcpy (p, "&a=");
		p = armor (p + strlen (p), z->application);
		strcpy (p, "&c=");
		p = armor (p + strlen (p), z->cpu_model);
		strcpy (p, "&f=");
		p = armor (p + strlen (p), z->cpu_features);
		sprintf (p, "&L1=%d&L2=%d&np=%d&hp=%d&m=%d&s=%d&h=%d&r=%d",
			 z->L1_cache_size, z->L2_cache_size, z->num_cpus,
			 z->num_hyperthread, z->mem_installed, z->cpu_speed,
			 z->hours_per_day, z->rolling_average);
		p += strlen (p);
		if (z->L3_cache_size > 0) {
			sprintf (p, "&L3=%d", z->L3_cache_size);
			p += strlen (p);
		}
		if (z->user_id[0]) {
			strcpy (p, "&u=");
			p = armor (p + strlen (p), z->user_id);
		}
		if (z->computer_name[0]) {
			strcpy (p, "&cn=");
			p = armor (p + strlen (p), z->computer_name);
		}
		break;
		}
	case PRIMENET_PROGRAM_OPTIONS:
		{
		struct primenetProgramOptions *z;

		z = (struct primenetProgramOptions *) pkt;
		strcpy (p, "&t=po&g=");
		p = armor (p + strlen (p), z->computer_guid);
		if (z->cpu_num != -1) {
			sprintf (p, "&c=%d", z->cpu_num);
			p += strlen (p);
		}
		if (z->num_workers != -1) {
			sprintf (p, "&nw=%d", z->num_workers);
			p += strlen (p);
		}
		if (z->work_preference != -1) {
			sprintf (p, "&w=%d", z->work_preference);
			p += strlen (p);
		}
		if (z->priority != -1) {
			sprintf (p, "&Priority=%d", z->priority);
			p += strlen (p);
		}
		if (z->daysOfWork != -1) {
			sprintf (p, "&DaysOfWork=%d", z->daysOfWork);
			p += strlen (p);
		}
		if (z->dayMemory != -1) {
			sprintf (p, "&DayMemory=%d", z->dayMemory);
			p += strlen (p);
		}
		if (z->nightMemory != -1) {
			sprintf (p, "&NightMemory=%d", z->nightMemory);
			p += strlen (p);
		}
		if (z->dayStartTime != -1) {
			sprintf (p, "&DayStartTime=%d", z->dayStartTime);
			p += strlen (p);
		}
		if (z->nightStartTime != -1) {
			sprintf (p, "&NightStartTime=%d", z->nightStartTime);
			p += strlen (p);
		}
		if (z->runOnBattery != -1) {
			sprintf (p, "&RunOnBattery=%d", z->runOnBattery);
			p += strlen (p);
		}
		break;
		}
	case PRIMENET_REGISTER_ASSIGNMENT:	/* register assignment */
		{
		struct primenetRegisterAssignment *z;

		z = (struct primenetRegisterAssignment *) pkt;
		strcpy (p, "&t=ra&g=");
		p = armor (p + strlen (p), z->computer_guid);
		sprintf (p, "&c=%d&w=%d", z->cpu_num, z->work_type);
		p = p + strlen (p);
		if (z->work_type == PRIMENET_WORK_TYPE_FACTOR) {
			sprintf (p, "&n=%d&sf=%g",
				 z->n, z->how_far_factored);
			p = p + strlen (p);
			if (z->factor_to != 0.0) {
				sprintf (p, "&ef=%g", z->factor_to);
				p = p + strlen (p);
			}
		}
		if (z->work_type == PRIMENET_WORK_TYPE_PFACTOR) {
			sprintf (p, "&A=%.0f&b=%d&n=%d&C=%d&sf=%g&saved=%g",
				 z->k, z->b, z->n, z->c, z->how_far_factored,
				 z->tests_saved);
			p = p + strlen (p);
		}
		if (z->work_type == PRIMENET_WORK_TYPE_FIRST_LL ||
		    z->work_type == PRIMENET_WORK_TYPE_DBLCHK) {
			sprintf (p, "&n=%d&sf=%g&p1=%d",
				 z->n, z->how_far_factored,
				 z->has_been_pminus1ed);
			p = p + strlen (p);
		}
		if (z->work_type == PRIMENET_WORK_TYPE_PMINUS1) {
			sprintf (p, "&A=%.0f&b=%d&n=%d&C=%d&B1=%.0f",
				 z->k, z->b, z->n, z->c, z->B1);
			p = p + strlen (p);
			if (z->B2 != 0.0) {
				sprintf (p, "&B2=%.0f", z->B2);
				p = p + strlen (p);
			}
		}
		if (z->work_type == PRIMENET_WORK_TYPE_ECM) {
			sprintf (p, "&A=%.0f&b=%d&n=%d&C=%d&B1=%.0f",
				 z->k, z->b, z->n, z->c, z->B1);
			p = p + strlen (p);
			if (z->B2 != 0.0) {
				sprintf (p, "&B2=%.0f", z->B2);
				p = p + strlen (p);
			}
			sprintf (p, "&CR=%d", z->curves);
			p = p + strlen (p);
		}
		if (z->work_type == PRIMENET_WORK_TYPE_PRP) {
			sprintf (p, "&A=%.0f&b=%d&n=%d&C=%d&sf=%g&saved=%g",
				 z->k, z->b, z->n, z->c, z->how_far_factored,
				 z->tests_saved);
			p = p + strlen (p);
		}
		break;
		}
	case PRIMENET_GET_ASSIGNMENT:		/* get assignment */
		{
		struct primenetGetAssignment *z;

		z = (struct primenetGetAssignment *) pkt;
		strcpy (p, "&t=ga&g=");
		p = armor (p + strlen (p), z->computer_guid);
		sprintf (p, "&c=%d", z->cpu_num);
		p += strlen (p);
		break;
		}
	case PRIMENET_ASSIGNMENT_PROGRESS:
		{
		struct primenetAssignmentProgress *z;

		z = (struct primenetAssignmentProgress *) pkt;
		strcpy (p, "&t=ap&g=");
		p = armor (p + strlen (p), z->computer_guid);
		strcpy (p, "&k=");
		p = armor (p + strlen (p), z->assignment_uid);
		if (z->stage[0]) {
			strcpy (p, "&stage=");
			p = armor (p + strlen (p), z->stage);
		}
		/* Server does not like a pcercent complete of 100%. */
		/* Just in case caller passes that value in, convert it */
		sprintf (p, "&c=%lu&p=%f&d=%lu&e=%lu",
			 z->cpu_num,
			 z->pct_complete < 100.0 ? z->pct_complete : 99.99,
			 z->next_update, z->end_date);
		p += strlen (p);
		if (z->fftlen) {
			sprintf (p, "&fftlen=%d", z->fftlen);
			p += strlen (p);
		}
		break;
		}
	case PRIMENET_ASSIGNMENT_RESULT:
		{
		struct primenetAssignmentResult *z;

		z = (struct primenetAssignmentResult *) pkt;
		strcpy (p, "&t=ar&g=");
		p = armor (p + strlen (p), z->computer_guid);
		if (z->assignment_uid[0]) {
			strcpy (p, "&k=");
			p = armor (p + strlen (p), z->assignment_uid);
		} else {
			strcpy (p, "&k=0");
			p += strlen (p);
		}
		if (z->message[0]) {
			strcpy (p, "&m=");
			p = armor (p + strlen (p), z->message);
		}
		sprintf (p, "&r=%d&d=%d", z->result_type, z->done);
		p += strlen (p);
		if (z->result_type == PRIMENET_AR_LL_RESULT) {
			sprintf (p, "&n=%d&sc=%d", z->n, z->shift_count);
			p += strlen (p);
			strcpy (p, "&rd=");
			p = armor (p + strlen (p), z->residue);
			strcpy (p, "&ec=");
			p = armor (p + strlen (p), z->error_count);
		}
		if (z->result_type == PRIMENET_AR_LL_PRIME) {
			sprintf (p, "&n=%d&sc=%d", z->n, z->shift_count);
			p += strlen (p);
			strcpy (p, "&ec=");
			p = armor (p + strlen (p), z->error_count);
		}
		if (z->result_type == PRIMENET_AR_TF_FACTOR) {
			sprintf (p, "&n=%d&sf=%g", z->n, z->start_bits);
			p += strlen (p);
			strcpy (p, "&f=");
			p = armor (p + strlen (p), z->factor);
		}
		if (z->result_type == PRIMENET_AR_P1_FACTOR) {
			sprintf (p, "&A=%.0f&b=%d&n=%d&c=%d&B1=%.0f",
				 z->k, z->b, z->n, z->c, z->B1);
			p += strlen (p);
			if (z->B2) {
				sprintf (p, "&B2=%.0f", z->B2);
				p += strlen (p);
			}
			strcpy (p, "&f=");
			p = armor (p + strlen (p), z->factor);
		}
		if (z->result_type == PRIMENET_AR_ECM_FACTOR) {
			sprintf (p, "&A=%.0f&b=%d&n=%d&c=%d&CR=%d&B1=%.0f&stage=%d",
				 z->k, z->b, z->n, z->c, z->curves, z->B1, z->stage);
			p += strlen (p);
			if (z->B2) {
				sprintf (p, "&B2=%.0f", z->B2);
				p += strlen (p);
			}
			strcpy (p, "&f=");
			p = armor (p + strlen (p), z->factor);
		}

		if (z->result_type == PRIMENET_AR_TF_NOFACTOR) {
			sprintf (p, "&n=%d&sf=%g&ef=%g",
				 z->n, z->start_bits, z->end_bits);
			p += strlen (p);
		}
		if (z->result_type == PRIMENET_AR_P1_NOFACTOR) {
			sprintf (p, "&A=%.0f&b=%d&n=%d&c=%d&B1=%.0f",
				 z->k, z->b, z->n, z->c, z->B1);
			p += strlen (p);
			if (z->B2) {
				sprintf (p, "&B2=%.0f", z->B2);
				p += strlen (p);
			}
		}
		if (z->result_type == PRIMENET_AR_ECM_NOFACTOR) {
			sprintf (p, "&A=%.0f&b=%d&n=%d&c=%d&CR=%d&B1=%.0f",
				 z->k, z->b, z->n, z->c, z->curves, z->B1);
			p += strlen (p);
			if (z->B2) {
				sprintf (p, "&B2=%.0f", z->B2);
				p += strlen (p);
			}
		}
		if (z->result_type == PRIMENET_AR_PRP_RESULT) {
			sprintf (p, "&A=%.0f&b=%d&n=%d&c=%d",
				 z->k, z->b, z->n, z->c);
			p += strlen (p);
			strcpy (p, "&rd=");
			p = armor (p + strlen (p), z->residue);
			strcpy (p, "&ec=");
			p = armor (p + strlen (p), z->error_count);
		}
		if (z->result_type == PRIMENET_AR_PRP_PRIME) {
			sprintf (p, "&A=%.0f&b=%d&n=%d&c=%d",
				 z->k, z->b, z->n, z->c);
			p += strlen (p);
			strcpy (p, "&ec=");
			p = armor (p + strlen (p), z->error_count);
		}
		if (z->fftlen) {
			sprintf (p, "&fftlen=%d", z->fftlen);
			p += strlen (p);
		}
//bug - can we support a 0 result_type that only sends a message.
//might need this for sending old results file
		break;
		}
	case PRIMENET_ASSIGNMENT_UNRESERVE:	/* assignment unreserve */
		{
		struct primenetAssignmentUnreserve *z;

		z = (struct primenetAssignmentUnreserve *) pkt;
		strcpy (p, "&t=au&g=");
		p = armor (p + strlen (p), z->computer_guid);
		strcpy (p, "&k=");
		p = armor (p + strlen (p), z->assignment_uid);
		break;
		}
	case PRIMENET_BENCHMARK_DATA:
		{
		struct primenetBenchmarkData *z;
		unsigned int i;

		z = (struct primenetBenchmarkData *) pkt;
		strcpy (p, "&t=bd&g=");
		p = armor (p + strlen (p), z->computer_guid);
		if (z->user_comment[0]) {
			strcpy (p, "&c=");
			p = armor (p + strlen (p), z->user_comment);
		}
		for (i = 0; i < z->num_data_points; i++) {
			*p++ = '&';
			p = armor (p, z->data_points[i].bench);
			sprintf (p, "=%g", z->data_points[i].timing);
			p += strlen (p);
		}
		break;
		}
	case PRIMENET_PING_SERVER:
		{
		struct primenetPingServer *z;

		z = (struct primenetPingServer *) pkt;
		sprintf (p, "&t=ps&q=%d", z->ping_type);
		p += strlen (p);
		break;
		}
	}

/* Append the security string */

#ifdef _V5_SECURITY_MODULE_PRESENT_
	{
		char	p_v5key[33];
		make_v5_client_key (p_v5key, COMPUTER_GUID);
		secure_v5_url (args, p_v5key);
	}
#endif
	return (0);
}


/*////////////////////////////////////////////////////////////////////////////
//
// HTTP downloaded response page parsing procedures
//
/////////////////////////////////////////////////////////////////////////////*/

/* skip over the token name and point to the data string */

char* skip_token (char *s)
{
	while (*s && *s != '=') s++;
	if (*s == '=') s++;
	return (s);
}


/* copy the data string up to the next '\r' delimiter character */

char* copy_value (char *buf, char *s)
{
	while (*s && *s != '\r') *buf++ = *s++;
	if (*s == '\r') s++;
	*buf = 0;
	return (s);
}

/* parse various tokens from the response */

char *find_id (
	char	*buf,
	char	*id)
{
	unsigned int idlen;
	idlen = (unsigned int) strlen (id);
	while (*buf) {
		if (memcmp (buf, id, idlen) == 0 && buf[idlen] == '=')
			return (buf + idlen + 1);
		while (*buf && *buf != '>' && *buf != '\n') buf++;
		while (*buf && (*buf == '>' || *buf == '\r' || *buf == '\n' || *buf == ' ')) buf++;
	}
	return (NULL);
}

int parse_string (
	char	*buf,
	char	*id,
	char	*result_buf,
	unsigned int result_buflen)
{
	unsigned int i;
	buf = find_id (buf, id);
	if (buf == NULL) return (FALSE);
	for  (i = 0; i < result_buflen-1 && buf[i] && buf[i] != '\n'; i++)
		result_buf[i] = (buf[i] == '\r' ? '\n' : buf[i]);
	result_buf[i] = 0;
//bug -raise error if string too long
	return (TRUE);
}

int parse_int (
	char	*buf,
	char	*id,
	int32_t	*result)
{
	buf = find_id (buf, id);
	if (buf == NULL) return (FALSE);
	if (buf[0] != '-' && (buf[0] < '0' || buf[0] > '9')) return (FALSE);
	*result = atoi (buf);
//bug -raise error if not integer
	return (TRUE);
}

int parse_uint (
	char	*buf,
	char	*id,
	uint32_t *result)
{
	buf = find_id (buf, id);
	if (buf == NULL) return (FALSE);
	if (buf[0] < '0' || buf[0] > '9') return (FALSE);
	*result = atoi (buf);
//bug -raise error if not integer
	return (TRUE);
}

void parse_double (
	char	*buf,
	char	*id,
	double	*result)
{
	buf = find_id (buf, id);
	if (buf == NULL) return;
	*result = atof (buf);
//bug -raise error if not number
}
	


/*
// parse_page: reads the server response page tokens and values
//             and converts these back into a C structure
*/

int parse_page (char *response_buf, short operation, void *pkt)

{
	char	*s;
	char	buf[400], errtxt[200];
	int32_t	res;

/* Get result code, which is always first */

	s = response_buf;
	if (!parse_int (s, "pnErrorResult", &res)) {
		LogMsg ("PnErrorResult value missing.  Full response was:\n");
		LogMsg (response_buf);

		/* Look for PHP timeout response and convert it to */
		/* a server busy error code. */

		if (strstr (s, "execution time") != NULL &&
		    strstr (s, "exceeded") != NULL)
			return (PRIMENET_ERROR_SERVER_BUSY);
			
		return (PRIMENET_ERROR_PNERRORRESULT);
	}
	if (!parse_string (s, "pnErrorDetail", errtxt, sizeof (errtxt))) {
		LogMsg ("PnErrorDetail string missing\n");
		return (PRIMENET_ERROR_PNERRORDETAIL);
	}

/* If result is non-zero print out an error message */

	if (res) {
		char	buf[2000];
		char	*resmsg;

/* Convert the error number to text */

		switch (res) {
		case PRIMENET_ERROR_SERVER_BUSY:
			resmsg = "Server busy";
			break;
		case PRIMENET_ERROR_INVALID_VERSION:
			resmsg = "Invalid version";
			break;
		case PRIMENET_ERROR_INVALID_TRANSACTION:
			resmsg = "Invalid transaction";
			break;
		case PRIMENET_ERROR_INVALID_PARAMETER:
			resmsg = "Invalid parameter";
			break;
		case PRIMENET_ERROR_ACCESS_DENIED:
			resmsg = "Access denied";
			break;
		case PRIMENET_ERROR_DATABASE_CORRUPT:
			resmsg = "Server database malfunction";
			break;
		case PRIMENET_ERROR_DATABASE_FULL_OR_BROKEN:
			resmsg = "Server database full or broken";
			break;
		case PRIMENET_ERROR_INVALID_USER:
			resmsg = "Invalid user";
			break;
		case PRIMENET_ERROR_UNREGISTERED_CPU:
			resmsg = "CPU not registered";
			break;
		case PRIMENET_ERROR_OBSOLETE_CLIENT:
			resmsg = "Obsolete client, please upgrade";
			break;
		case PRIMENET_ERROR_STALE_CPU_INFO:
			resmsg = "Stale cpu info";
			break;
		case PRIMENET_ERROR_CPU_IDENTITY_MISMATCH:
			resmsg = "CPU identity mismatch";
			break;
		case PRIMENET_ERROR_CPU_CONFIGURATION_MISMATCH:
			resmsg = "CPU configuration mismatch";
			break;
		case PRIMENET_ERROR_NO_ASSIGNMENT:
			resmsg = "No assignment";
			break;
		case PRIMENET_ERROR_INVALID_ASSIGNMENT_KEY:
			resmsg = "Invalid assignment key";
			break;
		case PRIMENET_ERROR_INVALID_ASSIGNMENT_TYPE:
			resmsg = "Invalid assignment type";
			break;
		case PRIMENET_ERROR_INVALID_RESULT_TYPE:
			resmsg = "Invalid result type";
			break;
		default:
			resmsg = "Unknown error code";
			break;
		}

/* Print out the error code, text, and details */
		
		sprintf (buf, "PrimeNet error %d: %s\n", res, resmsg);
		LogMsg (buf);
		sprintf (buf, "%s\n", errtxt);
		LogMsg (buf);
	}

/* If there was no error code but there was some error text, then print */
/* the error text. */

	else if (strcmp (errtxt, "SUCCESS")) {
		LogMsg ("PrimeNet success code with additional info:\n");
		sprintf (buf, "%s\n", errtxt);
		LogMsg (buf);
	}

/* Parse remaining response parameters */

	switch (operation) {
	case PRIMENET_UPDATE_COMPUTER_INFO:	/* update computer info */
		{
		struct primenetUpdateComputerInfo *z;

		z = (struct primenetUpdateComputerInfo *) pkt;
		parse_string (s, "g", z->computer_guid, sizeof (z->computer_guid));
		parse_string (s, "u", z->user_id, sizeof (z->user_id));
		parse_string (s, "un", z->user_name, sizeof (z->user_name));
		parse_string (s, "cn", z->computer_name, sizeof (z->computer_name));
		parse_uint (s, "od", &z->options_counter);

		if (!strcmp (z->user_id, "admin_user_anon"))
			strcpy (z->user_id, "ANONYMOUS");

		break;
		}
	case PRIMENET_PROGRAM_OPTIONS:
		{
		struct primenetProgramOptions *z;

		z = (struct primenetProgramOptions *) pkt;
		z->num_workers = -1;
		parse_int (s, "nw", &z->num_workers);
		z->work_preference = -1;
		parse_int (s, "w", &z->work_preference);
		z->priority = -1;
		parse_int (s, "Priority", &z->priority);
		z->daysOfWork = -1;
		parse_int (s, "DaysOfWork", &z->daysOfWork);
		z->dayMemory = -1;
		parse_int (s, "DayMemory", &z->dayMemory);
		z->nightMemory = -1;
		parse_int (s, "NightMemory", &z->nightMemory);
		z->dayStartTime = -1;
		parse_int (s, "DayStartTime", &z->dayStartTime);
		z->nightStartTime = -1;
		parse_int (s, "NightStartTime", &z->nightStartTime);
		z->runOnBattery = -1;
		parse_int (s, "RunOnBattery", &z->runOnBattery);
		parse_uint (s, "od", &z->options_counter);
		break;
		}
	case PRIMENET_GET_ASSIGNMENT:
		{
		struct primenetGetAssignment *z;

		z = (struct primenetGetAssignment *) pkt;
		parse_string (s, "k", z->assignment_uid, sizeof (z->assignment_uid));
		parse_uint (s, "w", &z->work_type);
		parse_double (s, "A", &z->k);
		parse_uint (s, "b", &z->b);
		parse_uint (s, "n", &z->n);
		parse_int (s, "c", &z->c);
		parse_uint (s, "p1", &z->has_been_pminus1ed);
		parse_double (s, "sf", &z->how_far_factored);
		parse_double (s, "ef", &z->factor_to);
		parse_double (s, "B1", &z->B1);
		parse_double (s, "B2", &z->B2);
		parse_uint (s, "CR", &z->curves);
		parse_double (s, "saved", &z->tests_saved);
		parse_string (s, "kf", z->known_factors, sizeof (z->known_factors));
		break;
		}
	case PRIMENET_REGISTER_ASSIGNMENT:
		{
		struct primenetRegisterAssignment *z;

		z = (struct primenetRegisterAssignment *) pkt;
		parse_string (s, "k", z->assignment_uid, sizeof (z->assignment_uid));
		break;
		}
	case PRIMENET_ASSIGNMENT_PROGRESS:
		break;
	case PRIMENET_ASSIGNMENT_RESULT:
		break;
	case PRIMENET_ASSIGNMENT_UNRESERVE:
		break;
	case PRIMENET_BENCHMARK_DATA:
		break;
	case PRIMENET_PING_SERVER:
		{
		struct primenetPingServer *z;

		z = (struct primenetPingServer *) pkt;
		parse_string (s, "r", z->ping_response, sizeof (z->ping_response));
		break;
		}
	}
	return (res);
}

/*
// LoadPrimenet: call from main thread at startup
*/

void LoadPrimenet ()
{

/* The cURL documentation strongly recommends calling curl_global_init */
/* from the main thread rather than letting curl_easy_init do the */
/* initialization when other threads are running. */

	if (IniSectionGetInt (INI_FILE, iniSection, "UseCURL", 1))
		curl_global_init (CURL_GLOBAL_ALL);
}


/*
// Primenet: main interface to Prime95.exe
*/

int PRIMENET (short operation, void *pkt)
{
	int	status;
	char args[4096];		/* formatted arguments buffer */
	char pbuf[4096];		/* return page buffer */

/* Check if we should use old Primenet server */

	if (USE_V4) return (v4_PRIMENET (operation, pkt));

/* Assemble GET/POST arguments */

	status = format_args (args, operation, pkt);
	if (status) return (status);

/* Send arguments, read back resulting page */

	if (IniSectionGetInt (INI_FILE, iniSection, "UseCURL", 1))
		status = pnHttpServerCURL (pbuf, sizeof (pbuf), args);
	else
		status = pnHttpServer (pbuf, sizeof (pbuf), args);
	if (status) return (status);

/* Extract results from returned page into packet */

	return (parse_page (pbuf, operation, pkt));
}










#include "security.c"
#ifndef IPSHASH
#define v4_hash_packet(a,b)
#endif

/*
// Primenet: main interface to Prime95.exe, accepts
//           and returns PrimeNet 3.0 packets
*/

int v4_PRIMENET (short operation, void *pkt)
{
	int	status, status2;
	char args[4096];		/* formatted arguments buffer */
	char pbuf[4096];		/* return page buffer */

/* Assemble GET/POST arguments */

	status = v4_format_args (args, operation, pkt);
	if (status == 999) { status = 0; goto skipit; }
	if (status) return (status);

/* Send arguments, read back resulting page */

	if (IniSectionGetInt (INI_FILE, iniSection, "UseCURL", 1))
		status = pnHttpServerCURL (pbuf, sizeof (pbuf), args);
	else
		status = pnHttpServer (pbuf, sizeof (pbuf), args);
	if (status) return (status);

/* Extract results from returned page into packet */

	status = v4_parse_page (pbuf, operation, pkt);

/* Some v5 messages require two v4 messages */

skipit:
	if (operation != PRIMENET_ASSIGNMENT_RESULT &&
	    operation != PRIMENET_UPDATE_COMPUTER_INFO)
		return (status);

/* Assemble GET/POST arguments */

	status2 = v4_format_args (args, operation + 1000, pkt);
	if (status2) return (status2);

/* Send arguments, read back resulting page */

	if (IniSectionGetInt (INI_FILE, iniSection, "UseCURL", 1))
		status2 = pnHttpServerCURL (pbuf, sizeof (pbuf), args);
	else
		status2 = pnHttpServer (pbuf, sizeof (pbuf), args);
	if (status2) return (status2);

/* Extract results from returned page into packet */

	status2 = v4_parse_page (pbuf, operation + 1000, pkt);
	if (status) return (status);
	return (status2);
}

/* Utility routines */

unsigned long exp_from_aid (char id[33])
{
	unsigned long x;

	sscanf (id+24, "%x", &x);
	return (x);
}

short type_from_aid (char id[33])
{
	char	buf[9];
	unsigned long x;

	memcpy (buf, id+16, 8); buf[8] = 0;
	sscanf (buf, "%x", &x);
	return ((short) x);
}

void gen_aid (char id[33], short type, unsigned long exp)
{
	sprintf (id, "0000000000000000%08X%08X", (long) type, exp);
}

int faclim (unsigned long n)
{
	return ((n > FAC71) ? 71 :
		(n > FAC70) ? 70 :
		(n > FAC69) ? 69 :
		(n > FAC68) ? 68 :
		(n > FAC67) ? 67 :
		(n > FAC66) ? 66 : 65);
}

/*
// format_args: format a HTTP argument string from a PrimeNet v4 packet
*/

int v4_format_args (char* args, short operation, void* pkt)
{
	char	V4_COMPID[40];

	strcpy (V4_COMPID, COMPID);	// Truncate COMPID to 12 chars
	V4_COMPID[12] = 0;

	*args = 0;

	switch (operation) {
	case PRIMENET_UPDATE_COMPUTER_INFO:	/* update computer info */
		{
		struct primenetUpdateComputerInfo *z;
		struct {
			short	versionNumber;
			short	structSize;
			short	hashFunction;
			unsigned short hash;
			unsigned short salt;
			char	userID[15];
			char	userPW[9];
			char	computerID[13];
			char	userName[80];
			char	userEmailAddr[80];
			char	oldUserID[15];
			char	oldUserPW[9];
			char	bUserOptions;
		} v4pkt;

		z = (struct primenetUpdateComputerInfo *) pkt;

		memset (&v4pkt, 0, sizeof (v4pkt));
		v4pkt.versionNumber = 4;
		v4pkt.structSize = sizeof (v4pkt);
		v4pkt.hashFunction = 4;
		strcpy (v4pkt.userID, V4_USERID);
		strcpy (v4pkt.userPW, V4_USERPWD);
		strcpy (v4pkt.computerID, V4_COMPID);
		strcpy (v4pkt.userName,	V4_USERNAME);
		IniGetString (INI_FILE, "UserEmailAddr", v4pkt.userEmailAddr,
			sizeof (v4pkt.userEmailAddr), NULL);

		v4_hash_packet (0, &v4pkt);

		sprintf (args,"uu&4&%d&4&%u&%u&%s&%s&%s",
			 v4pkt.structSize, v4pkt.hash, v4pkt.salt,
			 v4pkt.userID, v4pkt.userPW, v4pkt.computerID);
		args += strlen (args);
		*args++ = '&';
		args = armor (args, v4pkt.userName);
		*args++ = '&';
		args = armor (args, v4pkt.userEmailAddr);
		sprintf (args,"&%s&%s&%d",
			 v4pkt.oldUserID, v4pkt.oldUserPW, v4pkt.bUserOptions);
		args += strlen (args);
		break;
		}
	case PRIMENET_UPDATE_COMPUTER_INFO+1000:  /* update computer info */
		{
		struct primenetUpdateComputerInfo *z;
		struct {
			short	versionNumber;
			short	structSize;
			short	hashFunction;
			unsigned short hash;
			unsigned short salt;
			char	userID[15];
			char	userPW[9];
			char	computerID[13];
			short	cpu_type;
			short	speed;
			short	hours_per_day;
		} v4pkt;

		z = (struct primenetUpdateComputerInfo *) pkt;

		memset (&v4pkt, 0, sizeof (v4pkt));
		v4pkt.versionNumber = 4;
		v4pkt.structSize = 54;
		v4pkt.hashFunction = 4;
		strcpy (v4pkt.userID, V4_USERID);
		strcpy (v4pkt.userPW, V4_USERPWD);
		strcpy (v4pkt.computerID, V4_COMPID);
		if (strstr (z->cpu_model, "Intel") != NULL)
			v4pkt.cpu_type = 12;	// Pentium 4
		else
			v4pkt.cpu_type = 11;	// AMD Athlon
		v4pkt.speed = z->cpu_speed;
		v4pkt.hours_per_day = CPU_HOURS;

		v4_hash_packet (6, &v4pkt);

		sprintf (args,"mi&4&54&4&%u&%u&%s&%s&%s&%d&%d&%d",
			 v4pkt.hash, v4pkt.salt, v4pkt.userID,
			 v4pkt.userPW, v4pkt.computerID, v4pkt.cpu_type,
			 v4pkt.speed, v4pkt.hours_per_day);
		args += strlen (args);
		break;
		}
	case PRIMENET_PROGRAM_OPTIONS:
		{
		return (999);
		}
	case PRIMENET_REGISTER_ASSIGNMENT:	/* register assignment */
		{
		struct primenetRegisterAssignment *z;
		int	work_type;

		z = (struct primenetRegisterAssignment *) pkt;

		switch (z->work_type) {
		case PRIMENET_WORK_TYPE_FACTOR:
			work_type = 1;  // V4_PRIMENET_ASSIGN_FACTOR
			break;
		case PRIMENET_WORK_TYPE_DBLCHK:
			work_type = 4;  // V4_PRIMENET_ASSIGN_DBLCHK
			break;
		case PRIMENET_WORK_TYPE_FIRST_LL:
			work_type = 2;  // V4_PRIMENET_ASSIGN_TEST
			break;
		}
		gen_aid (z->assignment_uid, work_type, z->n);

		return (999);
		}
	case PRIMENET_GET_ASSIGNMENT:		/* get assignment */
		{
		struct primenetGetAssignment *z;
		struct {
			short	versionNumber;
			short	structSize;
			short	hashFunction;
			unsigned short hash;
			unsigned short salt;
			char	userID[15];
			char	userPW[9];
			char	computerID[13];
			short	requestType;
			short	programType;
			unsigned long exponent;
			double	how_far_factored;
		} v4pkt;

		z = (struct primenetGetAssignment *) pkt;

		memset (&v4pkt, 0, sizeof (v4pkt));
		v4pkt.versionNumber = 4;
		v4pkt.structSize = 64;
		v4pkt.hashFunction = 4;
		strcpy (v4pkt.userID, V4_USERID);
		strcpy (v4pkt.userPW, V4_USERPWD);
		strcpy (v4pkt.computerID, V4_COMPID);
		v4pkt.programType = 2;
		v4pkt.how_far_factored = 64.0;
		switch (WORK_PREFERENCE[z->cpu_num]) {
		case PRIMENET_WP_WHATEVER:
			v4pkt.requestType =
				(CPU_SPEED < 500) ? 1 :
				(CPU_SPEED < 1500) ? 4 : 2;
			break;
		case PRIMENET_WP_FACTOR_LMH:
		case PRIMENET_WP_FACTOR:
		case PRIMENET_WP_PMINUS1:
		case PRIMENET_WP_ECM_SMALL:
		case PRIMENET_WP_ECM_FERMAT:
		case PRIMENET_WP_ECM_CUNNINGHAM:
		default:
			v4pkt.requestType = 1;  // V4_PRIMENET_ASSIGN_FACTOR
			break;
		case PRIMENET_WP_PFACTOR:
		case PRIMENET_WP_LL_DBLCHK:
			v4pkt.requestType = 4;  // V4_PRIMENET_ASSIGN_DBLCHK
			break;
		case PRIMENET_WP_LL_FIRST:
		case PRIMENET_WP_LL_FIRST_NOFAC:
		case PRIMENET_WP_PRP:
			v4pkt.requestType = 2;  // V4_PRIMENET_ASSIGN_TEST
			break;
		case PRIMENET_WP_WORLD_RECORD:
		case PRIMENET_WP_LL_10M:
		case PRIMENET_WP_LL_100M:
			v4pkt.requestType = 2;  // V4_PRIMENET_ASSIGN_TEST
			v4pkt.requestType |= 16;  // V4_PRIMENET_ASSIGN_BIGONES
			break;
		}
		v4_hash_packet (1, &v4pkt);

		sprintf (args,"ga&4&64&4&%u&%u&%s&%s&%s&%d&%d&%04.1f",
			 v4pkt.hash, v4pkt.salt, v4pkt.userID, v4pkt.userPW,
			 v4pkt.computerID, v4pkt.requestType,
			 v4pkt.programType, v4pkt.how_far_factored);
		args += strlen (args);
		break;
		}
	case PRIMENET_ASSIGNMENT_PROGRESS:
		{
		struct primenetAssignmentProgress *z;
		struct {
			short	versionNumber;
			short	structSize;
			short	hashFunction;
			unsigned short hash;
			unsigned short salt;
			char	userID[15];
			char	userPW[9];
			char	computerID[13];
			unsigned long exponent;
			unsigned long days;
			short	requestType;
			short	programType;
			unsigned long iteration;
			unsigned long nextMsg;
		} v4pkt;

		z = (struct primenetAssignmentProgress *) pkt;

		memset (&v4pkt, 0, sizeof (v4pkt));
		v4pkt.versionNumber = 4;
		v4pkt.structSize = 68;
		v4pkt.hashFunction = 4;
		strcpy (v4pkt.userID, V4_USERID);
		strcpy (v4pkt.userPW, V4_USERPWD);
		strcpy (v4pkt.computerID, V4_COMPID);
		v4pkt.exponent = exp_from_aid (z->assignment_uid);
		v4pkt.days = z->end_date / 86400;
		v4pkt.requestType = type_from_aid (z->assignment_uid);
		v4pkt.programType = 2;
		v4pkt.iteration = (unsigned long) (z->pct_complete / 100.0 * v4pkt.exponent);
		v4pkt.nextMsg = z->next_update / 86400;

		v4_hash_packet (2, &v4pkt);

		sprintf (args,"cd&4&68&4&%u&%u&%s&%s&%s&%ld&%ld&%d&%d&%ld&%ld",
			 v4pkt.hash, v4pkt.salt, v4pkt.userID, v4pkt.userPW,
			 v4pkt.computerID, v4pkt.exponent, v4pkt.days,
			 v4pkt.requestType, v4pkt.programType,
			 v4pkt.iteration, v4pkt.nextMsg);
		args += strlen (args);
		break;
		}
	case PRIMENET_ASSIGNMENT_RESULT:
		{
		struct primenetAssignmentResult *z;
		struct {
			short	versionNumber;
			short	structSize;
			short	hashFunction;
			unsigned short hash;
			unsigned short salt;
			char	userID[15];
			char	userPW[9];
			char	computerID[13];
			char	pad1[1];
			unsigned long exponent;
			short	resultType;
			char	pad2[2];
			union resinf {
				double	how_far_factored;
				char	residue[17];
				char	factor[32];
			} r;
		} v4pkt;

		z = (struct primenetAssignmentResult *) pkt;

		if (z->result_type == 0 ||
		    z->result_type == PRIMENET_AR_ECM_NOFACTOR ||
		    z->result_type == PRIMENET_AR_PRP_RESULT ||
		    z->result_type == PRIMENET_AR_PRP_PRIME ||
		    !z->assignment_uid[0])
			return (999);
		if (z->result_type == PRIMENET_AR_TF_NOFACTOR &&
		    z->end_bits != faclim (z->n))
			return (999);

		memset (&v4pkt, 0, sizeof (v4pkt));
		v4pkt.versionNumber = 4;
		v4pkt.structSize = 88;
		v4pkt.hashFunction = 4;
		strcpy (v4pkt.userID, V4_USERID);
		strcpy (v4pkt.userPW, V4_USERPWD);
		strcpy (v4pkt.computerID, V4_COMPID);
		v4pkt.exponent = exp_from_aid (z->assignment_uid);

		if (z->result_type == PRIMENET_AR_LL_RESULT) {
			v4pkt.resultType = 1;
			strcpy (v4pkt.r.residue, z->residue);
		}
		if (z->result_type == PRIMENET_AR_LL_PRIME) {
			v4pkt.resultType = 3;
			strcpy (v4pkt.r.residue, "0000000000000000");
		}
		if (z->result_type == PRIMENET_AR_TF_FACTOR ||
		    z->result_type == PRIMENET_AR_P1_FACTOR ||
		    z->result_type == PRIMENET_AR_ECM_FACTOR) {
			v4pkt.resultType = 2;
			if (strlen (z->factor) <= 31)
				strcpy (v4pkt.r.factor, z->factor);
			else
				memcpy (v4pkt.r.factor, z->factor, 31);
		}
		if (z->result_type == PRIMENET_AR_TF_NOFACTOR) {
			v4pkt.resultType = 0;
			v4pkt.r.how_far_factored = z->end_bits;
		}
		if (z->result_type == PRIMENET_AR_P1_NOFACTOR) {
			v4pkt.resultType = 0;
			v4pkt.r.how_far_factored = 74.5;
		}
		v4_hash_packet (4, &v4pkt);

		if (z->result_type == PRIMENET_AR_TF_NOFACTOR ||
		    z->result_type == PRIMENET_AR_P1_NOFACTOR)
			sprintf (args, "ar&4&88&4&%u&%u&%s&%s&%s&%ld&0&%04.1f",
				 v4pkt.hash, v4pkt.salt, v4pkt.userID,
				 v4pkt.userPW, v4pkt.computerID, v4pkt.exponent,
				 v4pkt.r.how_far_factored);
		else
			sprintf (args, "ar&4&88&4&%u&%u&%s&%s&%s&%ld&%d&%s",
				 v4pkt.hash, v4pkt.salt, v4pkt.userID,
				 v4pkt.userPW, v4pkt.computerID, v4pkt.exponent,
				 v4pkt.resultType, v4pkt.r.factor);
		args += strlen (args);
		break;
		}
	case PRIMENET_ASSIGNMENT_RESULT+1000:
		{
		struct primenetAssignmentResult *z;
		struct {
			short	versionNumber;
			short	structSize;
			short	hashFunction;
			unsigned short hash;
			unsigned short salt;
			char	userID[15];
			char	userPW[9];
			char	computerID[13];
			char	message[200];
			char	pad1[1];
		} v4pkt;

		z = (struct primenetAssignmentResult *) pkt;

		memset (&v4pkt, 0, sizeof (v4pkt));
		v4pkt.versionNumber = 4;
		v4pkt.hashFunction = 4;
		strcpy (v4pkt.userID, V4_USERID);
		strcpy (v4pkt.userPW, V4_USERPWD);
		strcpy (v4pkt.computerID, V4_COMPID);
		strcpy (v4pkt.message, z->message);
		v4pkt.structSize = (short) (sizeof (v4pkt) - 200 + strlen (v4pkt.message) + 1);

		v4_hash_packet (3, &v4pkt);

		sprintf (args,"rm&4&%d&4&%u&%u&%s&%s&%s",
			 v4pkt.structSize, v4pkt.hash, v4pkt.salt,
			 v4pkt.userID, v4pkt.userPW, v4pkt.computerID);
		args += strlen (args);
		*args++ = '&';
		args = armor (args, v4pkt.message);
		break;
		}
	case PRIMENET_ASSIGNMENT_UNRESERVE:	/* assignment unreserve */
		{
		struct primenetAssignmentUnreserve *z;
		struct {
			short	versionNumber;
			short	structSize;
			short	hashFunction;
			unsigned short hash;
			unsigned short salt;
			char	userID[15];
			char	userPW[9];
			char	computerID[13];
			char	pad1[1];
			unsigned long exponent;
			short	resultType;
			char	pad2[2];
			char	factor[32];
		} v4pkt;

		z = (struct primenetAssignmentUnreserve *) pkt;

		memset (&v4pkt, 0, sizeof (v4pkt));
		v4pkt.versionNumber = 4;
		v4pkt.structSize = 88;
		v4pkt.hashFunction = 4;
		strcpy (v4pkt.userID, V4_USERID);
		strcpy (v4pkt.userPW, V4_USERPWD);
		strcpy (v4pkt.computerID, V4_COMPID);
		v4pkt.exponent = exp_from_aid (z->assignment_uid);
		v4pkt.resultType = 4;

		v4_hash_packet (4, &v4pkt);

		sprintf (args, "ar&4&88&4&%u&%u&%s&%s&%s&%ld&4",
			 v4pkt.hash, v4pkt.salt, v4pkt.userID,
			 v4pkt.userPW, v4pkt.computerID, v4pkt.exponent);
		args += strlen (args);
		break;
		}
	case PRIMENET_BENCHMARK_DATA:
		{
		return (999);
		break;
		}
	case PRIMENET_PING_SERVER:
		{
		struct primenetPingServer *z;

		z = (struct primenetPingServer *) pkt;
		sprintf (args,"ps&4&.&.");
		args += strlen (args);
		break;
		}
	}
	return (0);
}

/*
// parse_page: reads the server response page tokens and values
//             and converts these back into a C structure
*/

int v4_parse_page (char *response_buf, short operation, void *pkt)
{
	char	*s;
	int32_t	res;

/* get result code, which is always first */

	s = response_buf;
	if (operation == PRIMENET_REGISTER_ASSIGNMENT)
		res = 0;
	else if (!parse_int (s, "pnResult", &res)) {
		LogMsg ("PnErrorResult value missing\n");
		return (PRIMENET_ERROR_PNERRORRESULT);
	}

/* Ignore error codes that cause much consternation in v24 */

	if (res == 11) res = 0;  // exponent already tested
	if (res == 17) res = 0;  // exponent assigned to someone else

/* If result is non-zero print out an error message */

	if (res) {
		char	buf[2000];
		char	*resmsg;

		if (res == 23) res = PRIMENET_ERROR_SERVER_BUSY;
		else if (res == 5) res = PRIMENET_ERROR_ACCESS_DENIED;
		else res += 20000;

/* Convert the error number to text */

		switch (res) {
		case PRIMENET_ERROR_SERVER_BUSY:
			resmsg = "Server busy";
			break;
		case PRIMENET_ERROR_ACCESS_DENIED:
			resmsg = "Access denied";
			break;
		default:
			resmsg = "Unknown error code";
			break;
		}

/* Print out the error code, text, and details */

		sprintf (buf, "PrimeNet error %d: %s\n", res, resmsg);
		LogMsg (buf);
	}

/* Generate a new user id on access denied error */

	if (res == PRIMENET_ERROR_ACCESS_DENIED) {
		LogMsg ("A new userID and password will be generated.\n");
		V4_USERID[0] = 0;
		V4_USERPWD[0] = 0;
		spoolMessage (PRIMENET_UPDATE_COMPUTER_INFO, NULL);
	}

/* Parse remaining response parameters */

	switch (operation) {
	case PRIMENET_UPDATE_COMPUTER_INFO:	/* update computer info */
		{
		struct primenetUpdateComputerInfo *z;

		z = (struct primenetUpdateComputerInfo *) pkt;
		parse_string (s, "g", z->computer_guid, sizeof (z->computer_guid));
		parse_string (s, "u", z->user_id, sizeof (z->user_id));
		parse_string (s, "un", z->user_name, sizeof (z->user_name));
		parse_string (s, "cn", z->computer_name, sizeof (z->computer_name));
		z->options_counter = 
			IniGetInt (LOCALINI_FILE, "SrvrP00", 1);

		if (V4_USERID[0] == 0) {
			parse_string (s, "userID", V4_USERID, sizeof (V4_USERID));
			parse_string (s, "userPW", V4_USERPWD, sizeof (V4_USERPWD));
			if (V4_USERID[strlen(V4_USERID)-1] == '\n')
				V4_USERID[strlen(V4_USERID)-1] = 0;
			if (V4_USERPWD[strlen(V4_USERPWD)-1] == '\n')
				V4_USERPWD[strlen(V4_USERPWD)-1] = 0;
			IniWriteString (INI_FILE, "UserID", V4_USERID);
			IniWriteString (INI_FILE, "UserPWD", V4_USERPWD);
			strcpy (z->user_id, V4_USERID);
		}

		break;
		}
	case PRIMENET_PROGRAM_OPTIONS:
		{
		struct primenetProgramOptions *z;

		z = (struct primenetProgramOptions *) pkt;
		z->num_workers = -1;
		z->work_preference = -1;
		z->priority = -1;
		z->daysOfWork = -1;
		z->dayMemory = -1;
		z->nightMemory = -1;
		z->dayStartTime = -1;
		z->nightStartTime = -1;
		z->runOnBattery = -1;
		z->options_counter = IniGetInt (LOCALINI_FILE, "SrvrP00", 1);
		break;
		}
	case PRIMENET_GET_ASSIGNMENT:
		{
		struct primenetGetAssignment *z;
		uint32_t v4_work_type;

		z = (struct primenetGetAssignment *) pkt;
		parse_uint (s, "requestType", &v4_work_type);
		parse_uint (s, "exponent", &z->n);
		parse_double (s, "how_factored", &z->how_far_factored);

		gen_aid (z->assignment_uid, z->work_type, z->n);

		switch (v4_work_type) {
		case 1:  // V4_PRIMENET_ASSIGN_FACTOR
			z->work_type = PRIMENET_WORK_TYPE_FACTOR;
			break;
		case 4:  // V4_PRIMENET_ASSIGN_DBLCHK
			z->work_type = PRIMENET_WORK_TYPE_DBLCHK;
			break;
		case 2:  // V4_PRIMENET_ASSIGN_TEST
			z->work_type = PRIMENET_WORK_TYPE_FIRST_LL;
			break;
		}

		z->k = 1.0;
		z->b = 2;
		z->c = -1;

		if (z->how_far_factored != (double) (int) z->how_far_factored) {
			z->how_far_factored = z->how_far_factored - 0.5;
			z->has_been_pminus1ed = 1;
		} else
			z->has_been_pminus1ed = 0;
		z->factor_to = faclim (z->n);
		z->B1 = 0.0;
		z->B2 = 0.0;
		z->curves = 0;
		z->tests_saved = 0;
		z->known_factors[0] = 0; 
		break;
		}
	case PRIMENET_REGISTER_ASSIGNMENT:
		break;
	case PRIMENET_ASSIGNMENT_PROGRESS:
		break;
	case PRIMENET_ASSIGNMENT_RESULT:
		break;
	case PRIMENET_ASSIGNMENT_UNRESERVE:
		break;
	case PRIMENET_BENCHMARK_DATA:
		break;
	case PRIMENET_PING_SERVER:
		{
		struct primenetPingServer *z;

		z = (struct primenetPingServer *) pkt;
		parse_string (s, "primenetServerName", z->ping_response, sizeof (z->ping_response));
		break;
		}
	}
	return (res);
}
