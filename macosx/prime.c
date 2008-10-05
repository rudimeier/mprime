/* Copyright 1995-2008 Mersenne Research, Inc. */
/* Author:  George Woltman */
/* Email: woltman@alum.mit.edu */

/* Include files needed by all ports */
#include "prime.h"
#include <ctype.h>
#include <fcntl.h>
#include <math.h>
#include <memory.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

/* Required Linux files */
#ifdef __linux__
#include <dirent.h>
#include <unistd.h>
#include <linux/unistd.h>
#include <asm/unistd.h>
#include <sys/resource.h>
#include <sys/sysinfo.h>
#include <sys/time.h>
#include <sys/timeb.h>
#endif

/* Required Mac OS X files */
#ifdef __APPLE__
#include <dirent.h>
#include <pthread.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/resource.h>
#include <sys/sysctl.h>
#include <sys/time.h>
#include <sys/timeb.h>
#define PTHREAD_MIN_PRIORITY 0		/* Missing #defines from pthreads.h */
#define PTHREAD_MAX_PRIORITY 31		/* Missing #defines from pthreads.h */
#endif

/* Required FreeBSD files */
#ifdef __FreeBSD__
#include <dirent.h>
#include <pthread.h>
#include <unistd.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/timeb.h>
#endif

/* Required OS/2 header files */
#ifdef __IBMC__
#define INCL_DOS
#define INCL_DOSPROFILE
#include <os2.h>
#include <direct.h>
#include <io.h>
#include <process.h>
#include <sys/timeb.h>
typedef int pid_t;
#include "dosqss.h"
#endif

/* Required Watcom C (is this defunct DOS port???) header files */
#ifdef __WATCOMC__
#include <dirent.h>
#include <unistd.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <sys/timeb.h>
#endif

/* Globals */

#ifndef __WATCOMC__
#ifndef __APPLE__
#define OPEN_MAX 20
#endif
#endif

#ifdef MPRIME_LOADAVG
#define LINUX_LDAV_FILE "/proc/loadavg"
int volatile SLEEP_STOP = 0;
long LOAD_CHECK_TIME = 0;
double HI_LOAD = 0.0;
double LO_LOAD = 0.0;
#endif

int volatile THREAD_KILL = 0;
int NO_GUI = 1;
int VERBOSE = 0;
int MENUING = 0;

/* Common code */

#include "gwutil.h"
#include "commona.c"
#include "commonb.c"
#include "commonc.c"
#include "ecm.c"
#include "primenet.c"
#include "gwtest.c"

/* Signal handlers */

void sigterm_handler(int signo)
{
	stop_workers_for_escape ();
	if (signo != SIGINT) THREAD_KILL = TRUE;
	(void)signal(signo, sigterm_handler);
}

#ifdef MPRIME_LOADAVG

/* Routine to get the current load average */
double get_load_average (void)
{
#ifdef __linux__
	char	ldavgbuf[40];
	double	load_avg;
	int	fd, count;

	fd = open (LINUX_LDAV_FILE, O_RDONLY);
	if (fd == -1) return (-1.0);
	count = read (fd, ldavgbuf, 40);
	(void) close (fd);
	if (count <= 0) return (-1.0);
	count = sscanf (ldavgbuf, "%lf", &load_avg);
	if (count < 1) return (-1.0);
	return (load_avg);
#endif
#if defined (__FreeBSD__) || defined (__APPLE__)
	double load[3];

	if (getloadavg (load, sizeof(load)/sizeof(load[0])) < 0) return (-1.0);
	return (load[0]);
#endif
}

/* load_handler: call by signal routine,
   sets SLEEP_STOP to TRUE if load is too high */
void load_handler (
	int	sig)
{
	double  load_avg;

	load_avg = get_load_average ();
	if (load_avg < 0.0) return;
  
	if (SLEEP_STOP) {
		if (load_avg < LO_LOAD)
			SLEEP_STOP = FALSE;
	} else {
		if (load_avg > HI_LOAD)
			SLEEP_STOP = TRUE;
	}
}

/* init_load_check: initialises timer that calls load_handler
   every LOAD_CHECK_TIME seconds */
void init_load_check (void)
{
	struct itimerval timer, otimer;
	struct sigaction sigact;
	int	ret;

	timer.it_interval.tv_sec  =  LOAD_CHECK_TIME;
	timer.it_interval.tv_usec =  0;
	timer.it_value.tv_sec     =  LOAD_CHECK_TIME;
	timer.it_value.tv_usec    =  0;

	ret = setitimer (ITIMER_REAL, &timer, &otimer);
	if (ret < 0) return;
  
	sigact.sa_handler = &load_handler;
	sigemptyset(&sigact.sa_mask);
	sigact.sa_flags =  SA_RESTART;
	ret = sigaction(SIGALRM, &sigact, NULL);
	if (ret < 0) { /* clean up after ourselves */
		setitimer (ITIMER_REAL, &otimer, NULL);
	}
}

/* test_sleep: tests if SLEEP_STOP is set and sleeps until load is normal
   again or WORKER_THREADS_STOPPING is set
*/
void test_sleep (void) 
{
	sigset_t newmask;

	while (SLEEP_STOP && !WORKER_THREADS_STOPPING) {
		sigemptyset (&newmask);
		sigsuspend (&newmask);
	}
}
#endif

/* Main entry point! */

int main (
	int	argc,
	char	*argv[])
{
	char	buf[256];
	int	named_ini_files = -1;
	int	contact_server = 0;
	int	torture_test = 0;
	int	i;
	char	*p;

/* catch termination signals */

	(void)signal(SIGTERM, sigterm_handler);
	(void)signal(SIGINT, sigterm_handler);

/* No buffering of output */

	setvbuf (stdout, NULL, _IONBF, 0);

/* Change to the executable's directory */
/* NOTE:  This only change's the working directory if the user typed */
/* in a full path to the executable (as opposed to finding it on the PATH) */

	strcpy (buf, argv[0]);
	p = strrchr (buf, '/');
	if (p != NULL) {
		*p = 0;
		_chdir (buf);
	}

/* Initialize gwnum call back routines.  Using callback routines lets the */
/* gwnum library have a nice clean interface for users that do not need */
/* additional functionality that only prime95 uses. */

	StopCheckRoutine = stopCheck;
	OutputBothRoutine = OutputBoth;

/* Process command line switches */

	for (i = 1; i < argc; i++) {
		p = argv[i];

		if (*p++ != '-') break;
		switch (*p++) {

/* Accept a -A switch indicating an alternate set of INI files */
/* are to be used. */

		case 'A':
		case 'a':
			named_ini_files = 0;
			while (isspace (*p)) p++;
			while (isdigit (*p)) {
				named_ini_files = named_ini_files * 10 + (*p - '0');
				p++;
			}
			break;

/* -C - contact the server now, then exit */

		case 'C':
		case 'c':
			contact_server = 1;
			VERBOSE = TRUE;
			NO_GUI = FALSE;
			break;
			
/* -D - debug */

		case 'D':
		case 'd':
			VERBOSE = TRUE;
			NO_GUI = FALSE;
			break;

/* -H - help */

		case 'H':
		case 'h':
		case '?':
			goto usage;

/* -M - Menu */

		case 'M':
		case 'm':
			MENUING = 1;
			NO_GUI = FALSE;
			break;

/* -S - status */

		case 'S':
		case 's':
			MENUING = 2;
			NO_GUI = FALSE;
			break;
		  

/* -T - Torture test */

		case 'T':
		case 't':
			torture_test = TRUE;
			break;

/* -V - version number */

		case 'V':
		case 'v':
			printf ("Mersenne Prime Test Program, Version %s.%d\n", VERSION, PORT);
			return (0); 

/* -W - use a different working directory */

		case 'W':
		case 'w':
			_chdir (p);
			break; 

/* Otherwise unknown switch */

		default:
			printf ("Invalid switch\n");
			goto usage;
		}
	}

/* Determine the names of the INI files, read them, do other initialization. */

	nameAndReadIniFiles (named_ini_files);

/* Read load averaging settings from INI files */

#ifdef MPRIME_LOADAVG
	IniGetString (INI_FILE, "MaxLoad", buf, sizeof (buf), "0");
	HI_LOAD = atof (buf);
	IniGetString (INI_FILE, "MinLoad", buf, sizeof (buf), "0");
	LO_LOAD = atof (buf);
	IniGetString (INI_FILE, "PauseTime", buf, sizeof (buf), "0");
	LOAD_CHECK_TIME = atol (buf);

/* Initialise load checking */

	if (HI_LOAD > 0.0 && LOAD_CHECK_TIME > 0)
		init_load_check ();
#endif

/* If running the torture test, do so now. */

	if (torture_test) {
		int	num_threads;

		VERBOSE = TRUE;
		NO_GUI = FALSE;
		num_threads = IniGetInt (INI_FILE, "TortureThreads",
					 NUM_CPUS * CPU_HYPERTHREADS);
		LaunchTortureTest (num_threads, TRUE);
	}

/* If this is a stress tester, then turn on menuing. */

	else if (IniGetInt (INI_FILE, "StressTester", 99) == 1) {
		MENUING = 1;
		VERBOSE = TRUE;
		NO_GUI = FALSE;
		main_menu ();
	}

/* On first run, get user id before contacting server */
/* for a work assignment.  To make first time user more comfortable, we will */
/* display data to the screen, rather than running silently. */

	else if (IniGetInt (INI_FILE, "StressTester", 99) == 99) {
		VERBOSE = TRUE;
		NO_GUI = FALSE;
		test_welcome ();
	}

/* If we are to contact the server, do so now.  This option lets the */
/* user create a batch file that contacts the server at regular intervals */
/* or when the ISP is contacted, etc. */

	else if (contact_server) {
		do_manual_comm_now ();
	}

/* Bring up the main menu */

	else if (MENUING == 1)
		main_menu ();
	else if (MENUING == 2)
		test_status();

/* Continue testing, return when worker threads exit. */

	else {
		linuxContinue ("Another mprime is already running!\n", TRUE);
	}

/* All done */

	return (0);

/* Invalid args message */

usage:	printf ("Usage: mprime [-cdhmstv] [-aN] [-wDIR]\n");
	printf ("-c\tContact the PrimeNet server, then exit.\n");
	printf ("-d\tPrint detailed information to stdout.\n");
	printf ("-h\tPrint this.\n");
	printf ("-m\tMenu to configure mprime.\n");
	printf ("-s\tDisplay status.\n");
	printf ("-t\tRun the torture test.\n");
	printf ("-v\tPrint the version number.\n");
	printf ("-aN\tUse an alternate set of INI and output files.\n");
	printf ("-wDIR\tRun from a different working directory.\n");
	printf ("\n");
	return (1);
}

/* Create an MDI output window for the thread -- unless we created */
/* one earlier and the user has not closed it */

void create_window (
	int	thread_num)
{
}

/* Set the title prefix for this MDI window - only called once */

void base_title (
	int	thread_num,
	char	*str)
{
}

/* Put a title on the MDI window */

void title (int thread_num, char *msg)
{
}

void flashWindowAndBeep (void)
{
	printf ("\007");
}

/* Do some work prior to launching worker threads */
/* Windows uses this to implement boot delay. */

void PreLaunchCallback (
	int	launch_type)
{
}

/* Do some work after worker threads have terminated */

void PostLaunchCallback (
	int	launch_type)
{
}

/* OSes that must poll for whether the ESC key was hit do it here. */
/* We use this opportunity to perform other miscellaneous tasks that */
/* can't be done any other way. */

void stopCheckCallback (
	int	thread_num)
{
#ifdef MPRIME_LOADAVG
	test_sleep ();
#endif
}

void RealOutputStr (int thread_num, char *buf)
{
static	int	last_char_out_was_newline = TRUE;
	if (VERBOSE || MENUING) {
		if (last_char_out_was_newline && !(MERGE_WINDOWS & MERGE_NO_PREFIX)) {
			if (thread_num == MAIN_THREAD_NUM)
				printf ("[Main thread");
			else if (thread_num == COMM_THREAD_NUM)
				printf ("[Comm thread");
			else if (NUM_WORKER_THREADS == 1)
				printf ("[Work thread");
			else
				printf ("[Worker #%d", thread_num+1);
			if (buf[0] == '[')
				printf (" %s", buf+1);
			else
				printf ("] %s", buf);
		} else
			printf ("%s", buf);
		last_char_out_was_newline = (buf[strlen(buf)-1] == '\n');
	}
}

/* Return TRUE if we are on battery power. */

int OnBattery (void)
{
}

/* The current implementation comes courtesy of Tim Wood and Dennis Gregorovic */

unsigned long physical_memory (void)
{
#ifdef __APPLE__
	int	mib[2];
	union {
		uint32_t ui32;
		uint64_t ui64;
	} value;
	size_t	len;

	mib[0] = CTL_HW;
	mib[1] = HW_MEMSIZE;
	len = sizeof (value);
	if (sysctl (mib, 2, &value, &len, NULL, 0) < 0)
		return (1024);		/* On error, guess 1GB */
	if (len == sizeof (uint32_t))
		return (value.ui32 >> 20);
	else
		return ((unsigned long) (value.ui64 >> 20));
#else
        struct sysinfo sys_info;

        if (sysinfo(&sys_info) != 0) return (1024);  /* Guess 1GB */

	return ((unsigned long)
		((double) sys_info.totalram *
		 (double) sys_info.mem_unit / 1048576.0));
#endif
}

unsigned long num_cpus (void)
{
#ifdef __APPLE__
	int	mib[2];
	int	ncpus;
	size_t	len;

	mib[0] = CTL_HW;
	mib[1] = HW_NCPU;
	len = sizeof (ncpus);
	sysctl (mib, 2, &ncpus, &len, NULL, 0);
	return (ncpus);
#else
	FILE	*fd;
	char	buf[80];
	int	count;

	count = 0;
	fd = fopen ("/proc/cpuinfo", "r");
	if (fd != NULL) {
		for ( ; ; ) {
			if (fscanf (fd, "%s", buf) == EOF) break;
			buf[9] = 0;
			if (strcmp (buf, "processor") == 0) count++;
		}
		fclose (fd);
	}
	if (count == 0) count = 1;
	return (count);
#endif
}

/* Return a better guess for amount of memory to use in a torture test. */
/* Caller passes in its guess for amount of memory to use, but this routine */
/* can reduce that guess based on OS-specific code that looks at amount */
/* of available physical memory. */
/* This code was written by an anonymous GIMPS user. */

unsigned long GetSuggestedMemory (unsigned long nDesiredMemory)
{
	return (nDesiredMemory);
}

int getDefaultTimeFormat (void)
{
	return (2);
}

void Sleep (
	long	ms) 
{
#ifdef __IBMC__
	DosSleep(ms);
#else
	usleep (ms * 1000);
#endif
}

/* Clear the array of active thread handles */

void clearThreadHandleArray (void)
{
}

/* Register a thread termination.  We remove the thread handle from the */
/* list of active worker threads. */

void registerThreadTermination (void)
{
}

/* When stopping or exiting we raise the priority of all worker threads */
/* so that they can terminate in a timely fashion even if there are other */
/* CPU bound tasks running. */

void raiseAllWorkerThreadPriority (void)
{
}


/* Set priority.  Map one (prime95's lowest priority) to 20 */
/* (linux's lowest priority).  Map eight (prime95's normal priority) to */
/* 0 (linux's normal priority). */

void setThreadPriorityAndAffinity (
	int	priority,		/* Priority, 1=low, 9=high */
	int	mask)			/* Affinity mask */
{
#ifdef __IBMC__
	DosSetPriority(PRTYS_PROCESS,
		(priority < 6) ? PRTYC_IDLETIME : PRTYC_REGULAR,
		(priority == 1 || priority == 6) ? PRTYD_MINIMUM :
		(priority == 2 || priority == 7) ? -10 :
		(priority == 3 || priority == 8) ? 0 :
		(priority == 4 || priority == 9) ? 10 :
		PRTYD_MAXIMUM,
		0);
#endif
#ifdef __linux__
	int	linux_priority, errcode;
	pid_t	thread_id;

/* I couldn't get the standard syscall0 declaration of gettid to */
/* work in my Linux distro.  Use the direct system call instead. */
	thread_id = (pid_t) syscall (__NR_gettid);

/* Set priority.  Map one (prime95's lowest priority) to 19 */
/* (linux's lowest priority).  Map eight (prime95's normal priority) to */
/* 0 (linux's normal priority). */

	linux_priority = (8 - (int) priority) * 19 / 7;
	errcode = setpriority (PRIO_PROCESS, thread_id, linux_priority);

/* Set affinity for this thread.  We assume the processor mask is the same */
/* as in Windows with physical CPUs representing the least significant bits */
/* and hyperthreaded logical CPUs as the next set of more significant bits */

	errcode = sched_setaffinity (thread_id, sizeof (mask), &mask);
#endif

#if defined (__APPLE__) || defined (__FreeBSD__)
static	int	default_priority = 0;
static	int	default_policy = 0;
	struct sched_param sp;

/* Get the default thread priority when a thread is first launched */
	if (default_priority == 0) {
		memset (&sp, 0, sizeof(struct sched_param));
	        if (pthread_getschedparam (pthread_self(),
					   &default_policy, &sp) >= 0) {
			default_priority = sp.sched_priority;
		} else {
			default_policy = SCHED_RR;
			default_priority = PTHREAD_MIN_PRIORITY +
					   (PTHREAD_MAX_PRIORITY -
					    PTHREAD_MIN_PRIORITY) / 2;
		}
	}

/* Map one (prime95's lowest priority) to PTHREAD_MIN_PRIORITY */
/* (pthread's lowest priority).  Map eight (prime95's normal priority) to */
/* pthread's default priority. */

	memset (&sp, 0, sizeof(struct sched_param));
	sp.sched_priority = PTHREAD_MIN_PRIORITY +
			    (priority - 1) *
				(default_priority - PTHREAD_MIN_PRIORITY) / 7;
	pthread_setschedparam (pthread_self(), default_policy, &sp);
#endif
}

void BlinkIcon (int thread_num, int x)
{
}

void ChangeIcon (int thread_num, int x)
{
}


/* This routine calls primeContinue unless there is another copy of mprime */
/* already running.  In that case, it outputs an optional error message. */

void linuxContinue (
	char	*error_message,
	int	wait_flag)
{
#ifdef __linux__
#define PROCNAME	"/proc/%d/exe"
#endif
#ifdef __FreeBSD__
#define PROCNAME	"/proc/%d/file"
#endif
	pid_t	my_pid, running_pid;
	char	filename[30];
	int	fd;
	struct stat filedata;
	ino_t	inode1, inode2;

/* Compare this process' ID and the pid from the INI file */

	my_pid = getpid ();
	openIniFile (LOCALINI_FILE, 1);
	running_pid = IniGetInt (LOCALINI_FILE, "Pid", 0);
	if (running_pid == 0 || my_pid == running_pid) goto ok;

#ifdef __APPLE__
	goto ok;
#else
#ifdef __OS2__

        {
            USHORT handle1 = 0, handle2 = 0;
            unsigned char buf[0x2000];
            if( !DosQuerySysState(0x01, 0, 0, 0, (PCHAR)buf, 0x2000) ) {
                PQPROCESS p = ((PQTOPLEVEL)buf)->procdata;
                while(p && p->rectype == 1) {
                    if( p->pid == running_pid ) handle1 = p->hndmod;
                    if( p->pid == my_pid ) handle2 = p->hndmod;
                    p = (PQPROCESS)(p->threads + p->threadcnt);
                }
                if( handle1 != handle2 ) goto ok;
            }
        }

#else

/* See if the two pids are running the same executable */

	sprintf (filename, PROCNAME, my_pid);
	fd = _open (filename, _O_RDONLY);
	if (fd < 0) goto ok;
	fstat (fd, &filedata);
	inode1 = filedata.st_ino;
	_close (fd);
	sprintf (filename, PROCNAME, running_pid);
	fd = _open (filename, _O_RDONLY);
	if (fd < 0) goto ok;
	fstat (fd, &filedata);
	inode2 = filedata.st_ino;
	_close (fd);
	if (inode1 != inode2) goto ok;
#endif
#endif

/* The two pids are running the same executable, raise an error and return */

	if (error_message != NULL) printf ("%s", error_message);
	return;

/* All is OK, save our pid, run, then delete our pid */

ok:
	IniWriteInt (LOCALINI_FILE, "Pid", my_pid);
	LaunchWorkerThreads (ALL_WORKERS, wait_flag);
	IniWriteInt (LOCALINI_FILE, "Pid", 0);
}

/* Load the PrimeNet DLL, make sure an internet connection is active */

int LoadPrimeNet (void)
{
	/* Init stuff */
	/* Set PRIMENET procedure pointer */
	/* return false if not connected to internet */

	int lines = 0;
#ifndef AOUT
	FILE* fd;
	char buffer[4096];
#ifdef __EMX__
	char command[128];
	char szProxyHost[120], *con_host;
	char *colon;

	IniSectionGetString (INI_FILE, "PrimeNet", "ProxyHost",
			     szProxyHost, 120, NULL);
	if (*szProxyHost) {
		if ((colon = strchr(szProxyHost, ':'))) {
			*colon = 0;
		}
		con_host = szProxyHost;
	} else {
		con_host = szSITE;
	}

	sprintf(command,"host %s",con_host);
#ifdef __DEBUG
	fprintf(stderr,"Command = %s\n",command);
#endif
	fd = popen(command,"r");
	if (fd != NULL) {
	  fgets(buffer, 199, fd);
#ifdef __DEBUG
	  fprintf(stderr,"Response = %s\n",buffer);
#endif
	  if (strncmp(buffer,"host:",5) != 0) {
	    fclose(fd);
	    return TRUE;
	  }
	  fclose(fd);
	}
#else
#ifdef __linux__
	/* Open file that will hopefully tell us if we are connected to */
	/* the Internet.  There are four possible settings for RouteRequired */
	/* 0:	Always return TRUE */
	/* 1:   Use old version 19 code */
	/* 2:   Use new code supplied by Matthew Ashton. */
	/* 99:	Default.  Use case 2 above but if cannot open /proc/net/route*/
	/*	then assume you are connected (we probably do not have read */
	/*	permission or this is a funny Linux setup). */
	{
	  int RtReq = IniSectionGetInt (INI_FILE, "PrimeNet", "RouteRequired", 99);
	  if (RtReq == 0) return (TRUE);
	  fd = fopen("/proc/net/route","r");
	  if (fd == NULL) return (RtReq == 99);
	/* We have a readable /proc/net/route file.  Use the new check */
	/* for an Internet connection written by Matthew Ashton. However, */
	/* we still support the old style check (just in case) by setting */
	/* RouteRequired to 1. */
	  if (RtReq >= 2) {
	    while (fgets(buffer, sizeof(buffer), fd)) {
	      int dest;
	      if(sscanf(buffer, "%*s %x", &dest) == 1 && dest == 0) {
		fclose (fd);
		return (TRUE);
	      }
	    }
	  }
	/* The old code for testing an Internet connection is below */
	  else {
	    fgets(buffer, 199, fd);
	    fgets(buffer, 199, fd);
	    while (!feof(fd)) {
	      if (strncmp(buffer, "lo", 2)) {
	        fclose(fd);
	        return TRUE;
	      }
	      fgets(buffer, 199, fd);
	    }
	  }
	  fclose(fd);
	}
#endif
#if defined (__FreeBSD__) || defined (__APPLE__) || defined (__WATCOMC__)
	/* The /proc/net/route test is not really meaningful under FreeBSD */
	/* There doesn't seem to be any meaningful test to see whether the */
	/* computer is connected to the Internet at the time using a non- */
	/* invasive test (which wouldn't, say, activate diald or ppp or */
	/* something else */
	return TRUE;
#endif                /* __FreeBSD__ */
#endif
#endif
	OutputStr (COMM_THREAD_NUM, "You are not connected to the Internet.\n");
	return FALSE;
}

/* Unload the PrimeNet DLL */

void UnloadPrimeNet (void)
{
}

/* Check if a program is currently running - not implemented for OS/2 */

void checkPauseListCallback (void)
{
#ifndef __OS2__
	FILE	*fd;
	char	buf[80];

	fd = popen ("ps -A -o ucomm", "r");
	if (fd != NULL) {
		while (fgets (buf, sizeof (buf), fd) != NULL) {
			int	len = strlen (buf);
			while (len && isspace (buf[len-1])) buf[--len] = 0;
			isInPauseList (buf);
		}
		fclose (fd);
	}
#endif
}
