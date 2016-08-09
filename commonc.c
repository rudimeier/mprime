/*----------------------------------------------------------------------
| Copyright 1995-2015 Mersenne Research, Inc.  All rights reserved
|
| This file contains routines and global variables that are common for
| all operating systems the program has been ported to.  It is included
| in one of the source code files of each port.  See common.h for the
| common #defines and common routine definitions.
|
| Commona contains information used only during setup
| Commonb contains information used only during execution
| Commonc contains information used during setup and execution
+---------------------------------------------------------------------*/

static const char JUNK[]="Copyright 1996-2015 Mersenne Research, Inc. All rights reserved";

char	INI_FILE[80] = {0};
char	LOCALINI_FILE[80] = {0};
char	WORKTODO_FILE[80] = {0};
char	RESFILE[80] = {0};
char	SPOOL_FILE[80] = {0};
char	LOGFILE[80] = {0};

char	USERID[21] = {0};
char	COMPID[21] = {0};
char	COMPUTER_GUID[33] = {0};
char	HARDWARE_GUID[33] = {0};
char	WINDOWS_GUID[33] = {0};
int	FIXED_GUIDS = 0;
int	USE_PRIMENET = 0;
int	DIAL_UP = 0;
unsigned int DAYS_OF_WORK = 5;
int	STRESS_TESTER = 0;
int volatile ERRCHK = 0;
int volatile SUM_INPUTS_ERRCHK = 0;	/* 1 to turn on sum(inputs) != sum(outputs) error checking */
unsigned int PRIORITY = 1;
unsigned int NUM_WORKER_THREADS = 1; /* Number of work threads to launch */
unsigned int WORK_PREFERENCE[MAX_NUM_WORKER_THREADS] = {0};
unsigned int CPU_AFFINITY[MAX_NUM_WORKER_THREADS] = {100};
unsigned int THREADS_PER_TEST[MAX_NUM_WORKER_THREADS] = {1};
				/* Number of threads gwnum can use in */
				/* computations. */
int	MANUAL_COMM = 0;
unsigned int volatile CPU_HOURS = 0;
int	CLASSIC_OUTPUT = 0;
int	OUTPUT_ROUNDOFF = 0;
unsigned long volatile ITER_OUTPUT = 0;
unsigned long volatile ITER_OUTPUT_RES = 999999999;
unsigned long volatile DISK_WRITE_TIME = 30;
unsigned int MODEM_RETRY_TIME = 2;
unsigned int NETWORK_RETRY_TIME = 70;
float	DAYS_BETWEEN_CHECKINS = 1.0;
int	NUM_BACKUP_FILES = 3;
int	SILENT_VICTORY = 0;
int	SILENT_VICTORY_PRP = 1;
int	RUN_ON_BATTERY = 1;
int	BATTERY_PERCENT = 0;
int	DEFEAT_POWER_SAVE = 1;
int	TRAY_ICON = TRUE;
int	HIDE_ICON = FALSE;
int	MERGE_WINDOWS = 0;		/* Flags indicating which MDI */
					/* windows to merge together */
double UNOFFICIAL_CPU_SPEED = 0.0;
unsigned int ROLLING_AVERAGE = 0;
unsigned int PRECISION = 2;
int	RDTSC_TIMING = 1;
int	TIMESTAMPING = 1;
int	CUMULATIVE_TIMING = 0;
int	SEQUENTIAL_WORK = 1;
int	WELL_BEHAVED_WORK = 0;
unsigned long INTERIM_FILES = 0;
unsigned long INTERIM_RESIDUES = 0;
unsigned long HYPERTHREADING_BACKOFF = 0;
int	THROTTLE_PCT = 0;

int	STARTUP_IN_PROGRESS = 0;/* True if displaying startup dialog boxes */

unsigned long NUM_CPUS = 1;	/* Number of CPUs/Cores in the computer */

gwmutex	INI_MUTEX;		/* Lock for accessing INI files */
gwmutex	INI_ADD_MUTEX;		/* Lock for accessing INI add-in files */
gwmutex	OUTPUT_MUTEX;		/* Lock for screen and results file access */
gwmutex	LOG_MUTEX;		/* Lock for prime.log access */
gwmutex	WORKTODO_MUTEX;		/* Lock for accessing worktodo structure */

int	LAUNCH_TYPE = 0;	/* Type of worker threads launched */
unsigned int WORKER_THREADS_ACTIVE = 0; /* Number of worker threads running */
int	WORKER_THREADS_STOPPING = 0; /* TRUE iff worker threads are stopping */

struct work_unit_array WORK_UNITS[MAX_NUM_WORKER_THREADS] = {0};
				/* An array of work units for each */
				/* worker thread */
unsigned int WORKTODO_COUNT = 0;/* Count of valid work lines */
unsigned int WORKTODO_IN_USE_COUNT = 0;/* Count of work units in use */
int	WORKTODO_CHANGED = 0;	/* Flag indicating worktodo file needs */
				/* writing */

#include "md5.c"

/* Generate the application string.  This is sent to the server in a */
/* UC (Update Computer info) call.  It is also displayed in the */
/* Help/About dialog box. */

void generate_application_string (
	char	*app_string)
{
#ifdef SECURITY_MODULES_PRESENT
	sprintf (app_string, "%s,Prime95,v%s,build %s",
#else
	sprintf (app_string, "%s,Untrusted Prime95,v%s,build %s",
#endif
		 PORT == 1 ? "Windows" :
		 PORT == 2 ? "Linux" :
		 PORT == 4 ? "Windows64" :
		 PORT == 5 ? "WindowsService" :
		 PORT == 6 ? "FreeBSD" :
		 PORT == 7 ? "OS/2" :
		 PORT == 8 ? "Linux64" :
		 PORT == 9 ? "Mac OS X" :
		 PORT == 10 ? "Mac OS X 64-bit" :
		 PORT == 11 ? "Haiku" :
		 PORT == 12 ? "FreeBSD 64-bit" : "Unknown",
		 VERSION, BUILD_NUM);
}

/* Calculate the 32-byte hex string for the hardware GUID.  We use the */
/* output of the CPUID function to generate this.  We don't include the */
/* cache size information because when a new processor comes out the CPUID */
/* does not recognize the cache size.  When a new version of prime95 is */
/* released that does recognize the cache size a different hardware GUID */
/* would be generated. */

void calc_hardware_guid (void)
{
	char	buf[500];

	sprintf (buf, "%s%d", CPU_BRAND, CPU_SIGNATURE);
	md5 (HARDWARE_GUID, buf);

/* Sometimes a user might want to run the program on several machines. */
/* Typically this is done by carrying the program and files around on a */
/* portable media such as floppy or USB memory stick.  In this case, */
/* we need to defeat the code that automatically detects hardware changes. */
/* The FixedHardwareUID INI option tells us to get the Windows and */
/* hardware hashes from the INI file rather than calculating them. */

	if (FIXED_GUIDS) {
		IniGetString (INI_FILE, "HardwareGUID", HARDWARE_GUID,
			      sizeof (HARDWARE_GUID), HARDWARE_GUID);
		IniWriteString (INI_FILE, "HardwareGUID", HARDWARE_GUID);
	}
}

/* Calculate the 32-byte hex string for the Windows-only GUID.  For */
/* non-Windows systems, we set WINDOWS_GUID to "". */
/* NOTE: In version 25.7 and earlier we used our first attempt at */
/* generating a unique ID.  In version 25.8 and later we use a more */
/* robust method of getting the serial number and SID.  We call the */
/* older code for all computer GUIDs that were generated by 25.7 */

void calc_windows_guid (void)
{
#ifdef _WINDOWS_
	int	algorithm_version;
	char	buf[500];

	algorithm_version = IniGetInt (INI_FILE, "WGUID_version", 1);
	if (algorithm_version == 1) {
		getWindowsSerialNumber (buf);
		getWindowsSID (buf + strlen (buf));
	} else {
		getWindowsSerialNumber_2 (buf);
		getWindowsSID_2 (buf + strlen (buf));
	}
	md5 (WINDOWS_GUID, buf);
#else
	WINDOWS_GUID[0] = 0;
#endif

/* Sometimes a user might want to run the program on several machines. */
/* Typically this is done by carrying the program and files around on a */
/* portable media such as floppy or USB memory stick.  In this case, */
/* we need to defeat the code that automatically detects hardware changes. */
/* The FixedHardwareUID INI option tells us to get the Windows and */
/* hardware hashes from the INI file rather than calculating them. */
/* NOTE: In a dual boot situation where Linux has already written out */
/* an empty WindowsGUID, then this code will write out the non-empty */
/* WindowsGUID when run on a Windows machine.  The server must not */
/* generate a CPU_IDENTITY_MISMATCH error in this case. */

	if (FIXED_GUIDS) {
		IniGetString (INI_FILE, "WindowsGUID", WINDOWS_GUID,
			      sizeof (WINDOWS_GUID), WINDOWS_GUID);
		IniWriteString (INI_FILE, "WindowsGUID", WINDOWS_GUID);
	}
}

/* Clear cached program options.  This is done when the server requests */
/* all program options to be resent.  This is a fail-safe that lets the */
/* client and server resync if the server detects an inconsistency (like */
/* getting an assignment for worker #2 with num_workers = 1 */

void clearCachedProgramOptions (void)
{
	int	tnum;
	char	section_name[32];

	IniWriteString (LOCALINI_FILE, "SrvrPO1", NULL);
	for (tnum = 0; tnum < MAX_NUM_WORKER_THREADS; tnum++) {
		sprintf (section_name, "Worker #%d", tnum+1);
		IniSectionWriteString (LOCALINI_FILE, section_name, "SrvrPO1", NULL);
	}
	IniWriteString (LOCALINI_FILE, "SrvrPO2", NULL);
	IniWriteString (LOCALINI_FILE, "SrvrPO3", NULL);
	IniWriteString (LOCALINI_FILE, "SrvrPO4", NULL);
	IniWriteString (LOCALINI_FILE, "SrvrPO5", NULL);
	IniWriteString (LOCALINI_FILE, "SrvrPO6", NULL);
	IniWriteString (LOCALINI_FILE, "SrvrPO7", NULL);
	IniWriteString (LOCALINI_FILE, "SrvrPO8", NULL);
	IniWriteString (LOCALINI_FILE, "SrvrPO9", NULL);
}

/* Generate a globally unique ID for this computer.  All Primenet */
/* communications are based on this value. */

void generate_computer_guid (void)
{
	char	buf[500];
	time_t	current_time;

	time (&current_time);
	sprintf (buf, "%s%d%f%d",
		 CPU_BRAND, CPU_SIGNATURE, CPU_SPEED, (int) current_time);
	md5 (COMPUTER_GUID, buf);
	IniWriteString (LOCALINI_FILE, "ComputerGUID", COMPUTER_GUID);

/* Clear out local copies of what we think the server knows about this computer */
/* The server now knows nothing about this computer because of the newly generated computer ID */

	IniWriteString (LOCALINI_FILE, "SrvrUID", NULL);
	IniWriteString (LOCALINI_FILE, "SrvrComputerName", NULL);
	clearCachedProgramOptions ();

/* Since we're generating a new computer GUID, we can use the latest, */
/* most robust version of caluculating the Windows GUID. */

	IniWriteInt (INI_FILE, "WGUID_version", 2);
	calc_windows_guid ();
}

/* Determine the CPU speed either empirically or by user overrides. */
/* getCpuType must be called prior to calling this routine. */

void getCpuSpeed (void)
{
	int	temp, old_cpu_speed, report_new_cpu_speed;

/* Guess the CPU speed using the RDTSC instruction */

	guessCpuSpeed ();

/* Now let the user override the cpu speed from the local.ini file */

	temp = IniGetInt (LOCALINI_FILE, "CpuSpeed", 99);
	if (temp != 99) CPU_SPEED = temp;

/* Make sure the cpu speed is reasonable */

	if (CPU_SPEED > 50000) CPU_SPEED = 50000;
	if (CPU_SPEED < 25) CPU_SPEED = 25;

/* Set the unofficial CPU speed.  The unofficial CPU speed is the */
/* last CPU speed measurement.  The official CPU speed is the one */
/* reported to the server. */

	UNOFFICIAL_CPU_SPEED = CPU_SPEED;

/* If CPU speed is much less than the official CPU speed, then set a new */
/* official CPU speed only after several slower measurements. */
/* The reason for this is that erroneously (due to one aberrant CPU speed */
/* calculation) reducing the speed we report to the server may result */
/* in erroneously unreserving exponents. */

	report_new_cpu_speed = FALSE;
	old_cpu_speed = IniGetInt (LOCALINI_FILE, "OldCpuSpeed", 0);
	if (CPU_SPEED < (double) old_cpu_speed * 0.97) {
		if (IniGetInt (LOCALINI_FILE, "NewCpuSpeedCount", 0) <= 5) {
			if (CPU_SPEED > (double) IniGetInt (LOCALINI_FILE, "NewCpuSpeed", 0))
				IniWriteInt (LOCALINI_FILE, "NewCpuSpeed", (int) (CPU_SPEED + 0.5));
			IniWriteInt (LOCALINI_FILE, "NewCpuSpeedCount", IniGetInt (LOCALINI_FILE, "NewCpuSpeedCount", 0) + 1);
			CPU_SPEED = old_cpu_speed;
		} else {
			if (CPU_SPEED < (double) IniGetInt (LOCALINI_FILE, "NewCpuSpeed", 0))
				CPU_SPEED = (double) IniGetInt (LOCALINI_FILE, "NewCpuSpeed", 0);
			report_new_cpu_speed = TRUE;
		}
	}

/* If CPU speed is close to last reported CPU speed, then use it. */
/* tell the server, recalculate new completion dates, and reset the */
/* rolling average.  Don't do this on the first run (before the Welcome */
/* dialog has been displayed). */

	else if (CPU_SPEED < (double) old_cpu_speed * 1.03) {
		IniWriteInt (LOCALINI_FILE, "NewCpuSpeedCount", 0);
		IniWriteInt (LOCALINI_FILE, "NewCpuSpeed", 0);
	}

/* If CPU speed is much larger than the speed reported to the server, then */
/* use this new speed and tell the server. */

	else {
		report_new_cpu_speed = TRUE;
	}

/* Report a new CPU speed.  Remember the new CPU speed, tell the server, */
/* recalculate new completion dates, and reset the rolling average in */
/* such a way as to reduce the chance of spurious unreserves.  Don't */
/* do this on the first run (before the Welcome dialog has been displayed). */

	if (report_new_cpu_speed) {
		IniWriteInt (LOCALINI_FILE, "OldCpuSpeed", (int) (CPU_SPEED + 0.5));
		IniWriteInt (LOCALINI_FILE, "NewCpuSpeedCount", 0);
		IniWriteInt (LOCALINI_FILE, "NewCpuSpeed", 0);
		if (old_cpu_speed) {
			if (WORKTODO_COUNT) {
				ROLLING_AVERAGE = (int) (ROLLING_AVERAGE * old_cpu_speed / CPU_SPEED);
				if (ROLLING_AVERAGE < 1000) ROLLING_AVERAGE = 1000;
			}
			else
				ROLLING_AVERAGE = 1000;
			IniWriteInt (LOCALINI_FILE, "RollingAverage", ROLLING_AVERAGE);
			IniWriteInt (LOCALINI_FILE, "RollingStartTime", 0);
			spoolMessage (PRIMENET_UPDATE_COMPUTER_INFO, NULL);
			UpdateEndDates ();
		}
	}
}

/* Set the CPU flags based on the CPUID instruction.  Also, the advanced */
/* user can override our guesses. */

void getCpuInfo (void)
{
	int	temp;

/* Get the CPU info using CPUID instruction */	

	guessCpuType ();

/* Allow overriding the number of physical processors.  For historical */
/* reasons, this code uses a variable called NUM_CPUS rather */
/* than the CPU_CORES value set by guessCpuType. */

	NUM_CPUS = IniGetInt (LOCALINI_FILE, "NumCPUs", CPU_CORES);
	temp = IniGetInt (LOCALINI_FILE, "NumPhysicalCores", 9999);
	if (temp != 9999) {
		CPU_HYPERTHREADS = NUM_CPUS / temp;
		if (CPU_HYPERTHREADS < 1) CPU_HYPERTHREADS = 1;
		NUM_CPUS = temp;
	}

/* Calculate hardware GUID (global unique identifier) using the CPUID info. */
/* Well, it isn't unique but it is about as good as we can do and still have */
/* portable code.  Do this calculation before user overrides values */
/* derived from CPUID results. */

	calc_hardware_guid ();

/* Let the user override the cpu flags from the local.ini file */

	temp = IniGetInt (LOCALINI_FILE, "CpuSupportsRDTSC", 99);
	if (temp == 0) CPU_FLAGS &= ~CPU_RDTSC;
	if (temp == 1) CPU_FLAGS |= CPU_RDTSC;
	temp = IniGetInt (LOCALINI_FILE, "CpuSupportsCMOV", 99);
	if (temp == 0) CPU_FLAGS &= ~CPU_CMOV;
	if (temp == 1) CPU_FLAGS |= CPU_CMOV;
	temp = IniGetInt (LOCALINI_FILE, "CpuSupportsPrefetch", 99);
	if (temp == 0) CPU_FLAGS &= ~CPU_PREFETCH;
	if (temp == 1) CPU_FLAGS |= CPU_PREFETCH;
	temp = IniGetInt (LOCALINI_FILE, "CpuSupportsSSE", 99);
	if (temp == 0) CPU_FLAGS &= ~CPU_SSE;
	if (temp == 1) CPU_FLAGS |= CPU_SSE;
#ifndef X86_64
	temp = IniGetInt (LOCALINI_FILE, "CpuSupportsSSE2", 99);
	if (temp == 0) CPU_FLAGS &= ~CPU_SSE2;
	if (temp == 1) CPU_FLAGS |= CPU_SSE2;
#else
	CPU_FLAGS |= CPU_SSE2;
#endif
	temp = IniGetInt (LOCALINI_FILE, "CpuSupportsSSE4", 99);
	if (temp == 0) CPU_FLAGS &= ~CPU_SSE41;
	if (temp == 1) CPU_FLAGS |= CPU_SSE41;
	temp = IniGetInt (LOCALINI_FILE, "CpuSupports3DNow", 99);
	if (temp == 0) CPU_FLAGS &= ~(CPU_3DNOW + CPU_3DNOW_PREFETCH);
	if (temp == 1) CPU_FLAGS |= (CPU_3DNOW + CPU_3DNOW_PREFETCH);
	temp = IniGetInt (LOCALINI_FILE, "CpuSupportsAVX", 99);
	if (temp == 0) CPU_FLAGS &= ~CPU_AVX;
	if (temp == 1) CPU_FLAGS |= CPU_AVX;
	temp = IniGetInt (LOCALINI_FILE, "CpuSupportsFMA3", 99);
	if (temp == 0) CPU_FLAGS &= ~CPU_FMA3;
	if (temp == 1) CPU_FLAGS |= CPU_FMA3;
	temp = IniGetInt (LOCALINI_FILE, "CpuSupportsFMA4", 99);
	if (temp == 0) CPU_FLAGS &= ~CPU_FMA4;
	if (temp == 1) CPU_FLAGS |= CPU_FMA4;
	temp = IniGetInt (LOCALINI_FILE, "CpuSupportsAVX2", 99);
	if (temp == 0) CPU_FLAGS &= ~CPU_AVX2;
	if (temp == 1) CPU_FLAGS |= CPU_AVX2;

/* Let the user override the L2 cache size in local.ini file */

	CPU_L2_CACHE_SIZE =
		IniGetInt (LOCALINI_FILE, "CpuL2CacheSize", CPU_L2_CACHE_SIZE);
	CPU_L2_CACHE_LINE_SIZE =
		IniGetInt (LOCALINI_FILE, "CpuL2CacheLineSize", CPU_L2_CACHE_LINE_SIZE);
	CPU_L2_SET_ASSOCIATIVE =
		IniGetInt (LOCALINI_FILE, "CpuL2SetAssociative", CPU_L2_SET_ASSOCIATIVE);

/* Let the user override the L3 cache size in local.ini file */

	CPU_L3_CACHE_SIZE =
		IniGetInt (LOCALINI_FILE, "CpuL3CacheSize", CPU_L3_CACHE_SIZE);
	CPU_L3_CACHE_LINE_SIZE =
		IniGetInt (LOCALINI_FILE, "CpuL3CacheLineSize", CPU_L3_CACHE_LINE_SIZE);
	CPU_L3_SET_ASSOCIATIVE =
		IniGetInt (LOCALINI_FILE, "CpuL3SetAssociative", CPU_L3_SET_ASSOCIATIVE);

/* Let the user override the CPUID brand string.  It should never be necessary. */
/* However, one Athlon owner's brand string became corrupted with illegal characters. */
	
	IniGetString (LOCALINI_FILE, "CpuBrand", CPU_BRAND, sizeof(CPU_BRAND), CPU_BRAND);

/* Let user override the number of hyperthreads */

	CPU_HYPERTHREADS =
		IniGetInt (LOCALINI_FILE, "CpuNumHyperthreads", CPU_HYPERTHREADS);
	if (CPU_HYPERTHREADS == 0) CPU_HYPERTHREADS = 1;

/* Let user override the CPU architecture */

	CPU_ARCHITECTURE = IniGetInt (LOCALINI_FILE, "CpuArchitecture", CPU_ARCHITECTURE);

/* Now get the CPU speed */

	getCpuSpeed ();
}

/* Format a long or very long textual cpu description */

void getCpuDescription (
	char	*buf,			/* A 512 byte buffer */
	int	long_desc)		/* True for a very long description */
{

/* Recalculate the CPU speed in case speed step has changed the original */
/* settings. */

	getCpuSpeed ();

/* Now format a pretty CPU description */

	sprintf (buf, "%s\nCPU speed: %.2f MHz", CPU_BRAND, UNOFFICIAL_CPU_SPEED);
	if (NUM_CPUS > 1 && CPU_HYPERTHREADS > 1)
		sprintf (buf + strlen (buf), ", %lu hyperthreaded cores", NUM_CPUS);
	else if (NUM_CPUS > 1)
		sprintf (buf + strlen (buf), ", %lu cores", NUM_CPUS);
	else if (CPU_HYPERTHREADS > 1)
		sprintf (buf + strlen (buf), ", with hyperthreading");
	strcat (buf, "\n");
	if (CPU_FLAGS) {
		strcat (buf, "CPU features: ");
//		if (CPU_FLAGS & CPU_RDTSC) strcat (buf, "RDTSC, ");
//		if (CPU_FLAGS & CPU_CMOV) strcat (buf, "CMOV, ");
//		if (CPU_FLAGS & CPU_MMX) strcat (buf, "MMX, ");
		if (CPU_FLAGS & CPU_3DNOW) strcat (buf, "3DNow!, ");
		else if (CPU_FLAGS & CPU_3DNOW_PREFETCH) strcat (buf, "3DNow! Prefetch, ");
		else if (CPU_FLAGS & CPU_PREFETCH) strcat (buf, "Prefetch, ");
		if (CPU_FLAGS & CPU_SSE) strcat (buf, "SSE, ");
		if (CPU_FLAGS & CPU_SSE2) strcat (buf, "SSE2, ");
		if (CPU_FLAGS & CPU_SSE41) strcat (buf, "SSE4, ");
		if (CPU_FLAGS & CPU_AVX) strcat (buf, "AVX, ");
		if (CPU_FLAGS & CPU_AVX2) strcat (buf, "AVX2, ");
		if (CPU_FLAGS & (CPU_FMA3 | CPU_FMA4)) strcat (buf, "FMA, ");
		strcpy (buf + strlen (buf) - 2, "\n");
	}
	strcat (buf, "L1 cache size: ");
	if (CPU_L1_CACHE_SIZE < 0) strcat (buf, "unknown\n");
	else sprintf (buf + strlen (buf), "%d KB\n", CPU_L1_CACHE_SIZE);
	strcat (buf, "L2 cache size: ");
	if (CPU_L2_CACHE_SIZE < 0) strcat (buf, "unknown\n");
	else {
		if (CPU_L2_CACHE_SIZE & 0x3FF)
			sprintf (buf + strlen (buf), "%d KB\n", CPU_L2_CACHE_SIZE);
		else
			sprintf (buf + strlen (buf), "%d MB\n", CPU_L2_CACHE_SIZE >> 10);
	}
	if (CPU_L3_CACHE_SIZE > 0) {
		if (CPU_L3_CACHE_SIZE & 0x3FF)
			sprintf (buf + strlen (buf) - 1, ", L3 cache size: %d KB\n", CPU_L3_CACHE_SIZE);
		else
			sprintf (buf + strlen (buf) - 1, ", L3 cache size: %d MB\n", CPU_L3_CACHE_SIZE >> 10);
	}
	if (! long_desc) return;
	strcat (buf, "L1 cache line size: ");
	if (CPU_L1_CACHE_LINE_SIZE < 0) strcat (buf, "unknown\n");
	else sprintf (buf+strlen(buf), "%d bytes\n", CPU_L1_CACHE_LINE_SIZE);
	strcat (buf, "L2 cache line size: ");
	if (CPU_L2_CACHE_LINE_SIZE < 0) strcat (buf, "unknown\n");
	else sprintf (buf+strlen(buf), "%d bytes\n", CPU_L2_CACHE_LINE_SIZE);
	if (CPU_L1_DATA_TLBS > 0)
		sprintf (buf + strlen (buf), "L1 TLBS: %d\n", CPU_L1_DATA_TLBS);
	if (CPU_L2_DATA_TLBS > 0)
		sprintf (buf + strlen (buf), "%sTLBS: %d\n",
			 CPU_L1_DATA_TLBS > 0 ? "L2 " : "",
			 CPU_L2_DATA_TLBS);
}

/* Determine if a number is prime */

int isPrime (
	unsigned long p)
{
	unsigned long i;
	for (i = 2; i * i <= p; i = (i + 1) | 1)
		if (p % i == 0) return (FALSE);
	return (TRUE);
}

/* Upper case a string */

void strupper (
	char	*p)
{
	for ( ; *p; p++) if (*p >= 'a' && *p <= 'z') *p = *p - 'a' + 'A';
}

/* Convert a string (e.g "11:30 AM") to minutes since midnight */

unsigned int strToMinutes (
	const char *buf)
{
	unsigned int hours, minutes, pm;

	pm = (strchr (buf, 'P') != NULL || strchr (buf, 'p') != NULL);
	hours = atoi (buf);
	while (isdigit (*buf)) buf++;
	while (*buf && ! isdigit (*buf)) buf++;
	minutes = atoi (buf);
	if (hours == 12) hours -= 12;
	if (pm) hours += 12;
	minutes = hours * 60 + minutes;
	minutes %= 1440;  // In case user entered 24:00 or 11:60 PM
	return (minutes);
}

/* Convert minutes since midnight to a string (e.g "11:30 AM") */

void minutesToStr (
	unsigned int minutes,
	char	*buf)
{
	unsigned int fmt_type, hours, pm;

	hours = minutes / 60;
	fmt_type = IniGetInt (INI_FILE, "AMPM", 0);
	if (fmt_type == 0) fmt_type = getDefaultTimeFormat ();
	if (fmt_type != 1) {
		sprintf (buf, "%d:%02d", hours, minutes % 60);
	} else {
		if (hours >= 12) hours -= 12, pm = 1;
		else pm = 0;
		if (hours == 0) hours = 12;
		sprintf (buf, "%d:%02d %s", hours, minutes % 60, pm ? "PM" : "AM");
	}
}

/* Convert the old DAY_MEMORY, NIGHT_MEMORY settings into the simpler */
/* and more powerful MEMORY setting. */

void write_memory_settings (
	unsigned int day_memory,
	unsigned int night_memory,
	unsigned int day_start_time,
	unsigned int day_end_time)
{
	char	buf[100];

	sprintf (buf, "%d during %d:%02d-%d:%02d else %d",
		 day_memory, day_start_time / 60, day_start_time % 60,
		 day_end_time / 60, day_end_time % 60, night_memory);
	IniWriteString (LOCALINI_FILE, "Memory", buf);
	IniWriteString (LOCALINI_FILE, "DayMemory", NULL);
	IniWriteString (LOCALINI_FILE, "NightMemory", NULL);
	IniWriteString (LOCALINI_FILE, "DayStartTime", NULL);
	IniWriteString (LOCALINI_FILE, "DayEndTime", NULL);
}

/* Convert the new MEMORY setting into the old DAY_MEMORY, NIGHT_MEMORY settings. */
/* Returns FALSE if the MEMORY setting cannot be converted. */

int read_memory_settings (
	unsigned int *day_memory,
	unsigned int *night_memory,
	unsigned int *day_start_time,
	unsigned int *day_end_time)
{
	const char *p, *during_clause, *else_clause;

/* Set up some default values */

	*day_memory = 8;
	*night_memory = 8;
	*day_start_time = 450;
	*day_end_time = 1410;

/* Get the memory settings.  If not found, return some defaults */

	p = IniSectionGetStringRaw (LOCALINI_FILE, NULL, "Memory");
	if (p == NULL) return (TRUE);

/* Find the time specifiers.  If none found, return the same value for */
/* day and night memory. */

	*day_memory = atol (p);
	during_clause = strstr (p, " during ");
	if (during_clause == NULL) {
		*night_memory = *day_memory;
		return (TRUE);
	}

/* Parse the time value to get the day start time and day end time */
/* It must be the exact same syntax as written by write_memory_settings */

	p = during_clause + 8;
	*day_start_time = atol (p) * 60;
	while (isdigit (*p)) p++;
	if (*p++ != ':') return (FALSE);
	*day_start_time += atol (p);
	while (isdigit (*p)) p++;
	if (*p++ != '-') return (FALSE);
	*day_end_time = atol (p) * 60;
	while (isdigit (*p)) p++;
	if (*p++ != ':') return (FALSE);
	*day_end_time += atol (p);
	while (isdigit (*p)) p++;

/* If user hand-edited the prime.txt file and entered an illegal time */
/* value, then arbitrarily set the value within the proper range */

	*day_start_time %= 1440;
	*day_end_time %= 1440;

/* Handle the else clause */

	else_clause = strstr (p, " else ");
	if (p != else_clause) return (FALSE);
	*night_memory = atol (p + 6);

	return (strstr (p, " during ") == NULL);
}

/* Determine the names of the INI files, then read them.  This is also the */
/* perfect time to initialize mutexes and do other initializations. */

void nameAndReadIniFiles (
	int	named_ini_files)
{
	char	buf[513];

/* Initialize mutexes */

	gwmutex_init (&INI_MUTEX);
	gwmutex_init (&INI_ADD_MUTEX);
	gwmutex_init (&OUTPUT_MUTEX);
	gwmutex_init (&LOG_MUTEX);
	gwmutex_init (&WORKTODO_MUTEX);

/* Figure out the names of the INI files */

	if (named_ini_files < 0) {
		strcpy (INI_FILE, "prime.ini");
		strcpy (LOCALINI_FILE, "local.ini");
		strcpy (SPOOL_FILE, "prime.spl");
		strcpy (WORKTODO_FILE, "worktodo.ini");
		strcpy (RESFILE, "results.txt");
		strcpy (LOGFILE, "prime.log");
	} else {
		sprintf (INI_FILE, "prim%04d.ini", named_ini_files);
		sprintf (LOCALINI_FILE, "loca%04d.ini", named_ini_files);
		sprintf (SPOOL_FILE, "prim%04d.spl", named_ini_files);
		sprintf (WORKTODO_FILE, "work%04d.ini", named_ini_files);
		sprintf (RESFILE, "resu%04d.txt", named_ini_files);
		sprintf (LOGFILE, "prim%04d.log", named_ini_files);
	}

/* Let the user rename these files and pick a different working directory */

	IniGetString (INI_FILE, "WorkingDir", buf, sizeof(buf), NULL);
	IniGetString (INI_FILE, "local.ini", LOCALINI_FILE, 80, LOCALINI_FILE);
	IniGetString (INI_FILE, "prime.spl", SPOOL_FILE, 80, SPOOL_FILE);
	IniGetString (INI_FILE, "worktodo.ini", WORKTODO_FILE, 80, WORKTODO_FILE);
	IniGetString (INI_FILE, "results.txt", RESFILE, 80, RESFILE);
	IniGetString (INI_FILE, "prime.log", LOGFILE, 80, LOGFILE);
	IniGetString (INI_FILE, "prime.ini", INI_FILE, 80, INI_FILE);
	if (buf[0]) {
		(void) _chdir (buf);
		IniFileReread (INI_FILE);
	}

/* Merge an old primenet.ini file into a special section of prime.ini */

	if (fileExists ("primenet.ini"))
		iniAddFileMerge (INI_FILE, "primenet.ini", "PrimeNet");

/* Perform other one-time initializations */

	LoadPrimenet ();
	init_spool_file_and_comm_code ();
	init_timed_event_handler ();

/* Create and name the main window */

	MERGE_WINDOWS = (int) IniGetInt (INI_FILE, "MergeWindows", MERGE_MAINCOMM_WINDOWS);
	create_window (MAIN_THREAD_NUM);
	base_title (MAIN_THREAD_NUM, "Main thread");

/* Output a startup message */

	sprintf (buf, "Mersenne number primality test program version %s\n", VERSION);
	OutputStr (MAIN_THREAD_NUM, buf);

/* Now that we have the proper file name, read the ini file into */
/* global variables. */

	readIniFiles ();

/* Generate the computer's UID if none exists */

	if (!COMPUTER_GUID[0]) generate_computer_guid ();

/* Calculate the hopefully-never-changing Windows GUID */

	calc_windows_guid ();

/* A stress tester ceases to be a stress tester if he ever turns on */
/* primenet or has work in his worktodo.ini file */

	if (STRESS_TESTER == 1 && (USE_PRIMENET || WORKTODO_COUNT)) {
		STRESS_TESTER = 0;
		IniWriteInt (INI_FILE, "StressTester", 0);
	}

/* Output our calculated CPU architecture and characteristics */

	sprintf (buf, "Optimizing for CPU architecture: %s, ",
		 CPU_ARCHITECTURE == CPU_ARCHITECTURE_PRE_SSE2 ? "Pre-SSE2" :
		 CPU_ARCHITECTURE == CPU_ARCHITECTURE_PENTIUM_4 ? "Pentium 4" :
		 CPU_ARCHITECTURE == CPU_ARCHITECTURE_PENTIUM_M ? "Pentium M" :
		 CPU_ARCHITECTURE == CPU_ARCHITECTURE_CORE ? "Core Solo/Duo" :
		 CPU_ARCHITECTURE == CPU_ARCHITECTURE_CORE_2 ? "Core 2" :
		 CPU_ARCHITECTURE == CPU_ARCHITECTURE_CORE_I7 ? "Core i3/i5/i7" :
		 CPU_ARCHITECTURE == CPU_ARCHITECTURE_ATOM ? "Atom" :
		 CPU_ARCHITECTURE == CPU_ARCHITECTURE_INTEL_OTHER ? "Unknown Intel" :
		 CPU_ARCHITECTURE == CPU_ARCHITECTURE_AMD_K8 ? "AMD K8" :
		 CPU_ARCHITECTURE == CPU_ARCHITECTURE_AMD_K10 ? "AMD K10" :
		 CPU_ARCHITECTURE == CPU_ARCHITECTURE_AMD_BULLDOZER ? "AMD Bulldozer" :
		 CPU_ARCHITECTURE == CPU_ARCHITECTURE_OTHER ? "Not Intel and not AMD" : "Undefined");
	strcat (buf, "L2 cache size: ");
	if (CPU_L2_CACHE_SIZE < 0) strcat (buf, "unknown");
	else if (CPU_L2_CACHE_SIZE & 0x3FF)
		sprintf (buf + strlen (buf), "%d KB", CPU_L2_CACHE_SIZE);
	else
		sprintf (buf + strlen (buf), "%d MB", CPU_L2_CACHE_SIZE >> 10);
	if (CPU_L3_CACHE_SIZE > 0) {
		if (CPU_L3_CACHE_SIZE & 0x3FF)
			sprintf (buf + strlen (buf), ", L3 cache size: %d KB", CPU_L3_CACHE_SIZE);
		else
			sprintf (buf + strlen (buf), ", L3 cache size: %d MB", CPU_L3_CACHE_SIZE >> 10);
	}
	strcat (buf, "\n");
	OutputStr (MAIN_THREAD_NUM, buf);

/* Dynamically determine which logical hyperthreaded CPUs map to physical CPUs */

	generate_affinity_scramble ();

/* Start some initial timers */

	add_timed_event (TE_ROLLING_AVERAGE, 6*60*60);
}

/* Init the communications code.  We used to do this at the end of */
/* nameAndReadIniFiles, but running mprime with the -s or -t argument */
/* caused spurious creation of a prime.spl file. */

void initCommCode (void) {

/* Start or stop the communication timers.  This needs to be called */
/* every time the INI file is read in case there have been changes to */
/* the USE_PRIMENET or MANUAL_COMM variables. */

//bug - do this on all rereads of prime.ini?  If so, only after comm windows
//bug - is created (and I'm not sure the comm window should be recreated
//bug - on every prime.ini reread (and what of the windows bugs where
//bug - calling create_window from other than the main thread can hang
	set_comm_timers ();

/* Tell the server about any changed program options */

	spoolMessage (PRIMENET_PROGRAM_OPTIONS, NULL);

/* If we're using primenet and we don't know the userid, then send an */
/* update computer message to get the userid. */

	if (USE_PRIMENET && USERID[0] == 0)
		spoolMessage (PRIMENET_UPDATE_COMPUTER_INFO, NULL);
}

/* Read or re-read the INI files & and do other initialization */

int readIniFiles (void)
{
	int	rc, temp;
	int	day_memory;
	int	night_memory;
	char	buf[80];

/* Force the INI files to be reread, just in case they were hand edited. */
/* Incorporate any additions from .add files */

	IniFileReread (INI_FILE);
	IniFileReread (LOCALINI_FILE);
	incorporateIniAddFiles ();

/* Get the CPU type, speed, and capabilities. */

	getCpuInfo ();

/* Convert V4 work preferences to v5 work preferences */

	if (!IniGetInt (INI_FILE, "V24OptionsConverted", 0)) {
		IniWriteInt (INI_FILE, "V24OptionsConverted", 1);
		if (IniGetInt (INI_FILE, "WorkPreference", 0) == 16) // Ten million digit
			IniWriteInt (INI_FILE, "WorkPreference", PRIMENET_WP_WORLD_RECORD);
		else if (IniGetInt (INI_FILE, "WorkPreference", 0) == 2)
			IniWriteInt (INI_FILE, "WorkPreference", PRIMENET_WP_LL_FIRST);
		else if (IniGetInt (INI_FILE, "WorkPreference", 0) == 4)
			IniWriteInt (INI_FILE, "WorkPreference", PRIMENET_WP_LL_DBLCHK);
		else if (IniGetInt (INI_FILE, "WorkPreference", 0) == 1)
			IniWriteInt (INI_FILE, "WorkPreference", PRIMENET_WP_FACTOR);
		IniSectionWriteInt (INI_FILE, "PrimeNet", "Debug", 0);
	}

/* Put commonly used options in global variables */

	IniGetString (LOCALINI_FILE, "ComputerGUID", COMPUTER_GUID, sizeof (COMPUTER_GUID), NULL);
	FIXED_GUIDS = IniGetInt (INI_FILE, "FixedHardwareUID", 0);

	IniGetString (INI_FILE, "V5UserID", USERID, sizeof (USERID), NULL);
	IniGetString (LOCALINI_FILE, "ComputerID", COMPID, sizeof (COMPID), NULL);
	sanitizeString (COMPID);
	USE_PRIMENET = (int) IniGetInt (INI_FILE, "UsePrimenet", 0);
	DIAL_UP = (int) IniGetInt (INI_FILE, "DialUp", 0);
	DAYS_OF_WORK = (unsigned int) IniGetInt (INI_FILE, "DaysOfWork", 5);
	if (DAYS_OF_WORK > 180) DAYS_OF_WORK = 180;

	CPU_HOURS = (unsigned int) IniGetInt (LOCALINI_FILE, "CPUHours", 24);
	if (CPU_HOURS < 1) CPU_HOURS = 1;
	if (CPU_HOURS > 24) CPU_HOURS = 24;

	ROLLING_AVERAGE = (unsigned int) IniGetInt (LOCALINI_FILE, "RollingAverage", 1000);
	if (ROLLING_AVERAGE < 10) ROLLING_AVERAGE = 10;
	if (ROLLING_AVERAGE > 4000) ROLLING_AVERAGE = 4000;
	/* Rolling average needs to be reset when AVX capable machines upgrade from v26 to v27. */
	if (CPU_FLAGS & CPU_AVX && ! IniGetInt (LOCALINI_FILE, "RollingAverageIsFromV27", 0)) {
		ROLLING_AVERAGE = 1000;
		IniWriteInt (LOCALINI_FILE, "RollingAverage", 1000);
		IniWriteInt (LOCALINI_FILE, "RollingAverageIsFromV27", 1);
	}

	PRECISION = (unsigned int) IniGetInt (INI_FILE, "PercentPrecision", 2);
	if (PRECISION > 6) PRECISION = 6;

	CLASSIC_OUTPUT = IniGetInt (INI_FILE, "ClassicOutput", 0);
	OUTPUT_ROUNDOFF = IniGetInt (INI_FILE, "OutputRoundoff", 0);
	ITER_OUTPUT = IniGetInt (INI_FILE, "OutputIterations", 10000);
	if (ITER_OUTPUT > 999999999) ITER_OUTPUT = 999999999;
	if (ITER_OUTPUT <= 0) ITER_OUTPUT = 1;
	ITER_OUTPUT_RES = IniGetInt (INI_FILE, "ResultsFileIterations", 999999999);
	if (ITER_OUTPUT_RES > 999999999) ITER_OUTPUT_RES = 999999999;
	if (ITER_OUTPUT_RES < 1000) ITER_OUTPUT_RES = 1000;
	DISK_WRITE_TIME = IniGetInt (INI_FILE, "DiskWriteTime", 30);
	MODEM_RETRY_TIME = (unsigned int) IniGetInt (INI_FILE, "NetworkRetryTime", 2);
	if (MODEM_RETRY_TIME < 1) MODEM_RETRY_TIME = 1;
	if (MODEM_RETRY_TIME > 300) MODEM_RETRY_TIME = 300;
	NETWORK_RETRY_TIME = (unsigned int)
		IniGetInt (INI_FILE, "NetworkRetryTime2",
			   MODEM_RETRY_TIME > 70 ? MODEM_RETRY_TIME : 70);
	if (NETWORK_RETRY_TIME < 1) NETWORK_RETRY_TIME = 1;
	if (NETWORK_RETRY_TIME > 300) NETWORK_RETRY_TIME = 300;
	
	IniGetString (INI_FILE, "DaysBetweenCheckins", buf, sizeof (buf), "1");
	DAYS_BETWEEN_CHECKINS = (float) atof (buf);
	if (DAYS_BETWEEN_CHECKINS > 7.0) DAYS_BETWEEN_CHECKINS = 7.0;				/* 7 day maximum */
	if (DAYS_BETWEEN_CHECKINS * 24.0 < 1.0) DAYS_BETWEEN_CHECKINS = (float) (1.0 / 24.0);	/* 1 hour minimum */
	SILENT_VICTORY = (int) IniGetInt (INI_FILE, "SilentVictory", 0);
	SILENT_VICTORY_PRP = (int) IniGetInt (INI_FILE, "SilentVictoryPRP", 1);
	RUN_ON_BATTERY = (int) IniGetInt (LOCALINI_FILE, "RunOnBattery", 1);
	BATTERY_PERCENT = (int) IniGetInt (INI_FILE, "BatteryPercent", 0);
	DEFEAT_POWER_SAVE = (int) IniGetInt (LOCALINI_FILE, "DefeatPowerSave", 1);

	STRESS_TESTER = (int) IniGetInt (INI_FILE, "StressTester", 99);
	temp = (int) IniGetInt (INI_FILE, "ErrorCheck", 0);
	ERRCHK = (temp != 0);
	temp = (int) IniGetInt (INI_FILE, "SumInputsErrorCheck", 0);
	SUM_INPUTS_ERRCHK = (temp != 0);
	NUM_WORKER_THREADS = IniGetInt (LOCALINI_FILE, "WorkerThreads", NUM_CPUS);
	if (NUM_WORKER_THREADS < 1) NUM_WORKER_THREADS = 1;
	if (NUM_WORKER_THREADS > MAX_NUM_WORKER_THREADS) NUM_WORKER_THREADS = MAX_NUM_WORKER_THREADS;
	PRIORITY = (unsigned int) IniGetInt (INI_FILE, "Priority", 1);
	if (PRIORITY < 1) PRIORITY = 1;
	if (PRIORITY > 10) PRIORITY = 10;
	PTOGetAll (INI_FILE, "WorkPreference", WORK_PREFERENCE, 0);
	PTOGetAll (LOCALINI_FILE, "Affinity", CPU_AFFINITY, 100);
	PTOGetAll (LOCALINI_FILE, "ThreadsPerTest", THREADS_PER_TEST, 1);
	MANUAL_COMM = (int) IniGetInt (INI_FILE, "ManualComm", 0);
	HIDE_ICON = (int) IniGetInt (INI_FILE, "HideIcon", 0);
	TRAY_ICON = (int) IniGetInt (INI_FILE, "TrayIcon", 1);
	MERGE_WINDOWS = (int) IniGetInt (INI_FILE, "MergeWindows", MERGE_MAINCOMM_WINDOWS);

/* Convert old TwoBackupFiles boolean to new NumBackupFiles integer.  Old default */
/* was 2 save files, new default is 3 save files. */	

	temp = (int) IniGetInt (INI_FILE, "TwoBackupFiles", 2);
	NUM_BACKUP_FILES = (int) IniGetInt (INI_FILE, "NumBackupFiles", temp+1);

/* Convert the old DAY_MEMORY, NIGHT_MEMORY settings into the simpler, all-inclusive */
/* and more powerful MEMORY setting. */

	day_memory = IniGetInt (LOCALINI_FILE, "DayMemory", 0);
	night_memory = IniGetInt (LOCALINI_FILE, "NightMemory", 0);
	if (day_memory || night_memory) {
		int	day_start_time, day_end_time;
		day_memory = (day_memory < 0) ? physical_memory () + day_memory : day_memory;
		if (day_memory < 8) day_memory = 8;
		if ((unsigned long) day_memory > physical_memory () - 8) day_memory = physical_memory () - 8;
		night_memory = (night_memory < 0) ? physical_memory () + night_memory : night_memory;
		if (night_memory < 8) night_memory = 8;
		if ((unsigned long) night_memory > physical_memory () - 8) night_memory = physical_memory () - 8;
		day_start_time = IniGetInt (LOCALINI_FILE, "DayStartTime", 450);
		day_end_time = IniGetInt (LOCALINI_FILE, "DayEndTime", 1410);
		write_memory_settings (day_memory, night_memory, day_start_time, day_end_time);
	}
	read_mem_info ();

/* Get the option controlling which timer to use.  If the high resolution */
/* performance counter is not available on this machine, then add 10 to */
/* the RDTSC_TIMING value. */

	RDTSC_TIMING = IniGetInt (INI_FILE, "RdtscTiming", 1);
	if (RDTSC_TIMING < 10 && ! isHighResTimerAvailable ())
		RDTSC_TIMING += 10;

/* Other oddball options */

	TIMESTAMPING = IniGetInt (INI_FILE, "TimeStamp", 1);
	CUMULATIVE_TIMING = IniGetInt (INI_FILE, "CumulativeTiming", 0);
	SEQUENTIAL_WORK = IniGetInt (INI_FILE, "SequentialWorkToDo", 1);
	WELL_BEHAVED_WORK = IniGetInt (INI_FILE, "WellBehavedWork", 0);

	read_pause_info ();
	read_load_average_info ();

	INTERIM_FILES = IniGetInt (INI_FILE, "InterimFiles", 0);
	INTERIM_RESIDUES = IniGetInt (INI_FILE, "InterimResidues", INTERIM_FILES);
	HYPERTHREADING_BACKOFF =
		IniGetInt (INI_FILE, "HyperthreadingBackoff",
0);//bug		   CPU_HYPERTHREADS <= 1 ? 0 : 30);

/* Option to slow down the program by sleeping after every iteration.  You */
/* might use this on a laptop or a computer running in a hot room to keep */
/* temperatures down and thus reduce the chance of a hardware error.  */

	THROTTLE_PCT = IniGetInt (INI_FILE, "Throttle", 0);

/* Now read the work-to-do file */

	rc = readWorkToDoFile ();
	if (rc) {
		OutputStr (MAIN_THREAD_NUM, "Error reading worktodo.txt file\n");
		return (rc);
	}

/* Return success */

	return (0);
}

/****************************************************************************/
/*          Portable routines to read and write INI files.                  */
/****************************************************************************/

/*----------------------------------------------------------------------
| NOTE:  These only work if you open no more than 10 ini files.  Also you must
| not change the working directory at any time during program execution.
+---------------------------------------------------------------------*/

#define INI_LINE_NORMAL		0	/* A normal keyword=value INI line */
#define INI_LINE_COMMENT	2	/* A comment line */
#define INI_LINE_HEADER		3	/* A section header line */

struct IniLine {
	char	*keyword;
	char	*value;
	int	line_type;
};
struct IniCache {
	char	*filename;
	int	immediate_writes;
	int	dirty;
	unsigned int num_lines;
	unsigned int array_size;
	struct IniLine **lines;
};

void growIniLineArray (
	struct IniCache *p)
{
	struct IniLine **newlines;

	if (p->num_lines != p->array_size) return;

	newlines = (struct IniLine **)
		malloc ((p->num_lines + 100) * sizeof (struct IniLine *));
	if (p->num_lines) {
		memcpy (newlines, p->lines, p->num_lines * sizeof (struct IniLine *));
		free (p->lines);
	}
	p->lines = newlines;
	p->array_size = p->num_lines + 100;
}

/* We used to name files xxxx.ini.  unfortunately, Windows backup/restore */
/* thought these files should be restored to there old values when you */
/* rollback a driver and/or some other reasons.  Now we name our files */
/* xxxx.txt.  This routine converts an old ini name to a new txt name. */

void mangleIniFileName (
	const char *inputFileName,
	char	*outputFileName,
	int	*mangled)	/* Return TRUE if output name != input name */
{
	int	len;
	strcpy (outputFileName, inputFileName);
	len = (int) strlen (outputFileName);
	if (len >= 4 &&
	    outputFileName[len-4] == '.' &&
	    (outputFileName[len-3] == 'I' || outputFileName[len-3] == 'i') &&
	    (outputFileName[len-2] == 'N' || outputFileName[len-2] == 'n') &&
	    (outputFileName[len-1] == 'I' || outputFileName[len-1] == 'i')) {
		strcpy (outputFileName+len-3, "txt");
		*mangled = TRUE;
	} else
		*mangled = FALSE;
}

struct IniCache *openIniFile (
	const char *filename,
	int	forced_read)
{
static	struct IniCache *cache[10] = {0};
	struct IniCache *p;
	FILE	*fd;
	unsigned int i;
	char	line[1024], newFileName[80];
	char	*val;
	int	mangled;

/* See if file is cached */

	for (i = 0; i < 10; i++) {
		p = cache[i];
		if (p == NULL) {
			p = (struct IniCache *) malloc (sizeof (struct IniCache));
			p->filename = (char *) malloc (strlen (filename) + 1);
			strcpy (p->filename, filename);
			p->immediate_writes = 1;
			p->dirty = 0;
			p->num_lines = 0;
			p->array_size = 0;
			p->lines = NULL;
			forced_read = 1;
			cache[i] = p;
			break;
		}
		if (strcmp (filename, p->filename) == 0)
			break;
	}

/* Skip reading the ini file if appropriate */

	if (!forced_read) return (p);
	if (p->dirty) return (p);

/* Free the data if we've already read some in */

	for (i = 0; i < p->num_lines; i++) {
		free (p->lines[i]->keyword);
		free (p->lines[i]->value);
		free (p->lines[i]);
	}
	p->num_lines = 0;

/* Read the IniFile */

	mangleIniFileName (filename, newFileName, &mangled);
	fd = fopen (newFileName, "r");
	if (fd == NULL && mangled) fd = fopen (filename, "r");
	if (fd == NULL) return (p);

	while (fgets (line, sizeof (line), fd)) {
		if (line[strlen(line)-1] == '\n') line[strlen(line)-1] = 0;
		if (line[0] && line[strlen(line)-1] == '\r') line[strlen(line)-1] = 0;

/* Allocate and fill in a new line structure */

		growIniLineArray (p);
		i = p->num_lines++;
		p->lines[i] = (struct IniLine *) malloc (sizeof (struct IniLine));

/* Flag section headers */

		if (line[0] == '[') {
			char	*q;

			p->lines[i]->keyword = (char *) malloc (strlen (line) + 1);
			p->lines[i]->value = (char *) malloc (strlen (line) + 1);
			p->lines[i]->line_type = INI_LINE_HEADER;
			strcpy (p->lines[i]->value, line);
			strcpy (p->lines[i]->keyword, line+1);
			for (q = p->lines[i]->keyword; *q; q++)
				if (*q == ']') {
					*q = 0;
					break;
				}
		}

/* Save comment lines - any line that doesn't begin with a letter */

		else if ((line[0] < 'A' || line[0] > 'Z') &&
			 (line[0] < 'a' || line[0] > 'z')) {
			p->lines[i]->keyword = NULL;
			p->lines[i]->value = (char *) malloc (strlen (line) + 1);
			p->lines[i]->line_type = INI_LINE_COMMENT;
			strcpy (p->lines[i]->value, line);
		}

/* Otherwise, parse keyword=value lines */

		else {
			val = strchr (line, '=');
			if (val == NULL) {
				char	buf[1200];
				sprintf (buf, "Illegal line in %s: %s\n", newFileName, line);
				OutputSomewhere (MAIN_THREAD_NUM, buf);
				p->lines[i]->keyword = NULL;
				p->lines[i]->value = (char *) malloc (strlen (line) + 1);
				p->lines[i]->line_type = INI_LINE_COMMENT;
				strcpy (p->lines[i]->value, line);
			} else {
				*val++ = 0;
				p->lines[i]->keyword = (char *) malloc (strlen (line) + 1);
				p->lines[i]->value = (char *) malloc (strlen (val) + 1);
				p->lines[i]->line_type = INI_LINE_NORMAL;
				strcpy (p->lines[i]->keyword, line);
				strcpy (p->lines[i]->value, val);
			}
		}
	}
	fclose (fd);

	return (p);
}

/* Write a changed INI file to disk */

void writeIniFile (
	struct IniCache *p)
{
	int	fd, mangled;
	unsigned int j;
	char	buf[2000], newFileName[80];

/* Delay writing the file unless this INI file is written */
/* to immediately */

	if (!p->immediate_writes) {
		p->dirty = 1;
		return;
	}

/* Create and write out the INI file */

	mangleIniFileName (p->filename, newFileName, &mangled);
	fd = _open (newFileName, _O_CREAT | _O_TRUNC | _O_WRONLY | _O_TEXT, CREATE_FILE_ACCESS);
	if (fd < 0) return;
	for (j = 0; j < p->num_lines; j++) {
		if (p->lines[j]->line_type == INI_LINE_COMMENT) {
			strcpy (buf, p->lines[j]->value);
		} else if (p->lines[j]->line_type == INI_LINE_HEADER) {
			strcpy (buf, p->lines[j]->value);
		} else {
			strcpy (buf, p->lines[j]->keyword);
			strcat (buf, "=");
			strcat (buf, p->lines[j]->value);
		}
		strcat (buf, "\n");
		(void) _write (fd, buf, (unsigned int) strlen (buf));
	}
	p->dirty = 0;
	_close (fd);
	if (mangled) _unlink (p->filename);
}

/* Routines to help analyze a timed line in an INI file */

void parseTimeLine (
	const char **line,
	int	*start_day,
	int	*end_day,
	int	*start_time,
	int	*end_time)
{
	const char *p;

/* Get the days of the week, e.g. 1-5 */

	p = *line;
	*start_day = atoi (p); while (isdigit (*p)) p++;
	if (*p == '-') {
		p++;
		*end_day = atoi (p); while (isdigit (*p)) p++;
	} else
		*end_day = *start_day;

/* Now do time portion.  If none present, then assume the numbers we */
/* parsed above were times, not days of the week. */

	if (*p == '/')
		p++;
	else {
		p = *line;
		*start_day = 1;
		*end_day = 7;
	} 
	*start_time = atoi (p) * 60; while (isdigit (*p)) p++;
	if (*p == ':') {
		p++;
		*start_time += atoi (p); while (isdigit (*p)) p++;
	}
	if (*p == '-') p++;			/* Skip '-' */
	*end_time = atoi (p) * 60; while (isdigit (*p)) p++;
	if (*p == ':') {
		p++;
		*end_time += atoi (p); while (isdigit (*p)) p++;
	}

/* Return ptr to next time interval on the line */

	if (*p++ == ',') *line = p;
	else *line = NULL;
}

int analyzeTimeLine (
	const char *line,
	time_t	current_t,
	unsigned int *wakeup_time)
{
	struct tm *x;
	int	current_time;
	const char *p;
	int	day, start_day, end_day, start_time, end_time;
	int	full_start_time, full_end_time;
	int	wakeup_t, min_wakeup_t;

/* Break current time into a more easily maniupulated form */

	x = localtime (&current_t);
	current_time = (x->tm_wday ? x->tm_wday : 7) * 24 * 60;
	current_time += x->tm_hour * 60 + x->tm_min;

/* Process each interval on the line */

	p = line;
	min_wakeup_t = 0;
	while (p != NULL) {
		parseTimeLine (&p, &start_day, &end_day, &start_time, &end_time);

/* Treat each day in the range as a separate time interval to process */

		for (day = start_day; day <= end_day; day++) {

/* We allow end_time to be less than start_time.  We treat this as */
/* the next day.  Thus 15:00-01:00 means 3PM to 1AM the next day. */

			full_start_time = day * 24 * 60 + start_time;
			full_end_time = day * 24 * 60 + end_time;
			if (end_time < start_time) full_end_time += 24 * 60;

/* Is the current time in this interval? */

			if (current_time >= full_start_time &&
			    current_time < full_end_time)
				goto winner;

/* Now check for the really sick case, where end_time was less than */
/* start_time and we've wrapped from day 7 back to day 1 */

			if (end_time < start_time && day == 7 &&
			    current_time < full_end_time - 7 * 24 * 60)
				goto winner;

/* No, see if this start time should be our new wakeup time. */

			if (full_start_time >= current_time)
				wakeup_t = (full_start_time - current_time) * 60;
			else
				wakeup_t = (full_start_time + 7 * 24 * 60 - current_time) * 60;
			if (min_wakeup_t == 0 || min_wakeup_t > wakeup_t)
				min_wakeup_t = wakeup_t;
		}
	}

/* Current time was not in any of the intervals */

	*wakeup_time = min_wakeup_t;
	return (FALSE);

/* Current time is in this interval, compute the wakeup time */

winner:	wakeup_t = (full_end_time - current_time) * 60;

/* Also, look for a start time that matches the end time and replace */
/* the end time.  For example, if current time is 18:00 and the */
/* Time= entry is 0:00-8:00,17:00-24:00, then the */
/* end time of 24:00 should be replaced with 8:00 of the next day. */
/* Be sure not to infinite loop in this time entry: 0:00-8:00,8:00-24:00 */

	p = line;
	while (p != NULL && wakeup_t < 10 * 24 * 60) {
		parseTimeLine (&p, &start_day, &end_day, &start_time, &end_time);

/* Treat each day in the range as a separate time interval to process */

		for (day = start_day; day <= end_day; day++) {
			int	this_full_start_time, this_full_end_time;

/* If this start time is the same as the winning end time, then set the new */
/* wakeup time to be the end of this interval.  Be sure to handle the tricky */
/* wrap around that occurs when end_time < start_time. */

			this_full_start_time = day * 24 * 60 + start_time;
			if (this_full_start_time != full_end_time &&
			    this_full_start_time != full_end_time - 7 * 24 * 60) continue;

			this_full_end_time = day * 24 * 60 + end_time;
			if (end_time < start_time) this_full_end_time += 24 * 60;
			wakeup_t += (this_full_end_time - this_full_start_time) * 60;
			full_end_time = this_full_end_time;
			p = line;
			break;
		}
	}

/* Return indicator that current time was covered by one of the intervals */

	*wakeup_time = wakeup_t + 1;
	return (TRUE);
}

/* INI file values can contain be conditional based on the day of the */
/* week and the time of day.  For example, this INI file value is */
/* file has different properties during the work week and weekend. */
/*	Priority=1 during 1-5/8:30-17:30 else 5			*/

void parse_timed_ini_value (
	const char *line,		/* INI value line to analyze */
	unsigned int *start_offset,	/* Returned start offset */
	unsigned int *len,		/* Returned length */
	unsigned int *seconds_valid)	/* Returned length of time */
					/* value is good for */
{
	time_t	current_time;
	const char *rest_of_line, *during_clause, *else_clause;
	unsigned int min_wakeup_time, wakeup_time;

/* Get the current time - so that we compare each timed section */
/* with the same current_time value */

	time (&current_time);

/* Loop processing each timed section in the line */

	rest_of_line = line;
	min_wakeup_time = 0;
	for ( ; ; ) {

/* If we don't see a "during" clause, then either there are no timed sections */
/* or we've reached the final else clause.  Return the else clause value. */

		during_clause = strstr (rest_of_line, " during ");
		if (during_clause == NULL) {
			*start_offset = (unsigned int) (rest_of_line - line);
			*len = (unsigned int) strlen (rest_of_line);
			*seconds_valid = min_wakeup_time;
			break;
		}

/* We've got a timed section, see if the current time is */
/* within this timed section. */

		if (analyzeTimeLine (during_clause+8, current_time, &wakeup_time)) {
			*start_offset = (unsigned int) (rest_of_line - line);
			*len = (unsigned int) (during_clause - rest_of_line);
			*seconds_valid = wakeup_time;
			break;
		}

/* We're not in this timed section, remember which timed section */
/* will come into effect first.  This will be the end time of the "else" */
/* section. */

		if (min_wakeup_time == 0 || wakeup_time < min_wakeup_time)
			min_wakeup_time = wakeup_time;

/* Move on to the next timed section. */

		else_clause = strstr (during_clause, " else ");
		if (else_clause != NULL) rest_of_line = else_clause + 6;
		else rest_of_line += strlen (rest_of_line);
	}
}

/* Utility routines used in copying an INI setting value */

void truncated_strcpy_with_len (
	char	*buf,
	unsigned int bufsize,
	const char *val,
	unsigned int valsize)
{
	if (valsize >= bufsize) valsize = bufsize - 1;
	memcpy (buf, val, valsize);
	buf[valsize] = 0;
}

void truncated_strcpy (
	char	*buf,
	unsigned int bufsize,
	const char *val)
{
	truncated_strcpy_with_len (buf, bufsize, val, (unsigned int) strlen (val));
}

/* Get a keyword's value from a specific section of the INI file. */
/* Do not process any timed sections. */

const char *IniSectionGetStringRaw (
	const char *filename,
	const char *section,
	const char *keyword)
{
	struct IniCache *p;
	unsigned int i;
	const char *retval;

/* Open ini file */

	gwmutex_lock (&INI_MUTEX);
	p = openIniFile (filename, 0);

/* Skip to the correct section */

	i = 0;
	if (section != NULL) {
		for ( ; i < p->num_lines; i++) {
			if (p->lines[i]->line_type == INI_LINE_HEADER &&
			    _stricmp (section, p->lines[i]->keyword) == 0) {
				i++;
				break;
			}
		}
	}

/* Look for the keyword within this section */

	for ( ; ; i++) {
		if (i == p->num_lines ||
		    p->lines[i]->line_type == INI_LINE_HEADER) {
			retval = NULL;
			break;
		}
		if (p->lines[i]->line_type == INI_LINE_NORMAL &&
		    _stricmp (keyword, p->lines[i]->keyword) == 0) {
			retval = p->lines[i]->value;
			break;
		}
	}

/* Unlock and return */

	gwmutex_unlock (&INI_MUTEX);
	return (retval);
}

/* Get a keyword's value from a specific section of the INI file. */
/* Return length of time this timed INI setting is good for. */

void IniSectionGetTimedString (
	const char *filename,
	const char *section,
	const char *keyword,
	char	*val,
	unsigned int val_bufsize,
	const char *default_val,
	unsigned int *seconds)
{
	const char *p;
	unsigned int start, len;

/* Lookup the keyword */

	p = IniSectionGetStringRaw (filename, section, keyword);

/* If we found the keyword in the INI file, then */
/* support different return values based on the time of day. */

	if (p != NULL) {
		parse_timed_ini_value (p, &start, &len, seconds);
		if (len) {
			truncated_strcpy_with_len (val, val_bufsize, p+start, len);
			return;
		}
	} else {
		*seconds = 0;
	}

/* Copy the default value to the caller's buffer */

	if (default_val)
		truncated_strcpy (val, val_bufsize, default_val);
	else
		val[0] = 0;
}

/* Get a keyword's value from a specific section of the INI file. */

void IniSectionGetString (
	const char *filename,
	const char *section,
	const char *keyword,
	char	*val,
	unsigned int val_bufsize,
	const char *default_val)
{
	unsigned int seconds;
	IniSectionGetTimedString (filename, section, keyword, val, val_bufsize, default_val, &seconds);
}

long IniSectionGetTimedInt (
	const char *filename,
	const char *section,
	const char *keyword,
	long	default_val,
	unsigned int *seconds)
{
	char	buf[20], defval[20];
	sprintf (defval, "%ld", default_val);
	IniSectionGetTimedString (filename, section, keyword, buf, 20, defval, seconds);
	return (atol (buf));
}

long IniSectionGetInt (
	const char *filename,
	const char *section,
	const char *keyword,
	long	default_val)
{
	unsigned int seconds;
	return (IniSectionGetTimedInt (filename, section, keyword, default_val, &seconds));
}

void IniSectionWriteString (
	const char *filename,
	const char *section,
	const char *keyword,
	const char *val)
{
	struct IniCache *p;
	unsigned int i, j, insertion_point;

/* Open ini file */

	gwmutex_lock (&INI_MUTEX);
	p = openIniFile (filename, 0);

/* Skip to the correct section.  If the section does not exist, create it */

	i = 0;
	if (section != NULL) {
		for ( ; ; i++) {
			if (i == p->num_lines) {
				if (val == NULL) goto done;
				growIniLineArray (p);
				p->lines[i] = (struct IniLine *)
					malloc (sizeof (struct IniLine));
				p->lines[i]->line_type = INI_LINE_COMMENT;
				p->lines[i]->keyword = NULL;
				p->lines[i]->value = (char *) malloc (1);
				p->lines[i]->value[0] = 0;
				p->num_lines++;
				i++;
				growIniLineArray (p);
				p->lines[i] = (struct IniLine *)
					malloc (sizeof (struct IniLine));
				p->lines[i]->line_type = INI_LINE_HEADER;
				p->lines[i]->keyword = (char *)
					malloc (strlen (section) + 1);
				strcpy (p->lines[i]->keyword, section);
				p->lines[i]->value = (char *)
					malloc (strlen (section) + 3);
				sprintf (p->lines[i]->value, "[%s]", section);
				p->num_lines++;
				i++;
				break;
			}
			if (p->lines[i]->line_type == INI_LINE_HEADER &&
			    _stricmp (section, p->lines[i]->keyword) == 0) {
				i++;
				break;
			}
		}
	}

/* Look for the keyword within this section */

	insertion_point = i;
	for ( ; ; i++) {
		if (i == p->num_lines ||
		    p->lines[i]->line_type == INI_LINE_HEADER ||
		    (p->lines[i]->line_type != INI_LINE_COMMENT &&
		     _stricmp (p->lines[i]->keyword, "Time") == 0)) {

/* Ignore request if we are deleting line */

			if (val == NULL) goto done;

/* Make sure the line array has room for the new line */

			growIniLineArray (p);

/* Shuffle entries down to make room for this entry */

			i = insertion_point;
			for (j = p->num_lines; j > i; j--)
				p->lines[j] = p->lines[j-1];

/* Allocate and fill in a new line structure */

			p->lines[i] = (struct IniLine *) malloc (sizeof (struct IniLine));
			p->lines[i]->line_type = INI_LINE_NORMAL;
			p->lines[i]->keyword = (char *) malloc (strlen (keyword) + 1);
			strcpy (p->lines[i]->keyword, keyword);
			p->lines[i]->value = NULL;
			p->num_lines++;
			break;
		}

/* If this is not a blank line, then if we need to insert a new line, */
/* insert it after this line.  In other words, insert new entries before */
/* any blank lines at the end of a section */

		if (p->lines[i]->line_type != INI_LINE_COMMENT ||
		    p->lines[i]->value[0]) {
			insertion_point = i + 1;
		}

/* If this is the keyword we are looking for, then we will replace the */
/* value if it has changed. */

		if (p->lines[i]->line_type == INI_LINE_NORMAL &&
		    _stricmp (keyword, p->lines[i]->keyword) == 0) {
			if (val != NULL &&
			    strcmp (val, p->lines[i]->value) == 0) goto done;
			break;
		}
	}

/* Delete the line if requested */

	if (val == NULL) {

/* Free the data associated with the given line */

		free (p->lines[i]->keyword);
		free (p->lines[i]->value);
		free (p->lines[i]);

/* Delete the line from the lines array */

		for (i++; i < p->num_lines; i++) p->lines[i-1] = p->lines[i];
		p->num_lines--;
	}

/* Replace the value associated with the keyword */

	else {
		free (p->lines[i]->value);
		p->lines[i]->value = (char *) malloc (strlen (val) + 1);
		strcpy (p->lines[i]->value, val);
	}

/* Write the INI file back to disk */

	writeIniFile (p);

/* Unlock and return */

done:	gwmutex_unlock (&INI_MUTEX);
}

void IniSectionWriteInt (
	const char *filename,
	const char *section,
	const char *keyword,
	long	val)
{
	char	buf[20];
	sprintf (buf, "%ld", val);
	IniSectionWriteString (filename, section, keyword, buf);
}

void IniSectionWriteFloat (
	const char *filename,
	const char *section,
	const char *keyword,
	float	val)
{
	char	buf[20];
	sprintf (buf, "%f", val);
	IniSectionWriteString (filename, section, keyword, buf);
}

/* Shorthand routines for reading and writing from the global section */

void IniGetTimedString (
	const char *filename,
	const char *keyword,
	char	*val,
	unsigned int val_bufsize,
	const char *default_val,
	unsigned int *seconds)
{
	IniSectionGetTimedString (filename, NULL, keyword, val, val_bufsize, default_val, seconds);
}

void IniGetString (
	const char *filename,
	const char *keyword,
	char	*val,
	unsigned int val_bufsize,
	const char *default_val)
{
	IniSectionGetString (filename, NULL, keyword, val, val_bufsize, default_val);
}

long IniGetTimedInt (
	const char *filename,
	const char *keyword,
	long	default_val,
	unsigned int *seconds)	     
{
	return (IniSectionGetTimedInt (filename, NULL, keyword, default_val, seconds));
}

long IniGetInt (
	const char *filename,
	const char *keyword,
	long	default_val)
{
	return (IniSectionGetInt (filename, NULL, keyword, default_val));
}

/* Write a string to the INI file. */

void IniWriteString (
	const char *filename,
	const char *keyword,
	const char *val)
{
	IniSectionWriteString (filename, NULL, keyword, val);
}

/* Write an integer to the INI file. */

void IniWriteInt (
	const char *filename,
	const char *keyword,
	long	val)
{
	IniSectionWriteInt (filename, NULL, keyword, val);
}

/* Write a float to the INI file. */

void IniWriteFloat (
	const char *filename,
	const char *keyword,
	float	val)
{
	IniSectionWriteFloat (filename, NULL, keyword, val);
}

/* Reread the INI file. */

void IniFileReread (
	const char *filename)
{
	struct IniCache *p;
	gwmutex_lock (&INI_MUTEX);
	p = openIniFile (filename, 1);
	gwmutex_unlock (&INI_MUTEX);
}

//
//bug - do we want to offer two new routines:  immediate_writes_on and
// immediate_writes_off?  This would speed up bulk changes to the INI file.
//

/****************************************************************************/
/*               Utility routines to work with ".add" files                 */
/****************************************************************************/

/* See if an "add file" file exists.  An add file lets the user create a */
/* prime.add, local.add, or worktodo.add file to overwrite/append options */
/* or append work while the program is running.  This is especially handy */
/* for workstations that are not physically accessible but there file */
/* systems are by network.  If this feature was not available, the only */
/* safe method for updating would be to stop the program, edit the .ini */
/* file, and restart the program. */

int addFileExists (void)
{
	char	filename[80];
	char	*dot;

	strcpy (filename, INI_FILE);
	dot = strrchr (filename, '.');
	if (dot != NULL) {
		strcpy (dot, ".add");
		if (fileExists (filename)) return (TRUE);
		strcpy (dot, ".add.txt");
		if (fileExists (filename)) return (TRUE);
	}
	strcpy (filename, LOCALINI_FILE);
	dot = strrchr (filename, '.');
	if (dot != NULL) {
		strcpy (dot, ".add");
		if (fileExists (filename)) return (TRUE);
		strcpy (dot, ".add.txt");
		if (fileExists (filename)) return (TRUE);
	}
	strcpy (filename, WORKTODO_FILE);
	dot = strrchr (filename, '.');
	if (dot != NULL) {
		strcpy (dot, ".add");
		if (fileExists (filename)) return (TRUE);
		strcpy (dot, ".add.txt");
		if (fileExists (filename)) return (TRUE);
	}
	return (FALSE);
}

/* Merge one "add file" into an ini file.  Assumes the ini file has been */
/* freshly re-read from disk.  This also is used to copy an old primenet.ini */
/* to a [PrimeNet] section of prime.ini */

void iniAddFileMerge (
	char	*ini_filename,
	char	*add_filename,
	char	*section_to_copy_to)
{
	struct IniCache *p, *q;
	char	*section;
	unsigned int j;

/* Obtain a lock so that only one thread adds to the INI file.  We will release */
/* lock once we've deleted the add file. */
	
	gwmutex_lock (&INI_ADD_MUTEX);

/* Open ini files */

	p = openIniFile (ini_filename, 0);
	q = openIniFile (add_filename, 1);

/* Save up all the writes */

	p->immediate_writes = FALSE;

/* Loop through all the lines in the add file, adding them to the */
/* base ini file */

	section = section_to_copy_to;
	for (j = 0; j < q->num_lines; j++) {
		if (q->lines[j]->line_type == INI_LINE_HEADER) {
			if (section_to_copy_to == NULL)
				section = q->lines[j]->keyword;
		}
		else if (q->lines[j]->line_type != INI_LINE_COMMENT)
			IniSectionWriteString (ini_filename, section,
						q->lines[j]->keyword,
						q->lines[j]->value);
	}

/* Output all the saved up writes */

	p->immediate_writes = TRUE;
	writeIniFile (p);
	
/* Delete the add file */

	_unlink (add_filename);

/* Unlock and return */

	gwmutex_unlock (&INI_ADD_MUTEX);
}

/* Merge all INI ".add files" into their corresponding base files */

void incorporateIniAddFiles (void)
{
	char	filename[80];
	char	*dot;

/* Merge additions to prime.ini */

	strcpy (filename, INI_FILE);
	dot = strrchr (filename, '.');
	if (dot != NULL) {
		strcpy (dot, ".add");
		if (fileExists (filename))
			iniAddFileMerge (INI_FILE, filename, NULL);
		strcpy (dot, ".add.txt");
		if (fileExists (filename))
			iniAddFileMerge (INI_FILE, filename, NULL);
	}

/* Merge additions to local.ini */

	strcpy (filename, LOCALINI_FILE);
	dot = strrchr (filename, '.');
	if (dot != NULL) {
		strcpy (dot, ".add");
		if (fileExists (filename))
			iniAddFileMerge (LOCALINI_FILE, filename, NULL);
		strcpy (dot, ".add.txt");
		if (fileExists (filename))
			iniAddFileMerge (LOCALINI_FILE, filename, NULL);
	}
}

/* Merge optional worktodo.add file into their worktodo.ini file */

int incorporateWorkToDoAddFile (void)
{
static	int	worktodo_add_disabled = FALSE;
	char	filename[80];
	char	*dot;
	int	rc;
	FILE	*fd;
	unsigned int tnum;
	char	line[2048];

/* If add files have been disabled (see below) then we're all done */

	if (worktodo_add_disabled) return (0);

/* Merge additions to worktodo.ini */

	strcpy (filename, WORKTODO_FILE);
	dot = strrchr (filename, '.');
	if (dot == NULL) return (0);

/* Open the worktodo.add file, it is OK if this file does not exist. */

	strcpy (dot, ".add");
	fd = fopen (filename, "r");
	if (fd == NULL) {
		strcpy (dot, ".add.txt");
		fd = fopen (filename, "r");
		if (fd == NULL) return (0);
	}

/* As an ugly kludge, we append lines from worktodo.add as comments in the */
/* in-memory version of worktodo.ini.  Later we will write worktodo.ini to */
/* disk and reprocess it entirely.  Loop processing each worktodo.add line */

	gwmutex_lock (&WORKTODO_MUTEX);
	tnum = 0;
	while (fgets (line, sizeof (line), fd)) {
		struct work_unit *w;

/* Remove trailing CRLFs */

		if (line[strlen(line)-1] == '\n')
			line[strlen(line)-1] = 0;
		if (line[0] && line[strlen(line)-1] == '\r')
			line[strlen(line)-1] = 0;
		if (line[0] == 0) continue;

/* If this is a section header find the matching section header in */
/* worktodo.ini.  If no match is found, add this to the first empty thread */
/* or the very last thread */

		if (line[0] == '[') {
			struct work_unit *w;
			for (tnum = 0; ; tnum++) {
				w = WORK_UNITS[tnum].first;
				if (w == NULL) {
					if (tnum) break;
				} else {
					if (w->work_type == WORK_NONE &&
					    _stricmp (w->comment, line) == 0)
						break;
				}
				if (tnum == NUM_WORKER_THREADS - 1) {
					w = NULL;
					break;
				}
			}
			if (w != NULL) continue;
		}

/* Allocate a work unit structure */

		w = (struct work_unit *) malloc (sizeof (struct work_unit));
		if (w == NULL) goto nomem;
		memset (w, 0, sizeof (struct work_unit));

/* Save new line as a comment.  It will be properly parsed when we re-read */
/* the worktodo.ini file. */

		w->work_type = WORK_NONE;
		w->comment = (char *) malloc (strlen (line) + 1);
		if (w->comment == NULL) goto nomem;
		strcpy (w->comment, line);

/* Grow the work_unit array if necessary and add this entry */

		rc = addToWorkUnitArray (tnum, w, FALSE);
		if (rc) goto retrc;
	}

/* Close the file, free the lock and return success */

	fclose (fd);
	gwmutex_unlock (&WORKTODO_MUTEX);

/* Write the combined worktodo.ini file */

	WORKTODO_CHANGED = TRUE;
	rc = writeWorkToDoFile (TRUE);
	if (rc) return (rc);

/* Delete the worktodo.add file.  If file exists after we tried to delete it */
/* then permanently disable worktodo.add processing (to avoid an infinite */
/* loop growing and growing the worktodo file with redundant work! */

	_unlink (filename);
	if (fileExists (filename)) {
		OutputBoth (MAIN_THREAD_NUM,
			    "ERROR:  Can't delete worktodo.add file\n");
		worktodo_add_disabled = TRUE;
	}

/* Now reprocess the combined and freshly written worktodo.ini file */

	return (readWorkToDoFile ());

/* Handle an error during the reading of the add file */

nomem:	rc = OutOfMemory (MAIN_THREAD_NUM);
retrc:	fclose (fd);
	gwmutex_unlock (&WORKTODO_MUTEX);
	return (rc);
}

/****************************************************************************/
/*          Utility routines to work with per-thread options (PTO)          */
/****************************************************************************/

void PTOGetAll (
	char	*ini_filename,		/* Ini file containing the options */
	char	*keyword,		/* Ini file keyword */
	unsigned int *array,		/* Options array */
	unsigned int def_val)		/* Default value */
{
	int	i, global_val;
	char	section_name[32];

/* Copy the global option setting to the entire array */

	global_val = IniGetInt (ini_filename, keyword, -1);
	for (i = 0; i < MAX_NUM_WORKER_THREADS; i++) {
		if (global_val == -1) array[i] = def_val;
		else array[i] = global_val;
	}

/* Now look for any section specific overrides */

	for (i = 0; i < MAX_NUM_WORKER_THREADS; i++) {
		sprintf (section_name, "Worker #%d", i+1);
		array[i] = IniSectionGetInt (ini_filename, section_name,
					     keyword, array[i]);
	}
}

void PTOSetAll (
	char	*ini_filename,		/* Ini file containing the options */
	char	*keyword,		/* Ini file keyword */
	char	*shadow_keyword,	/* Ini file keyword for value the */
					/* server has stored for this option */
	unsigned int *array,		/* Options array */
	unsigned int new_val)		/* New option value */
{
	int	i;
	char	section_name[32];

/* Copy the global option setting to the entire array */

	IniWriteInt (ini_filename, keyword, new_val);
	for (i = 0; i < MAX_NUM_WORKER_THREADS; i++) {
		array[i] = new_val;
	}
	if (shadow_keyword != NULL)
		IniWriteInt (LOCALINI_FILE, shadow_keyword, new_val);

/* Now clear all section specific overrides */

	for (i = 0; i < MAX_NUM_WORKER_THREADS; i++) {
		sprintf (section_name, "Worker #%d", i+1);
		IniSectionWriteString (ini_filename, section_name, keyword, NULL);
		if (shadow_keyword != NULL)
			IniSectionWriteString (LOCALINI_FILE, section_name,
					       shadow_keyword, NULL);
	}
}

void PTOSetOne (
	char	*ini_filename,		/* Ini file containing the options */
	char	*keyword,		/* Ini file keyword */
	char	*shadow_keyword,	/* Ini file keyword for value the */
					/* server has stored for this option */
	unsigned int *array,		/* Options array */
	int	tnum,			/* Thread number */
	unsigned int new_val)		/* New option value */
{
	char	section_name[32];

/* Will changing this option cause us to switch from one global setting */
/* to individual settings on each thread? */

	if (PTOIsGlobalOption (array)) {
		int	i;

/* If option has not changed, then do nothing */

		if (array[tnum] == new_val) return;

/* Delete the global option and set the thread specific option for */
/* each thread. */

		IniWriteString (ini_filename, keyword, NULL);
		if (shadow_keyword != NULL)
			IniWriteString (LOCALINI_FILE, shadow_keyword, NULL);
		for (i = 0; i < (int) NUM_WORKER_THREADS; i++) {
			sprintf (section_name, "Worker #%d", i+1);
			IniSectionWriteInt (ini_filename, section_name,
					    keyword, array[i]);
			if (shadow_keyword != NULL)
				IniSectionWriteInt (LOCALINI_FILE, section_name,
						    shadow_keyword, array[i]);
		}
	}

/* Set the option for just this one thread */

	array[tnum] = new_val;
	sprintf (section_name, "Worker #%d", tnum+1);
	IniSectionWriteInt (ini_filename, section_name, keyword, new_val);
	if (shadow_keyword != NULL)
		IniSectionWriteInt (LOCALINI_FILE, section_name,
				    shadow_keyword, new_val);
}

int PTOIsGlobalOption (
	unsigned int *array)		/* Options array */
{
	int	i;

	for (i = 1; i < (int) NUM_WORKER_THREADS; i++)
		if (array[i-1] != array[i]) return (FALSE);
	return (TRUE);
}

int PTOHasOptionChanged (
	char	*shadow_keyword,	/* Ini file keyword for value the */
					/* server has stored for this option */
	unsigned int *array,		/* Options array */
	int	tnum)			/* Thread number */
{
	char	section_name[32];

/* If this is the thread number that tells us we are dealing global options */
/* then return TRUE if this is a globally set option and it has changed. */

	if (tnum == -1)
		return (PTOIsGlobalOption (array) &&
			array[0] != IniGetInt (LOCALINI_FILE, shadow_keyword, -1));

/* Otherwise, return TRUE if this is not a globally set option and it has changed. */

	else {
		sprintf (section_name, "Worker #%d", tnum+1);
		return (!PTOIsGlobalOption (array) &&
			array[tnum] != IniSectionGetInt (LOCALINI_FILE, section_name, shadow_keyword, -1));
	}
}


/****************************************************************************/
/*                Utility routines to output messages                       */
/****************************************************************************/

/* Output string to screen or results file */

void OutputSomewhere (
	int	thread_num,
	const char *buf)
{
	if (NO_GUI) writeResults (buf);
	else OutputStr (thread_num, buf);
}

/* Output string to both the screen and results file */

EXTERNC
void OutputBoth (
	int	thread_num,
	const char *buf)
{
	OutputStr (thread_num, buf);
	writeResults (buf);
}

/* Output string to screen.  Prefix it with an optional timestamp. */

void OutputStr (
	int	thread_num,
	const char *buf)
{

/* Grab a lock so that multiple threads cannot output at the same time. */
/* Simultaneous output would be real bad on OS's that output all this data */
/* to a single screen/window (such as a Linux command line implementation). */

	gwmutex_lock (&OUTPUT_MUTEX);

/* Format a timestamp.  Prefix every line in the input buffer with */
/* the timestamp. */

	if (TIMESTAMPING) {
		time_t	this_time;
		char	tmpbuf[200], fmtbuf[40];

		time (&this_time);
		strcpy (tmpbuf, ctime (&this_time)+4);

		/* Eliminate seconds and year or just year */
		if (TIMESTAMPING & 1) tmpbuf[12] = 0;
		else tmpbuf[15] = 0;

		/* Eliminate date or zero-suppress day */
		if (TIMESTAMPING >= 3)
			safe_strcpy (tmpbuf, tmpbuf+7);
		else if (tmpbuf[4] == '0' || tmpbuf[4] == ' ')
			safe_strcpy (tmpbuf+4, tmpbuf+5);

		sprintf (fmtbuf, "[%s] ", tmpbuf);

/* Output the prefix for every line in the buffer */ 

		do {
			const char *eol;

			eol = strchr (buf, '\n');
			if (eol != NULL) eol++;
			else eol = buf + strlen (buf);
			RealOutputStr (thread_num, fmtbuf);
			while (buf != eol) {
				int	len;
				len = (int) (eol - buf);
				if (len >= sizeof (tmpbuf))
					len = sizeof (tmpbuf) - 1;
				memcpy (tmpbuf, buf, len);
				tmpbuf[len] = 0;
				RealOutputStr (thread_num, tmpbuf);
				buf += len;
			}
		} while (*buf);
	}

/* No timestamp - just output the possibly multi-line buffer. */

	else
		RealOutputStr (thread_num, buf);

/* Free the lock and return */

	gwmutex_unlock (&OUTPUT_MUTEX);
}

/* Output string to screen without prefixing a timestamp.  Only used when */
/* outputting a line in chunks (we do this when benchmarking). */

void OutputStrNoTimeStamp (
	int	thread_num,
	const char *buf)
{
	gwmutex_lock (&OUTPUT_MUTEX);
	RealOutputStr (thread_num, buf);
	gwmutex_unlock (&OUTPUT_MUTEX);
}

/* Output an out-of-memory error */

int OutOfMemory (
	int thread_num)
{
	OutputStr (thread_num, "Out of memory!\n");
	return (STOP_OUT_OF_MEM);
}

/****************************************************************************/
/*               Routines to process worktodo.ini files                     */
/****************************************************************************/

/* The worktodo structures are accessed by worker threds, the communication */
/* thread, the timer threads, and the GUI thread.  We must be VERY CAREFUL */
/* with our locking scheme to make sure there are no crashes or hangs. */

/* The simplest scheme just locks on entry to each of these routines. */
/* This solves some basic problems such as two threads writing the worktodo */
/* file at the same time.  However, more difficult problems still exist. */
/* For example, a worker thread could delete the work unit while the */
/* GUI (Test/Status) or comm thread (sending completion dates) try to */
/* has a pointer to the work unit structure.  A dereference after the */
/* delete will crash. */

/* The scheme I came up with increments a per work unit counter that */
/* indicates the work unit is in use.  Deletes are not allowed while the */
/* work unit is in use. */

/* We complicate this scheme by differentiating between in-use by a worker */
/* thread (it could take a very long time to decrement the counter) and */
/* access by the comm or GUI threads (these should be very quick accesses). */

/* Finally, the comm thread has a work-unit in use when it sends a message */
/* to the server.  If this hangs, the timer code must decrement the in-use */
/* counter and kill the thread (or leave it in a hung state). */


/* Count commas in a string */

unsigned int countCommas (
	char *p)
{
	unsigned int cnt;
	for (cnt = 0; ; cnt++) {
		p = strchr (p, ',');
		if (p == NULL) break;
		p++;
	}
	return (cnt);
}

/* Do some more initialization of work_unit fields.  These values do */
/* not appear in the worktodo.ini file, but need initializing in a */
/* common place. */

void auxiliaryWorkUnitInit (
	struct work_unit *w)
{

/* Compute factor_to if not set already */

	if ((w->work_type == WORK_FACTOR ||
	     w->work_type == WORK_TEST ||
	     w->work_type == WORK_DBLCHK) &&
	    w->factor_to == 0.0)
		w->factor_to = factorLimit (w);

/* Initialize the number of LL tests saved */

	if (w->work_type == WORK_TEST) w->tests_saved = 2.0;
	if (w->work_type == WORK_DBLCHK) w->tests_saved = 1.0;

/* Guard against wild tests_saved values.  Huge values will cause guess_pminus1_bounds */
/* to run for a very long time. */

	if (w->tests_saved > 10) w->tests_saved = 10;
}

/* Fill in a work unit's stage and percentage complete based on any */
/* save files. */

void pct_complete_from_savefile (
	struct work_unit *w)
{
	int	fd, res;
	unsigned long version;
	char	filename[32];

/* Generate the save file name */

	tempFileName (w, filename);

/* See if there is an intermediate file.  If there is read it to get our */
/* stage and percent complete.  This is a little tricky for LL and PRP tests */
/* which could have 3 different types of save files (LL/PRP, P-1, or trial */
/* factoring). */

	for ( ; ; ) {
		fd = _open (filename, _O_BINARY | _O_RDONLY);
		if (fd <= 0) {
			if ((w->work_type == WORK_TEST ||
			     w->work_type == WORK_DBLCHK ||
			     w->work_type == WORK_PRP) && filename[0] == 'p') {
				filename[0] = 'm';
				continue;
			}
			if ((w->work_type == WORK_TEST ||
			     w->work_type == WORK_DBLCHK) && filename[0] == 'm') {
				filename[0] = 'f';
				continue;
			}
			break;
		}

/* Read the header */

		res = read_header (fd, &version, w, NULL);
		_close (fd);
		if (!res) break;

/* We've successfully read a save file header, return now that the work */
/* unit's stage and pct_complete fields have been filled in. */

		return;
	}

/* We have not started working on this item */

	w->stage[0] = 0;
	w->pct_complete = 0.0;
}

/* Add a work_unit to the work_unit array.  Grow the work_unit */
/* array if necessary */

int addToWorkUnitArray (
	unsigned int tnum,	/* Thread number that will run work unit */
	struct work_unit *w,	/* Work unit to add to array */
	int	add_to_end)	/* TRUE if we are to add to end of array */
{
	struct work_unit *insertion_point;

/* If add_to_end is set, then we are reading the worktodo.ini file. */
/* Add the entry to the end of the array */

	if (add_to_end)
		insertion_point = WORK_UNITS[tnum].last;

/* Otherwise, add ADVANCEDTEST work units to the front of the array after */
/* any comment lines. */

	else if (w->work_type == WORK_ADVANCEDTEST) {
		if (WORK_UNITS[tnum].first != NULL &&
		    WORK_UNITS[tnum].first->work_type == WORK_NONE)
			for (insertion_point = WORK_UNITS[tnum].first;
			     insertion_point->next != NULL;
			     insertion_point = insertion_point->next) {
				if (insertion_point->next->work_type != WORK_NONE)
					break;
			}
		else
			insertion_point = NULL;
	}

/* Add all other work units before the last blank line. */

	else {
		for (insertion_point = WORK_UNITS[tnum].last;
		     insertion_point != NULL;
		     insertion_point = insertion_point->prev) {
			if (insertion_point->work_type != WORK_NONE ||
			    insertion_point->comment[0]) break;
		}
	}

/* Now insert the work unit after the insertion point */

	w->prev = insertion_point;
	if (insertion_point == NULL) {
		w->next = WORK_UNITS[tnum].first;
		WORK_UNITS[tnum].first = w;
	} else {
		w->next = insertion_point->next;
		insertion_point->next = w;
	}
	if (w->next == NULL)
		WORK_UNITS[tnum].last = w;
	else
		w->next->prev = w;

/* Bump count of valid work lines and wake up thread if it is waiting */
/* for work to do. */

	if (w->work_type != WORK_NONE) {
		WORKTODO_COUNT++;
		restart_one_waiting_worker (tnum, RESTART_WORK_AVAILABLE);
	}

/* Return success */

	return (0);
}

/* Read the entire worktodo.ini file into memory.  Return error_code */
/* if we have a memory or file I/O error. */

int readWorkToDoFile (void)
{
	FILE	*fd;
	unsigned int tnum, i, linenum;
	int	rc, mangled;
	char	line[16384], newFileName[80];

/* Grab the lock so that comm thread cannot try to add work units while */
/* file is being read in. */

	for (i = 1; ; i++) {
		gwmutex_lock (&WORKTODO_MUTEX);

/* Make sure no other threads are accessing work units right now. */
/* There should be no worker threads active so any use should be short-lived. */

		if (WORKTODO_IN_USE_COUNT == 0 && !WORKTODO_CHANGED) break;
		gwmutex_unlock (&WORKTODO_MUTEX);
		if (i <= 10) {
			Sleep (50);
			continue;
		}

/* Uh oh, the lock hasn't been released after half-a-second.  This happens processing large */
/* worktodo.txt files in communicateWithServer (see James Heinrich's complaints in 26.4 thread). */
/* As a workaround, we'll simply not re-read the worktodo.txt file now.  We only reread the file */
/* to pick up any manual edits that may have taken place since the last time worktodo.txt was */
/* read in (and to process worktodo.add).  Hopefully the comm-with-server thread will finish up */
/* and we can successfully re-read the worktodo.txt file at a later time. */

		return (0);
	}

/* Clear file needs writing flag and count of worktodo lines */

	WORKTODO_CHANGED = FALSE;
	WORKTODO_COUNT = 0;

/* Free old work_units for each worker thread. */
/* We sometimes reread the worktodo.ini file in case the user */
/* manually edits the file while the program is running. */

	for (tnum = 0; tnum < MAX_NUM_WORKER_THREADS; tnum++) {
		struct work_unit *w, *next_w;
		for (w = WORK_UNITS[tnum].first; w != NULL; w = next_w) {
			next_w = w->next;
			free (w->known_factors);
			free (w->comment);
			free (w);
		}
		WORK_UNITS[tnum].first = NULL;
		WORK_UNITS[tnum].last = NULL;
	}

/* Read the lines of the work file.  It is OK if the worktodo.ini file */
/* does not exist. */

	mangleIniFileName (WORKTODO_FILE, newFileName, &mangled);
	fd = fopen (newFileName, "r");
	if (fd == NULL && mangled) fd = fopen (WORKTODO_FILE, "r");
	if (fd == NULL) goto done;

	tnum = 0;
	linenum = 0;
	while (fgets (line, sizeof (line), fd)) {
	    struct work_unit *w;
	    char keyword[20];
	    char *value;

/* Remove trailing CRLFs */

	    if (line[strlen(line)-1] == '\n') line[strlen(line)-1] = 0;
	    if (line[0] && line[strlen(line)-1] == '\r') line[strlen(line)-1] = 0;
	    linenum++;

/* Allocate a work unit structure */

	    w = (struct work_unit *) malloc (sizeof (struct work_unit));
	    if (w == NULL) goto nomem;
	    memset (w, 0, sizeof (struct work_unit));

/* A section header precedes each worker thread's work units.  The first */
/* section need not be preceeded by a section header. */

	    if (line[0] == '[' && linenum > 1) {
		tnum++;
		if (tnum == NUM_WORKER_THREADS) {
		    char	buf[100];
		    sprintf (buf,
			     "Too many sections in worktodo.txt.  Line #%u\n",
			     linenum);
		    OutputSomewhere (MAIN_THREAD_NUM, buf);
		}
		if (tnum >= MAX_NUM_WORKER_THREADS)
			tnum = MAX_NUM_WORKER_THREADS - 1;
	    }

/* All lines other than keyword=value are saved as comment lines. */

	    if (((line[0] < 'A' || line[0] > 'Z') &&
		 (line[0] < 'a' || line[0] > 'z'))) {
comment:	w->work_type = WORK_NONE;
		w->comment = (char *) malloc (strlen (line) + 1);
		if (w->comment == NULL) goto nomem;
		strcpy (w->comment, line);
		goto wdone;
	    }

/* Otherwise, parse keyword=value lines */

	    value = strchr (line, '=');
	    if (value == NULL || (int) (value - (char *) line) >= sizeof (keyword) - 1) {
		char	buf[2100];
illegal_line:	sprintf (buf, "Illegal line in worktodo.txt file: %s\n", line);
		OutputSomewhere (MAIN_THREAD_NUM, buf);
		goto comment;
	    }
	    *value = 0;
	    strcpy (keyword, line);
	    *value++ = '=';

/* Set some default values.  Historically, this program worked on */
/* Mersenne numbers only.  Default to an FFT length chosen by gwnum library. */

	    w->k = 1.0;
	    w->b = 2;
	    w->c = -1;
	    w->forced_fftlen = 0;
	    w->extension[0] = 0;

/* Parse the optional assignment_uid */

	    if ((value[0] == 'N' || value[0] == 'n') &&
		(value[1] == '/') &&
		(value[2] == 'A' || value[2] == 'a') &&
		(value[3] == ',')) {
		w->ra_failed = TRUE;
		safe_strcpy (value, value+4);
	    }
	    for (i = 0; ; i++) {
		if (!(value[i] >= '0' && value[i] <= '9') &&
		    !(value[i] >= 'A' && value[i] <= 'F') &&
		    !(value[i] >= 'a' && value[i] <= 'f')) break;
		if (i == 31) {
			if (value[32] != ',') break;
			value[32] = 0;
			strcpy (w->assignment_uid, value);
			safe_strcpy (value, value+33);
			break;
		}
	    }

/* Parse the FFT length to use.  The syntax is FFT_length for x87 cpus and */
/* FFT2_length for SSE2 machines.  We support two syntaxes so that an */
/* assignment moved from an x87 to-or-from an SSE2 machine will recalculate */
/* the soft FFT crossover. */

	    if ((value[0] == 'F' || value[0] == 'f') &&
	        (value[1] == 'F' || value[1] == 'f') &&
	        (value[2] == 'T' || value[2] == 't')) {
		int	sse2;
		unsigned long fftlen;
		char	*p;

		if (value[3] == '2') {
			sse2 = TRUE;
			p = value+5;
		} else {
			sse2 = FALSE;
			p = value+4;
		}
		fftlen = atoi (p);
		while (isdigit (*p)) p++;
		if (*p == 'K' || *p == 'k') fftlen <<= 10, p++;
		if (*p == 'M' || *p == 'm') fftlen <<= 20, p++;
		if (*p == ',') p++;
		safe_strcpy (value, p);
		if ((sse2 && (CPU_FLAGS & CPU_SSE2)) ||
		    (!sse2 && ! (CPU_FLAGS & CPU_SSE2)))
			w->forced_fftlen = fftlen;
	    }

/* Parse the optional file extension to use on save files (no good use */
/* right now, was formerly used for multiple workers ECMing the same number) */

	    if ((value[0] == 'E' || value[0] == 'e') &&
	        (value[1] == 'X' || value[1] == 'x') &&
	        (value[2] == 'T' || value[2] == 't') &&
		value[3] == '=') {
		char	*comma, *p;

		p = value+4;
		comma = strchr (p, ',');
		if (comma != NULL) {
			*comma = 0;
			if (strlen (p) > 8) p[8] = 0;
			strcpy (w->extension, p);
			safe_strcpy (value, comma+1);
		}
	    }

/* Handle Test= and DoubleCheck= lines.					*/
/*	Test=exponent,how_far_factored,has_been_pminus1ed		*/
/*	DoubleCheck=exponent,how_far_factored,has_been_pminus1ed	*/

	    if (_stricmp (keyword, "Test") == 0) {
		float	sieve_depth;
		w->work_type = WORK_TEST;
		sieve_depth = 0.0;
		sscanf (value, "%lu,%f,%d",
				&w->n, &sieve_depth, &w->pminus1ed);
		w->sieve_depth = sieve_depth;
		w->tests_saved = 2.0;
	    }
	    else if (_stricmp (keyword, "DoubleCheck") == 0) {
		float	sieve_depth;
		w->work_type = WORK_DBLCHK;
		sieve_depth = 0.0;
		sscanf (value, "%lu,%f,%d",
				&w->n, &sieve_depth, &w->pminus1ed);
		w->sieve_depth = sieve_depth;
		w->tests_saved = 1.0;
	    }

/* Handle AdvancedTest= lines. */
/*	AdvancedTest=exponent */

	    else if (_stricmp (keyword, "AdvancedTest") == 0) {
		w->work_type = WORK_ADVANCEDTEST;
		sscanf (value, "%lu", &w->n);
	    }

/* Handle Factor= lines.  Old style is:					*/
/*	Factor=exponent,how_far_factored				*/
/* New style is:							*/
/*	Factor=exponent,how_far_factored,how_far_to_factor_to		*/

	    else if (_stricmp (keyword, "Factor") == 0) {
		float	sieve_depth, factor_to;
		w->work_type = WORK_FACTOR;
		sieve_depth = 0.0;
		factor_to = 0.0;
		sscanf (value, "%lu,%f,%f",
				&w->n, &sieve_depth, &factor_to);
		w->sieve_depth = sieve_depth;
		w->factor_to = factor_to;
	    }

/* Handle Pfactor= lines.  Old style is:				*/
/*	Pfactor=exponent,how_far_factored,double_check_flag		*/
/* New style is:							*/
/*	Pfactor=k,b,n,c,how_far_factored,ll_tests_saved_if_factor_found	*/

	    else if (_stricmp (keyword, "PFactor") == 0) {
		float	sieve_depth;
		w->work_type = WORK_PFACTOR;
		sieve_depth = 0.0;
		if (countCommas (value) > 3) {		/* New style */
			char	*q;
			float	tests_saved;
			tests_saved = 0.0;
			q = strchr (value, ','); *q = 0; w->k = atof (value);
			sscanf (q+1, "%lu,%lu,%ld,%f,%f",
				&w->b, &w->n, &w->c, &sieve_depth,
				&tests_saved);
			w->sieve_depth = sieve_depth;
			w->tests_saved = tests_saved;
		} else {				/* Old style */
			int	dblchk;
			sscanf (value, "%lu,%f,%d",
				&w->n, &sieve_depth, &dblchk);
			w->sieve_depth = sieve_depth;
			w->tests_saved = dblchk ? 1.0 : 2.0;
		}
	    }

/* Handle ECM= lines.  Old style is: */
/*   ECM=exponent,B1,B2,curves_to_do,unused[,specific_sigma,plus1,B2_start] */
/* New style is: */
/*   ECM2=k,b,n,c,B1,B2,curves_to_do[,specific_sigma,B2_start][,"factors"] */

	    else if (_stricmp (keyword, "ECM") == 0) {
		char	*q;
		w->work_type = WORK_ECM;
		sscanf (value, "%ld", &w->n);
		if ((q = strchr (value, ',')) == NULL) goto illegal_line;
		w->B1 = atof (q+1);
		if ((q = strchr (q+1, ',')) == NULL) goto illegal_line;
		w->B2 = atof (q+1);
		if ((q = strchr (q+1, ',')) == NULL) goto illegal_line;
		w->curves_to_do = atoi (q+1);
		if ((q = strchr (q+1, ',')) == NULL) goto illegal_line;
		q = strchr (q+1, ',');
		w->curve = 0;
		if (q != NULL) {
			w->curve = atof (q+1);
			q = strchr (q+1, ',');
		}
		if (q != NULL) {
			w->c = atoi (q+1);
			if (w->c == 0) w->c = -1; /* old plus1 arg */
			q = strchr (q+1, ',');
		}
		w->B2_start = w->B1;
		if (q != NULL) {
			double j;
			j = atof (q+1);
			if (j > w->B1) w->B2_start = j;
		}
	    } else if (_stricmp (keyword, "ECM2") == 0) {
		int	i;
		char	*q;
		w->work_type = WORK_ECM;
		w->k = atof (value);
		if ((q = strchr (value, ',')) == NULL) goto illegal_line;
		sscanf (q+1, "%lu,%lu,%ld", &w->b, &w->n, &w->c);
		for (i = 1; i <= 3; i++)
			if ((q = strchr (q+1, ',')) == NULL) goto illegal_line;
		w->B1 = atof (q+1);
		if ((q = strchr (q+1, ',')) == NULL) goto illegal_line;
		w->B2 = atof (q+1);
		if ((q = strchr (q+1, ',')) == NULL) goto illegal_line;
		w->curves_to_do = atoi (q+1);
		q = strchr (q+1, ',');
		w->curve = 0;
		if (q != NULL && q[1] != '"') {
			w->curve = atof (q+1);
			q = strchr (q+1, ',');
		}
		w->B2_start = w->B1;
		if (q != NULL && q[1] != '"') {
			double j;
			j = atof (q+1);
			if (j > w->B1) w->B2_start = j;
			q = strchr (q+1, ',');
		}
		if (q != NULL && q[1] == '"') {
			w->known_factors = (char *) malloc (strlen (q));
			if (w->known_factors == NULL) goto nomem;
			strcpy (w->known_factors, q+2);
		}
	    }

/* Handle Pminus1 lines:  Old style:				*/
/*	Pminus1=exponent,B1,B2,plus1[,B2_start]			*/
/* New style is:						*/
/*	Pminus1=k,b,n,c,B1,B2[,how_far_factored][,B2_start][,"factors"] */

	    else if (_stricmp (keyword, "Pminus1") == 0) {
		char	*q;
		w->work_type = WORK_PMINUS1;
		if (countCommas (value) <= 4) {
			sscanf (value, "%ld", &w->n);
			if ((q = strchr (value, ',')) == NULL)
				goto illegal_line;
			w->B1 = atof (q+1);
			if ((q = strchr (q+1, ',')) == NULL) goto illegal_line;
			w->B2 = atof (q+1);
			if ((q = strchr (q+1, ',')) == NULL) goto illegal_line;
			sscanf (q+1, "%ld", &w->c);
			q = strchr (q+1, ',');
			if (w->c == 0) w->c = -1; /* old plus1 arg */
			if (q != NULL) {
				double j;
				j = atof (q+1);
				if (j > w->B1) w->B2_start = j;
			}
		} else {
			w->k = atof (value);
			if ((q = strchr (value, ',')) == NULL)
				goto illegal_line;
			sscanf (q+1, "%lu,%lu,%ld", &w->b, &w->n, &w->c);
			for (i = 1; i <= 3; i++)
				if ((q = strchr (q+1, ',')) == NULL)
					goto illegal_line;
			w->B1 = atof (q+1);
			if ((q = strchr (q+1, ',')) == NULL) goto illegal_line;
			w->B2 = atof (q+1);
			q = strchr (q+1, ',');
			w->sieve_depth = 0.0;
			if (q != NULL && q[1] != '"') {
				double	j;
				j = atof (q+1);
				if (j < 100.0) {
					w->sieve_depth = j;
					q = strchr (q+1, ',');
				}
			}
			w->B2_start = 0;
			if (q != NULL && q[1] != '"') {
				double	j;
				j = atof (q+1);
				if (j > w->B1) w->B2_start = j;
				q = strchr (q+1, ',');
			}
			if (q != NULL && q[1] == '"') {
				w->known_factors = (char *) malloc (strlen (q));
				if (w->known_factors == NULL) goto nomem;
				strcpy (w->known_factors, q+2);
			}
		}
	    }

/* Handle PRP= lines.							*/
/*	PRP=k,b,n,c[,how_far_factored,tests_saved][,known_factors]	*/
/* A tests_saved value of 0.0 will bypass any P-1 factoring		*/

	    else if (_stricmp (keyword, "PRP") == 0) {
		char	*q;

		w->work_type = WORK_PRP;
		w->k = atof (value);
		if ((q = strchr (value, ',')) == NULL) goto illegal_line;
		sscanf (q+1, "%lu,%lu,%ld", &w->b, &w->n, &w->c);
		for (i = 1; i <= 2; i++)
			if ((q = strchr (q+1, ',')) == NULL) goto illegal_line;
		q = strchr (q+1, ',');

		if (q != NULL && q[1] != '"') {
			w->sieve_depth = atof (q+1);
			if ((q = strchr (q+1, ',')) == NULL) goto illegal_line;
			w->tests_saved = atof (q+1);
			q = strchr (q+1, ',');
		} else {
			w->sieve_depth = 0.0;
			w->tests_saved = 0.0;
		}
		if (q != NULL && q[1] == '"') {
			w->known_factors = (char *) malloc (strlen (q));
			if (w->known_factors == NULL) goto nomem;
			strcpy (w->known_factors, q+2);
		}
	    }

/* Uh oh.  We have a worktodo.ini line we cannot process. */

	    else if (_stricmp (keyword, "AdvancedFactor") == 0) {
		OutputSomewhere (MAIN_THREAD_NUM, "Worktodo error: AdvancedFactor no longer supported\n");
		goto comment;
	    } else {
		goto illegal_line;
	    }

/* Trim trailing non-digit characters from known factors list (this should be the closing double quote) */
/* Turn all non-digit characters into commas (they should be anyway) */

	    if (w->known_factors != NULL) {
		for (i = (unsigned int) strlen (w->known_factors);
		     i > 0 && !isdigit (w->known_factors[i-1]);
		     i--);
		w->known_factors[i] = 0;
		for (i = 0; i < (unsigned int) strlen (w->known_factors); i++)
			if (!isdigit (w->known_factors[i])) w->known_factors[i] = ',';
	    }

/* Make sure this line of work from the file makes sense. The exponent */
/* should be a prime number, bounded by values we can handle, and we */
/* should never be asked to factor a number more than we are capable of. */

	    if (w->k == 1.0 && w->b == 2 && !isPrime (w->n) && w->c == -1 &&
	        w->work_type != WORK_ECM && w->work_type != WORK_PMINUS1) {
		char	buf[80];
		sprintf (buf, "Error: Worktodo.txt file contained composite exponent: %ld\n", w->n);
		OutputBoth (MAIN_THREAD_NUM, buf);
		goto illegal_line;
	    }
	    if ((w->work_type == WORK_TEST ||
	         w->work_type == WORK_DBLCHK ||
	         w->work_type == WORK_ADVANCEDTEST) &&
	        (w->n < MIN_PRIME ||
		 (w->forced_fftlen == 0 &&
		  w->n > (unsigned long) (CPU_FLAGS & (CPU_AVX | CPU_SSE2) ? MAX_PRIME_SSE2 : MAX_PRIME)))) {
		char	buf[80];
		sprintf (buf, "Error: Worktodo.txt file contained bad LL exponent: %ld\n", w->n);
		OutputBoth (MAIN_THREAD_NUM, buf);
		goto illegal_line;
	    }
	    if (w->work_type == WORK_FACTOR && w->n < 20000) {
		char	buf[100];
		sprintf (buf, "Error: Use ECM instead of trial factoring for exponent: %ld\n", w->n);
		OutputBoth (MAIN_THREAD_NUM, buf);
		goto illegal_line;
	    }
	    if (w->work_type == WORK_FACTOR && w->n > MAX_FACTOR) {
		char	buf[100];
		sprintf (buf, "Error: Worktodo.txt file contained bad factoring assignment: %ld\n", w->n);
		OutputBoth (MAIN_THREAD_NUM, buf);
		goto illegal_line;
	    }

/* A user discovered a case where a computer that dual boots between 32-bit prime95 */
/* and 64-bit prime95 can run into problems.  If near the FFT limit an FFT length is */
/* picked and written to worktodo.ini.  When running the other executable, that FFT */
/* length may not be supported leading to a "cannot initialize FFT error".  For */
/* example, the 2800K FFT length is implemented in 64-bit prime95, but not 32-bit prime95. */
/* The quick workaround here is to ignore FFT lengths from the worktodo file if that FFT */
/* length is not supported.  This is non-optimal because the proper FFT size will */
/* have to be recalculated. */

	    if (w->forced_fftlen && gwmap_fftlen_to_max_exponent (w->forced_fftlen) == 0) {
		    char	buf[100];
		    sprintf (buf, "Warning: Ignoring unsupported FFT length, %ld, on line %u of worktodo.txt.\n",
			     w->forced_fftlen, linenum);
		    OutputBoth (MAIN_THREAD_NUM, buf);
		    w->forced_fftlen = 0;
	    }

/* Do more initialization of the work_unit structure */

	    auxiliaryWorkUnitInit (w);

/* Grow the work_unit array if necessary and add this entry */

wdone:	    rc = addToWorkUnitArray (tnum, w, TRUE);
	    if (rc) goto retrc;
	}

/* Now that we've finished reading the worktodo file, set stage */
/* and pct_complete based on existing save files. */

	for (tnum = 0; tnum < MAX_NUM_WORKER_THREADS; tnum++) {
	    struct work_unit *w;
	    int first_real_work_line;

	    first_real_work_line = TRUE;
	    for (w = WORK_UNITS[tnum].first; w != NULL; w = w->next) {

/* Init assuming we won't find a save file. */

		w->stage[0] = 0;
		w->pct_complete = 0.0;

/* Skip comment lines */

		if (w->work_type == WORK_NONE) continue;

/* If well behaved work is set, only first work lines can have a save file */

		if (WELL_BEHAVED_WORK && !first_real_work_line) goto next_wu;

/* Set stage and pct_complete for work units that have already begun */
/* based on data in the save files.  Only do this for the first appearance */
/* of a number for a worker.  For example, if a worker has several entries */
/* ECMing the same number, only the first entry will have the pct_complete set. */
/* We also assume an existing save file is used for the first worker of a thread */
/* rather than a non-first work unit in an earlier thread. */

		if (!first_real_work_line) {
			int	tnum2;
			struct work_unit *w2;

/* See if any other worker's first work unit is testing the same number. */
/* If so, assume any existing save files are for that worker */

			for (tnum2 = 0; tnum2 < MAX_NUM_WORKER_THREADS; tnum2++) {
				for (w2 = WORK_UNITS[tnum2].first; w2 != NULL; w2 = w2->next) {
					if (w2->work_type == WORK_NONE) continue;
					if (w2->work_type == w->work_type &&
					    w2->k == w->k &&
					    w2->b == w->b &&
					    w2->n == w->n &&
					    w2->c == w->c) goto next_wu;
					break;
				}
			}

/* See if any earlier work units in this worker are testing the same number. */
/* If so, assume any existing save files are for that work unit. */

			for (w2 = WORK_UNITS[tnum].first; w2 != w; w2 = w2->next) {
				if (w2->work_type == w->work_type &&
				    w2->k == w->k &&
				    w2->b == w->b &&
				    w2->n == w->n &&
				    w2->c == w->c) goto next_wu;
			}
		}

/* Now see if an existing save file can be used to set stage and pct_complete */

		pct_complete_from_savefile (w);

/* Progress to the next work unit */

next_wu:	first_real_work_line = FALSE;
	    }
	}

/* Close the file, free the lock and return success */

	fclose (fd);
done:	gwmutex_unlock (&WORKTODO_MUTEX);

/* Almost done.  Incorporate the optional worktodo.add file. */

	return (incorporateWorkToDoAddFile ());
	
/* Close the file, free the lock and return error from routine we called */

retrc:	fclose (fd);
	gwmutex_unlock (&WORKTODO_MUTEX);
	return (rc);

/* Free the lock and return out of memory error code */

nomem:	fclose (fd);
	gwmutex_unlock (&WORKTODO_MUTEX);
	return (OutOfMemory (MAIN_THREAD_NUM));
}

/* Write the updated worktodo.ini to disk */

int writeWorkToDoFile (
	int	force)		/* Force writing file even if WELL_BEHAVED */
{
	char	newFileName[80];
	int	fd, mangled, last_line_was_blank;
	unsigned int tnum;

/* If work to do hasn't changed, then don't write the file */

	if (!WORKTODO_CHANGED) return (0);

/* If the well-behaved-work-option is on, then only write the file every */
/* half hour.  The user should set this option when the worktodo file is */
/* long and the work units complete quickly.  This commonly happens when */
/* trial factoring a large number of exponents to a low limit. */

	if (WELL_BEHAVED_WORK && !force) {
		static time_t last_time_written = 0;
		time_t	current_time;
		time (&current_time);
		if (current_time < last_time_written + 1800) return (0);
		last_time_written = current_time;
	}

/* Grab the lock so that comm thread cannot try to add work units while */
/* file is being written. */

	gwmutex_lock (&WORKTODO_MUTEX);

/* Create the WORKTODO.INI file */

	mangleIniFileName (WORKTODO_FILE, newFileName, &mangled);
	fd = _open (newFileName, _O_CREAT | _O_TRUNC | _O_WRONLY | _O_TEXT, CREATE_FILE_ACCESS);
	if (fd < 0) {
		OutputBoth (MAIN_THREAD_NUM, "Error creating worktodo.txt file\n");
		return (STOP_FILE_IO_ERROR);
	}

/* Loop over all worker threads */

	last_line_was_blank = FALSE;
	for (tnum = 0; tnum < MAX_NUM_WORKER_THREADS; tnum++) {
	    struct work_unit *w;
	    unsigned int len;

/* If we've processed all the worker threads and there is nothing left */
/* to output, then we are done. */

	    if (tnum >= NUM_WORKER_THREADS &&
		WORK_UNITS[tnum].first == NULL) break;

/* Output a standardized section header */

	    if (tnum || NUM_WORKER_THREADS > 0) {
		char	buf[40];
		if (tnum && !last_line_was_blank && _write (fd, "\n", 1) != 1)
			goto write_error;
		sprintf (buf, "[Worker #%d]\n", tnum+1);
		len = (unsigned int) strlen (buf);
		if (_write (fd, buf, len) != len) goto write_error;
		last_line_was_blank = FALSE;
	    }

/* Loop over each assignment for this worker thread */

	    for (w = WORK_UNITS[tnum].first; w != NULL; w = w->next) {
		char	idbuf[100];
		char	buf[4096];

/* Do not output deleted lines */

		if (w->work_type == WORK_DELETED) continue;
		
/* Do not output section header */

		if (w == WORK_UNITS[tnum].first &&
		    w->work_type == WORK_NONE &&
		    w->comment[0] == '[') continue;
		
/* Format the optional assignment id */

		idbuf[0] = 0;
		if (w->assignment_uid[0])
			sprintf (idbuf, "%s,", w->assignment_uid);
		else if (w->ra_failed)
			sprintf (idbuf, "%s,", "N/A");

/* Format the FFT length */

		if (w->forced_fftlen) {
			strcat (idbuf, "FFT");
			if (CPU_FLAGS & CPU_SSE2) strcat (idbuf, "2");
			if ((w->forced_fftlen & 0xFFFFF) == 0)
				sprintf (idbuf+strlen(idbuf), "=%luM,", w->forced_fftlen >> 20);
			else if ((w->forced_fftlen & 0x3FF) == 0)
				sprintf (idbuf+strlen(idbuf), "=%luK,", w->forced_fftlen >> 10);
			else
				sprintf (idbuf+strlen(idbuf), "=%lu,", w->forced_fftlen);
		}

/* Output the optional file name extension (no good use right now, */
/* was formerly used for multiple workers ECMing the same number) */

		if (w->extension[0]) {
			sprintf (idbuf+strlen(idbuf), "EXT=%s,", w->extension);
		}

/* Write out comment lines just as we read them in */
/* Format normal work unit lines */

		switch (w->work_type) {

		case WORK_NONE:
			strcpy (buf, w->comment);
			break;

		case WORK_TEST:
			sprintf (buf, "Test=%s%lu,%.0f,%d",
				 idbuf, w->n, w->sieve_depth, w->pminus1ed);
			break;

		case WORK_DBLCHK:
			sprintf (buf, "DoubleCheck=%s%lu,%.0f,%d",
				 idbuf, w->n, w->sieve_depth, w->pminus1ed);
			break;

		case WORK_ADVANCEDTEST:
			sprintf (buf, "AdvancedTest=%lu", w->n);
			break;

		case WORK_FACTOR:
			sprintf (buf, "Factor=%s%ld,%.0f,%.0f",
				 idbuf, w->n, w->sieve_depth, w->factor_to);
			break;

		case WORK_PFACTOR:
			sprintf (buf, "Pfactor=%s%.0f,%lu,%lu,%ld,%g,%g",
				 idbuf, w->k, w->b, w->n, w->c, w->sieve_depth,
				 w->tests_saved);
			break;

		case WORK_ECM:
			sprintf (buf,
				 "ECM2=%s%.0f,%lu,%lu,%ld,%.0f,%.0f,%u",
				 idbuf, w->k, w->b, w->n, w->c, w->B1, w->B2,
				 w->curves_to_do);
			if (w->B2_start > w->B1)
				sprintf (buf + strlen (buf),
					 ",%.0f,%.0f",
					 w->curve, w->B2_start);
			else if (w->curve)
				sprintf (buf + strlen (buf),
					 ",%.0f",
					 w->curve);
			if (w->known_factors != NULL)
				sprintf (buf + strlen (buf),
					 ",\"%s\"",
					 w->known_factors);
			break;

		case WORK_PMINUS1:
			sprintf (buf, "Pminus1=%s%.0f,%lu,%lu,%ld,%.0f,%.0f",
				 idbuf, w->k, w->b, w->n, w->c, w->B1, w->B2);
			if (w->sieve_depth > 0.0)
				sprintf (buf + strlen (buf), ",%.0f", w->sieve_depth);
			if (w->B2_start > w->B1)
				sprintf (buf + strlen (buf), ",%.0f", w->B2_start);
			if (w->known_factors != NULL)
				sprintf (buf + strlen (buf), ",\"%s\"", w->known_factors);
			break;

		case WORK_PRP:
			sprintf (buf, "PRP=%s%.0f,%lu,%lu,%ld", idbuf, w->k, w->b, w->n, w->c);
			if (w->tests_saved > 0.0)
				sprintf (buf + strlen (buf), ",%g,%g", w->sieve_depth, w->tests_saved);
			if (w->known_factors != NULL)
				sprintf (buf + strlen (buf), ",\"%s\"", w->known_factors);
			break;
		}

/* Write out the formatted line */

		strcat (buf, "\n");
		len = (unsigned int) strlen (buf);
		if (_write (fd, buf, len) != len) {
write_error:		OutputBoth (MAIN_THREAD_NUM,
				    "Error writing worktodo.txt file\n");
			_close (fd);
			return (STOP_FILE_IO_ERROR);
		}
		last_line_was_blank = (len == 1);
	    }
	}

/* Close file, unlock, and return success */

	_close (fd);
	if (mangled) _unlink (WORKTODO_FILE);
	WORKTODO_CHANGED = FALSE;
	gwmutex_unlock (&WORKTODO_MUTEX);
	return (0);
}

/* Return a worktodo.ini entry for the given worker thread */

struct work_unit *getNextWorkToDoLine (
	int	thread_num,		/* Thread number starting from 0 */
	struct work_unit *w,		/* Current WorkToDo entry (or NULL) */
	int	usage)			/* Short vs. long term usage */
{
	struct work_unit *next;		/* Next WorkToDo entry (or NULL) */

	ASSERTG (thread_num < (int) NUM_WORKER_THREADS);

/* Grab the lock so that other threads do not add or delete lines */
/* while we are finding the next entry. */

	gwmutex_lock (&WORKTODO_MUTEX);

/* Get pointer to the first or next work unit */

	next = (w == NULL) ? WORK_UNITS[thread_num].first : w->next;
	while (next != NULL && next->work_type == WORK_DELETED)
		next = next->next;

/* Decrement the use count */

	if (w != NULL) {
		w->in_use_count--;
		WORKTODO_IN_USE_COUNT--;
		if (usage == LONG_TERM_USE) w->in_use_count &= ~0x80000000;

/* Free the work unit if it has been deleted and use count is now zero */

		if (w->work_type == WORK_DELETED && w->in_use_count == 0) {

/* Unlink the work unit from the list */

			if (w->prev == NULL)
				WORK_UNITS[thread_num].first = w->next;
			else
				w->prev->next = w->next;
			if (w->next == NULL)
				WORK_UNITS[thread_num].last = w->prev;
			else
				w->next->prev = w->prev;

/* Free memory allocated for this work unit */

			free (w->known_factors);
			free (w->comment);
			free (w);
		}
	}

/* Increment the in-use count.  If this is a long-term usage, remember */
/* the fact so that deleting the work unit can be prohibited. */

	if (next != NULL) {
		next->in_use_count++;
		WORKTODO_IN_USE_COUNT++;
		if (usage == LONG_TERM_USE) next->in_use_count |= 0x80000000;
	}

/* Unlock and return */

	gwmutex_unlock (&WORKTODO_MUTEX);
	return (next);
}

/* Return a worktodo.ini entry for the given worker thread */

void decrementWorkUnitUseCount (
	struct work_unit *w,		/* WorkToDo entry */
	int	usage)			/* Short vs. long term usage */
{
	ASSERTG (w->in_use_count != 0);

/* Grab the lock */

	gwmutex_lock (&WORKTODO_MUTEX);

/* Decrement the reader count before getting the next entry */

	w->in_use_count--;
	WORKTODO_IN_USE_COUNT--;
	if (usage == LONG_TERM_USE) w->in_use_count &= ~0x80000000;

/* Unlock and return */

	gwmutex_unlock (&WORKTODO_MUTEX);
}

/* Add a line of work to the work-to-do INI file. */

int addWorkToDoLine (
	int	tnum,		/* Worker thread to process this work unit */
	struct work_unit *w)	/* Type of work */
{
	struct work_unit *malloc_w;
	int	rc;

/* Do more initialization of the work_unit structure */

	auxiliaryWorkUnitInit (w);

/* Copy work unit from stack to malloc'ed area */

	malloc_w = (struct work_unit *) malloc (sizeof (struct work_unit));
	if (malloc_w == NULL) return (OutOfMemory (MAIN_THREAD_NUM));
	memcpy (malloc_w, w, sizeof (struct work_unit));

/* Grab the lock so that comm thread and/or worker threads do not */
/* access structure while the other is adding/deleting lines. */

	gwmutex_lock (&WORKTODO_MUTEX);

/* Add the work unit to the end of the array.  Well, actually before the */
/* last set of blank lines. */

	rc = addToWorkUnitArray (tnum, malloc_w, FALSE);
	if (rc) goto retrc;

/* Set flag indicating worktodo needs writing */

	WORKTODO_CHANGED = TRUE;

/* Unlock and write the worktodo.ini file to disk */

	gwmutex_unlock (&WORKTODO_MUTEX);
	return (writeWorkToDoFile (FALSE));

/* Unlock and return error code from called routine */

retrc:	gwmutex_unlock (&WORKTODO_MUTEX);
	return (rc);
}

/* Caller has updated a work unit structure such that the work-to-do */
/* INI file needs to be written.  */

int updateWorkToDoLine (
	int	tnum,		/* Worker thread to process this work unit */
	struct work_unit *w)	/* Type of work */
{

/* Set flag indicating worktodo needs writing */

	WORKTODO_CHANGED = TRUE;

/* Write the worktodo.ini file to disk */

	return (writeWorkToDoFile (FALSE));
}

/* Delete a line of work from the work-to-do INI file.  Even if an error */
/* occurs the work unit has been deleted and the use_count decremented. */
/* The work_unit pointer (w) will be set to the previous work unit so that */
/* getNextWorkToDoLine works properly. */

int deleteWorkToDoLine (
	int	tnum,		/* Worker thread to process this work unit */
	struct work_unit *w,	/* Work unit pointer */
	int	stop_if_in_progress) /* Stop thread processing work unit */
{

/* If this work unit has already been deleted, then ignore this delete */
/* request.  This should only happen in a bizarre race condition. */

	if (w->work_type == WORK_DELETED) return (0);

/* Grab the lock so that comm thread and/or worker threads do not */
/* access structure while the other is adding/deleting lines. */

	gwmutex_lock (&WORKTODO_MUTEX);

/* If we are deleting a work unit that is currently being worked on */
/* then set a flag so that the worker thread will abort the work unit */
/* at its first opportunity.  There is a race condition whereby the worker */
/* thread may move on to a different work unit before testing the */
/* abort-work-unit flag.  This is harmless as the new work unit will be */
/* aborted and then restarted. */

	if ((w->in_use_count & 0x80000000) && stop_if_in_progress)
		stop_worker_for_abort (tnum);

/* Decrement count of valid work lines */

	if (w->work_type != WORK_NONE) WORKTODO_COUNT--;

/* Mark this work unit deleted */

	w->work_type = WORK_DELETED;

/* Set flag indicating worktodo needs writing */

	WORKTODO_CHANGED = TRUE;

/* Unlock and write the worktodo.ini file to disk */

	gwmutex_unlock (&WORKTODO_MUTEX);
	return (writeWorkToDoFile (FALSE));
}

/* Return TRUE if a work unit is currently being worked on. */

int isWorkUnitActive (
	struct work_unit *w)	/* work_unit pointer */
{

/* At the present time only a worker thread actively working on a work unit */
/* sets the LONG_TERM_USE flag of in_use_count.  Use that fact, to decide */
/* if the work unit is active. */

	if (w->in_use_count & 0x80000000) return (TRUE);
	return (FALSE);
}

/* Make a guess as to how much longer a chore should take. */

double work_estimate (
	int	thread_num,		/* Thread number doing the work */
	struct work_unit *w)
{
	double	timing, est, pct_complete;
	int	can_use_multiple_threads;
	unsigned int i, total_threads;

/* I suppose there are race conditions where a deleted work unit could */
/* get here.  Return an estimate of 0.0. */

	est = 0.0;

/* Make sure the pct_complete is between 0.0 and 1.0.  There is presently */
/* a bug in P-1 code that is setting this value to more than 1.0 */

	pct_complete = w->pct_complete;
	if (pct_complete < 0.0) pct_complete = 0.0;
	if (pct_complete > 1.0) pct_complete = 1.0;

/* Only large SSE2 FFTs can use multiple threads. */

	can_use_multiple_threads = (CPU_FLAGS & CPU_SSE2 && w->n > 172700);

/* For ECM, estimating time is very difficult because it depends on how */
/* much memory is available for temporaries in stage 2.  There are also */
/* significant increases in time when numbers no longer fit in the L2 */
/* or L3 cache.  */
/* For small numbers, we do about 12.86 * B1 squarings in stage 1 and */
/* 0.0617 * (B2 - B1) squarings in stage 2.  The extra overhead due to */
/* adds/subs and gwfftmul, gwfftfftmul being less efficient than gwsquare */
/* adds overhead, I'm guessing 10% for tiny numbers that likely fit in the */
/* cache to 20% for numbers that don't fit in the cache. */
/* For larger numbers, fewer temporaries and greater costs for modinv cause */
/* us to eventually switch from a 2 FFTs/prime to 4 FFTs/prime strategy in */
/* stage 2.  This switch occurs around 5,000,000 bits on a machine using */
/* modest amounts of memory.  We'll be doing 0.1261 * (B2 - B1) stage 2 */
/* squarings then.  Between 100,000 bits and 5,000,000 bits we'll gradually */
/* increase the stage 2 cost to account for the fewer temporaries resulting */
/* in more modular inverses combined with modular inverses getting more */
/* and more expensive. */
/* Also note: the stage text is C<curve#>S<stage#>. */

	if (w->work_type == WORK_ECM) {
		int	full_curves_to_do, stage, bits;
		double	stage1_time, stage2_time, overhead, B2_minus_B1;

		full_curves_to_do = w->curves_to_do;
		stage = 0;
		if (w->stage[0]) {
			full_curves_to_do -= atoi (&w->stage[1]);
			stage = atoi (&w->stage[strlen(w->stage)-1]);
		}

		timing = gwmap_to_timing (w->k, w->b, w->n, w->c);
		bits = (int) (w->n * log ((double) w->b) / log (2.0));
		if (bits <= 80000) overhead = 1.10;
		else if (bits >= 1500000) overhead = 1.20;
		else overhead = 1.10 + ((double) bits - 80000.0) / 1420000.0 * (1.20 - 1.10);

		stage1_time = (12.86 * w->B1) * timing * overhead;

		if (bits <= 100000) stage2_time = 0.0617;
		else if (bits >= 5000000) stage2_time = 0.1261;
		else stage2_time = 0.0617 + ((double) bits - 100000.0) / 4900000.0 * (0.1261 - 0.0617);
		B2_minus_B1 = (w->B2 > 0.0) ? w->B2 - w->B1 : 99.0 * w->B1;
		stage2_time = stage2_time * B2_minus_B1 * timing * overhead;

		est = (double) full_curves_to_do * (stage1_time + stage2_time);
		if (stage == 1)
			est += stage1_time * (1.0 - pct_complete) + stage2_time;
		if (stage == 2)
			est += stage2_time * (1.0 - pct_complete);
	}

/* For P-1, estimate about 1.4545 * B1 squarings in stage 1 and 0.06154 * B2 */
/* squarings in stage 2.  Note that the stage 2 estimate is quite */
/* optimistic for large numbers as fewer temporaries will result in nearly */
/* double the number of squarings.  Also, pass 2 squarings are 28.5% slower */
/* (due to all the adds). */ 

	if (w->work_type == WORK_PMINUS1 || w->work_type == WORK_PFACTOR) {
		int	stage;
		double	B1, B2;
		double	stage1_time, stage2_time;

		if (w->work_type == WORK_PFACTOR) {
			unsigned long guess_B1, guess_B2;
			unsigned long squarings;
			double	prob;
			guess_pminus1_bounds (thread_num, w->k, w->b, w->n, w->c,
					      w->sieve_depth, w->tests_saved,
					      &guess_B1, &guess_B2,
					      &squarings, &prob);
			B1 = guess_B1;
			B2 = guess_B2;
		} else {
			B1 = w->B1;
			B2 = w->B2;
		}

		if (w->stage[0]) stage = atoi (&w->stage[1]);
		else stage = 0;

		timing = gwmap_to_timing (w->k, w->b, w->n, w->c);
		stage1_time = timing * (1.4545 * B1);
		if (B2)
			stage2_time = timing * (0.06154 * (B2 - B1)) * 1.285;
		else
			stage2_time = timing * (0.06154 * 99.0 * B1) * 1.285;

		if (stage == 0)
			est = stage1_time + stage2_time;
		if (stage == 1)
			est = stage1_time * (1.0 - pct_complete) + stage2_time;
		if (stage == 2)
			est = stage2_time * (1.0 - pct_complete);
	}

/* If factoring, guess how long that will take.  Timings are based on */
/* the factoring benchmark for my 2 GHz P4. */
/*	Best time for 60 bit trial factors: 15.123 ms. */
/*	Best time for 61 bit trial factors: 15.021 ms. */
/*	Best time for 62 bit trial factors: 15.080 ms. */
/*	Best time for 63 bit trial factors: 16.127 ms. */
/*	Best time for 64 bit trial factors: 16.143 ms. */
/*	Best time for 65 bit trial factors: 20.230 ms. */
/*	Best time for 66 bit trial factors: 20.212 ms. */
/*	Best time for 67 bit trial factors: 20.244 ms. */
/* Factoring M35000011 from 2^60 to 2^61 takes 513 seconds.  Solve for */
/* constant C in this formula:  15.1 ms * 2^61 * C = 513 seconds */
/*	C = 513 sec / 2^61 / 15.1 ms */
/*	C = 513000 ms / 2^61 / 15.1 ms = 33974 / 2^61 */
/* Our estimate for factoring Mp to 2^i is then: */
/*	time_in_ms = benchmark_in_ms * 2^i * C * (35,000,011 / p) */
/* Which simplifies to: */
/*	time_in_seconds = benchmark_in_ms * 2^(i-48) * 145000 / p */

	if (w->work_type == WORK_FACTOR) {
		int	i, tf_level;

		can_use_multiple_threads = FALSE;

		if (w->stage[0]) tf_level = atoi (&w->stage[2]);
		else tf_level = 0;

		est = 0.0;
		for (i = (int) w->sieve_depth+1; i <= (int) w->factor_to; i++) {
			if (i < 48) continue;
			timing = (i > 64) ? 20.2 : (i > 62) ? 16.1 : 15.1;
			if (i <= 72) timing *= 145000.0 * (1L << (i - 48)) / w->n;
			else timing *= 145000.0 * 16777216.0 * (1L << (i - 72)) / w->n;
			if (i == tf_level)
				est += timing * (1.0 - pct_complete);
			else
				est += timing;
		}
		est = est * 2000.0 / CPU_SPEED;
	}

/* If testing add in the Lucas-Lehmer testing time */

	if (w->work_type == WORK_TEST ||
	    w->work_type == WORK_ADVANCEDTEST ||
	    w->work_type == WORK_DBLCHK) {
		est = w->n * gwmap_to_timing (w->k, w->b, w->n, w->c);
		if (w->stage[0] == 'L') est *= (1.0 - pct_complete);
	}

/* If PRPing add in the PRP testing time */

	if (w->work_type == WORK_PRP) {
		est = w->n * log ((double) w->b) / log (2.0) * gwmap_to_timing (w->k, w->b, w->n, w->c);
		if (w->stage[0] == 'P') est *= (1.0 - pct_complete);
	}

/* Factor in the hours per day the computer is running and the */
/* rolling average */

	est *= (24.0 / CPU_HOURS) * (1000.0 / ROLLING_AVERAGE);

/* If the worker uses multiple CPUs to execute the FFT, then adjust the */
/* time estimate.  As a rough estimate, assume the first additional CPU */
/* reduces the time by a factor of 1.7.  Also assume each additional CPU */
/* is less beneficial. */

	if (can_use_multiple_threads &&
	    THREADS_PER_TEST[thread_num] > CPU_HYPERTHREADS) {
		double	effective_num_cpus, cpu_value;
		effective_num_cpus = 0.0;
		cpu_value = 1.0;
		for (i = 0; i < THREADS_PER_TEST[thread_num]; i += CPU_HYPERTHREADS) {
			effective_num_cpus += cpu_value;
			cpu_value *= 0.7;
		}
		est = est / effective_num_cpus;
	}

/* If the user is unwisely running more threads than there are logical */
/* CPUs, then increase the time estimate */

	total_threads = 0;
	for (i = 0; i < NUM_WORKER_THREADS; i++)
		total_threads += THREADS_PER_TEST[i];
	if (total_threads > NUM_CPUS * CPU_HYPERTHREADS)
		est *= (double) total_threads /
		       (double) (NUM_CPUS * CPU_HYPERTHREADS);

/* Return the total estimated time in seconds */

	return (est);
}


/* Determine how much we should factor (in bits) */

unsigned int factorLimit (
	struct work_unit *w)
{
	unsigned long p;
	unsigned int test;

/* If this is trial factoring work with a specified end point, then */
/* return that end_point. */

	if (w->factor_to != 0.0) return ((unsigned int) w->factor_to);

/* For LL tests, determine the optimal trial factoring end point. */
/* This is based on timings from my 2GHz P4. */

	p = w->n;
	if (p > FAC80) test = 80;	/* Test all 80 bit factors */
	else if (p > FAC79) test = 79;	/* Test all 79 bit factors */
	else if (p > FAC78) test = 78;	/* Test all 78 bit factors */
	else if (p > FAC77) test = 77;	/* Test all 77 bit factors */
	else if (p > FAC76) test = 76;	/* Test all 76 bit factors */
	else if (p > FAC75) test = 75;	/* Test all 75 bit factors */
	else if (p > FAC74) test = 74;	/* Test all 74 bit factors */
	else if (p > FAC73) test = 73;	/* Test all 73 bit factors */
	else if (p > FAC72) test = 72;	/* Test all 72 bit factors */
	else if (p > FAC71) test = 71;	/* Test all 71 bit factors */
	else if (p > FAC70) test = 70;	/* Test all 70 bit factors */
	else if (p > FAC69) test = 69;	/* Test all 69 bit factors */
	else if (p > FAC68) test = 68;	/* Test all 68 bit factors */
	else if (p > FAC67) test = 67;	/* Test all 67 bit factors */
	else if (p > FAC66) test = 66;	/* Test all 66 bit factors */
	else if (p > FAC65) test = 65;	/* Test all 65 bit factors */
	else if (p > FAC64) test = 64;	/* Test all 64 bit factors */
	else if (p > FAC63) test = 63;	/* Test all 63 bit factors */
	else if (p > FAC62) test = 62;	/* Test all 62 bit factors */
	else if (p > FAC61) test = 61;	/* Test all 61 bit factors */
	else if (p > FAC60) test = 60;	/* Test all 60 bit factors */
	else if (p > FAC59) test = 59;	/* Test all 59 bit factors */
	else if (p > FAC58) test = 58;	/* Test all 58 bit factors */
	else if (p > FAC57) test = 57;	/* Test all 57 bit factors */
	else if (p > FAC56) test = 56;	/* Test all 56 bit factors */
	else test = 40;			/* Test all 40 bit factors */

/* If double-checking, then trial factor to one less bit depth because */
/* a found factor will only save one LL test, not two. */

	if (w->work_type == WORK_DBLCHK) test--;

/* Return the computed end point. */

	return (test);
}

/**************************************************************/
/*            Routines to compute the rolling average         */
/**************************************************************/

//bug - someway to avoid local.ini updates (to disk) if well_behaved_work
// is set!!

/* Convert a string to a hash value */

unsigned long string_to_hash (
	char	*str)		/* String to hash */
{
	char	md5val[33];
	int	i, j;
	char	*p;
	unsigned long hash, val;

/* Use md5 to generate a number to hash */

	md5 (md5val, str);
	for (i = 0, p = md5val, hash = 0; i < 4; i++) {
		for (j = 0, val = 0; j < 8; j++, p++) {
			val = (val << 4) +
			      (*p < 'A' ? *p - '0' : *p - 'A' + 10);
		}
		hash += val;
	}
	return (hash & 0x7FFFFFFF);
}

/* Build the hash value as we go */

unsigned long build_rolling_hash (
	struct work_unit *w)		/* Work unit to add in to hash */
{
	char	buf[80];

/* Use md5 on a character rep of the number */

	gw_as_string (buf, w->k, w->b, w->n, w->c);
	return (string_to_hash (buf));
}

/* Adjust rolling average computation variables when a work unit completes */

void rolling_average_work_unit_complete (
	int	thread_num,		/* Thread number that completed work */
	struct work_unit *completed_w)
{
	unsigned long hash, time_to_complete;
	struct work_unit *first_w, *second_w;

/* Get current rolling average computation variables */

	hash = IniGetInt (LOCALINI_FILE, "RollingHash", 0);
	time_to_complete = IniGetInt (LOCALINI_FILE, "RollingCompleteTime", 0);

/* Grab the lock so that other threads do not add or delete lines */
/* while we are making this calculation. */

	gwmutex_lock (&WORKTODO_MUTEX);

/* Get first and second work unit pointers */

	for (first_w = WORK_UNITS[thread_num].first;
	     first_w != NULL && first_w->work_type == WORK_NONE;
	     first_w = first_w->next);
	for (second_w = first_w->next;
	     second_w != NULL && second_w->work_type == WORK_NONE;
	     second_w = second_w->next);

/* Only modify the values if we completed the first work_unit and there */
/* is a second work_unit */

	if (first_w == completed_w && second_w != NULL) {

/* Adjust the hash value.  Subtract out the completed work unit and add in */
/* the next work unit. */

		hash -= build_rolling_hash (first_w);
		hash += build_rolling_hash (second_w);

/* Adjust the estimated time to complete the first work unit */

		time_to_complete += (unsigned long)
			work_estimate (thread_num, second_w);
	}

/* Unlock access to the worktodo.ini structures */

	gwmutex_unlock (&WORKTODO_MUTEX);

/* Update the changed rolling average computation values */

	IniWriteInt (LOCALINI_FILE, "RollingHash", hash);
	IniWriteInt (LOCALINI_FILE, "RollingCompleteTime", time_to_complete);
}

/* Adjust the rolling average */

void adjust_rolling_average (void)
{
	unsigned int tnum;
	unsigned long time_to_complete, starting_time_to_complete;
	time_t	current_time, starting_time;
	unsigned long hash, time_in_this_period;
	double	rolling_average_this_period, pct;

/* Grab the lock so that other threads do not add or delete lines */
/* while we are making this calculation. */

	gwmutex_lock (&WORKTODO_MUTEX);

/* Look at the first work unit in each thread */

	hash = 0;
	time_to_complete = 0;
	for (tnum = 0; tnum < NUM_WORKER_THREADS; tnum++) {
		struct work_unit *w;

		for (w = WORK_UNITS[tnum].first; w != NULL; w = w->next)
			if (w->work_type != WORK_NONE &&
			    w->work_type != WORK_DELETED) break;
		if (w == NULL) continue;
	   
/* Build the hash value */

		hash += build_rolling_hash (w);

/* Get the new estimated time to complete this work unit */

		time_to_complete += (unsigned long) work_estimate (tnum, w);
	}

/* Unlock access to the worktodo.ini structures */

	gwmutex_unlock (&WORKTODO_MUTEX);

/* Get the current time and starting time */

	starting_time = IniGetInt (LOCALINI_FILE, "RollingStartTime", 0);
	time (&current_time);

/* Make sure hash codes match.  This protects us against making incorrect */
/* rolling average adjustments when user manually edits worktodo.ini. */

	if (hash != IniGetInt (LOCALINI_FILE, "RollingHash", 0))
		goto no_update;

/* Now compute the rolling average for this time period.  For example, if */
/* the time-to-complete decreased 5 hours over the last 6 hours of elapsed */
/* time, then the rolling average for the last 6 hours was 5/6 of the */
/* current rolling average. */

	time_in_this_period = (unsigned long) (current_time - starting_time);
	starting_time_to_complete =
		IniGetInt (LOCALINI_FILE, "RollingCompleteTime", 0);

	/* Sanity checks for bogus time periods and completion estimates*/

	if (starting_time == 0 ||
	    current_time <= starting_time ||
	    time_in_this_period > 30 * 86400 ||
	    starting_time_to_complete <= time_to_complete)
		goto no_update;

	rolling_average_this_period =
		(double) (starting_time_to_complete - time_to_complete) /
		(double) NUM_WORKER_THREADS /
		(double) time_in_this_period *
		(double) ROLLING_AVERAGE;

	/* If the user is running more worker threads than there are */
	/* CPUs, then adjust the rolling average upward accordingly. */

	if (NUM_WORKER_THREADS > NUM_CPUS)
		rolling_average_this_period *=
			(double) NUM_WORKER_THREADS / (double) NUM_CPUS;

	/* More safeguard against bogus or abruptly changing data */

	if (rolling_average_this_period > 50000.0) goto no_update;
	if (rolling_average_this_period < 0.5 * ROLLING_AVERAGE)
		rolling_average_this_period = 0.5 * ROLLING_AVERAGE;
	if (rolling_average_this_period > 2.0 * ROLLING_AVERAGE)
		rolling_average_this_period = 2.0 * ROLLING_AVERAGE;

/* Calculate the new rolling average - we use a 30-day rolling average */

	pct = (double) time_in_this_period / (30.0 * 86400.0);
	ROLLING_AVERAGE = (unsigned long)
		((1.0 - pct) * ROLLING_AVERAGE +
		 pct * rolling_average_this_period + 0.5);

	/* Safeguards against excessive rolling average values */

	if (ROLLING_AVERAGE < 20) ROLLING_AVERAGE = 20;
	if (ROLLING_AVERAGE > 4000) ROLLING_AVERAGE = 4000;

/* Update rolling average data in the local.ini file */

no_update:
	IniWriteInt (LOCALINI_FILE, "RollingHash", hash);
	IniWriteInt (LOCALINI_FILE, "RollingStartTime", (unsigned long) current_time);
	IniWriteInt (LOCALINI_FILE, "RollingCompleteTime", time_to_complete);
	IniWriteInt (LOCALINI_FILE, "RollingAverage", ROLLING_AVERAGE);
}

/* If a work_unit performed less work than estimated (say by unexpectedly */
/* finding a factor) then do not update the rolling average this period. */

void invalidateNextRollingAverageUpdate (void)
{
	IniWriteInt (LOCALINI_FILE, "RollingStartTime", 0);
	adjust_rolling_average ();
}

/**************************************************************/
/*      Routines to aid in reading and writing save files     */
/**************************************************************/

/* Generate temporary file name */

void tempFileName (
	struct work_unit *w,
	char	*buf)
{
	char	c;
	unsigned long p;

/* WARNING:  Version 24 had a bug where exponents between 61 million and */
/* 70 million used funky characters in the file name.  I've fixed the bug */
/* here, but users will have to rename some intermediate files. */

	p = w->n;
	if (p < 80000000) {
		sprintf (buf, "p%07li", p % 10000000);
		if (p >= 10000000)	/* buf[1] ranges from A-Y */
			buf[1] = (char) ('A' + (p / 1000000) - 10);
		if (p >= 35000000)	/* buf[2] ranges from A-Z */
			c = buf[1], buf[1] = buf[2], buf[2] = (char)(c - 25);
		if (p >= 61000000)	/* buf[3] ranges from B-T */
			c = buf[2], buf[2] = buf[3], buf[3] = (char)(c - 25);
	} else
		sprintf (buf, "p%ld", p);

/* Use different first letters for different work types.  This isn't */
/* completely compatible with v24 (the -An extension and P-1 and ECM when c=1 */

	if (w->work_type == WORK_FACTOR) buf[0] = 'f';
	if (w->work_type == WORK_ECM) buf[0] = 'e';
	if (w->work_type == WORK_PMINUS1 || w->work_type == WORK_PFACTOR) buf[0] = 'm';

/* Prior to version 25.9 build 5, the pfactor work type used p as the */
/* first letter, we now use m.  To reduce upgrading problems, old save */
/* file names are renamed. */

	if (w->work_type == WORK_PFACTOR) {
		char	v258_filename[32];
		sprintf (v258_filename, "p%s", buf+1);
		rename (v258_filename, buf);
	}

/* Prior to version 25.9 build 4, if c was 1 then P-1 and ECM used */
/* a different first letter in the filename.  From now on, we will no */
/* longer do this.  To reduce upgrading problems, old save file names */
/* are renamed. */

	if (w->c == 1 && buf[0] == 'm') {
		char	v258_filename[32];
		sprintf (v258_filename, "l%s", buf+1);
		rename (v258_filename, buf);
	}
	if (w->c == 1 && buf[0] == 'e') {
		char	v258_filename[32];
		sprintf (v258_filename, "d%s", buf+1);
		rename (v258_filename, buf);
	}

/* Prior to version 25.9 build 4, we did not use k or c in generating the */
/* filename.  Thus, 10223*2^11111111+1 and 67607*2^11111111+1 would both */
/* use the save save file -- a definite problem for Seventeen or Bust. */
/* From now on, we will use k and c to generate the filename.  To reduce */
/* upgrading problems, old save file names are renamed. */

	if (w->k != 1.0 || abs(w->c) != 1) {
		char	v258_filename[32];
		strcpy (v258_filename, buf);
		buf[1] = 0;
		if (w->k != 1.0) sprintf (buf+strlen(buf), "%g", fmod (w->k, 1000000.0));
		sprintf (buf+strlen(buf), "_%ld", p);
		if (abs(w->c) != 1) sprintf (buf+strlen(buf), "_%d", abs(w->c) % 1000);
		rename (v258_filename, buf);
		if (buf[0] == 'p') {
			v258_filename[0] = buf[0] = 'q';
			rename (v258_filename, buf);
			buf[0] = 'p';
		}
	}

/* Append extension */

	if (w->extension[0]) {
		strcat (buf, ".");
		strcat (buf, w->extension);
	}
}

/* See if the given file exists */

int fileExists (
	char	*filename)
{
	int	fd;
	fd = _open (filename, _O_RDONLY | _O_BINARY);
	if (fd < 0) return (0);
	_close (fd);
	return (1);
}

/* Routines to read and write a byte array from and to a save file */

int read_array (
	int	fd,
	char	*buf,
	unsigned long len,
	unsigned long *sum)
{
	unsigned long i;
	unsigned char *ubuf;

	if (_read (fd, buf, len) != len) return (FALSE);
	ubuf = (unsigned char *) buf;
	if (sum != NULL)
		for (i = 0; i < len; i++)
			*sum = (uint32_t) (*sum + ubuf[i]);
	return (TRUE);
}

int write_array (
	int	fd,
	char	*buf,
	unsigned long len,
	unsigned long *sum)
{
	unsigned long i;
	unsigned char *ubuf;

	if (len == 0) return (TRUE);
	if (_write (fd, buf, len) != len) return (FALSE);
	ubuf = (unsigned char *) buf;
	if (sum != NULL)
		for (i = 0; i < len; i++)
			*sum = (uint32_t) (*sum + ubuf[i]);
	return (TRUE);
}

/* Routines to read and write a gwnum from and to a save file */

int read_gwnum (
	int	fd,
	gwhandle *gwdata,
	gwnum	g,
	unsigned long *sum)
{
	giant	tmp;
	unsigned long i, len, giantlen, bytes;

	if (!read_long (fd, &len, sum)) return (FALSE);
	if (len == 0) return (FALSE);

	giantlen = ((int) gwdata->bit_length >> 5) + 10;
	if (len > giantlen) return (FALSE);
	tmp = popg (&gwdata->gdata, giantlen);
	if (tmp == NULL) return (FALSE);	// BUG - we should return some other error code
						// otherwise caller will likely delete save file.

	bytes = len * sizeof (uint32_t);
	if (_read (fd, tmp->n, bytes) != bytes) goto errexit;
	if (len && tmp->n[len-1] == 0) goto errexit;
	tmp->sign = len;
	*sum = (uint32_t) (*sum + len);
	for (i = 0; i < len; i++) *sum = (uint32_t) (*sum + tmp->n[i]);
	gianttogw (gwdata, tmp, g);
	pushg (&gwdata->gdata, 1);
	return (TRUE);

// Free memory and return failure

errexit:
	pushg (&gwdata->gdata, 1);
	return (FALSE);
}

int write_gwnum (
	int	fd,
	gwhandle *gwdata,
	gwnum	g,
	unsigned long *sum)
{
	giant	tmp;
	unsigned long i, len, bytes;

	tmp = popg (&gwdata->gdata, ((int) gwdata->bit_length >> 5) + 10);
	if (tmp == NULL) return (FALSE);
	gwtogiant (gwdata, g, tmp);
	len = tmp->sign;
	if (len == 0) return (FALSE);
	if (!write_long (fd, len, sum)) return (FALSE);
	bytes = len * sizeof (uint32_t);
	if (_write (fd, tmp->n, bytes) != bytes) return (FALSE);
	*sum = (uint32_t) (*sum + len);
	for (i = 0; i < len; i++) *sum = (uint32_t) (*sum + tmp->n[i]);
	pushg (&gwdata->gdata, 1);
	return (TRUE);
}

/* Routines to read and write values from and to a save file */

int read_short (			/* Used for old-style save files */
	int	fd,
	short	*val)
{
	if (_read (fd, val, sizeof (short)) != sizeof (short)) return (FALSE);
	return (TRUE);
}

int read_long (
	int	fd,
	unsigned long *val,
	unsigned long *sum)
{
	uint32_t tmp;

	if (_read (fd, &tmp, sizeof (uint32_t)) != sizeof (uint32_t))
		return (FALSE);
	if (sum != NULL) *sum = (uint32_t) (*sum + tmp);
	*val = tmp;
	return (TRUE);
}

int write_long (
	int	fd,
	unsigned long val,
	unsigned long *sum)
{
	uint32_t tmp;

	tmp = (uint32_t) val;
	if (_write (fd, &tmp, sizeof (uint32_t)) != sizeof (uint32_t))
		return (FALSE);
	if (sum != NULL) *sum = (uint32_t) (*sum + tmp);
	return (TRUE);
}

int read_slong (
	int	fd,
	long	*val,
	unsigned long *sum)
{
	int32_t tmp;

	if (_read (fd, &tmp, sizeof (int32_t)) != sizeof (int32_t))
		return (FALSE);
	if (sum != NULL) *sum = (uint32_t) (*sum + (uint32_t) tmp);
	*val = tmp;
	return (TRUE);
}

int write_slong (
	int	fd,
	long	val,
	unsigned long *sum)
{
	int32_t tmp;

	tmp = (int32_t) val;
	if (_write (fd, &tmp, sizeof (int32_t)) != sizeof (int32_t))
		return (FALSE);
	if (sum != NULL) *sum = (uint32_t) (*sum + (uint32_t) tmp);
	return (TRUE);
}

int read_longlong (
	int	fd,
	uint64_t *val,
	unsigned long *sum)
{
	if (_read (fd, val, sizeof (uint64_t)) != sizeof (uint64_t))
		return (FALSE);
	if (sum != NULL) *sum = (uint32_t) (*sum + (*val >> 32) + *val);
	return (TRUE);
}

int write_longlong (
	int	fd,
	uint64_t val,
	unsigned long *sum)
{
	if (_write (fd, &val, sizeof (uint64_t)) != sizeof (uint64_t))
		return (FALSE);
	if (sum != NULL) *sum = (uint32_t) (*sum + (val >> 32) + val);
	return (TRUE);
}

int read_double (
	int	fd,
	double	*val,
	unsigned long *sum)
{
	if (_read (fd, val, sizeof (double)) != sizeof (double))
		return (FALSE);
	if (sum != NULL) *sum = (uint32_t) (*sum + (uint32_t) *val);
	return (TRUE);
}

int write_double (
	int	fd,
	double	val,
	unsigned long *sum)
{
	if (_write (fd, &val, sizeof (double)) != sizeof (double))
		return (FALSE);
	if (sum != NULL) *sum += (uint32_t) (*sum + (uint32_t) val);
	return (TRUE);
}

/* Routines to read and write the common header portion of all save files */
/* The save file format is: */
/*	u32		magic number  (different for ll, p-1, prp, tf, ecm) */
/*	u32		version number */
/*	double		k in k*b^n+c */
/*	u32		b in k*b^n+c */
/*	u32		n in k*b^n+c */
/*	s32		c in k*b^n+c */
/*	double		pct complete */
/*	char(11)	stage */
/*	char(1)		pad */
/*	u32		checksum of all following data */

/* Read and test the "magic number" in the first 4 bytes of the file. */
/* We use this to detect save files created by this program prior to the */
/* common header format. */

int read_magicnum (
	int	fd,
	unsigned long magicnum)
{
	unsigned long filenum;

/* Read the magic number from the first 4 bytes */

	_lseek (fd, 0, SEEK_SET);
	if (!read_long (fd, &filenum, NULL)) return (FALSE);

/* Return TRUE if the magic number matches the caller's desired magic number */

	return (filenum == magicnum);
}

/* Read the rest of the common header */

int read_header (
	int	fd,
	unsigned long *version,
	struct work_unit *w,
	unsigned long *sum)
{
	double	k;
	unsigned long b, n;
	long	c;
	char	pad;
	char	stage[11];
	double	pct_complete;
	unsigned long trash_sum;

/* Skip past the magic number in the first 4 bytes */

	_lseek (fd, sizeof (uint32_t), SEEK_SET);

/* Read the header */

	if (!read_long (fd, version, NULL)) return (FALSE);
	if (!read_double (fd, &k, NULL)) return (FALSE);
	if (!read_long (fd, &b, NULL)) return (FALSE);
	if (!read_long (fd, &n, NULL)) return (FALSE);
	if (!read_slong (fd, &c, NULL)) return (FALSE);
	if (!read_array (fd, stage, 11, NULL)) return (FALSE);
	if (!read_array (fd, &pad, 1, NULL)) return (FALSE);
	if (!read_double (fd, &pct_complete, NULL)) return (FALSE);
	if (sum == NULL) sum = &trash_sum;
	if (!read_long (fd, sum, NULL)) return (FALSE);

/* Validate the k,b,n,c values */

	if (k != w->k || b != w->b || n != w->n || c != w->c) return (FALSE);

/* Set the work unit's stage and pct_complete fields */

	stage[10] = 0;
	strcpy (w->stage, stage);
	if (pct_complete < 0.0) pct_complete = 0.0;
	if (pct_complete > 1.0) pct_complete = 1.0;
	w->pct_complete = pct_complete;

/* Return success */

	return (TRUE);
}

int write_header (
	int	fd,
	unsigned long magicnum,
	unsigned long version,
	struct work_unit *w)
{
	char	pad = 0;
	uint32_t sum = 0;

	if (!write_long (fd, magicnum, NULL)) return (FALSE);
	if (!write_long (fd, version, NULL)) return (FALSE);
	if (!write_double (fd, w->k, NULL)) return (FALSE);
	if (!write_long (fd, w->b, NULL)) return (FALSE);
	if (!write_long (fd, w->n, NULL)) return (FALSE);
	if (!write_slong (fd, w->c, NULL)) return (FALSE);
	if (!write_array (fd, w->stage, 11, NULL)) return (FALSE);
	if (!write_array (fd, &pad, 1, NULL)) return (FALSE);
	if (!write_double (fd, w->pct_complete, NULL)) return (FALSE);
	if (!write_long (fd, sum, NULL)) return (FALSE);
	return (TRUE);
}

#define CHECKSUM_OFFSET	48

int read_checksum (
	int	fd,
	unsigned long *sum)
{
	_lseek (fd, CHECKSUM_OFFSET, SEEK_SET);
	if (!read_long (fd, sum, NULL)) return (FALSE);
	return (TRUE);
}

int write_checksum (
	int	fd,
	unsigned long sum)
{
	_lseek (fd, CHECKSUM_OFFSET, SEEK_SET);
	if (!write_long (fd, sum, NULL)) return (FALSE);
	return (TRUE);
}



/* Format a message for writing to the results file and sending to the */
/* server.  We prepend the assignment ID to the message.  If operating */
/* offline, then prepend other information. */

void formatMsgForResultsFile (
	char	*buf,		/* Msg to prepend to, 200 characters max */
	struct work_unit *w)
{
	char	newbuf[2000];

/* Output a USERID/COMPID prefix for result messages */

	if (!USERID[0])
		strcpy (newbuf, buf);
	else if (!COMPID[0])
		sprintf (newbuf, "UID: %s, %s", USERID, buf);
	else
		sprintf (newbuf, "UID: %s/%s, %s", USERID, COMPID, buf);

/* Output the assignment ID too.  That alone should be enough info to */
/* credit the correct userID.  However, we still output the user ID as it */
/* is far more human friendly than an assignment ID. */

	if (w->assignment_uid[0])
		sprintf (newbuf + strlen (newbuf) - 1,
			 ", AID: %s\n", w->assignment_uid);

/* Now truncate the message to 200 characters */

	newbuf[200] = 0;
	strcpy (buf, newbuf);
}

/* Open the results file and write a line to the end of it. */

int writeResults (
	const char *msg)
{
static	time_t	last_time = 0;
	time_t	this_time;
	int	fd;

/* Open file, position to end */

	gwmutex_lock (&OUTPUT_MUTEX);
	fd = _open (RESFILE, _O_TEXT | _O_RDWR | _O_CREAT | _O_APPEND, CREATE_FILE_ACCESS);
	if (fd < 0) {
		gwmutex_unlock (&OUTPUT_MUTEX);
		LogMsg ("Error opening results file to output this message:\n");
		LogMsg (msg);
		return (FALSE);
	}

/* If it has been at least 5 minutes since the last time stamp */
/* was output, then output a new timestamp */

	time (&this_time);
	if (this_time - last_time > 300) {
		char	buf[32];
		last_time = this_time;
		buf[0] = '[';
		strcpy (buf+1, ctime (&this_time));
		buf[25] = ']';
		buf[26] = '\n';
		(void) _write (fd, buf, 27);
	}

/* Output the message */

	if (_write (fd, msg, (unsigned int) strlen (msg)) < 0) goto fail;
	_close (fd);
	gwmutex_unlock (&OUTPUT_MUTEX);
	return (TRUE);

/* On a write error, close file and return error flag */

fail:	_close (fd);
	gwmutex_unlock (&OUTPUT_MUTEX);
	LogMsg ("Error writing message to results file:\n");
	LogMsg (msg);
	return (FALSE);
}


/****************************************************************************/
/*               Spool File and Server Communication Code                   */
/****************************************************************************/

#define	SPOOL_FILE_MAGICNUM	0x73d392ac
#define SPOOL_FILE_VERSION	1
/* Offset to the header words (just past the magicnum and version num) */
#define SPOOL_FILE_HEADER_OFFSET (2 * sizeof (uint32_t))
/* Offset to messages (past magicnum, version num, and two header words) */
#define SPOOL_FILE_MSG_OFFSET (4 * sizeof (uint32_t))

gwthread COMMUNICATION_THREAD = 0;	/* Handle for comm thread */
gwmutex	SPOOL_FILE_MUTEX;		/* Lock governing spool file access */
int	GET_PING_INFO = 0;		/* Flag to get ping info */
int	GLOBAL_SEND_MSG_COUNT = 0;	/* Used to detect hung comm threads */
struct work_unit *LOCKED_WORK_UNIT = NULL; /* Work unit to unlock if comm */
					/* thread hangs. */

/* Spool file header word flags */

#define HEADER_FLAG_MSGS	0x0001	/* informational msgs in file */
#define HEADER_FLAG_PO		0x0004	/* exchange program options */
#define HEADER_FLAG_END_DATES	0x0008	/* completion dates need sending */
#define HEADER_FLAG_QUIT_GIMPS	0x0010	/* return all work (quit gimps) */
#define HEADER_FLAG_UC		0x0020	/* computer info has changed */
#define HEADER_FLAG_WORK_QUEUE	0x0040	/* check if enough work is queued */

void communicateWithServer (void *arg);

/* Init the spool file and communication code */

void init_spool_file_and_comm_code (void)
{
	gwmutex_init (&SPOOL_FILE_MUTEX);
	GET_PING_INFO = 0;
}

/* Add or delete the comm timers based on the communication global */
/* variables.  This is called at start up and whenever the user toggles */
/* the USE_PRIMENET or MANUAL_COMM global variables. */

void set_comm_timers (void)
{
	time_t	last_time, current_time;

/* If we are not doing automatic server communications, then make sure */
/* no communication timers are active.  If, by chance, we are in the */
/* middle of communicating with the server, then let the thread complete. */

//bug - should we destroy the comm window if !USE_PRIMENET??? User may
//have interesting data there, however it is in prime.log too.
//Maybe just change the window title.
// or do we need to leave it active until a Quit GIMPS succeeds??
	if (!USE_PRIMENET) {
		delete_timed_event (TE_COMM_SERVER);
		delete_timed_event (TE_COMPLETION_DATES);
		delete_timed_event (TE_WORK_QUEUE_CHECK);
		return;
	}
	if (MANUAL_COMM)
		delete_timed_event (TE_COMM_SERVER);

/* Create and name all the communication window. */

	create_window (COMM_THREAD_NUM);
	base_title (COMM_THREAD_NUM, "Communication thread");
	ChangeIcon (COMM_THREAD_NUM, IDLE_ICON);
	title (COMM_THREAD_NUM, "Inactive");

/* Update completion dates on the server if it has been a month since we */
/* last updated the server.  Otherwise, start a timer to make this happen */
/* at the appropriate time. */

/* Get the current time and when the completion dates were last sent */

	time (&current_time);
	last_time = IniGetInt (LOCALINI_FILE, "LastEndDatesSent", 0);

/* If it's been the correct number of days, then update the end dates */

	if (current_time < last_time ||
	    current_time > (time_t)(last_time + DAYS_BETWEEN_CHECKINS * 86400.0))
		UpdateEndDates ();
	else
		add_timed_event (TE_COMPLETION_DATES,
				 (int) (last_time +
					DAYS_BETWEEN_CHECKINS * 86400.0 -
					current_time));

/* Add the event that checks if the work queue has enough work every 6 */
/* hours.  As a side effect, this will also start the comm thread in case */
/* there is an old spool file hanging around. */

	add_timed_event (TE_WORK_QUEUE_CHECK, 5);  /* Start in 5 seconds */
}

/* Routine to fire up the communication thread in response to a user */
/* request to communicate with the server now. */

void do_manual_comm_now (void)
{
	gwmutex_lock (&SPOOL_FILE_MUTEX);
	if (!COMMUNICATION_THREAD)
		gwthread_create (&COMMUNICATION_THREAD,
				 &communicateWithServer, NULL);
	gwmutex_unlock (&SPOOL_FILE_MUTEX);
}

/* Clear rate-limiting counters and timers.  Any time the user explicitly */
/* chooses Test/Continue, we reset the rate limits.  We do this so that */
/* if there is an explicit need to communicate with the server frequently */
/* in the short term, the user has a way to do it.  The rate limits are */
/* here to guard against runaway clients from pummeling the server. */

void clear_comm_rate_limits (void)
{
//bug - clear other rate limiters here....

/* If we are paused for an hour because of a failed connection attempt, */
/* then kill the timer and comm now.  The user may have edited the proxy */
/* info in the INI file so that comm will now work. */

	gwmutex_lock (&SPOOL_FILE_MUTEX);
	if (!MANUAL_COMM &&
	    !COMMUNICATION_THREAD &&
	    is_timed_event_active (TE_COMM_SERVER)) {
		delete_timed_event (TE_COMM_SERVER);
		gwthread_create (&COMMUNICATION_THREAD,
				 &communicateWithServer, NULL);
	}
	gwmutex_unlock (&SPOOL_FILE_MUTEX);
}

/* Ping the server and retrieve info for an About box. */
/* The communication thread will call pingServerResponse with the results. */

void pingServer (void)
{

/* Set global variable indicating we want to get ping information. */
/* Then fire up the communication thread. */

	gwmutex_lock (&SPOOL_FILE_MUTEX);
	GET_PING_INFO = 1;
	if (!COMMUNICATION_THREAD)
		gwthread_create (&COMMUNICATION_THREAD,
				 &communicateWithServer, NULL);
	gwmutex_unlock (&SPOOL_FILE_MUTEX);
}

/* Update completion dates on the server.  Set a flag in */
/* the spool file saying this is necessary. */

void UpdateEndDates (void)
{
	spoolMessage (PRIMENET_ASSIGNMENT_PROGRESS, NULL);
}

/* Write a message to the spool file */

void spoolMessage (
	short	msgType,
	void	*msg)
{
	int	fd;
	unsigned long magicnum, version;
	unsigned long header_word;

/* If we're not using primenet, ignore this call */

	if (!USE_PRIMENET) return;

/* Obtain mutex before accessing spool file */

	gwmutex_lock (&SPOOL_FILE_MUTEX);

/* Open the spool file */

	fd = _open (SPOOL_FILE, _O_RDWR | _O_BINARY | _O_CREAT, CREATE_FILE_ACCESS);
	if (fd < 0) {
		LogMsg ("ERROR: Unable to open spool file.\n");
		gwmutex_unlock (&SPOOL_FILE_MUTEX);
		return;
	}

/* If the file is empty, write the spool file header */

	if (!read_long (fd, &magicnum, NULL)) {
		write_long (fd, SPOOL_FILE_MAGICNUM, NULL);
		write_long (fd, SPOOL_FILE_VERSION, NULL);
		write_long (fd, 0, NULL);
		write_long (fd, 0, NULL);
		header_word = 0;
	}

/* Otherwise, read and validate header.  If it is bad try to salvage */
/* the spool file data. */

	else if (magicnum != SPOOL_FILE_MAGICNUM ||
		 !read_long (fd, &version, NULL) ||
		 version != SPOOL_FILE_VERSION ||
		 !read_long (fd, &header_word, NULL)) {
		_close (fd);
		gwmutex_unlock (&SPOOL_FILE_MUTEX);
		salvageCorruptSpoolFile ();
		spoolMessage (msgType, msg);
		return;
	}

/* If this is a message telling us to check if enough work is queued up, */
/* then set the proper bit in the header word. */

	if (msgType == MSG_CHECK_WORK_QUEUE)
		header_word |= HEADER_FLAG_WORK_QUEUE;

/* If this is a maintain user info message, then set the header */
/* word appropriately.  At Scott's request also send computer info. */

	else if (msgType == PRIMENET_UPDATE_COMPUTER_INFO) {
		header_word |= HEADER_FLAG_UC;
	}

/* If this is a exchange program options message, then set the header */
/* word appropriately. */

	else if (msgType == PRIMENET_PROGRAM_OPTIONS) {
		header_word |= HEADER_FLAG_PO;
	}

/* If this is an update completion dates message, then set the header word */
/* appropriately.  At Scott's request also send computer info. */

	else if (msgType == PRIMENET_ASSIGNMENT_PROGRESS)
		header_word |= HEADER_FLAG_END_DATES + HEADER_FLAG_UC;

/* Ugly little hack when quitting GIMPS */

	else if (msgType == MSG_QUIT_GIMPS)
		header_word |= HEADER_FLAG_QUIT_GIMPS;

/* Otherwise this is a result, unreserve, or benchmark data */

	else
		header_word |= HEADER_FLAG_MSGS;

/* Write the new header word */

	_lseek (fd, SPOOL_FILE_HEADER_OFFSET, SEEK_SET);
	write_long (fd, header_word, NULL);

/* Write out a full message */

	if (msgType == PRIMENET_ASSIGNMENT_RESULT ||
	    msgType == PRIMENET_ASSIGNMENT_UNRESERVE ||
	    msgType == PRIMENET_BENCHMARK_DATA) {
		char	buf[1024];
		short	datalen;

/* Skip the remaining messages */

		while (_read (fd, buf, sizeof (buf)));

/* Append the latest message */

		datalen = (msgType == PRIMENET_ASSIGNMENT_RESULT) ?
			sizeof (struct primenetAssignmentResult) :
			(msgType == PRIMENET_ASSIGNMENT_UNRESERVE) ?
			sizeof (struct primenetAssignmentUnreserve) :
			sizeof (struct primenetBenchmarkData);
		(void) _write (fd, &msgType, sizeof (short));
		(void) _write (fd, &datalen, sizeof (short));
		(void) _write (fd, msg, datalen);
	}

/* Close the spool file */

	_close (fd);

/* Fire up the communication thread if we are not waiting for the user to */
/* complete the initial startup dialog boxes and we are not doing manual */
/* communication and the communication thread is not already active and we */
/* are not waiting some time to retry after a failed communication attempt. */

	if (!STARTUP_IN_PROGRESS &&
	    !MANUAL_COMM &&
	    !COMMUNICATION_THREAD &&
	    !is_timed_event_active (TE_COMM_SERVER))
		gwthread_create (&COMMUNICATION_THREAD,
				 &communicateWithServer, NULL);

/* Release mutex before accessing spool file */

	gwmutex_unlock (&SPOOL_FILE_MUTEX);
}

/* Copy an existing results file to the spool file */
/* This is used when converting from manual to automatic mode */
/* It provides an extra chance that the existing results file */
/* will get sent to us. */

void spoolExistingResultsFile (void)
{
	int	i;
	char	*filename;
	char	line[256];
	FILE	*fd;

	for (i = 1; i <= 2; i++) {
		if (i == 1) filename = "results.txt";
		if (i == 2) {
			if (!strcmp (RESFILE, "results.txt")) continue;
			filename = RESFILE;
		}
		fd = fopen (filename, "r");
		if (fd == NULL) continue;
		while (fgets (line, sizeof (line) - 1, fd)) {
			if (line[0] == '[') continue;
			if (strstr (line, "Res64") == NULL &&
			    strstr (line, "factor:") == NULL &&
			    strstr (line, "completed P-1") == NULL &&
			    strstr (line, "no factor") == NULL) continue;
			if (line[0] == 'U' && strchr (line, ',') != NULL)
				safe_strcpy (line, strchr (line, ',') + 2);
//BUG - what to do with this??? Pass an unknown result type which
//forces it to be treated like a manual result?  Does v5 understand
//v4 result strings?
//			spoolMessage (PRIMENET_RESULT_MESSAGE, line);
		}
		fclose (fd);
	}
}

/* Unreserve an exponent */

int unreserve (
	unsigned long p)
{
	unsigned int tnum;
	int	rc, found_one;

/* Find exponent in worktodo.ini and delete it if present */

	found_one = FALSE;
	for (tnum = 0; tnum < NUM_WORKER_THREADS; tnum++) {
	    struct work_unit *w;

	    w = NULL;
	    for ( ; ; ) {

/* Read the line of the work file */

		w = getNextWorkToDoLine (tnum, w, SHORT_TERM_USE);
		if (w == NULL) break;
		if (w->work_type == WORK_NONE) continue;

/* Skip the line if exponent is not a match */

		if (w->n != p) continue;
		found_one = TRUE;

/* Build a packet and spool message */

		if (w->assignment_uid[0]) {
			struct primenetAssignmentUnreserve pkt;
			memset (&pkt, 0, sizeof (pkt));
			strcpy (pkt.computer_guid, COMPUTER_GUID);
			strcpy (pkt.assignment_uid, w->assignment_uid);
			spoolMessage (PRIMENET_ASSIGNMENT_UNRESERVE, &pkt);
		}

/* Delete the line.  The work unit will immediately be removed from the */
/* worktodo.ini file and will be deleted from the in-memory structures */
/* once the last in-use lock is released. */

		rc = deleteWorkToDoLine (tnum, w, TRUE);
		if (rc) return (rc);
	    }
	}

/* If we didn't find a match then output a message to the main window */

	if (!found_one) {
		char	buf[90];
		sprintf (buf, "Error unreserving exponent: %lu not found in worktodo.txt\n", p);
		OutputStr (MAIN_THREAD_NUM, buf);
	}

/* Return successfully */

	return (0);
}

/* Output message to screen and prime.log file */

void LogMsg (
	const char *str)
{
	int	fd;
	unsigned long filelen;
static	time_t	last_time = 0;
	time_t	this_time;

/* Output it to the screen */

	OutputStr (COMM_THREAD_NUM, str);

/* Open the log file and position to the end */

	gwmutex_lock (&LOG_MUTEX);
	fd = _open (LOGFILE, _O_TEXT | _O_RDWR | _O_CREAT, CREATE_FILE_ACCESS);
	if (fd < 0) {
		gwmutex_unlock (&LOG_MUTEX);
		OutputStr (COMM_THREAD_NUM, "Unable to open log file.\n");
		return;
	}
	filelen = _lseek (fd, 0L, SEEK_END);

/* If the log file has grown too big, lose the first 100,000 bytes */

	if (filelen > (unsigned long) IniGetInt (INI_FILE, "MaxLogFileSize", 2000000)) {
		char	*buf, *p;
		int	bytes_read;

		buf = (char *) malloc (filelen);
		if (buf != NULL) {
			_lseek (fd, 100000L, SEEK_SET);
			strcpy (buf, "Prior log file entries removed.\n");
			for (p = buf + strlen (buf); (bytes_read = _read (fd, p, 50000)) != 0; p += bytes_read)
				/*do nothing*/;
		       	_close (fd);
			fd = _open (LOGFILE, _O_TEXT | _O_RDWR | _O_TRUNC, CREATE_FILE_ACCESS);
			if (fd < 0) {
				free (buf);
				gwmutex_unlock (&LOG_MUTEX);
				OutputStr (COMM_THREAD_NUM, "Unable to truncate log file.\n");
				return;
			}
			(void) _write (fd, buf, (unsigned int) (p - buf));
			free (buf);
		}
	}

/* If it has been at least 5 minutes since the last time stamp */
/* was output, then output a new timestamp */

	time (&this_time);
	if (this_time - last_time > 300) {
		char	buf[48];
		last_time = this_time;
		buf[0] = '[';
		strcpy (buf+1, ctime (&this_time));
		sprintf (buf+25, " - ver %s]\n", VERSION);
		(void) _write (fd, buf, (unsigned int) strlen (buf));
	}

/* Output the message */

	(void) _write (fd, str, (unsigned int) strlen (str));
	_close (fd);
	gwmutex_unlock (&LOG_MUTEX);
}

/* Format text for prime.log */

void kbnc_to_text (
	char	*buf,
	int	primenet_work_type,
	double	k,
	unsigned long b,
	unsigned long n,
	long	c)
{
	char	num[80], *work_type_str;

	switch (primenet_work_type) {
	case PRIMENET_WORK_TYPE_FACTOR:
		work_type_str = "Trial factor";
		k = 1.0; b = 2; c = -1;
		break;
	case PRIMENET_WORK_TYPE_FIRST_LL:
		work_type_str = "LL";
		k = 1.0; b = 2; c = -1;
		break;
	case PRIMENET_WORK_TYPE_DBLCHK:
		work_type_str = "Double check";
		k = 1.0; b = 2; c = -1;
		break;
	case PRIMENET_WORK_TYPE_ECM:
		work_type_str = "ECM";
		break;
	case PRIMENET_WORK_TYPE_PFACTOR:
		work_type_str = "P-1";
		k = 1.0; b = 2; c = -1;
		break;
	case PRIMENET_WORK_TYPE_PMINUS1:
		work_type_str = "P-1";
		break;
	case PRIMENET_WORK_TYPE_PRP:
		work_type_str = "PRP";
		break;
	default:
		work_type_str = "Unknown work type";
		break;
	}
	gw_as_string (num, k, b, n, c);
	sprintf (buf, "%s %s", work_type_str, num);
}

void aid_to_text (
	char	*buf,
	char	*aid)
{
	unsigned int tnum;
	struct work_unit *w;

/* Scan all work units until we find a matching assignment id */

	for (tnum = 0; tnum < NUM_WORKER_THREADS; tnum++) {
	    w = NULL;
	    for ( ; ; ) {
		w = getNextWorkToDoLine (tnum, w, SHORT_TERM_USE);
		if (w == NULL) break;
		if (w->work_type == WORK_NONE) continue;
		if (strcmp (aid, w->assignment_uid)) continue;
		gw_as_string (buf, w->k, w->b, w->n, w->c);
		decrementWorkUnitUseCount (w, SHORT_TERM_USE);
		return;
	    }
	}
	sprintf (buf, "assignment %s", aid);
}


/* Read a spooled message */

void readMessage (
	int	fd,
	long	*offset,	/* Offset of the message */
	short	*msgType,	/* 0 = no message */
	void	*msg)
{
	short	datalen;

/* Loop until a message that hasn't already been sent is found */

	for ( ; ; ) {
		*offset = _lseek (fd, 0, SEEK_CUR);
		if (_read (fd, msgType, sizeof (short)) != sizeof (short))
			break;
		if (_read (fd, &datalen, sizeof (short)) != sizeof (short))
			break;

/* Read the body of the message */

		if (_read (fd, msg, datalen) != datalen) break;

/* Loop if message has already been sent */

		if (*msgType == -1) continue;

/* Return if msgType is one we expected. */

		if (*msgType == PRIMENET_ASSIGNMENT_RESULT ||
		    *msgType == PRIMENET_ASSIGNMENT_UNRESERVE ||
		    *msgType == PRIMENET_BENCHMARK_DATA)
			return;

/* MsgType was unexpected.  This indicates a corrupt spool file. */

		LogMsg ("Corrupt spool file.  Message ignored.\n");
	}

/* On file read errors (or EOF), return code indicating we're done reading */
/* the spool file. */
	
	*msgType = 0;
}


/* Send a message that was read from the spool file */

int sendMessage (
	short	msgType,
	void	*msg)
{
	struct primenetAssignmentResult *pkt;
	int	local_send_msg_count;
	int	return_code;
	char	buf[400], info[200];

/* Load the primenet library and do any special handling */
/* required prior to calling primenet */

	if (!LoadPrimeNet ()) return (PRIMENET_ERROR_MODEM_OFF);

/* Prepend all messages with the userid and computer id */

	pkt = (struct primenetAssignmentResult *) msg;

/* Print a message on the screen and in the log file */

	switch (msgType) {
	case PRIMENET_UPDATE_COMPUTER_INFO:
		LogMsg ("Updating computer information on the server\n");
		break;
	case PRIMENET_PROGRAM_OPTIONS:
		LogMsg ("Exchanging program options with server\n");
		break;
	case PRIMENET_PING_SERVER:
		OutputStr (COMM_THREAD_NUM, "Contacting PrimeNet Server.\n");
		break;
	case PRIMENET_GET_ASSIGNMENT:
		LogMsg ("Getting assignment from server\n");
		break;
	case PRIMENET_REGISTER_ASSIGNMENT:
		kbnc_to_text (info,
			      ((struct primenetRegisterAssignment *)pkt)->work_type,
			      ((struct primenetRegisterAssignment *)pkt)->k,
			      ((struct primenetRegisterAssignment *)pkt)->b,
			      ((struct primenetRegisterAssignment *)pkt)->n,
			      ((struct primenetRegisterAssignment *)pkt)->c);
		sprintf (buf, "Registering assignment: %s\n", info);
		LogMsg (buf);
		break;
	case PRIMENET_ASSIGNMENT_PROGRESS:
		{
			time_t	this_time;
			char	timebuf[30];
			time (&this_time);
			this_time += ((struct primenetAssignmentProgress *)pkt)->end_date;
			strcpy (timebuf, ctime (&this_time)+4);
			safe_strcpy (timebuf+6, timebuf+15);
			aid_to_text (info, ((struct primenetAssignmentProgress *)pkt)->assignment_uid);
			sprintf (buf, "Sending expected completion date for %s: %s",
				 info, timebuf);
			LogMsg (buf);
		}
		break;
	case PRIMENET_ASSIGNMENT_UNRESERVE:
		aid_to_text (info, ((struct primenetAssignmentUnreserve *)pkt)->assignment_uid);
		sprintf (buf, "Unreserving %s\n", info);
		LogMsg (buf);
		break;
	case PRIMENET_ASSIGNMENT_RESULT:
		sprintf (buf, "Sending result to server: %s\n",
			 ((struct primenetAssignmentResult *)pkt)->message);
		LogMsg (buf);
		break;
	case PRIMENET_BENCHMARK_DATA:
		LogMsg ("Sending benchmark data to server\n");
		break;
	}

/* Fill in the common header fields */

	pkt->versionNumber = PRIMENET_VERSION;

/* Send the message.  Kill the comm thread if server hasn't responded */
/* in 15 minutes. */

	local_send_msg_count = ++GLOBAL_SEND_MSG_COUNT;
	add_timed_event (TE_COMM_KILL, 15*60);
	return_code = PRIMENET (msgType, pkt);
	delete_timed_event (TE_COMM_KILL);

/* If the kill thread timer fired because our communication with the */
/* server hung, yet somehow the thread magically got unstuck, then return */
/* a dummy error code.  This is necessary because we have decremented the */
/* use count of the work_unit we are sending a completion date on. */

	if (GLOBAL_SEND_MSG_COUNT != local_send_msg_count) return (9999);

/* Print a result message on the screen and in the log file */

	if (return_code == 0)
	switch (msgType) {
	case PRIMENET_GET_ASSIGNMENT:
		kbnc_to_text (info,
			      ((struct primenetGetAssignment *)pkt)->work_type,
			      ((struct primenetGetAssignment *)pkt)->k,
			      ((struct primenetGetAssignment *)pkt)->b,
			      ((struct primenetGetAssignment *)pkt)->n,
			      ((struct primenetGetAssignment *)pkt)->c);
		sprintf (buf, "Got assignment %s: %s\n",
			 ((struct primenetGetAssignment *)pkt)->assignment_uid,
			 info);
		LogMsg (buf);
		break;
	case PRIMENET_REGISTER_ASSIGNMENT:
		sprintf (buf, "Assignment registered as: %s\n",
			 ((struct primenetRegisterAssignment *)pkt)->assignment_uid);
		LogMsg (buf);
		break;
	}

/* Return the return code */

	return (return_code);
}


/* Get program options from the server. */

int getProgramOptions (void)
{
	struct primenetProgramOptions pkt;
	int	rc, tnum, mem_changed, restart, original_rob, mem_readable;
	unsigned int day_memory, night_memory, day_start_time, day_end_time;

/* Get the old-style memory settings */

	mem_readable = read_memory_settings (&day_memory, &night_memory,
					     &day_start_time, &day_end_time);

/* Loop once for global options (tnum = -1) and once for each worker */
/* thread (to get from the server the thread-specific options). */

	restart = FALSE;
	mem_changed = FALSE;
	original_rob = RUN_ON_BATTERY;
	for (tnum = -1; tnum < (int) NUM_WORKER_THREADS; tnum++) {

/* Init the packet.  A packet where no options are sent to the server tells */
/* the server to send us all the options saved on the server. */

		memset (&pkt, 0, sizeof (pkt));
		strcpy (pkt.computer_guid, COMPUTER_GUID);
		pkt.cpu_num = tnum;
		pkt.work_preference = -1;
		pkt.priority = -1;
		pkt.daysOfWork = -1;
		pkt.dayMemory = -1;
		pkt.nightMemory = -1;
		pkt.dayStartTime = -1;
		pkt.nightStartTime = -1;
		pkt.runOnBattery = -1;
		pkt.num_workers = -1;

/* Get the options from the server */

		LOCKED_WORK_UNIT = NULL;
		rc = sendMessage (PRIMENET_PROGRAM_OPTIONS, &pkt);
		if (rc) return (rc);

/* The server will send all the options it has stored in its database. */
/* Copy these to our global variables and INI files. */

		if (pkt.work_preference != -1) {
			if (tnum == -1) {
				PTOSetAll (INI_FILE, "WorkPreference", "SrvrPO1", WORK_PREFERENCE, pkt.work_preference);
			} else {
				PTOSetOne (INI_FILE, "WorkPreference", "SrvrPO1", WORK_PREFERENCE, tnum, pkt.work_preference);
			}
		}

/* These options cannot be set for each thread.  In theory, the server */
/* should only send these options when tnum is -1 (i.e. a global option). */
/* However, for robustness, we accept the option even if the server sends */
/* as only for this thread. */

		if (pkt.priority != -1) {
			PRIORITY = pkt.priority;
			IniWriteInt (INI_FILE, "Priority", PRIORITY);
			IniWriteInt (LOCALINI_FILE, "SrvrPO2", PRIORITY);
			restart = TRUE;
		}

		if (pkt.daysOfWork != -1) {
			DAYS_OF_WORK = pkt.daysOfWork;
			IniWriteInt (INI_FILE, "DaysOfWork", DAYS_OF_WORK);
			IniWriteInt (LOCALINI_FILE, "SrvrPO3", DAYS_OF_WORK);
		}

		if (pkt.dayMemory != -1) {
			if (day_memory != pkt.dayMemory) {
				day_memory = pkt.dayMemory;
				mem_changed = TRUE;
			}
			IniWriteInt (LOCALINI_FILE, "SrvrPO4", day_memory);
		}

		if (pkt.nightMemory != -1) {
			if (night_memory != pkt.nightMemory) {
				night_memory = pkt.nightMemory;
				mem_changed = TRUE;
			}
			IniWriteInt (LOCALINI_FILE, "SrvrPO5", night_memory);
		}

		if (pkt.dayStartTime != -1) {
			if (day_start_time != pkt.dayStartTime) {
				day_start_time = pkt.dayStartTime;
				mem_changed = TRUE;
			}
			IniWriteInt (LOCALINI_FILE, "SrvrPO6", day_start_time);
		}

		if (pkt.nightStartTime != -1) {
			if (day_end_time != pkt.nightStartTime) {
				day_end_time = pkt.nightStartTime;
				mem_changed = TRUE;
			}
			IniWriteInt (LOCALINI_FILE, "SrvrPO7", day_end_time);
		}

		if (pkt.runOnBattery != -1) {
			RUN_ON_BATTERY = pkt.runOnBattery;
			IniWriteInt (INI_FILE, "RunOnBattery", RUN_ON_BATTERY);
			IniWriteInt (LOCALINI_FILE, "SrvrPO8", RUN_ON_BATTERY);
		}

		if (pkt.num_workers != -1) {
			if (pkt.num_workers > (int) (NUM_CPUS * CPU_HYPERTHREADS))
				pkt.num_workers = NUM_CPUS * CPU_HYPERTHREADS;
			NUM_WORKER_THREADS = pkt.num_workers;
			IniWriteInt (LOCALINI_FILE, "WorkerThreads", NUM_WORKER_THREADS);
			IniWriteInt (LOCALINI_FILE, "SrvrPO9", NUM_WORKER_THREADS);
			restart = TRUE;
		}
		if (mem_readable && mem_changed)
			write_memory_settings (day_memory, night_memory, day_start_time, day_end_time);
	}

/* When we have finished getting the options for every CPU, then update */
/* the counter in the INI file that tracks the options counter.  The server */
/* increments the counter whenever the options are edited on the server. */
/* This updated count is sent in a uc pkt and compared to the value saved */
/* in the INI file.  If it is different we know we need to get the program */
/* options from the server. */

	IniWriteInt (LOCALINI_FILE, "SrvrP00", pkt.options_counter);

/* If memory settings, priority, num_workers, or run-on-battery changed, */
/* then restart threads that may be affected by the change. */

	if (mem_readable && mem_changed) mem_settings_have_changed ();
	if (restart) stop_workers_for_restart ();
	if (original_rob != RUN_ON_BATTERY) run_on_battery_changed ();
	return (0);
}


/* Send program options to the server if they've been changed locally. */

int sendProgramOptions (
	int	*talked_to_server)
{
	struct primenetProgramOptions pkt;
	int	rc, tnum, options_changed, work_pref_changed, mem_readable;
	unsigned int local_od;
	unsigned int day_memory, night_memory, day_start_time, day_end_time;

/* Get the old-style memory settings */

	mem_readable = read_memory_settings (&day_memory, &night_memory,
					     &day_start_time, &day_end_time);

/* Loop once for global options (tnum = -1) and once for each thread */
/* to send thread-specific options. */

	for (tnum = -1; tnum < (int) NUM_WORKER_THREADS; tnum++) {

		memset (&pkt, 0, sizeof (pkt));
		strcpy (pkt.computer_guid, COMPUTER_GUID);
		pkt.cpu_num = tnum;

		work_pref_changed = FALSE;
		options_changed = FALSE;

		pkt.work_preference = -1;
		if (PTOHasOptionChanged ("SrvrPO1", WORK_PREFERENCE, tnum)) {
			if (tnum == -1)
				pkt.work_preference = WORK_PREFERENCE[0];
			else
				pkt.work_preference = WORK_PREFERENCE[tnum];
			work_pref_changed = TRUE;
			options_changed = TRUE;
		}

		pkt.priority = -1;
		if (tnum == -1 &&
		    PRIORITY != IniGetInt (LOCALINI_FILE, "SrvrPO2", -1)) {
			pkt.priority = PRIORITY;
			options_changed = TRUE;
		}

		pkt.daysOfWork = -1;
		if (tnum == -1 &&
		    DAYS_OF_WORK != IniGetInt (LOCALINI_FILE, "SrvrPO3", -1)) {
			pkt.daysOfWork = DAYS_OF_WORK;
			options_changed = TRUE;
		}

		pkt.dayMemory = -1;
		if (tnum == -1 && mem_readable &&
		    day_memory != IniGetInt (LOCALINI_FILE, "SrvrPO4", -1)) {
			pkt.dayMemory = day_memory;
			options_changed = TRUE;
		}

		pkt.nightMemory = -1;
		if (tnum == -1 && mem_readable &&
		    night_memory != IniGetInt (LOCALINI_FILE, "SrvrPO5", -1)) {
			pkt.nightMemory = night_memory;
			options_changed = TRUE;
		}

		pkt.dayStartTime = -1;
		if (tnum == -1 && mem_readable &&
		    day_start_time != IniGetInt (LOCALINI_FILE, "SrvrPO6", -1)) {
			pkt.dayStartTime = day_start_time;
			options_changed = TRUE;
		}

		pkt.nightStartTime = -1;
		if (tnum == -1 && mem_readable &&
		    day_end_time != IniGetInt (LOCALINI_FILE, "SrvrPO7", -1)) {
			pkt.nightStartTime = day_end_time;
			options_changed = TRUE;
		}

		pkt.runOnBattery = -1;
		if (tnum == -1 &&
		    RUN_ON_BATTERY != IniGetInt (LOCALINI_FILE, "SrvrPO8", -1)) {
			pkt.runOnBattery = RUN_ON_BATTERY;
			options_changed = TRUE;
		}

		pkt.num_workers = -1;
		if (tnum == -1 &&
		    NUM_WORKER_THREADS != IniGetInt (LOCALINI_FILE, "SrvrPO9", -1)) {
			pkt.num_workers = NUM_WORKER_THREADS;
			options_changed = TRUE;
		}

/* If options haven't changed, we're done */

		if (!options_changed) continue;

/* Send the changed options */
	
		LOCKED_WORK_UNIT = NULL;
		rc = sendMessage (PRIMENET_PROGRAM_OPTIONS, &pkt);
		if (rc) return (rc);
		*talked_to_server = TRUE;

/* Now save the shadow copy of the changed options in our LOCALINI file */
/* so we can detect future changes to the program options (even if user */
/* hand edits prime.ini!) */

		if (work_pref_changed) {
			if (tnum == -1)
				PTOSetAll (INI_FILE, "WorkPreference", "SrvrPO1", WORK_PREFERENCE, WORK_PREFERENCE[0]);
			else
				PTOSetOne (INI_FILE, "WorkPreference", "SrvrPO1", WORK_PREFERENCE, tnum, WORK_PREFERENCE[tnum]);
        }
		if (tnum == -1) {
			IniWriteInt (LOCALINI_FILE, "SrvrPO2", PRIORITY);
			IniWriteInt (LOCALINI_FILE, "SrvrPO3", DAYS_OF_WORK);
			if (mem_readable) {
				IniWriteInt (LOCALINI_FILE, "SrvrPO4", day_memory);
				IniWriteInt (LOCALINI_FILE, "SrvrPO5", night_memory);
				IniWriteInt (LOCALINI_FILE, "SrvrPO6", day_start_time);
				IniWriteInt (LOCALINI_FILE, "SrvrPO7", day_end_time);
			}
			IniWriteInt (LOCALINI_FILE, "SrvrPO8", RUN_ON_BATTERY);
			IniWriteInt (LOCALINI_FILE, "SrvrPO9", NUM_WORKER_THREADS);
		}

/* Increment the options counter so that we can detect when the options */
/* are changed on the server.  The update we just did incremented the */
/* server's option counter.  We don't want to simply replace the options */
/* counter with the one in the pkt because of this scenario:  All synced at */
/* od=25, web page options change makes od=26, local change causes us to */
/* send new option above, the od value returned in the pkt is now 27. */
/* If we write 27 to our ini file then we will miss downloading the web */
/* change that caused the counter to become 26.  We blend the packet od */
/* value with the INI file value for maximum robustness (provides some */
/* protection from the INI file value somehow getting larger than the */
/* returned od value - causing us to miss future web page updates). */

		local_od = IniGetInt (LOCALINI_FILE, "SrvrP00", -1) + 1;
		if (local_od > pkt.options_counter)
			local_od = pkt.options_counter;
		IniWriteInt (LOCALINI_FILE, "SrvrP00", local_od);
	}
	return (0);
}


/* Send any queued up messages to the server.  See if we have enough */
/* work queued up.  If we have too much work, give some back */

#define RETRY_EXCEEDED "Retry count exceeded.\n"

void communicateWithServer (void *arg)
{
static	int	obsolete_client = FALSE;
static	int	send_message_retry_count = 0;
	unsigned long magicnum, version;
	unsigned long header_words[2];/* Flag words from spool file */
				/* We copy the header word to detect */
				/* any changes to the header word while */
				/* we are communicating with the server */
	int	fd;		/* Spool file handle */
	long	msg_offset;	/* File offset of current message */
	unsigned int tnum;
	double	est, work_to_get, unreserve_threshold;
	int	rc, stop_reason;
	int	talked_to_server = FALSE;
	int	server_options_counter, retry_count;
	char	buf[1000];

/* If we got an obsolete client error code earlier, then do not attempt */
/* any more communication with the server. */

	if (obsolete_client) goto leave;

/* Change title to show comm thread is active */

	ChangeIcon (COMM_THREAD_NUM, WORKING_ICON);
	title (COMM_THREAD_NUM, "Active");

/* There have been reports of computers losing their names. */
/* I can only see this happening if COMPID gets clobbered.  To combat this */
/* possibility, reread the computer name from the INI file. */

	if (COMPID[0] == 0) {
		IniGetString (LOCALINI_FILE, "ComputerID", COMPID,
			      sizeof (COMPID), NULL);
		sanitizeString (COMPID);
	}

/* This is the retry entry point.  Used when an error causes us to go back */
/* and start reprocessing from the header words.  Also used when we get all */
/* done and find that someone wrote to the spool file while this thread was */
/* running. */

	retry_count = 0;
retry:

/* Ping the server if so requested. */

	if (GET_PING_INFO == 1) {
		struct primenetPingServer pkt;
		memset (&pkt, 0, sizeof (pkt));
//bug		pkt.versionNumber = PRIMENET_VERSION;
		pkt.ping_type = 0;
		LOCKED_WORK_UNIT = NULL;
		rc = sendMessage (PRIMENET_PING_SERVER, &pkt);
		OutputStr (MAIN_THREAD_NUM, "\n");
		if (rc)
			OutputStr (MAIN_THREAD_NUM, "Failure pinging server");
		else
			OutputStr (MAIN_THREAD_NUM, pkt.ping_response);
		OutputStr (MAIN_THREAD_NUM, "\n");
		GET_PING_INFO = 0;
	}

/* Obtain the lock controlling spool file access.  Open the spool file. */

	gwmutex_lock (&SPOOL_FILE_MUTEX);
	fd = _open (SPOOL_FILE, _O_RDWR | _O_BINARY);
	if (fd < 0) goto locked_leave;

/* Read and validate the spool file header */

	if (!read_long (fd, &magicnum, NULL) ||
	    magicnum != SPOOL_FILE_MAGICNUM ||
	    !read_long (fd, &version, NULL) ||
	    version != SPOOL_FILE_VERSION ||
	    !read_long (fd, &header_words[0], NULL) ||
	    !read_long (fd, &header_words[1], NULL)) {
		_close (fd);
		gwmutex_unlock (&SPOOL_FILE_MUTEX);
		salvageCorruptSpoolFile ();
		goto retry;
	}

/* We must perform some complicated shenanigans to detect any changes to */
/* header word while this thread is running.  The header word is copied */
/* to the second word and the first word cleared.  That way any new calls */
/* to spoolMessage will set the first header word.  If we happen to */
/* crash then we must be careful not to lose the unprocessed bits in */
/* the second header word. */

	header_words[1] |= header_words[0];
	header_words[0] = 0;
	_lseek (fd, SPOOL_FILE_HEADER_OFFSET, SEEK_SET);
	write_long (fd, header_words[0], NULL);
	write_long (fd, header_words[1], NULL);

/* Close and unlock the spool file.  This allows worker threads to write */
/* messages to the spool file while we are sending messages. */

	_close (fd);
	gwmutex_unlock (&SPOOL_FILE_MUTEX);
	
/* Make sure we don't pummel the server with data.  Suppose user uses */
/* Advanced/Factor to find tiny factors again.  At least, make sure */
/* he only sends the data once every 5 minutes. */

//bug	next_comm_time = this_time + 300;  //bug - manual_comm overrides this
// move this test to SendMessage?  It could spawn add_timed_event
// rather than creating this thread.  However, does spreading out these
// messages help?  No.  It only helps if prime95 is buggy (like the spool file
// can't be deleted or modified so it just keeps being resent?)

//bug - check this global before each message is sent... in case ping is requested while
// a comm is in progress

/* Send computer info first */

	server_options_counter = -1;
	if (header_words[1] & HEADER_FLAG_UC) {
		struct primenetUpdateComputerInfo pkt;
		unsigned long server_uid, server_computer_name;

		memset (&pkt, 0, sizeof (pkt));
		strcpy (pkt.computer_guid, COMPUTER_GUID);
		strcpy (pkt.hardware_guid, HARDWARE_GUID);
		strcpy (pkt.windows_guid, WINDOWS_GUID);
		generate_application_string (pkt.application);
		strcpy (pkt.cpu_model, CPU_BRAND);
		pkt.cpu_features[0] = 0;
//		if (CPU_FLAGS & CPU_RDTSC) strcat (pkt.cpu_features, "RDTSC,");
//		if (CPU_FLAGS & CPU_CMOV) strcat (pkt.cpu_features, "CMOV,");
		if (CPU_FLAGS & CPU_PREFETCH) strcat (pkt.cpu_features, "Prefetch,");
		if (CPU_FLAGS & CPU_3DNOW) strcat (pkt.cpu_features, "3DNow!,");
//		if (CPU_FLAGS & CPU_MMX) strcat (pkt.cpu_features, "MMX,");
		if (CPU_FLAGS & CPU_SSE) strcat (pkt.cpu_features, "SSE,");
		if (CPU_FLAGS & CPU_SSE2) strcat (pkt.cpu_features, "SSE2,");
		if (CPU_FLAGS & CPU_SSE41) strcat (pkt.cpu_features, "SSE4,");
		if (CPU_FLAGS & CPU_AVX) strcat (pkt.cpu_features, "AVX,");
		if (CPU_FLAGS & CPU_AVX2) strcat (pkt.cpu_features, "AVX2,");
		if (CPU_FLAGS & (CPU_FMA3 | CPU_FMA4)) strcat (pkt.cpu_features, "FMA, ");
		if (pkt.cpu_features[0])
			pkt.cpu_features[strlen (pkt.cpu_features) - 1] = 0;
		pkt.L1_cache_size = CPU_L1_CACHE_SIZE;
		pkt.L2_cache_size = CPU_L2_CACHE_SIZE;
		pkt.L3_cache_size = CPU_L3_CACHE_SIZE;
		pkt.num_cpus = NUM_CPUS;
		pkt.num_hyperthread = CPU_HYPERTHREADS;
		pkt.mem_installed = physical_memory ();
		pkt.cpu_speed = (int) CPU_SPEED;
		pkt.hours_per_day = CPU_HOURS;
		pkt.rolling_average = ROLLING_AVERAGE;
		/* The spec says only send the UserID if it has changed. */
		/* Rather than store the UID that we last sent to the server */
		/* we store a hash value so the user is not tempted to hand edit it. */
		server_uid = IniGetInt (LOCALINI_FILE, "SrvrUID", 0);
		if (string_to_hash (USERID) != server_uid)
			strcpy (pkt.user_id, USERID);
		/* The spec says only send the computer name if it changed */
		/* Again store a hash value so the user won't hand edit the value */
		server_computer_name = IniGetInt (LOCALINI_FILE, "SrvrComputerName", 0);
		if (string_to_hash (COMPID) != server_computer_name)
			strcpy (pkt.computer_name, COMPID);
		LOCKED_WORK_UNIT = NULL;
		rc = sendMessage (PRIMENET_UPDATE_COMPUTER_INFO, &pkt);
		if (rc) goto error_exit;
		talked_to_server = TRUE;
		strcpy (USERID, pkt.user_id);
		IniWriteString (INI_FILE, "V5UserID", USERID);
		IniWriteInt (LOCALINI_FILE, "SrvrUID", string_to_hash (USERID));
		strcpy (COMPID, pkt.computer_name);
		IniWriteString (LOCALINI_FILE, "ComputerID", COMPID);
		IniWriteInt (LOCALINI_FILE, "SrvrComputerName", string_to_hash (COMPID));
		server_options_counter = pkt.options_counter;
	}

/* When we first register a computer, the server may have initialized */
/* some default program options for us to download.  Get them before */
/* we overwrite them with our default option values. */

//bug - who wins when user has also set some options locally during the
//startup_in_progress code path?  What about options that came from a
// copied prime.ini file - which take precedence?
	if (server_options_counter == 1 &&
	    server_options_counter > IniGetInt (LOCALINI_FILE, "SrvrP00", 0)) {
		rc = getProgramOptions ();
		if (rc) goto error_exit;
		talked_to_server = TRUE;
	}

/* Send our changed program options.  Do this before they could get */
/* overwritten by downloading changed options on the server (since the */
/* server sends all options, not just the changed options). */

	rc = sendProgramOptions (&talked_to_server);
	if (rc) goto error_exit;

/* Finally get program options if they've been changed on the server. */

	if (server_options_counter > IniGetInt (LOCALINI_FILE, "SrvrP00", 0)) {
		rc = getProgramOptions ();
		if (rc) goto error_exit;
		talked_to_server = TRUE;
	}

/* Send the messages */

	msg_offset = SPOOL_FILE_MSG_OFFSET;
	for ( ; ; ) {
		short	msgType;
		union {
			struct primenetAssignmentResult ar;
			struct primenetAssignmentUnreserve au;
			struct primenetBenchmarkData bd;
		} msg;
		long	new_offset;

/* Read a message */

		gwmutex_lock (&SPOOL_FILE_MUTEX);
		fd = _open (SPOOL_FILE, _O_RDONLY | _O_BINARY);
		if (fd < 0) goto locked_leave;
		_lseek (fd, msg_offset, SEEK_SET);
		readMessage (fd, &msg_offset, &msgType, &msg);
		new_offset = _lseek (fd, 0, SEEK_CUR);
		_close (fd);
		gwmutex_unlock (&SPOOL_FILE_MUTEX);

/* If there was a message, send it now.  Ignore most errors.  We do this */
/* so that a corrupt spool file will not "get stuck" trying to send the */
/* same corrupt message over and over again. */

		if (msgType == 0) break;
		for ( ; ; ) {
			LOCKED_WORK_UNIT = NULL;
			rc = sendMessage (msgType, &msg);
			/* If the computer ID is bad and not ours (a damaged */
			/* spool file or we were forced to generate a new */
			/* computer GUID?) then try again using our computer ID */
			if ((rc == PRIMENET_ERROR_UNREGISTERED_CPU ||
			     rc == PRIMENET_ERROR_CPU_CONFIGURATION_MISMATCH ||
			     rc == PRIMENET_ERROR_STALE_CPU_INFO) &&
			    strcmp (msg.ar.computer_guid, COMPUTER_GUID)) {
				OutputStr (COMM_THREAD_NUM, "Retrying message with this computer's GUID.\n");
				strcpy (msg.ar.computer_guid, COMPUTER_GUID);
				continue;
			}
			break;
		}
		if (rc >= PRIMENET_FIRST_INTERNAL_ERROR ||
		    rc == PRIMENET_ERROR_SERVER_BUSY ||
		    rc == PRIMENET_ERROR_OBSOLETE_CLIENT ||
		    rc == PRIMENET_ERROR_UNREGISTERED_CPU ||
		    rc == PRIMENET_ERROR_CPU_CONFIGURATION_MISMATCH ||
		    rc == PRIMENET_ERROR_STALE_CPU_INFO)
			goto error_exit;
		talked_to_server = TRUE;

/* Handle errors that show the server processed the message properly. */
/* These errors can happen with unwanted results (retesting known primes, */
/* PRP results, ECM or P-1 on non-Mersennes, etc.) */ 

		if (rc == PRIMENET_ERROR_INVALID_ASSIGNMENT_KEY ||
		    rc == PRIMENET_ERROR_INVALID_RESULT_TYPE ||
		    rc == PRIMENET_ERROR_NO_ASSIGNMENT ||
		    rc == PRIMENET_ERROR_WORK_NO_LONGER_NEEDED ||
		    rc == PRIMENET_ERROR_INVALID_PARAMETER)
			rc = 0;

/* Just in case the error was casued by some kind of unexpected */
/* and unrepeatable server problem, retry sending the message 5 times */
/* before giving up and moving on to the next one. */

		if (rc) {
			if (send_message_retry_count++ < 5) {
				rc = PRIMENET_ERROR_SERVER_BUSY;
				goto error_exit;
			}
			LogMsg ("Deleting unprocessed message from spool file.\n");
		}
		send_message_retry_count = 0;

/* Flag the message as successfully sent.  Even if there was an error, the */
/* error code is such that resending the message will not be helpful. */

		gwmutex_lock (&SPOOL_FILE_MUTEX);
		fd = _open (SPOOL_FILE, _O_RDWR | _O_BINARY);
		if (fd < 0) goto locked_leave;
		_lseek (fd, msg_offset, SEEK_SET);
		msgType = -1;
		(void) _write (fd, &msgType, sizeof (short));
		_close (fd);
		gwmutex_unlock (&SPOOL_FILE_MUTEX);
		msg_offset = new_offset;
	}

/* Loop over all worker threads to get enough work for each thread */

	unreserve_threshold = IniGetInt (INI_FILE, "UnreserveDays", 30) * 86400.0;
	work_to_get = DAYS_OF_WORK * 86400.0;
	for (tnum = 0; tnum < NUM_WORKER_THREADS; tnum++) {
	    struct work_unit *w;
	    int	num_work_units;

/* Get work to do until we've accumulated enough to keep us busy for */
/* a while.  If we have too much work to do, lets give some back. */

	    num_work_units = 0;
	    est = 0.0;
	    for (w = NULL; ; ) {
		int	registered_assignment;

/* Read the line of the work file */

		w = getNextWorkToDoLine (tnum, w, SHORT_TERM_USE);
		if (w == NULL) break;
		if (w->work_type == WORK_NONE) continue;

/* If we are quitting GIMPS or we have too much work queued up, */
/* then return the assignment. */

		if (header_words[1] & HEADER_FLAG_QUIT_GIMPS ||
		    (est >= work_to_get + unreserve_threshold &&
		     ! isWorkUnitActive (w) &&
		     w->pct_complete == 0.0)) {

			if (w->assignment_uid[0]) {
				struct primenetAssignmentUnreserve pkt;
				memset (&pkt, 0, sizeof (pkt));
				strcpy (pkt.computer_guid, COMPUTER_GUID);
				strcpy (pkt.assignment_uid, w->assignment_uid);
				LOCKED_WORK_UNIT = w;
				rc = sendMessage (PRIMENET_ASSIGNMENT_UNRESERVE, &pkt);
				if (rc && rc != PRIMENET_ERROR_INVALID_ASSIGNMENT_KEY) {
					decrementWorkUnitUseCount (w, SHORT_TERM_USE);
					goto error_exit;
				}
				talked_to_server = TRUE;
			}
			stop_reason = deleteWorkToDoLine (tnum, w, TRUE);
			if (stop_reason) goto error_exit;
			continue;
		}

/* Adjust our time estimate */

		num_work_units++;
		est += work_estimate (tnum, w);

/* Register assignments that were not issued by the server */

		registered_assignment = FALSE;
		if (!w->assignment_uid[0] && !w->ra_failed) {
			struct primenetRegisterAssignment pkt;
			memset (&pkt, 0, sizeof (pkt));
			strcpy (pkt.computer_guid, COMPUTER_GUID);
			pkt.cpu_num = tnum;
			if (w->work_type == WORK_FACTOR)
				pkt.work_type = PRIMENET_WORK_TYPE_FACTOR;
			if (w->work_type == WORK_PMINUS1)
				pkt.work_type = PRIMENET_WORK_TYPE_PMINUS1;
			if (w->work_type == WORK_PFACTOR)
				pkt.work_type = PRIMENET_WORK_TYPE_PFACTOR;
			if (w->work_type == WORK_ECM)
				pkt.work_type = PRIMENET_WORK_TYPE_ECM;
			if (w->work_type == WORK_TEST ||
			    w->work_type == WORK_ADVANCEDTEST)
				pkt.work_type = PRIMENET_WORK_TYPE_FIRST_LL;
			if (w->work_type == WORK_DBLCHK)
				pkt.work_type = PRIMENET_WORK_TYPE_DBLCHK;
			if (w->work_type == WORK_PRP)
				pkt.work_type = PRIMENET_WORK_TYPE_PRP;
			pkt.k = w->k;
			pkt.b = w->b;
			pkt.n = w->n;
			pkt.c = w->c;
			pkt.how_far_factored = w->sieve_depth;
			pkt.factor_to = w->factor_to;
			pkt.has_been_pminus1ed = w->pminus1ed;
			pkt.B1 = w->B1;
			pkt.B2 = w->B2;
			pkt.curves = w->curves_to_do;
			pkt.tests_saved = w->tests_saved;
			LOCKED_WORK_UNIT = w;
			rc = sendMessage (PRIMENET_REGISTER_ASSIGNMENT, &pkt);
			if (rc &&
			    rc != PRIMENET_ERROR_NO_ASSIGNMENT &&
			    rc != PRIMENET_ERROR_INVALID_ASSIGNMENT_TYPE &&
			    rc != PRIMENET_ERROR_INVALID_PARAMETER) {
				decrementWorkUnitUseCount (w, SHORT_TERM_USE);
				goto error_exit;
			}
			talked_to_server = TRUE;
			if (rc)
				w->ra_failed = TRUE;
			else {
				strcpy (w->assignment_uid, pkt.assignment_uid);
				registered_assignment = TRUE;
			}
			updateWorkToDoLine (tnum, w);
		}

/* Update the server on the work unit's projected completion date */
/* If we get an invalid assignment key, then the user probably unreserved */
/* the exponent using the web forms - delete it from our work to do file. */

		if ((header_words[1] & HEADER_FLAG_END_DATES || registered_assignment) &&
		    w->assignment_uid[0]) {
			struct primenetAssignmentProgress pkt2;
			memset (&pkt2, 0, sizeof (pkt2));
			strcpy (pkt2.computer_guid, COMPUTER_GUID);
			pkt2.cpu_num = tnum;
			strcpy (pkt2.assignment_uid, w->assignment_uid);
			strcpy (pkt2.stage, w->stage);
			pkt2.pct_complete = w->pct_complete * 100.0;
			pkt2.end_date = (unsigned long) est;
			pkt2.next_update = (uint32_t) (DAYS_BETWEEN_CHECKINS * 86400.0);
			pkt2.fftlen = w->fftlen;
			LOCKED_WORK_UNIT = w;
			rc = sendMessage (PRIMENET_ASSIGNMENT_PROGRESS, &pkt2);
			if (rc == PRIMENET_ERROR_INVALID_ASSIGNMENT_KEY ||
			    rc == PRIMENET_ERROR_WORK_NO_LONGER_NEEDED) {
				est -= work_estimate (tnum, w);
				rc = deleteWorkToDoLine (tnum, w, TRUE);
			} else if (rc) {
				decrementWorkUnitUseCount (w, SHORT_TERM_USE);
				goto error_exit;
			}
			talked_to_server = TRUE;
		}
	    }

/* If we are quitting gimps, do not get more work.  Note that there is a */
/* race condition when quitting gimps that makes it possible to get here */
/* with a request to get exponents and USE_PRIMENET not set. */

	    if (header_words[1] & HEADER_FLAG_QUIT_GIMPS ||
		IniGetInt (INI_FILE, "NoMoreWork", 0) ||
		!USE_PRIMENET)
		    continue;

/* If we don't have enough work to do, get more work from the server. */

	    while (est < work_to_get &&
		   ! IniGetInt (INI_FILE, "NoMoreWork", 0) &&
		   num_work_units < IniGetInt (INI_FILE, "MaxExponents", 15)) {
		struct primenetGetAssignment pkt1;
		struct primenetAssignmentProgress pkt2;
		struct work_unit w;

/* Get a work unit to process */

		memset (&pkt1, 0, sizeof (pkt1));
		strcpy (pkt1.computer_guid, COMPUTER_GUID);
		pkt1.cpu_num = tnum;
		LOCKED_WORK_UNIT = NULL;
		rc = sendMessage (PRIMENET_GET_ASSIGNMENT, &pkt1);
		if (rc) goto error_exit;
		talked_to_server = TRUE;

/* Sanity check that the server hasn't sent us bogus (too small) exponents to test */

		if (pkt1.n < 15000000 &&
		    (pkt1.work_type == PRIMENET_WORK_TYPE_FACTOR ||
		     pkt1.work_type == PRIMENET_WORK_TYPE_PFACTOR ||
		     pkt1.work_type == PRIMENET_WORK_TYPE_FIRST_LL ||
		     pkt1.work_type == PRIMENET_WORK_TYPE_DBLCHK)) {
			sprintf (buf, "Server sent bad exponent: %lu.\n", (unsigned long) pkt1.n);
			LogMsg (buf);
			goto error_exit;
		}

/* Format the work_unit structure based on the work_type */

		memset (&w, 0, sizeof (w));
		strcpy (w.assignment_uid, pkt1.assignment_uid);
		w.k = 1.0;	/* Set some default values */
		w.b = 2;
		w.n = pkt1.n;
		w.c = -1;
		switch (pkt1.work_type) {
		case PRIMENET_WORK_TYPE_FACTOR:
			w.work_type = WORK_FACTOR;
			w.sieve_depth = pkt1.how_far_factored;
			w.factor_to = pkt1.factor_to;
			break;
		case PRIMENET_WORK_TYPE_PFACTOR:
			w.work_type = WORK_PFACTOR;
			w.sieve_depth = pkt1.how_far_factored;
			w.pminus1ed = pkt1.has_been_pminus1ed;
			w.tests_saved = pkt1.tests_saved;
			break;
		case PRIMENET_WORK_TYPE_FIRST_LL:
			w.work_type = WORK_TEST;
			w.sieve_depth = pkt1.how_far_factored;
			w.pminus1ed = pkt1.has_been_pminus1ed;
			break;
		case PRIMENET_WORK_TYPE_DBLCHK:
			w.work_type = WORK_DBLCHK;
			w.sieve_depth = pkt1.how_far_factored;
			w.pminus1ed = pkt1.has_been_pminus1ed;
			break;
		case PRIMENET_WORK_TYPE_PMINUS1:
			w.work_type = WORK_PMINUS1;
			w.k = pkt1.k;
			w.b = pkt1.b;
			w.c = pkt1.c;
			w.B1 = pkt1.B1;
			w.B2 = pkt1.B2;
			break;
		case PRIMENET_WORK_TYPE_ECM:
			w.work_type = WORK_ECM;
			w.k = pkt1.k;
			w.b = pkt1.b;
			w.c = pkt1.c;
			w.B1 = pkt1.B1;
			w.B2 = pkt1.B2;
			w.curves_to_do = pkt1.curves;
			break;
		case PRIMENET_WORK_TYPE_PRP:
			w.work_type = WORK_PRP;
			w.k = pkt1.k;
			w.b = pkt1.b;
			w.c = pkt1.c;
			w.sieve_depth = pkt1.how_far_factored;
			w.tests_saved = pkt1.tests_saved;
			break;
		default:
			sprintf (buf, "Received unknown work type: %lu.\n",
				 (unsigned long) pkt1.work_type);
			LogMsg (buf);
			goto error_exit;
		}
		if (pkt1.known_factors[0]) { /* ECM, P-1, PRP may have this */
			w.known_factors = (char *)
				malloc (strlen (pkt1.known_factors) + 1);
			if (w.known_factors == NULL) {
				LogMsg ("Memory allocation error\n");
				goto error_exit;
			}
			strcpy (w.known_factors, pkt1.known_factors);
		}

/* Write the exponent to our worktodo file, before acknowledging the */
/* assignment with a projected completion date. */

		stop_reason = addWorkToDoLine (tnum, &w);
		if (stop_reason) goto error_exit;

//bug Can we somehow verify that the worktodo.ini line got written???
//bug (The rogue 'cat' problem)  If not, turn off USE_PRIMENET.
//bug Or rate limit get assignments to N per day (where N takes into account unreserves?)
//bug or have uc return the number of active assignments and ga the capability
//bug to retrieve them

/* Add work unit to our time estimate */

		num_work_units++;
		est = est + work_estimate (tnum, &w);

/* Acknowledge the exponent by sending a projected completion date */

		memset (&pkt2, 0, sizeof (pkt2));
		strcpy (pkt2.computer_guid, COMPUTER_GUID);
		pkt2.cpu_num = tnum;
		strcpy (pkt2.assignment_uid, pkt1.assignment_uid);
		pkt2.pct_complete = 0.0;
		pkt2.end_date = (unsigned long) est;
		pkt2.next_update = (uint32_t) (DAYS_BETWEEN_CHECKINS * 86400.0);
		LOCKED_WORK_UNIT = NULL;
		rc = sendMessage (PRIMENET_ASSIGNMENT_PROGRESS, &pkt2);
		if (rc == PRIMENET_ERROR_INVALID_ASSIGNMENT_KEY) {
			w.assignment_uid[0] = 0;
			updateWorkToDoLine (tnum, &w);
		} else if (rc) {
			goto error_exit;
		}
		talked_to_server = TRUE;
	    }
	}

/* Set some ini flags after we've successfully quit gimps */

	if ((header_words[1] & HEADER_FLAG_QUIT_GIMPS && USE_PRIMENET) ||
	    (IniGetInt (INI_FILE, "NoMoreWork", 0) && WORKTODO_COUNT == 0)) {
		USE_PRIMENET = 0;
		IniWriteInt (INI_FILE, "UsePrimenet", 0);
		IniWriteInt (INI_FILE, "NoMoreWork", 0);
		OutputSomewhere (COMM_THREAD_NUM, "Successfully quit GIMPS.\n");
	}

/* After sending new completion dates remember the current time */
/* so that we can send new completion dates in a month.  Set a timer */
/* to send them again. */

	else if (header_words[1] & HEADER_FLAG_END_DATES) {
		time_t current_time;
		time (&current_time);
		IniWriteInt (LOCALINI_FILE, "LastEndDatesSent", (long) current_time);
		if (!MANUAL_COMM)
			add_timed_event (TE_COMPLETION_DATES,
					 (int) (DAYS_BETWEEN_CHECKINS * 86400.0));
	}

/* Delete the spool file. However, we don't delete the file if any writes */
/* took place after we read the header words.  We detect writes by examining */
/* the first header word. */

	gwmutex_lock (&SPOOL_FILE_MUTEX);
	fd = _open (SPOOL_FILE, _O_RDWR | _O_BINARY);
	if (fd < 0) goto locked_leave;
	_lseek (fd, SPOOL_FILE_HEADER_OFFSET, SEEK_SET);
	read_long (fd, &header_words[0], NULL);
	read_long (fd, &header_words[1], NULL);
	if (header_words[0]) {
		if (header_words[1] & HEADER_FLAG_QUIT_GIMPS)
			header_words[0] |= HEADER_FLAG_QUIT_GIMPS;
		header_words[1] = 0;
		_lseek (fd, SPOOL_FILE_HEADER_OFFSET, SEEK_SET);
		write_long (fd, header_words[0], NULL);
		write_long (fd, header_words[1], NULL);
		_close (fd);
		gwmutex_unlock (&SPOOL_FILE_MUTEX);
		goto retry;
	}
	_close (fd);
	_unlink (SPOOL_FILE);

/* Tell user we're done communicating, then exit this communication thread */

	if (talked_to_server)
		OutputStr (COMM_THREAD_NUM, "Done communicating with server.\n");
	goto locked_leave;

/* We got an error communicating with the server.  Check for error codes */
/* that require special handling. */

error_exit:

/* If an UPDATE_COMPUTER_INFO command gets an invalid user error, then */
/* the user entered a non-existant userid.  Reset the userid field to the */
/* last value sent to or sent by the server.  This will cause the next uc */
/* command to get the current userid value from the server. */

	if (rc == PRIMENET_ERROR_INVALID_USER) {
		IniGetString (LOCALINI_FILE, "SrvrUID",
			      USERID, sizeof (USERID), NULL);
		spoolMessage (PRIMENET_UPDATE_COMPUTER_INFO, NULL);
		if (retry_count++ < 10) goto retry;
		OutputStr (COMM_THREAD_NUM, RETRY_EXCEEDED);
	}

/* If an UPDATE_COMPUTER_INFO command gets an cpu identity mismatch error */
/* then the previously registered hardware and windows hashes for this */
/* computer id don't match.  Either the machine hardware has been upgraded */
/* or Windows reinstalled OR the user moved the entire directory to another */
/* machine (a very common way that users install on multiple machines). */
/* First, as a safety measure, reget the computer UID from the INI file in */
/* case it became corrupted in memory.  If no corruption detected then */
/* generate a new computer uid and do an update computer info command again. */

	if (rc == PRIMENET_ERROR_CPU_IDENTITY_MISMATCH) {
		generate_computer_guid ();
//bug - ask user before changing uid?  Delete worktodo? Manipulate computer
// name?
//bug - log a useful message describing error and our remedial action?
		spoolMessage (PRIMENET_UPDATE_COMPUTER_INFO, NULL);
		if (retry_count++ < 10) goto retry;
		OutputStr (COMM_THREAD_NUM, RETRY_EXCEEDED);
	}

/* If the error is STALE_CPU_INFO, then the server is simply */
/* requesting us to do an update computer info again.  Spool the message */
/* then loop to send that message. */

	if (rc == PRIMENET_ERROR_STALE_CPU_INFO) {
		spoolMessage (PRIMENET_UPDATE_COMPUTER_INFO, NULL);
		if (retry_count++ < 10) goto retry;
		OutputStr (COMM_THREAD_NUM, RETRY_EXCEEDED);
	}

/* If the error is UNREGISTERED_CPU, then the COMPUTER_GUID */
/* value is corrupt.  Reget it from the INI file just in case the value */
/* became corrupt in memory.  Otherwise, generate and register a new computer uid */

	if (rc == PRIMENET_ERROR_UNREGISTERED_CPU) {
		// bug
//bug - log a useful message describing error and our remedial action?
		generate_computer_guid ();
		spoolMessage (PRIMENET_UPDATE_COMPUTER_INFO, NULL);
		if (retry_count++ < 10) goto retry;
		OutputStr (COMM_THREAD_NUM, RETRY_EXCEEDED);
	}

/* If the error is CPU_CONFIGURATION_MISMATCH, then the server has */
/* detected an inconsistency in the program options.  Resend them all. */

	if (rc == PRIMENET_ERROR_CPU_CONFIGURATION_MISMATCH) {
		clearCachedProgramOptions ();
		spoolMessage (PRIMENET_UPDATE_COMPUTER_INFO, NULL);
		if (retry_count++ < 10) goto retry;
		OutputStr (COMM_THREAD_NUM, RETRY_EXCEEDED);
	}

/* If the error indicates this client is too old, then turn off the */
/* "use primenet" setting until the user upgrades. */

	if (rc == PRIMENET_ERROR_OBSOLETE_CLIENT) {
//bug - delete spool file?NO  write a messsage to screen/log/etc???
//bug - log a useful message describing error and our remedial action?
		OutputStr (MAIN_THREAD_NUM, "Client is obsolete.  Please upgrage.\n");
		obsolete_client = TRUE;
	}

/* If an invalid work preference error is generated (can happen if user */
/* hand edits the prime.ini file) then set the work preference to zero */
/* (get work that makes the most sense) and resend the program options. */
/* As usual, reget the work preference from the ini file just in case */
/* the in-memory copy was inexplicably corrupted. */

	if (rc == PRIMENET_ERROR_INVALID_WORK_TYPE) {
//bug		short	ini_work_preference;
//bug - log a useful message describing error and our remedial action?
//bug		ini_work_preference = (short)
//bug			IniGetInt (INI_FILE, "WorkPreference", 0);
//bug		if (WORK_PREFERENCE != ini_work_preference)
//bug			WORK_PREFERENCE = ini_work_preference;
//bug		else
//bug			WORK_PREFERENCE = 0;
//bug		IniWriteInt (INI_FILE, "WorkPreference", WORK_PREFERENCE);
		spoolMessage (PRIMENET_PROGRAM_OPTIONS, NULL);
		if (retry_count++ < 10) goto retry;
		OutputStr (COMM_THREAD_NUM, RETRY_EXCEEDED);
	}

/* Otherwise, print out message saying where more help may be found. */

	OutputStr (COMM_THREAD_NUM, "Visit http://mersenneforum.org for help.\n");

/* We could not contact the server.  Set up a timer to relaunch the */
/* communication thread. */

	if (!MANUAL_COMM) {
		unsigned int retry_time;
		retry_time = (rc == PRIMENET_ERROR_MODEM_OFF) ?
			MODEM_RETRY_TIME : NETWORK_RETRY_TIME;
		sprintf (buf, "Will try contacting server again in %d %s.\n",
			 retry_time, retry_time == 1 ? "minute" : "minutes");
		OutputStr (COMM_THREAD_NUM, buf);
		ChangeIcon (COMM_THREAD_NUM, IDLE_ICON);
		sprintf (buf, "Waiting %d %s",
			 retry_time, retry_time == 1 ? "minute" : "minutes");
		title (COMM_THREAD_NUM, buf);
		add_timed_event (TE_COMM_SERVER, retry_time * 60);
	}

/* Clear handle indicating communication thread is active.  Return.  This */
/* will end the communication thread. */

leave:
	gwmutex_lock (&SPOOL_FILE_MUTEX);
locked_leave:
	COMMUNICATION_THREAD = 0;
	gwmutex_unlock (&SPOOL_FILE_MUTEX);
	if (!is_timed_event_active (TE_COMM_SERVER)) {
		ChangeIcon (COMM_THREAD_NUM, IDLE_ICON);
		title (COMM_THREAD_NUM, "Inactive");
	}
}

/* This routine tries to salvage information from a corrupt spool file */

void salvageCorruptSpoolFile (void)
{
	int	fd;
	char	filename[128];
	char	inbuf[1000];
	struct primenetAssignmentResult pkt;
	char	*in, *out;
	int	i, j, len, pkts_sent;

/* Output a message */

	OutputBoth (MAIN_THREAD_NUM,
		    "Spool file is corrupt.  Attempting to salvage data.\n");

/* Rename corrupt spool file before extracting data from it */

	strcpy (filename, SPOOL_FILE);
	filename[strlen(filename)-1]++;
	gwmutex_lock (&SPOOL_FILE_MUTEX);
	_unlink (filename);
	rename (SPOOL_FILE, filename);
	gwmutex_unlock (&SPOOL_FILE_MUTEX);

/* Read corrupt file.  Find ASCII characters and send them to the server */
/* in primenetAssignmentResult packets.  Process up to 100000 bytes of */
/* the spool file or 100 lines of output. */

	fd = _open (filename, _O_RDONLY | _O_BINARY);
	if (fd >= 0) {
		memset (&pkt, 0, sizeof (pkt));
		strcpy (pkt.computer_guid, COMPUTER_GUID);
		pkt.result_type = PRIMENET_AR_NO_RESULT;
		pkts_sent = 0;
		strcpy (pkt.message, "RECOVERY: ");
		out = pkt.message + 10;
		for (i = 0; i < 100 && pkts_sent < 100; i++) {
			len = _read (fd, &inbuf, sizeof (inbuf));
			if (len <= 0) break;
			for (j = 0, in = inbuf; j < len; j++, in++) {
				if (! isprint (*in) && *in != '\n') continue;
				*out++ = *in;
				if (*in == '\n');
				else if ((int) (out - pkt.message) < sizeof (pkt.message) - 1) continue;
				else *out++ = '\n';
				*out++ = 0;
				if ((int) (out - pkt.message) > 20) {
					spoolMessage (PRIMENET_ASSIGNMENT_RESULT, &pkt);
					pkts_sent++;
				}
				out = pkt.message + 10;
			}
		}
		_close (fd);
		*out++ = '\n';
		*out++ = 0;
		if ((int) (out - pkt.message) > 20)
			spoolMessage (PRIMENET_ASSIGNMENT_RESULT, &pkt);
	}

/* Delete the corrupt spool file */

	_unlink (filename);
}

/****************************************************************************/
/*                        Timed Events Handler                              */
/****************************************************************************/

gwmutex	TIMED_EVENTS_MUTEX;		/* Lock for timed events data */
gwevent TIMED_EVENTS_CHANGED;		/* Signal for telling the scheuler */
					/* there has been a change in the */
					/* array of timed events */
gwthread TIMED_EVENTS_THREAD = 0;	/* Thread for timed events */

struct {
	int	active;			/* TRUE if event is active */
	time_t	time_to_fire;		/* When to start this event */
} timed_events[MAX_TIMED_EVENTS];	/* Array of active timed events */

void timed_events_scheduler (void *arg);


void init_timed_event_handler (void)
{
	int	i;
	gwmutex_init (&TIMED_EVENTS_MUTEX);
	gwevent_init (&TIMED_EVENTS_CHANGED);
	TIMED_EVENTS_THREAD = 0;
	for (i = 0; i < MAX_TIMED_EVENTS; i++)
		timed_events[i].active = FALSE;
}

void add_timed_event (
	int	event_number,		/* Which event to add */
	int	time_to_fire)		/* When to start event (seconds from now) */
{
	time_t	this_time;

	gwmutex_lock (&TIMED_EVENTS_MUTEX);

	time (&this_time);
	timed_events[event_number].active = TRUE;
	timed_events[event_number].time_to_fire = this_time + time_to_fire;

	if (!TIMED_EVENTS_THREAD)
		gwthread_create (&TIMED_EVENTS_THREAD, &timed_events_scheduler, NULL);
	else
		gwevent_signal (&TIMED_EVENTS_CHANGED);

	gwmutex_unlock (&TIMED_EVENTS_MUTEX);
}

void delete_timed_event (
	int	event_number)		/* Which event to delete */
{
	gwmutex_lock (&TIMED_EVENTS_MUTEX);
	if (timed_events[event_number].active) {
		timed_events[event_number].active = FALSE;
		gwevent_signal (&TIMED_EVENTS_CHANGED);
	}
	gwmutex_unlock (&TIMED_EVENTS_MUTEX);
}

int is_timed_event_active (
	int	event_number)		/* Which event to test */
{
	return (timed_events[event_number].active);
}

time_t timed_event_fire_time (
	int	event_number)		/* Which event to get fire time of */
{
	return (timed_events[event_number].time_to_fire);
}

void timed_events_scheduler (void *arg)
{
	time_t	this_time;
	int	i;
	time_t	wake_up_time;
	int	there_are_active_events;

/* Loop forever, sleeping until next event fires up */

	for ( ; ; ) {

/* Determine how long until the next timed event */

		gwmutex_lock (&TIMED_EVENTS_MUTEX);
		there_are_active_events = FALSE;
		for (i = 0; i < MAX_TIMED_EVENTS; i++) {
			if (!timed_events[i].active) continue;
			if (!there_are_active_events ||
			    wake_up_time > timed_events[i].time_to_fire)
				wake_up_time = timed_events[i].time_to_fire;
			there_are_active_events = TRUE;
		}
		gwevent_reset (&TIMED_EVENTS_CHANGED);
		gwmutex_unlock (&TIMED_EVENTS_MUTEX);

/* Sleep until we have work to do.  If there are no active events, then */
/* sleep for a million seconds. */

		time (&this_time);
		if (!there_are_active_events) wake_up_time = this_time + 1000000;
		if (wake_up_time > this_time) {
			int	rc;
			rc = gwevent_wait (&TIMED_EVENTS_CHANGED,
					   (int) (wake_up_time - this_time));
			if (rc != GWEVENT_TIMED_OUT) continue;
		}

/* Do action associated with any timed events that have triggered */

		time (&this_time);
		for (i = 0; i < MAX_TIMED_EVENTS; i++) {
			int	fire;
			gwmutex_lock (&TIMED_EVENTS_MUTEX);
			fire = (timed_events[i].active && this_time >= timed_events[i].time_to_fire);
			gwmutex_unlock (&TIMED_EVENTS_MUTEX);
			if (!fire) continue;
			switch (i) {
			case TE_MEM_CHANGE:	/* Night/day memory change event */
				timed_events[i].active = FALSE;
				mem_settings_have_changed ();
				break;
			case TE_PAUSE_WHILE:	/* Check pause_while_running event */
				timed_events[i].active = FALSE;
				checkPauseWhileRunning ();
				break;
			case TE_WORK_QUEUE_CHECK:	/* Check work queue event */
				timed_events[i].time_to_fire = this_time + TE_WORK_QUEUE_CHECK_FREQ;
				spoolMessage (MSG_CHECK_WORK_QUEUE, NULL);
				break;
			case TE_COMM_SERVER:	/* Retry communication with server event */
				timed_events[i].active = FALSE;
				gwthread_create (&COMMUNICATION_THREAD, &communicateWithServer, NULL);
				break;
			case TE_COMM_KILL:	/* Kill hung communication thread */
				timed_events[i].active = FALSE;
				GLOBAL_SEND_MSG_COUNT++;
				if (LOCKED_WORK_UNIT != NULL)
					decrementWorkUnitUseCount (LOCKED_WORK_UNIT, SHORT_TERM_USE);
//bug				gwthread_kill (&COMMUNICATION_THREAD);
				break;
			case TE_PRIORITY_WORK:	/* Check for priority work event */
				timed_events[i].time_to_fire = this_time + TE_PRIORITY_WORK_FREQ;
				check_for_priority_work ();
				break;
			case TE_COMPLETION_DATES:	/* Send expected completion dates event */
				timed_events[i].active = FALSE;
				UpdateEndDates ();
				break;
			case TE_THROTTLE:	/* Sleep due to Throttle=n event */
				timed_events[i].time_to_fire =
					this_time + handleThrottleTimerEvent ();
				break;
			case TE_SAVE_FILES:	/* Timer to trigger writing save files */
						/* Also check for add files */
				timed_events[i].active = FALSE;
				if (addFileExists ()) stop_workers_for_add_files ();
				saveFilesTimer ();
				break;
			case TE_BATTERY_CHECK:	/* Test battery status */
				timed_events[i].time_to_fire = this_time + TE_BATTERY_CHECK_FREQ;
				test_battery ();
				break;
			case TE_ROLLING_AVERAGE: /* Adjust rolling average event */
				timed_events[i].time_to_fire = this_time + TE_ROLLING_AVERAGE_FREQ;
				adjust_rolling_average ();
				break;
			case TE_READ_PAUSE_DATA: /* Reread PauseWhileRunning info */
				timed_events[i].active = FALSE;
				read_pause_info ();
				break;
			case TE_READ_INI_FILE: /* Reread Ini files */
				timed_events[i].active = FALSE;
				stop_workers_for_reread_ini ();
				break;
			case TE_LOAD_AVERAGE:	/* Check load average event */
				timed_events[i].active = FALSE;
				checkLoadAverage ();
				break;
			}
		}

/* Loop to calculate sleep time until next event */

	}
}
