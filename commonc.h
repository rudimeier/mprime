/* Copyright 1995-2015 Mersenne Research, Inc.  All rights reserved */

/* Constants */

#define VERSION		"28.7"
#define BUILD_NUM	"1"
/* The list of assigned OS ports follows: */
/* Win9x (prime95) #1 */
/* Linux (mprime)  #2 */
/* Solaris	   #3 (never happened) */
/* Win 64-bit	   #4 */
/* WinNT (ntprime) #5 */
/* FreeBSD(mprime) #6 */
/* OS/2		   #7 */
/* Linux x86-64	   #8 */
/* Mac OS X	   #9 */
/* Mac OS X 64-bit #10 */
/* Haiku	   #11 */
/* FreeBSD 64-bit  #12 */

#define MIN_PRIME	5L		/* Smallest testable prime */
#define MAX_FACTOR	2000000000	/* Largest factorable Mersenne number*/
#define ERROR_RATE	0.018		/* Estimated error rate on clean run */

/* We want to increase the max number of worker threads to 64, but 32-bit Windows */
/* SetThreadAffinityMask only supports 32-bit masks.  So we'll assume any user that */
/* has a machine powerful enough to need more than 32 workers will run a 64-bit OS. */
/* 2014 UPDATE: Machines with more than 64 logical CPUs are in the works.  Windows has */
/* a bizarre affinity mechanism so we'll keep the 64 worker restriction there. */
/* On other OSes, we'll increase the limit to 512. */

#ifdef X86_64
#ifdef _WINDOWS_
#define MAX_NUM_WORKER_THREADS 64	/* Number of launchable work threads */
#else
#define MAX_NUM_WORKER_THREADS 512	/* Number of launchable work threads */
#endif
#else
#define MAX_NUM_WORKER_THREADS 32	/* Number of launchable work threads */
#endif

/* Factoring limits based on complex formulas given the speed of the */
/* factoring code vs. the speed of the Lucas-Lehmer code */
/* As an example, examine factoring to 2^68 (finding all 68-bit factors). */
/* First benchmark a machine to get LL iteration times and trial factoring */
/* times for a (16KB sieve of p=35000011). */
/*	We want to find when time spend eliminating an exponent with */
/* trial factoring equals time saved running 2 LL tests. */

/*	runs to find a factor (68) *
	#16KB sections (2^68-2^67)/p/(120/16)/(16*1024*8) *
	factoring_benchmark = 2.0 * LL test time (p * ll_benchmark)

	simplifying:

	68 * (2^68-2^67)/p/(120/16)/(16*1024*8) * facbench = 2 * p * llbench
	68 * 2^67 / p / (120/16) / 2^17 * facbench = 2 * p * lltime
	68 * 2^49 / p / (120/16) * facbench = p * lltime
	68 * 2^49 / (120/16) * facbench = p^2 * lltime
	68 * 2^53 / 120 * facbench = p^2 * lltime
	68 * 2^53 / 120 * facbench / lltime = p^2
	sqrt (68 * 2^53 / 120 * facbench / lltime) = p
*/

/* Now lets assume 30% of these factors would have been found by P-1.  So
   we only save a relatively quick P-1 test instead 2 LL tests.  Thus:
	sqrt (68 / 0.7 * 2^53 / 120 * facbench / lltime) = p
*/

/* Now factor in that 35000000 does 19 squarings, but 70000000 requires 20.
   Thus, if maxp is the maximum exponent that can be handled by an FFT size:
	sqrt (68 / 0.7 * 2^53 / 120 *
	      facbench * (1 + LOG2 (maxp/35000000) / 19) / lltime) = p
*/

/* Now factor in that errors sometimes force us to run more than 2 LL tests.
   Assume, 2.04 on average:
	sqrt (68 / 0.7 * 2^53 / 120 *
	      facbench * (1 + LOG2 (maxp/35000000) / 19) / lltime / 1.02) = p
*/

/* These breakeven points we're calculated on a 2.5 GHz Core 2 using 64-bit prime95 v26.6: */

#define FAC80	600000000L	/* The actual FAC80 breakeven point is even higher! */
#define FAC79	516800000L
#define FAC78	408400000L
#define FAC77	322100000L
#define FAC76	253500000L
#define FAC75	199500000L
#define FAC74	153400000L
#define FAC73	120000000L
#define FAC72	96830000L
#define FAC71	77910000L
#define FAC70	60940000L
#define FAC69	48800000L
#define FAC68	38300000L
#define FAC67	29690000L	/* We didn't bother calculating any smaller breakevens */

/* These breakeven points we're calculated on a 2.0 GHz P4 Northwood (using v24?): */

//#define FAC80	516000000L
//#define FAC79	420400000L
//#define FAC78	337400000L
//#define FAC77	264600000L
//#define FAC76	227300000L
//#define FAC75	186400000L
//#define FAC74	147500000L
//#define FAC73	115300000L
//#define FAC72	96830000L
//#define FAC71	75670000L
//#define FAC70	58520000L
//#define FAC69	47450000L
//#define FAC68	37800000L
//#define FAC67	29690000L
#define FAC66	23390000L

/* These breakevens we're calculated a long time ago on unknown hardware: */

#define FAC65	13380000L
#define FAC64	8250000L
#define FAC63	6515000L
#define FAC62	5160000L
#define FAC61	3960000L
#define FAC60	2950000L
#define FAC59	2360000L
#define FAC58	1930000L
#define FAC57	1480000L
#define FAC56	1000000L

/* Global variables */

extern int USE_V4;
extern char V4_USERID[15];
extern char V4_USERPWD[9];
extern char V4_USERNAME[80];

extern char INI_FILE[80];		/* Name of the prime INI file */
extern char LOCALINI_FILE[80];		/* Name of the local INI file */
extern char WORKTODO_FILE[80];		/* Name of the work-to-do INI file */
extern char RESFILE[80];		/* Name of the results file */
extern char SPOOL_FILE[80];		/* Name of the spool file */
extern char LOGFILE[80];		/* Name of the server log file */

extern char USERID[21];			/* User's ID */
extern char COMPID[21];			/* Computer name */
extern char COMPUTER_GUID[33];		/* Global unique computer ID */
extern int USE_PRIMENET;		/* TRUE if we're using PrimeNet */
extern int DIAL_UP;			/* TRUE if we're dialing into */
					/* PrimeNet server */
extern unsigned int NUM_WORKER_THREADS; /* Number of work threads to launch */
extern unsigned int WORK_PREFERENCE[MAX_NUM_WORKER_THREADS];
					/* Type of work (factoring, testing, */
					/* etc.) to get from the server. */
extern unsigned int CPU_AFFINITY[MAX_NUM_WORKER_THREADS];
					/* NT Processor affinity */
extern unsigned int THREADS_PER_TEST[MAX_NUM_WORKER_THREADS];
					/* Number of threads gwnum can use */
					/* in computations. */
extern unsigned int DAYS_OF_WORK;	/* How much work to retrieve from */
					/* the primenet server */
extern int STRESS_TESTER;		/* 1 if stress testing */
extern int volatile ERRCHK;		/* 1 to turn on roundoff error checking */
extern int volatile SUM_INPUTS_ERRCHK;	/* 1 to turn on sum(inputs) != sum(outputs) error checking */
extern unsigned int PRIORITY;		/* Desired priority level */
extern int MANUAL_COMM;			/* Set on if user explicitly starts */
					/* all communication with the server */
extern unsigned int volatile CPU_HOURS;	/* Hours per day program will run */
extern int CLASSIC_OUTPUT;		/* LL and PRP output to worker windows should use the pre-v28.5 classic style */
extern int OUTPUT_ROUNDOFF;		/* LL and PRP output to worker windows shound include the roundoff error */
extern unsigned long volatile ITER_OUTPUT;/* Iterations between outputs */
extern unsigned long volatile ITER_OUTPUT_RES;/* Iterations between results */
					/* file outputs */
extern unsigned long volatile DISK_WRITE_TIME;
					/* Number of minutes between writing */
					/* intermediate results to disk */
extern unsigned int MODEM_RETRY_TIME;	/* How often to try sending msgs */
					/* to primenet server whem modem off */
extern unsigned int NETWORK_RETRY_TIME;	/* How often to try sending msgs */
					/* to primenet server */
extern float DAYS_BETWEEN_CHECKINS;	/* Days between sending updated */
					/* completion dates to the server */
extern int NUM_BACKUP_FILES;		/* Between 1 and 3 backup files (or 99 */
					/* for overwrite) */
extern int SILENT_VICTORY;		/* Quiet find of new Mersenne prime */
extern int SILENT_VICTORY_PRP;		/* Quiet find of new PRP */
extern int RUN_ON_BATTERY;		/* Run program even on battery power */
extern int BATTERY_PERCENT;		/* Pause if battery below this percent charged */
extern int DEFEAT_POWER_SAVE;		/* If possible, call OS to force computer to run at full speed */
extern int TRAY_ICON;			/* Display tiny tray icon */
extern int HIDE_ICON;			/* Display no icon */
extern int MERGE_WINDOWS;		/* Flags indicating which MDI */
					/* windows to merge together */
#define MERGE_MAIN_WINDOW	0x1	/* Merge main into first worker */
#define MERGE_COMM_WINDOW	0x2	/* Merge comm into first worker */
#define MERGE_WORKER_WINDOWS	0x4	/* Merge all workers into one window */
#define MERGE_MAINCOMM_WINDOWS	0x8	/* Merge main and comm windows */
#define MERGE_NO_PREFIX		0x20	/* Output thread prefix on each line flag */

extern double UNOFFICIAL_CPU_SPEED;	/* Last calculated CPU Speed in MHz.  We only */
					/* lower the official CPU_SPEED (the one reported to */
					/* the server) after several lower measurements */
extern unsigned int ROLLING_AVERAGE;	/* Ratio of this computer's speed */
					/* compared to the expected speed */
					/* for this CPU */
extern unsigned int PRECISION;		/* Number of decimal places to output*/
					/* in percent complete lines */
extern int RDTSC_TIMING;		/* True if RDTSC is used to time */
extern int TIMESTAMPING;		/* True is timestamps to be output */
extern int CUMULATIVE_TIMING;		/* True if outputting cumulative time*/
extern int SEQUENTIAL_WORK;		/* TRUE (the  default) if undocumented */
					/* SequentialWorkToDo is set */
extern int WELL_BEHAVED_WORK;		/* TRUE if undocumented feature */
					/* "well behaved worktodo file" */
					/* is on.  This reduces the number */
					/* of times worktodo.ini is read */
					/* and written. */
extern unsigned long INTERIM_FILES;	/* Create save file every N iters */  
extern unsigned long INTERIM_RESIDUES;	/* Print residue every N iterations */
extern unsigned long HYPERTHREADING_BACKOFF; /* Pause prime95 if iterations */
					/* get too slow. */
extern int THROTTLE_PCT;		/* Percent CPU time prog should run */

extern int STARTUP_IN_PROGRESS;		/* TRUE if startup dialogs are up */

extern unsigned long NUM_CPUS;		/* Number of CPUs/Cores in computer */

extern int LAUNCH_TYPE;			/* Type of worker threads launched */
extern unsigned int WORKER_THREADS_ACTIVE;/* Num worker threads running */
extern int WORKER_THREADS_STOPPING;	/* TRUE iff worker threads stopping */

extern unsigned int WORKTODO_COUNT;	/* Count of valid work lines */

extern int GIMPS_QUIT;			/* TRUE if we just successfully */
					/* quit the GIMPS project */

extern gwthread COMMUNICATION_THREAD;	/* Handle for comm thread.  Set when */
					/* comm thread is active. */

/* Common routines */

void generate_application_string (char *);
void getCpuInfo (void);
void getCpuDescription (char *, int);

int isPrime (unsigned long p);
unsigned int strToMinutes (const char	*);
void minutesToStr (unsigned int, char *);
void write_memory_settings (unsigned int, unsigned int, unsigned int, unsigned int);
int read_memory_settings (unsigned int *, unsigned int *, unsigned int *, unsigned int *);

void nameAndReadIniFiles (int named_ini_files);
void initCommCode (void);
int readIniFiles (void);

void IniSectionGetString (const char *, const char *, const char *, char *, unsigned int, const char *);
void IniSectionGetTimedString (const char *, const char *, const char *, char *, unsigned int, const char *, unsigned int *);
long IniSectionGetInt (const char *, const char *, const char *, long);
long IniSectionGetTimedInt (const char *, const char *, const char *, long, unsigned int *);
void IniSectionWriteString (const char *, const char *, const char *, const char *);
void IniSectionWriteInt (const char *, const char *, const char *, long);
void IniSectionWriteFloat (const char *, const char *, const char *, float);
void IniGetString (const char *, const char *, char *, unsigned int, const char *);
void IniGetTimedString (const char *, const char *, char *, unsigned int, const char *, unsigned int *);
long IniGetInt (const char *, const char *, long);
long IniGetTimedInt (const char *, const char *, long, unsigned int *);
void IniWriteString (const char *, const char *, const char *);
void IniWriteInt (const char *, const char *, long);
void IniWriteFloat (const char *, const char *, float);
void IniFileReread (const char *);
void processTimedIniFile (const char *);

void parse_timed_ini_value (const char *, unsigned int *, unsigned int *, unsigned int *);
const char *IniSectionGetStringRaw (const char *, const char *, const char *);

int addFileExists (void);
void incorporateIniAddFiles (void);
int incorporateWorkToDoAddFile (void);
void iniAddFileMerge (char *, char *, char *);

void PTOGetAll (char *ini_filename, char *keyword, unsigned int *array,
		unsigned int def_val);
void PTOSetAll (char *ini_filename, char *keyword, char	*shadow_keyword,
		unsigned int *array, unsigned int new_val);
void PTOSetOne (char *ini_filename, char *keyword, char *shadow_keyword,
		unsigned int *array, int tnum, unsigned int new_val);
int PTOIsGlobalOption (unsigned int *array);
int PTOHasOptionChanged (char *shadow_keyword, unsigned int *array, int tnum);


#define MAIN_THREAD_NUM		-2
#define COMM_THREAD_NUM		-1
void create_window (int thread_num);
void destroy_window (int thread_num);
void TileViews (void);
void base_title (int, const char *);
void title (int, const char *);
#define	WORKING_ICON	0
#define	IDLE_ICON	1
void ChangeIcon (int, int);
void BlinkIcon (int, int);
EXTERNC void OutputBoth (int, const char *);
void OutputStr (int, const char *);
void OutputStrNoTimeStamp (int, const char *);
void RealOutputStr (int, const char *);
void OutputSomewhere (int, const char *);
void LogMsg (const char *);
int OutOfMemory (int);

/* Structures and definitions dealing with the worktodo.ini file */

#define WORK_FACTOR		0
#define WORK_TEST		1
#define WORK_DBLCHK		2
#define WORK_ADVANCEDTEST	3
#define WORK_ECM		4
#define WORK_PMINUS1		5
#define WORK_PFACTOR		6
#define WORK_PRP		10
#define WORK_NONE		100	/* Comment line in worktodo.ini */
#define WORK_DELETED		101	/* Deleted work_unit */

struct work_unit {		/* One line from the worktodo file */
	int	work_type;	/* Type of work to do */
	char	assignment_uid[33]; /* Primenet assignment ID */
	char	extension[9];	/* Optional save file extension */
	double	k;		/* K in k*b^n+c */
	unsigned long b;	/* B in k*b^n+c */
	unsigned long n;	/* N in k*b^n+c */
	signed long c;		/* C in k*b^n+c */
	unsigned long forced_fftlen;/* Forced FFT length to use.  Primarily */
				/* used for implementing soft FFT */
				/* crossovers.  Zero means default fftlen */
	double	sieve_depth;	/* How far it has been trial factored */
	double	factor_to;	/* How far we should trial factor to */
	int	pminus1ed;	/* TRUE if has been P-1 factored */
	double	B1;		/* ECM and P-1 - Stage 1 bound */
	double	B2_start;	/* ECM and P-1 - Stage #2 start */
	double	B2;		/* ECM and P-1 - Stage #2 end */
	unsigned int curves_to_do; /* ECM - curves to try */
	double	curve;		/* ECM - Specific curve to test (debug tool) */
	double	tests_saved;	/* Pfactor - primality tests saved if */
				/* a factor is found */
	char	*known_factors;	/* ECM, P-1, PRP - list of known factors */
	char	*comment;	/* Comment line in worktodo.ini */
		/* Runtime variables */
	struct work_unit *next; /* Next in doubly-linked list */
	struct work_unit *prev; /* Previous in doubly-linked list */
	int	in_use_count;	/* Count of threads accessing this work unit */
	int	high_memory_usage;/* Set if we are using a lot of memory */
				/* If user changes the available memory */
				/* settings, then we should stop and */
				/* restart our computations */
	char	stage[10];	/* Test stage (e.g. TF,P-1,LL) */
	double	pct_complete;	/* Percent complete (misnomer as value is */
				/* between 0.0 and 1.0) */
	unsigned long fftlen;	/* FFT length in use */
	int	ra_failed;	/* Set when register assignment fails, tells */
				/* us not to try registering it again. */
};
struct work_unit_array {	/* All the lines for one worker thread */
	struct work_unit *first; /* First work unit */
	struct work_unit *last;	/* Last work unit */
};

int readWorkToDoFile (void);
int writeWorkToDoFile (int);
#define SHORT_TERM_USE		0
#define LONG_TERM_USE		1
struct work_unit *getNextWorkToDoLine (int, struct work_unit *, int);
void decrementWorkUnitUseCount (struct work_unit *, int);
int addWorkToDoLine (int, struct work_unit *);
int updateWorkToDoLine (int, struct work_unit *);
int deleteWorkToDoLine (int, struct work_unit *, int);
int isWorkUnitActive (struct work_unit *);
int addToWorkUnitArray (unsigned int, struct work_unit *, int);

void rolling_average_work_unit_complete (int, struct work_unit *);
void invalidateNextRollingAverageUpdate (void);

/* More miscellaneous routines */

double work_estimate (int thread_num, struct work_unit *);
unsigned int factorLimit (struct work_unit *);
void guess_pminus1_bounds (int, double, unsigned long, unsigned long, signed long,
			   double, double, unsigned long *,
			   unsigned long *, unsigned long *, double *);

void strupper (char *);
void tempFileName (struct work_unit *, char *);
int fileExists (char *);

int read_array (int fd, char *buf, unsigned long len, unsigned long *sum);
int write_array (int fd, char *buf, unsigned long len, unsigned long *sum);
int read_gwnum (int fd, gwhandle *gwdata, gwnum g, unsigned long *sum);
int write_gwnum (int fd, gwhandle *gwdata, gwnum g, unsigned long *sum);
int read_short (int fd, short *val);
int read_long (int fd, unsigned long *val, unsigned long *sum);
int write_long (int fd, unsigned long val, unsigned long *sum);
int read_slong (int fd, long *val, unsigned long *sum);
int write_slong (int fd, long val, unsigned long *sum);
int read_longlong (int fd, uint64_t *val, unsigned long *sum);
int write_longlong (int fd, uint64_t val, unsigned long *sum);
int read_double (int fd, double *val, unsigned long *sum);
int write_double (int fd, double dbl, unsigned long *sum);
int read_magicnum (int fd, unsigned long magicnum);
int read_header (int fd, unsigned long *version,
		 struct work_unit *w, unsigned long *sum);
int write_header (int fd, unsigned long magicnum, unsigned long version,
		  struct work_unit *w);
int read_checksum (int fd, unsigned long *sum);
int write_checksum (int fd, unsigned long sum);

void formatMsgForResultsFile (char *, struct work_unit *);
int writeResults (const char *);

/* Routines called by common routines */

unsigned long physical_memory (void);
unsigned long GetSuggestedMemory (unsigned long nDesiredMemory);
int getDefaultTimeFormat (void);

/******************************************************************************
*                 Spool File and Server Communication Code                    *
******************************************************************************/

void init_spool_file_and_comm_code (void);
void set_comm_timers (void);
void clear_comm_rate_limits (void);
void do_manual_comm_now (void);
void pingServer (void);
void UpdateEndDates (void);
void ConditionallyUpdateEndDates (void);
#define MSG_CHECK_WORK_QUEUE		998
#define MSG_QUIT_GIMPS			999
void spoolMessage (short, void *);
void spoolExistingResultsFile (void);
int unreserve (unsigned long);
void salvageCorruptSpoolFile (void);

int LoadPrimeNet (void);
void UnloadPrimeNet (void);
int PRIMENET (short, void *);

/******************************************************************************
*                           Timed Events Handler                              *
******************************************************************************/

#define TE_MEM_CHANGE		0	/* Night/day memory change event */
#define TE_PAUSE_WHILE		1	/* Check pause_while_running event */
#define TE_WORK_QUEUE_CHECK	2	/* Check work queue event */
#define TE_COMM_SERVER		3	/* Retry communication with server event */
#define TE_COMM_KILL		4	/* Kill hung communication thread event */
#define TE_PRIORITY_WORK	5	/* Check for priority work event */
#define TE_COMPLETION_DATES	6	/* Send expected completion dates event */
#define TE_THROTTLE		7	/* Sleep due to Throttle=n event */
#define TE_SAVE_FILES		8	/* Trigger the writing of save files */
#define TE_BATTERY_CHECK	9	/* Check battery status frequently */
#define TE_ROLLING_AVERAGE	10	/* Adjust rolling average */
#define TE_READ_PAUSE_DATA	11	/* Reread PauseWhileRunning info */
#define TE_READ_INI_FILE	12	/* Reread prime.txt settings because */
					/* a during/else time period has ended */
#define TE_LOAD_AVERAGE		13	/* Linux/FreeBSD/Apple load average check */

#define MAX_TIMED_EVENTS	14	/* Maximum number of timed events */

void init_timed_event_handler (void);

void add_timed_event (
	int	event_number,		/* Which event to add */
	int	time_to_fire);		/* When to start event (seconds from now) */

void delete_timed_event (
	int	event_number);		/* Which event to delete */

int is_timed_event_active (
	int	event_number);		/* Which event to test */

time_t timed_event_fire_time (
	int	event_number);		/* Which event to get fire time of */

#define TE_WORK_QUEUE_CHECK_FREQ 6*60*60 /* Check work queue every 6 hrs. */
#define TE_PRIORITY_WORK_FREQ	 6*60*60 /* Check priority work every 6 hrs. */
#define TE_BATTERY_CHECK_FREQ	 15	/* Check battery every 15 sec. */
#define TE_THROTTLE_FREQ	 5	/* Throttle every 5 sec. */
#define TE_ROLLING_AVERAGE_FREQ	 12*60*60 /* Adjust rolling every 12 hr. */
