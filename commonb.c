/*----------------------------------------------------------------------
| This file contains routines and global variables that are common for
| all operating systems the program has been ported to.  It is included
| in one of the source code files of each port.
|
| Commona contains information used only during setup
| Commonb contains information used only during execution
| Commonc contains information used during setup and execution
+---------------------------------------------------------------------*/
 
/* Globals for error messages */

char ERRMSG0[] = "Iteration: %ld/%ld, %s";
char ERRMSG1A[] = "ERROR: ILLEGAL SUMOUT\n";
char ERRMSG1B[] = "ERROR: SUM(INPUTS) != SUM(OUTPUTS), %.16g != %.16g\n";
char ERRMSG1C[] = "ERROR: ROUND OFF (%.10g) > 0.40\n";
char ERRMSG1D[] = "ERROR: Shift counter corrupt.\n";
char ERRMSG1E[] = "ERROR: Illegal double encountered.\n";
char ERRMSG1F[] = "ERROR: FFT data has been zeroed!\n";
char ERRMSG2[] = "Possible hardware failure, consult the readme.txt file.\n";
char ERRMSG3[] = "Continuing from last save file.\n";
char ERRMSG4[] = "Waiting five minutes before restarting.\n";
char ERRMSG5[] = "For added safety, redoing iteration using a slower, more reliable method.\n";
char ERROK[] = "Disregard last error.  Result is reproducible and thus not a hardware problem.\n";
char READFILEERR[] = "Error reading intermediate file: %s\n";
char WRITEFILEERR[] = "Error writing intermediate file: %s\n";
char RENAME_MSG[] = "Renaming intermediate file %s to %s.\n";

/* PauseWhileRunning globals */

struct pause_info {
	int	thread_num;		/* Worker thread to pause */
	int	low_mem;		/* Flag set for LowMemWhileRunning entries */
	int	workers_affected;	/* Number of workers affected */
	char	*program_name;		/* Pause if running this program */
	char	matching_program[80];	/* Running program that matched this entry */
	struct pause_info *next;	/* Next in linked list of program names */
};

int	PAUSE_MUTEX_INITIALIZED = 0;
gwmutex	PAUSE_MUTEX;		/* Lock for accessing pause globals */
struct pause_info *PAUSE_DATA = NULL;
int	PAUSE_WHILE_RUNNING_FREQ = 10;

/* Globals for stopping and starting worker threads */

/* Note that we have one flag byte for each worker thread.  We could */
/* use one bit per worker thread, but then we need to have locks around */
/* updates so that 2 worker threads don't interleave a read-modify-write */
/* operation. */

int	STOP_FOR_RESTART = FALSE;/* Flag indicating we should stop and */
				/* restart all worker threads because an */
				/* important option changed in the GUI. */
				/* One example is changing the priority */
				/* for worker threads. */
int	STOP_FOR_REREAD_INI = FALSE;/* Flag indicating all workers must */
				/* stop because a during/else time period */
				/* has ended and INI file must be reread. */
char	STOP_FOR_MEM_CHANGED[MAX_NUM_WORKER_THREADS] = {0};
				/* Flags indicating it is time to stop */
				/* workers due to day/night memory change. */
int	STOP_FOR_BATTERY = FALSE;/* Flag indicating it is time to stop */
				/* workers due to running on battery. */
char	STOP_FOR_PRIORITY_WORK[MAX_NUM_WORKER_THREADS] = {0};
				/* Flags indicating it is time to switch */
				/* a worker to high priority work. */
struct pause_info *STOP_FOR_PAUSE[MAX_NUM_WORKER_THREADS] = {NULL};
				/* Flags saying worker thread should */
				/* pause while another program runs */
struct pause_info *STOP_FOR_LOW_MEMORY = NULL; /* Set when LowMemWhileRunning active */
				/* Workers using lots of memory will be stopped */
char	STOP_FOR_THROTTLE[MAX_NUM_WORKER_THREADS] = {0};
				/* Flags indicating it is time to pause */
				/* a worker for throttling. */
char	STOP_FOR_ABORT[MAX_NUM_WORKER_THREADS] = {0};
				/* Abort work unit due to unreserve, factor */
				/* found in a different thread, server */
				/* request, or any other reason. */
char	ACTIVE_WORKERS[MAX_NUM_WORKER_THREADS] = {0};
				/* Flags indicating which worker threads */
				/* are active. */
char	WRITE_SAVE_FILES[MAX_NUM_WORKER_THREADS] = {0};
				/* Flags indicating it is time to write */
				/* a save file. */

char	WORK_AVAILABLE_OR_STOP_INITIALIZED[MAX_NUM_WORKER_THREADS] = {0};
gwevent WORK_AVAILABLE_OR_STOP[MAX_NUM_WORKER_THREADS] = {0};
				/* Signal for telling primeContinue that */
				/* work is now available or all threads */
				/* are stopping */
char	USER_START_OR_STOP_INITIALIZED[MAX_NUM_WORKER_THREADS] = {0};
gwevent	USER_START_OR_STOP[MAX_NUM_WORKER_THREADS] = {0};
				/* Signal for telling implement_stop_one_worker */
				/* that the user wants this worker to start or */
				/* all threads are stopping */
char	END_PAUSE_OR_STOP_INITIALIZED[MAX_NUM_WORKER_THREADS] = {0};
gwevent END_PAUSE_OR_STOP[MAX_NUM_WORKER_THREADS] = {0};
				/* Signal for telling implement_pause */
				/* that the pause has ended or */
				/* all threads are stopping */
char	OFF_BATTERY_OR_STOP_INITIALIZED[MAX_NUM_WORKER_THREADS] = {0};
gwevent OFF_BATTERY_OR_STOP[MAX_NUM_WORKER_THREADS] = {0};
				/* Signal for telling implement_stop_battery */
				/* that AC power has been restored or */
				/* all threads are stopping */
char	MEM_WAIT_OR_STOP_INITIALIZED[MAX_NUM_WORKER_THREADS] = {0};
gwevent MEM_WAIT_OR_STOP[MAX_NUM_WORKER_THREADS] = {0};
				/* Signal for telling avail_mem that */
				/* it can now determine the available memory */

/* Globals for memory manager */

#define DEFAULT_MEM_USAGE 24	/* 24MB default */
unsigned long AVAIL_MEM = 0;	/* Memory available now */
unsigned long MAX_MEM = 0;	/* Max memory available */
unsigned long AVAIL_MEM_PER_WORKER[MAX_NUM_WORKER_THREADS] = {0};
				/* Maximum memory each worker can use */
unsigned long MAX_HIGH_MEM_WORKERS =  0; /* Maximum number of workers */
				/* allowed to use lots of memory */
#define MEM_USAGE_NOT_SET 0x1	/* The mem_in_use value is just a guess */
				/* as the work unit for the thread has not */
				/* started yet or is restarting. */
#define MEM_RESTARTING	0x2	/* The mem_in_use value will be changing */
				/* soon as the thread is being restarted */
				/* because it was using too much memory. */
#define MEM_WILL_BE_VARIABLE_USAGE 0x4
				/* The current work unit will be a */
				/* variable memory user.  We just don't */
				/* know how much it will use yet. */
#define	MEM_VARIABLE_USAGE 0x8	/* The current work unit is using a */
				/* lot of memory now and if needed could */
				/* use less if restarted. */
#define MEM_WAITING	0x10	/* Set if thread is waiting for another thread */
				/* to stop before returning from set_memory_usage */
char	MEM_FLAGS[MAX_NUM_WORKER_THREADS] = {0};
				/* Flags indicating which threads will be */
				/* affected by a change in memory settings. */
unsigned int MEM_IN_USE[MAX_NUM_WORKER_THREADS] = {0};
				/* Array containing memory in use by each */
				/* worker thread */

#define MEM_RESTART_LOWMEM_ENDS 0x1
				/* Worker needs to restart when */
				/* the LowMemWhileRunning program ends. */
#define MEM_RESTART_MAX_MEM_AVAILABLE 0x2
				/* Worker needs to restart when available memory */
				/* equals maximum memory.  This happens when */
				/* stage 2 is delayed until max memory is available. */
#define MEM_RESTART_MAX_MEM_CHANGE 0x4
				/* Current work unit needs to restart if */
				/* max mem changes.  P-1 may choose different */
				/* bounds because of the change */
#define MEM_RESTART_TOO_MANY_HIGHMEM 0x8
				/* Worker needs to restart because */
				/* MAX_HIGH_MEM_WORKERS exceeded. */
#define MEM_RESTART_MORE_AVAIL 0x10 /* One of the worker's work units did not */
				/* have enough memory to run.  If memory */
				/* becomes available restart the worker. */
#define MEM_RESTART_IF_MORE 0x20 /* The current work unit could use more memory */
				/* and should be restarted if more becomes */
				/* available. */

char	MEM_RESTART_FLAGS[MAX_NUM_WORKER_THREADS] = {0};
unsigned int MEM_RESTART_MIN_AMOUNT[MAX_NUM_WORKER_THREADS] = {0};
unsigned int MEM_RESTART_DESIRED_AMOUNT[MAX_NUM_WORKER_THREADS] = {0};
				/* Only restart if this amount of memory  */
				/* is available */
unsigned int MEM_RESTART_IF_MORE_AMOUNT[MAX_NUM_WORKER_THREADS] = {0};

int	MEM_MUTEX_INITIALIZED = FALSE;
gwmutex	MEM_MUTEX;		/* Lock for accessing mem globals */

/**************************************************************/
/*    Routines dealing with thread priority and affinity      */
/**************************************************************/

/* Set the thread priority correctly.  Most screen savers run at priority 4. */
/* Most application's run at priority 9 when in foreground, 7 when in */
/* background.  In selecting the proper thread priority I've assumed the */
/* program usually runs in the background. */ 

/* This routine is also responsible for setting the thread's CPU affinity. */
/* If there are N cpus with hyperthreading, then physical cpu 0 is logical */
/* cpu 0 and N, physical cpu 1 is logical cpu 1 and N+1, etc. */

void SetPriority (
	struct PriorityInfo *info)
{
	unsigned int i;
	int	mask;

/* Benchmarking affinity.  For hyperthreaded CPUs, we put auxillary */
/* threads onto the logical CPUs before moving onto the next physical CPU. */  

	if (info->type == SET_PRIORITY_BENCHMARKING) {
		mask = 1 << (info->aux_thread_num / CPU_HYPERTHREADS);
		for (i = 1; i < CPU_HYPERTHREADS; i++)
			mask |= (mask << NUM_CPUS);
	}

/* Torture test affinity.  If we're running the same number of torture */
/* tests as CPUs if the system then set affinity.  Otherwise, let the */
/* threads run on any CPU. */

	else if (info->type == SET_PRIORITY_TORTURE) {
		if (info->num_threads == NUM_CPUS * CPU_HYPERTHREADS) {
			mask = 1 << info->thread_num;
		} else if (info->num_threads == NUM_CPUS) {
			mask = 1 << (info->thread_num / CPU_HYPERTHREADS);
			for (i = 1; i < CPU_HYPERTHREADS; i++)
				mask |= (mask << NUM_CPUS);
		} else
			mask = -1;
	}

/* QA affinity.  Just let the threads run on any CPU. */

	else if (info->type == SET_PRIORITY_QA) {
		mask = -1;
	}

/* Normal worker threads.  Pay attention to the affinity option set */
/* by the user. */

/* A CPU_AFFINITY setting of 99 means "run on any CPU". */

	else if (CPU_AFFINITY[info->thread_num] == 99) {
 		mask = -1;
	}

/* A small CPU_AFFINITY setting means run only on that CPU.  Since there is */
/* no way to explicitly tell us which CPU to run an auxillary thread on, */
/* we put the first auxillary thread on any hyperthreaded logical CPUs. */
/* The next auxillary thread goes on the next physical CPU. */

	else if (CPU_AFFINITY[info->thread_num] < 100) {
		mask = 1 << CPU_AFFINITY[info->thread_num];
		mask <<= ((info->aux_thread_num % CPU_HYPERTHREADS) *
							CPU_HYPERTHREADS);
		mask <<= (info->aux_thread_num / CPU_HYPERTHREADS);
	}

/* A CPU_AFFINITY setting of 100 means "smart affinity assignments". */
/* We've now reached that case. */

/* If all worker threadts are not set to smart affinity, then it is */
/* just too hard to figure out what is best.  Just let the OS run the */
/* threads on any CPU. */

	else if (! PTOIsGlobalOption (CPU_AFFINITY))
	 	mask = -1;

/* If number of worker threads equals number of logical cpus then run each */
/* worker thread on its own logical CPU.  If the user also has us running */
/* auxillary threads, then the user has made a really bad decision and a */
/* performance hit will occur. */

	else if (NUM_WORKER_THREADS == NUM_CPUS * CPU_HYPERTHREADS) { 
		mask = 1 << info->thread_num;
		if (info->aux_thread_num) mask = -1;
	}

/* If number of worker threads equals number of physical cpus then run each */
/* worker thread on its own physical CPU.  Run auxillary threads on the same */
/* physical CPU.  This should be advantageous on hyperthreaded CPUs.  The */
/* should be careful to not run more auxillary threads than available */
/* logical CPUs created by hyperthreading. */ 

	else if (NUM_WORKER_THREADS == NUM_CPUS) {
		mask = 1 << (info->thread_num / CPU_HYPERTHREADS);
		for (i = 1; i < CPU_HYPERTHREADS; i++)
			mask |= (mask << NUM_CPUS);
	}

/* Otherwise, just run on any CPU. */

	else
		mask = -1;

/* Output an informative message */

	if (NUM_CPUS > 1 || CPU_HYPERTHREADS > 1) {
		char	buf[120];

		if (info->aux_thread_num == 0)
			strcpy (buf, "Setting affinity to run worker on ");
		else
			sprintf (buf, "Setting affinity to run helper thread %d on ", info->aux_thread_num);
		if (mask == -1) {
			strcat (buf, "any logical CPU.\n");
		} else {
			int	i, count;
			char	cpu_list[80];
			for (i = count = 0; i < 32; i++) {
				if (! (mask & (1 << i))) continue;
				count++;
				if (count == 1) sprintf (cpu_list, "%d", i);
				else sprintf (cpu_list + strlen(cpu_list), ",%d", i);
			}
			sprintf (buf + strlen(buf), "logical CPU%s%s\n",
				 (count == 1) ? " #" : "s ", cpu_list);
		}
		OutputStr (info->thread_num, buf);
	}

/* Call OS-specific routine to set the priority and affinity */

	setThreadPriorityAndAffinity (PRIORITY, mask);
}

/* Gwnum thread callback routine */

void SetAuxThreadPriority (int aux_thread_num, int action, void *data)
{
	struct PriorityInfo sp_info;

/* Handle thread start action.  Set the thread priority. */

	if (action == 0) {
		memcpy (&sp_info, data, sizeof (struct PriorityInfo));
		sp_info.aux_thread_num = aux_thread_num;
		SetPriority (&sp_info);
	}

/* Handle thread terminate action.  Remove thread handle from list */
/* of active worker threads. */

	if (action == 1) {
		registerThreadTermination ();
	}
}

/**************************************************************/
/*       Routines and globals dealing with stop codes         */
/*             and the write save files timer                 */
/**************************************************************/

/* This routine checks if the worker thread needs to be stopped for any */
/* reason whatsoever.  If the worker thread should stop, a stop reason */
/* is returned.  The routine is declared EXTERNC becasue it can be called */
/* by the C code in giants that does GCD. */

EXTERNC int stopCheck (
	int	thread_num)	/* Worker thread number */
{

/* Call an OS-specific callback routine.  This gives OSes that poll for */
/* the ESC key a hook.  They can perform any necessary miscellaneous */
/* functions and check for the ESC key to call stop_workers_for_escape. */

	stopCheckCallback (thread_num);

/* If the ESC key was hit, stop processing.  This takes precedence over */
/* all other stop reasons.  This also happens when the program is exiting. */

	if (WORKER_THREADS_STOPPING) return (STOP_ESCAPE);

/* If an important option changed in the GUI, restart all threads. */
/* For example, the user changes the priority of all worker threads. */	

	if (STOP_FOR_RESTART) return (STOP_RESTART);

/* If the during/else time period has ended, stop processing all worker */
/* threads so prime.txt and local.txt can be reprocessed. */

	if (STOP_FOR_REREAD_INI) return (STOP_REREAD_INI);

/* If the memory settings have changed, stop processing affected worker */
/* threads so they can allocate more or less memory. */

	if (STOP_FOR_MEM_CHANGED[thread_num]) {
		STOP_FOR_MEM_CHANGED[thread_num] = 0;
		return (STOP_MEM_CHANGED);
	}

/* If we are on battery power, stop processing all worker */
/* threads until we cease running on the battery. */

	if (STOP_FOR_BATTERY) return (STOP_BATTERY);

/* If the thread needs to go do some higher priority work, then stop */
/* processing this work_unit and reprocess the worktodo file. */

	if (STOP_FOR_PRIORITY_WORK[thread_num]) {
		STOP_FOR_PRIORITY_WORK[thread_num] = 0;
		return (STOP_PRIORITY_WORK);
	}

/* If the thread needs to abort the current work unit, then return */
/* that stop code. */

	if (STOP_FOR_ABORT[thread_num]) {
		STOP_FOR_ABORT[thread_num] = 0;
		return (STOP_ABORT);
	}

/* If the thread needs to stop because the user has explicitly stopped (or never */
/* started) this worker, then return the proper stop code. */

	if (!ACTIVE_WORKERS[thread_num]) return (STOP_WORKER);

/* Check if thread should pause because another process is running. */
/* When pause completes, check stop codes again.  We may have been paused */
/* a long time during which other stop timers may have fired. */

	if (STOP_FOR_PAUSE[thread_num] != NULL) {
		return (STOP_PAUSE);
// Do we want to offer INI option to do an immediate pause (next 2 lines) instead???
//		implement_pause (thread_num);
//		return (stopCheck (thread_num));
	}

/* If the thread needs to pause because of the throttle option, then */
/* do so now. */

	if (STOP_FOR_THROTTLE[thread_num]) {
		STOP_FOR_THROTTLE[thread_num] = 0;
		implementThrottle (thread_num);
	}

/* No need to stop */

	return (0);
}

/* Clear flags controlling the stopping of worker threads. */

int stop_events_initialized = FALSE;

void init_stop_code (void)
{
	STOP_FOR_RESTART = FALSE;
	STOP_FOR_REREAD_INI = FALSE;
	STOP_FOR_BATTERY = FALSE;
	memset (STOP_FOR_MEM_CHANGED, 0, sizeof (STOP_FOR_MEM_CHANGED));
	memset (STOP_FOR_PRIORITY_WORK, 0, sizeof (STOP_FOR_PRIORITY_WORK));
	memset (STOP_FOR_PAUSE, 0, sizeof (STOP_FOR_PAUSE));
	memset (STOP_FOR_THROTTLE, 0, sizeof (STOP_FOR_THROTTLE));
	memset (STOP_FOR_ABORT, 0, sizeof (STOP_FOR_ABORT));
	memset (WRITE_SAVE_FILES, 0, sizeof (WRITE_SAVE_FILES));
}

/* Signal threads waiting for work to do */

void restart_waiting_workers (
	int	restart_flags)
{
	int	thread_num;
	for (thread_num = 0; thread_num < MAX_NUM_WORKER_THREADS; thread_num++)
		restart_one_waiting_worker (thread_num, restart_flags);
}

void restart_one_waiting_worker (
	int	thread_num,
	int	restart_flags)
{
	if (restart_flags & RESTART_USER_START &&
	    USER_START_OR_STOP_INITIALIZED[thread_num]) {
		gwevent_signal (&USER_START_OR_STOP[thread_num]);
	}
	if (restart_flags & RESTART_WORK_AVAILABLE &&
	    WORK_AVAILABLE_OR_STOP_INITIALIZED[thread_num]) {
		gwevent_signal (&WORK_AVAILABLE_OR_STOP[thread_num]);
	}
	if (restart_flags & RESTART_END_PAUSE &&
	    END_PAUSE_OR_STOP_INITIALIZED[thread_num]) {
		gwevent_signal (&END_PAUSE_OR_STOP[thread_num]);
	}
	if (restart_flags & RESTART_MEM_WAIT &&
	    MEM_WAIT_OR_STOP_INITIALIZED[thread_num]) {
		gwevent_signal (&MEM_WAIT_OR_STOP[thread_num]);
	}
	if (restart_flags & RESTART_BATTERY &&
	    OFF_BATTERY_OR_STOP_INITIALIZED[thread_num]) {
		gwevent_signal (&OFF_BATTERY_OR_STOP[thread_num]);
	}
}

/* Set flags so that worker threads will stop due to ESC key being pressed. */

void stop_workers_for_escape (void)
{
	if (WORKER_THREADS_ACTIVE) {
		OutputStr (MAIN_THREAD_NUM, "Stopping all worker threads.\n");
		WORKER_THREADS_STOPPING = TRUE;
		restart_waiting_workers (RESTART_ALL);
		raiseAllWorkerThreadPriority ();
	}
}

/* Set flag so that all worker threads stop and restart because an */
/* important INI option changed (like thread priority).  This routine only */
/* restarts "genuine" work threads - not benchmarking and torture test */
/* work threads. */

void stop_workers_for_restart (void)
{
	if (WORKER_THREADS_ACTIVE &&
	    LAUNCH_TYPE == LD_CONTINUE &&
	    ! STOP_FOR_RESTART) {
		OutputStr (MAIN_THREAD_NUM, "Restarting all worker threads.\n");
		STOP_FOR_RESTART = TRUE;
		restart_waiting_workers (RESTART_ALL);
		if (NUM_WORKER_THREADS > WORKER_THREADS_ACTIVE)
			create_worker_windows (NUM_WORKER_THREADS);
	}
}

/* Set flag so that all worker threads stop and restart because an */
/* important INI option changed (like thread priority). */

void stop_workers_for_add_files (void)
{
	if (WORKER_THREADS_ACTIVE && ! STOP_FOR_RESTART) {
		OutputStr (MAIN_THREAD_NUM, "Restarting all worker threads to process .add file.\n");
		STOP_FOR_RESTART = TRUE;
		restart_waiting_workers (RESTART_ALL);
	}
}

/* Set flag so that worker threads will stop due to Time= end time being */
/* reached.  We need to stop all worker threads, reprocess prime.ini, and */
/* restart the worker threads. */

void stop_workers_for_reread_ini (void)
{
	if (WORKER_THREADS_ACTIVE && ! STOP_FOR_REREAD_INI) {
		OutputStr (MAIN_THREAD_NUM, "Restarting all worker threads with new timed INI settings.\n");
		STOP_FOR_REREAD_INI = TRUE;
		restart_waiting_workers (RESTART_ALL);
	}
}

/* Set flags so that worker threads will stop due to day/night memory */
/* changeover. */

void stop_worker_for_mem_changed (
	int	thread_num)
{
	if (WORKER_THREADS_ACTIVE && ! STOP_FOR_MEM_CHANGED[thread_num]) {
		OutputStr (thread_num, "Restarting worker with new memory settings.\n");
		MEM_FLAGS[thread_num] |= MEM_RESTARTING;
		STOP_FOR_MEM_CHANGED[thread_num] = 1;
		restart_one_waiting_worker (thread_num, RESTART_ALL);
	}
}

/* Set flag so that worker thread will stop to do priority work. */

void stop_worker_for_advanced_test (
	int	thread_num)
{
	if (WORKER_THREADS_ACTIVE && ! STOP_FOR_PRIORITY_WORK[thread_num]) {
		OutputStr (thread_num, "Restarting worker to do LL test.\n");
		STOP_FOR_PRIORITY_WORK[thread_num] = 1;
	}
}

/* Set flag so that worker thread will stop to do priority work. */

void stop_worker_for_priority_work (
	int	thread_num)
{
	if (WORKER_THREADS_ACTIVE && ! STOP_FOR_PRIORITY_WORK[thread_num]) {
		OutputStr (thread_num, "Restarting worker to do factoring prior to an LL test.\n");
		STOP_FOR_PRIORITY_WORK[thread_num] = 1;
	}
}

/* Set flags so that worker threads will stop for throttling. */

void stop_workers_for_throttle (void)
{
	if (WORKER_THREADS_ACTIVE)
		memset (STOP_FOR_THROTTLE, 1, sizeof (STOP_FOR_THROTTLE));
}

/* Set flags so that worker thread will abort processing its current */
/* work unit.  There are many reasons to do this: unreserve, factor found */
/* in another thread, server request, etc. */

void stop_worker_for_abort (
	int	thread_num)
{
	if (WORKER_THREADS_ACTIVE)
		STOP_FOR_ABORT[thread_num] = 1;
}

/* Start save files timer */

void start_save_files_timer ()
{
	add_timed_event (TE_SAVE_FILES, DISK_WRITE_TIME * 60);
}

/* Stop save files timer */

void stop_save_files_timer ()
{
	delete_timed_event (TE_SAVE_FILES);
}

/* Set flags so that worker threads will write save files */
/* at next convenient opportunity. */

void saveFilesTimer ()
{
	memset (WRITE_SAVE_FILES, 1, sizeof (WRITE_SAVE_FILES));
	start_save_files_timer ();
}

/* Return TRUE if it is time to write a save file. */

int testSaveFilesFlag (
	int	thread_num)
{
	if (WRITE_SAVE_FILES[thread_num]) {
		WRITE_SAVE_FILES[thread_num] = 0;
		return (TRUE);
	}
	return (FALSE);
}

/**************************************************************/
/*      Routines dealing with Day/Night memory settings       */
/**************************************************************/

/* Read the Memory settings from INI file */

void read_mem_info (void)
{
	const char *p;
	int	tnum;
	unsigned int seconds, seconds_until_reread;

/* Initalize the memory mutex and other memory related events */

	if (!MEM_MUTEX_INITIALIZED) {
		MEM_MUTEX_INITIALIZED = 1;
		gwmutex_init (&MEM_MUTEX);
	}

/* Lock just in case memory routines are accessing this data */

	gwmutex_lock (&MEM_MUTEX);

/* Kill the timer that triggers rereading the memory info */

	delete_timed_event (TE_MEM_CHANGE);

/* Read and parse the Memory data from the INI file */

	seconds_until_reread = 0;
	AVAIL_MEM = IniGetTimedInt (LOCALINI_FILE, "Memory", 8, &seconds);
	if (seconds && (seconds_until_reread == 0 || seconds < seconds_until_reread))
		seconds_until_reread = seconds;
	for (tnum = 0; tnum < (int) MAX_NUM_WORKER_THREADS; tnum++) {
		AVAIL_MEM_PER_WORKER[tnum] = IniGetTimedInt (LOCALINI_FILE, "Memory", AVAIL_MEM, &seconds);
		if (seconds && (seconds_until_reread == 0 || seconds < seconds_until_reread))
			seconds_until_reread = seconds;
	}

/* Compute the maximum memory setting.  If not found, assume 8MB. */

	MAX_MEM = 8;
	p = IniSectionGetStringRaw (LOCALINI_FILE, NULL, "Memory");
	if (p != NULL) for ( ; ; ) {
		unsigned long mem = atol (p);
		if (mem > MAX_MEM) MAX_MEM = mem;
		p = strstr (p, " else ");
		if (p == NULL) break;
		p = p + 6;
	}

/* Get the maximum number of workers that can use lots of memory */
/* Default is AVAIL_MEM / 200MB rounded off. */
	
	MAX_HIGH_MEM_WORKERS = IniGetTimedInt (LOCALINI_FILE, "MaxHighMemWorkers",
					       (AVAIL_MEM + 100) / 200, &seconds);
	if (seconds && (seconds_until_reread == 0 || seconds < seconds_until_reread))
		seconds_until_reread = seconds;
	if (MAX_HIGH_MEM_WORKERS < 1) MAX_HIGH_MEM_WORKERS = 1;

/* Add the event that fires when the memory settings expire. */

	if (seconds_until_reread)
		add_timed_event (TE_MEM_CHANGE, seconds_until_reread);

/* Unlock */

	gwmutex_unlock (&MEM_MUTEX);
}

/* This routine initializes mem_changed globals.  It must be called prior */
/* to launching the worker threads. */

void init_mem_state (void)
{
	int	i;

/* Clear flags saying thread is affected by changes in the memory settings. */
/* Assume each worker thread will use a default amount of memory. */

	for (i = 0; i < MAX_NUM_WORKER_THREADS; i++) {
		MEM_FLAGS[i] = MEM_USAGE_NOT_SET;
		MEM_IN_USE[i] = DEFAULT_MEM_USAGE;
		MEM_RESTART_FLAGS[i] = 0;
	}
}

/* Clear flags that keep track if the thread needs restarting */
/* on available memory changes. */

void clear_memory_restart_flags (
	int	thread_num)
{
	MEM_RESTART_FLAGS[thread_num] = 0;
}

/* Set thread to default memory usage.  For now, this is 24MB -- roughly */
/* the amount of memory used by LL test using a 2.5M FFT. */

void set_default_memory_usage (
	int	thread_num)
{
	MEM_FLAGS[thread_num] = MEM_USAGE_NOT_SET;
	MEM_IN_USE[thread_num] = DEFAULT_MEM_USAGE;

/* Clear restart flags that only apply to current work unit as opposed to */
/* most flags which are not reset until primeContinue reprocesses the worker's */
/* complete list of work. */

	MEM_RESTART_FLAGS[thread_num] &= ~MEM_RESTART_MAX_MEM_CHANGE;
	MEM_RESTART_FLAGS[thread_num] &= ~MEM_RESTART_IF_MORE;
}

/* Set flag that restarts worker if max mem changes. */
/* Needed for Pfactor so that we can compute new bounds */
/* should max mem change. */

void set_restart_if_max_memory_change (
	int	thread_num)
{
	MEM_RESTART_FLAGS[thread_num] |= MEM_RESTART_MAX_MEM_CHANGE;
}

/* Set flag that restarts worker if max mem changes. */
/* Needed for Pfactor so that we can compute new bounds */
/* should max mem change. */

void clear_restart_if_max_memory_change (
	int	thread_num)
{
	MEM_RESTART_FLAGS[thread_num] &= ~MEM_RESTART_MAX_MEM_CHANGE;
}

/* Set flag that restarts current work unit if more memory */
/* becmes available.  Used when stage 2 got far less memory than */
/* it wanted and significantly more memory whould speed up stage 2. */

void set_restart_if_more_memory_available (
	int	thread_num,
	unsigned int memory)		/* Memory needed for a restart in MB */
{
	MEM_RESTART_FLAGS[thread_num] |= MEM_RESTART_IF_MORE;
	MEM_RESTART_IF_MORE_AMOUNT[thread_num] = memory;
}

/* If the caller of avail_mem wasn't happy with the amount of memory */
/* returned, he can call this routine to set flags so that worker will be */
/* restarted when more memory becomes available. */

int avail_mem_not_sufficient (
	int	thread_num,
	unsigned long min_memory,	/* Minumum memory in MB */
	unsigned long desired_memory)	/* Desired memory in MB */
{
	OutputStr (thread_num, "Other workers are using lots of memory now.\n");
	if (MEM_RESTART_FLAGS[thread_num] & MEM_RESTART_MORE_AVAIL) {
		if (min_memory < MEM_RESTART_MIN_AMOUNT[thread_num])
			MEM_RESTART_MIN_AMOUNT[thread_num] = min_memory;
		if (desired_memory < MEM_RESTART_DESIRED_AMOUNT[thread_num])
			MEM_RESTART_DESIRED_AMOUNT[thread_num] = desired_memory;
	} else {
		MEM_RESTART_FLAGS[thread_num] |= MEM_RESTART_MORE_AVAIL;
		MEM_RESTART_MIN_AMOUNT[thread_num] = min_memory;
		MEM_RESTART_DESIRED_AMOUNT[thread_num] = desired_memory;
	}
	return (STOP_NOT_ENOUGH_MEM);
}

/* Internal routine that returns TRUE if other threads are using lots of */
/* the available memory.  We use this to delay ECM and P-1 stage 2 while other */
/* stage 2's are running. */

int are_threads_using_lots_of_memory (
	int	thread_num)
{
	int	max_high_mem, i;

/* Get the user configurable count of workers that are allowed to use */
/* lots of memory.  If this equals the number of workers (default) then */
/* there is no need to scan the workers */

	max_high_mem = MAX_HIGH_MEM_WORKERS;
	if (max_high_mem >= (int) NUM_WORKER_THREADS) return (FALSE);

/* If there are enough threads with variable memory usage, then return TRUE. */
/* To guard against an ECM stage 2 that really isn't using a whole lot of */
/* memory, also require the thread to be using 50MB. */

	for (i = 0; i < (int) NUM_WORKER_THREADS; i++)
		if (i != thread_num &&
		    (MEM_FLAGS[i] & MEM_VARIABLE_USAGE ||
		     MEM_FLAGS[i] & MEM_WILL_BE_VARIABLE_USAGE) &&
		    MEM_IN_USE[i] >= 50) {
			max_high_mem--;
			if (max_high_mem == 0) return (TRUE);
		}
	return (FALSE);
}

/* Each worker thread tells us how much memory it will be using.  This may */
/* cause other worker threads to restart if they are using more than their */
/* fair share. NOTE: This mem routine takes its argument in bytes */
/* rather than megabytes. */
/* Variable usage callers must examine the return code!  During startup */
/* all threads may not have determined their memory needs.  This routine */
/* returns TRUE if caller should recalculate the amount of memory available */
/* for use because we previously overestimated the amount of memory available */
/* to the thread. */

int set_memory_usage (
	int	thread_num,
	int	flags,		/* Valid values are MEM_VARIABLE_USAGE */
				/* and MEM_USAGE_NOT_SET */
	unsigned long memory)	/* Memory in use (in bytes!!) */
{
	int	i, best_thread, worst_thread, all_threads_set;
	unsigned long mem_usage;

/* Obtain lock before accessing memory global variables */

	gwmutex_lock (&MEM_MUTEX);

/* Set or clear flag indicating thread is executing code that can choose a */
/* different amount of memory to use. */

	if (flags & MEM_VARIABLE_USAGE)
		MEM_FLAGS[thread_num] |= MEM_VARIABLE_USAGE;
	else
		MEM_FLAGS[thread_num] &= ~MEM_VARIABLE_USAGE;
	MEM_FLAGS[thread_num] &= ~MEM_WILL_BE_VARIABLE_USAGE;

/* Set flag indicating we are guessing how much memory this thread */
/* will use because the thread has not started its work unit. */

	if (flags & MEM_USAGE_NOT_SET)
		MEM_FLAGS[thread_num] |= MEM_USAGE_NOT_SET;
	else
		MEM_FLAGS[thread_num] &= ~MEM_USAGE_NOT_SET;
	MEM_FLAGS[thread_num] &= ~MEM_RESTARTING;

/* Record the amount of memory being used */

	MEM_IN_USE[thread_num] = (memory >> 20) + 1;

/* Sum up the amount of memory used by all threads.  In case we've allocated */
/* too much memory, select a variable thread to restart.  We do this to make */
/* the thread reduce its memory usage so that the other threads will be OK. */
/* We'll restart the variable thread using the most memory. */

	mem_usage = 0;
	worst_thread = -1;
	for (i = 0; i < (int) NUM_WORKER_THREADS; i++) {
		mem_usage += MEM_IN_USE[i];
		if ((MEM_FLAGS[i] & MEM_VARIABLE_USAGE ||
		     MEM_FLAGS[i] & MEM_WILL_BE_VARIABLE_USAGE) &&
		    (worst_thread == -1 ||
		     MEM_IN_USE[i] > MEM_IN_USE[worst_thread]))
			worst_thread = i;
	}

/* If we have allocated more than the maximum allowable, then stop a */
/* thread to free up some memory.  We also make sure we are using significantly */
/* more memory than we should be so that minor fluctuations in memory */
/* usage by the fixed threads do not cause needless restarts.  The 32MB */
/* threshold is arbitrary. */

	if (mem_usage > AVAIL_MEM + 32) {

/* If the current thread is the worst thread (should only happen if there has */
/* been a wild change in other thread's memory usage between the call to */
/* avail_mem and the call to set_memory_usage), then return to caller and */
/* tell it to try again.  WARNING:  this could cause an infinite */
/* loop if caller misbehaves and tries to use the same amount of memory. */

		if (worst_thread == thread_num) {
			set_default_memory_usage (thread_num);
			gwmutex_unlock (&MEM_MUTEX);
			return (TRUE);
		}

/* If we found a worst thread and that thread has actually allocated */
/* memory (MEM_VARIABLE_USAGE), as opposed to being in the process of */
/* figuring out its memory needs (MEM_WILL_BE_VARIABLE_USAGE), then */
/* stop the offending thread. */

		if (worst_thread >= 0 && MEM_FLAGS[worst_thread] & MEM_VARIABLE_USAGE) {
			stop_worker_for_mem_changed (worst_thread);

/* Wait for the stop to take effect so that we don't briefly over-allocate memory. */

			MEM_FLAGS[thread_num] |= MEM_WAITING;
			gwmutex_unlock (&MEM_MUTEX);
			MEM_WAIT_OR_STOP_INITIALIZED[thread_num] = 1;
			gwevent_init (&MEM_WAIT_OR_STOP[thread_num]);
			gwevent_reset (&MEM_WAIT_OR_STOP[thread_num]);
			gwevent_wait (&MEM_WAIT_OR_STOP[thread_num], 20);
			gwevent_destroy (&MEM_WAIT_OR_STOP[thread_num]);
			MEM_WAIT_OR_STOP_INITIALIZED[thread_num] = 0;
			gwmutex_lock (&MEM_MUTEX);
			MEM_FLAGS[thread_num] &= ~MEM_WAITING;
		}
	}

/* See if all fixed usage threads have set their memory usage */

	all_threads_set = TRUE;
	for (i = 0; i < (int) NUM_WORKER_THREADS; i++) {
		if (MEM_FLAGS[i] & MEM_WILL_BE_VARIABLE_USAGE) continue;
		if (MEM_FLAGS[i] & MEM_USAGE_NOT_SET ||
		    MEM_FLAGS[i] & MEM_RESTARTING) {
			all_threads_set = FALSE;
			break;
		}
	}

/* If all fixed usage threads have called this routine setting their memory */
/* usage, then signal an event to wake up one of variable usage workers */
/* that is waiting for all fixed usage workers to compute their memory usage. */

	if (all_threads_set) {
		best_thread = -1;
		for (i = 0; i < (int) NUM_WORKER_THREADS; i++) {
			if (MEM_FLAGS[i] & MEM_WAITING) {
				best_thread = i;
				break;
			}
			if (MEM_FLAGS[i] & MEM_WILL_BE_VARIABLE_USAGE &&
			    (best_thread == -1 ||
			     MEM_IN_USE[i] < MEM_IN_USE[best_thread]))
				best_thread = i;
		}
		if (best_thread >= 0) {
			restart_one_waiting_worker (best_thread, RESTART_MEM_WAIT);
			all_threads_set = FALSE;
		}
	}

/* If a worker is waiting for a reduction in the number of workers */
/* using lots of memory, then check to see if it can run now. */

	if (all_threads_set && mem_usage < AVAIL_MEM - 32) {
		for (i = 0; i < (int) NUM_WORKER_THREADS; i++) {
			if (! (MEM_RESTART_FLAGS[i] & MEM_RESTART_TOO_MANY_HIGHMEM)) continue;
			if (are_threads_using_lots_of_memory (i)) continue;
			stop_worker_for_mem_changed (i);
			all_threads_set = FALSE;
			break;
		}
	}

/* If all fixed and variable usage threads have set their memory usage, */
/* then if we have enough free memory restart a work unit that could use */
/* more memory. */

	if (all_threads_set && mem_usage < AVAIL_MEM - 32) {
		best_thread = -1;
		for (i = 0; i < (int) NUM_WORKER_THREADS; i++) {
			if (MEM_RESTART_FLAGS[i] & MEM_RESTART_IF_MORE &&
			    MEM_RESTART_IF_MORE_AMOUNT[i] < AVAIL_MEM - mem_usage)
				best_thread = i;
		}
		if (best_thread >= 0) {
			stop_worker_for_mem_changed (best_thread);
			all_threads_set = FALSE;
		}
	}
		
/* If all fixed and variable usage threads have set their memory usage, */
/* then if we have enough free memory restart a thread that couldn't */
/* run a work unit due to lack of available memory. */

	if (all_threads_set && mem_usage < AVAIL_MEM - 32) {
		best_thread = -1;
		for (i = 0; i < (int) NUM_WORKER_THREADS; i++) {
			if (MEM_RESTART_FLAGS[i] & MEM_RESTART_MORE_AVAIL &&
			    MEM_RESTART_MIN_AMOUNT[i] < AVAIL_MEM - mem_usage)
				best_thread = i;
		}
		if (best_thread >= 0) {
			stop_worker_for_mem_changed (best_thread);
			all_threads_set = FALSE;
		}
	}

/* All done */

	gwmutex_unlock (&MEM_MUTEX);
	return (FALSE);
}

/* Return maximum memory (in MB) that will ever be available for a variable usage thread. */

unsigned int max_mem (void)
{
	return (MAX_MEM);
}

/* Return memory (in MB) now available for a variable usage thread. */
/* This routine takes into account the memory used by other worker threads. */
/* NOTE: caller is expected to have called are_threads_using_lots_of_memory */
/* to make sure too many workers don't become high memory users. */

int avail_mem (
	int	thread_num,
	unsigned long minimum_memory,	/* If this much memory (in bytes!) */
					/* can be returned without restarting other */
					/* workers, then do so */
	unsigned long desired_memory,	/* If this much memory (in bytes!) */
					/* can be returned without restarting other */
					/* workers, then do so */
	unsigned int *memory)		/* Returned available memory, in MB */
{
	int	i, fixed_threads[MAX_NUM_WORKER_THREADS];
	unsigned long fixed_usage, variable_usage, num_variable_threads, avail, diff;

/* Check if we are in a period of forced low memory usage */

	if (is_LowMemWhileRunning_active (thread_num)) {
		MEM_RESTART_FLAGS[thread_num] |= MEM_RESTART_LOWMEM_ENDS;
		return (STOP_NOT_ENOUGH_MEM);
	}

/* Check if we are only supposed to run high memory workers when the maximum */
/* amount memory is available. */

	if (IniGetInt (INI_FILE, "OnlyRunStage2WithMaxMemory", 0) &&
	    AVAIL_MEM != MAX_MEM) {
		OutputStr (thread_num, "Waiting for maximum available memory to run stage 2.\n");
		MEM_RESTART_FLAGS[thread_num] |= MEM_RESTART_MAX_MEM_AVAILABLE;
		return (STOP_NOT_ENOUGH_MEM);
	}

/* Check if we must wait for more memory to become available.  This */
/* happens when we reach the maximum allowable number of threads using a lot */
/* of memory. */

	if (are_threads_using_lots_of_memory (thread_num)) {
		OutputStr (thread_num, "Exceeded limit on number of workers that can use lots of memory.\n");
		MEM_RESTART_FLAGS[thread_num] |= MEM_RESTART_TOO_MANY_HIGHMEM;
		return (STOP_NOT_ENOUGH_MEM);
	}

/* Obtain lock before accessing memory global variables */

	gwmutex_lock (&MEM_MUTEX);

/* Set flag saying this will be a variable usage thread.  Remember the */
/* "good enough" value as it will be helpful in determining the best */
/* value this routine should return (for this thread and other threads) */

	MEM_FLAGS[thread_num] |= MEM_WILL_BE_VARIABLE_USAGE;
	desired_memory = (desired_memory >> 20) + 1;
	MEM_IN_USE[thread_num] = desired_memory;

/* If any workers have not yet set their memory usage, then wait for them */
/* to do so.  This allows us to accurately guage how much fixed memory */
/* is consumed and how many variable usage workers there are. */
/* Just in case we wake up from the timeout (should rarely happen), we try*/
/* to stagger the timeouts by adding the thread number. */

	for (i = 0; i < (int) NUM_WORKER_THREADS; i++) {
		if (i == thread_num) continue;
		if (MEM_FLAGS[i] & MEM_USAGE_NOT_SET ||
		    MEM_FLAGS[i] & MEM_RESTARTING) {
			gwmutex_unlock (&MEM_MUTEX);
			MEM_WAIT_OR_STOP_INITIALIZED[thread_num] = 1;
			gwevent_init (&MEM_WAIT_OR_STOP[thread_num]);
			gwevent_reset (&MEM_WAIT_OR_STOP[thread_num]);
			gwevent_wait (&MEM_WAIT_OR_STOP[thread_num], 20 + thread_num);
			gwevent_destroy (&MEM_WAIT_OR_STOP[thread_num]);
			MEM_WAIT_OR_STOP_INITIALIZED[thread_num] = 0;
			gwmutex_lock (&MEM_MUTEX);
		}
	}

/* Sum up the amount of memory used by threads that cannot adjust their */
/* memory usage.  Also count how many threads (including this one) can */
/* adjust their memory usage. */

	fixed_usage = 0;
	variable_usage = 0;
	num_variable_threads = 0;
	for (i = 0; i < (int) NUM_WORKER_THREADS; i++) {
		if (MEM_FLAGS[i] & MEM_VARIABLE_USAGE ||
		    MEM_FLAGS[i] & MEM_WILL_BE_VARIABLE_USAGE) {
			num_variable_threads++;
			variable_usage += MEM_IN_USE[i];
			fixed_threads[i] = FALSE;
		} else {
			fixed_usage += MEM_IN_USE[i];
			fixed_threads[i] = TRUE;
		}
	}

/* We can now calculate how much memory is available for the threads */
/* that are using a variable amount of memory.  */

	avail = (AVAIL_MEM > fixed_usage) ? AVAIL_MEM - fixed_usage : 0;
	for (i = 0; i < (int) NUM_WORKER_THREADS; i++) {
		unsigned long avail_per_worker;
		if (i == thread_num) continue;
		if (fixed_threads[i]) continue;
		avail_per_worker = avail / num_variable_threads;

/* If any variable threads are either using less memory than they are */
/* allowed to use or all variable threads can fit in available memory, */
/* then treat this worker like a fixed memory user. */

		if (MEM_FLAGS[i] & MEM_VARIABLE_USAGE &&
		    (MEM_IN_USE[i] < avail_per_worker ||
		     fixed_usage + variable_usage <= AVAIL_MEM)) {
			avail -= MEM_IN_USE[i];
			fixed_threads[i] = TRUE;
			num_variable_threads--;
			i = -1;	    /* Restart loop */
			continue;
		}

/* If any variable thread is prohibited from using its full share */
/* of the remaining available pool, then distribute the excess among */
/* the other variable usage threads. */

		if (MEM_FLAGS[i] & MEM_WILL_BE_VARIABLE_USAGE &&
		    AVAIL_MEM_PER_WORKER[i] < avail_per_worker) {
			avail -= AVAIL_MEM_PER_WORKER[i];
			fixed_threads[i] = TRUE;
			num_variable_threads--;
			i = -1;	    /* Restart loop */
			continue;
		}
	}
	avail = avail / num_variable_threads;

/* Free lock after accessing memory global variables */

	gwmutex_unlock (&MEM_MUTEX);

/* If there weren't enough memory available, try again later */

	minimum_memory = (minimum_memory >> 20) + 1;
	if (avail < minimum_memory)
		return (avail_mem_not_sufficient (thread_num, minimum_memory, desired_memory));

/* Return the amount of memory this thread can use.  If all variable */
/* threads can obtain their desired memory, then distribute the excess */
/* among all the variable threads.  Otherwise, return my pro-rata share */
/* of variable memory, any overcommitted workers will be restarted once this */
/* thread calls set_memory_usage letting us know how much of the available */
/* memory it actually used. */

	if (fixed_usage + variable_usage <= AVAIL_MEM)
		*memory = desired_memory +
			  (AVAIL_MEM - (fixed_usage + variable_usage)) / num_variable_threads;
	else
		*memory = avail;

/* If memory exceeds this worker's maximum, then only return */
/* this worker's maximum. */

	if (*memory > AVAIL_MEM_PER_WORKER[thread_num])
		*memory = AVAIL_MEM_PER_WORKER[thread_num];

/* As a first approximation, mark the work unit as available for restart */
/* if more memory is available whenever we are near minimum_memory.  The caller */
/* can override our guess if he so desires */	

	diff = desired_memory - minimum_memory;
	if (*memory <= minimum_memory + diff / 4)
		set_restart_if_more_memory_available (thread_num, diff / 4);
	else if (*memory <= minimum_memory + diff / 2)
		set_restart_if_more_memory_available (thread_num, diff / 2);

/* Return clean stop code */

	return (0);
}

/* Routine to notify all worker threads the day/night memory settings */
/* have changed.  This is called when the memory change timer fires OR */
/* when memory settings are changed by the GUI. */

void mem_settings_have_changed (void)
{
	unsigned int old_avail_mem, old_max_mem;
	int	tnum;

/* Recompute the available memory and restart the memory changed timer */

	old_avail_mem = AVAIL_MEM;
	old_max_mem = MAX_MEM;
	read_mem_info ();

/* If the worker threads are not active then no workers need restarting */

	if (! WORKER_THREADS_ACTIVE) return;

/* If maximum memory has changed see which threads need restarting. */
/* Those threads that are in stage 1 of pfactor work will want to compute */
/* new bounds. */

	if (MAX_MEM != old_max_mem)
		for (tnum = 0; tnum < (int) NUM_WORKER_THREADS; tnum++)
			if (MEM_RESTART_FLAGS[tnum] & MEM_RESTART_MAX_MEM_CHANGE)
				stop_worker_for_mem_changed (tnum);

/* If available memory is now equal to maximum memory see which threads */
/* need restarting. Those threads that postponed work because they only */
/* run during memory need restarting. */

	if (AVAIL_MEM == MAX_MEM)
		for (tnum = 0; tnum < (int) NUM_WORKER_THREADS; tnum++)
			if (MEM_RESTART_FLAGS[tnum] & MEM_RESTART_MAX_MEM_AVAILABLE)
				stop_worker_for_mem_changed (tnum);

/* If available memory has increased we may pick a thread to restart. */
/* Those threads that postponed work because there wasn't enough memory */
/* need restarting. */

	if (AVAIL_MEM > old_avail_mem)
		for (tnum = 0; tnum < (int) NUM_WORKER_THREADS; tnum++)
			if (MEM_RESTART_FLAGS[tnum] & MEM_RESTART_MORE_AVAIL)
				stop_worker_for_mem_changed (tnum);

/* If any worker now exceeds (by 10MB) the per-worker maximum, then restart. */

	for (tnum = 0; tnum < (int) NUM_WORKER_THREADS; tnum++)
		if (MEM_FLAGS[tnum] & MEM_VARIABLE_USAGE &&
		    MEM_IN_USE[tnum] > AVAIL_MEM_PER_WORKER[tnum] + 10)
			stop_worker_for_mem_changed (tnum);

/* If available memory has decreased we may pick a thread to restart. */
/* If total memory in use is greater than the new available, then pick */
/* one of the variable threads to restart.  Note that if any threads */
/* haven't yet set their memory usage, then when they do set their memory */
/* usage this overcommitment will be sorted out then. */

	if (AVAIL_MEM < old_avail_mem) {
		unsigned long mem_usage;
		int	worst_thread;

		mem_usage = 0;
		worst_thread = -1;
		for (tnum = 0; tnum < (int) NUM_WORKER_THREADS; tnum++) {
			mem_usage += MEM_IN_USE[tnum];
			if (MEM_FLAGS[tnum] & MEM_USAGE_NOT_SET ||
			    MEM_FLAGS[tnum] & MEM_RESTARTING ||
			    MEM_FLAGS[tnum] & MEM_WILL_BE_VARIABLE_USAGE) {
				worst_thread = -1;
				break;
			}
			if (MEM_FLAGS[tnum] & MEM_VARIABLE_USAGE &&
			    (worst_thread == -1 ||
			     MEM_IN_USE[tnum] > MEM_IN_USE[worst_thread]))
				worst_thread = tnum;
		}
		if (mem_usage > AVAIL_MEM + 32 && worst_thread != -1)
			stop_worker_for_mem_changed (worst_thread);
	}
}

/* Routine to force any workers that are using lots of memory to stop */
/* and restart.  This happens when LowMemWhileRunnning is activated. */

void stop_high_memory_workers (void)
{
	int	i;

/* If the worker threads are not active then no workers need restarting */

	if (! WORKER_THREADS_ACTIVE) return;

/* Obtain lock before accessing memory global variables */

	gwmutex_lock (&MEM_MUTEX);

/* Look for workers marked with variable usage */

	for (i = 0; i < (int) NUM_WORKER_THREADS; i++) {
		if (MEM_FLAGS[i] & MEM_VARIABLE_USAGE ||
		    MEM_FLAGS[i] & MEM_WILL_BE_VARIABLE_USAGE)
			stop_worker_for_mem_changed (i);
	}

/* All done */

	gwmutex_unlock (&MEM_MUTEX);
}

/* Routine to restart workers that were stopped due to LowMemWhileRunning */

void restart_high_memory_workers (void)
{
	int	tnum;

/* If the worker threads are not active then no workers need restarting */

	if (! WORKER_THREADS_ACTIVE) return;

/* Restart the workers affected by LowMemWhileRunning */

	for (tnum = 0; tnum < (int) NUM_WORKER_THREADS; tnum++)
		if (MEM_RESTART_FLAGS[tnum] & MEM_RESTART_LOWMEM_ENDS) {
			MEM_RESTART_FLAGS[tnum] &= ~MEM_RESTART_LOWMEM_ENDS;
			stop_worker_for_mem_changed (tnum);
		}
}


/**************************************************************/
/*           Routines dealing running on battery              */
/**************************************************************/

void start_battery_timer (void)
{
	if (RUN_ON_BATTERY) return;
	add_timed_event (TE_BATTERY_CHECK, TE_BATTERY_CHECK_FREQ);
}

void stop_battery_timer (void)
{
	delete_timed_event (TE_BATTERY_CHECK);
}

/* This routine is called if the user changes the RUN_ON_BATTERY setting */
/* from the GUI or it is changed by talking to the server. */

void run_on_battery_changed (void)
{
	if (WORKER_THREADS_ACTIVE) {
		stop_battery_timer ();
		start_battery_timer ();
	}
}

void test_battery (void)
{
	if (OnBattery ()) STOP_FOR_BATTERY = TRUE;
	else if (STOP_FOR_BATTERY) {
		STOP_FOR_BATTERY = FALSE;
		restart_waiting_workers (RESTART_BATTERY);
	}
}

/* Stopping while on battery power, restart thread only when AC power */
/* is restored. */

void implement_stop_battery (
	int	thread_num)
{

/* Output message, change title and icon */

	title (thread_num, "Battery Pause");
	OutputStr (thread_num, "Worker stopped while on battery power.\n");
	ChangeIcon (thread_num, IDLE_ICON);

/* Wait for AC power */

	OFF_BATTERY_OR_STOP_INITIALIZED[thread_num] = 1;
	gwevent_init (&OFF_BATTERY_OR_STOP[thread_num]);
	gwevent_reset (&OFF_BATTERY_OR_STOP[thread_num]);
	gwevent_wait (&OFF_BATTERY_OR_STOP[thread_num], 0);
	gwevent_destroy (&OFF_BATTERY_OR_STOP[thread_num]);
	OFF_BATTERY_OR_STOP_INITIALIZED[thread_num] = 0;

/* Output message, change title and icon */

	title (thread_num, "Working");
	OutputStr (thread_num, "AC power restored, restarting worker.\n");
	ChangeIcon (thread_num, WORKING_ICON);
}

/**************************************************************/
/*           Routines dealing with priority work              */
/**************************************************************/

void start_priority_work_timer (void)
{
	if (IniGetInt (INI_FILE, "SequentialWorkToDo", 1)) return;
	add_timed_event (TE_PRIORITY_WORK, TE_PRIORITY_WORK_FREQ);
}

void stop_priority_work_timer (void)
{
	delete_timed_event (TE_PRIORITY_WORK);
}

/* Returns true if this is a priority work item */

int isPriorityWork (
	struct work_unit *w)
{
	if (isWorkUnitActive (w)) return (FALSE);
	if (w->work_type == WORK_ADVANCEDTEST) return (TRUE);
	if ((w->work_type == WORK_TEST || w->work_type == WORK_DBLCHK) &&
	    (w->sieve_depth < w->factor_to || !w->pminus1ed))
		return (TRUE);
	if (w->work_type == WORK_PRP && !w->pminus1ed)
		return (TRUE);
	return (FALSE);
}

/* For all threads, check if any of the Lucas-Lehmer test lines also */
/* require factoring. This will force factoring to be done first - giving */
/* us more accurate estimates of how much work is queued up. */

void check_for_priority_work (void)
{
	int	tnum;
	struct work_unit *w;

	for (tnum = 0; tnum < (int) NUM_WORKER_THREADS; tnum++) {
		w = NULL;
		for ( ; ; ) {
			w = getNextWorkToDoLine (tnum, w, SHORT_TERM_USE);
			if (w == NULL) break;
			if (isPriorityWork (w))
				stop_worker_for_priority_work (tnum);
		}
       }
}

/**************************************************************/
/*     Routines dealing with stopping specific workers        */
/**************************************************************/

void mark_workers_active (
	int	thread_num)	/* Number of workers to mark active */
				/* or (<= 0) the only worker to mark */
{
	int	i;

	memset (ACTIVE_WORKERS, 0, sizeof (ACTIVE_WORKERS));
	for (i = 0; i < thread_num; i++) ACTIVE_WORKERS[i] = 1;
	if (thread_num <= 0) ACTIVE_WORKERS[-thread_num] = 1;
}

void start_one_worker (
	int	thread_num)
{
	if (thread_num < 0 || thread_num > (int) NUM_WORKER_THREADS) {
		OutputStr (MAIN_THREAD_NUM, "Worker thread number out of range.\n");
		return;
	}
	if (ACTIVE_WORKERS[thread_num]) {
		OutputStr (MAIN_THREAD_NUM, "Worker thread is already running.\n");
		return;
	}
	ACTIVE_WORKERS[thread_num] = 1;

	// Restart the worker
	restart_one_waiting_worker (thread_num, RESTART_USER_START);
}

void stop_one_worker (
	int	thread_num)
{
	if (thread_num < 0 || thread_num > (int) NUM_WORKER_THREADS) {
		OutputStr (MAIN_THREAD_NUM, "Worker thread number out of range.\n");
		return;
	}
	if (!ACTIVE_WORKERS[thread_num]) {
		OutputStr (MAIN_THREAD_NUM, "Worker thread is already stopped.\n");
		return;
	}
	ACTIVE_WORKERS[thread_num] = 0;
}

void implement_stop_one_worker (
	int	thread_num)
{

/* If some race condition has caused the worker active flag */
/* to be set, then do not wait for an event. */

	if (ACTIVE_WORKERS[thread_num]) return;

/* Output a worker stopping message and change the icon */

//// bug - this message will be output even if worker never started
	OutputStr (thread_num, "Stopping worker at user request.\n");
	ChangeIcon (thread_num, IDLE_ICON);	/* Idle icon while stopped */

/* Set memory usage to zero */

	set_memory_usage (thread_num, 0, 0);

/* Initialize and then wait for the event */

	USER_START_OR_STOP_INITIALIZED[thread_num] = 1;
	gwevent_init (&USER_START_OR_STOP[thread_num]);
	gwevent_reset (&USER_START_OR_STOP[thread_num]);
	gwevent_wait (&USER_START_OR_STOP[thread_num], 0);

/* Destroy the event */

	USER_START_OR_STOP_INITIALIZED[thread_num] = 0;
	gwevent_destroy (&USER_START_OR_STOP[thread_num]);

/* Output a worker starting message and change the icon */

////bug - why output a restart message if we are only resuming to exit?
	OutputStr (thread_num, "Resuming worker at user request.\n");
	ChangeIcon (thread_num, WORKING_ICON);
}

unsigned int active_workers_count (void)
{
	unsigned int i, count;

	for (i = count = 0; i < WORKER_THREADS_ACTIVE; i++)
		if (ACTIVE_WORKERS[i]) count++;
	return (count);
}

/**************************************************************/
/*        Routines dealing with "pause while running"         */
/**************************************************************/

/* Internal routine to parse a PauseWhileRunning or LowMemWhileRunning entry */

void parse_pause_info (
       char	*buf,		/* Comma separated list of program names */
       int	thread_num,	/* Worker thread to pause */
       int	low_mem)	/* Flag for LowMemWhileRunning */
{
	struct pause_info *data;
	char	*p, *bracket, *comma;

	if (*buf == 0) return;

	for (p = buf; ; p = comma + 1) {
		comma = strchr (p, ',');
		if (comma != NULL) *comma = 0;

		data = (struct pause_info *) malloc (sizeof (struct pause_info));
		if (data == NULL) return;
		data->next = PAUSE_DATA;
		PAUSE_DATA = data;

		data->thread_num = thread_num;
		data->low_mem = low_mem;
		bracket = strchr (p, '[');
		if (bracket != NULL) {
			*bracket = 0;
			data->workers_affected = atoi (bracket+1);
		} else
			data->workers_affected = MAX_NUM_WORKER_THREADS;

		if (!low_mem && *p == '*')
			data->program_name = NULL;
		else {
			data->program_name = (char *) malloc (strlen (p) + 1);
			if (data->program_name == NULL) return;
			strupper (p);
			strcpy (data->program_name, p);
		}

		if (comma == NULL) break;
	}
}

/* Read the PauseWhileRunning and LowMemWhileRunning settings */

void read_pause_info (void)
{
	int	tnum;
	char	buf[250];
	unsigned int seconds, seconds_until_reread;

/* Initalize the mutex */

	if (!PAUSE_MUTEX_INITIALIZED) {
		PAUSE_MUTEX_INITIALIZED = 1;
		gwmutex_init (&PAUSE_MUTEX);
	}

/* Lock just in case implement_pause is accessing this data */

	gwmutex_lock (&PAUSE_MUTEX);

/* Kill the timer that triggers rereading the pause info */

	delete_timed_event (TE_READ_PAUSE_DATA);

/* Free the previous pause data */

	while (PAUSE_DATA != NULL) {
		struct pause_info *p;
		p = PAUSE_DATA;
		PAUSE_DATA = p->next;
		if (p->program_name != NULL) free (p->program_name);
		free (p);
	}

/* Read and parse the PauseWhileRunning data from the ini file */

	seconds_until_reread = 0;
	IniGetTimedString (INI_FILE, "PauseWhileRunning", buf, sizeof (buf), NULL, &seconds);
	if (seconds && (seconds_until_reread == 0 || seconds < seconds_until_reread))
		seconds_until_reread = seconds;
	parse_pause_info (buf, ALL_WORKERS, FALSE);
	for (tnum = 0; tnum < (int) NUM_WORKER_THREADS; tnum++) {
		char	section_name[32];
		sprintf (section_name, "Worker #%d", tnum+1);
		IniSectionGetTimedString (INI_FILE, section_name, "PauseWhileRunning", buf, sizeof (buf), NULL, &seconds);
		if (seconds && (seconds_until_reread == 0 || seconds < seconds_until_reread))
			seconds_until_reread = seconds;
		parse_pause_info (buf, tnum, FALSE);
	}
	PAUSE_WHILE_RUNNING_FREQ = IniGetTimedInt (INI_FILE, "PauseCheckInterval", 10, &seconds);
	if (seconds && (seconds_until_reread == 0 || seconds < seconds_until_reread))
		seconds_until_reread = seconds;

/* Also read in the LowMemWhileRunning program list */

	IniGetTimedString (INI_FILE, "LowMemWhileRunning", buf, sizeof (buf), NULL, &seconds);
	if (seconds && (seconds_until_reread == 0 || seconds < seconds_until_reread))
		seconds_until_reread = seconds;
	parse_pause_info (buf, ALL_WORKERS, TRUE);

/* Add the event that fires when the memory settings expire. */

	if (seconds_until_reread)
		add_timed_event (TE_READ_PAUSE_DATA, seconds_until_reread);

/* If the pause timer is active, then call checkPauseWhileRunning so that */
/* we can decide which workers need to be paused based on this new pause info. */

	if (is_timed_event_active (TE_PAUSE_WHILE)) {
		delete_timed_event (TE_PAUSE_WHILE);
		checkPauseWhileRunning ();
	}

/* Unlock */

	gwmutex_unlock (&PAUSE_MUTEX);
}

void start_pause_while_running_timer (void)
{
	if (PAUSE_DATA == NULL) return;
	add_timed_event (TE_PAUSE_WHILE, PAUSE_WHILE_RUNNING_FREQ);
}

void stop_pause_while_running_timer (void)
{
	delete_timed_event (TE_PAUSE_WHILE);
}

/* Internal routine to pick the "best" worker to pause */

int best_pause_candidate (
	struct pause_info **workers_to_pause)
{
	int	i;

/* Loop through all the workers.  Give preference to any worker that */
/* is paused waiting for work or is already paused. */

	for (i = 0; i < (int) NUM_WORKER_THREADS; i++) {
		if (!ACTIVE_WORKERS[i]) continue;
		if (workers_to_pause[i] != NULL) continue;
		if (WORK_AVAILABLE_OR_STOP_INITIALIZED[i]) return (i);
		if (STOP_FOR_PAUSE[i] != NULL) return (i);
	}

/* Loop through all the workers.  Give preference to any worker that */
/* hasn't gotten started yet or is in low mem state but would rather */
/* be doing high mem work.  Note the MEM_RESTART_MAX_MEM_CHANGE flag */
/* is not checked because that is the flag that recomputes pfactor */
/* bounds on change in max mem. */

	for (i = 0; i < (int) NUM_WORKER_THREADS; i++) {
		if (!ACTIVE_WORKERS[i]) continue;
		if (workers_to_pause[i] != NULL) continue;
		if (MEM_FLAGS[i] & MEM_USAGE_NOT_SET) return (i);
		if (MEM_FLAGS[i] & MEM_RESTARTING) return (i);
		if (MEM_RESTART_FLAGS[i] & ~MEM_RESTART_MAX_MEM_CHANGE) return (i);
	}

/* Loop through all the workers.  Return first one we find. */

	for (i = 0; i < (int) NUM_WORKER_THREADS; i++) {
		if (!ACTIVE_WORKERS[i]) continue;
		if (workers_to_pause[i] != NULL) continue;
		return (i);
	}

/* Return 0 if we've paused all the workers */

	return (0);
}

/* Every time the pause-while-running timer fires, this routine is called */

void checkPauseWhileRunning (void)
{
	struct pause_info *p, *lowmem;
	struct pause_info *workers_to_pause[MAX_NUM_WORKER_THREADS];
	int	i, named_program_entries;

/* Clear flag indicating a running program matched a pause_info entry */

	named_program_entries = FALSE;
	for (p = PAUSE_DATA; p != NULL; p = p->next) {
		p->matching_program[0] = 0;
		if (p->program_name != NULL) named_program_entries = TRUE;
	}

/* Call OS-specific routine to see if a process is running such that */
/* we should pause.  This OS-specific routine must get the list of active */
/* processes and call isInPauseList for each one. */

	checkPauseListCallback ();

/* Examine pause info entries to see if a period of forced low memory usage */
/* should be in effect. */

	lowmem = NULL;
	for (p = PAUSE_DATA; p != NULL; p = p->next) {
		if (!p->low_mem) continue;
		if (p->matching_program[0]) lowmem = p;
	}
	p = STOP_FOR_LOW_MEMORY;
	STOP_FOR_LOW_MEMORY = lowmem;
	if (p == NULL && STOP_FOR_LOW_MEMORY != NULL) {
		char	buf[150];
		sprintf (buf, "Entering a period of low memory usage because %s is running.\n",
			 lowmem->matching_program);
		OutputStr (MAIN_THREAD_NUM, buf);
		stop_high_memory_workers ();
	}
	if (p != NULL && STOP_FOR_LOW_MEMORY == NULL) {
		restart_high_memory_workers ();
	}

/* Examine pause info entries to see which ones matched.  In this pass */
/* we are looking for pause_info entries that pause a specific worker. */

	memset (workers_to_pause, 0, sizeof (workers_to_pause));
	for (p = PAUSE_DATA; p != NULL; p = p->next) {
		if (p->thread_num == ALL_WORKERS) continue;
		if (p->low_mem) continue;
		if (p->program_name == NULL || p->matching_program[0])
			workers_to_pause[p->thread_num] = p;
	}

/* Examine pause info entries to see which ones matched.  In this pass */
/* we are looking for pause_info entries that let us choose which worker */
/* we want to pause.  Choose the "best" worker to pause. */

	for (p = PAUSE_DATA; p != NULL; p = p->next) {
		if (p->low_mem) continue;
		if (p->program_name == NULL || p->matching_program[0]) {
			int	count;
			count = p->workers_affected;
			if (p->thread_num != ALL_WORKERS) count--;
			for (i = 0; i < count; i++)
				workers_to_pause[best_pause_candidate (workers_to_pause)] = p;
		}
	}

/* We have now determined which workers we want to pause.  Compare that */
/* to the workers that are currently paused.  Pause more workers or */
/* resume workers as appropriate. */

	for (i = 0; i < (int) NUM_WORKER_THREADS; i++) {
		p = STOP_FOR_PAUSE[i];
		STOP_FOR_PAUSE[i] = workers_to_pause[i];
		if (p != NULL && STOP_FOR_PAUSE[i] == NULL)
			restart_one_waiting_worker (i, RESTART_END_PAUSE);
	}

/* If there are any pause-for-specific program entries, then we must reset */
/* the timer to check the pause list in a few seconds.  If there are only */
/* star (match any program) entries, then we don't need to check the pause */
/* list until new PauseWhileRunning info is read from the INI file.  However, */
/* do not delete the timer as read_pause_info checks for the timer being active! */	

	add_timed_event (TE_PAUSE_WHILE, named_program_entries ? PAUSE_WHILE_RUNNING_FREQ : 1000000);
}

/* This routine is called by the OS-specific routine that gets the process */
/* list.  It returns TRUE if an active process is in the pause-while-running */
/* list. */

void isInPauseList (
	char	*program_name)
{
	struct pause_info *p;
	char	buf[512];

	strcpy (buf, program_name);
	strupper (buf);
	for (p = PAUSE_DATA; p != NULL; p = p->next) {
		if (p->program_name != NULL &&
		    strstr (buf, p->program_name) != NULL) {
			buf[sizeof(p->matching_program)-1] = 0;
			strcpy (p->matching_program, buf);
		}
	}
}

/* This routine implements a pause for one worker thread */

void implement_pause (
	int	thread_num)
{
	struct pause_info *p;
	char	buf[140];

/* Lock so that read_pause_info cannot free our structure while we are accessing it */

	gwmutex_lock (&PAUSE_MUTEX);

/* Get the pause_info struct that is causing us to pause. */
/* Return quickly if the pause has already been cancelled. */	

	p = STOP_FOR_PAUSE[thread_num];
	if (p == NULL) {
		gwmutex_unlock (&PAUSE_MUTEX);
		return;
	}

/* Output an informative message.  If we are in a sleep time period */
/* (a "*" PauseWhileRunning entry) then output a different message than */
/* if we are pausing because a specific program is running. */

	if (p->program_name == NULL) {
		time_t	sleep_time;
		char	*time_as_string;

		sleep_time = timed_event_fire_time (TE_READ_PAUSE_DATA);
		time_as_string = sleep_time ? ctime (&sleep_time) : "forever";
		if (NUM_WORKER_THREADS == 1)
			sprintf (buf, "Sleeping until %s\n", time_as_string);
		else if (p->workers_affected == 1)
			sprintf (buf, "Sleeping one worker until %s\n", time_as_string);
		else
			sprintf (buf, "Sleeping %d workers until %s\n", p->workers_affected, time_as_string);
		OutputStr (thread_num, buf);
		title (thread_num, "Sleeping");
	} else {
		sprintf (buf, "Pausing because %s is running.\n", p->matching_program);
		OutputStr (thread_num, buf);
		title (thread_num, "Paused");
	}

/* Unlock */

	gwmutex_unlock (&PAUSE_MUTEX);

/* Wait for the end of the pause */

	END_PAUSE_OR_STOP_INITIALIZED[thread_num] = 1;
	gwevent_init (&END_PAUSE_OR_STOP[thread_num]);
	gwevent_reset (&END_PAUSE_OR_STOP[thread_num]);
	gwevent_wait (&END_PAUSE_OR_STOP[thread_num], 0);
	gwevent_destroy (&END_PAUSE_OR_STOP[thread_num]);
	END_PAUSE_OR_STOP_INITIALIZED[thread_num] = 0;

/* Output another informative message */

	OutputStr (thread_num, "Resuming processing.\n");
	title (thread_num, "Resuming");
}

/* This routine checks for a forced low memory situation */

int is_LowMemWhileRunning_active (
	int	thread_num)
{
	struct pause_info *p;
	char	buf[140];

/* Lock so that read_pause_info cannot free our structure while we are accessing it */

	gwmutex_lock (&PAUSE_MUTEX);

/* Get the pause_info struct that is causing the low memory situation. */
/* Return quickly if not in a low memory situation. */

	p = STOP_FOR_LOW_MEMORY;
	if (p == NULL) {
		gwmutex_unlock (&PAUSE_MUTEX);
		return (FALSE);
	}

/* Output an informative message.  If we are in a sleep time period */
/* (a "*" PauseWhileRunning entry) then output a different message than */
/* if we are pausing because a specific program is running. */

	sprintf (buf, "Cannot use lots of memory because %s is running.\n", p->matching_program);
	OutputStr (thread_num, buf);

/* Unlock */

	gwmutex_unlock (&PAUSE_MUTEX);
	return (TRUE);
}

/**************************************************************/
/*             Routines dealing with throttling               */
/**************************************************************/

int	THROTTLE_SLEEP_TIME_IN_SEC = 0;
int	THROTTLE_SLEEP_TIME_IN_MS = 0;

void start_throttle_timer (void)
{
	if (THROTTLE_PCT <= 0 || THROTTLE_PCT >= 100) return;
	THROTTLE_SLEEP_TIME_IN_SEC = (int)
		((double) TE_THROTTLE_FREQ * (100.0 / (double) THROTTLE_PCT - 1.0));
	THROTTLE_SLEEP_TIME_IN_MS = (int)
		((double) (TE_THROTTLE_FREQ + THROTTLE_SLEEP_TIME_IN_SEC) *
			  (100.0 - (double) THROTTLE_PCT) * 10.0);
	add_timed_event (TE_THROTTLE, TE_THROTTLE_FREQ);
}

void stop_throttle_timer (void)
{
	delete_timed_event (TE_THROTTLE);
}

/* Every time the throttle timer fires, this routine is called */

int handleThrottleTimerEvent (void)
{
	stop_workers_for_throttle ();

/* Assume most threads will pause very soon.  Set timer to fire again */
/* after the idle time plus TE_THROTTLE_FREQ time.  That way each thread */
/* will run approximately TE_THROTTLE_FREQ seconds and idle */
/* THROTTLE_SLEEP_TIME_IN_SEC seconds for a CPU usage of THROTTLE_PCT. */

	return (TE_THROTTLE_FREQ + THROTTLE_SLEEP_TIME_IN_SEC);
}

/* This routine implements a throttle for one worker thread */

void implementThrottle (
	int	thread_num)
{
	int	totaltime;

/* Every 0.1 seconds see if we should resume processing.  We check */
/* frequently so that we can be responsive to an ESC or terminate command. */

	for (totaltime = 0; totaltime < THROTTLE_SLEEP_TIME_IN_MS; totaltime += 100) {
		if (WORKER_THREADS_STOPPING) return;
		Sleep (100);
	}
}

/**************************************************************/
/*                     Utility Routines                       */
/**************************************************************/

/* Return true is exponent yields a known Mersenne prime */

int isKnownMersennePrime (
	unsigned long p)
{
	return (p == 2 || p == 3 || p == 5 || p == 7 || p == 13 || p == 17 ||
		p == 19 || p == 31 || p == 61 || p == 89 || p == 107 ||
		p == 127 || p == 521 || p == 607 || p == 1279 || p == 2203 ||
		p == 2281 || p == 3217 || p == 4253 || p == 4423 ||
		p == 9689 || p == 9941 || p == 11213 || p == 19937 ||
		p == 21701 || p == 23209 || p == 44497 || p == 86243 ||
		p == 110503 || p == 132049 || p == 216091 || p == 756839 ||
		p == 859433 || p == 1257787 || p == 1398269 || p == 2976221 ||
		p == 3021377 || p == 6972593 || p == 13466917 ||
		p == 20996011 || p == 24036583 || p == 25964951 ||
		p == 30402457 || p == 32582657 || p == 37156667 || p == 43112609);
}

/* Make a string out of a 96-bit value (a found factor) */

void makestr (
	unsigned long hsw,
	unsigned long msw,
	unsigned long lsw,
	char	*buf)			/* An 80 character output buffer */
{
	int	i, j, k, carry;
	unsigned long x[3];
	char	pow[80];

	x[0] = hsw; x[1] = msw; x[2] = lsw;
	for (i = 0; i < 79; i++) pow[i] = '0', buf[i] = '0';
	pow[78] = '1';
	pow[79] = buf[79] = 0;

	for (i = 3; i--; ) {
		for (j = 0; j < 32; j++) {
			if (x[i] & 1) {
				carry = 0;
				for (k = 79; k--; ) {
					buf[k] = buf[k] - '0' +
						pow[k] - '0' + carry;
					carry = buf[k] / 10;
					buf[k] %= 10;
					buf[k] += '0';
				}
			}
			carry = 0;
			for (k = 79; k--; ) {
				pow[k] = (pow[k] - '0') * 2 + carry;
				carry = pow[k] / 10;
				pow[k] %= 10;
				pow[k] += '0';
			}
			x[i] >>= 1;
		}
	}
	while (buf[0] == '0') strcpy (buf, buf+1);
}

/* Sleep five minutes before restarting */

int SleepFive (
	int	thread_num)
{
	int	i;

	OutputStr (thread_num, ERRMSG4);
	BlinkIcon (thread_num, 10);		/* Blink icon for 10 seconds */
	for (i = 0; i < 100; i++) {
		Sleep (100);
		if (WORKER_THREADS_STOPPING) return (STOP_ESCAPE);
	}
	ChangeIcon (thread_num, IDLE_ICON);	/* Idle icon while stopped */
	for (i = 0; i < 2900; i++) {
		Sleep (100);
		if (WORKER_THREADS_STOPPING) return (STOP_ESCAPE);
	}
	ChangeIcon (thread_num, WORKING_ICON);	/* Back to the working icon */
	return (0);
}

/* Truncate a percentage to the requested number of digits. */
/* Truncating prevents 99.5% from showing up as 100% complete. */

double trunc_percent (
	double	percent)
{
	percent *= 100.0;
	if (percent > 100.0) percent = 100.0;
	percent -= 0.5 * pow (10.0, - (double) PRECISION);
	if (percent < 0.0) return (0.0);
	return (percent);
}

/*************************************/
/* Routines used to time code chunks */
/*************************************/

void clear_timers (double *timers, int num_timers) {
	int	i;
	for (i = 0; i < num_timers; i++) timers[i] = 0.0;
}

void clear_timer (
	double	*timers,
	int	i)
{
	timers[i] = 0.0;
}

void start_timer (
	double	*timers,
	int	i)
{
	if (RDTSC_TIMING < 10) {
		timers[i] -= getHighResTimer ();
	} else if (RDTSC_TIMING > 10 && (CPU_FLAGS & CPU_RDTSC)) {
		uint32_t hi, lo;
		rdtsc (&hi, &lo);
		timers[i] -= (double) hi * 4294967296.0 + lo;
	} else {
		struct _timeb timeval;
		_ftime (&timeval);
		timers[i] -= (double) timeval.time * 1000.0 + timeval.millitm;
	}
}

void end_timer (
	double	*timers,
	int	i)
{
	if (RDTSC_TIMING < 10) {
		timers[i] += getHighResTimer ();
	} else if (RDTSC_TIMING > 10 && (CPU_FLAGS & CPU_RDTSC)) {
		uint32_t hi, lo;
		rdtsc (&hi, &lo);
		timers[i] += (double) hi * 4294967296.0 + lo;
	} else {
		struct _timeb timeval;
		_ftime (&timeval);
		timers[i] += (double) timeval.time * 1000.0 + timeval.millitm;
	}
}

void divide_timer (
	double	*timers,
	int	i,
	int	j)
{
	timers[i] = timers[i] / (double) j;
}

double timer_value (
	double	*timers,
	int	i)
{
	if (RDTSC_TIMING < 10)
		return (timers[i] / getHighResTimerFrequency ());
	else if (RDTSC_TIMING > 10 && (CPU_FLAGS & CPU_RDTSC))
		return (timers[i] / CPU_SPEED / 1000000.0);
	else
		return (timers[i] / 1000.0);
}

#define TIMER_NL	0x1
#define TIMER_CLR	0x2
#define TIMER_OPT_CLR	0x4
#define TIMER_MS	0x8

void print_timer (
	double	*timers,
	int	i,
	char	*buf,
	int	flags)
{
	double	t;

/* The timer could be less than zero if the computer went into hibernation. */
/* Hibernation is where the memory image is saved to disk and the computer */
/* shut off.  Upon power up the memory image is restored but the RDTSC */
/* timestamp counter has been reset to zero. */

	buf += strlen (buf);
	t = timer_value (timers, i);
	if (t < 0.0) {
		strcpy (buf, "Unknown");
		timers[i] = 0.0;
	}

/* Format the timer value in one of several styles */

	else {
		int	style;

		style = IniGetInt (INI_FILE, "TimingOutput", 0);
		if (style == 0) {
			if (flags & TIMER_MS) style = 4;
			else style = 1;
		}

		if (style == 1)
			sprintf (buf, "%.3f sec.", t);
		else if (style == 2)
			sprintf (buf, "%.1f ms.", t * 1000.0);
		else if (style == 3)
			sprintf (buf, "%.2f ms.", t * 1000.0);
		else
			sprintf (buf, "%.3f ms.", t * 1000.0);
		if (RDTSC_TIMING == 12 && (CPU_FLAGS & CPU_RDTSC)) {
			sprintf (buf+strlen(buf), " (%.0f clocks)", timers[i]);
		}
	}

/* Append optional newline */

	if (flags & TIMER_NL) strcat (buf, "\n");

/* Clear the timer */

	if (flags & TIMER_CLR) timers[i] = 0.0;
	if ((flags & TIMER_OPT_CLR) && !CUMULATIVE_TIMING) timers[i] = 0.0;
}

/****************************************************************************/
/*             Portable routines to launch worker threads                   */
/****************************************************************************/

/* Structure used in launching one worker thread. */

struct LaunchData {
	int	thread_num;		/* This thread number */
	unsigned int num_threads;	/* Num threads to run */
	unsigned long p;		/* Exponent to time */
	unsigned long iters;		/* Iterations to time */
	int	stop_reason;		/* Returned stop reason */
};

/* Create windows for the worker threads.  Windows REALLY prefers this be */
/* done in the main thread.  Otherwise, deadlocks can occur. */

void create_worker_windows (
	int	num_threads)
{
	int	tnum;
	char	buf[80];

/* Make sure each worker thread has a window to output to */

	for (tnum = 0; tnum < num_threads; tnum++) {
		create_window (tnum);
		if (NUM_CPUS * CPU_HYPERTHREADS > 1)
			sprintf (buf, "Worker #%d", tnum+1);
		else
			strcpy (buf, "Worker");
		base_title (tnum, buf);
	}
}

/* Launch the worker threads to process work units */

int LaunchWorkerThreads (
	int	thread_num,		/* Specific worker to launch or */
					/* special value ALL_WORKERS */
	int	wait_flag)		/* TRUE if we wait for workers to */
					/* end before returning. */
{
	struct LaunchData *ld;
	gwthread thread_handle;

/* If workers are already active, then call routine that restarts */
/* individual workers. */

	if (WORKER_THREADS_ACTIVE && (LAUNCH_TYPE == LD_CONTINUE || LAUNCH_TYPE == LD_TORTURE)) {
		if (thread_num == ALL_WORKERS) {
			for (thread_num = 0; thread_num < (int) NUM_WORKER_THREADS; thread_num++)
				if (! ACTIVE_WORKERS[thread_num])
					start_one_worker (thread_num);
		} else
			start_one_worker (thread_num);
		return (0);
	}

/* Create the launcher data structure, create the windows, then launch */

	ld = (struct LaunchData *) malloc (sizeof (struct LaunchData));
	if (ld == NULL) return (OutOfMemory (MAIN_THREAD_NUM));
	ld->num_threads = NUM_WORKER_THREADS;
	LAUNCH_TYPE = LD_CONTINUE;
	create_worker_windows (NUM_WORKER_THREADS);
	mark_workers_active (thread_num == ALL_WORKERS ? NUM_WORKER_THREADS : -thread_num);
	if (wait_flag) {
		gwthread_create_waitable (&thread_handle, &Launcher, ld);
		gwthread_wait_for_exit (&thread_handle);
	} else
		gwthread_create (&thread_handle, &Launcher, ld);
	return (0);
}

/* Launch threads to do a torture test */

int LaunchTortureTest (
	unsigned long num_threads,	/* Number of torture tests to run */
	int	wait_flag)		/* TRUE if we wait for workers to */
					/* end before returning. */
{
	struct LaunchData *ld;
	gwthread thread_handle;

	ld = (struct LaunchData *) malloc (sizeof (struct LaunchData));
	if (ld == NULL) return (OutOfMemory (MAIN_THREAD_NUM));
	ld->num_threads = num_threads;
	LAUNCH_TYPE = LD_TORTURE;
	create_worker_windows (num_threads);
	mark_workers_active (num_threads);
	if (wait_flag) {
		gwthread_create_waitable (&thread_handle, &Launcher, ld);
		gwthread_wait_for_exit (&thread_handle);
	} else
		gwthread_create (&thread_handle, &Launcher, ld);
	return (0);
}

/* Launch a thread to do a benchmark */

int LaunchBench (void)
{
	struct LaunchData *ld;
	gwthread thread_handle;

	ld = (struct LaunchData *) malloc (sizeof (struct LaunchData));
	if (ld == NULL) return (OutOfMemory (MAIN_THREAD_NUM));
	ld->num_threads = 1;
	LAUNCH_TYPE = LD_BENCH;
	create_worker_windows (1);
	mark_workers_active (1);
	gwthread_create (&thread_handle, &Launcher, ld);
	return (0);
}

/* Launch the worker thread(s) to process Advanced/Time */

int LaunchAdvancedTime (
	unsigned long p,		/* Exponent to time */
	unsigned long iters)		/* Iterations to time */
{
	struct LaunchData *ld;
	gwthread thread_handle;

	ld = (struct LaunchData *) malloc (sizeof (struct LaunchData));
	if (ld == NULL) return (OutOfMemory (MAIN_THREAD_NUM));
	if (p >= 9900 && p <= 9919) ld->num_threads = p - 9900 + 1;
	else if (p >= 9920 && p <= 9939) ld->num_threads = p - 9920 + 1;
	else ld->num_threads = 1;
	ld->p = p;
	ld->iters = iters;
	LAUNCH_TYPE = LD_TIME;
	create_worker_windows (ld->num_threads);
	mark_workers_active (ld->num_threads);
	gwthread_create (&thread_handle, &Launcher, ld);
	return (0);
}

/* Launch all worker threads */

void Launcher (void *arg)
{
	struct LaunchData *ld;
	unsigned int tnum;
	int	stop_reason;
	gwthread handles[MAX_NUM_WORKER_THREADS];
	struct LaunchData *ldwork[MAX_NUM_WORKER_THREADS];

/* This thread will create more worker threads if necessary and */
/* then become thread number 0. */

	ld = (struct LaunchData *) arg;

/* If worker threads are active then stop them all.  This can */
/* happen when we choose Torture Test, Benchmark, or Advanced/Time from */
/* the menus while the worker threads are running */

	if (WORKER_THREADS_ACTIVE) {
		stop_workers_for_escape ();
		while (WORKER_THREADS_STOPPING) Sleep (50);
	}

/* Set flags so that GUI knows worker threads are active */

	WORKER_THREADS_ACTIVE = ld->num_threads;
	WORKER_THREADS_STOPPING = FALSE;

/* Output a starting worker threads message */

	if (ld->num_threads > 1)
		OutputStr (MAIN_THREAD_NUM, "Starting workers.\n");
	else
		OutputStr (MAIN_THREAD_NUM, "Starting worker.\n");

/* Every time the user chooses Test/Continue, clear any timers that */
/* prevents communication for a period of time.  This allows the user */
/* to try something and if it doesn't work, ESC and choose Test/Continue */
/* to try some other system settings (without waiting an hour). */

	clear_comm_rate_limits ();

/* Clear array of active thread handles */

again:	clearThreadHandleArray ();

/* Reread prime.ini, local.ini, and worktodo.ini files just in case user */
/* hand edited it.  We don't officially support this, but we'll do it */
/* anyway.  Also, check for a .add file, which we do officially support. */
/* If the user edited the ini files changing the number of worker threads */
/* then handle that here.  We also jump here if the threads were restarted */
/* because the user changed the number of worker threads using dialog boxes. */
/* NOTE: If the user increases the number of threads, then he will not see */
/* worker windows until he does a stop and restart. */

	stop_reason = readIniFiles ();
	if (stop_reason) {
		OutputStr (MAIN_THREAD_NUM, "Error rereading INI files.\n");
		return;
	}
	if (LAUNCH_TYPE == LD_CONTINUE) ld->num_threads = NUM_WORKER_THREADS;

/* Initialize flags that cause the worker threads to stop at the */
/* appropriate time */

	init_stop_code ();

/* Init the code that keeps track of the memory used by each worker thread */

	init_mem_state ();

/* Run OS-specific code prior to launching the worker threads */

	PreLaunchCallback (LAUNCH_TYPE);

/* Change the icon */

	ChangeIcon (MAIN_THREAD_NUM, WORKING_ICON);

/* Start all appropriate timers */

	if (LAUNCH_TYPE == LD_CONTINUE) {

/* Start timer that tells us to write save files every so often */

		start_save_files_timer ();

/* Start the timer that checks battery status */

		start_battery_timer ();

/* Start the timer that checks for priority work */

		start_priority_work_timer ();

/* Start the timer that checks the pause-while-running list */

		start_pause_while_running_timer ();

/* Start the throttle timer */

		start_throttle_timer ();
	}

/* Launch more worker threads if needed */

	for (tnum = 1; tnum < ld->num_threads; tnum++) {
		ldwork[tnum] = (struct LaunchData *) malloc (sizeof (struct LaunchData));
		if (ldwork[tnum] == NULL) {
			OutOfMemory (MAIN_THREAD_NUM);
			return;
		}
		memcpy (ldwork[tnum], ld, sizeof (struct LaunchData));
		ldwork[tnum]->thread_num = tnum;
		gwthread_create_waitable (&handles[tnum], &LauncherDispatch, ldwork[tnum]);
	}

/* This thread is a worker thread too.  Call dispatching routine. */

	ld->thread_num = 0;
	LauncherDispatch (ld);
	stop_reason = ld->stop_reason;

/* Wait for other threads to finish */
/* Combine the stop reason with the stop reason returned by other threads */

	for (tnum = 1; tnum < ld->num_threads; tnum++) {
		gwthread_wait_for_exit (&handles[tnum]);
		if (stop_reason == 0)
			stop_reason = ldwork[tnum]->stop_reason;
		else if (stop_reason == STOP_ESCAPE ||
			 ldwork[tnum]->stop_reason == STOP_ESCAPE)
			stop_reason = STOP_ESCAPE;
		free (ldwork[tnum]);
	}

/* Write the worktodo file in case the WELL_BEHAVED_WORK flag caused us */
/* to delay writing the file. */

	if (LAUNCH_TYPE == LD_CONTINUE) {
		writeWorkToDoFile (TRUE);

/* Clear timers we started earlier */

		stop_save_files_timer ();
		stop_battery_timer ();
		stop_priority_work_timer ();
		stop_pause_while_running_timer ();
		stop_throttle_timer ();
	}

/* Change the icon */

	ChangeIcon (MAIN_THREAD_NUM, IDLE_ICON);

/* Run OS-specific code after worker threads terminate */

	PostLaunchCallback (LAUNCH_TYPE);

/* Restart all worker threads if the stop reason tells us to.  Make sure */
/* we set num_threads in case the reason for the restart is a change to */
/* NUM_WORKER_THREADS. */

	if (stop_reason == STOP_RESTART) {
		OutputStr (MAIN_THREAD_NUM, "Restarting all worker threads using new settings.\n");
		goto again;
	}

/* Restart all worker threads if the stop reason tells us to reread the */
/* INI file.  Make sure we set num_threads in case the reason for the restart */
/* is a change to NUM_WORKER_THREADS. */

	if (stop_reason == STOP_REREAD_INI) {
		OutputStr (MAIN_THREAD_NUM, "Restarting all worker threads using new timed prime.txt settings.\n");
		goto again;
	}

/* Output informative message */

	if (LAUNCH_TYPE == LD_CONTINUE || LAUNCH_TYPE == LD_TORTURE)
		OutputStr (MAIN_THREAD_NUM, "Execution halted.\n");
	if (LAUNCH_TYPE == LD_CONTINUE)
		OutputStr (MAIN_THREAD_NUM, "Choose Test/Continue to restart.\n");

/* Clear flags so that GUI knows worker threads are not active */

	WORKER_THREADS_ACTIVE = 0;
	WORKER_THREADS_STOPPING = FALSE;

/* Free the ld structure and exit the first worker thread */

	free (ld);
}

/* Now that the worker thread has been created, call the correct routine */
/* to do some work. */

void LauncherDispatch (void *arg)
{
	struct LaunchData *ld;
	int	stop_reason;

	ld = (struct LaunchData *) arg;

/* Change the title bar */

	title (ld->thread_num, "Starting");
	OutputStr (ld->thread_num, "Worker starting\n");
	ChangeIcon (ld->thread_num, WORKING_ICON);

/* Dispatch to the correct code */

	switch (LAUNCH_TYPE) {
	case LD_CONTINUE:
		stop_reason = primeContinue (ld->thread_num);
		break;
	case LD_TIME:
		stop_reason = primeTime (ld->thread_num, ld->p, ld->iters);
		break;
	case LD_BENCH:
		stop_reason = primeBench (ld->thread_num);
		break;
	case LD_TORTURE:
		stop_reason = tortureTest (ld->thread_num, ld->num_threads);
		break;
	}

/* Change the title bar and output a line to the window */

	title (ld->thread_num, "Not running");
	OutputStr (ld->thread_num, "Worker stopped.\n");
	ChangeIcon (ld->thread_num, IDLE_ICON);

/* Set the return code and exit this worker thread */

	ld->stop_reason = stop_reason;
}

/****************************************************************************/
/*                       Process the work units                             */
/****************************************************************************/

/* Continue factoring/testing Mersenne numbers */

int primeContinue (
	int	thread_num)
{
	struct PriorityInfo sp_info;
	struct work_unit *w;
	unsigned int pass, sequential;
	int	stop_reason;

/* Set the process/thread priority */

	sp_info.type = SET_PRIORITY_NORMAL_WORK;
	sp_info.thread_num = thread_num;
	sp_info.aux_thread_num = 0;
	SetPriority (&sp_info);

/* Loop until the ESC key is hit or the entire work-to-do INI file */
/* is processed and we are not connected to the server. */

	sequential = IniGetInt (INI_FILE, "SequentialWorkToDo", 1);
	for ( ; ; ) {

/* Check for a stop code.  We do this here in case the work-to-do file */
/* is empty (this call will be our only chance to check for a stop code). */

	stop_reason = stopCheck (thread_num);
	if (stop_reason) goto check_stop_code;

/* Clear flags that says we need to restart this thread if memory settings */
/* change.  If a work_unit cannot be processed because of a lack of */
/* available memory, then we will set these flags. */

	clear_memory_restart_flags (thread_num);

/* Make three passes over the worktodo.ini file looking for the ideal */
/* piece of work to do.  In pass 1, we look for high-priority work.  This */
/* includes trial and P-1 factoring prior to an LL test.  If a factor is */
/* found, it can reduce the amount of work we have queued up, requiring */
/* us to ask the server for more.  We also do AdvancedTest= lines in */
/* pass 1.  In pass 2, we process the file in order (except for LL tests */
/* that are not yet ready because the P-1 factoring has not completed). */
/* In pass 3, as a last resort we start P-1 stage 2 even if they will share */
/* memory with another P-1 in stage 2 and we start LL tests where P-1 */
/* factoring is stalled because of low memory. */
/* Skip first pass on large well-behaved work files. */

	for (pass = (WELL_BEHAVED_WORK || sequential) ? 2 : 1;
	     pass <= 3;
	     pass++) {

/* Examine each line in the worktodo.ini file */

	    for (w = NULL; ; ) {

/* Read the line from the work file, break when out of lines */
/* Skip comment lines from worktodo.ini */
		
		w = getNextWorkToDoLine (thread_num, w, LONG_TERM_USE);
		if (w == NULL) break;
		if (w->work_type == WORK_NONE) continue;

/* Clear flags indicating this work_unit is using a lot of memory */

		set_default_memory_usage (thread_num);

/* Handle a factoring assignment */

		if (w->work_type == WORK_FACTOR && pass == 2) {
			stop_reason = primeFactor (thread_num, &sp_info, w, 0);
		}

/* Do special P-1 factoring work. */

		if (w->work_type == WORK_PFACTOR && pass == 2) {
			stop_reason = pfactor (thread_num, &sp_info, w, pass);
		}

/* Run the LL test */

		if (w->work_type == WORK_ADVANCEDTEST ||
		    w->work_type == WORK_TEST ||
		    w->work_type == WORK_DBLCHK) {
			stop_reason = prime (thread_num, &sp_info, w, pass);
		}

/* See if this is an ECM factoring line */

		if (w->work_type == WORK_ECM && pass == 2) {
			stop_reason = ecm (thread_num, &sp_info, w, pass);
		}

/* See if this is an P-1 factoring line */

		if (w->work_type == WORK_PMINUS1 && pass == 2) {
			stop_reason = pminus1 (thread_num, &sp_info, w, pass);
		}

/* Run a PRP test */

		if (w->work_type == WORK_PRP) {
			stop_reason = prp (thread_num, &sp_info, w, pass);
		}

/* Set us back to default memory usage */

		set_default_memory_usage (thread_num);

/* If the work unit completed remove it from the worktodo.ini file and */
/* move on to the next entry */

		if (stop_reason == STOP_WORK_UNIT_COMPLETE) {
			rolling_average_work_unit_complete (thread_num, w);
			stop_reason = deleteWorkToDoLine (thread_num, w, FALSE);
		}

/* If a work unit could not be processed because there isn't enough memory, */
/* then move on to the next worktodo entry while we wait for more memory. */

		if (stop_reason == STOP_NOT_ENOUGH_MEM) {
			OutputStr (thread_num, "Looking for work that uses less memory.\n");
			stop_reason = 0;
		}

/* If we are aborting this work unit (probably because it is being deleted) */
/* then print a message. */

		if (stop_reason == STOP_ABORT)
			OutputStr (thread_num, "Aborting processing of this work unit.\n");

/* If stop reason is set then unlock this work unit and go process the */
/* stop reason.  Otherwise, no work was done, move on to the next entry */
/* in the worktodo.ini file. */

		if (stop_reason) {
			decrementWorkUnitUseCount (w, LONG_TERM_USE);
			goto check_stop_code;
		}

/* Process next work unit in the current pass */

	    }

/* Make another pass over the worktodo.ini file */

	}

/* Check for all the possible stop codes we must handle here.  Those */
/* that terminate the worker thread are not handled here. */

check_stop_code:

/* If we are aborted a work unit (probably because it is being deleted) */
/* then start again. */

	if (stop_reason == STOP_ABORT) continue;

/* If we need to do priority work then reprocess the entire worktodo.ini. */

	if (stop_reason == STOP_PRIORITY_WORK) continue;

/* If we need to restart with the new memory settings, do so. */

	if (stop_reason == STOP_MEM_CHANGED) continue;

/* If the user is specifically stopping this worker, then stop until */
/* the user restarts the worker. */

	if (stop_reason == STOP_WORKER) {
		implement_stop_one_worker (thread_num);
		continue;
	}

/* If the worker is pausing because another program is running */
/* then implement that now. */

	if (stop_reason == STOP_PAUSE) {
		implement_pause (thread_num);
		continue;
	}

/* If the worker is pausing because we are now on battery power, then */
/* implement that now. */

	if (stop_reason == STOP_BATTERY) {
		implement_stop_battery (thread_num);
		continue;
	}

/* The stop reason was not caught above.  It must be a fatal error or a */
/* stop code that causes the worker thread to terminate. */

	if (stop_reason) return (stop_reason);

/* Ugh, we made three passes over the worktodo file and couldn't find */
/* any work to do.  I think this can only happen if we are low on memory */
/* or the worktodo file is empty. */

//bug? - only do this if two attempts are made at executing work?  Because
// work might have been added to the front of the file???

/* Output a message saying this worker thread is waiting for work */

	title (thread_num, "Waiting for work");
	OutputStr (thread_num, "No work to do at the present time.  Waiting.\n");
	ChangeIcon (thread_num, IDLE_ICON);

/* Set memory usage to zero */

	set_memory_usage (thread_num, 0, 0);

/* Spool a message to check the work queue.  Since we have no work queued */
/* up, this should cause us to get some work from the server. */

	spoolMessage (MSG_CHECK_WORK_QUEUE, NULL);

/* Wait for a mem-changed event OR communication attempt (it might get work) */
/* OR user entering new work via the dialog boxes OR the discovery of a .add */
/* file OR wait for a thread stop event (like ESC or shutdown). */

	WORK_AVAILABLE_OR_STOP_INITIALIZED[thread_num] = 1;
	gwevent_init (&WORK_AVAILABLE_OR_STOP[thread_num]);
	gwevent_reset (&WORK_AVAILABLE_OR_STOP[thread_num]);
	gwevent_wait (&WORK_AVAILABLE_OR_STOP[thread_num], 3600);
	gwevent_destroy (&WORK_AVAILABLE_OR_STOP[thread_num]);
	WORK_AVAILABLE_OR_STOP_INITIALIZED[thread_num] = 0;
	OutputStr (thread_num, "Resuming.\n");
	ChangeIcon (thread_num, WORKING_ICON);

/* Loop scanning the work-to-do file.  Hopefully the event triggered */
/* because we now have work to do. */

	}
}

/************************/
/* Trial Factoring code */
/************************/

/* This defines the C / assembly language communication structure */

#define NEW_STACK_SIZE	(4096+256)
struct facasm_data {
	uint32_t EXPONENT;		/* Mersenne number to factor */
	uint32_t FACPASS;		/* Which of 16 factoring passes */
	uint32_t FACHSW;		/* High word of found factor */
	uint32_t FACMSW;		/* Middle word of found factor */
	uint32_t FACLSW;		/* Low word of found factor */
	uint32_t cpu_flags;		/* Copy of CPU_FLAGS */
	uint32_t firstcall;		/* Flag set on first facpasssetup */
	uint32_t pad[5];
	uint32_t xmm_data[100];		/* XMM data initialized in C code */
};

/* This defines the factoring data handled in C code.  The handle */
/* abstracts all the internal details from callers of the factoring code. */

typedef struct {
	struct	facasm_data *asm_data;	/* Memory for factoring code */
} fachandle;

EXTERNC void setupf (struct facasm_data *);	/* Assembly code, setup */
EXTERNC int factor64 (struct facasm_data *);	/* Assembly code, do work */

/* Prepare a factoring run */

int factorSetup (
	int	thread_num,
	unsigned long p,
	fachandle *facdata)
{
	void	*asm_data_alloc;
	struct facasm_data *asm_data;

/* Allocate 1MB for the assembly code global data.  This area is preceded */
/* by a temporary stack.  This allows the assembly code to access the global */
/* data using offsets from the stack pointer.  We zero the first 64KB, */
/* asm code requires this (such as XMM_COMPARE_VALn). */

	asm_data_alloc = aligned_malloc (1000000, 4096);
	if (asm_data_alloc == NULL) {
		OutputStr (thread_num, "Error allocating memory for trial factoring.\n");
		return (STOP_OUT_OF_MEM);
	}
	facdata->asm_data = asm_data = (struct facasm_data *)
		((char *) asm_data_alloc + NEW_STACK_SIZE);
	memset (asm_data, 0, 65536);

/* Init */

	asm_data->EXPONENT = p;
	asm_data->cpu_flags = CPU_FLAGS;
	asm_data->firstcall = 0;

/* Setup complete */

	return (0);
}

/* Prepare for one of the 16 factoring passes */

int factorPassSetup (
	int	thread_num,
	unsigned long pass,
	fachandle *facdata)		/* Handle returned by factorSetup */
{
	struct facasm_data *asm_data;

/* Call the factoring setup assembly code */

	asm_data = (struct facasm_data *) facdata->asm_data;
	asm_data->FACPASS = pass;
	setupf (asm_data);

/* If using the SSE2 factoring code, do more initialization */
/* We need to initialize much of the following data: */
/*	XMM_INITVAL		DD	0,0,0,0
	XMM_INVFAC		DD	0,0,0,0
	XMM_I1			DD	0,0,0,0
	XMM_I2			DD	0,0,0,0
	XMM_F1			DD	0,0,0,0
	XMM_F2			DD	0,0,0,0
	XMM_F3			DD	0,0,0,0
	XMM_TWO_120_MODF1	DD	0,0,0,0
	XMM_TWO_120_MODF2	DD	0,0,0,0
	XMM_TWO_120_MODF3	DD	0,0,0,0
	XMM_INIT120BS		DD	0,0
	XMM_INITBS		DD	0,0
	XMM_BS			DD	0,0
	XMM_SHIFTER		DD	64 DUP (0)
	TWO_TO_FACSIZE_PLUS_62	DQ	0.0
	SSE2_LOOP_COUNTER	DD	0 */

#ifndef X86_64
	if (asm_data->cpu_flags & CPU_SSE2) {
		unsigned long i, p, bits_in_factor;
		uint32_t *xmm_data;

/* Compute the number of bits in the factors we will be testing */

		if (asm_data->FACHSW)
			bits_in_factor = 64, i = asm_data->FACHSW;
		else if (asm_data->FACMSW)
			bits_in_factor = 32, i = asm_data->FACMSW;
		else return (0);
		while (i) bits_in_factor++, i >>= 1;

/* Factors 63 bits and below use the non-SSE2 code */

		if (bits_in_factor <= 63) return (0);

/* Set XMM_SHIFTER values (the first shifter value is not used). */
/* Also compute the initial value. */

		xmm_data = asm_data->xmm_data;
		p = asm_data->EXPONENT;
		for (i = 0; p > bits_in_factor + 59; i++) {
			xmm_data[48+i*2] = (p & 1) ? 1 : 0;
			p >>= 1;
		}
		xmm_data[0] =			/* XMM_INITVAL */
		xmm_data[2] = p >= 90 ? 0 : (1 << (p - 60));
		xmm_data[40] = 62 - (120 - bits_in_factor);/* XMM_INIT120BS */
		xmm_data[42] = 62 - (p - bits_in_factor);/* XMM_INITBS */
		xmm_data[112] = i;		/* SSE2_LOOP_COUNTER */
		*(double *)(&xmm_data[110]) =	/* TWO_TO_FACSIZE_PLUS_62 */
			pow ((double) 2.0, (int) (bits_in_factor + 62));

/* Set XMM_BS to 60 - (120 - fac_size + 1) as defined in factor64.mac */

		xmm_data[44] = bits_in_factor - 61;
	}
#endif

/* Setup complete */

	return (0);
}

/* Factor one "chunk".  The assembly code decides how big a chunk is. */

int factorChunk (
	fachandle *facdata)		/* Handle returned by factorSetup */
{
	return (factor64 (facdata->asm_data));
}

/* Cleanup after making a factoring run */

void factorDone (
	fachandle *facdata)		/* Handle returned by factorSetup */
{

/* Free assembly code work area */

	if (facdata->asm_data != NULL) {
		aligned_free ((char *) facdata->asm_data - NEW_STACK_SIZE);
		facdata->asm_data = NULL;
	}
}

/* Wrapper code that verifies any factors found by the assembly code */
/* res is set to 2 if a factor was not found, 1 otherwise */

int factorAndVerify (
	int	thread_num,
	unsigned long p,
	fachandle *facdata,
	int	*res)
{
	unsigned long hsw, msw;
	int	stop_reason;

/* Remember starting point in case of an error */

	hsw = facdata->asm_data->FACHSW;
	msw = facdata->asm_data->FACMSW;

/* Call assembly code */


loop:	*res = factorChunk (facdata);

/* If a factor was not found, return. */

	if (*res == 2) return (stopCheck (thread_num));

/* Otherwise verify the factor. */

	if (facdata->asm_data->FACHSW ||
	    facdata->asm_data->FACMSW ||
	    facdata->asm_data->FACLSW > 1) {
		giant	f, x;

		f = newgiant (100);
		itog ((int) facdata->asm_data->FACHSW, f);
		gshiftleft (32, f);
		uladdg (facdata->asm_data->FACMSW, f);
		gshiftleft (32, f);
		uladdg (facdata->asm_data->FACLSW, f);

		x = newgiant (100);
		itog (2, x);
		powermod (x, p, f);
		*res = isone (x);

		free (f);
		free (x);

		if (*res) return (0);
	}

/* If factor is no good, print an error message, sleep, and */
/* restart the factoring code. */

	OutputBoth (thread_num, "ERROR: Incorrect factor found.\n");
	facdata->asm_data->FACHSW = hsw;
	facdata->asm_data->FACMSW = msw;
	stop_reason = SleepFive (thread_num);
	if (stop_reason) return (stop_reason);
	stop_reason = factorSetup (thread_num, p, facdata);
	if (stop_reason) return (stop_reason);
	stop_reason = factorPassSetup (thread_num, p, facdata);
	if (stop_reason) return (stop_reason);
	goto loop;
}

/* Trial factor a Mersenne number prior to running a Lucas-Lehmer test */

char FACMSG[] = "Trial factoring M%%ld to 2^%%d is %%.%df%%%% complete.";
char SHORT_FACMSG[] = "Trial factoring M%ld to 2^%d.";

#define FACTOR_MAGICNUM		0x1567234D
#define FACTOR_VERSION		1

int primeFactor (
	int	thread_num,
	struct PriorityInfo *sp_info,	/* SetPriority information */
	struct work_unit *w,
	unsigned int factor_limit_adjustment)
{
	fachandle facdata;		/* Handle to the factoring data */
	unsigned long p;		/* Exponent to factor */
	unsigned long bits;		/* How far already factored in bits */
	unsigned long test_bits;	/* How far to factor to */
	long	factor_found;		/* Returns true if factor found */
	int	fd;			/* Continuation file handle or zero */
	int	first_iter_msg, continuation, stop_reason, find_smaller_factor;
	unsigned long endpthi, endptlo;
	double	endpt, startpt;		/* For computing percent complete */
	unsigned long pass;		/* Factoring pass 0 through 15 */
	unsigned long report_bits;	/* When to report results one bit */
					/* at a time */
	char	filename[32];
	char	buf[200], str[80];
	double	timers[2];

/* Init */

	factor_found = 0;
	p = w->n;
	bits = (unsigned int) w->sieve_depth;

/* Determine how much we should factor (in bits) */

	test_bits = (unsigned int) w->factor_to - factor_limit_adjustment;

/* Is exponent already factored enough? This should never happen with */
/* WORK_FACTOR work units.  However, I suppose the user could have */
/* manually changed the line in worktodo.ini.  So send a message to */
/* server saying we didn't do any factoring but we are done with */
/* this work unit.  Then delete the work unit. */

	if (bits >= test_bits) {
		if (w->work_type == WORK_FACTOR) {
			struct primenetAssignmentResult pkt;
			memset (&pkt, 0, sizeof (pkt));
			strcpy (pkt.computer_guid, COMPUTER_GUID);
			strcpy (pkt.assignment_uid, w->assignment_uid);
			pkt.result_type = PRIMENET_AR_TF_NOFACTOR;
			pkt.n = p;
			pkt.done = TRUE;
			spoolMessage (PRIMENET_ASSIGNMENT_RESULT, &pkt);
			return (STOP_WORK_UNIT_COMPLETE);
		}
		return (0);
	}

/* Setup the factoring code */

	stop_reason = factorSetup (thread_num, p, &facdata);
	if (stop_reason) return (stop_reason);

/* Record the amount of memory being used by this thread (1MB). */

	set_memory_usage (thread_num, 0, 1 << 20);

/* Check for a v24 continuation file.  These were named pXXXXXXX.  The */
/* first 16 bits contained a 2 to distinguish it from a LL save file. */
/* In v25, we name the file fXXXXXXX and use the common header format */
/* to make Test/Status and computing completion dates easier. */

	continuation = FALSE;
	tempFileName (w, filename);
	filename[0] = 'p';
	fd = _open (filename, _O_BINARY | _O_RDONLY);
	if (fd > 0) {
		short	type;
		short	shortdummy;
		unsigned long longdummy, fachsw, facmsw;
		short	file_factor_found, file_bits, file_pass;

		if (read_short (fd, &type) &&
		    type == 2 &&
		    read_long (fd, &longdummy, NULL) &&
		    read_short (fd, &shortdummy) &&
		    read_short (fd, &file_factor_found) &&
		    read_short (fd, &shortdummy) &&
		    read_short (fd, &file_bits) &&
		    read_short (fd, &file_pass) &&
		    read_long (fd, &fachsw, NULL) &&
		    read_long (fd, &facmsw, NULL) &&
		    read_long (fd, &endpthi, NULL) &&
		    read_long (fd, &endptlo, NULL)) {
			OutputBoth (thread_num, "Using old-style factoring save file.\n");
			facdata.asm_data->FACHSW = fachsw;
			facdata.asm_data->FACMSW = facmsw;
			factor_found = file_factor_found;
			bits = file_bits;
			pass = file_pass;
			continuation = TRUE;
			_close (fd);
			_unlink (filename);
		} else {
			_close (fd);
		}
	}
	
/* Read v25+ continuation file */

	filename[0] = 'f';
	if (!continuation &&
	    (fd = _open (filename, _O_BINARY | _O_RDONLY)) > 0) {
		unsigned long version, sum, fachsw, facmsw;

		if (read_magicnum (fd, FACTOR_MAGICNUM) &&
		    read_header (fd, &version, w, &sum) &&
		    version == FACTOR_VERSION &&
		    read_long (fd, (unsigned long *) &factor_found, NULL) &&
		    read_long (fd, &bits, NULL) &&
		    read_long (fd, &pass, NULL) &&
		    read_long (fd, &fachsw, NULL) &&
		    read_long (fd, &facmsw, NULL) &&
		    read_long (fd, &endpthi, NULL) &&
		    read_long (fd, &endptlo, NULL)) {
			_close (fd);
			facdata.asm_data->FACHSW = fachsw;
			facdata.asm_data->FACMSW = facmsw;
			continuation = TRUE;
		} else {
			sprintf (buf, READFILEERR, filename);
			OutputBoth (thread_num, buf);
			_close (fd);
			_unlink (filename);
		}
	}

/* Init the title */

	sprintf (buf, "Factoring M%ld", p);
	title (thread_num, buf);
	sprintf (buf, "%s trial factoring of M%ld to 2^%d\n",
		 fd > 0 ? "Resuming" : "Starting", p, test_bits);
	OutputStr (thread_num, buf);

/* Clear all timers */

	clear_timers (timers, sizeof (timers) / sizeof (timers[0]));

/* When doing easy and quick trial factoring on a Mersenne number, */
/* do not send a message to the server for every bit level we complete. */
/* If we did, the client would spend more CPU time sending messages to the */
/* server than actually factoring numbers.  Here we calculate the threshold */
/* where we'll start reporting results one bit at time.  We've arbitrarily */
/* chosen the difficulty in trial factoring M80000000 to 2^60 as the */
/* point where it is worthwhile to report results one bit at a time. */

	report_bits = (unsigned long)
		(60.0 + log ((double) p / 80000000.0) / log (2.0));
	if (report_bits >= test_bits) report_bits = test_bits;

/* Loop testing larger and larger factors until we've tested to the */
/* appropriate number of bits.  Advance one bit at a time because it */
/* is faster to look for factors at lower bit levels first. */
/* We always enter this loop if there is a continuation file because v23 */
/* had higher factoring limits and if we upgrade to v25 midstream, we */
/* might not send a factoring complete message to the server if we don't */
/* finish off the current bit level. */

	while (test_bits > bits || continuation) {
	    unsigned int end_bits;
	    unsigned long iters, iters_r;

/* Advance one bit at a time to minimize wasted time looking for a */
/* second factor after a first factor is found. */

	    end_bits = (bits < 50) ? 50 : bits + 1;
	    if (end_bits > test_bits) end_bits = test_bits;
	    sprintf (w->stage, "TF%d", end_bits);

/* Compute the ending point for each pass */

	    if (!continuation) {
		if (end_bits < 64) {
			endpthi = 0;
			endptlo = 1L << (end_bits-32);
		} else {
			endpthi = 1L << (end_bits-64);
			endptlo = 0;
		}
	    }

/* Precompute some constant for calculating percent complete */

	    if (bits < 32) startpt = 0.0;
	    else startpt = pow ((double) 2.0, (int) (bits-32));
	    endpt = endpthi * 4294967296.0 + endptlo;

/* Sixteen passes.  Two for the 1 or 7 mod 8 factors times two for the */
/* 1 or 2 mod 3 factors times four for the 1, 2, 3, or 4 mod 5 factors. */

	    iters_r = 0;
	    iters = 0;
	    first_iter_msg = (continuation ? 1 : 2);
	    if (! continuation) pass = 0;
	    for ( ; pass < 16; pass++) {

/* Set the starting point only if we are not resuming from */
/* a continuation file.  For no particularly good reason we */
/* quickly redo trial factoring for factors below 2^50. */

		if (continuation)
			continuation = FALSE;
		else {
			if (bits < 50) {
				facdata.asm_data->FACHSW = 0;
				facdata.asm_data->FACMSW = 0;
			} else if (bits < 64) {
				facdata.asm_data->FACHSW = 0;
				facdata.asm_data->FACMSW = 1L << (bits-32);
			} else {
				facdata.asm_data->FACHSW = 1L << (bits-64);
				facdata.asm_data->FACMSW = 0;
			}
		}

/* Only test for factors less than 2^32 on the first pass */

		if (facdata.asm_data->FACHSW == 0 &&
		    facdata.asm_data->FACMSW == 0 && pass != 0)
			facdata.asm_data->FACMSW = 1;

/* Setup the factoring program */

		stop_reason = factorPassSetup (thread_num, pass, &facdata);
		if (stop_reason) {
			factorDone (&facdata);
			return (stop_reason);
		}

/* Loop until all factors tested or factor found */

		for ( ; ; ) {
			int	res;
			double	currentpt;

/* Do a chunk of factoring */

			start_timer (timers, 0);
#ifdef SERVER_TESTING
			if (facdata.asm_data->FACMSW >= 0xFFF00000) {
				facdata.asm_data->FACHSW++;
				facdata.asm_data->FACMSW = 0;
			} else
				facdata.asm_data->FACMSW += 0x100000;
			stop_reason = stopCheck (thread_num);
			if (rand () == 1234 && rand () < 3000) res = 0;
			else res = 2;
#else
			stop_reason = factorAndVerify (thread_num, p, &facdata, &res);
#endif
			end_timer (timers, 0);
			if (res != 2) break;

/* Compute new percentage complete (of this bit level) */

			currentpt = facdata.asm_data->FACHSW * 4294967296.0 +
				    facdata.asm_data->FACMSW;
			if (currentpt > endpt) currentpt = endpt;
			w->pct_complete =
				(pass + (currentpt - startpt) /
					(endpt - startpt)) / 16.0;

/* Output informative message.  Usually this includes a percent complete, however, */
/* when just beginning a bit level (first_iter_msg == 2) we don't as the percentage */
/* is close to zero. */

			if (++iters >= ITER_OUTPUT || first_iter_msg) {
				char	fmt_mask[80];
				double	pct;
				pct = trunc_percent (w->pct_complete);
				if (first_iter_msg == 2) {
					sprintf (buf, "M%ld to 2^%d", p, end_bits);
				} else {
					sprintf (fmt_mask,
						 "%%.%df%%%% of M%%ld to 2^%%d",
						 PRECISION);
					sprintf (buf, fmt_mask, pct, p, end_bits);
				}
				title (thread_num, buf);
				if (first_iter_msg == 2) {
					sprintf (buf, SHORT_FACMSG, p, end_bits);
				} else {
					sprintf (fmt_mask, FACMSG, PRECISION);
					sprintf (buf, fmt_mask, p, end_bits, pct);
				}
				if (first_iter_msg) {
					strcat (buf, "\n");
					clear_timer (timers, 0);
				} else {
					strcat (buf, "  Time: ");
					print_timer (timers, 0, buf,
						     TIMER_NL | TIMER_OPT_CLR);
				}
				OutputStr (thread_num, buf);
				iters = 0;
				first_iter_msg = FALSE;
			}

/* Output informative message */

			if (++iters_r >= ITER_OUTPUT_RES ||
			    (NO_GUI && stop_reason)) {
				char	fmt_mask[80];
				double	pct;
				pct = trunc_percent (w->pct_complete);
				sprintf (fmt_mask, FACMSG, PRECISION);
				sprintf (buf, fmt_mask, p, end_bits, pct);
				strcat (buf, "\n");
				writeResults (buf);
				iters_r = 0;
			}

/* If an escape key was hit, write out the results and return */

			if (stop_reason || testSaveFilesFlag (thread_num)) {
				fd = _open (filename, _O_BINARY | _O_WRONLY | _O_CREAT, CREATE_FILE_ACCESS);
				write_header (fd, FACTOR_MAGICNUM, FACTOR_VERSION, w);
				write_long (fd, factor_found, NULL);
				write_long (fd, bits, NULL);
				write_long (fd, pass, NULL);
				write_long (fd, facdata.asm_data->FACHSW, NULL);
				write_long (fd, facdata.asm_data->FACMSW, NULL);
				write_long (fd, endpthi, NULL);
				write_long (fd, endptlo, NULL);
				_commit (fd);
				_close (fd);
				if (stop_reason) {
					factorDone (&facdata);
					return (stop_reason);
				}
			}

/* Test for completion */

			if (facdata.asm_data->FACHSW > endpthi ||
			    (facdata.asm_data->FACHSW == endpthi &&
			     facdata.asm_data->FACMSW >= endptlo))
				goto nextpass;
		}

/* Set flag indicating a factor has been found! */

		factor_found = TRUE;

/* We used to continue factoring to find a smaller factor in a later pass. */
/* We'll continue to do this if the found factor is really small (less than */
/* 2^56) or if the user sets FindSmallestFactor in prime.ini. */

		find_smaller_factor =
			(end_bits <= (unsigned int) IniGetInt (INI_FILE, "FindSmallestFactor", 56));

/* Format and output a message */

		makestr (facdata.asm_data->FACHSW,
			 facdata.asm_data->FACMSW,
			 facdata.asm_data->FACLSW, str);
		sprintf (buf, "M%ld has a factor: %s\n", p, str);
		OutputStr (thread_num, buf);
		formatMsgForResultsFile (buf, w);
		writeResults (buf);

/* Send assignment result to the server.  To avoid flooding the server */
/* with small factors from users needlessly redoing factoring work, make */
/* sure the factor is more than 50 bits or so. */

		if (strlen (str) >= 15 ||
		    IniGetInt (INI_FILE, "SendAllFactorData", 0)) {
			struct primenetAssignmentResult pkt;
			memset (&pkt, 0, sizeof (pkt));
			strcpy (pkt.computer_guid, COMPUTER_GUID);
			strcpy (pkt.assignment_uid, w->assignment_uid);
			strcpy (pkt.message, buf);
			pkt.result_type = PRIMENET_AR_TF_FACTOR;
			pkt.n = p;
			strcpy (pkt.factor, str);
			pkt.start_bits = (bits < report_bits) ?
				     (unsigned int) w->sieve_depth : bits;
			pkt.done = !find_smaller_factor;
			spoolMessage (PRIMENET_ASSIGNMENT_RESULT, &pkt);
		}

/* If we're looking for smaller factors, set a new end point.  Otherwise, */
/* skip all remaining passes. */

		if (!find_smaller_factor) break;

		if (facdata.asm_data->FACMSW != 0xFFFFFFFF) {
			endpthi = facdata.asm_data->FACHSW;
			endptlo = facdata.asm_data->FACMSW+1;
		} else {
			endpthi = facdata.asm_data->FACHSW+1;
			endptlo = 0;
		}
	        endpt = endpthi * 4294967296.0 + endptlo;

/* Do next of the 16 passes */

nextpass:	;
	    }

/* If we've found a factor then we need to send an assignment done */
/* message if we continued to look for a smaller factor. */

	    if (factor_found) {
		if (w->assignment_uid[0] && find_smaller_factor) {
			struct primenetAssignmentResult pkt;
			memset (&pkt, 0, sizeof (pkt));
			strcpy (pkt.computer_guid, COMPUTER_GUID);
			strcpy (pkt.assignment_uid, w->assignment_uid);
			pkt.result_type = PRIMENET_AR_NO_RESULT;
			pkt.done = TRUE;
			spoolMessage (PRIMENET_ASSIGNMENT_RESULT, &pkt);
		}
		break;
	    }

/* Output a no factor found message */

	    if (end_bits >= report_bits) {
		unsigned int start_bits;

		start_bits = (end_bits == report_bits) ?
				(unsigned int) w->sieve_depth : bits;
		if (start_bits < 32)
		    sprintf (buf,
			     "M%ld no factor to 2^%d, Wd%d: %08lX\n",
			     p, end_bits, PORT, SEC3 (p));
		else
		    sprintf (buf,
			     "M%ld no factor from 2^%d to 2^%d, Wd%d: %08lX\n",
			     p, start_bits, end_bits, PORT, SEC3 (p));
		OutputStr (thread_num, buf);
		formatMsgForResultsFile (buf, w);
		writeResults (buf);

/* Send no factor found message to the server for each bit */
/* level (i.e. one bit at a time).  As always to avoid swamping */
/* the server with needless data, do not send small bit level */
/* messages - that work has already been done. */

		if (end_bits >= 50 ||
		    IniGetInt (INI_FILE, "SendAllFactorData", 0)) {
			struct primenetAssignmentResult pkt;
			memset (&pkt, 0, sizeof (pkt));
			strcpy (pkt.computer_guid, COMPUTER_GUID);
			strcpy (pkt.assignment_uid, w->assignment_uid);
			strcpy (pkt.message, buf);
			pkt.result_type = PRIMENET_AR_TF_NOFACTOR;
			pkt.n = p;
			pkt.start_bits = start_bits;
			pkt.end_bits = end_bits;
			pkt.done = (w->work_type == WORK_FACTOR) &&
				   (end_bits >= test_bits);
			spoolMessage (PRIMENET_ASSIGNMENT_RESULT, &pkt);
		}
	    }

/* Advance the how far factored variable */

	    bits = end_bits;
	}

/* Clean up allocated factoring data */

	factorDone (&facdata);

/* Delete the continuation file */

	_unlink (filename);

/* If we found a factor, then we likely performed much less work than */
/* we estimated.  Make sure we do not update the rolling average with */
/* this inaccurate data. */

	if (factor_found) invalidateNextRollingAverageUpdate ();

/* If we finished this work unit, return the happy news */

	if (factor_found || w->work_type == WORK_FACTOR)
		return (STOP_WORK_UNIT_COMPLETE);

/* Update the worktodo file */

	w->sieve_depth = bits;
	stop_reason = updateWorkToDoLine (thread_num, w);
	if (stop_reason) return (stop_reason);

/* All done */

	return (0);
}

/***************************************/
/* Routines to run a Lucas-Lehmer test */
/***************************************/

/* Structure for holding lucas setup data */

typedef struct {		/* Some of the data kept during LL test */
	gwhandle gwdata;	/* When we multithread the gwnum code, */
				/* gwsetup will return a handle */
	gwnum	lldata;		/* Number in the lucas sequence */
	unsigned long units_bit; /* Shift count */
} llhandle;

/* Prepare for running a Lucas-Lehmer test.  Caller must have already */
/* called gwinit. */

int lucasSetup (
	int	thread_num,	/* Worker thread number */
	unsigned long p,	/* Exponent to test */
	unsigned long fftlen,	/* Specific FFT length to use, or zero */
	llhandle *lldata)	/* Common LL data structure */
{
	int	res;

/* Init LL data structure */

	lldata->lldata = NULL;
	lldata->units_bit = 0;

/* Init the FFT code for squaring modulo 1.0*2^p-1.  NOTE: As a kludge for */
/* the benchmarking code, an odd FFTlen sets up the 1.0*2^p+1 FFT code. */

	gwset_specific_fftlen (&lldata->gwdata, fftlen & ~1);
	if (fftlen & 1)
		res = gwsetup (&lldata->gwdata, 1.0, 2, p, 1);
	else
		res = gwsetup (&lldata->gwdata, 1.0, 2, p, -1);

/* If we were unable to init the FFT code, then print an error message */
/* and return an error code.  There is one exception, when we are doing */
/* a benchmark of all possible FFT implementations, do not print an error */
/* message. */

	if (res) {
		if (!lldata->gwdata.bench_pick_nth_fft) {
			char	buf[80];
			sprintf (buf, "Cannot initialize FFT code, errcode=%d\n", res);
			OutputBoth (thread_num, buf);
		}
		return (STOP_FATAL_ERROR);
	}
		
/* Allocate memory for the Lucas-Lehmer data (the number to square) */

	lldata->lldata = gwalloc (&lldata->gwdata);
	if (lldata->lldata == NULL) {
		gwdone (&lldata->gwdata);
		OutputStr (thread_num, "Error allocating memory for FFT data.\n");
		return (STOP_OUT_OF_MEM);
	}
	return (0);
}

/* Clean up after running a Lucas-Lehmer test */

void lucasDone (
	llhandle *lldata)	/* Common LL data structure */
{

/* Free memory for the Lucas-Lehmer data */

	gwfree (&lldata->gwdata, lldata->lldata);

/* Cleanup the FFT code */

	gwdone (&lldata->gwdata);
}

/* Generate the 64-bit residue of a Lucas-Lehmer test.  Returns -1 for an */
/* illegal result, 0 for a zero result, 1 for a non-zero result. */

int generateResidue64 (
	llhandle *lldata,
	unsigned long *reshi,
	unsigned long *reslo)
{
	giant	tmp;
	int	err_code;

	*reshi = *reslo = 0;
	tmp = popg (&lldata->gwdata.gdata, (lldata->gwdata.n >> 5) + 5);
	err_code = gwtogiant (&lldata->gwdata, lldata->lldata, tmp);
	if (err_code < 0) return (err_code);
	if (tmp->sign == 0) return (0);
	gshiftright (lldata->units_bit, tmp);
	if (tmp->sign > 0) *reslo = tmp->n[0];
	if (tmp->sign > 1) *reshi = tmp->n[1];
	pushg (&lldata->gwdata.gdata, 1);
	return (1);
}

/* Return TRUE if a continuation file exists.  If one does exist, */
/* make sure it is named pXXXXXXX. */

int continuationFileExists (
	int	thread_num,
	char	*filename)
{
	char	backupname[32];
	char	buf[80];

	if (fileExists (filename)) return (TRUE);
	strcpy (backupname, filename);
	backupname[0] = 'q';
	if (fileExists (backupname)) {
		sprintf (buf, RENAME_MSG, backupname, filename);
		OutputBoth (thread_num, buf);
		rename (backupname, filename);
		return (TRUE);
	}
	backupname[0] = 'r';
	if (fileExists (backupname)) {
		sprintf (buf, RENAME_MSG, backupname, filename);
		OutputBoth (thread_num, buf);
		rename (backupname, filename);
		return (TRUE);
	}
	return (FALSE);
}

/* Read the data portion of an intermediate Lucas-Lehmer results file */

int convertOldStyleLLSaveFile (
	llhandle *lldata,
	int	fd,
	unsigned long *counter,
	unsigned long *error_count)
{
	unsigned long i, buggy_error_count;
	unsigned long fftlen;
	unsigned long sum, filesum;
	int	bits, zero;		/* Guard against a zeroed out file */
	unsigned short type;

	_lseek (fd, 0, SEEK_SET);
	if (! read_short (fd, (short *) &type)) goto err;
	if (! read_long (fd, counter, NULL)) goto err;

/* Check for corrupt LL continuation files. */
/* Type 2 files are factoring continuation files */
/* Type 3 files are obsolete Advanced / Factoring continuation files */
/* Type 4 files are Advanced / Factoring continuation files */

	if (type <= 7) goto err;

/* Deduce the fftlen from the type field */

	if (type & 1)
		fftlen = (unsigned long) type - 1;
	else	
		fftlen = (unsigned long) type * 1024;

/* Handle case where the save file was for a different FFT length than */
/* we would prefer to use. */

	if (fftlen != gwfftlen (&lldata->gwdata)) {
		OutputStr (MAIN_THREAD_NUM, "FFT length mismatch in old LL save file\n");
		goto err;
	}

/* Read the fft data */

	bits = (int) lldata->gwdata.BITS_PER_WORD + 1;
	sum = 0;
	zero = TRUE;
	for (i = 0; i < gwfftlen (&lldata->gwdata); i++) {
		long	x;
		if (bits <= 15) {
			short y;
			if (! read_short (fd, &y)) goto err;
			x = (long) y;
		} else {
			if (! read_slong (fd, &x, NULL)) goto err;
		}
		if ((x & 0xFF000000) != 0 && (x & 0xFF000000) != 0xFF000000)
			goto err;
		sum = (uint32_t) (sum + x);
		if (x) zero = FALSE;
		set_fft_value (&lldata->gwdata, lldata->lldata, i, x);
	}
	if (!read_long (fd, &filesum, NULL)) goto err;
	if (!read_long (fd, &lldata->units_bit, &sum)) goto err;
	if (!read_long (fd, &buggy_error_count, &sum)) goto err;
	/* Now read in the correct V19 error count */
	if (!read_long (fd, error_count, &sum)) goto err;
	/* Kludge so that buggy v17 save files are rejected */
	/* V18 and later flip the bottom checksum bit */
	if (lldata->units_bit != 0) sum ^= 0x1;
	if (filesum != sum) goto err;
	if (zero) goto err;
	if (lldata->units_bit >= lldata->gwdata.n) goto err;
	_close (fd);
	return (TRUE);
err:	_close (fd);
	return (FALSE);
}

/* Write intermediate Lucas-Lehmer results to a file */
/* The LL save file format is: */
/*	52-bytes	standard header for all work types */
/*	u32		error_count */
/*	u32		iteration counter */
/*	u32		shift_count */
/*	gwnum		FFT data (u32 len, array u32s) */

#define LL_MAGICNUM		0x2c7330a8
#define LL_VERSION		1
#define LL_ERROR_COUNT_OFFSET	52

int writeLLSaveFile (
	llhandle *lldata,
	char	*filename,
	struct work_unit *w,
	unsigned long counter,
	unsigned long error_count)
{
	int	fd;
	unsigned long sum = 0;

/* If we are allowed to create multiple intermediate files, then */
/* write to a file called rXXXXXXX. */

	if (TWO_BACKUP_FILES && strlen (filename) == 8)
		filename[0] = 'r';

/* Now save to the intermediate file */

	fd = _open (filename, _O_BINARY | _O_WRONLY | _O_TRUNC | _O_CREAT, CREATE_FILE_ACCESS);
	if (fd < 0) return (FALSE);

	if (!write_header (fd, LL_MAGICNUM, LL_VERSION, w)) goto err;

	if (!write_long (fd, error_count, &sum)) goto err;
	if (!write_long (fd, counter, &sum)) goto err;
	if (!write_long (fd, lldata->units_bit, &sum)) goto err;
	if (!write_gwnum (fd, &lldata->gwdata, lldata->lldata, &sum)) goto err;

	if (!write_checksum (fd, sum)) goto err;

	_commit (fd);
	_close (fd);

/* Now rename the intermediate files */

	if (TWO_BACKUP_FILES && strlen (filename) == 8) {
		char	backupname[16];
		strcpy (backupname, filename);
		backupname[0] = 'q'; filename[0] = 'p';
		_unlink (backupname);
		 rename (filename, backupname);
		backupname[0] = 'r';
		rename (backupname, filename);
	}

	return (TRUE);

/* An error occured.  Delete the current file. */

err:	_close (fd);
	_unlink (filename);
	return (FALSE);
}

/* Update the error count in an intermediate file */

void writeNewErrorCount (
	char	*filename,
	unsigned long new_error_count)
{
	int	fd;
	unsigned long sum, old_error_count;

/* Open the intermediate file, position past the FFT data */

	fd = _open (filename, _O_BINARY | _O_RDWR);
	if (fd < 0) return;

/* Read in the checksum and old error count */

	if (!read_checksum (fd, &sum)) goto err;
	_lseek (fd, LL_ERROR_COUNT_OFFSET, SEEK_SET);
	if (!read_long (fd, &old_error_count, NULL)) goto err;

/* Update the checksum */

	sum = sum - old_error_count + new_error_count;

/* Write out the checksum and new error count */

	if (!write_checksum (fd, sum)) goto err;
	_lseek (fd, LL_ERROR_COUNT_OFFSET, SEEK_SET);
	if (!write_long (fd, new_error_count, NULL)) goto err;

/* Close file and return */

err:	_close (fd);
}

/* Read the data portion of an intermediate Lucas-Lehmer results file */

int readLLSaveFile (
	llhandle *lldata,
	char	*filename,
	struct work_unit *w,
	unsigned long *counter,
	unsigned long *error_count)
{
	int	fd;
	unsigned long sum, filesum, version;

	fd = _open (filename, _O_BINARY | _O_RDONLY);
	if (fd <= 0) return (FALSE);

	if (!read_magicnum (fd, LL_MAGICNUM))
		return (convertOldStyleLLSaveFile (lldata, fd, counter,
						   error_count));
	if (!read_header (fd, &version, w, &filesum)) goto err;
	if (version != LL_VERSION) goto err;

	sum = 0;
	if (!read_long (fd, error_count, &sum)) goto err;
	if (!read_long (fd, counter, &sum)) goto err;
	if (!read_long (fd, &lldata->units_bit, &sum)) goto err;
	if (lldata->units_bit >= lldata->gwdata.n) goto err;

	if (!read_gwnum (fd, &lldata->gwdata, lldata->lldata, &sum)) goto err;

	if (filesum != sum) goto err;
	_close (fd);
	return (TRUE);
err:	_close (fd);
	return (FALSE);
}

/* Increment the error counter.  The error counter is one 32-bit */
/* field that contains 5 values - a flag if this is a contiuation */
/* from a save file that did not track error counts, a count of */
/* errors that were reproducible, a count of ILLEAL SUMOUTs, */
/* a count of convolution errors above 0.4, and a count of */
/* SUMOUTs not close enough to SUMINPs. */

void inc_error_count (
	int	type,
	unsigned long *error_count)
{
	unsigned long addin, maxval, temp;
	
	addin = 1 << (type * 8);
	maxval = ((type == 3) ? 127 : 255) * addin;
	temp = *error_count & maxval;
	if (temp != maxval) temp += addin;
	*error_count = (*error_count & ~maxval) + temp;
}

/* Prepare for subtracting 2 from the squared result.  Also keep track */
/* of the location of the ever changing units bit. */

void lucas_fixup (
	llhandle *lldata,
	unsigned long p)	/* Exponent being tested */
{

/* We are about to square the number, the units bit position will double */

	lldata->units_bit <<= 1;
	if (lldata->units_bit >= p) lldata->units_bit -= p;

/* Tell gwnum code the value to subtract 2 from the squared result. */

	gwsetaddinatbit (&lldata->gwdata, -2, lldata->units_bit);
}

/* Generate random FFT data for timing the Lucas-Lehmer code */

void generateRandomData (
	llhandle *lldata)
{
	unsigned long i;

/* Fill data space with random values. */

	srand ((unsigned) time (NULL));
	for (i = 0; i < gwfftlen (&lldata->gwdata); i++) {
		set_fft_value (&lldata->gwdata, lldata->lldata, i, rand() & 0xFF);
	}
}

/* For exponents that are near an FFT limit, do 1000 sample iterations */
/* to see if we should use the smaller or larger FFT size.  We examine */
/* the average roundoff error to determine which FFT size to use. */

int pick_fft_size (
	int	thread_num,
	struct work_unit *w)
{
	llhandle lldata;
	char	buf[120];
	double	softpct, total_error, avg_error, max_avg_error;
	unsigned long small_fftlen, large_fftlen;
	int	i, stop_reason;

/* We only do this for Mersenne numbers */

	if (w->k != 1.0 || w->b != 2 || w->c != -1) return (0);

/* We don't do this for small exponents.  We've not studied the average */
/* error enough on smaller FFT sizes to intelligently pick the FFT size. */
/* Also, for really large exponents there is no larger FFT size to use! */

	if (w->n <= 5000000) return (0);

/* If we've already calculated the best FFT size, then return */

	if (w->forced_fftlen) return (0);

/* If we have an old-style (v24) entry for this exponent in local.ini, then */
/* use it.  V25 and later store the fft length in the worktodo.ini file. */

	IniGetString (LOCALINI_FILE, "SoftCrossoverData", buf, sizeof (buf), "0");
	if (w->n == (unsigned long) atol (buf)) {
		char	*comma;
		unsigned long fftlen;
		comma = strchr (buf, ',');
		if (comma != NULL) {
			*comma++ = 0;
			fftlen = atol (comma);
			comma = strchr (comma, ',');
			if (comma != NULL) {
				unsigned long sse2;
				*comma++ = 0;
				sse2 = atol (comma);
				if ((sse2 && (CPU_FLAGS & CPU_SSE2)) ||
				    (!sse2 && !(CPU_FLAGS & CPU_SSE2))) {
					w->forced_fftlen = fftlen;
					stop_reason = updateWorkToDoLine (thread_num, w);
					if (stop_reason) return (stop_reason);
				}
			}
		}
	}

/* Get the info on how what percentage of exponents on either side of */
/* an FFT crossover we will do this 1000 iteration test. */

	IniGetString (INI_FILE, "SoftCrossover", buf, sizeof (buf), "0.2");
	softpct = atof (buf) / 100.0;

/* If this exponent is not close to an FFT crossover, then we are done */

	small_fftlen = gwmap_to_fftlen (1.0, 2,
			(unsigned long) ((1.0 - softpct) * w->n), -1);
	large_fftlen = gwmap_to_fftlen (1.0, 2,
			(unsigned long) ((1.0 + softpct) * w->n), -1);
	if (small_fftlen == large_fftlen || large_fftlen == 0) return (0);

/* Let the user be more conservative or more aggressive in picking the */
/* acceptable average error.  By default, we accept an average error */
/* between 0.241 and 0.243 depending on the FFT size. */

	max_avg_error = 0.241 + 0.002 *
		(log ((double) small_fftlen) - log ((double) 262144.0)) /
		(log ((double) 4194304.0) - log ((double) 262144.0));
	IniGetString (INI_FILE, "SoftCrossoverAdjust", buf, sizeof (buf), "0");
	max_avg_error += atof (buf);

/* Print message to let user know what is going on */

	sprintf (buf,
		 "Trying 1000 iterations for exponent %ld using %dK FFT.\n",
		 w->n, small_fftlen / 1024);
	OutputBoth (thread_num, buf);
	sprintf (buf,
		 "If average roundoff error is above %.5g, then a larger FFT will be used.\n",
		 max_avg_error);
	OutputBoth (thread_num, buf);

/* Init the FFT code using the smaller FFT size */

	gwinit (&lldata.gwdata);
	stop_reason = lucasSetup (thread_num, w->n, small_fftlen, &lldata);
	if (stop_reason) return (stop_reason);

/* Fill data space with random values then do one squaring to make */
/* the data truly random. */

	generateRandomData (&lldata);
	gwsetnormroutine (&lldata.gwdata, 0, TRUE, 0);
	gwstartnextfft (&lldata.gwdata, TRUE);
	gwsquare (&lldata.gwdata, lldata.lldata);

/* Average the roundoff error over a 1000 iterations. */

	for (i = 0, total_error = 0.0; ; ) {
		gw_clear_maxerr (&lldata.gwdata);
		gwsquare (&lldata.gwdata, lldata.lldata);
		total_error += gw_get_maxerr (&lldata.gwdata);
		stop_reason = stopCheck (thread_num);
		if (stop_reason) {
			lucasDone (&lldata);
			return (stop_reason);
		}
		if (++i == 1000) break;
		if (i % 100 == 0) {
			sprintf (buf,
				 "After %d iterations average roundoff error is %.5g.\n",
				 i, total_error / (double) i);
			OutputStr (thread_num, buf);
		}
	}
	avg_error = total_error / 1000.0;
	lucasDone (&lldata);

/* Now decide which FFT size to use based on the average error. */
/* Save this info in worktodo.ini so that we don't need to do this again. */

	w->forced_fftlen = (avg_error <= max_avg_error) ? small_fftlen : large_fftlen;
	stop_reason = updateWorkToDoLine (thread_num, w);
	if (stop_reason) return (stop_reason);

/* Output message to user informing him of the outcome. */

	sprintf (buf,
		 "Final average roundoff error is %.5g, using %dK FFT for exponent %ld.\n",
		 avg_error, w->forced_fftlen / 1024, w->n);
	OutputBoth (thread_num, buf);
	return (0);
}

/* Test if we are near the maximum exponent this fft length can test */
/* We only support this (careful iterations when near fft limit) for */
/* Mersenne numbers. */

int exponent_near_fft_limit (
	gwhandle *gwdata)		/* Handle returned by gwsetup */
{
	char	pct[30];
	IniGetString (INI_FILE, "NearFFTLimitPct", pct, sizeof(pct), "0.5");
	return (gwnear_fft_limit (gwdata, atof (pct)));
}

/* Do an LL iteration very carefully.  This is done after a normal */
/* iteration gets a roundoff error above 0.40.  This careful iteration */
/* will not generate a roundoff error. */

void careful_iteration (
	llhandle *lldata,		/* Handle from lucasSetup */
	unsigned long p)		/* Exponent being tested */
{
	gwnum	hi, lo;
	unsigned long i;

/* Copy the data to hi and lo.  Zero out half the FFT data in each. */

	hi = gwalloc (&lldata->gwdata);
	lo = gwalloc (&lldata->gwdata);
	gwcopy (&lldata->gwdata, lldata->lldata, hi);
	gwcopy (&lldata->gwdata, lldata->lldata, lo);
	for (i = 0; i < gwfftlen (&lldata->gwdata)/2; i++)
		set_fft_value (&lldata->gwdata, hi, i, 0);
	for ( ; i < gwfftlen (&lldata->gwdata); i++)
		set_fft_value (&lldata->gwdata, lo, i, 0);

/* Now do the squaring using three multiplies and adds */

	gwsetnormroutine (&lldata->gwdata, 0, 0, 0);
	gwstartnextfft (&lldata->gwdata, FALSE);
	gwsetaddin (&lldata->gwdata, 0);
	gwfft (&lldata->gwdata, hi, hi);
	gwfft (&lldata->gwdata, lo, lo);
	gwfftfftmul (&lldata->gwdata, lo, hi, lldata->lldata);
	gwfftfftmul (&lldata->gwdata, hi, hi, hi);
	lucas_fixup (lldata, p);
	gwfftfftmul (&lldata->gwdata, lo, lo, lo);
	gwaddquick (&lldata->gwdata, lldata->lldata, lldata->lldata);
	gwaddquick (&lldata->gwdata, hi, lldata->lldata);
	gwadd (&lldata->gwdata, lo, lldata->lldata);

/* Since our error recovery code cannot cope with an error during a careful */
/* iteration, make sure the error variable is cleared.  This shouldn't */
/* ever happen, but two users inexplicably ran into this problem. */

	gw_clear_error (&lldata->gwdata);

/* Free memory and return */

	gwfree (&lldata->gwdata, hi);
	gwfree (&lldata->gwdata, lo);
}

/* Output the good news of a new prime to the screen in an infinite loop */

void good_news (void *arg)
{
	char	buf[80];

	title (MAIN_THREAD_NUM, "New Prime!!!");
	sprintf (buf, "New Mersenne Prime!!!!  M%d is prime!\n", (int) (intptr_t) arg);
	while (WORKER_THREADS_ACTIVE && ! WORKER_THREADS_STOPPING) {
		OutputStr (MAIN_THREAD_NUM, buf);
		flashWindowAndBeep ();
		Sleep (50);
	}
}

/* Do the Lucas-Lehmer test */

int prime (
	int	thread_num,		/* Worker thread number */
	struct PriorityInfo *sp_info,	/* SetPriority information */
	struct work_unit *w,		/* Worktodo entry */
	int	pass)			/* PrimeContinue pass */
{
	llhandle lldata;
	unsigned long p;
	unsigned long counter;
	unsigned long error_count;
	unsigned long iters;
	char	filename[32];
	double	timers[2];
	double	inverse_p;
	double	reallyminerr = 1.0;
	double	reallymaxerr = 0.0;
	double	*addr1;
	int	first_iter_msg, saving, near_fft_limit, sleep5;
	unsigned long high32, low32;
	int	rc, isPrime, stop_reason;
	char	buf[160], fft_desc[100];
	int	slow_iteration_count;
	double	best_iteration_time;
	unsigned long last_counter = 0;		/* Iteration of last error */
	int	maxerr_recovery_mode = 0;	/* Big roundoff err rerun */
	double	last_suminp = 0.0;
	double	last_sumout = 0.0;
	double	last_maxerr = 0.0;

/* Initialize */

	p = w->n;

/* Do some of the trial factoring.  We treat factoring that is part of a */
/* LL test as priority work (done in pass 1).  We don't do all the trial */
/* factoring as the last 2 bit levels take a lot of time and are unlikely */
/* to find a factor.  The P-1 test will probably be needed anyway and */
/* may find a factor thus saving us from doing the last 2 bit levels. */

	if ((w->work_type == WORK_TEST || w->work_type == WORK_DBLCHK) &&
	    ! IniGetInt (INI_FILE, "SkipTrialFactoring", 0)) {
		stop_reason = primeFactor (thread_num, sp_info, w, 2);
		if (stop_reason) return (stop_reason);
	}

/* See if this exponent needs P-1 factoring.  We treat P-1 factoring */
/* that is part of an LL test as priority work (done in pass 1). */

	if ((w->work_type == WORK_TEST || w->work_type == WORK_DBLCHK) &&
	    ! w->pminus1ed) {
		stop_reason = pfactor (thread_num, sp_info, w, pass);
		if (stop_reason) {
			if (pass == 3 && stop_reason == STOP_NOT_ENOUGH_MEM)
				stop_reason = 0;
			else
				return (stop_reason);
		}
	}

/* Do the rest of the trial factoring. */

	if ((w->work_type == WORK_TEST || w->work_type == WORK_DBLCHK) &&
	    ! IniGetInt (INI_FILE, "SkipTrialFactoring", 0)) {
		stop_reason = primeFactor (thread_num, sp_info, w, 0);
		if (stop_reason) return (stop_reason);
	}

/* Done with pass 1 priority work.  Return to do more priority work. */

	if (pass == 1 && w->work_type != WORK_ADVANCEDTEST) return (0);

/* Figure out which FFT size we should use */

	stop_reason = pick_fft_size (thread_num, w);
	if (stop_reason) return (stop_reason);

/* Make sure the first-time user runs a successful self-test. */
/* The one-hour self-test may have been useful when it was first introduced */
/* but I think it now does little to catch buggy machines (they eventually */
/* work OK for an hour) and does create user confusion and annoyance. */

#ifdef ONE_HOUR_SELF_TEST
	if (w->work_type == WORK_TEST || w->work_type == WORK_DBLCHK) {
		stop_reason = selfTest (thread_num, sp_info, w);
		if (stop_reason) return (stop_reason);
	}
#endif

/* Loop reading from save files (and backup save files) */

readloop:
	tempFileName (w, filename);

/* Setup the LL test */

	gwinit (&lldata.gwdata);
	gwset_num_threads (&lldata.gwdata, THREADS_PER_TEST[thread_num]);
	gwset_thread_callback (&lldata.gwdata, SetAuxThreadPriority);
	gwset_thread_callback_data (&lldata.gwdata, sp_info);
	stop_reason = lucasSetup (thread_num, p, w->forced_fftlen, &lldata);
	if (stop_reason) return (stop_reason);

/* Record the amount of memory being used by this thread. */

	set_memory_usage (thread_num, 0,
			  gwmemused (&lldata.gwdata) + gwnum_size (&lldata.gwdata));

/* Read an LL save file.  On error try the backup intermediate file. */

	if (continuationFileExists (thread_num, filename)) {
		if (! readLLSaveFile (&lldata, filename, w, &counter, &error_count) ||
		    counter > w->n) {
			lucasDone (&lldata);
			sprintf (buf, READFILEERR, filename);
			OutputBoth (thread_num, buf);
			_unlink (filename);
			goto readloop;
		}
	}

/* Start off with the 1st Lucas number */

	else {
		counter = 2;
		error_count = 0;
	}

/* Hyperthreading backoff is an option to pause the program when iterations */
/* take longer than usual.  This is useful on hyperthreaded machines so */
/* that prime95 doesn't steal cycles from a foreground task, thus hurting */
/* the computers responsiveness. */

	best_iteration_time = 1.0e50;
	slow_iteration_count = 0;

/* Clear all timers */

	clear_timers (timers, sizeof (timers) / sizeof (timers[0]));

/* Init the title */

	sprintf (buf, "%ld / %ld", counter, p);
	title (thread_num, buf);

/* Init vars for Test/Status and CommunicateWithServer */

	strcpy (w->stage, "LL");
	inverse_p = 1.0 / (double) p;
	w->pct_complete = (double) counter * inverse_p;

/* Start off with the 1st Lucas number - four */
/* Note we do something a little strange here.  We actually set the */
/* first number to 4 but shifted by a random amount.  This lets two */
/* different machines check the same Mersenne number and operate */
/* on different FFT data - thus greatly reducing the chance that */
/* a CPU or program error corrupts the results. */

	if (counter == 2) {
		unsigned long i, word, bit_in_word;
		uint32_t hi, lo;
		srand ((unsigned) time (NULL));
		lldata.units_bit = (rand () << 16) + rand ();
		if (CPU_FLAGS & CPU_RDTSC) { rdtsc(&hi,&lo); lldata.units_bit += lo; }
		lldata.units_bit = lldata.units_bit % p;
		bitaddr (&lldata.gwdata, (lldata.units_bit + 2) % p, &word, &bit_in_word);
		for (i = 0; i < gwfftlen (&lldata.gwdata); i++) {
			set_fft_value (&lldata.gwdata, lldata.lldata, i, (i == word) ? (1L << bit_in_word) : 0);
		}
		first_iter_msg = FALSE;
	} else
		first_iter_msg = TRUE;

/* Output a message indicating we are starting/resuming an LL test. */
/* Also tell user the FFT length. */

	gwfft_description (&lldata.gwdata, fft_desc);
	sprintf (buf, "%s primality test of M%ld using %s\n",
		 (counter == 2) ? "Starting" : "Resuming", p, fft_desc);
	OutputStr (thread_num, buf);

/* If we are near the maximum exponent this fft length can test, then we */
/* will error check all iterations */

	near_fft_limit = exponent_near_fft_limit (&lldata.gwdata);

/* Get address of second FFT data element.  We'll use this for very */
/* quickly checking for zeroed FFT data. */

	addr1 = addr (&lldata.gwdata, lldata.lldata, 1);

/* Compute numbers in the lucas series, write out every 30 minutes to a file */

	iters = 0;
	while (counter < p) {
		int	echk;

/* On first iteration create a save file so that writeNewErrorCount */
/* can properly keep track of error counts. */
/* Also save right after we pass an errored iteration and several */
/* iterations before retesting an errored iteration so that we don't */
/* have to backtrack very far to do a careful_iteration	(we don't do the */
/* iteration immediately before because on the P4 a save operation will */
/* change the FFT data and make the error non-reproducible. */
/* Error check the last 50 iterations, before writing an */
/* intermediate file (either user-requested stop or a */
/* 30 minute interval expired), and every 128th iteration. */

		stop_reason = stopCheck (thread_num);
		saving = stop_reason ||
			 (counter == 2 && p > 1500000) ||
			 counter == last_counter-8 ||
			 counter == last_counter ||
			 testSaveFilesFlag (thread_num);
		echk = saving || near_fft_limit || ERRCHK ||
			(counter >= p - 50) || ((counter & 127) == 0);
		gw_clear_maxerr (&lldata.gwdata);

/* Do a Lucas-Lehmer iteration */

		timers[1] = 0.0;
		start_timer (timers, 1);

/* If we are recovering from a big roundoff error, then run one */
/* iteration using three multiplies where half the data is zeroed. */
/* This won't run into any roundoff problems and will protect us from */
/* roundoff errors up to 0.6. */

		if (maxerr_recovery_mode && counter == last_counter) {
			careful_iteration (&lldata, p);
			maxerr_recovery_mode = 0;
			echk = 0;
		}

/* Otherwise, do a normal iteration */

#ifndef SERVER_TESTING
		else {
			gwsetnormroutine (&lldata.gwdata, 0, echk, 0);
			gwstartnextfft (&lldata.gwdata,
					!saving && !maxerr_recovery_mode &&
					counter+1 != p &&
					(INTERIM_FILES == 0 ||
					 (counter+1) % INTERIM_FILES > 0) &&
					(INTERIM_RESIDUES == 0 ||
					 (counter+1) % INTERIM_RESIDUES > 2));
			lucas_fixup (&lldata, p);
			gwsquare (&lldata.gwdata, lldata.lldata);
		}
#endif

/* End iteration timing and increase count of iterations completed */

		end_timer (timers, 1);
		timers[0] += timers[1];
		iters++;

/* If the sum of the output values is an error (such as infinity) */
/* then raise an error. */

		if (gw_test_illegal_sumout (&lldata.gwdata)) {
			sprintf (buf, ERRMSG0, counter, p, ERRMSG1A);
			OutputBoth (thread_num, buf);
			inc_error_count (2, &error_count);
			sleep5 = TRUE;
			goto restart;
		}

/* Check that the sum of the input numbers squared is approximately */
/* equal to the sum of unfft results.  Since this check may not */
/* be perfect, check for identical results after a restart. */

		if (gw_test_mismatched_sums (&lldata.gwdata)) {
			if (counter == last_counter &&
			    gwsuminp (&lldata.gwdata, lldata.lldata) == last_suminp &&
			    gwsumout (&lldata.gwdata, lldata.lldata) == last_sumout) {
				OutputBoth (thread_num, ERROK);
				inc_error_count (3, &error_count);
				gw_clear_error (&lldata.gwdata);
			} else {
				char	msg[100];
				sprintf (msg, ERRMSG1B,
					 gwsuminp (&lldata.gwdata, lldata.lldata),
					 gwsumout (&lldata.gwdata, lldata.lldata));
				sprintf (buf, ERRMSG0, counter, p, msg);
				OutputBoth (thread_num, buf);
				last_counter = counter;
				last_suminp = gwsuminp (&lldata.gwdata, lldata.lldata);
				last_sumout = gwsumout (&lldata.gwdata, lldata.lldata);
				inc_error_count (0, &error_count);
				sleep5 = TRUE;
				goto restart;
			}
		}

/* Check for excessive roundoff error.  If round off is too large, repeat */
/* the iteration to see if this was a hardware error.  If it was repeatable */
/* then repeat the iteration using a safer, slower method.  This can */
/* happen when operating near the limit of an FFT. */

		if (echk && gw_get_maxerr (&lldata.gwdata) >= 0.40625) {
			if (counter == last_counter &&
			    gw_get_maxerr (&lldata.gwdata) == last_maxerr) {
				OutputBoth (thread_num, ERROK);
				inc_error_count (3, &error_count);
				gw_clear_error (&lldata.gwdata);
				OutputBoth (thread_num, ERRMSG5);
				maxerr_recovery_mode = 1;
				sleep5 = FALSE;
				goto restart;
			} else {
				char	msg[100];
				sprintf (msg, ERRMSG1C, gw_get_maxerr (&lldata.gwdata));
				sprintf (buf, ERRMSG0, counter, p, msg);
				OutputBoth (thread_num, buf);
				last_counter = counter;
				last_maxerr = gw_get_maxerr (&lldata.gwdata);
				inc_error_count (1, &error_count);
				sleep5 = FALSE;
				goto restart;
			}
		}

/* Check if the units_bit is corrupt.  This will make sure we are always */
/* subtracting 2 from the FFT data.  If the FFT data was mysteriously zeroed */
/* and the units_bit value was corrupt then we could get a false positive */
/* result.  With this fix we should get into a safe -2, 2, 2, 2 loop. */

		if (lldata.units_bit >= p) {
			sprintf (buf, ERRMSG0, counter, p, ERRMSG1D);
			OutputBoth (thread_num, buf);
			inc_error_count (2, &error_count);
			sleep5 = TRUE;
			goto restart;
		}

/* Check if the FFT data has been zeroed. This will help reduce the chances */
/* of another false positive being reported. */

#ifndef SERVER_TESTING
		if (*addr1 == 0.0 && p > 1000 &&
		    counter > 50 && counter < p-2 && counter != last_counter) {
			unsigned long i;		
			for (i = 2; ; i++) {
				if (*addr (&lldata.gwdata, lldata.lldata, i) != 0.0) break;
				if (i == 50) {
					sprintf (buf, ERRMSG0, counter, p, ERRMSG1F);
					OutputBoth (thread_num, buf);
					inc_error_count (2, &error_count);
					last_counter = counter;
					sleep5 = TRUE;
					goto restart;
				}
			}
		}
#endif

/* Update counter, percentage complete, and maximum round-off error */

		counter++;
		w->pct_complete = (double) counter * inverse_p;
		if (ERRCHK) {
			if (gw_get_maxerr (&lldata.gwdata) < reallyminerr && counter > 30)
				reallyminerr = gw_get_maxerr (&lldata.gwdata);
			if (gw_get_maxerr (&lldata.gwdata) > reallymaxerr)
				reallymaxerr = gw_get_maxerr (&lldata.gwdata);
		}

/* Print a message every so often */

		if (counter % ITER_OUTPUT == 0 || first_iter_msg) {
			char	fmt_mask[80];
			double	pct;
			pct = trunc_percent (w->pct_complete);
			sprintf (fmt_mask, "%%.%df%%%% of M%%ld", PRECISION);
			sprintf (buf, fmt_mask, pct, p);
			title (thread_num, buf);
			sprintf (fmt_mask,
				 "Iteration: %%ld / %%ld [%%.%df%%%%]",
				 PRECISION);
			sprintf (buf, fmt_mask, counter, p, pct);
			if (ERRCHK && counter > 30) {
				sprintf (buf+strlen(buf),
					 ".  Round off: %10.10f to %10.10f",
					 reallyminerr, reallymaxerr);
			}
			if (first_iter_msg) {
				strcat (buf, ".\n");
				clear_timer (timers, 0);
			} else {
				strcat (buf, ".  Per iteration time: ");
				divide_timer (timers, 0, iters);
				print_timer (timers, 0, buf,
					     TIMER_NL | TIMER_OPT_CLR);
			}
			OutputStr (thread_num, buf);
			if (!CUMULATIVE_TIMING) iters = 0;
			first_iter_msg = FALSE;
		}

/* Print a results file message every so often */

		if (counter % ITER_OUTPUT_RES == 0 || (NO_GUI && stop_reason)) {
			sprintf (buf, "Iteration %ld / %ld\n", counter, p);
			writeResults (buf);
		}

/* Write results to a file every DISK_WRITE_TIME minutes */
/* On error, retry in 10 minutes (it could be a temporary */
/* disk-full situation) */

		if (saving) {
			if (! writeLLSaveFile (&lldata, filename, w, counter,
					       error_count)) {
				sprintf (buf, WRITEFILEERR, filename);
				OutputBoth (thread_num, buf);
			}
		}

/* If an escape key was hit, write out the results and return */

		if (stop_reason) {
			char	fmt_mask[80];
			sprintf (fmt_mask,
				 "Stopping primality test of M%%ld at iteration %%ld [%%.%df%%%%]\n",
				 PRECISION);
			sprintf (buf, fmt_mask, p, counter,
				 trunc_percent (w->pct_complete));
			OutputStr (thread_num, buf);
			lucasDone (&lldata);
			return (stop_reason);
		}

/* Output the 64-bit residue at specified interims.  Also output the */
/* residues for the next two iterations so that we can compare our */
/* residues to programs that start counter at zero or one. */

		if (INTERIM_RESIDUES && counter % INTERIM_RESIDUES <= 2) {
			generateResidue64 (&lldata, &high32, &low32);
			sprintf (buf, 
				 "M%ld interim Wd%d residue %08lX%08lX at iteration %ld\n",
				 p, PORT, high32, low32, counter);
			OutputBoth (thread_num, buf);
		}

/* Write a save file every INTERIM_FILES iterations. */

		if (INTERIM_FILES && counter % INTERIM_FILES == 0) {
			char	interimfile[20];
			sprintf (interimfile, "%.8s.%03d",
				 filename, counter / INTERIM_FILES);
			writeLLSaveFile (&lldata, interimfile, w, counter,
				         error_count);
		}

/* If ten iterations take 40% longer than a typical iteration, then */
/* assume a foreground process is running and sleep for a short time */
/* to give the foreground process more CPU time.  Even though a foreground */
/* process runs at higher priority, hyperthreading will cause this */
/* program to run at an equal priority, hurting responsiveness. */

		if (HYPERTHREADING_BACKOFF && p > 10000000) {
			if (timers[1] < best_iteration_time)
				best_iteration_time = timers[1];
			if (timers[1] > 1.40 * best_iteration_time) {
				if (slow_iteration_count == 10) {
					sprintf (buf, "Pausing %d seconds.\n",
						 HYPERTHREADING_BACKOFF);
					OutputStr (thread_num, buf);
					Sleep (HYPERTHREADING_BACKOFF * 1000);
				}
				slow_iteration_count++;
			} else
				slow_iteration_count = 0;
		}
	}

/* Check for a successful completion */
/* We found a prime if result is zero */
/* Note that all values of -1 is the same as zero */

	rc = generateResidue64 (&lldata, &high32, &low32);
	if (rc < 0) {
		sprintf (buf, ERRMSG0, counter, p, ERRMSG1E);
		OutputBoth (thread_num, buf);
		inc_error_count (2, &error_count);
		sleep5 = TRUE;
		goto restart;
	}
	isPrime = (rc == 0);

/* Format the output message */

	if (isPrime)
		sprintf (buf, "M%ld is prime! Wd%d: %08lX,%08lX\n",
			 p, PORT, SEC1 (p), error_count);
	else
		sprintf (buf,
			 "M%ld is not prime. Res64: %08lX%08lX. Wd%d: %08lX,%ld,%08lX\n",
			 p, high32, low32, PORT,
			 SEC2 (p, high32, low32, lldata.units_bit, error_count),
			 lldata.units_bit, error_count);
	OutputStr (thread_num, buf);
	formatMsgForResultsFile (buf, w);
	rc = writeResults (buf);

/* Output results to the screen, results file, and server */

	{
		struct primenetAssignmentResult pkt;
		memset (&pkt, 0, sizeof (pkt));
		strcpy (pkt.computer_guid, COMPUTER_GUID);
		strcpy (pkt.assignment_uid, w->assignment_uid);
		strcpy (pkt.message, buf);
		pkt.result_type =
			isPrime ? PRIMENET_AR_LL_PRIME : PRIMENET_AR_LL_RESULT;
		pkt.n = p;
		sprintf (pkt.residue, "%08lX%08lX", high32, low32);
		pkt.shift_count = lldata.units_bit;
		sprintf (pkt.error_count, "%08lX", error_count);
		pkt.fftlen = gwfftlen (&lldata.gwdata);
		pkt.done = TRUE;
		spoolMessage (PRIMENET_ASSIGNMENT_RESULT, &pkt);
	}

/* Delete the continuation files - assuming the results file write */
/* was successful. */

	if (!isPrime || isKnownMersennePrime (p)) {
		if (rc) _unlink (filename);
		filename[0] = 'q';
		_unlink (filename);
	}

/* Clean up */

	lucasDone (&lldata);

/* Output good news to the screen in an infinite loop */

	if (isPrime && !SILENT_VICTORY && !isKnownMersennePrime (p)) {
		gwthread thread_handle;
		gwthread_create (&thread_handle, &good_news, (void *) p);
	}

/* All done */

	return (STOP_WORK_UNIT_COMPLETE);

/* An error occured, output a message saying we are restarting, sleep, */
/* then try restarting at last save point. */

restart:if (sleep5) OutputBoth (thread_num, ERRMSG2);
	OutputBoth (thread_num, ERRMSG3);

/* Update the error count in the save file */

	writeNewErrorCount (filename, error_count);

/* Sleep five minutes before restarting */

	if (sleep5) {
		stop_reason = SleepFive (thread_num);
		if (stop_reason) return (stop_reason);
	}

/* Return so that last continuation file is read in */

	lucasDone (&lldata);
	goto readloop;
}

/*********************/
/* Torture test code */
/*********************/

#define TORTURE1 "Beginning a continuous self-test to check your computer.\n"
#if defined (__linux__) || defined (__FreeBSD__) || defined (__EMX__)
#define TORTURE2 "Please read stress.txt.  Hit ^C to end this test.\n"
#else
#define TORTURE2 "Please read stress.txt.  Choose Test/Stop to end this test.\n"
#endif
#define SELFMSG1A "The program will now perform a self-test to make sure the\n"
#define SELFMSG1B "Lucas-Lehmer code is working properly on your computer.\n"
#define SELFMSG1C "This will take about an hour.\n"
#define SELF1 "Test %i, %i Lucas-Lehmer iterations of M%ld using %s.\n"
#define SELFFAIL "FATAL ERROR: Final result was %08lX, expected: %08lX.\n"
char SELFFAIL1[] = "ERROR: ILLEGAL SUMOUT\n";
char SELFFAIL2[] = "FATAL ERROR: Resulting sum was %.16g, expected: %.16g\n";
char SELFFAIL3[] = "FATAL ERROR: Rounding was %.10g, expected less than 0.4\n";
char SELFFAIL4[] = "Possible hardware failure, consult readme.txt file, restarting test.\n";
char SELFFAIL5[] = "Hardware failure detected, consult stress.txt file.\n";
char SELFFAIL6[] = "Maximum number of warnings exceeded.\n";

#define SELFPASS "Self-test %iK passed!\n"
char SelfTestIniMask[] = "SelfTest%iPassed";

struct self_test_info {
	unsigned long p;
	unsigned long iters;
	unsigned long reshi;
};

#define MAX_SELF_TEST_ITERS	405
struct self_test_info SELF_TEST_DATA[MAX_SELF_TEST_ITERS] = {
{560000001, 100, 0x7F853A0A}, {420000001, 150, 0x89665E7E},
{280000001, 200, 0xC32CAD46}, {210000001, 300, 0x89823329},
{140000001, 400, 0x15EF4F24}, {110000001, 500, 0x893C9000},
{78643201, 400, 0x2D9C8904}, {78643199, 400, 0x7D469182},
{75497473, 400, 0x052C7FD8}, {75497471, 400, 0xCCE7495D},
{71303169, 400, 0x467A9338}, {71303167, 400, 0xBBF8B37D},
{68157441, 400, 0xBE71E616}, {68157439, 400, 0x93A71CC2},
{66060289, 400, 0xF296BB99}, {66060287, 400, 0x649EEF2A},
{62390273, 400, 0xBC8DFC27}, {62390271, 400, 0xDE7D5B5E},
{56623105, 400, 0x0AEBF972}, {56623103, 400, 0x1BA96297},
{53477377, 400, 0x5455F347}, {53477375, 400, 0xCE1C7F78},
{50331649, 400, 0x3D746AC8}, {50331647, 400, 0xE23F2DE6},
{49807361, 400, 0xB43EF4C5}, {49807359, 400, 0xA8BEB02D},
{47185921, 400, 0xD862563C}, {47185919, 400, 0x17281086},
{41943041, 400, 0x0EDA1F92}, {41943039, 400, 0xDE6911AE},
{39845889, 400, 0x43D8A96A}, {39845887, 400, 0x3D118E8F},
{37748737, 400, 0x38261154}, {37748735, 400, 0x22B34CD2},
{35651585, 400, 0xB0E48D2E}, {35651583, 400, 0xCC3340C6},
{34865153, 400, 0xD2C00E6C}, {34865151, 400, 0xFA644F69},
{33030145, 400, 0x83E5738D}, {33030143, 400, 0x6EDBC5B5},
{31195137, 400, 0xFF9591CF}, {31195135, 400, 0x04577C70},
{29884417, 400, 0xACC36457}, {29884415, 400, 0xC0FE7B1E},
{28311553, 400, 0x780EB8F5}, {28311551, 400, 0xE6D128C3},
{26738689, 400, 0x09DC45B0}, {26738687, 400, 0xDC7C074A},
{24903681, 400, 0xA482CF1E}, {24903679, 400, 0x4B3F5121},
{23592961, 400, 0xAFE3C198}, {23592959, 400, 0xCF9AD48C},
{20971521, 400, 0x304EC13B}, {20971519, 400, 0x9C4E157E},
{19922945, 400, 0x83FE36D9}, {19922943, 400, 0x9C60E7A2},
{18874369, 400, 0x83A9F8CB}, {18874367, 400, 0x5A6E22E0},
{17825793, 400, 0xF3A90A5E}, {17825791, 400, 0x6477CA76},
{17432577, 400, 0xCAB36E6A}, {17432575, 400, 0xB8F814C6},
{16515073, 400, 0x91EFCB1C}, {16515071, 400, 0xA0C35CD9},
{15597569, 400, 0x12E057AD}, {15597567, 400, 0xC4EFAEFD},
{14942209, 400, 0x1C912A7B}, {14942207, 400, 0xABA9EA6E},
{14155777, 400, 0x4A943A4E}, {14155775, 400, 0x00789FB9},
{13369345, 400, 0x27A041EE}, {13369343, 400, 0xA8B01A41},
{12451841, 400, 0x4DC891F6}, {12451839, 400, 0xA75BF824},
{11796481, 400, 0xFDD67368}, {11796479, 400, 0xE0237D19},
{10485761, 400, 0x15419597}, {10485759, 400, 0x154D473B},
{10223617, 400, 0x26039EB7}, {10223615, 400, 0xC9DFB1A4},
{9961473, 400, 0x3EB29644}, {9961471, 400, 0xE2AB9CB2},
{9437185, 400, 0x42609D65}, {9437183, 400, 0x77ED0792},
{8716289, 400, 0xCCA0C17B}, {8716287, 400, 0xD47E0E85},
{8257537, 400, 0x80B5C05F}, {8257535, 400, 0x278AE556},
{7798785, 400, 0x55A2468D}, {7798783, 400, 0xCF62032E},
{7471105, 400, 0x0AE03D3A}, {7471103, 400, 0xD8AB333B},
{7077889, 400, 0xC516359D}, {7077887, 400, 0xA23EA7B3},
{6684673, 400, 0xA7576F00}, {6684671, 400, 0x057E57F4},
{6422529, 400, 0xC779D2C3}, {6422527, 400, 0xA8263D37},
{6225921, 400, 0xB46AEB2F}, {6225919, 400, 0xD0A5FD5F},
{5898241, 400, 0xE46E76F9}, {5898239, 400, 0x29ED63B2},
{5505025, 400, 0x83566CC3}, {5505023, 400, 0x0B9CBE64},
{5242881, 400, 0x3CC408F6}, {5242879, 400, 0x0EA4D112},
{4980737, 400, 0x6A2056EF}, {4980735, 400, 0xE03CC669},
{4718593, 400, 0x87622D6B}, {4718591, 400, 0xF79922E2},
{4587521, 400, 0xE189A38A}, {4587519, 400, 0x930FF36C},
{4358145, 400, 0xDFEBF850}, {4358143, 400, 0xBB63D330},
{4128769, 400, 0xC0844AD1}, {4128767, 400, 0x25BDBFC3},
{3932161, 400, 0x7A525A7E}, {3932159, 400, 0xF30C9045},
{3735553, 400, 0xFAD79E97}, {3735551, 400, 0x005ED15A},
{3538945, 400, 0xDDE5BA46}, {3538943, 400, 0x15ED5982},
{3342337, 400, 0x1A6E87E9}, {3342335, 400, 0xECEEA390},
{3276801, 400, 0x3341C77F}, {3276799, 400, 0xACA2EE28},
{3112961, 400, 0x2BDF9D2B}, {3112959, 400, 0xA0AC8635},
{2949121, 400, 0x36EDB768}, {2949119, 400, 0x53FD5473},
{2785281, 400, 0x66816C94}, {2785279, 400, 0x059E8D6B},
{2654999, 400, 0x07EE900D}, {2621441, 400, 0x2BC1DACD},
{2621439, 400, 0xBCBA58F1}, {2653987, 400, 0xB005CACC},
{2651879, 400, 0x38DCD06B}, {2654003, 400, 0x1ED556E7},
{2620317, 400, 0x09DB64F8}, {2539613, 400, 0x4146EECA},
{2573917, 400, 0x939DA3B3}, {2359297, 400, 0x73A131F0},
{2359295, 400, 0x53A92203}, {2646917, 400, 0x71D4E5A2},
{2605473, 400, 0xE11637FC}, {2495213, 400, 0x89D80370},
{2540831, 400, 0x2CF01FBB}, {2654557, 400, 0x4106F46F},
{2388831, 400, 0xA508B5A7}, {2654777, 400, 0x9E744AA3},
{2584313, 400, 0x800E9A61}, {2408447, 400, 0x8C91E8AA},
{2408449, 400, 0x437ECC01}, {2345677, 400, 0x60AEE9C2},
{2332451, 400, 0xAB209667}, {2330097, 400, 0x3FB88055},
{2333851, 400, 0xFE4ECF19}, {2444819, 400, 0x56BF33C5},
{2555671, 400, 0x9DC03527}, {2654333, 400, 0xE81BCF40},
{2543123, 400, 0x379CA95D}, {2432123, 400, 0x5952676A},
{2321123, 400, 0x24DCD25F}, {2654227, 400, 0xAC3B7F2B},
{2329999, 400, 0xF5E902A5}, {2293761, 400, 0x9E4BBB8A},
{2293759, 400, 0x1901F07B}, {2236671, 400, 0x45EB162A},
{2193011, 400, 0x382B6E4B}, {2329001, 400, 0x4FF052BB},
{2327763, 400, 0x3B315213}, {2325483, 400, 0x0DC5165A},
{2323869, 400, 0xD220E27F}, {2315679, 400, 0xF650BE33},
{2004817, 400, 0xC2FF3440}, {2130357, 400, 0xC25804D8},
{2288753, 400, 0xA4DD9AAD}, {2266413, 400, 0x675257DB},
{2244765, 400, 0xC08FF487}, {2222517, 400, 0x1A128B22},
{2200339, 400, 0x0EB0E827}, {2328117, 400, 0x0A24673A},
{2329557, 400, 0x2E267692}, {2188001, 400, 0xD012AF6A},
{2166567, 400, 0x509BA41A}, {2144651, 400, 0x54CFC0E6},
{2122923, 400, 0xA47068E6}, {2100559, 400, 0xACFAB4E1},
{2088461, 400, 0xEA01E860}, {2066543, 400, 0x847DF0D0},
{2044767, 400, 0x04225888}, {2022823, 400, 0x6EA34B32},
{2328527, 400, 0xC55E3E05}, {2327441, 400, 0x207C8CEC},
{2326991, 400, 0x0A4F2ACD}, {2009987, 400, 0xE6A59DEF},
{1999999, 400, 0xD645A18F}, {1966081, 400, 0xB88828A1},
{1966079, 400, 0x5BD87C45}, {1998973, 400, 0xCBDD74F7},
{1997651, 400, 0x666B0CB1}, {1675001, 400, 0x50A94DB7},
{1977987, 400, 0x30D1CD1F}, {1955087, 400, 0x5B9426A4},
{1933071, 400, 0x23C1AF0B}, {1911957, 400, 0xF7699248},
{1899247, 400, 0x11C76E04}, {1877431, 400, 0xA3299B39},
{1855067, 400, 0x35243683}, {1833457, 400, 0xCF630DC0},
{1811987, 400, 0x7C7022EC}, {1799789, 400, 0xEFEC47B7},
{1777773, 400, 0x0F16E2D6}, {1755321, 400, 0x1AC5D492},
{1733333, 400, 0x5DA0555E}, {1711983, 400, 0xDC19DA8B},
{1699779, 400, 0x2B44914E}, {1677323, 400, 0x03D3980B},
{1995091, 400, 0x922E555B}, {1993041, 400, 0x0CA8451B},
{1991991, 400, 0xDFFB212D}, {1679779, 400, 0x51D75E0F},
{1684993, 400, 0x048BBCE8}, {1970009, 400, 0x646E0DFA},
{1957445, 400, 0xC8D244ED}, {1999997, 400, 0x5FC899D0},
{1998983, 400, 0x1CD518AA}, {1999007, 400, 0xA9DD8591},
{1674999, 400, 0xDB0169D8}, {1638401, 400, 0xD3F8A8C5},
{1638399, 400, 0xF270D8DD}, {1674997, 400, 0xC824EF15},
{1674551, 400, 0xD844AEAD}, {1674001, 400, 0x8F5EFA50},
{1345001, 400, 0x18EE2E2D}, {1655083, 400, 0x09B30DEE},
{1633941, 400, 0x0B87C8B1}, {1611557, 400, 0x6B57E48D},
{1599549, 400, 0x48EA38B2}, {1577771, 400, 0xCE84D9DC},
{1555947, 400, 0x6797EEF4}, {1533349, 400, 0xD6897409},
{1511861, 400, 0x8A8177AC}, {1499625, 400, 0x56BB6FB3},
{1477941, 400, 0xF3DD8ED3}, {1455931, 400, 0x31A222C7},
{1433069, 400, 0x28F01E1B}, {1411747, 400, 0x680C6E39},
{1399449, 400, 0xB7F01A54}, {1377247, 400, 0xE656F652},
{1355991, 400, 0xB2AA2819}, {1350061, 400, 0x31F9A728},
{1673881, 400, 0xA51D38E4}, {1672771, 400, 0x5474B6F9},
{1671221, 400, 0x2710DDEA}, {1670551, 400, 0x31FC3838},
{1660881, 400, 0x4C5B22C5}, {1650771, 400, 0x998F747B},
{1655001, 400, 0x164659A6}, {1674339, 400, 0xED2D23E2},
{1344999, 400, 0x158AA064}, {1310721, 400, 0x5694A427},
{1310719, 400, 0x258BDDE3}, {1344997, 400, 0x1D059D4F},
{1344551, 400, 0x60606AA3}, {1344001, 400, 0x9AC6AB36},
{1322851, 400, 0x3A000D0A}, {1300993, 400, 0x77CB0184},
{1288771, 400, 0x7431D9E2}, {1266711, 400, 0xB4BC4E8D},
{1244881, 400, 0x48BC9FF9}, {1222991, 400, 0x3F5FC39E},
{1200881, 400, 0xD5DF4944}, {1188441, 400, 0xD9D8968B},
{1166661, 400, 0xD4AB97F4}, {1144221, 400, 0x9940943B},
{1122001, 400, 0x647406B8}, {1100881, 400, 0x3AD40CE0},
{1088511, 400, 0xD578BB51}, {1066837, 400, 0x2F82BFBB},
{1044811, 400, 0x7C6EDDD1}, {1022991, 400, 0x6A1C2DD4},
{1000001, 400, 0x2879748F}, {1343881, 400, 0xB59E8006},
{1342771, 400, 0x87563FFE}, {1341221, 400, 0x29AD6127},
{1340551, 400, 0x17DB4ACB}, {1330881, 400, 0x9642F068},
{942079, 1000, 0xE528A9B0}, {974849, 1000, 0x79791EDB},
{983041, 1000, 0x29216C43}, {901121, 1000, 0x26C4E660},
{917503, 1000, 0x5F244685}, {933889, 1000, 0x62490F57},
{851967, 1000, 0x331AA906}, {860161, 1000, 0x41185F27},
{884735, 1000, 0x7BC7A661}, {802817, 1000, 0xA9645693},
{819199, 1000, 0x48AFB0A5}, {835585, 1000, 0x706437D3},
{753663, 1000, 0x99C43F31}, {778241, 1000, 0x1729A6C4},
{786431, 1000, 0x61080929}, {720897, 1000, 0x1E96863D},
{737279, 1000, 0x1B07A764}, {745473, 1000, 0x7BCE80AA},
{655359, 1000, 0x1107F161}, {659457, 1000, 0x589C16A4},
{688127, 1000, 0xD01E5A85}, {622593, 1000, 0x26F6FC8C},
{630783, 1000, 0x4DD2E603}, {638977, 1000, 0xC88F34B4},
{589823, 1000, 0x0290B60B}, {602113, 1000, 0xEFCD5BA8},
{614399, 1000, 0x6408F880}, {557057, 1000, 0xC30FE589},
{565247, 1000, 0xF4CA3679}, {573441, 1000, 0xF8F039AA},
{532479, 1000, 0x0072FE03}, {540673, 1000, 0xDA0E0D99},
{544767, 1000, 0x62443C6B}, {491521, 1000, 0x3F520DFA},
{516095, 1000, 0xA6BD9423}, {524289, 1000, 0xCD591388},
{466943, 1000, 0xE10EE929}, {471041, 1000, 0x18752F40},
{487423, 1000, 0x933FFF17}, {442369, 1000, 0xC22471C3},
{450559, 1000, 0x025B1320}, {458753, 1000, 0xE296CC00},
{417791, 1000, 0x080C803C}, {425985, 1000, 0xB2095F04},
{430079, 1000, 0x98B1EC61}, {393217, 1000, 0x26DD79ED},
{401407, 1000, 0x2F0F75F9}, {409601, 1000, 0xAEFAC2F8},
{372735, 1000, 0xCB6D00A2}, {376833, 1000, 0x915D5458},
{389119, 1000, 0x6188E38D}, {344065, 1000, 0x4D0C5089},
{360447, 1000, 0x84AC5CFD}, {368641, 1000, 0x72414364},
{319487, 1000, 0x24ED1BE9}, {327681, 1000, 0x3101106A},
{329727, 1000, 0x5BDB69AF}, {307201, 1000, 0x68536CD1},
{311295, 1000, 0x69778074}, {315393, 1000, 0x429D4950},
{286719, 1000, 0x1A31A686}, {294913, 1000, 0xF55727C6},
{301055, 1000, 0x33BDB242}, {272385, 1000, 0xEF6EC4B4},
{278527, 1000, 0x05530FD5}, {282625, 1000, 0x34A4E699},
{262143, 1000, 0xA9638844}, {266241, 1000, 0xE0969CED},
{270335, 1000, 0x14AD54BE}, {243713, 1000, 0xC19AEA91},
{245759, 1000, 0x7538BF0B}, {258049, 1000, 0x73F541AD},
{229375, 1000, 0x6E42B26A}, {233473, 1000, 0x1964F897},
{235519, 1000, 0x661BBC3F}, {215041, 1000, 0x04D5D2F0},
{221183, 1000, 0xA89E7764}, {225281, 1000, 0x20876BED},
{204799, 1000, 0xD20C2126}, {208897, 1000, 0x9D4DCF0E},
{212991, 1000, 0x1FF00E2A}, {194561, 1000, 0x6ED1CB70},
{196607, 1000, 0x3190D5F5}, {200705, 1000, 0xFAD28F5A},
{184319, 1000, 0x360EF08E}, {186369, 1000, 0x0F001482},
{188415, 1000, 0x86FCE4D6}, {164865, 1000, 0x4942B002},
{172031, 1000, 0xC5AF29DB}, {180225, 1000, 0x35D49D74},
{157695, 1000, 0x5422FACF}, {159745, 1000, 0xB5CD03A1},
{163839, 1000, 0x1CA6048E}, {150529, 1000, 0x7412F09C},
{153599, 1000, 0xA9FAAE69}, {155649, 1000, 0xA7B736AF},
{141311, 1000, 0x7A5D0730}, {143361, 1000, 0x580F4DC4},
{147455, 1000, 0x176B299A}, {135169, 1000, 0x65AC10A4},
{136191, 1000, 0xC4591D37}, {139265, 1000, 0xBCE1FC80},
{129023, 1000, 0xAFE1E7A8}, {131073, 1000, 0xC5AAB12F},
{133119, 1000, 0xDE51C35A}, {117761, 1000, 0x054A26F6},
{121855, 1000, 0x55AF2385}, {122881, 1000, 0x652827AC},
{112639, 1000, 0x6FA4DB24}, {114689, 1000, 0x0BBAF161},
{116735, 1000, 0xB85F0E8E}, {106497, 1000, 0xF833D925},
{107519, 1000, 0x80F177D8}, {110593, 1000, 0x1A56AA86},
{100351, 1000, 0x1DE12CE6}, {102401, 1000, 0x19F967B4},
{104447, 1000, 0xF9F3CDFD}
};

#define MAX_SELF_TEST_ITERS2	376
struct self_test_info SELF_TEST_DATA2[MAX_SELF_TEST_ITERS2] = {
{560000001, 100, 0x7F853A0A}, {420000001, 150, 0x89665E7E},
{280000001, 200, 0xC32CAD46}, {210000001, 300, 0x89823329},
{140000001, 400, 0x15EF4F24}, {110000001, 500, 0x893C9000},
{77497473, 900, 0xF0B43F54}, {76497471, 900, 0xF30AFA95},
{75497473, 900, 0x32D8D3A7}, {75497471, 900, 0x9E689331},
{74497473, 900, 0xD43166A4}, {73497471, 900, 0x639E4F0C},
{72303169, 900, 0x74BDED5C}, {71303169, 900, 0xA2147B5C},
{71303167, 900, 0x717525AB}, {70303167, 900, 0xD716B4F0},
{68060289, 1000, 0xF90C7BFF}, {67060287, 1000, 0xFE9BF47C},
{66060289, 1000, 0x057C60F5}, {66060287, 1000, 0x2ECC97CE},
{65390273, 1000, 0xC55C6369}, {64390271, 1000, 0x48552448},
{63390273, 1000, 0x6FF8CD84}, {62390273, 1000, 0x42ACEB15},
{62390271, 1000, 0x48764DF8}, {61390271, 1000, 0xD5408698},
{57623105, 1200, 0x098B4491}, {56623105, 1200, 0x5E720717},
{56623103, 1200, 0x1980D8BC}, {55623103, 1200, 0xEDD592B6},
{53477377, 1200, 0xBAEF5CCC}, {53477375, 1200, 0x2F296FC8},
{52331647, 1200, 0xA1EAE85D}, {51331649, 1200, 0xE3B39845},
{50331649, 1200, 0x53543DF2}, {50331647, 1200, 0x0049E54B},
{48185921, 1500, 0x78F4AEAA}, {47185921, 1500, 0x4D7FFDDC},
{47185919, 1500, 0x059D196F}, {46185919, 1500, 0x38B1D9AD},
{45943041, 1500, 0x7670FDDF}, {44943039, 1500, 0xA859BBD7},
{43943041, 1500, 0xD673E000}, {42943039, 1500, 0x6B69D8CE},
{41943041, 1500, 0x6E92CE47}, {41943039, 1500, 0x888BEE79},
{39151585, 1900, 0x3B06496C}, {38748737, 1900, 0x6429E0FD},
{38251583, 1900, 0x04AD7F99}, {37748737, 1900, 0x47659BC5},
{37748735, 1900, 0x2DFA41B0}, {36748735, 1900, 0x1A1DA557},
{36251585, 1900, 0x83F23FA8}, {35651585, 1900, 0x3598B4B9},
{35651583, 1900, 0x7E443962}, {35251583, 1900, 0x1CE4D084},
{34230145, 2100, 0x0FDE9717}, {33730143, 2100, 0x54EB5333},
{33030145, 2100, 0xF37897B8}, {33030143, 2100, 0x52B3981B},
{32595137, 2100, 0xA76D0805}, {32095135, 2100, 0xCF443ACD},
{31595137, 2100, 0xA6DEA70A}, {31195137, 2100, 0x0777442D},
{31195135, 2100, 0x9B265F8F}, {30695135, 2100, 0xA3BC760F},
{29311553, 2500, 0xFD1D6D74}, {28811551, 2500, 0xE720BFD3},
{28311553, 2500, 0xA11F75AB}, {28311551, 2500, 0x7E0471E5},
{27738689, 2500, 0xD246DC55}, {27238687, 2500, 0x806A3A62},
{26738689, 2500, 0x8E8450B1}, {26738687, 2500, 0xD4A0DBC9},
{26138689, 2500, 0x47C47755}, {25638687, 2500, 0x7E9C7E8E},
{24903681, 3100, 0x50835AB8}, {24903679, 3100, 0xAE3D2F94},
{24092961, 3100, 0x7B540B4D}, {23892959, 3100, 0xA0D4EC50},
{23592961, 3100, 0x47FBD6FE}, {23592959, 3100, 0x09FD89AB},
{22971521, 3100, 0x99DFEDB9}, {21871519, 3100, 0x35A8B46A},
{20971521, 3100, 0x94C12572}, {20971519, 3100, 0x1F6D3003},
{19922945, 4000, 0x86B106EB}, {19922943, 4000, 0xE1CE3C1A},
{19374367, 4000, 0xD1045A66}, {19174369, 4000, 0x3247CE82},
{18874369, 4000, 0x33BB2689}, {18874367, 4000, 0x6856F21F},
{18474367, 4000, 0x95E2F6FA}, {18274367, 4000, 0x61182009},
{18274369, 4000, 0xB2FD8175}, {18074369, 4000, 0x7F242A6E},
{17432577, 4500, 0x632CAD0B}, {17432575, 4500, 0xC9C79F07},
{17115073, 4500, 0xF2B70D4B}, {16815071, 4500, 0x71B22529},
{16515073, 4500, 0xAB1CC854}, {16515071, 4500, 0xF54D05D7},
{16297569, 4500, 0x6B5F72DA}, {15997567, 4500, 0x9669F188},
{15597569, 4500, 0x352BFCCF}, {15597567, 4500, 0x36B164ED},
{14942209, 5300, 0xEA5DB53B}, {14942207, 5300, 0x6CC650A2},
{14155777, 5300, 0xEB7C125D}, {14155775, 5300, 0xB4C8B09B},
{13969343, 5300, 0x832359A5}, {13669345, 5300, 0x7EE99140},
{13369345, 5300, 0xCDF43471}, {13369343, 5300, 0x343FEA12},
{13069345, 5300, 0x65B17A9B}, {12969343, 5300, 0x063F492B},
{12451841, 6500, 0xCB168E5D}, {12451839, 6500, 0xE91EEB5A},
{12196481, 6500, 0x0A261B7E}, {11796481, 6500, 0x38100A5F},
{11796479, 6500, 0x78FCF8C5}, {11596479, 6500, 0x8C481635},
{11285761, 6500, 0x2580BC8D}, {10885759, 6500, 0x54030992},
{10485761, 6500, 0x054660AA}, {10485759, 6500, 0x50F74AF0},
{9961473, 7800, 0x7991161C}, {9961471, 7800, 0x627F3BEE},
{9837183, 7800, 0xBC67A608}, {9737185, 7800, 0x9A0CBC59},
{9537183, 7800, 0xA6A509A6}, {9437185, 7800, 0x877C09B6},
{9437183, 7800, 0x1D259540}, {9337185, 7800, 0x5EF3F14C},
{9237183, 7800, 0x5780245F}, {9137185, 7800, 0x6C1162A9},
{8716289, 9000, 0x2011133F}, {8716287, 9000, 0xEEEC1181},
{8516289, 9000, 0xF1D93A69}, {8316287, 9000, 0x53D6E3CB},
{8257537, 9000, 0x38DB98D6}, {8257535, 9000, 0x7D1BECA7},
{8098785, 9000, 0x51E9FA27}, {7998783, 9000, 0xF7F14FF2},
{7798785, 9000, 0x8437BC4D}, {7798783, 9000, 0x9E28D8E1},
{7471105, 11000, 0xEFDA89EA}, {7471103, 11000, 0x4061C4BF},
{7377889, 11000, 0x65ABE846}, {7277887, 11000, 0x02B0EBD7},
{7077889, 11000, 0x336E1030}, {7077887, 11000, 0x685B792E},
{6984673, 11000, 0x3AE19FAF}, {6884671, 11000, 0x2A0ED16A},
{6684673, 11000, 0x206A3512}, {6684671, 11000, 0x4FD9980A},
{6225921, 13000, 0x1A922371}, {6225919, 13000, 0xC0F63BD8},
{6198241, 13000, 0xDA664501}, {6098239, 13000, 0xB92015CD},
{5898241, 13000, 0xDA384BD9}, {5898239, 13000, 0x20B59AC8},
{5705025, 13000, 0x941A2DA0}, {5605023, 13000, 0xCFDF5835},
{5505025, 13000, 0x37A6C972}, {5505023, 13000, 0x6252AB5C},
{5120737, 17000, 0x512705D0}, {5030735, 17000, 0x633E3E74},
{4980737, 17000, 0xD8245D49}, {4980735, 17000, 0xFB2C3530},
{4888593, 17000, 0xE3C6EDBC}, {4818591, 17000, 0x89E7FE48},
{4718593, 17000, 0xA23C713D}, {4718591, 17000, 0xC7BA41D6},
{4698593, 17000, 0xA0194103}, {4648591, 17000, 0xD5A50A23},
{4501145, 19000, 0x7BAF4344}, {4458143, 19000, 0x686F6B13},
{4358145, 19000, 0x682E6643}, {4358143, 19000, 0x974DA6CC},
{4298769, 19000, 0x1FC0E577}, {4228767, 19000, 0x46B5F3CD},
{4128769, 19000, 0x59332478}, {4128767, 19000, 0x4AF5C8B8},
{4028769, 19000, 0x542C17CB}, {3978767, 19000, 0x76E41351},
{3835553, 22000, 0x9058FE40}, {3785551, 22000, 0x45EF5C15},
{3735553, 22000, 0x2700B350}, {3735551, 22000, 0x09EDCEAD},
{3688945, 22000, 0x626C29D3}, {3618943, 22000, 0x82B1D4D1},
{3538945, 22000, 0x70331CC6}, {3538943, 22000, 0x00FEB746},
{3342337, 22000, 0x7CEE24AE}, {3342335, 22000, 0x1802D072},
{3242961, 27000, 0xE877F863}, {3172959, 27000, 0x04C9F1F7},
{3112961, 27000, 0x241E93DB}, {3112959, 27000, 0x8D359307},
{2949121, 27000, 0x6B545E09}, {2949119, 27000, 0xAFD6F417},
{2885281, 27000, 0x439E57E6}, {2785281, 27000, 0xB4E40DFE},
{2785279, 27000, 0x3787D3FA}, {2685279, 27000, 0x902967B7},
{2605473, 34000, 0xE21C344E}, {2584313, 34000, 0xFDBCFCB2},
{2573917, 34000, 0x89B5012C}, {2540831, 34000, 0x201BAA90},
{2539613, 34000, 0x2226BA6B}, {2495213, 34000, 0xE3577D9F},
{2408447, 34000, 0x594C9155}, {2388831, 34000, 0x55CE9F16},
{2359297, 34000, 0x09A72A40}, {2359295, 34000, 0x621E8BF9},
{2244765, 39000, 0xEC2F362D}, {2236671, 39000, 0x4B50CA20},
{2222517, 39000, 0x8DA427C0}, {2193011, 39000, 0xD1DE8993},
{2130357, 39000, 0x4B5EBB90}, {2122923, 39000, 0x5F9110FC},
{2100559, 39000, 0xE0CF8904}, {2088461, 39000, 0x26AD1DEA},
{2066543, 39000, 0xB78C9237}, {2004817, 39000, 0x3D7838F8},
{1933071, 46000, 0x86323D21}, {1911957, 46000, 0x500CFEAD},
{1899247, 46000, 0x128667DF}, {1877431, 46000, 0x2A59B6B5},
{1855067, 46000, 0xBE9AABF5}, {1833457, 46000, 0xB84D7929},
{1777773, 46000, 0x771E0A9D}, {1755321, 46000, 0xF93334E3},
{1699779, 46000, 0x07B46DEE}, {1677323, 46000, 0x910E0320},
{1633941, 56000, 0x455509CD}, {1611557, 56000, 0x0F51FA1E},
{1599549, 56000, 0x646A96B0}, {1577771, 56000, 0xA4A21303},
{1555947, 56000, 0x80B84725}, {1533349, 56000, 0x23E9F7B1},
{1477941, 56000, 0x593F208F}, {1455931, 56000, 0x11002C52},
{1433069, 56000, 0x5B641D8B}, {1411747, 56000, 0x5EAE18A8},
{1322851, 75000, 0xD5C50F2E}, {1310721, 75000, 0x855E44A2},
{1310719, 75000, 0xC0836C1F}, {1300993, 75000, 0xF62263D6},
{1288771, 75000, 0x867EBBAB}, {1266711, 75000, 0xBA1FF3BE},
{1244881, 75000, 0xCE8199EB}, {1222991, 75000, 0xCDE49EF5},
{1200881, 75000, 0xC8610F6C}, {1188441, 75000, 0xFC772495},
{1150221, 84000, 0xA3334541}, {1144221, 84000, 0x44307B03},
{1122001, 84000, 0x9B937DCF}, {1108511, 84000, 0x9F3D191E},
{1100881, 84000, 0xBAF4EA2D}, {1096837, 84000, 0xAA9396F1},
{1088511, 84000, 0xB0CB2704}, {1066837, 84000, 0x031F202C},
{1044811, 84000, 0x7EA89CFE}, {1022991, 84000, 0xD42294C8},
{983041, 100000, 0x4052BBC0}, {974849, 100000, 0xB0E9EB07},
{942079, 100000, 0xEE230987}, {933889, 100000, 0x58FA63B0},
{917503, 100000, 0x8B457209}, {901121, 100000, 0xD2325FC4},
{884735, 100000, 0xCBB5A603}, {860161, 100000, 0xBC240C77},
{854735, 100000, 0xE8BE766D}, {851967, 100000, 0x09AD9B74},
{827279, 120000, 0x64B01894}, {819199, 120000, 0xF97F1E2B},
{802817, 120000, 0xC4EDBC3C}, {795473, 120000, 0x046584E0},
{786431, 120000, 0xC6BA553D}, {778241, 120000, 0x856A5147},
{753663, 120000, 0xC7895B4A}, {745473, 120000, 0x42B47EA2},
{737279, 120000, 0x29E477B8}, {720897, 120000, 0x97111FA7},
{662593, 160000, 0x32472A99}, {659457, 160000, 0xEF49D340},
{655359, 160000, 0x75C12C38}, {644399, 160000, 0xDE632783},
{638977, 160000, 0xDCDB98B4}, {630783, 160000, 0x6B8F0706},
{622593, 160000, 0xD732286D}, {614399, 160000, 0x2489EFB3},
{612113, 160000, 0xCAE00EC6}, {602113, 160000, 0x792AD67D},
{580673, 180000, 0xC508CAFA}, {573441, 180000, 0xB0680C2B},
{565247, 180000, 0xF1DBB762}, {557057, 180000, 0x374F647B},
{544767, 180000, 0x3DC41F49}, {540673, 180000, 0x949A4CB7},
{532479, 180000, 0xEA06DC97}, {524289, 180000, 0xA76CE14A},
{522479, 180000, 0xAA8EAC14}, {516095, 180000, 0x04F0CC23},
{501041, 210000, 0xD9F72F62}, {496943, 210000, 0xD62D5380},
{487423, 210000, 0x55ACB2FD}, {471041, 210000, 0xB6AEAB0E},
{466943, 210000, 0x251CDE78}, {458753, 210000, 0xDC40CADB},
{450559, 210000, 0x2AD0CF72}, {442369, 210000, 0x5FF2E46E},
{441041, 210000, 0x1194CC23}, {436943, 210000, 0x0272AF35},
{420217, 270000, 0xD233852A}, {409601, 270000, 0x6F89825C},
{401407, 270000, 0x3D9DE818}, {393217, 270000, 0xDE8E6FF0},
{392119, 270000, 0x30CA58B7}, {389119, 270000, 0x80975797},
{376833, 270000, 0xC75824DB}, {372735, 270000, 0xF8BE0932},
{368641, 270000, 0xA48AC5E3}, {360447, 270000, 0x7DD29C13},
{339487, 340000, 0xA7311A6D}, {335393, 340000, 0xD9704DF2},
{331681, 340000, 0x3316A003}, {329727, 340000, 0xE46D5991},
{327681, 340000, 0xBEDA4A7B}, {319487, 340000, 0xB25C84FF},
{315393, 340000, 0xF5AD1DDA}, {311295, 340000, 0xFE41A12A},
{308295, 340000, 0x03AAC47E}, {307201, 340000, 0xFC08ACCC},
{291913, 380000, 0xC56AB884}, {286719, 380000, 0x248EF622},
{282625, 380000, 0x50A98488}, {280335, 380000, 0x9B64A843},
{278527, 380000, 0x39D5B7DB}, {274335, 380000, 0x48623B41},
{270335, 380000, 0xC04B857A}, {266241, 380000, 0xFE4475F6},
{262143, 380000, 0xADC3ECE9}, {260335, 380000, 0x15B8F9EF},
{250519, 460000, 0xA2FE3B50}, {245759, 460000, 0xC6D800D6},
{245281, 460000, 0x4F23AA34}, {243713, 460000, 0xB30EC823},
{235519, 460000, 0x31FD709E}, {233473, 460000, 0x8FCC69C2},
{231183, 460000, 0xD59255CC}, {229375, 460000, 0x788520D0},
{225281, 460000, 0xD669C8BC}, {221183, 460000, 0x9B915F4B},
{212991, 560000, 0x0555250D}, {210415, 560000, 0x3FC3CCD7},
{208897, 560000, 0x9FF8F462}, {204799, 560000, 0x294EB549},
{200705, 560000, 0x80B1222F}, {196607, 560000, 0x8AB8D945},
{194561, 560000, 0x4140E623}, {188415, 560000, 0xFA0A3453},
{186369, 560000, 0xAC17EAB6}, {184319, 560000, 0x835F341B},
{172031, 800000, 0xF6BD0728}, {163839, 800000, 0x26C78657},
{159745, 800000, 0x6ACBB961}, {157695, 800000, 0x3EA979F3},
{155649, 800000, 0x09C7ADE4}, {153599, 800000, 0xF601EB92},
{147455, 800000, 0x0AA97D21}, {143361, 800000, 0xEA6A01F1},
{141311, 800000, 0x9BB8A6A3}, {135169, 800000, 0xECA55A45}
};

#ifdef ONE_HOUR_SELF_TEST
int selfTest (
	int	thread_num,
	struct PriorityInfo *sp_info,
	struct work_unit *w)
{
	unsigned long fftlen;
	char	iniName[32];
	int	self_test_errors, self_test_warnings;

/* What fft length are we running? */

	if (w->forced_fftlen)
		fftlen = w->forced_fftlen;
	else
		fftlen = gwmap_to_fftlen (1.0, 2, w->n, -1);

/* If fftlength is less than 64K return (we don't have any small exponents */
/* in our self test data) */

	if (fftlen < 65536) return (0);

/* Make sure we haven't run this self-test already. */

	sprintf (iniName, SelfTestIniMask, (int) (fftlen/1024));
	if (IniGetInt (LOCALINI_FILE, iniName, 0)) return (0);
#ifdef SERVER_TESTING
	return (0);
#endif

/* Make sure the user really wants to spend an hour doing this now */

	OutputStr (thread_num, SELFMSG1A);
	OutputStr (thread_num, SELFMSG1B);
	OutputStr (thread_num, SELFMSG1C);

/* Do the self test */

	self_test_errors = 0;
	self_test_warnings = 0;
	return (selfTestInternal (thread_num, sp_info, fftlen, 60, NULL, 0, NULL,
				  &self_test_errors, &self_test_warnings));
}
#endif

int tortureTest (
	int	thread_num,
	int	num_threads)
{
	struct PriorityInfo sp_info;
#define SELF_FFT_LENGTHS	49
	int	lengths[SELF_FFT_LENGTHS] = {1024,8,10,896,768,12,14,640,512,16,20,448,384,24,28,320,256,32,40,224,192,48,56,160,128,64,80,112,96,1280,1536,1792,2048,2560,3072,3584,4096,5120,6144,7168,8192,10240,12288,14336,16384,20480,24576,28672,32768};
	int	data_index[SELF_FFT_LENGTHS] = {0};
	int	min_fft, max_fft, test_time;
	int	self_test_errors, self_test_warnings;
	int	i, run_indefinitely, stop_reason;
	time_t	start_time, current_time;
	unsigned int memory;	/* Memory to use during torture test */
	void	*bigbuf = NULL;

/* Set the process/thread priority */

	sp_info.type = SET_PRIORITY_TORTURE;
	sp_info.thread_num = thread_num;
	sp_info.aux_thread_num = 0;
	sp_info.num_threads = num_threads;
	SetPriority (&sp_info);

/* We used to support a menu option to run the self-test for an hour on */
/* each FFT length.  If we ever decide to resupport this option, change */
/* the run_indefiitely argument to an argument and change the output */
/* message below. */

loop:	run_indefinitely = TRUE;

/* Make sure the user really wants to spend many hours doing this now */

	if (run_indefinitely) {
		OutputStr (thread_num, TORTURE1);
		OutputStr (thread_num, TORTURE2);
		test_time = IniGetInt (INI_FILE, "TortureTime", 15);
	}

/* Determine fft lengths we should run and allocate a big block */
/* of memory to test. */

	min_fft = IniGetInt (INI_FILE, "MinTortureFFT", 8);
	max_fft = IniGetInt (INI_FILE, "MaxTortureFFT", 4096);
	memory = IniGetInt (INI_FILE, "TortureMem", 8);
	while (memory > 8 && bigbuf == NULL) {
		bigbuf = aligned_malloc (memory * 1000000, 128);
		if (bigbuf == NULL) memory--;
	}

/* Now self-test each fft length */

	stop_reason = 0;
	self_test_errors = 0;
	self_test_warnings = 0;
	time (&start_time);
	for ( ; ; ) {
	    for (i = 0; i < SELF_FFT_LENGTHS; i++) {
		if (lengths[i] < min_fft || lengths[i] > max_fft)
			continue;
		stop_reason =
			selfTestInternal (thread_num, &sp_info, lengths[i]*1024,
					  test_time, &data_index[i],
					  memory, bigbuf,
					  &self_test_errors,
					  &self_test_warnings);
		if (stop_reason) {
			char	buf[120];
			int	hours, minutes;
			time (&current_time);
			minutes = (int) (current_time - start_time) / 60;
			hours = minutes / 60;
			minutes = minutes % 60;
			strcpy (buf, "Torture Test ran ");
			if (hours)
				sprintf (buf+strlen(buf), "%d hours, ", hours);
			sprintf (buf+strlen(buf),
				 "%d minutes - %d errors, %d warnings.\n",
				 minutes, self_test_errors, self_test_warnings);
			OutputStr (thread_num, buf);
			run_indefinitely = FALSE;
			break;
		}
	    }
	    if (! run_indefinitely) break;
	}

/* Self test completed!  Free memory. */

	aligned_free (bigbuf);

/* If this was a user requested stop, then wait for a restart */
	
	while (stop_reason == STOP_WORKER) {
		implement_stop_one_worker (thread_num);
		stop_reason = stopCheck (thread_num);
		if (stop_reason == 0) goto loop;
	}

/* All done */

	return (stop_reason);
}

int selfTestInternal (
	int	thread_num,
	struct PriorityInfo *sp_info,
	unsigned long fftlen,
	unsigned int test_time,	/* Number of minutes to self-test */
	int	*torture_index,	/* Index into self test data array */
	unsigned int memory,	/* MB of memory the torture test can use */
	void	*bigbuf,	/* Memory block for the torture test */
	int	*errors,	/* Returned count of self test errors */
	int	*warnings)	/* Returned count of self test warnings */
{
	llhandle lldata;
	unsigned long k, limit, num_threads;
	unsigned int i, iter;
	char	buf[120];
	char	iniName[32];
	time_t	start_time, current_time;
	struct self_test_info *test_data;
	unsigned int test_data_count;
	int	cpu_supports_3dnow, stop_reason;

/* Set the title */

	title (thread_num, "Self-Test");

/* Pick which self test data array to use.  Machines are much faster now */
/* compared to when the torture test was introduced.  This new self test */
/* data will run more iterations and thus stress the cpu more by spending */
/* less time in the initialization code. */

	if (CPU_SPEED < 1000.0) {
		test_data = SELF_TEST_DATA;
		test_data_count = MAX_SELF_TEST_ITERS;
	} else {
		test_data = SELF_TEST_DATA2;
		test_data_count = MAX_SELF_TEST_ITERS2;
	}

/* Decide how many threads the torture test can use.  This should only */
/* really be needed for QA purposes as the user can probably create more */
/* more stress by running one torture test window for each CPU core. */

	num_threads = IniGetInt (INI_FILE, "TortureTestThreads", 1);

/* Determine the range from which we'll choose an exponent to test. */

	limit = gwmap_fftlen_to_max_exponent (fftlen);

/* Get the current time */

	time (&start_time);

/* Start in the self test data array where we left off the last time */
/* torture test executed this FFT length. */

	i = (torture_index == NULL) ? 0 : *torture_index;

/* Loop testing various exponents from self test data array until */
/* time runs out */

	for (iter = 1; ; iter++) {
		char	fft_desc[80];
		unsigned long p, reshi, reslo;
		unsigned int ll_iters, num_gwnums;
		gwnum	*gwarray, g;

/* Find next self test data entry to work on */

		for ( ; ; i++) {

/* Wrap in the self test data array */

			if (i == test_data_count) i = 0;

/* Now select the actual exponent */

			p = test_data[i].p;
			if (p > limit) continue;

/* The SSE2 carry propagation code gets into trouble if there are too */
/* few bits per FFT word!  Thus, we'll require at least 8 bits per */
/* word here.  Now that the number of iterations changes for each FFT */
/* length I'm raising the requirement to 10 bits to keep timings roughly */
/* equal. */

			if (p / fftlen < 10) continue;

/* We've found an exponent to test! */

			break;
		}

/* Anecdotal evidence suggests that some AMD chips struggle with FFTs */
/* that do not use the prefetchw instruction (they would fail a version */
/* 23.8 torture test, but pass a version 24.11 torture test).  Thus, I've */
/* made the torture test use both the prefetchw and non-prefetchw FFTs. */

		cpu_supports_3dnow = (CPU_FLAGS & CPU_3DNOW) ? TRUE : FALSE;

/* Now run Lucas setup, for extra safety double the maximum allowable */
/* sum(inputs) vs. sum(outputs) difference. */

		gwinit (&lldata.gwdata);
		gwset_num_threads (&lldata.gwdata, num_threads);
		lldata.gwdata.GW_BIGBUF = (char *) bigbuf;
		lldata.gwdata.GW_BIGBUF_SIZE = (bigbuf != NULL) ? memory * 1000000 : 0;
		if (cpu_supports_3dnow && p > 5000000 &&
		    torture_index != NULL &&
		    (test_data[i].reshi & 1) &&
		    IniGetInt (INI_FILE, "Special3DNowTorture", 1))
			CPU_FLAGS &= ~CPU_3DNOW;
		stop_reason = lucasSetup (thread_num, p, fftlen, &lldata);
		if (cpu_supports_3dnow)
			CPU_FLAGS |= CPU_3DNOW;
		if (stop_reason) return (stop_reason);
		lldata.gwdata.MAXDIFF *= 2.0;

/* Output start message */

		ll_iters = test_data[i].iters;
		gwfft_description (&lldata.gwdata, fft_desc);
		sprintf (buf, SELF1, iter, ll_iters, p, fft_desc);
		OutputStr (thread_num, buf);

/* Determine how many gwnums we can allocate in the memory we are given */

		if (memory <= 8 || (iter & 1) == 0)
			num_gwnums = 1;
		else {
			num_gwnums = (unsigned int)
				(((double) memory * 1000000.0 -
				  (double) gwmemused (&lldata.gwdata)) /
				 (double) (gwnum_size (&lldata.gwdata)));
			if (num_gwnums < 1) num_gwnums = 1;
			if (num_gwnums > ll_iters) num_gwnums = ll_iters;
		}

/* Allocate gwnums to eat up the available memory */

		gwarray = (gwnum *) malloc (num_gwnums * sizeof (gwnum));
		if (gwarray == NULL) {
			lucasDone (&lldata);
			return (OutOfMemory (thread_num));
		}
		gwarray[0] = lldata.lldata;
		for (k = 1; k < num_gwnums; k++) {
			gwarray[k] = gwalloc (&lldata.gwdata);
			if (gwarray[k] == NULL) {
				num_gwnums = k;
				break;
			}
		}

/* Init data area with a pre-determined value */

restart_test:	dbltogw (&lldata.gwdata, 4.0, lldata.lldata);
		g = lldata.lldata;

/* Do Lucas-Lehmer iterations */

		for (k = 0; k < ll_iters; k++) {

/* Copy previous squared value (so we plow through memory) */

			if (k && num_gwnums > 1) {
				gwnum	prev;
				prev = g;
				g = gwarray[k % num_gwnums];
				gwcopy (&lldata.gwdata, prev, g);
			}

/* One Lucas-Lehmer test with error checking */

			gwsetnormroutine (&lldata.gwdata, 0, 1, 0);
			gwstartnextfft (&lldata.gwdata, k != ll_iters - 1);
			lucas_fixup (&lldata, p);
			gwsquare (&lldata.gwdata, g);

/* If the sum of the output values is an error (such as infinity) */
/* then raise an error. */

			if (gw_test_illegal_sumout (&lldata.gwdata)) {
				OutputBoth (thread_num, SELFFAIL1);
				(*warnings)++;
				if (*warnings < 100) {
					OutputBoth (thread_num, SELFFAIL4);
					goto restart_test;
				} else {
					OutputBoth (thread_num, SELFFAIL6);
					lucasDone (&lldata);
					free (gwarray);
					return (STOP_FATAL_ERROR);
				}
			}

/* Check that the sum of the input numbers squared is approximately */
/* equal to the sum of unfft results. */

			if (gw_test_mismatched_sums (&lldata.gwdata)) {
				sprintf (buf, SELFFAIL2,
					 gwsumout (&lldata.gwdata, g),
					 gwsuminp (&lldata.gwdata, g));
				OutputBoth (thread_num, buf);
				OutputBoth (thread_num, SELFFAIL5);
				(*errors)++;
				lucasDone (&lldata);
				free (gwarray);
				return (STOP_FATAL_ERROR);
			}

/* Make sure round off error is tolerable */

			if (gw_get_maxerr (&lldata.gwdata) > 0.45) {
				sprintf (buf, SELFFAIL3, gw_get_maxerr (&lldata.gwdata));
				OutputBoth (thread_num, buf);
				OutputBoth (thread_num, SELFFAIL5);
				(*errors)++;
				lucasDone (&lldata);
				free (gwarray);
				return (STOP_FATAL_ERROR);
			}

/* Abort if user demands it */

			stop_reason = stopCheck (thread_num);
			if (stop_reason) {
				lucasDone (&lldata);
				free (gwarray);
				return (stop_reason);
			}
		}

/* Compare final 32 bits with the pre-computed array of correct residues */

		if (g != lldata.lldata)
			gwcopy (&lldata.gwdata, g, lldata.lldata);
		generateResidue64 (&lldata, &reshi, &reslo);
		lucasDone (&lldata);
		free (gwarray);
		if (reshi != test_data[i].reshi) {
			sprintf (buf, SELFFAIL, reshi, test_data[i].reshi);
			OutputBoth (thread_num, buf);
			OutputBoth (thread_num, SELFFAIL5);
			(*errors)++;
			return (STOP_FATAL_ERROR);
		}

/* Bump index into self test data array */

		i++;

/* Has time expired? */

		time (&current_time);
		if ((unsigned int) (current_time - start_time) >= test_time * 60) break;
	}

/* Save our position in self test data array for next time torture test */
/* executes this FFT length */

	if (torture_index != NULL) *torture_index = i;

/* We've passed the self-test.  Remember this in the .INI file */
/* so that we do not need to do this again. */

	sprintf (buf, SELFPASS, (int) (fftlen/1024));
	OutputBoth (thread_num, buf);
	sprintf (iniName, SelfTestIniMask, (int) (fftlen/1024));
	IniWriteInt (LOCALINI_FILE, iniName, 1);
	return (0);
}

/*******************************************/
/* Various QA and data analysis functions! */
/*******************************************/

/* Read a file of exponents to run LL iterations on as part of a QA process */
/* The format of this file is: */
/*	exponent,optional fft length,num iters,optional shift count,residue */
/* An Advanced/Time 9999 corresponds to type 0, Advanced/Time 9998 */
/* corresponds to type 1, etc. */
/* Type 0 executes much like an LL test, error checking and doing a */
/* careful iteration occasionally */
/* Type 1 does roundoff checking every iteration and accumulates */
/* statistics on the round off data. */
/* Type 2 and higher have not been used much and may not work */

int lucas_QA (
	int	thread_num,
	int	type)
{
	llhandle lldata;
	FILE	*fd;
	int	stop_reason;

/* Set the title, init random generator */

	title (thread_num, "QA");
	srand ((unsigned) time (NULL));

/* Open QA file */

	fd = fopen ("qa", "r");
	if (fd == NULL) {
		OutputStr (thread_num, "File named 'qa' could not be opened.\n");
		return (STOP_FILE_IO_ERROR);
	}

/* Loop until the entire file is processed */

	for ( ; ; ) {
		unsigned long p, fftlen, iters;
		char	buf[500], res[80];
		unsigned long reshi, reslo, units_bit;
		unsigned long i, word, bit_in_word, maxerrcnt, loops;
		double	maxsumdiff, maxerr, toterr, M, S;
		unsigned long ge_300, ge_325, ge_350, ge_375, ge_400;
		gwnum	t1, t2;
		unsigned int iters_unchecked;

/* Read a line from the file */

		p = 0;
		fscanf (fd, "%lu,%lu,%lu,%lu,%s\n",
			&p, &fftlen, &iters, &units_bit, &res);
		if (p == 0) break;

/* In a type 4 run, we decrement through exponents to find any with */
/* anamolously high average errors.  After selecting a tentative FFT */
/* crossover, we do a type 4 run looking for a higher average error */
/* below the crossover we selected. */

		for (loops = (type != 4 ? 1 : units_bit % 100); loops--; p-=2){

/* Now run Lucas setup */

		gwinit (&lldata.gwdata);
		stop_reason = lucasSetup (thread_num, p, fftlen, &lldata);
		lldata.units_bit = units_bit;
		if (stop_reason) return (stop_reason);
		maxsumdiff = 0.0;
		ge_300 = ge_325 = ge_350 = ge_375 = ge_400 = 0;
		maxerr = 0.0; maxerrcnt = 0; toterr = 0.0;
		iters_unchecked = (type > 3) ? 2 : 40;

/* Check for a randomized units bit */

		if (lldata.units_bit >= p) {
			uint32_t hi, lo;
			lldata.units_bit = (rand () << 16) + rand ();
			if (CPU_FLAGS & CPU_RDTSC) {
				rdtsc (&hi, &lo);
				lldata.units_bit += lo;
			}
			lldata.units_bit = lldata.units_bit % p;
			sprintf (buf, "Units bit = %lu\n", lldata.units_bit);
			OutputBoth (thread_num, buf);
		}

/* Init data area with a pre-determined value */

		bitaddr (&lldata.gwdata, (lldata.units_bit + 2) % p, &word, &bit_in_word);
		for (i = 0; i < gwfftlen (&lldata.gwdata); i++) {
			set_fft_value (&lldata.gwdata, lldata.lldata, i,
				       (type == 3 || type == 4) ?
					 (rand () & 1) ? rand () : -rand () :
				       (i == word) ? (1L << bit_in_word) : 0);
		}

/* The thorough, P-1, and ECM tests use more than one number */

		if (type == 2 || type == 3) {
			t1 = gwalloc (&lldata.gwdata);
			dbltogw (&lldata.gwdata, 234872639921.0, t1);
			gwfft (&lldata.gwdata, t1, t1);
			t2 = gwalloc (&lldata.gwdata);
			dbltogw (&lldata.gwdata, 1982387192367.0, t2);
			gwfft (&lldata.gwdata, t2, t2);
			lldata.gwdata.MAXDIFF *= 16;
		}

/* Do Lucas-Lehmer iterations */

		for (i = 0; i < iters; i++) {

/* One Lucas-Lehmer iteration with error checking */

			if (type == 0) {		/* Typical LL test */
				gwsetnormroutine (&lldata.gwdata, 0, (i & 63) == 37, 0);
				gwstartnextfft (&lldata.gwdata, i < iters / 2);
				if (i > iters / 2 && (i & 63) == 44)
					careful_iteration (&lldata, p);
				else {
					lucas_fixup (&lldata, p);
					gwsquare (&lldata.gwdata, lldata.lldata);
				}
			} else if (type == 1 || type == 4) { /* Gather stats */
				gwsetnormroutine (&lldata.gwdata, 0, 1, 0);
				gwstartnextfft (&lldata.gwdata, i < iters / 2);
				lucas_fixup (&lldata, p);
				gwsquare (&lldata.gwdata, lldata.lldata);
			} else if (type == 2) {		/* Thorough test */
				unsigned long j;
				for (j = 0; j < (i & 7); j++) {
					gwadd (&lldata.gwdata, lldata.lldata, lldata.lldata);
					lldata.units_bit = (lldata.units_bit+1) % p;
				}
				if ((i & 15) == 13) {
					gwadd3quick (&lldata.gwdata, lldata.lldata, lldata.lldata, t1);
					gwsub3quick (&lldata.gwdata, t1, lldata.lldata, lldata.lldata);
					gwadd3 (&lldata.gwdata, lldata.lldata, lldata.lldata, t1);
					gwsub3 (&lldata.gwdata, t1, lldata.lldata, lldata.lldata);
					gwaddsub4 (&lldata.gwdata, lldata.lldata, lldata.lldata, t1, t2);
					gwaddsub (&lldata.gwdata, t1, lldata.lldata);
					gwadd (&lldata.gwdata, t2, lldata.lldata);
				}
				lucas_fixup (&lldata, p);
				if ((i & 3) == 0) {
					gwsquare (&lldata.gwdata, lldata.lldata);
				} else if ((i & 3) == 1) {
					gwfft (&lldata.gwdata, lldata.lldata, lldata.lldata);
					gwfftfftmul (&lldata.gwdata, lldata.lldata, lldata.lldata, lldata.lldata);
				} else {
					gwfft (&lldata.gwdata, lldata.lldata, t1);
					gwfftmul (&lldata.gwdata, t1, lldata.lldata);
				}
			} else if (type == 3) {		/* Typical ECM run */
				lucas_fixup (&lldata, p);
				gwfftsub3 (&lldata.gwdata, t1, t2, t2);
				gwfft (&lldata.gwdata, lldata.lldata, lldata.lldata);
				gwfftfftmul (&lldata.gwdata, t2, lldata.lldata, t2);
				gwswap (t1, lldata.lldata);
				gwswap (t2, lldata.lldata);
			}

/* Keep track of the standard deviation - see Knuth vol 2 */

			if (i > iters_unchecked) {
				toterr += gw_get_maxerr (&lldata.gwdata);
				if (i == iters_unchecked + 1) {
					M = gw_get_maxerr (&lldata.gwdata);
					S = 0.0;
				} else {
					double	newM;
					newM = M + (gw_get_maxerr (&lldata.gwdata) - M) /
						   (i - iters_unchecked);
					S = S + (gw_get_maxerr (&lldata.gwdata) - M) * (gw_get_maxerr (&lldata.gwdata) - newM);
					M = newM;
				}
			}

/* Maintain range info */

			if (gw_get_maxerr (&lldata.gwdata) >= 0.300) ge_300++;
			if (gw_get_maxerr (&lldata.gwdata) >= 0.325) ge_325++;
			if (gw_get_maxerr (&lldata.gwdata) >= 0.350) ge_350++;
			if (gw_get_maxerr (&lldata.gwdata) >= 0.375) ge_375++;
			if (gw_get_maxerr (&lldata.gwdata) >= 0.400) ge_400++;

/* Maintain maximum error info */

			if (gw_get_maxerr (&lldata.gwdata) > maxerr) maxerr = gw_get_maxerr (&lldata.gwdata), maxerrcnt = 1;
			else if (gw_get_maxerr (&lldata.gwdata) == maxerr) maxerrcnt++;
			gw_clear_maxerr (&lldata.gwdata);

/* Maintain maximum suminp/sumout difference */

			if (fabs (gwsuminp (&lldata.gwdata, lldata.lldata) -
				  gwsumout (&lldata.gwdata, lldata.lldata)) > maxsumdiff) {
				maxsumdiff = fabs (gwsuminp (&lldata.gwdata, lldata.lldata) -
						   gwsumout (&lldata.gwdata, lldata.lldata));
			}

/* If the sum of the output values is an error (such as infinity) */
/* then raise an error.  For some reason these bad values are treated */
/* as zero by the C compiler.  There is probably a better way to */
/* check for this error condition. */

			if (gw_test_illegal_sumout (&lldata.gwdata)) {
				OutputBoth (thread_num, "Warning: ILLEGAL SUMOUT\n");
				dbltogw (&lldata.gwdata, 11.0, lldata.lldata);
				gw_clear_error (&lldata.gwdata);
			}

/* Check that the sum of the input numbers squared is approximately */
/* equal to the sum of unfft results. */

			if (gw_test_mismatched_sums (&lldata.gwdata)) {
				OutputBoth (thread_num, "Warning: SUMOUT MISMATCH\n");
				gw_clear_error (&lldata.gwdata);
			}

/* Abort if user demands it */

			stop_reason = stopCheck (thread_num);
			if (stop_reason) {
				lucasDone (&lldata);
				fclose (fd);
				return (stop_reason);
			}
		}

/* Generate residue and cleanup */

		generateResidue64 (&lldata, &reshi, &reslo);
		lucasDone (&lldata);

/* Output array of distributions of MAXERR */

		if (type == 1 || type == 3 || type == 4) {
			S = sqrt (S / (iters - iters_unchecked - 1));
			toterr /= iters - iters_unchecked;
			sprintf (buf, "avg: %6.6f, stddev: %6.6f, #stdev to 0.5: %6.6f\n",
				 toterr, S, (0.50 - toterr) / S);
			OutputBoth (thread_num, buf);
		}

/* Compare residue with correct residue from the input file */

		sprintf (buf, "%08X%08X", reshi, reslo);
		if (type <= 2 && _stricmp (res, buf)) {
			sprintf (buf, "Warning: Residue mismatch. Expected %s\n", res);
			OutputBoth (thread_num, buf);
		}

/* Output message */

		sprintf (buf, "Exp/iters: %lu/%lu, res: %08X%08X, maxerr: %6.6f/%lu, %lu/%lu/%lu/%lu/%lu, maxdiff: %9.9f/%9.9f\n",
			 p, iters, reshi, reslo, maxerr, maxerrcnt,
			 ge_300, ge_325, ge_350, ge_375, ge_400,
			 maxsumdiff, lldata.gwdata.MAXDIFF);
		OutputBoth (thread_num, buf);
		}
	}
	fclose (fd);

	return (0);
}

/* Test the factoring program */

int primeSieveTest (
	int	thread_num)
{
	fachandle facdata;
	char	buf[500];
	FILE	*fd;
	unsigned long p;
	int	stop_reason;
	uint32_t res, carryl, carryh;

/* Open factors file */

	fd = fopen ("factors", "r");

/* Loop until all the entire range is factored */

	while (fscanf (fd, "%ld", &p) && p) {
		unsigned long fachi, facmid, faclo;
		unsigned long i, pass;
		char fac[480];
		char *f;

/* What is the factor? */

		fscanf (fd, "%s", fac);
		fachi = facmid = faclo = 0;
		for (f = fac; *f; f++) {
			if (*f < '0' || *f > '9') continue;
			res = *f - '0';
			carryl = 0;
			muladdhlp (&res, &carryl, &carryh, faclo, 10);
			faclo = res;
			res = carryl;
			carryl = 0;
			muladdhlp (&res, &carryl, &carryh, facmid, 10);
			facmid = res;
			fachi = fachi * 10 + carryl;
			if (fachi >= 4194304 ||
			    (fachi >= 4096 && !(CPU_FLAGS & CPU_SSE2))) {
				sprintf (buf, "%ld%s factor too big.\n", p, fac);
				OutputBoth (thread_num, buf);
				goto nextp;
			}
		}

/* See if p is a prime */

		if (! isPrime (p)) {
			sprintf (buf, "%ld not a prime.\n", p);
			OutputBoth (thread_num, buf);
			goto nextp;
		}

/* Setup the factoring program */

		i = (fachi % 120 * 16 + facmid % 120 * 16 + faclo % 120) % 120;
		if (i == 1) pass = 0;
		else if (i == 7) pass = 1;
		else if (i == 17) pass = 2;
		else if (i == 23) pass = 3;
		else if (i == 31) pass = 4;
		else if (i == 41) pass = 5;
		else if (i == 47) pass = 6;
		else if (i == 49) pass = 7;
		else if (i == 71) pass = 8;
		else if (i == 73) pass = 9;
		else if (i == 79) pass = 10;
		else if (i == 89) pass = 11;
		else if (i == 97) pass = 12;
		else if (i == 103) pass = 13;
		else if (i == 113) pass = 14;
		else if (i == 119) pass = 15;
		else goto bad;
		stop_reason = factorSetup (thread_num, p, &facdata);
		if (stop_reason) {
			fclose (fd);
			return (stop_reason);
		}
		facdata.asm_data->FACHSW = fachi;
		facdata.asm_data->FACMSW = facmid;
		stop_reason = factorPassSetup (thread_num, pass, &facdata);
		if (stop_reason) {
			fclose (fd);
			factorDone (&facdata);
			return (stop_reason);
		}

/* Factor found, is it a match? */

		do {
			if (factorChunk (&facdata) != 2 &&
			    facdata.asm_data->FACHSW == fachi &&
			    facdata.asm_data->FACMSW == facmid &&
			    facdata.asm_data->FACLSW == faclo) {
				sprintf (buf, "%ld%s factored OK.\n", p, fac);
				OutputSomewhere (thread_num, buf);
				goto nextp;
			}
		} while (facdata.asm_data->FACMSW == facmid);

/* Uh oh. */

bad:		sprintf (buf, "%ld%s factor not found.\n", p, fac);
		OutputBoth (thread_num, buf);

/* If an escape key was hit, write out the results and return */

nextp:		stop_reason = stopCheck (thread_num);
		if (stop_reason) {
			fclose (fd);
			factorDone (&facdata);
			return (stop_reason);
		}
		p = 0;
	}

/* All done */

	fclose (fd);
	factorDone (&facdata);
	return (0);
}

/*********************/
/* Benchmarking code */
/*********************/

/* Time a few iterations of an LL test on a given exponent */

int primeTime (
	int	thread_num,
	unsigned long p,
	unsigned long iterations)
{
	struct PriorityInfo sp_info;
#define SAVED_LIMIT	10
	llhandle lldata;
	unsigned long i, j, saved, save_limit, num_threads;
	char	buf[80], fft_desc[100];
	double	time, saved_times[SAVED_LIMIT];
	int	days, hours, minutes, stop_reason;
	uint32_t *ASM_TIMERS;
	uint32_t best_asm_timers[32] = {0};
	double	timers[2];

/* Look for special values to run QA suites */

	if (p >= 9900 && p <= 9999) {
		sp_info.type = SET_PRIORITY_QA;
		sp_info.thread_num = thread_num;
		sp_info.aux_thread_num = 0;
		SetPriority (&sp_info);

		if (p >= 9994 && p <= 9999)
			return (lucas_QA (thread_num, 9999 - p));
		if (p == 9992)
			return (pminus1_QA (thread_num, &sp_info));
		if (p == 9991)
			return (ecm_QA (thread_num, &sp_info));
		if (p == 9990)
			return (primeSieveTest (thread_num));
		if (p >= 9900 && p <= 9919)
			return (test_randomly (thread_num, &sp_info));
		return (test_all_impl (thread_num, &sp_info));
	}

/* Set the process/thread priority */

	sp_info.type = SET_PRIORITY_BENCHMARKING;
	sp_info.thread_num = thread_num;
	sp_info.aux_thread_num = 0;
	SetPriority (&sp_info);

/* Loop through all possible num_thread values */

	for (num_threads = 1;
	     num_threads <= NUM_CPUS * CPU_HYPERTHREADS;
	     num_threads++) {

/* Clear all timers */

	clear_timers (timers, sizeof (timers) / sizeof (timers[0]));

/* Init the FFT code */

	gwinit (&lldata.gwdata);
	gwset_num_threads (&lldata.gwdata, num_threads);
	gwset_thread_callback (&lldata.gwdata, SetAuxThreadPriority);
	gwset_thread_callback_data (&lldata.gwdata, &sp_info);
	stop_reason = lucasSetup (thread_num, p, 0, &lldata);
	if (stop_reason) return (stop_reason);
	ASM_TIMERS = get_asm_timers (&lldata.gwdata);
	memset (ASM_TIMERS, 0, 32 * sizeof (uint32_t));

/* Output a message about the FFT length */

	gwfft_description (&lldata.gwdata, fft_desc);
	sprintf (buf, "Using %s\n", fft_desc);
	OutputStr (thread_num, buf);
	title (thread_num, "Timing");

/* Fill data space with random values. */

	generateRandomData (&lldata);

/* Do one squaring untimed, to prime the caches and start the */
/* post-FFT process going. */

	gwsetnormroutine (&lldata.gwdata, 0, ERRCHK != 0, 0);
	gwstartnextfft (&lldata.gwdata, TRUE);
	gwsquare (&lldata.gwdata, lldata.lldata);

/* Compute numbers in the lucas series */
/* Note that for reasons unknown, we've seen cases where printing out
/* the times on each iteration greatly impacts P4 timings. */

	save_limit = (p <= 4000000) ? SAVED_LIMIT : 1;
	for (i = 0, saved = 0; i < iterations; i++) {

/* Time a single squaring */

		start_timer (timers, 0);
		gwsquare (&lldata.gwdata, lldata.lldata);
		end_timer (timers, 0);
		timers[1] += timers[0];
		saved_times[saved++] = timers[0];
		timers[0] = 0;

/* Remember the best asm timers (used when I'm optimizing assembly code) */

		for (j = 0; j < 32; j++)
			if (i == 0 || ASM_TIMERS[j] < best_asm_timers[j])
				best_asm_timers[j] = ASM_TIMERS[j];

/* Output timer squaring times */

		if (saved == save_limit || i == iterations - 1) {
			for (j = 0; j < saved; j++) {
				sprintf (buf, "p: %u.  Time: ", p);
				timers[0] = saved_times[j];
				print_timer (timers, 0, buf, TIMER_MS | TIMER_NL | TIMER_CLR);
				OutputStr (thread_num, buf);
			}
			saved = 0;
		}

/* Abort early if so requested */

		stop_reason = stopCheck (thread_num);
		if (stop_reason) {
			lucasDone (&lldata);
			return (stop_reason);
		}
	}
	lucasDone (&lldata);
	time = timer_value (timers, 1);

/* Print an estimate for how long it would take to test this number */

	sprintf (buf, "Iterations: %d.  Total time: ", iterations);
	print_timer (timers, 1, buf, TIMER_NL | TIMER_CLR);
	OutputStr (thread_num, buf);
	time = time * p / iterations;
	days = (int) (time / 86400.0); time -= (double) days * 86400.0;
	hours = (int) (time / 3600.0); time -= (double) hours * 3600.0;
	minutes = (int) (time / 60.0);
	strcpy (buf, "Estimated time to complete this exponent: ");
	sprintf (buf+strlen(buf), days == 1 ? "%d day, " : "%d days, ", days);
	sprintf (buf+strlen(buf), hours == 1 ? "%d hour, " : "%d hours, ", hours);
	sprintf (buf+strlen(buf), minutes == 1 ? "%d minute.\n" : "%d minutes.\n", minutes);
	OutputStr (thread_num, buf);

/* I use these assembly language timers to time various chunks of */
/* assembly code.  Print these timers out. */

	for (i = 0; i < 32; i++) {
		sprintf (buf, "timer %d: %d\n", i, (int) best_asm_timers[i]);
		if (best_asm_timers[i]) OutputBoth (thread_num, buf);
	}

/* Loop through all possible thread counts */

	}

/* All done */

	return (0);
}

#define BENCH1 "Your timings will be written to the results.txt file.\n"
#define BENCH2 "Compare your results to other computers at http://www.mersenne.org/bench.htm\n"

int factorBench (
	int	thread_num,
	struct primenetBenchmarkData *pkt)
{
	fachandle facdata;
	unsigned long num_lengths, i, j;
	double	best_time;
	char	buf[512];
	int	bit_lengths[] = {58, 59, 60, 61, 62, 63, 64, 65, 66, 67};
	int	res, stop_reason;
	double	timers[2];

/* Loop over all trial factor lengths */

	num_lengths = sizeof (bit_lengths) / sizeof (int);
	for (i = 0; i < num_lengths; i++) {

/* Initialize for this bit length. */

		stop_reason = factorSetup (thread_num, 35000011, &facdata);
		if (stop_reason) return (stop_reason);
		if (bit_lengths[i] < 64) {
			facdata.asm_data->FACHSW = 0;
			facdata.asm_data->FACMSW = 1L << (bit_lengths[i]-32);
		} else {
			facdata.asm_data->FACHSW = 1L << (bit_lengths[i]-64);
			facdata.asm_data->FACMSW = 0;
		}
		stop_reason = factorPassSetup (thread_num, 0, &facdata);
		if (stop_reason) return (stop_reason);

/* Output start message for this bit length */

		sprintf (buf, "Timing trial factoring of M35000011 with %d bit length factors.  ", bit_lengths[i]);
		OutputStr (thread_num, buf);

/* Do one "iteration" untimed, to prime the caches. */

		res = factorChunk (&facdata);

/* Time 10 iterations. Take best time. */

		for (j = 0; j < 10; j++) {
			stop_reason = stopCheck (thread_num);
			if (stop_reason) {
				OutputStrNoTimeStamp (thread_num, "\n");
				OutputStr (thread_num, "Execution halted.\n");
				factorDone (&facdata);
				return (stop_reason);
			}
			clear_timers (timers, sizeof (timers) / sizeof (timers[0]));
			start_timer (timers, 0);
			res = factorChunk (&facdata);
			end_timer (timers, 0);
			if (j == 0 || timers[0] < best_time) best_time = timers[0];
		}
		factorDone (&facdata);

/* Print the best time for this bit length.  Take into account that */
/* X86_64 factoring code does 3 times as much work (bigger sieve). */

#ifdef X86_64
		best_time = best_time / 3;
#endif
		timers[0] = best_time;
		strcpy (buf, "Best time: ");
		print_timer (timers, 0, buf, TIMER_NL | TIMER_MS);
		OutputStrNoTimeStamp (thread_num, buf);
		sprintf (buf, "Best time for %d bit trial factors: ", bit_lengths[i]);
		print_timer (timers, 0, buf, TIMER_NL | TIMER_MS);
		writeResults (buf);

/* Accumulate best times to send to the server */

		if (pkt->num_data_points < PRIMENET_BENCH_MAX_DATAPOINTS) {
			sprintf (pkt->data_points[pkt->num_data_points].bench,
				 "TF%d", bit_lengths[i]);
			pkt->data_points[pkt->num_data_points].timing = timer_value (timers, 0);
			pkt->num_data_points++;
		}

	}
	return (0);
}

/* Time a few iterations of many FFT lengths */

int primeBench (
	int	thread_num)
{
	struct PriorityInfo sp_info;
	llhandle lldata;
	unsigned long num_lengths, i, ii, j, iterations;
	double	best_time;
	char	buf[512];
	unsigned int threads;
	int	fft_lengths[] = {4, 5, 6, 7, 8, 10, 12, 14, 16, 20, 24, 28, 32, 40, 48, 56, 64, 80, 96, 112, 128, 160, 192, 224, 256, 320, 384, 448, 512, 640, 768, 896, 1024, 1280, 1536, 1792, 2048, 2560, 3072, 3584, 4096, 5120, 6144, 7168, 8192, 10240, 12288, 14336, 16384, 20480, 24576, 28672, 32768};
	int	all_bench, plus1, stop_reason;
	double	timers[2];
	struct primenetBenchmarkData pkt;

/* Init */

	memset (&pkt, 0, sizeof (pkt));
	strcpy (pkt.computer_guid, COMPUTER_GUID);

/* Set the process/thread priority */

	sp_info.type = SET_PRIORITY_BENCHMARKING;
	sp_info.thread_num = thread_num;
	sp_info.aux_thread_num = 0;
	SetPriority (&sp_info);

/* Output startup message */

	title (thread_num, "Benchmarking");
	OutputStr (thread_num, BENCH1);
	OutputBoth (thread_num, BENCH2);

/* Output to the results file a full CPU description */

	getCpuDescription (buf, 1);
	writeResults (buf);
#ifdef X86_64
	sprintf (buf, "Prime95 64-bit version %s, RdtscTiming=%d\n",
		 VERSION, RDTSC_TIMING);
#else
	sprintf (buf, "Prime95 32-bit version %s, RdtscTiming=%d\n",
		 VERSION, RDTSC_TIMING);
#endif
	writeResults (buf);

/* Loop over all all possible multithread possibilities */

	for (threads = 1; threads <= NUM_CPUS * CPU_HYPERTHREADS; threads++) {
	  if (threads > CPU_HYPERTHREADS && threads % CPU_HYPERTHREADS != 0)
	    continue;
	  if (threads > 1) {
	    if (! (CPU_FLAGS & CPU_SSE2)) continue;  // Only SSE2 code supports multi-threaded FFTs
	    if (CPU_HYPERTHREADS == 1) {
	      sprintf (buf, "Timing FFTs using %d threads.\n", threads);
	      OutputBoth (thread_num, buf);
	    } else {
	      sprintf (buf, "Timing FFTs using %d threads on %d physical CPUs.\n", threads, (threads + CPU_HYPERTHREADS - 1) / CPU_HYPERTHREADS);
	      OutputBoth (thread_num, buf);
	    }
	  }

/* Loop over all FFT lengths */

	  all_bench = IniGetInt (INI_FILE, "AllBench", 0);
	  if (IniGetInt (INI_FILE, "FullBench", 0)) {
	    i = 0;
	    num_lengths = sizeof (fft_lengths) / sizeof (int);
	  } else {
	    i = 30;
	    num_lengths = sizeof (fft_lengths) / sizeof (int) - 8;
	    if (! (CPU_FLAGS & CPU_SSE2)) num_lengths -= 4;
	  }
	  for ( ; i < num_lengths; i++) {
	    for (plus1 = 0; plus1 <= 1; plus1++) {
	      for (ii = 1; ; ii++) {

/* Initialize for this FFT length.  Compute the number of iterations to */
/* time.  This is based on the fact that it doesn't take too long for */
/* my 1400 MHz P4 to run 10 iterations of a 1792K FFT. */

		gwinit (&lldata.gwdata);
		gwset_num_threads (&lldata.gwdata, threads);
		gwset_thread_callback (&lldata.gwdata, SetAuxThreadPriority);
		gwset_thread_callback_data (&lldata.gwdata, &sp_info);
		if (all_bench) lldata.gwdata.bench_pick_nth_fft = ii;
		stop_reason = lucasSetup (thread_num, fft_lengths[i] * 1024 * 17 + 1, fft_lengths[i] * 1024 + plus1, &lldata);
		if (stop_reason) break;
		iterations = (unsigned long) (10 * 1792 * CPU_SPEED / 1400 / fft_lengths[i]);
		if (iterations < 10) iterations = 10;
		if (iterations > 100) iterations = 100;

/* Output start message for this FFT length */

		sprintf (buf, "Timing %d iterations at %dK%s FFT length.  ",
			 iterations, fft_lengths[i],
			 plus1 ? " all-complex" : "");
		OutputStr (thread_num, buf);

/* Fill data space with random values. */

		generateRandomData (&lldata);

/* Do one squaring untimed, to prime the caches and start the */
/* POSTFFT optimization going. */

		gwsetnormroutine (&lldata.gwdata, 0, 0, 0);
		gwstartnextfft (&lldata.gwdata, TRUE);
		gwsquare (&lldata.gwdata, lldata.lldata);

/* Compute numbers in the lucas series */
/* Note that for reasons unknown, we've seen cases where printing out
/* the times on each iteration greatly impacts P4 timings. */

		for (j = 0; j < iterations; j++) {
			stop_reason = stopCheck (thread_num);
			if (stop_reason) {
				OutputStrNoTimeStamp (thread_num, "\n");
				OutputStr (thread_num, "Execution halted.\n");
				lucasDone (&lldata);
				return (stop_reason);
			}
			clear_timers (timers, sizeof (timers) / sizeof (timers[0]));
			start_timer (timers, 0);
			gwsquare (&lldata.gwdata, lldata.lldata);
			end_timer (timers, 0);
			if (j == 0 || timers[0] < best_time) best_time = timers[0];
		}
		lucasDone (&lldata);

/* Print the best time for this FFT length */

		timers[0] = best_time;
		strcpy (buf, "Best time: ");
		print_timer (timers, 0, buf, TIMER_NL | TIMER_MS);
		OutputStrNoTimeStamp (thread_num, buf);
		if (all_bench)
			sprintf (buf,
				 "Time FFTlen=%dK%s, Levels2=%d, clm=%d: ",
				 fft_lengths[i], plus1 ? " all-complex" : "",
				 lldata.gwdata.PASS2_LEVELS,
				 lldata.gwdata.PASS1_CACHE_LINES / 2);
		else
			sprintf (buf, "Best time for %dK%s FFT length: ",
				 fft_lengths[i], plus1 ? " all-complex" : "");
		print_timer (timers, 0, buf, TIMER_NL | TIMER_MS);
		writeResults (buf);

/* Accumulate best times to send to the server */

		if (!all_bench &&
		    pkt.num_data_points < PRIMENET_BENCH_MAX_DATAPOINTS) {
			if (threads == 1)
				sprintf (pkt.data_points[pkt.num_data_points].bench,
					 "FFT%dK", fft_lengths[i]);
			else
				sprintf (pkt.data_points[pkt.num_data_points].bench,
					 "FFT%dK %dT", fft_lengths[i], threads);
			pkt.data_points[pkt.num_data_points].timing = timer_value (timers, 0);
			pkt.num_data_points++;
		}

/* Do next FFT implementation (if all_bench set) or next FFT size (if */
/* all_bench is not set) */

		if (!all_bench) break;
	    }

/* Do all-complex implementation (if all_bench set) or next FFT size (if */
/* all_bench is not set) */

	    if (!all_bench) break;
	  }
	}
}

/* Now benchmark the trial factoring code */

	stop_reason = factorBench (thread_num, &pkt);
	if (stop_reason) return (stop_reason);
	OutputStr (thread_num, "Benchmark complete.\n");

/* Finally, send the benchmark data to the server. */

//bug - send bench data to server. (checkbox to allow sending data to server?)
//only do this if guid is registered? Or should we auto-register computer
//under ANONYMOUS userid for stress-testers.

	spoolMessage (PRIMENET_BENCHMARK_DATA, &pkt);
	return (0);
}

/****************************/
/* Probable Prime Test code */
/****************************/

/* Write intermediate PRP results to a file */
/* The PRP save file format is: */
/*	u32		magic number  (different for ll, p-1, prp, tf, ecm) */
/*	u32		version number */
/*	double		pct complete */
/*	char(11)	stage */
/*	char(1)		pad */
/*	u32		checksum of following data */
/*	u32		error_count */
/*	u32		iteration counter */
/*	gwnum		FFT data (u32 len, array u32s) */

#define PRP_MAGICNUM		0x87f2a91b
#define PRP_VERSION		1

int writePRPSaveFile (
	gwhandle *gwdata,
	gwnum	x,
	char	*filename,
	struct work_unit *w,
	unsigned long counter,
	unsigned long error_count)
{
	int	fd;
	unsigned long sum = 0;

/* If we are allowed to create multiple intermediate files, then */
/* write to a file called rXXXXXXX. */

	if (TWO_BACKUP_FILES && strlen (filename) == 8)
		filename[0] = 'r';

/* Now save to the intermediate file */

	fd = _open (filename, _O_BINARY | _O_WRONLY | _O_TRUNC | _O_CREAT, CREATE_FILE_ACCESS);
	if (fd < 0) return (FALSE);

	if (!write_header (fd, PRP_MAGICNUM, PRP_VERSION, w)) goto err;

	if (!write_long (fd, error_count, &sum)) goto err;
	if (!write_long (fd, counter, &sum)) goto err;
	if (!write_gwnum (fd, gwdata, x, &sum)) goto err;

	if (!write_checksum (fd, sum)) goto err;

	_commit (fd);
	_close (fd);

/* Now rename the intermediate files */

	if (TWO_BACKUP_FILES && strlen (filename) == 8) {
		char	backupname[16];
		strcpy (backupname, filename);
		backupname[0] = 'q'; filename[0] = 'p';
		_unlink (backupname);
		 rename (filename, backupname);
		backupname[0] = 'r';
		rename (backupname, filename);
	}

	return (TRUE);

/* An error occured.  Delete the current file. */

err:	_close (fd);
	_unlink (filename);
	return (FALSE);
}

/* Read the data portion of an intermediate PRP results file */

int readPRPSaveFile (
	gwhandle *gwdata,
	gwnum	x,
	char	*filename,
	struct work_unit *w,
	unsigned long *counter,
	unsigned long *error_count)
{
	int	fd;
	unsigned long sum, filesum, version;

	fd = _open (filename, _O_BINARY | _O_RDONLY);
	if (fd <= 0) return (FALSE);

	if (!read_magicnum (fd, PRP_MAGICNUM)) goto err;
	if (!read_header (fd, &version, w, &filesum)) goto err;
	if (version != PRP_VERSION) goto err;

	sum = 0;
	if (!read_long (fd, error_count, &sum)) goto err;
	if (!read_long (fd, counter, &sum)) goto err;

	if (!read_gwnum (fd, gwdata, x, &sum)) goto err;

	if (filesum != sum) goto err;
	_close (fd);
	return (TRUE);
err:	_close (fd);
	return (FALSE);
}

/* Do a PRP test */

int prp (
	int	thread_num,		/* Worker thread number */
	struct PriorityInfo *sp_info,	/* SetPriority information */
	struct work_unit *w,		/* Worktodo entry */
	int	pass)			/* PrimeContinue pass */
{
	gwhandle gwdata;
	gwnum	x;
	giant	N, tmp;
	int	first_iter_msg, res, stop_reason;
	int	echk, saving, near_fft_limit, sleep5, isProbablePrime;
	unsigned long Nlen, counter, iters, error_count;
	int	slow_iteration_count;
	double	timers[2];
	double	inverse_Nlen;
	double	reallyminerr = 1.0;
	double	reallymaxerr = 0.0;
	double	best_iteration_time;
	char	filename[32];
	char	buf[160], fft_desc[100], res64[17];
	unsigned long last_counter = 0;		/* Iteration of last error */
	int	maxerr_recovery_mode = 0;	/* Big roundoff err rerun */
	double	last_suminp = 0.0;
	double	last_sumout = 0.0;
	double	last_maxerr = 0.0;
	char	string_rep[80];
	int	string_rep_truncated;

/* See if this exponent needs P-1 factoring.  We treat P-1 factoring */
/* that is part of a PRP test as priority work (done in pass 1). */

	if (! w->pminus1ed) {
		stop_reason = pfactor (thread_num, sp_info, w, pass);
		if (stop_reason) {
			if (pass == 3 && stop_reason == STOP_NOT_ENOUGH_MEM)
				stop_reason = 0;
			else
				return (stop_reason);
		}
	}

/* Done with pass 1 priority work.  Return to do more priority work. */

	if (pass == 1) return (0);

/* Figure out which FFT size we should use */

	stop_reason = pick_fft_size (thread_num, w);
	if (stop_reason) return (stop_reason);

/* Make sure the first-time user runs a successful self-test. */
/* The one-hour self-test may have been useful when it was first introduced */
/* but I think it now does little to catch buggy machines (they eventually */
/* work OK for an hour) and does create user confusion and annoyance. */

#ifdef ONE_HOUR_SELF_TEST
	stop_reason = selfTest (thread_num, sp_info, w);
	if (stop_reason) return (stop_reason);
#endif

/* Init the FFT code for squaring modulo k*b^n+c. */

begin:	gwinit (&gwdata);
	gwset_num_threads (&gwdata, THREADS_PER_TEST[thread_num]);
	gwset_thread_callback (&gwdata, SetAuxThreadPriority);
	gwset_thread_callback_data (&gwdata, sp_info);
	gwset_specific_fftlen (&gwdata, w->forced_fftlen);
	res = gwsetup (&gwdata, w->k, w->b, w->n, w->c);

/* If we were unable to init the FFT code, then print an error message */
/* and return an error code. */

	if (res) {
		sprintf (buf, "PRP cannot initialize FFT code, errcode=%d\n", res);
		OutputBoth (thread_num, buf);
		return (STOP_FATAL_ERROR);
	}

/* Record the amount of memory being used by this thread. */

	set_memory_usage (thread_num, 0,
			  gwmemused (&gwdata) + gwnum_size (&gwdata));

/* Allocate memory for the PRP test */

	x = gwalloc (&gwdata);
	if (x == NULL) {
		gwdone (&gwdata);
		OutputStr (thread_num, "Error allocating memory for FFT data.\n");
		return (STOP_OUT_OF_MEM);
	}

/* Format the string representation of the test number */

	strcpy (string_rep, gwmodulo_as_string (&gwdata));
	if (w->known_factors == NULL)
		string_rep_truncated = FALSE;
	else if (strlen (w->known_factors) < 40) {
		char	*p;
		strcat (string_rep, "/");
		strcat (string_rep, w->known_factors);
		while ((p = strchr (string_rep, ',')) != NULL) *p = '/';
		string_rep_truncated = FALSE;
	} else {
		strcat (string_rep, "/known_factors");
		string_rep_truncated = TRUE;
	}

/* Init the title */

	sprintf (buf, "PRP %s", string_rep);
	title (thread_num, buf);

/* Loop reading from save files (and backup save files) */

readloop:
	tempFileName (w, filename);

/* Read a PRP save file.  On error try the backup save file. */

	if (continuationFileExists (thread_num, filename)) {
		if (! readPRPSaveFile (&gwdata, x, filename, w, &counter, &error_count)) {
			sprintf (buf, READFILEERR, filename);
			OutputBoth (thread_num, buf);
			_unlink (filename);
			goto readloop;
		}
		first_iter_msg = TRUE;
	}

/* Start off with the 1st PRP squaring */

	else {
		dbltogw (&gwdata, 3.0, x);
		counter = 0;
		error_count = 0;
		first_iter_msg = FALSE;
	}

/* Output a message saying we are starting/resuming the PRP test. */
/* Also output the FFT length. */

	gwfft_description (&gwdata, fft_desc);
	sprintf (buf, "%s PRP test of %s using %s\n",
		 (counter == 0) ? "Starting" : "Resuming",
		 string_rep, fft_desc);
	OutputStr (thread_num, buf);

/* Compute the number we are testing. */

	stop_reason = setN (&gwdata, thread_num, w, &N);
	if (stop_reason) goto exit;

/* Subtract 1 from N to compute a^(N-1) mod N.  Get the exact bit length */
/* of the number.  We will perform bitlen(N)-1 squarings for the PRP test. */

	iaddg (-1, N);
	Nlen = bitlen (N);

/* Hyperthreading backoff is an option to pause the program when iterations */
/* take longer than usual.  This is useful on hyperthreaded machines so */
/* that prime95 doesn't steal cycles from a foreground task, thus hurting */
/* the computers responsiveness. */

	best_iteration_time = 1.0e50;
	slow_iteration_count = 0;

/* Clear all timers */

	clear_timers (timers, sizeof (timers) / sizeof (timers[0]));

/* Init vars for Test/Status and CommunicateWithServer */

	strcpy (w->stage, "PRP");
	inverse_Nlen = 1.0 / (double) (Nlen - 1);
	w->pct_complete = (double) counter * inverse_Nlen;

/* If we are near the maximum exponent this fft length can test, then we */
/* will error check all iterations */

	near_fft_limit = exponent_near_fft_limit (&gwdata);

/* Do the PRP test */

//#define CHECK_ITER
#ifdef CHECK_ITER
{giant t1, t2;
t1 = popg (&gwdata.gdata, ((unsigned long) gwdata.bit_length >> 4) + 5);
t2 = popg (&gwdata.gdata, ((unsigned long) gwdata.bit_length >> 4) + 5);
gwtogiant (&gwdata, x, t1);
#endif
	gwsetmulbyconst (&gwdata, 3);
	iters = 0;
	while (counter < Nlen - 1) {

/* On first iteration create a save file so that writeNewErrorCount */
/* can properly keep track of error counts. */
/* Also save right after we pass an errored iteration and several */
/* iterations before retesting an errored iteration so that we don't */
/* have to backtrack very far to do a careful_iteration	(we don't do the */
/* iteration immediately before because on the P4 a save operation will */
/* change the FFT data and make the error non-reproducible. */
/* Error check the first and last 50 iterations, before writing an */
/* intermediate file (either user-requested stop or a */
/* 30 minute interval expired), and every 128th iteration. */

		stop_reason = stopCheck (thread_num);
		saving = stop_reason ||
			 (counter == 0 && Nlen > 1500000) ||
			 counter == last_counter-8 ||
			 counter == last_counter ||
			 testSaveFilesFlag (thread_num);
		echk = saving || near_fft_limit || ERRCHK ||
			counter < 50 || counter >= Nlen-51 ||
			((counter & 127) == 0);
		gw_clear_maxerr (&gwdata);

/* Do one PRP iteration */

		timers[1] = 0.0;
		start_timer (timers, 1);

/* Process this bit.  Use square carefully the first and last 30 iterations. */
/* This should avoid any pathological non-random bit pattterns.  Also square */
/* carefully during an error recovery. This will protect us from roundoff */
/* errors up to 0.6. */

		gwstartnextfft (&gwdata,
				!saving && !maxerr_recovery_mode &&
				counter > 35 && counter < Nlen-35 &&
				(INTERIM_FILES == 0 ||
				 (counter+1) % INTERIM_FILES > 0) &&
				(INTERIM_RESIDUES == 0 ||
				 (counter+1) % INTERIM_RESIDUES > 0));
#ifdef CHECK_ITER
squareg (t1);
if (bitval (N, Nlen-2-counter)) imulg (3, t1);
specialmodg (&gwdata, t1);
if (w->known_factors) {	iaddg (1, N); modg (N, t1); iaddg (-1, N); }
gwstartnextfft (&gwdata, 0);
echk=1;
#endif
		if (bitval (N, Nlen-2-counter)) {
			gwsetnormroutine (&gwdata, 0, echk, 1);
		} else {
			gwsetnormroutine (&gwdata, 0, echk, 0);
		}
		if (maxerr_recovery_mode && counter == last_counter) {
			gwsquare_carefully (&gwdata, x);
			maxerr_recovery_mode = 0;
			echk = 0;
		} else if (counter < 30 || counter > Nlen-32)
			gwsquare_carefully (&gwdata, x);
		else
			gwsquare (&gwdata, x);

#ifdef CHECK_ITER
gwtogiant (&gwdata, x, t2);
if (w->known_factors) {	iaddg (1, N); modg (N, t2); iaddg (-1, N); }
if (gcompg (t1, t2) != 0)
OutputStr (thread_num, "Iteration failed.\n");
//if (counter == 100) counter = Nlen-2;
#endif

/* End iteration timing and increase count of iterations completed */

		end_timer (timers, 1);
		timers[0] += timers[1];
		iters++;

/* If the sum of the output values is an error (such as infinity) */
/* then raise an error. */

		if (gw_test_illegal_sumout (&gwdata)) {
			sprintf (buf, ERRMSG0, counter+1, Nlen-1, ERRMSG1A);
			OutputBoth (thread_num, buf);
			inc_error_count (2, &error_count);
			sleep5 = TRUE;
			goto restart;
		}

/* Check that the sum of the input numbers squared is approximately */
/* equal to the sum of unfft results.  Since this check may not */
/* be perfect, check for identical results after a restart. */

		if (gw_test_mismatched_sums (&gwdata)) {
			if (counter == last_counter &&
			    gwsuminp (&gwdata, x) == last_suminp &&
			    gwsumout (&gwdata, x) == last_sumout) {
				OutputBoth (thread_num, ERROK);
				inc_error_count (3, &error_count);
				gw_clear_error (&gwdata);
			} else {
				char	msg[100];
				sprintf (msg, ERRMSG1B,
					 gwsuminp (&gwdata, x),
					 gwsumout (&gwdata, x));
				sprintf (buf, ERRMSG0, counter+1, Nlen-1, msg);
				OutputBoth (thread_num, buf);
				last_counter = counter;
				last_suminp = gwsuminp (&gwdata, x);
				last_sumout = gwsumout (&gwdata, x);
				inc_error_count (0, &error_count);
				sleep5 = TRUE;
				goto restart;
			}
		}

/* Check for excessive roundoff error.  If round off is too large, repeat */
/* the iteration to see if this was a hardware error.  If it was repeatable */
/* then repeat the iteration using a safer, slower method.  This can */
/* happen when operating near the limit of an FFT. */

		if (echk && gw_get_maxerr (&gwdata) >= 0.40625) {
			if (counter == last_counter &&
			    gw_get_maxerr (&gwdata) == last_maxerr) {
				OutputBoth (thread_num, ERROK);
				inc_error_count (3, &error_count);
				gw_clear_error (&gwdata);
				OutputBoth (thread_num, ERRMSG5);
				maxerr_recovery_mode = 1;
				sleep5 = FALSE;
				goto restart;
			} else {
				char	msg[100];
				sprintf (msg, ERRMSG1C, gw_get_maxerr (&gwdata));
				sprintf (buf, ERRMSG0, counter+1, Nlen-1, msg);
				OutputBoth (thread_num, buf);
				last_counter = counter;
				last_maxerr = gw_get_maxerr (&gwdata);
				inc_error_count (1, &error_count);
				sleep5 = FALSE;
				goto restart;
			}
		}

/* Update counter, percentage complete, and maximum round-off error */

		counter++;
		w->pct_complete = (double) counter * inverse_Nlen;
		if (ERRCHK) {
			if (counter > 30 &&
			    gw_get_maxerr (&gwdata) < reallyminerr)
				reallyminerr = gw_get_maxerr (&gwdata);
			if (gw_get_maxerr (&gwdata) > reallymaxerr)
				reallymaxerr = gw_get_maxerr (&gwdata);
		}

/* Print a message every so often */

		if (counter % ITER_OUTPUT == 0 || first_iter_msg) {
			char	fmt_mask[80];
			double	pct;
			pct = trunc_percent (w->pct_complete);
			sprintf (fmt_mask, "%%.%df%%%% of %%s", PRECISION);
			sprintf (buf, fmt_mask, pct, string_rep);
			title (thread_num, buf);
			sprintf (fmt_mask,
				 "Iteration: %%ld / %%ld [%%.%df%%%%]",
				 PRECISION);
			sprintf (buf, fmt_mask, counter, Nlen-1, pct);
			if (ERRCHK && counter > 30) {
				sprintf (buf+strlen(buf),
					 ".  Round off: %10.10f to %10.10f",
					 reallyminerr, reallymaxerr);
			}
			if (first_iter_msg) {
				strcat (buf, ".\n");
				clear_timer (timers, 0);
			} else {
				strcat (buf, ".  Per iteration time: ");
				divide_timer (timers, 0, iters);
				print_timer (timers, 0, buf,
					     TIMER_NL | TIMER_OPT_CLR);
			}
			OutputStr (thread_num, buf);
			if (!CUMULATIVE_TIMING) iters = 0;
			first_iter_msg = FALSE;
		}

/* Print a results file message every so often */

		if (counter % ITER_OUTPUT_RES == 0 || (NO_GUI && stop_reason)) {
			sprintf (buf, "Iteration %ld / %ld\n", counter, Nlen-1);
			writeResults (buf);
		}

/* Write results to a file every DISK_WRITE_TIME minutes */
/* On error, retry in 10 minutes (it could be a temporary */
/* disk-full situation) */

		if (saving) {
			if (! writePRPSaveFile (&gwdata, x, filename, w,
						counter, error_count)) {
				sprintf (buf, WRITEFILEERR, filename);
				OutputBoth (thread_num, buf);
			}
		}

/* If an escape key was hit, write out the results and return */

		if (stop_reason) {
			char	fmt_mask[80];
			sprintf (fmt_mask,
				 "Stopping PRP test of %%s at iteration %%ld [%%.%df%%%%]\n",
				 PRECISION);
			sprintf (buf, fmt_mask, string_rep,
				 counter, trunc_percent (w->pct_complete));
			OutputStr (thread_num, buf);
			goto exit;
		}

/* Output the 64-bit residue at specified interims.  Also output the */
/* residues for the next two iterations so that we can compare our */
/* residues to programs that start counter at zero or one. */

		if (INTERIM_RESIDUES && counter % INTERIM_RESIDUES == 0) {
			tmp = popg (&gwdata.gdata, ((unsigned long) gwdata.bit_length >> 5) + 5);
			gwtogiant (&gwdata, x, tmp);
			if (w->known_factors) {	iaddg (1, N); modg (N, tmp); iaddg (-1, N); }
			sprintf (buf, 
				 "%s interim Wd%d residue %08lX%08lX at iteration %ld\n",
				 string_rep, PORT, tmp->n[1], tmp->n[0], counter);
			OutputBoth (thread_num, buf);
			pushg (&gwdata.gdata, 1);
		}

/* Write a save file every INTERIM_FILES iterations. */

		if (INTERIM_FILES && counter % INTERIM_FILES == 0) {
			char	interimfile[20];
			sprintf (interimfile, "%.8s.%03d",
				 filename, counter / INTERIM_FILES);
			writePRPSaveFile (&gwdata, x, interimfile, w, counter,
					  error_count);
		}

/* If ten iterations take 40% longer than a typical iteration, then */
/* assume a foreground process is running and sleep for a short time */
/* to give the foreground process more CPU time.  Even though a foreground */
/* process runs at higher priority, hyperthreading will cause this */
/* program to run at an equal priority, hurting responsiveness. */

		if (HYPERTHREADING_BACKOFF && Nlen > 1500000) {
			if (timers[1] < best_iteration_time)
				best_iteration_time = timers[1];
			if (timers[1] > 1.40 * best_iteration_time) {
				if (slow_iteration_count == 10) {
					sprintf (buf, "Pausing %d seconds.\n",
						 HYPERTHREADING_BACKOFF);
					OutputStr (thread_num, buf);
					Sleep (HYPERTHREADING_BACKOFF * 1000);
				}
				slow_iteration_count++;
			} else
				slow_iteration_count = 0;
		}
	}
#ifdef CHECK_ITER
pushg(&gwdata.gdata, 2);}
#endif

/* See if we've found a probable prime.  If not, format a 64-bit residue. */

	tmp = popg (&gwdata.gdata, ((unsigned long) gwdata.bit_length >> 5) + 5);
	gwtogiant (&gwdata, x, tmp);
	if (w->known_factors) {	iaddg (1, N); modg (N, tmp); iaddg (-1, N); }
	isProbablePrime = isone (tmp);
	if (!isProbablePrime) {
		sprintf (res64, "%08lX%08lX", tmp->n[1], tmp->n[0]);
	}
	pushg (&gwdata.gdata, 1);
	gwfree (&gwdata, x);

/* Print results. */

	if (isProbablePrime)
		sprintf (buf, "%s is a probable prime! Wd%d: %08lX,%08lX\n",
			 string_rep, PORT, SEC1 (w->n), error_count);
	else
		sprintf (buf, "%s is not prime.  RES64: %s. Wd%d: %08lX,%08lX\n",
			 string_rep, res64, PORT, SEC1 (w->n), error_count);
	OutputStr (thread_num, buf);
	formatMsgForResultsFile (buf, w);

/* Update the output file */

	if ((isProbablePrime && IniGetInt (INI_FILE, "OutputPrimes", 1)) ||
	    (!isProbablePrime && IniGetInt (INI_FILE, "OutputComposites", 1)))
		writeResults (buf);

//if (ERRCHK) {
//	sprintf (buf, "Round off: %10.10f to %10.10f\n", reallyminerr, reallymaxerr);
//	OutputBoth (thread_num, buf);
//}

/* Output results to the server */

	{
		struct primenetAssignmentResult pkt;
		memset (&pkt, 0, sizeof (pkt));
		strcpy (pkt.computer_guid, COMPUTER_GUID);
		strcpy (pkt.assignment_uid, w->assignment_uid);
		strcpy (pkt.message, buf);
		pkt.result_type =
			isProbablePrime ? PRIMENET_AR_PRP_PRIME : PRIMENET_AR_PRP_RESULT;
		pkt.k = w->k;
		pkt.b = w->b;
		pkt.n = w->n;
		pkt.c = w->c;
		strcpy (pkt.residue, res64);
		sprintf (pkt.error_count, "%08lX", error_count);
		pkt.fftlen = gwfftlen (&gwdata);
		pkt.done = TRUE;
		spoolMessage (PRIMENET_ASSIGNMENT_RESULT, &pkt);
	}

/* Print known factors */

	if (string_rep_truncated) {
		char	*bigbuf;
		bigbuf = (char *) malloc (strlen (w->known_factors) + 100);
		if (bigbuf != NULL) {
			sprintf (bigbuf,
				 "Known factors used for PRP test were: %s\n",
				 w->known_factors);
			OutputBoth (thread_num, bigbuf);
			free (bigbuf);
		}
	}

/* Delete the continuation files. */

	_unlink (filename);
	filename[0] = 'q';
	_unlink (filename);

/* Return work unit completed stop reason */

	stop_reason = STOP_WORK_UNIT_COMPLETE;

/* Cleanup and exit */

exit:	gwdone (&gwdata);
	free (N);
	return (stop_reason);

/* An error occured, output a message saying we are restarting, sleep, */
/* then try restarting at last save point. */

restart:if (sleep5) OutputBoth (thread_num, ERRMSG2);
	OutputBoth (thread_num, ERRMSG3);

/* Update the error count in the save file */

	writeNewErrorCount (filename, error_count);

/* Sleep five minutes before restarting */

	if (sleep5) {
		stop_reason = SleepFive (thread_num);
		if (stop_reason) return (stop_reason);
	}

/* Return so that last continuation file is read in */

	gwdone (&gwdata);
	free (N);
	goto begin;
}
