/* TODO: test larger values of mul-by-const */
/*	all_impl should work on small ffts by using the next large fft size */
/*		to do the comparison (and/or x87) */

/* Generate a random k of the given size */

double gen_k (int bits)
{
	double	k;
	int	i;

/* To reduce number of generic reductions due to gcd (k,c) > 1, only */
/* return small k values not divisible by 3, 5, and 7 */

	do {
		k = 1.0;
		for (i = bits; --i; )
			k = k * 2.0 + ((i == 1) ? 1 : (rand () & 1));
	} while (k > 7.0 && k < 4000000000.0 &&
		 ((unsigned long) k % 3 == 0 ||
		  (unsigned long) k % 5 == 0 ||
		  (unsigned long) k % 7 == 0));
	return (k);
}

/* Generate a random c of the given size */

int gen_c (int bits)
{
	int	c, i;

/* To reduce number of generic reductions due to gcd (k,c) > 1, only */
/* return c values not divisible by 3, 5, and 7 */

	do {
		c = 1;
		for (i = bits; --i; )
			c = c * 2 + ((i == 1) ? 1 : (rand () & 1));
	} while (c > 7 && (c % 3 == 0 || c % 5 == 0 || c % 7 == 0));
	if (rand () & 1) c = -c;
	return (c);
}

/* Generate a random n of the given size */

unsigned int gen_n ()
{
	unsigned int MIN_N, MAX_N;
	unsigned int min_n_bits, max_n_bits, range_n_bits;
	unsigned int n, nbits;

/* Get control variables */

	MIN_N = IniSectionGetInt (INI_FILE, "QA", "MIN_N", 250);
	MAX_N = IniSectionGetInt (INI_FILE, "QA", "MAX_N", 70000000);

/* Compute number of bits to generate */
/* Note: 1/log(2) = 1.442695040888963407 */
	
	min_n_bits = (unsigned int)
			floor (log ((double) MIN_N) * 1.442695040888963407);
	max_n_bits = (unsigned int)
			ceil (log ((double) MAX_N) * 1.442695040888963407);
	range_n_bits = max_n_bits - min_n_bits + 1;

/* Loop until we successfully create an n value in the desired range */

	do {
		nbits = (unsigned int) rand () % range_n_bits + min_n_bits;
 		n = 1;
		while (--nbits) n = n * 2 + (rand () & 1);
	} while (n < MIN_N || n > MAX_N);

/* Return our random n value */

	return (n);
}


/* Compare gwnum to a giant */

void compare (int thread_num, gwhandle *gwdata, gwnum x, giant g)
{
	giant	tmp;

	tmp = popg (&gwdata->gdata, (gwdata->n >> 4) + 13);
	gwtogiant (gwdata, x, tmp);
	if (gcompg (g, tmp) != 0)
		OutputBoth (thread_num, "Test failed.\n");
	pushg (&gwdata->gdata, 1);
}

/* Generate random data to start the test */

void gen_data (gwhandle *gwdata, gwnum x, giant g)
{
	unsigned long i, len;
	int	seed;
	char	buf[100];

/* Set and output seed so that we can re-generate the random data */

	seed = IniSectionGetInt (INI_FILE, "QA", "SPECIFIC_SEED",
				 (int) time (NULL) + rand ());
	srand (seed);
	sprintf (buf, "Random seed is %d\n", seed);
	writeResults (buf);

/* Generate the random number */

	len = (((unsigned long) gwdata->bit_length) >> 5) + 1;
	for (i = 0; i < len; i++) {
		g->n[i] = ((unsigned long) rand() << 20) +
			  ((unsigned long) rand() << 10) +
			  (unsigned long) rand();
	}
	g->sign = len;
	specialmodg (gwdata, g);
	gianttogw (gwdata, g, x);
}

/* Set and print random seed */

void set_seed (
	int	thread_num)
{
	int	seed;
//	char	buf[100];

	seed = (int) time (NULL);
	srand (seed + thread_num);
//	sprintf (buf, "Random seed is %ld\n", seed);
//	OutputBoth (thread_num, buf);
}

/* Thoroughly test the current setup.  This is just like test_it except */
/* that rather than using giants code to test the results, we compare */
/* the final result to every possible FFT implementation for this FFT */
/* length.  This is useful in that more FFT code is tested, the tests are */
/* faster -- especially in the really large FFT cases where giants is not */
/* practical. */

void test_it_all (
	int	thread_num,		/* Worker thread number */
	struct PriorityInfo *sp_info,	/* SetPriority information */
	double	k,
	unsigned long b,
	unsigned long n,
	signed long c,
	int	threads)
{
	gwhandle gwdata;
	gwnum	x, x2, x3, x4;
	giant	g, g2, g3;
	int	i, ii, res, nth_fft;
	double	diff, maxdiff;
	char	buf[200], fft_desc[100];

/* Init */

	g = g2 = g3 = NULL;

/* Loop over both x87 and SSE2 implementations.  Pass 1 does x87 FFTs */
/* on SSE2 machines.  Pass 2 does the "natural" FFTs. */

	for (ii = 1; ii <= 2; ii++) {

	    if (ii == 1) {
#ifdef X86_64
		continue;
#else
		if (! (CPU_FLAGS & CPU_SSE2)) continue;
#endif
	    }

/* Loop over all possible FFT implementations */

	    nth_fft = 1;
	    for ( ; ; ) {
		gwinit (&gwdata);
		gwset_num_threads (&gwdata, threads);
		gwset_thread_callback (&gwdata, SetAuxThreadPriority);
		gwset_thread_callback_data (&gwdata, sp_info);
		if (ii == 1) gwdata.cpu_flags = CPU_FLAGS & ~CPU_SSE2;
		else gwdata.cpu_flags = CPU_FLAGS;
		gwdata.qa_pick_nth_fft = nth_fft;
		res = gwsetup (&gwdata, k, b, n, c);
		nth_fft = gwdata.qa_pick_nth_fft;
		if (res) break;

/* Output a startup message */

		gwfft_description (&gwdata, fft_desc);
		sprintf (buf, "QA of %s using %s\n",
			 gwmodulo_as_string (&gwdata), fft_desc);
		OutputBoth (thread_num, buf);

/* Alloc and init numbers */

		x = gwalloc (&gwdata);
		if (x == NULL) goto nomem;
		x2 = gwalloc (&gwdata);
		if (x2 == NULL) goto nomem;
		x3 = gwalloc (&gwdata);
		if (x3 == NULL) goto nomem;
		x4 = gwalloc (&gwdata);
		if (x4 == NULL) goto nomem;

		if (g == NULL) {
			g = newgiant (((unsigned long) gwdata.bit_length >> 4) + 20);
			if (g == NULL) goto nomem;
			gen_data (&gwdata, x, g);
		} else
			gianttogw (&gwdata, g, x);

/* Test 50 squarings */	

		gwcopy (&gwdata, x, x2);
		maxdiff = 0.0;
		gwsetnormroutine (&gwdata, 0, 1, 0); /* Enable error checking */
		for (i = 0; i < 50; i++) {
			/* Test POSTFFT sometimes */
			gwstartnextfft (&gwdata, (i & 3) == 2);

			/* Test gwsetaddin without and with POSTFFT set */
			if ((i == 45 || i == 46) && abs (gwdata.c) == 1)
				gwsetaddin (&gwdata, -31);

			/* Test several different ways to square a number */
			if (i >= 4 && i <= 7) {
				gwfft (&gwdata, x, x);
				gwfftfftmul (&gwdata, x, x, x);
			} else if (i >= 12 && i <= 15) {
				gwfft (&gwdata, x, x3);
				gwfftmul (&gwdata, x3, x);
			} else if (i >= 20 && i <= 23) {
				gwfft (&gwdata, x, x3);
				gwcopy (&gwdata, x3, x4);
				gwfftfftmul (&gwdata, x3, x4, x);
			} else
				gwsquare (&gwdata, x);

			/* Remember maximum difference */
			diff = fabs (gwsuminp (&gwdata, x) - gwsumout (&gwdata, x));
			if (diff > maxdiff) maxdiff = diff;
			if ((i == 45 || i == 46) && abs (gwdata.c) == 1)
				gwsetaddin (&gwdata, 0);
		}
		if (gwdata.MAXDIFF < 1e50)
			sprintf (buf, "Squares complete. MaxErr=%.8g, SumoutDiff=%.8g/%.8g(%d to 1)\n", gw_get_maxerr (&gwdata), maxdiff, gwdata.MAXDIFF, (int) (gwdata.MAXDIFF / maxdiff));
		else
			sprintf (buf, "Squares complete. MaxErr=%.10g\n", gw_get_maxerr (&gwdata));
		OutputBoth (thread_num, buf);

/* Test mul by const */

		gwsetmulbyconst (&gwdata, 3);
		gwsetnormroutine (&gwdata, 0, 1, 1);
		gwsquare (&gwdata, x);
		gwsetnormroutine (&gwdata, 0, 1, 0);
		diff = fabs (gwsuminp (&gwdata, x) - gwsumout (&gwdata, x));
		if (diff > maxdiff) maxdiff = diff;

/* Test square carefully */

		gwfree (&gwdata, x3); gwfree (&gwdata, x4);
		if (abs (gwdata.c) == 1) gwsetaddin (&gwdata, -42);
		gwsquare_carefully (&gwdata, x);
		gwfree (&gwdata, gwdata.GW_RANDOM); gwdata.GW_RANDOM = NULL;
		diff = fabs (gwsuminp (&gwdata, x) - gwsumout (&gwdata, x));
		if (diff > maxdiff) maxdiff = diff;
		if (abs (gwdata.c) == 1) gwsetaddin (&gwdata, 0);

/* Test gwaddquick, gwsubquick */

		x3 = gwalloc (&gwdata); if (x3 == NULL) goto nomem;
		x4 = gwalloc (&gwdata); if (x4 == NULL) goto nomem;
		gwadd3quick (&gwdata, x, x2, x3);
		gwsub3quick (&gwdata, x, x2, x4);

/* Test gwadd and gwsub */

		gwadd (&gwdata, x, x); gwadd (&gwdata, x, x); gwadd (&gwdata, x, x);
		gwsub (&gwdata, x3, x);
		gwadd (&gwdata, x4, x);
		gwadd3 (&gwdata, x3, x4, x2);
		gwsub3 (&gwdata, x3, x, x4);
		gwadd (&gwdata, x2, x);
		gwadd (&gwdata, x4, x);

/* Test gwaddsub */

		gwaddsub (&gwdata, x, x2);	// compute x+x2 and x-x2
		gwaddsub4 (&gwdata, x, x2, x3, x4); // compute x+x2 and x-x2
		gwadd (&gwdata, x2, x);
		gwadd (&gwdata, x3, x);
		gwadd (&gwdata, x4, x);

/* Test gwaddsmall */

		gwaddsmall (&gwdata, x, 111);

/* Do some multiplies to make sure that the adds and subtracts above */
/* normalized properly. */

		gwfft (&gwdata, x, x);
		gwfftfftmul (&gwdata, x, x, x);
		diff = fabs (gwsuminp (&gwdata, x) - gwsumout (&gwdata, x));
		if (diff > maxdiff) maxdiff = diff;

		gwfft (&gwdata, x, x2); gwcopy (&gwdata, x2, x); gwfftadd3 (&gwdata, x, x2, x4);
		gwfftmul (&gwdata, x4, x3);
		diff = fabs (gwsuminp (&gwdata, x3) - gwsumout (&gwdata, x3));
		if (diff > maxdiff) maxdiff = diff;
		gwfft (&gwdata, x3, x4);
		gwfftfftmul (&gwdata, x4, x2, x);
		diff = fabs (gwsuminp (&gwdata, x) - gwsumout (&gwdata, x));
		if (diff > maxdiff) maxdiff = diff;

/* Print final stats */

		if (gwdata.GWERROR) OutputBoth (thread_num, "GWERROR set during calculations.\n");
		if (maxdiff > gwdata.MAXDIFF) OutputBoth (thread_num, "Sumout failed during test.\n");
		if (gwdata.MAXDIFF < 1e50)
			sprintf (buf, "Test complete. MaxErr=%.8g, SumoutDiff=%.8g/%.8g(%d to 1)\n", gw_get_maxerr (&gwdata), maxdiff, gwdata.MAXDIFF, (int) (gwdata.MAXDIFF / maxdiff));
		else
			sprintf (buf, "Test complete. MaxErr=%.10g\n", gw_get_maxerr (&gwdata));
		OutputBoth (thread_num, buf);

/* Free some space (so that gwtogiant can use it for temporaries) */

		gwfree (&gwdata, x2);
		gwfree (&gwdata, x3);
		gwfree (&gwdata, x4);

/* Do the final compare */

		if (g2 == NULL) {
			g2 = newgiant (((unsigned long) gwdata.bit_length >> 4) + 20);
			if (g2 == NULL) goto nomem;
			gwtogiant (&gwdata, x, g2);
		} else {
			g3 = newgiant (((unsigned long) gwdata.bit_length >> 4) + 20);
			if (g3 == NULL) goto nomem;
			gwtogiant (&gwdata, x, g3);
			if (gcompg (g2, g3)) {
				strcpy (buf, "Mismatched result.\n");
				OutputBoth (thread_num, buf);
			} else {
				strcpy (buf, "Results match!\n");
				OutputBoth (thread_num, buf);
			}
			free (g3);
		}
		OutputBoth (thread_num, "\n");

/* Do next FFT implementation */

		gwdone (&gwdata);
	    }
	}
bye:	free (g);
	free (g2);
	return;

nomem:	OutputBoth (thread_num, "Out of memory\n");
	gwdone (&gwdata);
	goto bye;
}


/* Thoroughly test the current setup */

void test_it (
	int	thread_num,		/* Worker thread number */
	gwhandle *gwdata)
{
	gwnum	x, x2, x3, x4;
	giant	g, g2, g3, g4;
	int	i;
	double	diff, maxdiff = 0.0;
	char	buf[200];
	int	SQUARE_ONLY, CHECK_OFTEN;

/* Get control variables */

	SQUARE_ONLY = IniSectionGetInt (INI_FILE, "QA", "SQUARE_ONLY", 0);
	CHECK_OFTEN = IniSectionGetInt (INI_FILE, "QA", "CHECK_OFTEN", 0);

/* Alloc and init numbers */

	x = gwalloc (gwdata);
	x2 = gwalloc (gwdata);
	x3 = gwalloc (gwdata);
	x4 = gwalloc (gwdata);
	g = popg (&gwdata->gdata, ((unsigned long) gwdata->bit_length >> 4) + 20);
	g2 = popg (&gwdata->gdata, ((unsigned long) gwdata->bit_length >> 4) + 20);
	g3 = popg (&gwdata->gdata, ((unsigned long) gwdata->bit_length >> 4) + 20);
	g4 = popg (&gwdata->gdata, ((unsigned long) gwdata->bit_length >> 4) + 20);
	gen_data (gwdata, x, g);
	if (CHECK_OFTEN) compare (thread_num, gwdata, x, g);
	gwcopy (gwdata, x, x2); gtog (g, g2);

/* Test 50 squarings */	

	gwsetnormroutine (gwdata, 0, 1, 0);	/* Enable error checking */
	for (i = 0; i < 50; i++) {

		/* Test POSTFFT sometimes */
		gwstartnextfft (gwdata, (i & 3) == 2);

		/* Test gwsetaddin without and with POSTFFT set */
		if ((i == 45 || i == 46) && abs (gwdata->c) == 1)
			gwsetaddin (gwdata, -31);

		/* Test several different ways to square a number */
		if (i >= 4 && i <= 7) {
			gwfft (gwdata, x, x);
			gwfftfftmul (gwdata, x, x, x);
		} else if (i >= 12 && i <= 15) {
			gwfft (gwdata, x, x3);
			gwfftmul (gwdata, x3, x);
		} else if (i >= 20 && i <= 23) {
			gwfft (gwdata, x, x3);
			gwcopy (gwdata, x3, x4);
			gwfftfftmul (gwdata, x3, x4, x);
		} else
			gwsquare (gwdata, x);

		/* Remember maximum difference */
		diff = fabs (gwsuminp (gwdata, x) - gwsumout (gwdata, x));
		if (diff > maxdiff) maxdiff = diff;

		/* Square number (and do add-in) using giants code */
		squaregi (&gwdata->gdata, g);
		if ((i == 45 || i == 46) && abs (gwdata->c) == 1) {
			iaddg (-31, g);
			gwsetaddin (gwdata, 0);
		}
		specialmodg (gwdata, g);

		/* Compare results */
		if (CHECK_OFTEN && (i & 3) != 2) compare (thread_num, gwdata, x, g);
	}
	if (SQUARE_ONLY) goto done;

	/* Report interim results */
	if (gwdata->MAXDIFF < 1e50)
		sprintf (buf,
			 "Squares complete. MaxErr=%.8g, SumoutDiff=%.8g/%.8g(%d to 1)\n",
			 gw_get_maxerr (gwdata), maxdiff, gwdata->MAXDIFF,
			 (int) (gwdata->MAXDIFF / maxdiff));
	else
		sprintf (buf, "Squares complete. MaxErr=%.10g\n", gw_get_maxerr (gwdata));
	OutputBoth (thread_num, buf);

/* Test mul by const */

	gwsetmulbyconst (gwdata, 3);
	gwsetnormroutine (gwdata, 0, 1, 1);
	gwsquare (gwdata, x);
	gwsetnormroutine (gwdata, 0, 1, 0);
	diff = fabs (gwsuminp (gwdata, x) - gwsumout (gwdata, x));
	if (diff > maxdiff) maxdiff = diff;
	squaregi (&gwdata->gdata, g); imulg (3, g);
	specialmodg (gwdata, g);
	if (CHECK_OFTEN) compare (thread_num, gwdata, x, g);

/* Test square carefully */

	gwfree (gwdata, x3); gwfree (gwdata, x4);
	if (abs (gwdata->c) == 1) gwsetaddin (gwdata, -42);
	gwsquare_carefully (gwdata, x);
	gwfree (gwdata, gwdata->GW_RANDOM); gwdata->GW_RANDOM = NULL;
	diff = fabs (gwsuminp (gwdata, x) - gwsumout (gwdata, x));
	if (diff > maxdiff) maxdiff = diff;
	squaregi (&gwdata->gdata, g);
	if (abs (gwdata->c) == 1) { iaddg (-42, g); gwsetaddin (gwdata, 0); }
	specialmodg (gwdata, g);
	if (CHECK_OFTEN) compare (thread_num, gwdata, x, g);

/* Test gwaddquick, gwsubquick */

	x3 = gwalloc (gwdata);
	x4 = gwalloc (gwdata);
	gwadd3quick (gwdata, x, x2, x3); gtog (g, g3); addg (g2, g3);
	gwsub3quick (gwdata, x, x2, x4); gtog (g, g4); subg (g2, g4);
	if (CHECK_OFTEN) {
		specialmodg (gwdata, g3); compare (thread_num, gwdata, x3, g3);
		specialmodg (gwdata, g4); compare (thread_num, gwdata, x4, g4);
	}

/* Test gwadd and gwsub */

	gwadd (gwdata, x, x); gwadd (gwdata, x, x); gwadd (gwdata, x, x);
	imulg (8, g);
	gwsub (gwdata, x3, x); subg (g3, g);
	gwadd (gwdata, x4, x); addg (g4, g);
	gwadd3 (gwdata, x3, x4, x2); gtog (g3, g2); addg (g4, g2);
	gwsub3 (gwdata, x3, x, x4); gtog (g3, g4); subg (g, g4);
	if (CHECK_OFTEN) {
		specialmodg (gwdata, g); compare (thread_num, gwdata, x, g);
		specialmodg (gwdata, g2); compare (thread_num, gwdata, x2, g2);
		specialmodg (gwdata, g4); compare (thread_num, gwdata, x4, g4);
	}
	gwadd (gwdata, x2, x); addg (g2, g);
	gwadd (gwdata, x4, x); addg (g4, g);

/* Test gwaddsub */

	gwaddsub (gwdata, x, x2);	// compute x+x2 and x-x2
	subg (g, g2);		// x2-x
	addg (g, g);		// x+x
	addg (g2, g);		// x+x2
	negg (g2);		// x-x2
	gwaddsub4 (gwdata, x, x2, x3, x4);	// compute x+x2 and x-x2
	gtog (g, g3); addg (g2, g3);
	gtog (g, g4); subg (g2, g4);
	if (CHECK_OFTEN) {
		specialmodg (gwdata, g); compare (thread_num, gwdata, x, g);
		specialmodg (gwdata, g2); compare (thread_num, gwdata, x2, g2);
		specialmodg (gwdata, g3); compare (thread_num, gwdata, x3, g3);
		specialmodg (gwdata, g4); compare (thread_num, gwdata, x4, g4);
	}
	gwadd (gwdata, x2, x); addg (g2, g);
	gwadd (gwdata, x3, x); addg (g3, g);
	gwadd (gwdata, x4, x); addg (g4, g);

/* Test gwaddsmall */

	gwaddsmall (gwdata, x, 111);
	iaddg (111, g); specialmodg (gwdata, g);
	if (CHECK_OFTEN) compare (thread_num, gwdata, x, g);

/* Do some multiplies to make sure that the adds and subtracts above */
/* normalized properly. */

	gwfft (gwdata, x, x);
	gwfftfftmul (gwdata, x, x, x);
	diff = fabs (gwsuminp (gwdata, x) - gwsumout (gwdata, x));
	if (diff > maxdiff) maxdiff = diff;
	squaregi (&gwdata->gdata, g); specialmodg (gwdata, g);
	if (CHECK_OFTEN) compare (thread_num, gwdata, x, g);

	gwfft (gwdata, x, x2); gwcopy (gwdata, x2, x); gwfftadd3 (gwdata, x, x2, x4);
	gwfftmul (gwdata, x4, x3); addg (g3, g3); mulgi (&gwdata->gdata, g, g3); specialmodg (gwdata, g3);
	diff = fabs (gwsuminp (gwdata, x3) - gwsumout (gwdata, x3));
	if (diff > maxdiff) maxdiff = diff;
	if (CHECK_OFTEN) compare (thread_num, gwdata, x3, g3);
	gwfft (gwdata, x3, x4);
	gwfftfftmul (gwdata, x4, x2, x);
	diff = fabs (gwsuminp (gwdata, x) - gwsumout (gwdata, x));
	if (diff > maxdiff) maxdiff = diff;
	mulgi (&gwdata->gdata, g3, g); specialmodg (gwdata, g);
	if (CHECK_OFTEN) compare (thread_num, gwdata, x, g);

/* Do the final compare */

done:	if (!CHECK_OFTEN) compare (thread_num, gwdata, x, g);
	if (gwdata->GWERROR) OutputBoth (thread_num, "GWERROR set during calculations.\n");

/* Print final stats */

	if (maxdiff > gwdata->MAXDIFF) OutputBoth (thread_num, "Sumout failed during test.\n");
	if (gwdata->MAXDIFF < 1e50)
		sprintf (buf,
			 "Test complete. MaxErr=%.8g, SumoutDiff=%.8g/%.8g(%d to 1)\n",
			 gw_get_maxerr (gwdata), maxdiff, gwdata->MAXDIFF,
			 (int) (gwdata->MAXDIFF / maxdiff));
	else
		sprintf (buf, "Test complete. MaxErr=%.10g\n", gw_get_maxerr (gwdata));
	OutputBoth (thread_num, buf);
	OutputBoth (thread_num, "\n");

	pushg (&gwdata->gdata, 4);
}

/* Ramdomly generate k,b,n,c combinations to test.  Results are compared */
/* again giants code performing the same operations. */

int test_randomly (
	int	thread_num,		/* Worker thread number */
	struct PriorityInfo *sp_info)	/* SetPriority information */
{
	gwhandle gwdata;
	int	kbits, cbits, threads, stop_reason;
	double	k;
	unsigned int b;
	int	c;
	unsigned int n;
	char	buf[100], fft_desc[100], SPECIFIC_K[20], SPECIFIC_C[20];
	int	MAX_THREADS, MIN_K_BITS, MAX_K_BITS;
	int	MAX_C_BITS_FOR_SMALL_K, MAX_C_BITS_FOR_LARGE_K;
	int	SPECIFIC_B, SPECIFIC_N, SPECIFIC_L2_CACHE, SPECIFIC_THREADS;
	int	FFT_LIMIT_THRESHOLD;

/* Get control variables */

	MAX_THREADS = IniSectionGetInt (INI_FILE, "QA", "MAX_THREADS", NUM_CPUS * CPU_HYPERTHREADS);
	if (MAX_THREADS < 1) MAX_THREADS = 1;
	if (MAX_THREADS > (int) (NUM_CPUS * CPU_HYPERTHREADS))
		MAX_THREADS = NUM_CPUS * CPU_HYPERTHREADS;
	MIN_K_BITS = IniSectionGetInt (INI_FILE, "QA", "MIN_K_BITS", 1);
	MAX_K_BITS = IniSectionGetInt (INI_FILE, "QA", "MAX_K_BITS", 49);
	MAX_C_BITS_FOR_SMALL_K = IniSectionGetInt (INI_FILE, "QA", "MAX_C_BITS_FOR_SMALL_K", 4);
	MAX_C_BITS_FOR_LARGE_K = IniSectionGetInt (INI_FILE, "QA", "MAX_C_BITS_FOR_LARGE_K", 20);
	/* Only QA numbers that are within this percent of the FFT limit */
	/* The default is 100% of numbers are tested.  I set this to 1% to */
	/* test numbers near the FFT limits to look for excessive */
	/* round off problems. */
	FFT_LIMIT_THRESHOLD = IniSectionGetInt (INI_FILE, "QA", "FFT_LIMIT_THRESHOLD", 100);

	IniSectionGetString (INI_FILE, "QA", "SPECIFIC_K", SPECIFIC_K, 20, "1");
	SPECIFIC_B = IniSectionGetInt (INI_FILE, "QA", "SPECIFIC_B", 2);
	SPECIFIC_N = IniSectionGetInt (INI_FILE, "QA", "SPECIFIC_N", 0);
	IniSectionGetString (INI_FILE, "QA", "SPECIFIC_C", SPECIFIC_C, 20, "-1");
	SPECIFIC_L2_CACHE = IniSectionGetInt (INI_FILE, "QA", "SPECIFIC_L2_CACHE", CPU_L2_CACHE_SIZE);
	SPECIFIC_THREADS = IniSectionGetInt (INI_FILE, "QA", "SPECIFIC_THREADS", 1);

	set_seed (thread_num);
	b = 2;
	for ( ; ; ) {

/* Abort loop when requested */

		stop_reason = stopCheck (thread_num);
		if (stop_reason) break;

/* Generate the number to QA */

		if (MAX_K_BITS == MIN_K_BITS) kbits = 0;
		else kbits = (unsigned int) rand () % (MAX_K_BITS - MIN_K_BITS + 1);
		kbits += MIN_K_BITS;
		if (kbits <= 20)
			cbits = (unsigned int) rand () % MAX_C_BITS_FOR_SMALL_K + 1;
		else
			cbits = (unsigned int) rand () % MAX_C_BITS_FOR_LARGE_K + 1;
		k = gen_k (kbits);
		c = gen_c (cbits);
		n = gen_n ();
		threads = rand () % MAX_THREADS + 1;
		switch (rand () % 5) {
			case 0:	CPU_L2_CACHE_SIZE = 128; break;
			case 1:	CPU_L2_CACHE_SIZE = 256; break;
			case 2:	CPU_L2_CACHE_SIZE = 512; break;
			case 3:	CPU_L2_CACHE_SIZE = 1024; break;
			case 4:	CPU_L2_CACHE_SIZE = -1; break;
		}

/* Override the generated values and use the specified values.  This is */
/* primarily done for debugging. */

		if (SPECIFIC_N) {
			k = atof (SPECIFIC_K);
			b = SPECIFIC_B;
			n = SPECIFIC_N;
			c = atoi (SPECIFIC_C);
			CPU_L2_CACHE_SIZE = SPECIFIC_L2_CACHE;
			threads = SPECIFIC_THREADS;
			FFT_LIMIT_THRESHOLD = 100;
		}

/* Test using the default FFT implementation and comparing it to giants */

		gwinit (&gwdata);
		gwset_num_threads (&gwdata, threads);
		gwset_thread_callback (&gwdata, SetAuxThreadPriority);
		gwset_thread_callback_data (&gwdata, sp_info);
		gwsetup (&gwdata, k, b, n, c);
		if (! gwnear_fft_limit (&gwdata, FFT_LIMIT_THRESHOLD)) {
			gwdone (&gwdata);
			continue;
		}
//if (gwdata.ZERO_PADDED_FFT || gwdata.GENERAL_MOD) {
//	gwdone (&gwdata);
//	continue;
//}
		sprintf (buf, "Starting QA run on %s, kbits=%ld, cbits=%ld\n",
			 gwmodulo_as_string (&gwdata), kbits, cbits);
		OutputBoth (thread_num, buf);
		gwfft_description (&gwdata, fft_desc);
		if (gwfftlen (&gwdata) < 140000)
			sprintf (buf,
				 "Using %s, bits_per_word=%.5g/%.5g\n",
				 fft_desc,
				 virtual_bits_per_word (&gwdata),
				 gwdata.fft_max_bits_per_word);
		else
			sprintf (buf, "Using %s, bits_per_word=%.5g/%.5g, L2_cache_size=%d\n",
				 fft_desc,
				 virtual_bits_per_word (&gwdata),
				 gwdata.fft_max_bits_per_word,
				 CPU_L2_CACHE_SIZE);
		OutputBoth (thread_num, buf);
		test_it (thread_num, &gwdata);
		gwdone (&gwdata);
	}

/* All done */

	return (stop_reason);
}

/* Ramdomly generate k,b,n,c combinations to test.  Results are compared */
/* using all possible implementations of the FFT (different L2 cache sizes */
/* and non-SSE2). */

int test_all_impl (
	int	thread_num,		/* Worker thread number */
	struct PriorityInfo *sp_info)	/* SetPriority information */
{
	int	kbits, cbits, threads, stop_reason;
	double	k;
	unsigned int b;
	int	c;
	unsigned int n;
	int	MAX_THREADS, MIN_K_BITS, MAX_K_BITS;
	int	MAX_C_BITS_FOR_SMALL_K, MAX_C_BITS_FOR_LARGE_K;

/* Get control variables */

	MAX_THREADS = IniSectionGetInt (INI_FILE, "QA", "MAX_THREADS", NUM_CPUS * CPU_HYPERTHREADS);
	if (MAX_THREADS < 1) MAX_THREADS = 1;
	if (MAX_THREADS > (int) (NUM_CPUS * CPU_HYPERTHREADS))
		MAX_THREADS = NUM_CPUS * CPU_HYPERTHREADS;
	MIN_K_BITS = IniSectionGetInt (INI_FILE, "QA", "MIN_K_BITS", 1);
	MAX_K_BITS = IniSectionGetInt (INI_FILE, "QA", "MAX_K_BITS", 49);
	MAX_C_BITS_FOR_SMALL_K = IniSectionGetInt (INI_FILE, "QA", "MAX_C_BITS_FOR_SMALL_K", 4);
	MAX_C_BITS_FOR_LARGE_K = IniSectionGetInt (INI_FILE, "QA", "MAX_C_BITS_FOR_LARGE_K", 20);

	set_seed (thread_num);
	b = 2;
	for ( ; ; ) {

/* Abort loop when requested */

		stop_reason = stopCheck (thread_num);
		if (stop_reason) break;

/* Generate the number to QA */

		if (MAX_K_BITS == MIN_K_BITS) kbits = 0;
		else kbits = (unsigned int) rand () % (MAX_K_BITS - MIN_K_BITS + 1);
		kbits += MIN_K_BITS;
		if (kbits <= 20)
			cbits = (unsigned int) rand () % MAX_C_BITS_FOR_SMALL_K + 1;
		else
			cbits = (unsigned int) rand () % MAX_C_BITS_FOR_LARGE_K + 1;
		k = gen_k (kbits);
		c = gen_c (cbits);
		n = gen_n ();
		threads = rand () % MAX_THREADS + 1;

/* Test using multiple implementations of an FFT size */

		test_it_all (thread_num, sp_info, k, b, n, c, threads);
	}

/* All done */

	return (stop_reason);
}