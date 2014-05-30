/*----------------------------------------------------------------------
| gwtables.c
|
| This file contains the C routines to build sin/cos and weights tables
| that the FFT assembly code needs.
| 
|  Copyright 2002-2014 Mersenne Research, Inc.  All rights reserved.
+---------------------------------------------------------------------*/

/* Include files */

#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include "gwnum.h"
#include "gwtables.h"
#include "gwdbldbl.h"

#define USE_WPN4
#define USE_REDUCED_SINCOS_FFTS

/* Find the power of two greater than or equal to N. */

unsigned long pow_two_above_or_equal (
	unsigned long n)
{
	unsigned long result;

	result = 1;
	for (n = n - 1; n; n = n >> 1) result = result << 1;
	return (result);
}

/* This routine builds the sin/cos table used in a one pass AVX traditional */
/* DJB radix-4 FFT  - called by gwsetup.  If this is an all-complex FFT, */
/* then the root-of-minus-1 premultipliers are also built. */

double *yr4_build_onepass_sincos_table (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	double	*table)		/* Pointer to the table to fill in */
{
	unsigned long size, avx_increment, j, N, temp;
	int	pow2_count;

/* Count the power-of-two FFT levels after the initial FFT levels.  If odd, the */
/* the last levels will be a radix-8, if even the last levels will be a radix-4. */

	size = gwdata->FFTLEN / 8;
	for (pow2_count = 0; (size & 1) == 0; size /= 2) pow2_count++;

/* Init necessary variables */

	size = gwdata->FFTLEN / 8;		/* Complex values we're generating sin/cos data for */
	avx_increment = gwdata->FFTLEN / 4;

/* The first group has no sin/cos data! */

	if (pow2_count & 1) {
		N = 8;
		size /= 8;
	} else {
		N = 4;
		size /= 4;
	}

/* For the yr4_four_complex_djbfft and yr4_eight_reals_four_complex_djbfft */
/* building block levels, output the sin/cos values. */

	if ((size & 3) == 0) {
		while ((size & 3) == 0) {
			N = N * 4;
			size /= 4;
		}

/* For the yr4_eight_reals_four_complex_djbfft building block levels, output the */
/* sin/cos values needed.  The eight_reals doubles N because the real part of the FFT */
/* is one level behind the complex part of the FFT.  The four-complex sin/cos values */
/* are the same for all 3 of the upper YMM doubles. */		

		if (!gwdata->ALL_COMPLEX_FFT) {
			for (j = 0; j < N / 4; j++) {
				gwsincos125by4 (j, N*2, table);				/* For the eight_reals */
				gwsincos12by4 (j + avx_increment, N, table+1);		/* For the four-complex */
				table[3] = table[2] = table[1];
				table[7] = table[6] = table[5];
				table[11] = table[10] = table[9];
				table[15] = table[14] = table[13];
				table[19] = table[18] = table[17] = -table[1];
				table[23] = table[22] = table[21] = -table[5];
				table += 24;
			}
		}

/* Output the sin/cos values for the all complex FFTs, used by the yr4__b4cl_four_complex_djbfft macro. */
/* We only need one sin/cos value pair as we use the vbroadcastsd instruction to fill out the YMM register. */

		else {
			for (j = 0; j < N / 4; j++) {
				gwsincos12by1 (j, N, table);
				table += 4;
			}
		}
	}

/* For the yr5_five_complex_djbfft building block levels, output the sin/cos values. */

	while (size % 5 == 0) {
		while ((size % 5) == 0) {
			N = N * 5;
			size /= 5;
		}

/* For the yr5_ten_reals_five_complex_djbfft building block levels, output the */
/* sin/cos values needed.  The ten_reals doubles N because the real part of the FFT */
/* is one level behind the complex part of the FFT.  The five-complex sin/cos values */
/* are the same for all 3 of the upper YMM doubles. */

		if (!gwdata->ALL_COMPLEX_FFT) {
			for (j = 0; j < N / 5; j++) {
				gwsincos13by4 (j, N*2, table);				/* For the ten_reals */
				gwsincos12by4 (j + avx_increment, N, table+1);		/* For the five-complex */
				table[16] = table[8];					/* For the ten_reals */
				table[20] = table[12];
				table[8] = table[1];
				table[12] = table[5];
				table[24] = table[9];
				table[28] = table[13];
				table[3] = table[2] = table[1];				/* For the five-complex */
				table[7] = table[6] = table[5];
				table[11] = table[10] = table[9];
				table[15] = table[14] = table[13];
				table[19] = table[18] = table[17] = -table[9];
				table[23] = table[22] = table[21] = -table[13];
				table[27] = table[26] = table[25] = -table[1];
				table[31] = table[30] = table[29] = -table[5];
				table += 32;
			}
		}

/* Output the sin/cos data for the complex sections (used by the yr5_five_complex_djbfft building block). */

		else {
			for (j = 0; j < N / 5; j++) {
				gwsincos12by1 (j, N, table);
				table += 4;
			}
		}
	}

/* For the yr3_three_complex_djbfft building block levels, output the sin/cos values. */

	while (size % 3 == 0) {
		while ((size % 3) == 0) {
			N = N * 3;
			size /= 3;
		}

/* The yr3_six_reals building blocks require an extra sin/cos */
/* value.  The six_reals doubles N because the real part of the */
/* FFT is one level behind the complex part of the FFT. */

		if (!gwdata->ALL_COMPLEX_FFT) {
			for (j = 0; j < N / 3; j++) {
				gwsincos12by4 (j, N*2, table);				/* For the six-reals FFT */
				gwsincos1by4 (j + avx_increment, N, table+1);		/* For the three-complex FFT */
				table[3] = table[2] = table[1];
				table[7] = table[6] = table[5];
				table[11] = table[10] = table[9] = -table[1];
				table[15] = table[14] = table[13] = -table[5];
				table += 16;
			}
		}

/* Output the sin/cos data for the complex sections (used by the yr3_three_complex_djbfft building block). */

		else {
			for (j = 0; j < N / 3; j++) {
				gwsincos1by1 (j, N, table);
				table += 2;
			}
		}
	}
	ASSERTG (size == 1);
	if (size != 1) gwdata->GWERROR = GWERROR_INTERNAL + 1;

/* Real FFTs output one last set of sin/cos values for the first 8-reals FFT. */

	avx_increment = gwdata->FFTLEN / 32;
	if (! gwdata->ALL_COMPLEX_FFT) {
		N = gwdata->FFTLEN;
		for (j = 0; j < avx_increment; j++) {
			gwsincos125by4 (j, N, table);
			gwsincos125by4 (j + avx_increment, N, table+1);
			gwsincos125by4 (j + 2*avx_increment, N, table+2);
			gwsincos125by4 (j + 3*avx_increment, N, table+3);
			table += 24;
		}
	}

/* For all-complex FFTs, the first FFT level converts from real values to all complex */
/* values by multiplying by a root of -1 weight and doing a butterfly.  This is */
/* simplified because weights in the bottom half are sqrt(-1) times the */
/* matching weights in the upper half.  Thus, we butterfly upper_word * weight */
/* with bottom_word * i * weight.  That equals (upper_word + i * lower_word) * weight */
/* That is just a complex multiply with the half of the butterfly output values */
/* unneeded thanks to Hermetian symmetry. */

/* The second and third levels do a standard radix-4 four-complex-FFT */
/* building block with a post-multiply by a sin/cos root of unity. */

	else {
		N = gwdata->FFTLEN / 2;
		for (j = 0; j < avx_increment; j++) {

/* Here we compute the standard 0,1,2,3 * temp for the radix-4 sin/cos multipliers. */
/* Then we multiply in the roots-of-minus-one premultiplier.  The root-of-minus-one */
/* premultiplier was for 2N, and a root-of-minus-one-of-2N is the same as a root */
/* unity for 4N. */

			temp = j;
			gwsincos1plus0123by4 (temp, 4 * temp, N * 4, table); /* premult + temp*0-3 */
			temp += avx_increment;
			gwsincos1plus0123by4 (temp, 4 * temp, N * 4, table+1);
			temp += avx_increment;
			gwsincos1plus0123by4 (temp, 4 * temp, N * 4, table+2);
			temp += avx_increment;
			gwsincos1plus0123by4 (temp, 4 * temp, N * 4, table+3);
			table += 32;
		}
	}

/* Return address of the end of the table */

	return (table);
}

/* This routine builds a big/little flags table - used by one-pass AVX normalization */
/* routines */

double *yr4_build_onepass_biglit_table (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	double	*table)		/* Pointer to the table to fill in */
{
	struct gwasm_data *asm_data = (struct gwasm_data *) gwdata->asm_data;
	unsigned char *p = (unsigned char *) table;
	unsigned long i, j, top5bits;

/* Init table of first 8 big/lit values */

	memset (asm_data->u.ymm.YMM_FIRST_BIGLIT_VALUES, 0, sizeof (asm_data->u.ymm.YMM_FIRST_BIGLIT_VALUES));

/* Loop to build table for zero-padded case.  We only need half as much data */
/* because the upper half data are the same as the lower half data. */

	if (gwdata->ZERO_PADDED_FFT) {
		memset (p, 0, gwdata->FFTLEN / 8);
		for (i = 0; i < gwdata->FFTLEN / 2; i++) {

/* Only big words result in a bit being set in the biglit table */

			if (! is_big_word (gwdata, i)) continue;

/* Find where this data appears in the table we are building.  Use the same algorithm */
/* as addr_offset except no padding is necessary. */			

			top5bits = i / (gwdata->FFTLEN >> 5); j = i - top5bits * (gwdata->FFTLEN >> 5);
			j = ((top5bits >> 2) & 3) * (gwdata->FFTLEN >> 2) + (j << 3) + ((top5bits >> 4) << 2) + (top5bits & 3);

/* Add to the biglit table entry for each big word double in an AVX word */

			p[j >> 3] += 16 << (j & 3);

/* Create separate table for first 8 biglit values for carry propagation code */

			if (i < sizeof (asm_data->u.ymm.YMM_FIRST_BIGLIT_VALUES))
				asm_data->u.ymm.YMM_FIRST_BIGLIT_VALUES[i] = 16;
		}

/* Return pointer past the end of the biglit table */

		return ((double *) (p + gwdata->FFTLEN / 8));
	}

/* Loop to build table for the non-zero-padded case */

	memset (p, 0, gwdata->FFTLEN / 4);
	for (i = 0; i < gwdata->FFTLEN; i++) {

/* Only big words result in a bit being set in the biglit table */

		if (! is_big_word (gwdata, i)) continue;

/* Find where this data appears in the table we are building.  Use the same algorithm */
/* as addr_offset except no padding is necessary. */			

		top5bits = i / (gwdata->FFTLEN >> 5); j = i - top5bits * (gwdata->FFTLEN >> 5);
		j = ((top5bits >> 2) & 3) * (gwdata->FFTLEN >> 2) + (j << 3) + ((top5bits >> 4) << 2) + (top5bits & 3);

/* Add to the biglit table entry for each big word double in an AVX word */

		p[j >> 2] += 16 << (j & 3);

/* Create separate table for first 8 biglit values for carry propagation code */

		if (i < sizeof (asm_data->u.ymm.YMM_FIRST_BIGLIT_VALUES))
			asm_data->u.ymm.YMM_FIRST_BIGLIT_VALUES[i] = 16;
	}

/* Return pointer past the end of the biglit table */

	return ((double *) (p + gwdata->FFTLEN / 4));
}

/* This routine builds a normalization table - used by one-pass AVX radix-4 */
/* normalization routines. */

double *yr4_build_onepass_norm_table (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	double	*table)		/* Pointer to the table to fill in */
{
	struct gwasm_data *asm_data = (struct gwasm_data *) gwdata->asm_data;
	unsigned long i;
	intptr_t first_8_offsets[8];

/* Loop to build table for zero-padded case.  We only need half as much data */
/* because the upper half data are the same as the lower half data. */

	if (gwdata->ZERO_PADDED_FFT) {
		for (i = 0; i < gwdata->FFTLEN / 2; i++) {
			double	ttp, ttmp;
			unsigned long top5bits, j, table_entry;

/* Call double-precision routine to compute the two multipliers */

			gwfft_weights3 (gwdata->dd_data, i, &ttp, NULL, &ttmp);

/* Find where this data appears in the table we are building.  Use the same algorithm */
/* as addr_offset except no padding is necessary. */			

			top5bits = i / (gwdata->FFTLEN >> 5); j = i - top5bits * (gwdata->FFTLEN >> 5);
			j = ((top5bits >> 2) & 3) * (gwdata->FFTLEN >> 2) + (j << 3) + ((top5bits >> 4) << 2) + (top5bits & 3);
			table_entry = j >> 3;

/* Now set the entry for the proper double in an AVX word */

			table[table_entry*8+(j&3)] = ttmp;
			table[table_entry*8+4+(j&3)] = ttp;

/* Get offsets for carry propagation code to step through norm array */

			if (i <= 7) first_8_offsets[i] = (&table[table_entry*8+(j&3)] - table) * sizeof (double);
		}

/* Form pointer for next table */

		table += gwdata->FFTLEN;
	}

/* Loop to build table for non-zero-padded case */

	else {
		for (i = 0; i < gwdata->FFTLEN; i++) {
			double	ttp, ttmp;
			unsigned long top5bits, j, table_entry;

/* Call double-precision routine to compute the two multipliers */

			gwfft_weights3 (gwdata->dd_data, i, &ttp, NULL, &ttmp);

/* Find where this data appears in the table we are building.  Use the same algorithm */
/* as addr_offset except no padding is necessary. */			

			top5bits = i / (gwdata->FFTLEN >> 5); j = i - top5bits * (gwdata->FFTLEN >> 5);
			j = ((top5bits >> 2) & 3) * (gwdata->FFTLEN >> 2) + (j << 3) + ((top5bits >> 4) << 2) + (top5bits & 3);
			table_entry = j >> 2;

/* Now set the entry for the proper double in an AVX word */

			table[table_entry*8+(j&3)] = ttmp;
			table[table_entry*8+4+(j&3)] = ttp;

/* Get offsets for carry propagation code to step through norm array */

			if (i <= 7) first_8_offsets[i] = (&table[table_entry*8+(j&3)] - table) * sizeof (double);
		}

/* Create pointer for next table */

		table += gwdata->FFTLEN * 2;
	}

/* Make differences for carry propagation code to step through norm array */

	asm_data->u.ymm.YMM_NORM_INCR7 = first_8_offsets[7] - first_8_offsets[6];
	asm_data->u.ymm.YMM_NORM_INCR6 = first_8_offsets[6] - first_8_offsets[5];
	asm_data->u.ymm.YMM_NORM_INCR5 = first_8_offsets[5] - first_8_offsets[4];
	asm_data->u.ymm.YMM_NORM_INCR4 = first_8_offsets[4] - first_8_offsets[3];
	asm_data->u.ymm.YMM_NORM_INCR3 = first_8_offsets[3] - first_8_offsets[2];
	asm_data->u.ymm.YMM_NORM_INCR2 = first_8_offsets[2] - first_8_offsets[1];
	asm_data->u.ymm.YMM_NORM_INCR1 = first_8_offsets[1] - first_8_offsets[0];

/* Return pointer for next table */

	return (table);
}

/* This routine builds the sin/cos table used in pass 1 by the radix-4/8 DJB */
/* FFT with delayed sin/cos multiplies and with partial normalization. */

double *yr4dwpn_build_pass1_table (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	double	*table)		/* Pointer to the table to fill in */
{
	struct gwasm_data *asm_data = (struct gwasm_data *) gwdata->asm_data;
	unsigned long pass1_size, pass1_increment, delay_count;
	unsigned long group, i, j, k, N, temp, upper_avx_word;
	int	pow2_count;
	int	wpn4 = FALSE;		/* Flag indicating we are using wpn4 in pass 1 */
	int	rsc = FALSE;		/* Flag indicating we are using reduced sin/cos in pass 1 */

/* Initialize some needed constants */

	pass1_size = gwdata->FFTLEN / gwdata->PASS2_SIZE; /* Real values in a pass1 section */
	upper_avx_word = gwdata->PASS2_SIZE;
	pass1_increment = gwdata->PASS2_SIZE * 4;

/* Determine number of delay groups.  In a standard radix-4 FFT, there is only one sin/cos */
/* group in the last pass 1 level.  We reduce our memory usage by using a fixed sin/cos */
/* table in the first FFT levels and having multiple groups of sin/cos data in the last pass 1 level. */
/* I call these groups of sin/cos data in the last pass 1 level "delay groups". */

#ifdef USE_REDUCED_SINCOS_FFTS
	if (pass1_size % 7 == 0)
		delay_count = 14;
	else if (pass1_size == 384 || pass1_size == 768)
		delay_count = 12;
	else if ((pass1_size == 640 && gwdata->ALL_COMPLEX_FFT) ||
		 (pass1_size == 1280 && gwdata->ALL_COMPLEX_FFT))
		delay_count = 20;
	else if (pass1_size == 1280 && !gwdata->ALL_COMPLEX_FFT)
		delay_count = 40;
	else if (pass1_size % 5 == 0 && !gwdata->ALL_COMPLEX_FFT)
		delay_count = 10;
	else if ((pass1_size == 256 && !gwdata->ALL_COMPLEX_FFT) ||
		 (pass1_size == 512 && !gwdata->ALL_COMPLEX_FFT))
		delay_count = 8;
	else if ((pass1_size == 512 && gwdata->ALL_COMPLEX_FFT) ||
		 pass1_size == 1024)
		delay_count = 16;
	else
		delay_count = 4;
	if (pass1_size == 128 || pass1_size == 256 || pass1_size == 320 || pass1_size == 384 ||
	    pass1_size == 448 || pass1_size == 512 || pass1_size == 640 || pass1_size == 768 ||
	    pass1_size == 896 || pass1_size == 1024 || pass1_size == 1280) {
		rsc = TRUE;				// Someday we can convert the pass 1 sizes above 1280 
	}
#endif
	if (!rsc) {
		if (pass1_size % 7 == 0)
			delay_count = 14;
		else if ((pass1_size == 384 && !gwdata->ALL_COMPLEX_FFT) ||
			 (pass1_size == 768 && !gwdata->ALL_COMPLEX_FFT) ||
			 (pass1_size == 1536 && !gwdata->ALL_COMPLEX_FFT))
			delay_count = 12;
		else if (pass1_size % 5 == 0 && !gwdata->ALL_COMPLEX_FFT)
			delay_count = 10;
		else if (pass1_size == 256 && !gwdata->ALL_COMPLEX_FFT)
			delay_count = 8;
		else if ((pass1_size == 512 && !gwdata->ALL_COMPLEX_FFT) ||
			 (pass1_size == 1024 && !gwdata->ALL_COMPLEX_FFT) ||
			 (pass1_size == 1536 && gwdata->ALL_COMPLEX_FFT) ||
			 pass1_size == 2048)
			delay_count = 16;
		else
			delay_count = 4;
	}

/* Count the power-of-two FFT levels after the initial FFT levels.  If odd, the */
/* the last levels will be a radix-8, if even the last levels will be a radix-4. */

	pass1_size /= (delay_count * 2);
	for (pow2_count = 0; (pass1_size & 1) == 0; pass1_size /= 2) pow2_count++;

/* Set count of pass 1 blocks that share one set of two-to-phi grp multipliers */

	if (pow2_count & 1) gwdata->wpn_count = 8;
	else if ((gwdata->FFTLEN / gwdata->PASS2_SIZE == 1536 && !gwdata->ALL_COMPLEX_FFT) ||
		 gwdata->FFTLEN / gwdata->PASS2_SIZE == 1792 ||
		 gwdata->FFTLEN / gwdata->PASS2_SIZE == 2048) gwdata->wpn_count = 16;
	else gwdata->wpn_count = 4;
#ifdef USE_WPN4
	{
		gwdata->wpn_count *= 4;
		wpn4 = TRUE;
	}
#endif

/* Set counters for inorm, zpnorm and ygw_carries to use.  Remember that ygw_carries */
/* always works on data after it has been copied to the scratch area. */

	asm_data->count2 = gwdata->wpn_count / 4;
	asm_data->count3 = asm_data->addcount1 / asm_data->count2;
	if (asm_data->count2 == 1) {
		asm_data->count4 = 1;
		asm_data->count5 = asm_data->count3 / 2;
	} else {
		asm_data->count4 = asm_data->count2 / 2;
		asm_data->count5 = asm_data->count3;
	}

/* Set pointer to table of multipliers */

	gwdata->adjusted_pass1_premults = table;

/* Loop through all the pass 1 groups in the same order the assembly code will */
/* process the groups. */

	for (group = 0; group < upper_avx_word; group += gwdata->PASS1_CACHE_LINES) {

		pass1_size = gwdata->FFTLEN / gwdata->PASS2_SIZE; /* Real values in a pass1 section */
		pass1_size /= (delay_count * 2);	/* Complex values we're generating sin/cos data for */
		N = gwdata->PASS2_SIZE;

/* Output the sin/cos/premultiplier values for the radix-8 block that does the */
/* last 3 levels in pass 1.  NOTE:  We do not need the "j loop" (it would loop */
/* from zero to zero) when generating the sin/cos twiddle factors for the last */
/* levels of pass 1. */

		if (rsc && (pow2_count & 1)) {
			N = N * 8;

/* Output the complex sin/cos values needed for a standard yr8_8cl_eight_complex_djbfft */
/* on the last pass 1 level.  At runtime, we compute the actual sin/cos values from this. */

			for (i = 0; i < gwdata->PASS1_CACHE_LINES; i += 4) {
				// Asm code swizzled the input so that upper_avx_word is 1
				temp = group + i;
				gwsincos1234by4_raw (temp, N, table);
				gwsincos1234by4_raw (temp + 1, N, table+1);
				gwsincos1234by4_raw (temp + 2, N, table+2);
				gwsincos1234by4_raw (temp + 3, N, table+3);
				table += 32;
			}

/* For the yr8_sg8cl_sixteen_reals_fft8 building block, output the extra */
/* sin/cos values needed for the sixteen_reals. */

			if (!gwdata->ALL_COMPLEX_FFT) {
				for (i = 0; i < gwdata->PASS1_CACHE_LINES; i += 4) {
					// Asm code swizzled the input so that upper_avx_word is 1
					temp = group + i;
					gwsincos15913by4 (temp, N*2, table);
					gwsincos15913by4 (temp + 1, N*2, table+1);
					gwsincos15913by4 (temp + 2, N*2, table+2);
					gwsincos15913by4 (temp + 3, N*2, table+3);
					table += 32;
				}
			}

/* Output the sin/cos values for the complex delay groups -- specifically the yr8_rsc_sg8cl_eight_complex_fft8 macro. */

			for (k = 0; k < delay_count; k++) {
				if (k == 0 && !gwdata->ALL_COMPLEX_FFT) continue;
				for (i = 0; i < gwdata->PASS1_CACHE_LINES; i += 4) {
					unsigned long bigN, ktemp, actemp, avx_word;

/* Work on each AVX word.  Unlike the SSE2 build-table code, we must recalculate */
/* ktemp for each AVX word because the ASM code swizzles its inputs */

					for (avx_word = 0; avx_word < 4; avx_word++) {
						unsigned long final_group = group + i + avx_word;

/* If this is an all-complex FFT, the roots of minus 1 (same as roots of FFTLEN*2) are */
/* split to reduce memory requirements.  We apply the part of the all-complex premultiplier here. */

						if (gwdata->ALL_COMPLEX_FFT) {
							bigN = gwdata->FFTLEN * 2;
							actemp = final_group;
						} else {
							bigN = gwdata->FFTLEN;
							actemp = 0;
						}

/* Factor in the delayed part of the sin/cos multiplies from the first 2 levels.  In the first 2 levels */
/* we use a fixed sin/cos table based only on j, leaving the group+i part to be applied here by */
/* creating delay_count table entries.  For an all-complex FFT we multiply by 0,2,1,-1*(group+i) with N = FFTLEN/2. */
/* For an all-real FFT we multiply by 0,2,1,5*(group+i) with N = FFTLEN. */

						if (gwdata->ALL_COMPLEX_FFT && delay_count == 4) {
							if (k == 0)
								ktemp = 0;
							else if (k == 1)
								ktemp = 2 * final_group * 4;
							else if (k == 2)
								ktemp = 1 * final_group * 4;
							else
								ktemp = bigN - 1 * final_group * 4;
						} else if (gwdata->ALL_COMPLEX_FFT && delay_count == 12) {
							/* 0,2,1,-1 combined with 0,1,-1 */
							int	kmap[12] = {0,4,-4, 2,6,-2, 1,5,-3, -1,3,-5};
							if (kmap[k] >= 0)
								ktemp = kmap[k] * final_group * 4;
							else
								ktemp = bigN + kmap[k] * final_group * 4;
						} else if (gwdata->ALL_COMPLEX_FFT && delay_count == 20) {
							/* 0,2,1,-1 combined with 0,1,2,-2,-1 */
							int	kmap[20] = {0,4,8,-8,-4, 2,6,10,-6,-2, 1,5,9,-7,-3, -1,3,7,-9,-5};
							if (kmap[k] >= 0)
								ktemp = kmap[k] * final_group * 4;
							else
								ktemp = bigN + kmap[k] * final_group * 4;
						} else if (gwdata->ALL_COMPLEX_FFT) {
							/* 0,2,1,-1 combined with 0,2,1,-1 */
							int	kmap[16] = {0,8,4,-4, 2,10,6,-2, 1,9,5,-3, -1,7,3,-5};
							if (kmap[k] >= 0)
								ktemp = kmap[k] * final_group * 4;
							else
								ktemp = bigN + kmap[k] * final_group * 4;
						} else if (delay_count == 4) {
							if (k == 0)
								ktemp = 0;
							else if (k == 1)
								ktemp = 2 * final_group;
							else if (k == 2)
								ktemp = 1 * final_group;
							else
								ktemp = 5 * final_group;
						} else if (delay_count == 16) {
							/* 0,2,1,5 combined with one 0,2,1,5 and three 0,2,1,-1 */
							int	kmap[16] = {0,8,4,20, 2,18,10,-6, 1,17,9,-7, 5,21,13,-3};
							if (kmap[k] >= 0)
								ktemp = kmap[k] * final_group;
							else
								ktemp = bigN + kmap[k] * final_group;
						} else if (delay_count == 40) {
							/* 0...9 combined with one 0,2,1,5 and three 0,2,1,-1 */
							int	kmap[40] = {0,20,10,50,   1,41,21,-19,  2,42,22,-18,
									    3,43,23,-17,  4,44,24,-16,  5,45,25,-15,
									    6,46,26,-14,  7,47,27,-13,  8,48,28,-12,  9,49,29,-11};
							if (kmap[k] >= 0)
								ktemp = kmap[k] * final_group;
							else
								ktemp = bigN + kmap[k] * final_group;
						} else {			/* delay_count == 8, 10, 12, or 14 */
							/* Multipliers for the radix-16, radix-20, radix-24, or radix-28 step */
							ktemp = k * final_group;
						}

/* We now calculate the group all-complex premultiplier roots of minus 1 (same as roots of FFTLEN*2) */
/* combined with the delayed group multipliers. */

						gwsincos1by4_raw (actemp + ktemp, bigN, table + avx_word);
					}
					table += 8;
				}
			}
			pass1_size /= 8;
		}

/* Output the sin/cos/premultiplier values for the radix-4 block that does the */
/* last 2 levels in pass 1. */

		else if (rsc) {
			N = N * 4;

/* Output the complex sin/cos values needed for a standard yr4_4cl_four_complex_djbfft */
/* on the last pass 1 level.  At runtime, we compute the actual sin/cos values from this. */

			for (i = 0; i < gwdata->PASS1_CACHE_LINES; i += 4) {
				// Asm code swizzled the input so that upper_avx_word is 1
				temp = group + i;
				gwsincos12by4_raw (temp, N, table);
				gwsincos12by4_raw (temp + 1, N, table+1);
				gwsincos12by4_raw (temp + 2, N, table+2);
				gwsincos12by4_raw (temp + 3, N, table+3);
				table += 16;
			}

/* Output the extra sin/cos values needed for the eight_reals FFT work done */
/* on the last pass 1 level.  We double N because the real part of the FFT */
/* is one level behind the complex part of the FFT. */

			if (!gwdata->ALL_COMPLEX_FFT) {
				for (i = 0; i < gwdata->PASS1_CACHE_LINES; i += 4) {
					// Asm code swizzled the input so that upper_avx_word is 1
					temp = group + i;
					gwsincos15by4 (temp, N*2, table);
					gwsincos15by4 (temp + 1, N*2, table+1);
					gwsincos15by4 (temp + 2, N*2, table+2);
					gwsincos15by4 (temp + 3, N*2, table+3);
					table += 16;
				}
			}

/* Output the sin/cos values for the delay groups -- specifically the yr4_rsc_sg4cl_four_complex_fft4 macro. */

			for (k = 0; k < delay_count; k++) {
				if (k == 0 && !gwdata->ALL_COMPLEX_FFT) continue;
				for (i = 0; i < gwdata->PASS1_CACHE_LINES; i += 4) {
					unsigned long bigN, ktemp, actemp, avx_word;

/* Work on each AVX word.  Unlike the SSE2 build-table code, we must recalculate */
/* ktemp for each AVX word because the ASM code swizzles its inputs */

					for (avx_word = 0; avx_word < 4; avx_word++) {
						unsigned long final_group = group + i + avx_word;

/* If this is an all-complex FFT, the roots of minus 1 (same as roots of FFTLEN*2) are */
/* split to reduce memory requirements.  We apply part of the all-complex premultiplier here. */

						if (gwdata->ALL_COMPLEX_FFT) {
							bigN = gwdata->FFTLEN * 2;
							actemp = final_group;
						} else {
							bigN = gwdata->FFTLEN;
							actemp = 0;
						}

/* Factor in the delayed part of the sin/cos multiplies from the first 2 levels.  In the first 2 levels */
/* we use a fixed sin/cos table based only on j, leaving the group+i part to be applied here by */
/* creating delay_count table entries.  For an all-complex FFT we multiply by 0,2,1,-1*(group+i) with N = FFTLEN/2. */
/* For an all-real FFT we multiply by 0,2,1,5*(group+i) with N = FFTLEN. */

						if (gwdata->ALL_COMPLEX_FFT && delay_count == 4) {
							if (k == 0)
								ktemp = 0;
							else if	(k == 1)
								ktemp = 2 * final_group * 4;
							else if (k == 2)
								ktemp = 1 * final_group * 4;
							else
								ktemp = bigN - 1 * final_group * 4;
						} else if (gwdata->ALL_COMPLEX_FFT && delay_count == 12) {
							/* 0,2,1,-1 combined with 0,1,-1 */
							int	kmap[12] = {0,4,-4, 2,6,-2, 1,5,-3, -1,3,-5};
							if (kmap[k] >= 0)
								ktemp = kmap[k] * final_group * 4;
							else
								ktemp = bigN + kmap[k] * final_group * 4;
						} else if (gwdata->ALL_COMPLEX_FFT && delay_count == 20) {
							/* 0,2,1,-1 combined with 0,1,2,-2,-1 */
							int	kmap[20] = {0,4,8,-8,-4, 2,6,10,-6,-2, 1,5,9,-7,-3, -1,3,7,-9,-5};
							if (kmap[k] >= 0)
								ktemp = kmap[k] * final_group * 4;
							else
								ktemp = bigN + kmap[k] * final_group * 4;
						} else if (gwdata->ALL_COMPLEX_FFT) {
							/* 0,2,1,-1 combined with 0,2,1,-1 */
							int	kmap[16] = {0,8,4,-4, 2,10,6,-2, 1,9,5,-3, -1,7,3,-5};
							if (kmap[k] >= 0)
								ktemp = kmap[k] * final_group * 4;
							else
								ktemp = bigN + kmap[k] * final_group * 4;
						} else if (delay_count == 4) {
							if (k == 0)
								ktemp = 0;
							else if	(k == 1)
								ktemp = 2 * final_group;
							else if (k == 2)
								ktemp = 1 * final_group;
							else
								ktemp = 5 * final_group;
						} else if (delay_count == 16) {
							/* 0,2,1,5 combined with one 0,2,1,5 and three 0,2,1,-1 */
							int	kmap[16] = {0,8,4,20, 2,18,10,-6, 1,17,9,-7, 5,21,13,-3};
							if (kmap[k] >= 0)
								ktemp = kmap[k] * final_group;
							else
								ktemp = bigN + kmap[k] * final_group;
						} else if (delay_count == 40) {
							/* 0...9 combined with one 0,2,1,5 and three 0,2,1,-1 */
							int	kmap[40] = {0,20,10,50,   1,41,21,-19,  2,42,22,-18,
									    3,43,23,-17,  4,44,24,-16,  5,45,25,-15,
									    6,46,26,-14,  7,47,27,-13,  8,48,28,-12,  9,49,29,-11};
							if (kmap[k] >= 0)
								ktemp = kmap[k] * final_group;
							else
								ktemp = bigN + kmap[k] * final_group;
						} else {			/* delay_count == 8, 10, 12, or 14 */
							/* Multipliers for the radix-16, radix-20, radix-24, or radix-28 step */
							ktemp = k * final_group;
						}

/* We now calculate the group all-complex premultiplier roots of minus 1 (same as roots of FFTLEN*2) */
/* combined with the delayed group multipliers. */

						gwsincos1by4_raw (actemp + ktemp, bigN, table + avx_word);
					}
					table += 8;
				}
			}
			pass1_size /= 4;
		}

/* Output the sin/cos/premultiplier values for the radix-8 block that does the */
/* last 3 levels in pass 1.  NOTE:  We do not need the "j loop" (it would loop */
/* from zero to zero) when generating the sin/cos twiddle factors for the last */
/* levels of pass 1. */

		else if (pow2_count & 1) {
			N = N * 8;

/* For the yr8_sg8cl_sixteen_reals_fft8 building block, output the extra */
/* sin/cos values needed for the sixteen_reals. */

			if (!gwdata->ALL_COMPLEX_FFT) {
				for (i = 0; i < gwdata->PASS1_CACHE_LINES; i += 4) {
					// Asm code swizzled the input so that upper_avx_word is 1
					temp = group + i;
					gwsincos15913by4 (temp, N*2, table);
					gwsincos15913by4 (temp + 1, N*2, table+1);
					gwsincos15913by4 (temp + 2, N*2, table+2);
					gwsincos15913by4 (temp + 3, N*2, table+3);
					table += 32;
				}
			}

/* Output the sin/cos values for the complex delay groups -- specifically the yr8_sg8cl_eight_complex_fft8 macro. */

			for (k = 0; k < delay_count; k++) {
				for (i = 0; i < gwdata->PASS1_CACHE_LINES; i += 4) {
					unsigned long bigN, ktemp, actemp, avx_word;

/* Work on each AVX word.  Unlike the SSE2 build-table code, we must recalculate */
/* ktemp for each AVX word because the ASM code swizzles its inputs */

					for (avx_word = 0; avx_word < 4; avx_word++) {
						unsigned long final_group = group + i + avx_word;

/* If this is an all-complex FFT, the roots of minus 1 (same as roots of FFTLEN*2) are */
/* split to reduce memory requirements.  We apply the part of the all-complex premultiplier here. */

						if (gwdata->ALL_COMPLEX_FFT) {
							bigN = gwdata->FFTLEN * 2;
							actemp = final_group;
						} else {
							bigN = gwdata->FFTLEN;
							actemp = 0;
						}

/* Factor in the delayed part of the sin/cos multiplies from the first 2 levels.  In the first 2 levels */
/* we use a fixed sin/cos table based only on j, leaving the group+i part to be applied here by */
/* creating delay_count table entries.  For an all-complex FFT we multiply by 0,2,1,-1*(group+i) with N = FFTLEN/2. */
/* For an all-real FFT we multiply by 0,2,1,5*(group+i) with N = FFTLEN. */

						if (gwdata->ALL_COMPLEX_FFT && delay_count == 4) {
							if (k == 0)
								ktemp = 0;
							else if (k == 1)
								ktemp = 2 * final_group * 4;
							else if (k == 2)
								ktemp = 1 * final_group * 4;
							else
								ktemp = bigN - 1 * final_group * 4;
						} else if (gwdata->ALL_COMPLEX_FFT) {
							/* 0,2,1,-1 combined with 0,2,1,-1 */
							int	kmap[16] = {0,8,4,-4,2,10,6,-2,1,9,5,-3,-1,7,3,-5};
							if (kmap[k] >= 0)
								ktemp = kmap[k] * final_group * 4;
							else
								ktemp = bigN + kmap[k] * final_group * 4;
						} else if (delay_count == 4) {
							if (k == 0)
								ktemp = 0;
							else if (k == 1)
								ktemp = 2 * final_group;
							else if (k == 2)
								ktemp = 1 * final_group;
							else
								ktemp = 5 * final_group;
						} else if (delay_count == 16) {
							/* 0,2,1,5 combined with one 0,2,1,5 and three 0,2,1,-1 */
							int	kmap[16] = {0,8,4,20,2,18,10,-6,1,17,9,-7,5,21,13,-3};
							if (kmap[k] >= 0)
								ktemp = kmap[k] * final_group;
							else
								ktemp = bigN + kmap[k] * final_group;
						} else {			/* delay_count == 8, 10, 12, or 14 */
							/* Multipliers for the radix-16, radix-20, radix-24, or radix-28 step */
							ktemp = k * final_group;
						}

/* We now calculate the standard sin/cos twiddle factors (temp*0,1,2,3,4,5,6,7 roots of N) */
/* combined with the group all-complex premultiplier roots of minus 1 (same as roots of FFTLEN*2) */
/* combined with the delayed group multipliers. */

						temp = final_group * (bigN / N);
						gwsincos1plus01234567by4 (actemp + ktemp, temp, bigN, table + avx_word); /* premult,delay and temp*0-7 */
					}
					table += 64;
				}
			}
			pass1_size /= 8;
		}

/* Output the sin/cos/premultiplier values for the radix-4 block that does the */
/* last 2 levels in pass 1. */

		else {
			N = N * 4;

/* Output the extra sin/cos values needed for the eight_reals FFT work done */
/* on the last pass 1 level.  We double N because the real part of the FFT */
/* is one level behind the complex part of the FFT. */

			if (!gwdata->ALL_COMPLEX_FFT) {
				for (i = 0; i < gwdata->PASS1_CACHE_LINES; i += 4) {
					// Asm code swizzled the input so that upper_avx_word is 1
					temp = group + i;
					gwsincos15by4 (temp, N*2, table);
					gwsincos15by4 (temp + 1, N*2, table+1);
					gwsincos15by4 (temp + 2, N*2, table+2);
					gwsincos15by4 (temp + 3, N*2, table+3);
					table += 16;
				}
			}

/* Output the sin/cos values for the complex delay groups -- specifically the yr4_sg4cl_four_complex_fft4 macro. */

			for (k = 0; k < delay_count; k++) {
				for (i = 0; i < gwdata->PASS1_CACHE_LINES; i += 4) {
					unsigned long bigN, ktemp, actemp, avx_word;

/* Work on each AVX word.  Unlike the SSE2 build-table code, we must recalculate */
/* ktemp for each AVX word because the ASM code swizzles its inputs */

					for (avx_word = 0; avx_word < 4; avx_word++) {
						unsigned long final_group = group + i + avx_word;

/* If this is an all-complex FFT, the roots of minus 1 (same as roots of FFTLEN*2) are */
/* split to reduce memory requirements.  We apply part of the all-complex premultiplier here. */

						if (gwdata->ALL_COMPLEX_FFT) {
							bigN = gwdata->FFTLEN * 2;
							actemp = final_group;
						} else {
							bigN = gwdata->FFTLEN;
							actemp = 0;
						}

/* Factor in the delayed part of the sin/cos multiplies from the first 2 levels.  In the first 2 levels */
/* we use a fixed sin/cos table based only on j, leaving the group+i part to be applied here by */
/* creating delay_count table entries.  For an all-complex FFT we multiply by 0,2,1,-1*(group+i) with N = FFTLEN/2. */
/* For an all-real FFT we multiply by 0,2,1,5*(group+i) with N = FFTLEN. */

						if (gwdata->ALL_COMPLEX_FFT && delay_count == 4) {
							if (k == 0)
								ktemp = 0;
							else if	(k == 1)
								ktemp = 2 * final_group * 4;
							else if (k == 2)
								ktemp = 1 * final_group * 4;
							else
								ktemp = bigN - 1 * final_group * 4;
						} else if (gwdata->ALL_COMPLEX_FFT) {
							/* 0,2,1,-1 combined with 0,2,1,-1 */
							int	kmap[16] = {0,8,4,-4,2,10,6,-2,1,9,5,-3,-1,7,3,-5};
							if (kmap[k] >= 0)
								ktemp = kmap[k] * final_group * 4;
							else
								ktemp = bigN + kmap[k] * final_group * 4;
						} else if (delay_count == 4) {
							if (k == 0)
								ktemp = 0;
							else if	(k == 1)
								ktemp = 2 * final_group;
							else if (k == 2)
								ktemp = 1 * final_group;
							else
								ktemp = 5 * final_group;
						} else if (delay_count == 16) {
							/* 0,2,1,5 combined with one 0,2,1,5 and three 0,2,1,-1 */
							int	kmap[16] = {0,8,4,20,2,18,10,-6,1,17,9,-7,5,21,13,-3};
							if (kmap[k] >= 0)
								ktemp = kmap[k] * final_group;
							else
								ktemp = bigN + kmap[k] * final_group;
						} else {			/* delay_count == 8, 10, 12, or 14 */
							/* Multipliers for the radix-16, radix-20, radix-24, or radix-28 step */
							ktemp = k * final_group;
						}

/* We now calculate the standard sin/cos twiddle factors (temp*0,1,2,3 roots of N) */
/* combined with the group all-complex premultiplier roots of minus 1 (same as roots of FFTLEN*2) */
/* combined with the delayed group multipliers. */

						temp = final_group * (bigN / N);
						gwsincos1plus0123by4 (actemp + ktemp, temp, bigN, table + avx_word); /* premult,delay and temp*0-3 */
					}
					table += 32;
				}
			}
			pass1_size /= 4;
		}

/* Output multipliers for the four complex building blocks. */

		while ((pass1_size & 3) == 0) {

			N = N * 4;

/* For the wpn4 building block level, output a separate table of column normalization values before the sin/cos data. */

			if (wpn4 && N == gwdata->PASS2_SIZE * gwdata->wpn_count) {
				double *weights, *inv_weights;

/* The weights are output in separate tables before the sin/cos values.  This requires two registers */
/* to access the tables, but gains in that we can group data in cache lines better. */

				weights = table;
				table += N / pass1_increment * gwdata->PASS1_CACHE_LINES;
				inv_weights = table;
				table += N / pass1_increment * gwdata->PASS1_CACHE_LINES;

/* Output the weights before the sin/cos data, used by the yr4_4cl_wpn4_four_complex_djbfft macro. */
/* We apply the two-to-phi weight for the upper AVX words in the group multipliers.  There is a */
/* reason for doing it there rather than here (it reduces the number of valid fudge factor combinations */
/* for each AVX word from 16 to 5). */

				for (j = 0; j < N / 4; j += pass1_increment) {
				    for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
					double	not_used;
					temp = (group + j + i);
					gwfft_weights3 (gwdata->dd_data, temp, weights, &not_used, inv_weights);
					gwfft_weights3 (gwdata->dd_data, temp + N/4, weights+1, &not_used, inv_weights+1);
					gwfft_weights3 (gwdata->dd_data, temp + 2*N/4, weights+2, &not_used, inv_weights+2);
					gwfft_weights3 (gwdata->dd_data, temp + 3*N/4, weights+3, &not_used, inv_weights+3);
					weights += 4;
					inv_weights += 4;
				    }
				}
			}

/* For the non-wpn and wpn4 levels, output the sin/cos values. */

			if (wpn4 || N != gwdata->PASS2_SIZE * gwdata->wpn_count * 4) {

/* Output the sin/cos value for the complex sections, used by the yr4_4cl_four_complex_djbfft macro */

				for (j = 0; j < N / 4; j += pass1_increment) {
				    for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
					temp = (group + j + i);
					gwsincos12by4 (temp, N, table);
					gwsincos12by4 (temp + upper_avx_word, N, table+1);
					gwsincos12by4 (temp + 2 * upper_avx_word, N, table+2);
					gwsincos12by4 (temp + 3 * upper_avx_word, N, table+3);
					table += 16;

/* For the yr4_4cl_csc_eight_reals_fft building block levels, output the extra */
/* sin/cos values needed for the eight_reals.  The eight_reals doubles N because */
/* the real part of the FFT is one level behind the complex part of the FFT. */

					if (!gwdata->ALL_COMPLEX_FFT) {
						gwsincos15by4 (temp, N*2, table);
						gwsincos15by4 (temp + upper_avx_word, N*2, table+1);
						gwsincos15by4 (temp + 2 * upper_avx_word, N*2, table+2);
						gwsincos15by4 (temp + 3 * upper_avx_word, N*2, table+3);
						table += 16;
					}

				    }
				}
			}

/* For the wpn building block level, output the sin/cos and column normalization values. */

			if (!wpn4 && N == gwdata->PASS2_SIZE * gwdata->wpn_count * 4) {
				double *weights, *inv_weights;

/* The weights are output in separate tables before the sin/cos values.  This requires two registers */
/* to access the tables, but gains in that we can group data in cache lines better. */

				weights = table;
				table += N / 4 / pass1_increment * gwdata->PASS1_CACHE_LINES;
				inv_weights = table;
				table += N / 4 / pass1_increment * gwdata->PASS1_CACHE_LINES;

/* Output the sin/cos value for the complex sections, used by the yr4_4cl_wpn_four_complex_djbfft macro */
/* We apply the two-to-phi weight for the upper AVX words in the group multipliers.  There is a */
/* reason for doing it there rather than here (it reduces the number of valid fudge factor combinations */
/* for each AVX word from 16 to 5). */

				for (j = 0; j < N / 4; j += pass1_increment) {
				    for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
					temp = (group + j + i);
					gwsincos012by4_weighted (gwdata->dd_data, temp, upper_avx_word, N, temp, table);
					*weights++ = table[24];
					*inv_weights++ = table[25];
					table += 24;

/* For the yr4_4cl_csc_wpn_eight_reals_four_complex_djbfft building block levels, output the extra */
/* sin/cos values needed for the eight_reals.  The eight_reals doubles N because */
/* the real part of the FFT is one level behind the complex part of the FFT. */

					if (!gwdata->ALL_COMPLEX_FFT) {
						gwsincos15by4_weighted (gwdata->dd_data, temp, upper_avx_word, N*2, temp, table);
						table += 24;
					}
				    }
				}
			}

			pass1_size /= 4;
		}

/* For the yr5_five_complex_djbfft building block levels, output the sin/cos values. */

		while (pass1_size % 5 == 0) {
			N = N * 5;

/* Output the sin/cos data for the complex sections, (the yr5_5cl_five_complex_djbfft building block). */

			for (j = 0; j < N / 5; j += pass1_increment) {
				for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
					temp = (group + j + i);
					gwsincos12by4 (temp, N, table);
					gwsincos12by4 (temp + upper_avx_word, N, table+1);
					gwsincos12by4 (temp + 2 * upper_avx_word, N, table+2);
					gwsincos12by4 (temp + 3 * upper_avx_word, N, table+3);
					table += 16;

/* The yr5_5cl_csc_ten_reals building blocks require extra sin/cos values.  The ten_reals doubles N */
/* because the real part of the FFT is one level behind the complex part of the FFT. */

					if (!gwdata->ALL_COMPLEX_FFT) {
						gwsincos13by4 (temp, N*2, table);
						gwsincos13by4 (temp + upper_avx_word, N*2, table+1);
						gwsincos13by4 (temp + 2 * upper_avx_word, N*2, table+2);
						gwsincos13by4 (temp + 3 * upper_avx_word, N*2, table+3);
						table += 16;
					}
				}
			}
			pass1_size /= 5;
		}

/* For the yr3_3cl_three_complex_djbfft building block levels, output the sin/cos values. */

		while (pass1_size % 3 == 0) {
			N = N * 3;

/* Output the sin/cos data for the complex sections (used by the yr3_3cl_three_complex_djbfft building block). */

			for (j = 0; j < N / 3; j += pass1_increment) {
				for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
					temp = (group + j + i);
					gwsincos1by4 (temp, N, table);
					gwsincos1by4 (temp + upper_avx_word, N, table+1);
					gwsincos1by4 (temp + 2 * upper_avx_word, N, table+2);
					gwsincos1by4 (temp + 3 * upper_avx_word, N, table+3);
					table += 8;

/* The yr3_3cl_csc_six_reals building blocks require an extra sin/cos value.  The six_reals doubles N */
/* because the real part of the FFT is one level behind the complex part of the FFT. */

					if (!gwdata->ALL_COMPLEX_FFT) {
						gwsincos1by4 (temp, N*2, table);
						gwsincos1by4 (temp + upper_avx_word, N*2, table+1);
						gwsincos1by4 (temp + 2 * upper_avx_word, N*2, table+2);
						gwsincos1by4 (temp + 3 * upper_avx_word, N*2, table+3);
						table += 8;
					}
				}
			}
			pass1_size /= 3;
		}

/* Calculate the size of each group's sin/cos/premult data for pass1_get_next_block */

		if (group == 0)
			gwdata->pass1_premult_block_size = (unsigned long)
				((char *) table - (char *) gwdata->adjusted_pass1_premults) /
				gwdata->PASS1_CACHE_LINES;
	}

/* Return address of the end of the table */

	return (table);
}

/* This routine builds the fixed postmultiplier table used in pass 1 of the */
/* radix-4/8 delayed DJB FFT - called by gwsetup. */

double *yr4dwpn_build_fixed_pass1_table (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	double	*table)		/* Pointer to the table to fill in */
{
	unsigned long upper_avx_word, pass1_increment, pass1_size, i, j, N;

/* Initialize some needed constants */

	upper_avx_word = gwdata->PASS2_SIZE;
	pass1_increment = 4 * gwdata->PASS2_SIZE;
	pass1_size = gwdata->FFTLEN / gwdata->PASS2_SIZE; /* Real values in a pass1 section */

/* Real FFTs output one shared set of sin/cos values for the first 16-reals, 20-reals, 24-reals, 28-reals, or 8-reals FFT. */

	if (! gwdata->ALL_COMPLEX_FFT) {
		N = gwdata->FFTLEN;
		if (pass1_size == 256 || pass1_size == 512) {
			for (j = 0; j < N / 16; j += pass1_increment) {
				for (i = 1; i <= 7; i++) {	/* Create 7 twiddle factors */
					gwsincos1by4 (i * j, N, table);
					gwsincos1by4 (i * (j + upper_avx_word), N, table+1);
					gwsincos1by4 (i * (j + 2 * upper_avx_word), N, table+2);
					gwsincos1by4 (i * (j + 3 * upper_avx_word), N, table+3);
					table += 8;
				}
			}
			N = N / 16;
		}
		else if (pass1_size % 20 == 0) {
			for (j = 0; j < N / 20; j += pass1_increment) {
				for (i = 1; i <= 9; i++) {	/* Create 9 twiddle factors */
					gwsincos1by4 (i * j, N, table);
					gwsincos1by4 (i * (j + upper_avx_word), N, table+1);
					gwsincos1by4 (i * (j + 2 * upper_avx_word), N, table+2);
					gwsincos1by4 (i * (j + 3 * upper_avx_word), N, table+3);
					table += 8;
				}
			}
			N = N / 20;
			/* Sometimes we also use a fixed sin/cos table for */
			/* the next FFT levels to further reduce memory usage. */
#ifdef USE_REDUCED_SINCOS_FFTS
			if (pass1_size == 1280) {
				/* Output the sin/cos values for the complex data followed by sin/cos values for the real data */
				for (j = 0; j < N / 4; j += pass1_increment) {
					gwsincos12by4 (j, N, table);
					gwsincos12by4 (j + upper_avx_word, N, table+1);
					gwsincos12by4 (j + 2 * upper_avx_word, N, table+2);
					gwsincos12by4 (j + 3 * upper_avx_word, N, table+3);
					table += 16;
					gwsincos15by4 (j, N*2, table);
					gwsincos15by4 (j + upper_avx_word, N*2, table+1);
					gwsincos15by4 (j + 2 * upper_avx_word, N*2, table+2);
					gwsincos15by4 (j + 3 * upper_avx_word, N*2, table+3);
					table += 16;
				}
				N = N / 4;
			}
#endif
		}
		else if (pass1_size == 384 || pass1_size == 768 || pass1_size == 1536) {
			for (j = 0; j < N / 24; j += pass1_increment) {
				for (i = 1; i <= 11; i++) {	/* Create 11 twiddle factors */
					gwsincos1by4 (i * j, N, table);
					gwsincos1by4 (i * (j + upper_avx_word), N, table+1);
					gwsincos1by4 (i * (j + 2 * upper_avx_word), N, table+2);
					gwsincos1by4 (i * (j + 3 * upper_avx_word), N, table+3);
					table += 8;
				}
			}
			N = N / 24;
		}
		else if (pass1_size % 28 == 0) {
			for (j = 0; j < N / 28; j += pass1_increment) {
				for (i = 1; i <= 13; i++) {	/* Create 13 twiddle factors */
					gwsincos1by4 (i * j, N, table);
					gwsincos1by4 (i * (j + upper_avx_word), N, table+1);
					gwsincos1by4 (i * (j + 2 * upper_avx_word), N, table+2);
					gwsincos1by4 (i * (j + 3 * upper_avx_word), N, table+3);
					table += 8;
				}
			}
			N = N / 28;
		}
		else {
			for (j = 0; j < N / 8; j += pass1_increment) {
				gwsincos125by4 (j, N, table);
				gwsincos125by4 (j + upper_avx_word, N, table+1);
				gwsincos125by4 (j + 2 * upper_avx_word, N, table+2);
				gwsincos125by4 (j + 3 * upper_avx_word, N, table+3);
				table += 24;
			}
			N = N / 8;
			/* Sometimes we also use a fixed sin/cos table for */
			/* the next FFT levels to further reduce memory usage. */
			if (pass1_size == 1024 || pass1_size == 2048) {
				/* Output the sin/cos values for the complex data followed by sin/cos values for the real data */
				for (j = 0; j < N / 4; j += pass1_increment) {
					gwsincos12by4 (j, N, table);
					gwsincos12by4 (j + upper_avx_word, N, table+1);
					gwsincos12by4 (j + 2 * upper_avx_word, N, table+2);
					gwsincos12by4 (j + 3 * upper_avx_word, N, table+3);
					table += 16;
					gwsincos15by4 (j, N*2, table);
					gwsincos15by4 (j + upper_avx_word, N*2, table+1);
					gwsincos15by4 (j + 2 * upper_avx_word, N*2, table+2);
					gwsincos15by4 (j + 3 * upper_avx_word, N*2, table+3);
					table += 16;
				}
				N = N / 4;
			}
		}
	}

/* For all-complex FFTs, build the fixed roots-of-minus-one table and the */
/* DJB FFT sin/cos table.  Output these values in the same order they will */
/* be used in the first two levels of pass 1. */

	else {
		N = gwdata->FFTLEN / 2;
		for (j = 0; j < N / 4; j += pass1_increment) {
			/* Compute the roots-of-minus-one premultiplier.  The root-of-minus-one */
			/* premultiplier is for 2N, and a root-of-minus-one-of-2N is the same as */
			/* a root unity for 4N. */					
			gwsincos1plus0123by4 (j, N / 4, N * 4, table);
			gwsincos1plus0123by4 (j + upper_avx_word, N / 4, N * 4, table + 1);
			gwsincos1plus0123by4 (j + 2 * upper_avx_word, N / 4, N * 4, table + 2);
			gwsincos1plus0123by4 (j + 3 * upper_avx_word, N / 4, N * 4, table + 3);
			table += 32;
			/* Output the fixed sin/cos DJB FFT entry */
			gwsincos12by4 (j, N, table);
			gwsincos12by4 (j + upper_avx_word, N, table+1);
			gwsincos12by4 (j + 2 * upper_avx_word, N, table+2);
			gwsincos12by4 (j + 3 * upper_avx_word, N, table+3);
			table += 16;
		}
		N = N / 4;
		/* Sometimes we also use a fixed sin/cos table for */
		/* the next FFT levels to further reduce memory usage. */
#ifdef USE_REDUCED_SINCOS_FFTS
		if (pass1_size == 384 || pass1_size == 768) {
			for (j = 0; j < N / 3; j += pass1_increment) {
				gwsincos1by4 (j, N, table);
				gwsincos1by4 (j + upper_avx_word, N, table+1);
				gwsincos1by4 (j + 2 * upper_avx_word, N, table+2);
				gwsincos1by4 (j + 3 * upper_avx_word, N, table+3);
				table += 8;
			}
			N = N / 3;
		}
		if (pass1_size == 512 || pass1_size == 1024) {
			/* Output the sin/cos values for the complex data */
			for (j = 0; j < N / 4; j += pass1_increment) {
				gwsincos12by4 (j, N, table);
				gwsincos12by4 (j + upper_avx_word, N, table+1);
				gwsincos12by4 (j + 2 * upper_avx_word, N, table+2);
				gwsincos12by4 (j + 3 * upper_avx_word, N, table+3);
				table += 16;
			}
			N = N / 4;
		}
		if (pass1_size == 640 || pass1_size == 1280) {
			for (j = 0; j < N / 5; j += pass1_increment) {
				gwsincos12by4 (j, N, table);
				gwsincos12by4 (j + upper_avx_word, N, table+1);
				gwsincos12by4 (j + 2 * upper_avx_word, N, table+2);
				gwsincos12by4 (j + 3 * upper_avx_word, N, table+3);
				table += 16;
			}
			N = N / 5;
		}
#else
		if (pass1_size == 1536 || pass1_size == 2048) {
			/* Output the sin/cos values for the complex data */
			for (j = 0; j < N / 4; j += pass1_increment) {
				gwsincos12by4 (j, N, table);
				gwsincos12by4 (j + upper_avx_word, N, table+1);
				gwsincos12by4 (j + 2 * upper_avx_word, N, table+2);
				gwsincos12by4 (j + 3 * upper_avx_word, N, table+3);
				table += 16;
			}
			N = N / 4;
		}
#endif
	}

/* Return address of the end of the table */

	return (table);
}

/* This routine builds the sin/cos table used in all-complex pass 2 */
/* blocks in a traditional radix-4 FFT - called by gwsetup. */

double *yr4_build_pass2_complex_table (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	double	*table)		/* Pointer to the table to fill in */
{
	unsigned long i, j, N;
	unsigned int pow2_count;

/* If the pass 2 size is divisible by 3, then the initial levels do */
/* radix-3 steps which only requires one sin/cos value.  The first levels */
/* in pass 2 have an upper_avx_word of one. */

	N = gwdata->PASS2_SIZE;
	while (N % 3 == 0) {
		for (i = 0; i < N / 3; i += 4) {
			gwsincos1by4 (i, N, table);
			gwsincos1by4 (i+1, N, table+1);
			gwsincos1by4 (i+2, N, table+2);
			gwsincos1by4 (i+3, N, table+3);
			table += 8;
		}
		N = N / 3;
	}

/* For initial level radix-5 building block, output two sin/cos values. */

	while (N % 5 == 0) {
		for (i = 0; i < N / 5; i += 4) {
			gwsincos12by4 (i, N, table);
			gwsincos12by4 (i+1, N, table+1);
			gwsincos12by4 (i+2, N, table+2);
			gwsincos12by4 (i+3, N, table+3);
			table += 16;
		}
		N = N / 5;
	}

/* For the first level radix-4 block, output two sin/cos values. */

	for (i = 0; i < N / 4; i += 4) {
		gwsincos12by4 (i, N, table);
		gwsincos12by4 (i+1, N, table+1);
		gwsincos12by4 (i+2, N, table+2);
		gwsincos12by4 (i+3, N, table+3);
		table += 16;
	}
	N = N / 4;

/* Build a smaller table for the remaining FFT levels */

	for (i = N, pow2_count = 0; (i & 1) == 0; i >>= 1) pow2_count++;
	if (pow2_count & 1) N = 128;		/* Last seven levels done in L1 cache */
	else if (N > 256) N = 256;		/* Last eight levels done in L1 cache */
	for (i = 0; i < 4; i++) {
		for (j = 0; j < N / 4; j += 4) {
			gwsincos12by1 (i + j, N, table);
			table += 4;
		}
	}

/* Return address of the end of the table */

	return (table);
}

/* This routine builds the sin/cos table used in pass 2 of real FFTs, */

double *yr4_build_pass2_real_table (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	double	*table)		/* Pointer to the table to fill in */
{
	unsigned long avx_increment, j, N;
	int	swizzled;

/* All complex FFTs, don't need these tables */

	if (gwdata->ALL_COMPLEX_FFT) return (table);

/* Init */

	avx_increment = 1;
	N = gwdata->PASS2_SIZE;
	swizzled = FALSE;

/* Output sin/cos values for an initial 6-real, 10-real, or 8-real macro. */

	while (N % 3 == 0) {
		for (j = 0; j < N / 3; j += 4) {
			gwsincos1by4 (j, N*2, table);		/* For the six real (and three complex) */
			gwsincos1by4 (j + avx_increment, N*2, table+1);
			gwsincos1by4 (j + 2*avx_increment, N*2, table+2);
			gwsincos1by4 (j + 3*avx_increment, N*2, table+3);
			table += 8;
		}
		N = N / 3;
	}

	while (N % 5 == 0) {
		for (j = 0; j < N / 5; j += 4) {
			gwsincos13by4 (j, N*2, table);	/* For the ten real (and five complex) */
			gwsincos13by4 (j + avx_increment, N*2, table+1);
			gwsincos13by4 (j + 2*avx_increment, N*2, table+2);
			gwsincos13by4 (j + 3*avx_increment, N*2, table+3);
			table += 16;
		}
		N = N / 5;
	}

	for (j = 0; j < N / 4; j += 4) {
		gwsincos15by4 (j, N*2, table);
		gwsincos15by4 (j + avx_increment, N*2, table+1);
		gwsincos15by4 (j + 2*avx_increment, N*2, table+2);
		gwsincos15by4 (j + 3*avx_increment, N*2, table+3);
		table += 16;
	}
	N = N / 4;

/* Output one last sin/cos table for the remaining yr4_eight_reals_four_complex_djbfft */
/* building block levels.  The eight_reals doubles N because the real part of the FFT */
/* is one level behind the complex part of the FFT.  The four-complex sin/cos values */
/* are the same for all 3 of the upper YMM doubles. */		

	avx_increment = N * 2;
	for (j = 0; j < N / 4; j++) {
		gwsincos125by4 (j, N*2, table);				/* For the eight_reals */
		gwsincos12by4 (j + avx_increment, N, table+1);		/* For the four-complex */
		table[3] = table[2] = table[1];
		table[7] = table[6] = table[5];
		table[11] = table[10] = table[9];
		table[15] = table[14] = table[13];
		table[19] = table[18] = table[17] = -table[1];
		table[23] = table[22] = table[21] = -table[5];
		table += 24;
	}

/* Return address of the end of the table */

	return (table);
}

/* This routine builds a normalization table - used by radix-4 with partial */
/* normalization FFTs. */

double *yr4dwpn_build_norm_table (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	double	*table)		/* Pointer to the table to fill in */
{
	unsigned long i, j, num_cols, upper_avx_word;

/* Build the group multipliers table */

	upper_avx_word = gwdata->PASS2_SIZE;
	num_cols = gwdata->PASS2_SIZE * gwdata->wpn_count;
	for (i = 0; i < gwdata->FFTLEN / 2; i += num_cols) {
		for (j = 0; j < gwdata->FFTLEN; j += gwdata->FFTLEN / 2) {
			double	double1_weight, double2_weight, double3_weight, double4_weight;
			int	double1_sort, double2_sort, double3_sort, double4_sort;
			double	ttp, ttmp, ttp_over_b, ttmp_times_b, *tab20;

/* For zero-padded FFTs the upper half of the FFT has the exact same multipliers as the lower half. */
/* Thus we can cut the size of our group multiplier table in half. */

			if (gwdata->ZERO_PADDED_FFT && j >= gwdata->FFTLEN / 2) continue;

/* The sort order of the weights determines which fudge factor combination never occur. */

			double1_weight = gwfft_weight_exponent (gwdata->dd_data, i + j);
			double2_weight = gwfft_weight_exponent (gwdata->dd_data, i + j + upper_avx_word);
			double3_weight = gwfft_weight_exponent (gwdata->dd_data, i + j + 2 * upper_avx_word);
			double4_weight = gwfft_weight_exponent (gwdata->dd_data, i + j + 3 * upper_avx_word);

/* Now sort the four weights */

			double1_sort = (double1_weight > double2_weight) +
				       (double1_weight > double3_weight) +
				       (double1_weight > double4_weight);
			double2_sort = (double2_weight > double1_weight) +
				       (double2_weight > double3_weight) +
				       (double2_weight > double4_weight);
			double3_sort = (double3_weight > double1_weight) +
				       (double3_weight > double2_weight) +
				       (double3_weight > double4_weight);
			double4_sort = (double4_weight > double1_weight) +
				       (double4_weight > double2_weight) +
				       (double4_weight > double3_weight);

/* Call quad-precision routine to compute set of multipliers. */
/* We compute two-to-phi multiplied by the fudge factor so that normalize won't have to. */

			gwfft_weights_fudged (gwdata->dd_data, i + j, gwdata->b, &ttp, &ttmp, &ttp_over_b, &ttmp_times_b);

/* Set the LSW entries in an AVX word */

			table[0] = ttmp;
			table[4] = double1_sort < 3 ? ttmp : ttmp_times_b;
			table[8] = double1_sort < 2 ? ttmp : ttmp_times_b;
			table[12] = double1_sort < 1 ? ttmp : ttmp_times_b;
			table[16] = ttmp_times_b;
			tab20 = table + 20;
			tab20[0] = ttp;
			tab20[4] = double1_sort < 3 ? ttp : ttp_over_b;
			tab20[8] = double1_sort < 2 ? ttp : ttp_over_b;
			tab20[12] = double1_sort < 1 ? ttp : ttp_over_b;
			tab20[16] = ttp_over_b;

/* Set the remaining entries in an AVX word */

			gwfft_weights_fudged (gwdata->dd_data, i + j + upper_avx_word, gwdata->b, &ttp, &ttmp, &ttp_over_b, &ttmp_times_b);
			table[1] = ttmp;
			table[5] = double2_sort < 3 ? ttmp : ttmp_times_b;
			table[9] = double2_sort < 2 ? ttmp : ttmp_times_b;
			table[13] = double2_sort < 1 ? ttmp : ttmp_times_b;
			table[17] = ttmp_times_b;
			tab20[1] = ttp;
			tab20[5] = double2_sort < 3 ? ttp : ttp_over_b;
			tab20[9] = double2_sort < 2 ? ttp : ttp_over_b;
			tab20[13] = double2_sort < 1 ? ttp : ttp_over_b;
			tab20[17] = ttp_over_b;

			gwfft_weights_fudged (gwdata->dd_data, i + j + 2 * upper_avx_word, gwdata->b, &ttp, &ttmp, &ttp_over_b, &ttmp_times_b);
			table[2] = ttmp;
			table[6] = double3_sort < 3 ? ttmp : ttmp_times_b;
			table[10] = double3_sort < 2 ? ttmp : ttmp_times_b;
			table[14] = double3_sort < 1 ? ttmp : ttmp_times_b;
			table[18] = ttmp_times_b;
			tab20[2] = ttp;
			tab20[6] = double3_sort < 3 ? ttp : ttp_over_b;
			tab20[10] = double3_sort < 2 ? ttp : ttp_over_b;
			tab20[14] = double3_sort < 1 ? ttp : ttp_over_b;
			tab20[18] = ttp_over_b;

			gwfft_weights_fudged (gwdata->dd_data, i + j + 3 * upper_avx_word, gwdata->b, &ttp, &ttmp, &ttp_over_b, &ttmp_times_b);
			table[3] = ttmp;
			table[7] = double4_sort < 3 ? ttmp : ttmp_times_b;
			table[11] = double4_sort < 2 ? ttmp : ttmp_times_b;
			table[15] = double4_sort < 1 ? ttmp : ttmp_times_b;
			table[19] = ttmp_times_b;
			tab20[3] = ttp;
			tab20[7] = double4_sort < 3 ? ttp : ttp_over_b;
			tab20[11] = double4_sort < 2 ? ttp : ttp_over_b;
			tab20[15] = double4_sort < 1 ? ttp : ttp_over_b;
			tab20[19] = ttp_over_b;

			table += 40;
		}
	}
	return (table);
}

/* This routine builds the big/little flags table for an AVX r4dwpn (radix-4 with partial normalization) FFT */

double *yr4dwpn_build_biglit_table (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	double	*table)		/* Pointer to the table to fill in */
{
	struct gwasm_data *asm_data = (struct gwasm_data *) gwdata->asm_data;
const	int	START_OF_CHAIN = 0x8000;
const	int	END_OF_CHAIN = 0x4000;
	int	combos[16], next_combo[16];
	unsigned char combo_recorded[256];
	unsigned char *p;
	unsigned long upper_avx_word, fftlen_over_2, num_cols;
	unsigned long group, i, j, k, m, n, num_combos;

/* Init table of first 8 big/lit values */

	memset (asm_data->u.ymm.YMM_FIRST_BIGLIT_VALUES, 0, sizeof (asm_data->u.ymm.YMM_FIRST_BIGLIT_VALUES));

/* Big/lit flags form a very regular pattern.  For example, if there are 18.3 b's */
/* per FFT word then you get either a big word followed by two little words or a */
/* big  followed by three little words.  Here we determine which patterns of big/lit */
/* are possible in a ynorm_wpn macro which processes 2 AVX words. */

/* Generate all possible valid combinations of big/lit flags */

	memset (combo_recorded, 0, sizeof (combo_recorded));
	upper_avx_word = gwdata->PASS2_SIZE;
	fftlen_over_2 = gwdata->FFTLEN / 2;
	num_combos = 0;

	/* Loop over all cache lines */
	for (i = 0; i < upper_avx_word; i++) {
		for (j = 0; j < fftlen_over_2; j += 4*upper_avx_word) {
			int	combo;

			/* Generate combo for this cache line */
			combo = (is_big_word (gwdata, i + j + 3*upper_avx_word) << 7) +
				(is_big_word (gwdata, i + j + 2*upper_avx_word) << 6) +
				(is_big_word (gwdata, i + j + upper_avx_word) << 5) +
				(is_big_word (gwdata, i + j) << 4) +
				(is_big_word (gwdata, i + j + fftlen_over_2 + 3*upper_avx_word) << 3) +
				(is_big_word (gwdata, i + j + fftlen_over_2 + 2*upper_avx_word) << 2) +
				(is_big_word (gwdata, i + j + fftlen_over_2 + upper_avx_word) << 1) +
				(is_big_word (gwdata, i + j + fftlen_over_2));

			/* Ignore this combo if it is a duplicate.  Otherwise, add it to our combos collection. */
			if (! combo_recorded[combo]) {
				combo_recorded[combo] = 1;
				combos[num_combos++] = combo;
			}
		}
	}

/* Concatentate combos to save space.  Let's hope they fit in 48 entries. */

	/* Init the next-in-chain array */
	for (i = 0; i < num_combos; i++)
		next_combo[i] = START_OF_CHAIN + END_OF_CHAIN;

	/* Look for 2 chains where the end of one chain has elements in common with */
	/* the start of the other chain. */

	/* Examine all chain starts */
	for (i = 0; i < num_combos; i++) {
		int	chain_end;

		/* Skip if not the start of a chain */
		if (! (next_combo[i] & START_OF_CHAIN)) continue;

		/* Find end of chain */
		for (chain_end = i; ! (next_combo[chain_end] & END_OF_CHAIN); chain_end = next_combo[chain_end] & 0xFF);

		/* Now look at all chain ends */
		for (j = 0; j < num_combos; j++) {

			/* Skip if not a chain end */
			if (! (next_combo[j] & END_OF_CHAIN)) continue;
			/* Can't chain to ourselves! */
			if (j == chain_end) continue;

			/* See if chain end has common elements with the chain start */
			if ((combos[i] >> 4) == (combos[j] & 0x0F)) {
				next_combo[j] = (next_combo[j] & START_OF_CHAIN) + i;
				next_combo[i] &= ~START_OF_CHAIN;
				break;
			}
		}
	}

/* HACK: Remember the chains in ASM_TIMERS so that we can properly build LIMIT_INVERSE */
/* and LIMIT_BIGMAX at a later time. */

	n = 0;
	for (i = 0; i < num_combos; i++) {
		if (! (next_combo[i] & START_OF_CHAIN)) continue;

		/* Output up to 4 bytes for each entry in the chain */
		for (j = i; ; j = next_combo[j] & 0xFF) {
			((char *)gwdata->ASM_TIMERS)[n++] = (combos[j] >> 4);
			if (next_combo[j] & END_OF_CHAIN) {
				((char *)gwdata->ASM_TIMERS)[n++] = combos[j] & 0xF;
				break;
			}
		}
	}
GWASSERT (n <= 24);	// Lets see what the AVX limits really are!!!
	GWASSERT (n <= 48);
	if (n > 48) gwdata->GWERROR = GWERROR_INTERNAL + 2;

/* Determine the number of column two-to-phi multipliers */

	num_cols = gwdata->PASS2_SIZE * gwdata->wpn_count;

/* Loop to build table in exactly the same order that it will be */
/* used by the assembly code. */

	p = (unsigned char *) table;
	for (group = 0; group < upper_avx_word; group += gwdata->PASS1_CACHE_LINES) {
	    for (n = 0; n < gwdata->FFTLEN / 2; n += num_cols) {
	        for (j = 0; j < num_cols; j += gwdata->PASS2_SIZE * 4) {
		    for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
			for (k = 0; k < gwdata->FFTLEN; k += gwdata->FFTLEN / 2) {
				unsigned long word;

/* Now set the big/little flag for a LSW in an AVX pair */
/* Otherwise, set the big/little flag for a MSW in an AVX pair */

				word = group + j + n + i + k;
				*p = is_big_word (gwdata, word);
				if (is_big_word (gwdata, word + upper_avx_word)) *p += 2;
				if (is_big_word (gwdata, word + 2 * upper_avx_word)) *p += 4;
				if (is_big_word (gwdata, word + 3 * upper_avx_word)) *p += 8;
   
/* Set the ttp and ttmp fudge flags for two pass FFTs.  The fudge flag is */
/* set if the col mult * the grp mult is b times the correct fft_weight, */
/* meaning a mul by 1/b is required to generate the correct multiplier. */
/* Since we can't do equality compares on floats, this test is a little bit */
/* cryptic. */

				if (gwfft_weight_exponent (gwdata->dd_data, word) + 0.5 <
				    gwfft_weight_exponent (gwdata->dd_data, n + k) +
				    gwfft_weight_exponent (gwdata->dd_data, group + j + i))
					*p += 16;
				if (gwfft_weight_exponent (gwdata->dd_data, word + upper_avx_word) + 0.5 <
				    gwfft_weight_exponent (gwdata->dd_data, n + k + upper_avx_word) +
				    gwfft_weight_exponent (gwdata->dd_data, group + j + i))
					*p += 32;
				if (gwfft_weight_exponent (gwdata->dd_data, word + 2 * upper_avx_word) + 0.5 <
				    gwfft_weight_exponent (gwdata->dd_data, n + k + 2 * upper_avx_word) +
				    gwfft_weight_exponent (gwdata->dd_data, group + j + i))
					*p += 64;
				if (gwfft_weight_exponent (gwdata->dd_data, word + 3 * upper_avx_word) + 0.5 <
				    gwfft_weight_exponent (gwdata->dd_data, n + k + 3 * upper_avx_word) +
				    gwfft_weight_exponent (gwdata->dd_data, group + j + i))
					*p += 128;

/* Apply our method for reducing fudge factor data from 16 combinations down to 5 possibilities. */

				if (*p & 32) *p -= 16;
				if (*p & 128) *p -= 64;
				*p = ((*p >> 6) + ((*p >> 4) & 0x3)) * 16 + (*p & 0xF);
				p++;
			}

/* Combine last two big/lit 4-bit flag values into one 6-bit flags value. */

			p -= 2;
			for (m = 0; m <= 46; m++) {
				if (((char *)gwdata->ASM_TIMERS)[m]   == (p[0] & 0xF) &&
				    ((char *)gwdata->ASM_TIMERS)[m+1] == (p[1] & 0xF))
					break;
			}
			ASSERTG (m != 47);
			if (m == 47) gwdata->GWERROR = GWERROR_INTERNAL + 3;

/* Combine 1st and 2nd fudge factor flags into one byte. */

			p[0] = ((p[0] >> 4) << 5) + ((p[1] >> 4) << 2);

/* Output the 6-bit big/lit flags value */

			p[1] = (unsigned char) (m << 2);

/* Create separate table for first 8 biglit values for carry propagation code */

			if (group + n + j + i < sizeof (asm_data->u.ymm.YMM_FIRST_BIGLIT_VALUES))
				asm_data->u.ymm.YMM_FIRST_BIGLIT_VALUES[group + n + j + i] = p[1];

/* Move pointer to next big/lit table entry */

			p += 2;
		    }
		}
	    }
	}
	return ((double *) p);
}

/* This routine builds the sin/cos table used in pass 1 by a traditional */
/* DJB radix-4 FFT  - called by gwsetup.  If this is an all-complex FFT, */
/* then the root-of-minus-1 premultipliers are also built. */

double *r4_build_pass1_table (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	double	*table)		/* Pointer to the table to fill in */
{
	unsigned long pass1_size, pass1_increment, first_levels_size;
	unsigned long group, i, j, k, N, temp, upper_sse2_word;
	int	pow2_count;

/* Initialize some needed constants */

	pass1_increment = gwdata->PASS2_SIZE;
	upper_sse2_word = pass1_increment / 2;
	pass1_size = gwdata->FFTLEN / gwdata->PASS2_SIZE; /* Real values in a pass1 section */

/* Determine size of the first levels. */

	if (pass1_size % 7 == 0)
		first_levels_size = 28;
	else if (pass1_size % 5 == 0 && !gwdata->ALL_COMPLEX_FFT)
		first_levels_size = 20;
	else
		first_levels_size = 8;

/* Count the power-of-two FFT levels after the initial FFT levels.  If odd, the */
/* the last levels will be a radix-8, if even the last levels will be a radix-4. */

	pass1_size /= first_levels_size;
	for (pow2_count = 0; (pass1_size & 1) == 0; pass1_size /= 2) pow2_count++;

/* Set pointer to table of multipliers */

	gwdata->adjusted_pass1_premults = table;

/* Loop through all the pass 1 groups in the same order the assembly code will */
/* process the groups. */

	for (group = 0; group < upper_sse2_word; group += gwdata->PASS1_CACHE_LINES) {

		pass1_size = gwdata->FFTLEN / gwdata->PASS2_SIZE; /* Real values in a pass1 section */
		pass1_size /= first_levels_size; /* Complex values we're generating sin/cos data for */
		N = gwdata->PASS2_SIZE;

/* Output the sin/cos/premultiplier values for the radix-8 block that does the */
/* last 3 levels in pass 1.  NOTE:  We do not need the "j loop" (it would loop */
/* from zero to zero) when generating the sin/cos twiddle factors for the last */
/* levels of pass 1. */

		if (pow2_count & 1) {
			N = N * 8;

/* For the r4_g8cl_sixteen_reals_eight_complex_djbfft building block, output the extra */
/* sin/cos values needed for the sixteen_reals.  The sixteen_reals is done in three parts: */
/* 4 four_reals_fft, 1 eight_reals_fft, and 1 four-complex_fft.  The four_reals_fft needs */
/* four sin/cos twiddles using 2*N, the eight_reals_fft and four_complex_fft can use the same */
/* sin/cos twiddles that will be generated later for the r4r8_g8cl_eight_complex_fft8 macro. */

			if (!gwdata->ALL_COMPLEX_FFT) {
				for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
					temp = (group + i);
					gwsincos1plus0123by2 (temp, pass1_increment, N*2, table);
					gwsincos1plus0123by2 (temp + upper_sse2_word, pass1_increment, N*2, table+1);
					table += 16;
				}
			}

/* Output the sin/cos values for the complex group -- specifically the r8_sg8cl_eight_complex_djbfft macro. */

			for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {

/* We now calculate the standard sin/cos twiddle factors (temp*0,1,2,3,4,5,6,7 roots of N) */

				temp = group + i;
				gwsincos1plus01234567by2 (0, temp, N, table);
				gwsincos1plus01234567by2 (0, temp + upper_sse2_word, N, table+1);
				table += 32;
			}
			pass1_size /= 8;
		}

/* For the r4_four_complex_djbfft building block levels, output the sin/cos values. */

		while ((pass1_size & 3) == 0) {
			N = N * 4;

/* For the r4_eight_reals_four_complex_djbfft building block levels, output the extra */
/* sin/cos values needed for the eight_reals.  The eight_reals doubles N because */
/* the real part of the FFT is one level behind the complex part of the FFT. */

			if (!gwdata->ALL_COMPLEX_FFT) {
				for (j = 0; j < N*2 / 8; j += pass1_increment) {
					for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
						temp = (group + j + i);
						gwsincos15by2 (temp, N*2, table);
						gwsincos15by2 (temp + upper_sse2_word, N*2, table+1);
						table += 8;
					}
				}
			}

/* Output the sin/cos value for the complex sections, used by the r4_four_complex_djbfft macro */

			for (j = 0; j < N / 4; j += pass1_increment) {
				for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
					temp = (group + j + i);
					gwsincos12by2 (temp, N, table);
					gwsincos12by2 (temp + upper_sse2_word, N, table+1);
					table += 8;
				}
			}
			pass1_size /= 4;
		}

/* For the r5_five_complex_djbfft building block levels, output the sin/cos values. */

		while (pass1_size % 5 == 0) {
			N = N * 5;

/* The r5_ten_reals building blocks require extra sin/cos */
/* values.  The ten_reals doubles N because the real part of the */
/* FFT is one level behind the complex part of the FFT. */

			if (!gwdata->ALL_COMPLEX_FFT) {
				for (j = 0; j < N / 5; j += pass1_increment) {
					for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
						temp = (group + j + i);
						gwsincos13by2 (temp, N*2, table);
						gwsincos13by2 (temp + upper_sse2_word, N*2, table+1);
						table += 8;
					}
				}
			}

/* Output the sin/cos data for the complex sections, (the r5_five_complex_djbfft building block). */

			for (j = 0; j < N / 5; j += pass1_increment) {
				for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
					temp = (group + j + i);
					gwsincos12by2 (temp, N, table);
					gwsincos12by2 (temp + upper_sse2_word, N, table+1);
					table += 8;
				}
			}
			pass1_size /= 5;
		}

/* For the r3_three_complex_djbfft building block levels, output the sin/cos values. */

		while (pass1_size % 3 == 0) {
			N = N * 3;

/* The r3_six_reals building blocks require an extra sin/cos */
/* value.  The six_reals doubles N because the real part of the */
/* FFT is one level behind the complex part of the FFT. */

			if (!gwdata->ALL_COMPLEX_FFT) {
				for (j = 0; j < N / 3; j += pass1_increment) {
					for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
						temp = (group + j + i);
						gwsincos1by2 (temp, N*2, table);
						gwsincos1by2 (temp + upper_sse2_word, N*2, table+1);
						table += 4;
					}
				}
			}

/* Output the sin/cos data for the complex sections (used by the r3_three_complex_djbfft building block). */

			for (j = 0; j < N / 3; j += pass1_increment) {
				for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
					temp = (group + j + i);
					gwsincos1by2 (temp, N, table);
					gwsincos1by2 (temp + upper_sse2_word, N, table+1);
					table += 4;
				}
			}
			pass1_size /= 3;
		}
		ASSERTG (pass1_size == 1);
		if (pass1_size != 1) gwdata->GWERROR = GWERROR_INTERNAL + 4;

/* Real FFTs output one last set of sin/cos values for the first 20-reals, 28-reals, or 8-reals FFT. */

		if (! gwdata->ALL_COMPLEX_FFT) {
			N = gwdata->FFTLEN;
			pass1_size = gwdata->FFTLEN / gwdata->PASS2_SIZE; /* Real values in a pass1 section */
			if (pass1_size % 20 == 0) {
				for (j = 0; j < N / 20; j += pass1_increment) {
					for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
						temp = (group + j + i);
						for (k = 1; k <= 9; k++) {	/* Create 9 twiddle factors */
							gwsincos1by2 (k * temp, N, table);
							gwsincos1by2 (k * (temp + upper_sse2_word), N, table+1);
							table += 4;
						}
					}
				}
			}
			else if (pass1_size % 28 == 0) {
				for (j = 0; j < N / 28; j += pass1_increment) {
					for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
						temp = (group + j + i);
						for (k = 1; k <= 13; k++) {	/* Create 13 twiddle factors */
							gwsincos1by2 (k * temp, N, table);
							gwsincos1by2 (k * (temp + upper_sse2_word), N, table+1);
							table += 4;
						}
					}
				}
			}
			else {
				for (j = 0; j < N / 8; j += pass1_increment) {
					for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
						temp = (group + j + i);
						gwsincos125by2 (temp, N, table);
						gwsincos125by2 (temp + upper_sse2_word, N, table+1);
						table += 12;
					}
				}
			}
		}

/* For all-complex FFTs, the first FFT level converts from real values to all complex */
/* values by multiplying by a root of -1 weight and doing a butterfly.  This is */
/* simplified because weights in the bottom half are sqrt(-1) times the */
/* matching weights in the upper half.  Thus, we butterfly upper_word * weight */
/* with bottom_word * i * weight.  That equals (upper_word + i * lower_word) * weight */
/* That is just a complex multiply with the half of the butterfly output values */
/* unneeded thanks to Hermetian symmetry. */

/* The second and third levels do a standard radix-4 four-complex-FFT */
/* building block with a post-multiply by a sin/cos root of unity. */

		else {
			N = gwdata->FFTLEN / 2;
			for (j = 0; j < N / 4; j += pass1_increment) {
				for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
					temp = (group + j + i);

/* Here we compute the standard 0,1,2,3 * temp for the radix-4 sin/cos multipliers. */
/* Then we multiply in the roots-of-minus-one premultiplier.  The root-of-minus-one */
/* premultiplier was for 2N, and a root-of-minus-one-of-2N is the same as a root */
/* unity for 4N. */

					gwsincos1plus0123by2 (temp, 4 * temp, N * 4, table); /* premult + temp*0-3 */
					temp += upper_sse2_word;
					gwsincos1plus0123by2 (temp, 4 * temp, N * 4, table+1);
					table += 16;
				}
			}
		}

/* Calculate the size of each group's sin/cos/premult data for pass1_get_next_block */

		if (group == 0)
			gwdata->pass1_premult_block_size = (unsigned long)
				((char *) table - (char *) gwdata->adjusted_pass1_premults) /
				gwdata->PASS1_CACHE_LINES;
	}

/* Return address of the end of the table */

	return (table);
}

/* This routine builds the sin/cos table used in all-complex pass 2 */
/* blocks in a traditional radix-4 FFT - called by gwsetup. */

double *r4_build_pass2_complex_table (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	double	*table)		/* Pointer to the table to fill in */
{
	unsigned long i, N, limit, aux_table_size;

/* We also build a smaller auxiliary table so the final several levels aren't using */
/* cache-unfriendly large strides to access sin/cos data.  If the last levels use */
/* eight_complex macros set auxiliary table size to 512, otherwise the last levels */
/* use the four_complex macros and we'll build an auxiliary table size of 256. */

	if (gwdata->PASS2_SIZE <= 640 && gwdata->PASS2_SIZE % 3 != 0)
		aux_table_size = 0;
	else {
		for (i = 0, N = gwdata->PASS2_SIZE; (N & 1) == 0; N >>= 1) i++;
		if (i < 8)
			aux_table_size = 1 << i;
		else
			aux_table_size = (i & 1) ? 512 : 256;
	}
	
/* If the pass 2 size is divisible by 3, then the first level does a */
/* radix-3 step which only requires one sin/cos value.  Alas, rather than */
/* computing the needed N/3 sin/cos values we must compute 1/2 or 2/5 times N */
/* sin/cos values for use in later r4_x4cl_2sc_four_complex_djbfft or */
/* r5_x5cl_2sc_five_complex_djbfft levels. */

	N = gwdata->PASS2_SIZE;
	if (N % 3 == 0) {
		if (!aux_table_size ||  N / aux_table_size % 2 == 0) limit = N / 2;
		else if (N % 5 == 0) limit = N * 2 / 5;
		else limit = N / 3;
		for (i = 0; i < limit; i++) {
			gwsincos1by2 (i, N, table);
			table[1] = table[0];
			table[3] = table[2];
			table += 4;
		}
	}

/* For the radix-4 and radix-5 building blocks, output two sin/cos values. */
/* If the levels above the aux_table_size are all radix-5, then we can */
/* output a slightly smaller sin/cos table. */

	else {
		if (!aux_table_size ||  N / aux_table_size % 2 == 0) limit = N / 4;
		else limit = N / 5;
		for (i = 0; i < limit; i++) {
			gwsincos12by2 (i, N, table);
			table[1] = table[0];	/* temp */
			table[3] = table[2];
			table[5] = table[4];	/* 2 * temp */
			table[7] = table[6];
			table += 8;
		}
	}

/* Build the smaller auxiliary table */

	if (aux_table_size) {
		N = aux_table_size;
		for (i = 0; i < N / 4; i++) {
			gwsincos12by2 (i, N, table);
			table[1] = table[0];	/* temp */
			table[3] = table[2];
			table[5] = table[4];	/* 2 * temp */
			table[7] = table[6];
			table += 8;
		}
	}

/* Return address of the end of the table */

	return (table);
}

/* This routine builds the sin/cos table used in pass 2 of real FFTs, */

double *r4_build_pass2_real_table (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	double	*table)		/* Pointer to the table to fill in */
{
	unsigned long i, N;
	double	sincos[6];

/* All complex FFTs, don't need these tables */

	if (gwdata->ALL_COMPLEX_FFT) return (table);

/* Real FFTs with an initial radix-3 step needs w^n for 2N. */

	N = gwdata->PASS2_SIZE;
	if (N % 3 == 0) {
		for (i = 0; i < N / 3; i++) {
			gwsincos (i, N * 2, (double *) &sincos);
			table[0] = sincos[0];	/* temp for 2*N */
			table[1] = sincos[1];
			table += 2;
		}
		while (N % 3 == 0) N = N / 3;
	}

/* Real FFTs with a radix-5 step needs w^n and w^3n for 2N. */

	if (N % 5 == 0) {
		for (i = 0; i < N / 5; i++) {
			gwsincos3 (i, N * 2, (double *) &sincos);
			table[0] = sincos[0];	/* temp for 2*N */
			table[1] = sincos[1];
			table[2] = sincos[4];	/* 3 * temp for 2*N */
			table[3] = sincos[5];
			table += 4;
		}
		while (N % 5 == 0) N = N / 5;
	}

/* Real FFTs with eight_reals macros need w^n and w^5n for 2N. */

	for (i = 0; i < N / 4; i++) {
		gwsincos5 (i, N * 2, (double *) &sincos);
		table[0] = sincos[0];	/* temp for 2*N */
		table[1] = sincos[1];
		table[2] = sincos[2];	/* 5 * temp for 2*N */
		table[3] = sincos[3];
		table += 4;
	}

/* Return address of the end of the table */

	return (table);
}


/* This routine builds the fixed sin/cos premultiplier table used in pass 1 by the radix-4 */
/* delayed multiplier FFT - called by gwsetup.  It is like the traditional radix-4 FFT */
/* but with the all-complex premultipliers split into two parts to reduce memory. */

double *r4delay_build_fixed_premult_table (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	double	*table)		/* Pointer to the table to fill in */
{
	unsigned long pass1_increment, j, N;

/* For all-complex FFTs, build the fixed roots-of-minus-one table. */
/* Output these roots in the same order they will be used the first two */
/* levels of pass 1. */

	if (gwdata->ALL_COMPLEX_FFT) {
		N = gwdata->FFTLEN / 2;
		pass1_increment = gwdata->PASS2_SIZE;
		for (j = 0; j < N / 4; j += pass1_increment) {

/* Here we compute the roots-of-minus-one premultiplier.  The root-of-minus-one */
/* premultiplier was for 2N, and a root-of-minus-one-of-2N is the same as a root */
/* unity for 4N. */					

			gwsincos1plus0123by2 (j, N / 4, N * 4, table);
			table[1] = table[0];	/* premult  */
			table[3] = table[2];
			table[5] = table[4];	/* premult * -1 ^ 1/8 */
			table[7] = table[6];
			table[9] = table[8];	/* premult * -1 ^ 2/8 */
			table[11] = table[10];
			table[13] = table[12];	/* premult * -1 ^ 3/8 */
			table[15] = table[14];
			table += 16;
		}
	}

/* Return address of the end of the table */

	return (table);
}

/* This routine builds the sin/cos table used in pass 1 by the radix-4/8 DJB */
/* FFT with delayed sin/cos multiplies. */

double *r4delay_build_pass1_table (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	double	*table)		/* Pointer to the table to fill in */
{
	unsigned long pass1_size, pass1_increment, delay_count;
	unsigned long group, i, j, k, N, temp, upper_sse2_word;
	int	pow2_count;

/* Initialize some needed constants */

	pass1_increment = gwdata->PASS2_SIZE;
	upper_sse2_word = pass1_increment / 2;
	pass1_size = gwdata->FFTLEN / gwdata->PASS2_SIZE; /* Real values in a pass1 section */

/* Determine number of delay groups.  In a standard radix-4 FFT, there is only one sin/cos */
/* group in the last pass 1 level.  We reduce our memory usage by using just one fixed sin/cos */
/* table in the first FFT levels and having multiple groups of sin/cos data in the last pass 1 level. */
/* I call these groups of sin/cos data in the last pass 1 level "delay groups". */

	if (pass1_size % 7 == 0)
		delay_count = 14;
	else if (pass1_size % 5 == 0 && !gwdata->ALL_COMPLEX_FFT)
		delay_count = 10;
 	else if ((pass1_size == 512 && !gwdata->ALL_COMPLEX_FFT) ||
		 (pass1_size == 1024 && !gwdata->ALL_COMPLEX_FFT) ||
		 (pass1_size == 2560 && gwdata->ALL_COMPLEX_FFT) ||
		 (pass1_size == 5120 && gwdata->ALL_COMPLEX_FFT) ||
		 pass1_size == 1536 || pass1_size == 2048 || pass1_size == 3072 || pass1_size == 4096)
		delay_count = 16;
	else
		delay_count = 4;

/* Count the power-of-two FFT levels after the initial FFT levels.  If odd, the */
/* the last levels will be a radix-8, if even the last levels will be a radix-4. */

	pass1_size /= (delay_count * 2);
	for (pow2_count = 0; (pass1_size & 1) == 0; pass1_size /= 2) pow2_count++;

/* Set pointer to table of multipliers */

	gwdata->adjusted_pass1_premults = table;

/* Loop through all the pass 1 groups in the same order the assembly code will */
/* process the groups. */

	for (group = 0; group < upper_sse2_word; group += gwdata->PASS1_CACHE_LINES) {

		pass1_size = gwdata->FFTLEN / gwdata->PASS2_SIZE; /* Real values in a pass1 section */
		pass1_size /= (delay_count * 2);	/* Complex values we're generating sin/cos data for */
		N = gwdata->PASS2_SIZE;

/* Output the sin/cos/premultiplier values for the radix-8 block that does the */
/* last 3 levels in pass 1.  NOTE:  We do not need the "j loop" (it would loop */
/* from zero to zero) when generating the sin/cos twiddle factors for the last */
/* levels of pass 1. */

		if (pow2_count & 1) {
			N = N * 8;

/* For the r4_g8cl_sixteen_reals_eight_complex_djbfft building block, output the extra */
/* sin/cos values needed for the sixteen_reals.  The sixteen_reals is done in three parts: */
/* 4 four_reals_fft, 1 eight_reals_fft, and 1 four-complex_fft.  The four_reals_fft needs */
/* four sin/cos twiddles using 2*N, the eight_reals_fft and four_complex_fft can use the same */
/* sin/cos twiddles that will be generated later for the r4r8_g8cl_eight_complex_fft8 macro. */

			if (!gwdata->ALL_COMPLEX_FFT) {
				for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
					temp = (group + i);
					gwsincos1plus0123by2 (temp, pass1_increment, N*2, table);
					gwsincos1plus0123by2 (temp + upper_sse2_word, pass1_increment, N*2, table+1);
					table += 16;
				}
			}

/* Output the sin/cos values for the complex delay groups -- specifically the r8_sg8cl_eight_complex_fft8 macro. */

			for (k = 0; k < delay_count; k++) {
				for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
					unsigned long bigN, ktemp, actemp;

/* If this is an all-complex FFT, the roots of minus 1 (same as roots of FFTLEN*2) are */
/* split to reduce memory requirements.  We apply the part of the all-complex premultiplier here. */

					if (gwdata->ALL_COMPLEX_FFT) {
						bigN = gwdata->FFTLEN * 2;
						actemp = group + i;
					} else {
						bigN = gwdata->FFTLEN;
						actemp = 0;
					}

/* Factor in the delayed part of the sin/cos multiplies from the first 2 levels.  In the first 2 levels */
/* we use a fixed sin/cos table based only on j, leaving the group+i part to be applied here by */
/* creating delay_count table entries.  For an all-complex FFT we multiply by 0,2,1,-1*(group+i) with N = FFTLEN/2. */
/* For an all-real FFT we multiply by 0,2,1,5*(group+i) with N = FFTLEN. */

					if (gwdata->ALL_COMPLEX_FFT && delay_count == 4) {
						if (k == 0)
							ktemp = 0;
						else if (k == 1)
							ktemp = 2 * (group + i) * 4;
						else if (k == 2)
							ktemp = 1 * (group + i) * 4;
						else
							ktemp = bigN - 1 * (group + i) * 4;
					} else if (gwdata->ALL_COMPLEX_FFT) {
						/* 0,2,1,-1 combined with 0,2,1,-1 */
						int	kmap[16] = {0,8,4,-4,2,10,6,-2,1,9,5,-3,-1,7,3,-5};
						if (kmap[k] >= 0)
							ktemp = kmap[k] * (group + i) * 4;
						else
							ktemp = bigN + kmap[k] * (group + i) * 4;
					} else if (delay_count == 4) {
						if (k == 0)
							ktemp = 0;
						else if (k == 1)
							ktemp = 2 * (group + i);
						else if (k == 2)
							ktemp = 1 * (group + i);
						else
							ktemp = 5 * (group + i);
					} else if (delay_count == 16) {
						/* 0,2,1,5 combined with one 0,2,1,5 and three 0,2,1,-1 */
						int	kmap[16] = {0,8,4,20,2,18,10,-6,1,17,9,-7,5,21,13,-3};
						if (kmap[k] >= 0)
							ktemp = kmap[k] * (group + i);
						else
							ktemp = bigN + kmap[k] * (group + i);
					} else {			/* delay_count == 10 or 14 */
						/* Multipliers for the radix-20 or radix-28 step */
						ktemp = k * (group + i);
					}

/* We now calculate the standard sin/cos twiddle factors (temp*0,1,2,3,4,5,6,7 roots of N) */
/* combined with the group all-complex premultiplier roots of minus 1 (same as roots of FFTLEN*2) */
/* combined with the delayed group multipliers. */

					temp = (group + i) * (bigN / N);
					gwsincos1plus01234567by2 (actemp + ktemp, temp, bigN, table); /* premult,delay and temp*0-7 */
					temp = (group + i + upper_sse2_word) * (bigN / N);
					if (gwdata->ALL_COMPLEX_FFT) actemp += upper_sse2_word;
					gwsincos1plus01234567by2 (actemp + ktemp, temp, bigN, table+1);
					table += 32;
				}
			}
			pass1_size /= 8;
		}

/* Output the sin/cos/premultiplier values for the radix-4 block that does the */
/* last 2 levels in pass 1. */

		else {
			N = N * 4;

/* Output the extra sin/cos values needed for the eight_reals FFT work done */
/* on the last pass 1 level.  We double N because the real part of the FFT */
/* is one level behind the complex part of the FFT. */

			if (!gwdata->ALL_COMPLEX_FFT) {
				for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
					temp = (group + i);
					gwsincos15by2 (temp, N*2, table);
					gwsincos15by2 (temp + upper_sse2_word, N*2, table+1);
					table += 8;
				}
			}

/* Output the sin/cos values for the complex delay groups -- specifically the r4_sg4cl_eight_reals_four_complex_fft4 macro. */

			for (k = 0; k < delay_count; k++) {
				for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
					unsigned long bigN, ktemp, actemp;

/* If this is an all-complex FFT, the roots of minus 1 (same as roots of FFTLEN*2) are */
/* split to reduce memory requirements.  We apply the part of the all-complex premultiplier here. */

					if (gwdata->ALL_COMPLEX_FFT) {
						bigN = gwdata->FFTLEN * 2;
						actemp = group + i;
					} else {
						bigN = gwdata->FFTLEN;
						actemp = 0;
					}

/* Factor in the delayed part of the sin/cos multiplies from the first 2 levels.  In the first 2 levels */
/* we use a fixed sin/cos table based only on j, leaving the group+i part to be applied here by */
/* creating delay_count table entries.  For an all-complex FFT we multiply by 0,2,1,-1*(group+i) with N = FFTLEN/2. */
/* For an all-real FFT we multiply by 0,2,1,5*(group+i) with N = FFTLEN. */

					if (gwdata->ALL_COMPLEX_FFT && delay_count == 4) {
						if (k == 0)
							ktemp = 0;
						else if (k == 1)
							ktemp = 2 * (group + i) * 4;
						else if (k == 2)
							ktemp = 1 * (group + i) * 4;
						else
							ktemp = bigN - 1 * (group + i) * 4;
					} else if (gwdata->ALL_COMPLEX_FFT) {
						/* 0,2,1,-1 combined with 0,2,1,-1 */
							int	kmap[16] = {0,8,4,-4,2,10,6,-2,1,9,5,-3,-1,7,3,-5};
							if (kmap[k] >= 0)
								ktemp = kmap[k] * (group + i) * 4;
							else
								ktemp = bigN + kmap[k] * (group + i) * 4;
					} else if (delay_count == 4) {
						if (k == 0)
							ktemp = 0;
						else if (k == 1)
							ktemp = 2 * (group + i);
						else if (k == 2)
							ktemp = 1 * (group + i);
						else
							ktemp = 5 * (group + i);
					} else if (delay_count == 16) {
						/* 0,2,1,5 combined with one 0,2,1,5 and three 0,2,1,-1 */
						int	kmap[16] = {0,8,4,20,2,18,10,-6,1,17,9,-7,5,21,13,-3};
						if (kmap[k] >= 0)
							ktemp = kmap[k] * (group + i);
						else
							ktemp = bigN + kmap[k] * (group + i);
					} else {			/* delay_count == 10 or 14 */
						/* Multipliers for the radix-20 or radix-28 step */
						ktemp = k * (group + i);
					}

/* We now calculate the standard sin/cos twiddle factors (temp*0,1,2,3 roots of N) */
/* combined with the group all-complex premultiplier roots of minus 1 (same as roots of FFTLEN*2) */
/* combined with the delayed group multipliers. */

					temp = (group + i) * (bigN / N);
					gwsincos1plus0123by2 (actemp + ktemp, temp, bigN, table); /* premult,delay and temp*0-3 */
					temp = (group + i + upper_sse2_word) * (bigN / N);
					if (gwdata->ALL_COMPLEX_FFT) actemp += upper_sse2_word;
					gwsincos1plus0123by2 (actemp + ktemp, temp, bigN, table+1);
					table += 16;
				}
			}
			pass1_size /= 4;
		}

/* For the r4_four_complex_djbfft building block levels, output the sin/cos values. */

		while ((pass1_size & 3) == 0) {
			N = N * 4;

/* For the r4_eight_reals_four_complex_djbfft building block levels, output the extra */
/* sin/cos values needed for the eight_reals.  The eight_reals doubles N because */
/* the real part of the FFT is one level behind the complex part of the FFT. */

			if (!gwdata->ALL_COMPLEX_FFT) {
				for (j = 0; j < N*2 / 8; j += pass1_increment) {
					for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
						temp = (group + j + i);
						gwsincos15by2 (temp, N*2, table);
						gwsincos15by2 (temp + upper_sse2_word, N*2, table+1);
						table += 8;
					}
				}
			}

/* Output the sin/cos value for the complex sections, used by the r4_four_complex_djbfft macro */

			for (j = 0; j < N / 4; j += pass1_increment) {
				for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
					temp = (group + j + i);
					gwsincos12by2 (temp, N, table);
					gwsincos12by2 (temp + upper_sse2_word, N, table+1);
					table += 8;
				}
			}
			pass1_size /= 4;
		}

/* For the r5_five_complex_djbfft building block levels, output the sin/cos values. */

		while (pass1_size % 5 == 0) {
			N = N * 5;

/* The r5_ten_reals building blocks require extra sin/cos */
/* values.  The ten_reals doubles N because the real part of the */
/* FFT is one level behind the complex part of the FFT. */

			if (!gwdata->ALL_COMPLEX_FFT) {
				for (j = 0; j < N / 5; j += pass1_increment) {
					for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
						temp = (group + j + i);
						gwsincos13by2 (temp, N*2, table);
						gwsincos13by2 (temp + upper_sse2_word, N*2, table+1);
						table += 8;
					}
				}
			}

/* Output the sin/cos data for the complex sections, (the r5_five_complex_djbfft building block). */

			for (j = 0; j < N / 5; j += pass1_increment) {
				for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
					temp = (group + j + i);
					gwsincos12by2 (temp, N, table);
					gwsincos12by2 (temp + upper_sse2_word, N, table+1);
					table += 8;
				}
			}
			pass1_size /= 5;
		}

/* For the r3_three_complex_djbfft building block levels, output the sin/cos values. */

		while (pass1_size % 3 == 0) {
			N = N * 3;

/* The r3_six_reals building blocks require an extra sin/cos */
/* value.  The six_reals doubles N because the real part of the */
/* FFT is one level behind the complex part of the FFT. */

			if (!gwdata->ALL_COMPLEX_FFT) {
				for (j = 0; j < N / 3; j += pass1_increment) {
					for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
						temp = (group + j + i);
						gwsincos1by2 (temp, N*2, table);
						gwsincos1by2 (temp + upper_sse2_word, N*2, table+1);
						table += 4;
					}
				}
			}

/* Output the sin/cos data for the complex sections (used by the r3_three_complex_djbfft building block). */

			for (j = 0; j < N / 3; j += pass1_increment) {
				for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
					temp = (group + j + i);
					gwsincos1by2 (temp, N, table);
					gwsincos1by2 (temp + upper_sse2_word, N, table+1);
					table += 4;
				}
			}
			pass1_size /= 3;
		}
		ASSERTG (pass1_size == 1);
		if (pass1_size != 1) gwdata->GWERROR = GWERROR_INTERNAL + 5;

/* Calculate the size of each group's sin/cos/premult data for pass1_get_next_block */

		if (group == 0)
			gwdata->pass1_premult_block_size = (unsigned long)
				((char *) table - (char *) gwdata->adjusted_pass1_premults) /
				gwdata->PASS1_CACHE_LINES;
	}

/* Return address of the end of the table */

	return (table);
}

/* This routine builds the fixed postmultiplier table used in pass 1 of the */
/* radix-4/8 delayed DJB FFT - called by gwsetup. */

double *r4delay_build_fixed_pass1_table (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	double	*table)		/* Pointer to the table to fill in */
{
	unsigned long pass1_increment, upper_sse2_word, pass1_size, i, j, N;

/* Initialize some needed constants */

	pass1_increment = gwdata->PASS2_SIZE;
	upper_sse2_word = pass1_increment / 2;
	pass1_size = gwdata->FFTLEN / gwdata->PASS2_SIZE; /* Real values in a pass1 section */

/* Real FFTs output one shared set of sin/cos values for the first 20-reals, 28-reals, or 8-reals FFT. */

	if (! gwdata->ALL_COMPLEX_FFT) {
		N = gwdata->FFTLEN;
		if (pass1_size % 20 == 0) {
			for (j = 0; j < N / 20; j += pass1_increment) {
				for (i = 1; i <= 9; i++) {	/* Create 9 twiddle factors */
					gwsincos1by2 (i * j, N, table);
					gwsincos1by2 (i * (j + upper_sse2_word), N, table+1);
					table += 4;
				}
			}
			N = N / 20;
		}
		else if (pass1_size % 28 == 0) {
			for (j = 0; j < N / 28; j += pass1_increment) {
				for (i = 1; i <= 13; i++) {	/* Create 13 twiddle factors */
					gwsincos1by2 (i * j, N, table);
					gwsincos1by2 (i * (j + upper_sse2_word), N, table+1);
					table += 4;
				}
			}
			N = N / 28;
		}
		else {
			for (j = 0; j < N / 8; j += pass1_increment) {
				gwsincos125by2 (j, N, table);
				gwsincos125by2 (j + upper_sse2_word, N, table+1);
				table += 12;
			}
			N = N / 8;
			/* Sometimes we also use a fixed sin/cos table for */
			/* the next FFT levels to further reduce memory usage. */
			if (pass1_size == 512 || pass1_size == 1024 || pass1_size == 1536 ||
			    pass1_size == 2048 || pass1_size == 3072 || pass1_size == 4096) {
				/* Output the sin/cos values for the real data */
				for (j = 0; j < N*2 / 8; j += pass1_increment) {
					gwsincos15by2 (j, N*2, table);
					gwsincos15by2 (j + upper_sse2_word, N*2, table+1);
					table += 8;
				}
				/* Output the sin/cos values for the complex data */
				for (j = 0; j < N / 4; j += pass1_increment) {
					gwsincos12by2 (j, N, table);
					gwsincos12by2 (j + upper_sse2_word, N, table+1);
					table += 8;
				}
				N = N / 4;
			}
		}
	}

/* For all-complex FFTs, also build the fixed roots-of-minus-one table. */
/* Output these roots in the same order they will be used in the first two */
/* levels of pass 1. */

	else {
		N = gwdata->FFTLEN / 2;
		for (j = 0; j < N / 4; j += pass1_increment) {
			gwsincos12by2 (j, N, table);
			gwsincos12by2 (j + upper_sse2_word, N, table+1);
			table += 8;
		}
		N = N / 4;
		/* Sometimes we also use a fixed sin/cos table for */
		/* the next FFT levels to further reduce memory usage. */
		if (pass1_size == 1536 || pass1_size == 2048 || pass1_size == 2560 ||
		    pass1_size == 3072 || pass1_size == 4096 || pass1_size == 5120) {
			/* Output the sin/cos values for the complex data */
			for (j = 0; j < N / 4; j += pass1_increment) {
				gwsincos12by2 (j, N, table);
				gwsincos12by2 (j + upper_sse2_word, N, table+1);
				table += 8;
			}
			N = N / 4;
		}
	}

/* Return address of the end of the table */

	return (table);
}

/* This routine builds the sin/cos table used in pass 1 by the radix-4/8 DJB */
/* FFT with delayed sin/cos multiplies and with partial normalization. */

double *r4dwpn_build_pass1_table (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	double	*table)		/* Pointer to the table to fill in */
{
	struct gwasm_data *asm_data = (struct gwasm_data *) gwdata->asm_data;
	unsigned long pass1_size, pass1_increment, delay_count;
	unsigned long group, i, j, k, N, temp, upper_sse2_word;
	int	pow2_count;

/* Initialize some needed constants */

	pass1_increment = gwdata->PASS2_SIZE;
	upper_sse2_word = pass1_increment / 2;
	pass1_size = gwdata->FFTLEN / gwdata->PASS2_SIZE; /* Real values in a pass1 section */

/* Determine number of delay groups.  In a standard radix-4 FFT, there is only one sin/cos */
/* group in the last pass 1 level.  We reduce our memory usage by using just one fixed sin/cos */
/* table in the first FFT levels and having multiple groups of sin/cos data in the last pass 1 level. */
/* I call these groups of sin/cos data in the last pass 1 level "delay groups". */

	if (pass1_size % 7 == 0)
		delay_count = 14;
	else if (pass1_size % 5 == 0 && !gwdata->ALL_COMPLEX_FFT)
		delay_count = 10;
 	else if ((pass1_size == 512 && !gwdata->ALL_COMPLEX_FFT) ||
		 (pass1_size == 1024 && !gwdata->ALL_COMPLEX_FFT) ||
		 (pass1_size == 2560 && gwdata->ALL_COMPLEX_FFT) ||
		 (pass1_size == 5120 && gwdata->ALL_COMPLEX_FFT) ||
		 pass1_size == 1536 || pass1_size == 2048 || pass1_size == 3072 || pass1_size == 4096)
		delay_count = 16;
	else
		delay_count = 4;

/* Count the power-of-two FFT levels after the initial FFT levels.  If odd, the */
/* the last levels will be a radix-8, if even the last levels will be a radix-4. */

	pass1_size /= (delay_count * 2);
	for (pow2_count = 0; (pass1_size & 1) == 0; pass1_size /= 2) pow2_count++;

/* Set count of pass 1 blocks share one set of two-to-phi grp multipliers */

	if (pow2_count & 1) gwdata->wpn_count = 8;
	else if (gwdata->FFTLEN / gwdata->PASS2_SIZE == 1792 ||
		 gwdata->FFTLEN / gwdata->PASS2_SIZE == 2048) gwdata->wpn_count = 16;
	else gwdata->wpn_count = 4;
	asm_data->count2 = gwdata->wpn_count;
	asm_data->count3 = asm_data->addcount1 / gwdata->wpn_count;

/* Set pointer to table of multipliers */

	gwdata->adjusted_pass1_premults = table;

/* Loop through all the pass 1 groups in the same order the assembly code will */
/* process the groups. */

	for (group = 0; group < upper_sse2_word; group += gwdata->PASS1_CACHE_LINES) {

		pass1_size = gwdata->FFTLEN / gwdata->PASS2_SIZE; /* Real values in a pass1 section */
		pass1_size /= (delay_count * 2);	/* Complex values we're generating sin/cos data for */
		N = gwdata->PASS2_SIZE;

/* Output the sin/cos/premultiplier values for the radix-8 block that does the */
/* last 3 levels in pass 1.  NOTE:  We do not need the "j loop" (it would loop */
/* from zero to zero) when generating the sin/cos twiddle factors for the last */
/* levels of pass 1. */

		if (pow2_count & 1) {
			N = N * 8;

/* For the r4_g8cl_sixteen_reals_eight_complex_djbfft building block, output the extra */
/* sin/cos values needed for the sixteen_reals.  The sixteen_reals is done in three parts: */
/* 4 four_reals_fft, 1 eight_reals_fft, and 1 four-complex_fft.  The four_reals_fft needs */
/* four sin/cos twiddles using 2*N, the eight_reals_fft and four_complex_fft can use the same */
/* sin/cos twiddles that will be generated later for the r4r8_g8cl_eight_complex_fft8 macro. */

			if (!gwdata->ALL_COMPLEX_FFT) {
				for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
					temp = (group + i);
					gwsincos1plus0123by2 (temp, pass1_increment, N*2, table);
					gwsincos1plus0123by2 (temp + upper_sse2_word, pass1_increment, N*2, table+1);
					table += 16;
				}
			}

/* Output the sin/cos values for the complex delay groups -- specifically the r8_sg8cl_eight_complex_fft8 macro. */

			for (k = 0; k < delay_count; k++) {
				for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
					unsigned long bigN, ktemp, actemp;

/* If this is an all-complex FFT, the roots of minus 1 (same as roots of FFTLEN*2) are */
/* split to reduce memory requirements.  We apply the part of the all-complex premultiplier here. */

					if (gwdata->ALL_COMPLEX_FFT) {
						bigN = gwdata->FFTLEN * 2;
						actemp = group + i;
					} else {
						bigN = gwdata->FFTLEN;
						actemp = 0;
					}

/* Factor in the delayed part of the sin/cos multiplies from the first 2 levels.  In the first 2 levels */
/* we use a fixed sin/cos table based only on j, leaving the group+i part to be applied here by */
/* creating delay_count table entries.  For an all-complex FFT we multiply by 0,2,1,-1*(group+i) with N = FFTLEN/2. */
/* For an all-real FFT we multiply by 0,2,1,5*(group+i) with N = FFTLEN. */

					if (gwdata->ALL_COMPLEX_FFT && delay_count == 4) {
						if (k == 0)
							ktemp = 0;
						else if (k == 1)
							ktemp = 2 * (group + i) * 4;
						else if (k == 2)
							ktemp = 1 * (group + i) * 4;
						else
							ktemp = bigN - 1 * (group + i) * 4;
					} else if (gwdata->ALL_COMPLEX_FFT) {
						/* 0,2,1,-1 combined with 0,2,1,-1 */
						int	kmap[16] = {0,8,4,-4,2,10,6,-2,1,9,5,-3,-1,7,3,-5};
						if (kmap[k] >= 0)
							ktemp = kmap[k] * (group + i) * 4;
						else
							ktemp = bigN + kmap[k] * (group + i) * 4;
					} else if (delay_count == 4) {
						if (k == 0)
							ktemp = 0;
						else if (k == 1)
							ktemp = 2 * (group + i);
						else if (k == 2)
							ktemp = 1 * (group + i);
						else
							ktemp = 5 * (group + i);
					} else if (delay_count == 16) {
						/* 0,2,1,5 combined with one 0,2,1,5 and three 0,2,1,-1 */
						int	kmap[16] = {0,8,4,20,2,18,10,-6,1,17,9,-7,5,21,13,-3};
						if (kmap[k] >= 0)
							ktemp = kmap[k] * (group + i);
						else
							ktemp = bigN + kmap[k] * (group + i);
					} else {			/* delay_count == 10 or 14 */
						/* Multipliers for the radix-20 or radix-28 step */
						ktemp = k * (group + i);
					}

/* We now calculate the standard sin/cos twiddle factors (temp*0,1,2,3,4,5,6,7 roots of N) */
/* combined with the group all-complex premultiplier roots of minus 1 (same as roots of FFTLEN*2) */
/* combined with the delayed group multipliers. */

					temp = (group + i) * (bigN / N);
					gwsincos1plus01234567by2 (actemp + ktemp, temp, bigN, table); /* premult,delay and temp*0-7 */
					temp = (group + i + upper_sse2_word) * (bigN / N);
					if (gwdata->ALL_COMPLEX_FFT) actemp += upper_sse2_word;
					gwsincos1plus01234567by2 (actemp + ktemp, temp, bigN, table+1); /* premult,delay and temp*0-7 */
					table += 32;
				}
			}
			pass1_size /= 8;
		}

/* Output the sin/cos/premultiplier values for the radix-4 block that does the */
/* last 2 levels in pass 1. */

		else {
			N = N * 4;

/* Output the extra sin/cos values needed for the eight_reals FFT work done */
/* on the last pass 1 level.  We double N because the real part of the FFT */
/* is one level behind the complex part of the FFT. */

			if (!gwdata->ALL_COMPLEX_FFT) {
				for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
					temp = (group + i);
					gwsincos15by2 (temp, N*2, table);
					gwsincos15by2 (temp + upper_sse2_word, N*2, table+1);
					table += 8;
				}
			}

/* Output the sin/cos values for the complex delay groups -- specifically the r4_sg4cl_eight_reals_four_complex_fft4 macro. */

			for (k = 0; k < delay_count; k++) {
				for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
					unsigned long bigN, ktemp, actemp;

/* If this is an all-complex FFT, the roots of minus 1 (same as roots of FFTLEN*2) are */
/* split to reduce memory requirements.  We apply the part of the all-complex premultiplier here. */

					if (gwdata->ALL_COMPLEX_FFT) {
						bigN = gwdata->FFTLEN * 2;
						actemp = group + i;
					} else {
						bigN = gwdata->FFTLEN;
						actemp = 0;
					}

/* Factor in the delayed part of the sin/cos multiplies from the first 2 levels.  In the first 2 levels */
/* we use a fixed sin/cos table based only on j, leaving the group+i part to be applied here by */
/* creating delay_count table entries.  For an all-complex FFT we multiply by 0,2,1,-1*(group+i) with N = FFTLEN/2. */
/* For an all-real FFT we multiply by 0,2,1,5*(group+i) with N = FFTLEN. */

					if (gwdata->ALL_COMPLEX_FFT && delay_count == 4) {
						if (k == 0)
							ktemp = 0;
						else if (k == 1)
							ktemp = 2 * (group + i) * 4;
						else if (k == 2)
							ktemp = 1 * (group + i) * 4;
						else
							ktemp = bigN - 1 * (group + i) * 4;
					} else if (gwdata->ALL_COMPLEX_FFT) {
						/* 0,2,1,-1 combined with 0,2,1,-1 */
							int	kmap[16] = {0,8,4,-4,2,10,6,-2,1,9,5,-3,-1,7,3,-5};
							if (kmap[k] >= 0)
								ktemp = kmap[k] * (group + i) * 4;
							else
								ktemp = bigN + kmap[k] * (group + i) * 4;
					} else if (delay_count == 4) {
						if (k == 0)
							ktemp = 0;
						else if (k == 1)
							ktemp = 2 * (group + i);
						else if (k == 2)
							ktemp = 1 * (group + i);
						else
							ktemp = 5 * (group + i);
					} else if (delay_count == 16) {
						/* 0,2,1,5 combined with one 0,2,1,5 and three 0,2,1,-1 */
						int	kmap[16] = {0,8,4,20,2,18,10,-6,1,17,9,-7,5,21,13,-3};
						if (kmap[k] >= 0)
							ktemp = kmap[k] * (group + i);
						else
							ktemp = bigN + kmap[k] * (group + i);
					} else {			/* delay_count == 10 or 14 */
						/* Multipliers for the radix-20 or radix-28 step */
						ktemp = k * (group + i);
					}

/* We now calculate the standard sin/cos twiddle factors (temp*0,1,2,3 roots of N) */
/* combined with the group all-complex premultiplier roots of minus 1 (same as roots of FFTLEN*2) */
/* combined with the delayed group multipliers. */

					temp = (group + i) * (bigN / N);
					gwsincos1plus0123by2 (actemp + ktemp, temp, bigN, table); /* premult,delay and temp*0-3 */
					temp = (group + i + upper_sse2_word) * (bigN / N);
					if (gwdata->ALL_COMPLEX_FFT) actemp += upper_sse2_word;
					gwsincos1plus0123by2 (actemp + ktemp, temp, bigN, table+1);
					table += 16;
				}
			}
			pass1_size /= 4;
		}

/* Output multipliers for the four complex building blocks. */

		while ((pass1_size & 3) == 0) {

			N = N * 4;

/* For the non-wpn levels, output the sin/cos values. */

			if (N != gwdata->PASS2_SIZE * gwdata->wpn_count * 4) {

/* For the r4_eight_reals_four_complex_djbfft building block levels, output the extra */
/* sin/cos values needed for the eight_reals.  The eight_reals doubles N because */
/* the real part of the FFT is one level behind the complex part of the FFT. */

				if (!gwdata->ALL_COMPLEX_FFT) {
					for (j = 0; j < N*2 / 8; j += pass1_increment) {
					    for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
						temp = (group + j + i);
						gwsincos15by2 (temp, N*2, table);
						gwsincos15by2 (temp + upper_sse2_word, N*2, table+1);
						table += 8;
					    }
					}
				}

/* Output the sin/cos value for the complex sections, used by the r4_four_complex_djbfft macro */

				for (j = 0; j < N / 4; j += pass1_increment) {
				    for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
					temp = (group + j + i);
					gwsincos12by2 (temp, N, table);
					gwsincos12by2 (temp + upper_sse2_word, N, table+1);
					table += 8;
				    }
				}
			}

/* For the wpn building block level, output the sin/cos values. */

			else {

/* For the r4_wpn_eight_reals_four_complex_djbfft building block levels, output the extra */
/* sin/cos values needed for the eight_reals.  The eight_reals doubles N because */
/* the real part of the FFT is one level behind the complex part of the FFT. */

				if (!gwdata->ALL_COMPLEX_FFT) {
					for (j = 0; j < N*2 / 8; j += pass1_increment) {
					    for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
						temp = (group + j + i);
						gwsincos15by2_weighted (gwdata->dd_data, temp, upper_sse2_word, N*2, temp, table);
						table += 12;
					    }
					}
				}

/* Output the sin/cos value for the complex sections, used by the r4_wpn_four_complex_djbfft macro */
/* We apply the two-to-phi weight for the upper SSE2 word in the group multipliers.  There is a */
/* reason for doing it there rather than here (it reduces the number of valid fudge factor combinations */
/* for each SSE2 word from 4 to 3). */

				for (j = 0; j < N / 4; j += pass1_increment) {
				    for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
					temp = (group + j + i);
					gwsincos012by2_weighted (gwdata->dd_data, temp, upper_sse2_word, N, temp, table);
					table += 16;
				    }
				}
			}

			pass1_size /= 4;
		}

/* For the r5_five_complex_djbfft building block levels, output the sin/cos values. */

		while (pass1_size % 5 == 0) {
			N = N * 5;

/* The r5_ten_reals building blocks require extra sin/cos */
/* values.  The ten_reals doubles N because the real part of the */
/* FFT is one level behind the complex part of the FFT. */

			if (!gwdata->ALL_COMPLEX_FFT) {
				for (j = 0; j < N / 5; j += pass1_increment) {
					for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
						temp = (group + j + i);
						gwsincos13by2 (temp, N*2, table);
						gwsincos13by2 (temp + upper_sse2_word, N*2, table+1);
						table += 8;
					}
				}
			}

/* Output the sin/cos data for the complex sections, (the r5_five_complex_djbfft building block). */

			for (j = 0; j < N / 5; j += pass1_increment) {
				for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
					temp = (group + j + i);
					gwsincos12by2 (temp, N, table);
					gwsincos12by2 (temp + upper_sse2_word, N, table+1);
					table += 8;
				}
			}
			pass1_size /= 5;
		}

/* For the r3_three_complex_djbfft building block levels, output the sin/cos values. */

		while (pass1_size % 3 == 0) {
			N = N * 3;

/* The r3_six_reals building blocks require an extra sin/cos */
/* value.  The six_reals doubles N because the real part of the */
/* FFT is one level behind the complex part of the FFT. */

			if (!gwdata->ALL_COMPLEX_FFT) {
				for (j = 0; j < N / 3; j += pass1_increment) {
					for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
						temp = (group + j + i);
						gwsincos1by2 (temp, N*2, table);
						gwsincos1by2 (temp + upper_sse2_word, N*2, table+1);
						table += 4;
					}
				}
			}

/* Output the sin/cos data for the complex sections (used by the r3_three_complex_djbfft building block). */

			for (j = 0; j < N / 3; j += pass1_increment) {
				for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
					temp = (group + j + i);
					gwsincos1by2 (temp, N, table);
					gwsincos1by2 (temp + upper_sse2_word, N, table+1);
					table += 4;
				}
			}
			pass1_size /= 3;
		}

/* Calculate the size of each group's sin/cos/premult data for pass1_get_next_block */

		if (group == 0)
			gwdata->pass1_premult_block_size = (unsigned long)
				((char *) table - (char *) gwdata->adjusted_pass1_premults) /
				gwdata->PASS1_CACHE_LINES;
	}

/* Return address of the end of the table */

	return (table);
}

/* This routine builds a normalization table - used by radix-4 normalization */
/* routines.  It differs from the home-grown SSE2 routines in that the different */
/* memory layout for PFA makes this routine much simpler. */

double *r4_build_norm_table (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	double	*table,		/* Pointer to the table to fill in */
	int	col)		/* TRUE if building column, not group, table */
{
	unsigned long i, k, num_cols;

/* Handle one-pass FFTs first, there are no group multipliers */

	if (gwdata->PASS2_SIZE == 0) {
		if (!col) return (table);

/* Loop to build table */

		for (i = 0; i < gwdata->FFTLEN; i++) {
			unsigned long j, table_entry;
			double	ttp, ttmp;

/* Call double-precision routine to compute the two multipliers */

			gwfft_weights3 (gwdata->dd_data, i, &ttp, NULL, &ttmp);

/* Find where this data appears in the FFT array and in the table we are building. */

			j = addr_offset (gwdata, i) / sizeof (double);
			table_entry = j >> 1;

/* Now set the entry for the MSW or LSW in an SSE2 pair */

			table[table_entry*4+(j&1)] = ttmp;
			table[table_entry*4+2+(j&1)] = ttp;
		}
		return (table + gwdata->FFTLEN + gwdata->FFTLEN);
	}

/* Two pass FFTs are handled here */

	num_cols = gwdata->PASS2_SIZE / 2;
	if (col) {

/* Loop to build table */

		for (i = 0; i < num_cols; i++) {
			double	ttp, ttmp;

/* Call double-precision routine to compute the two multipliers */

			gwfft_weights3 (gwdata->dd_data, i, &ttp, NULL, &ttmp);

/* Now set the entry for BOTH the MSW and LSW in an SSE2 pair */

			table[i*4] = ttmp;
			table[i*4+1] = ttmp;
			table[i*4+2] = ttp;
			table[i*4+3] = ttp;
		}
		return (table + num_cols * 4);
	}

/* Build the group multipliers table */

	else {
		unsigned long h, hlimit, haddin, m, mmult, u, umult;

/* Loop to build table */

		umult = gwdata->FFTLEN / 2;
		hlimit = gwdata->FFTLEN / 4 / (2*num_cols);
		for (h = 0; h < hlimit; h++) {
			haddin = h * 2 * num_cols;
			mmult = gwdata->FFTLEN / 4;
			for (u = 0; u < 2; u++) {
				for (m = 0; m < 2; m++) {
					for (k = 0; k < 2; k++) {
						double	ttp, ttmp;
						long	n;

/* Call double-precision routine to compute the two multipliers */

						n = haddin + u * umult + m * mmult + k * num_cols;
						gwfft_weights3 (gwdata->dd_data, n, &ttp, &ttmp, NULL);

/* Now set the entry for BOTH the MSW and LSW in an SSE2 pair */

						table[k] = ttmp;
						table[2+k] = ttp;
					}
					table += 4;
				}
			}
		}
		return (table);
	}
}

/* This routine builds a normalization table - used by radix-4 with partial */
/* normalization FFTs. */

double *r4dwpn_build_norm_table (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	double	*table)		/* Pointer to the table to fill in */
{
	unsigned long i, j, num_cols, upper_sse2_word, second_word_in_cache_line;

/* Build the group multipliers table */

	upper_sse2_word = gwdata->PASS2_SIZE / 2;
	second_word_in_cache_line = gwdata->FFTLEN / 4;
	num_cols = gwdata->PASS2_SIZE * gwdata->wpn_count;
	for (i = 0; i < gwdata->FFTLEN / 4; i += num_cols) {
		for (j = 0; j < gwdata->FFTLEN; j += gwdata->FFTLEN / 2) {
			double	word1_lower_weight, word1_upper_weight, word2_lower_weight, word2_upper_weight;
			int	word1_lower_sort, word1_upper_sort, word2_lower_sort, word2_upper_sort;
			double	ttp, ttmp, ttp_over_b, ttmp_times_b, *tab10;

/* The sort order of the weights determines which fudge factor combination never occur. */

			word1_lower_weight = gwfft_weight_exponent (gwdata->dd_data, i + j);
			word1_upper_weight = gwfft_weight_exponent (gwdata->dd_data, i + j + upper_sse2_word);
			word2_lower_weight = gwfft_weight_exponent (gwdata->dd_data, i + j + second_word_in_cache_line);
			word2_upper_weight = gwfft_weight_exponent (gwdata->dd_data, i + j + second_word_in_cache_line + upper_sse2_word);

/* Now sort the four weights */

			word1_lower_sort = (word1_lower_weight > word1_upper_weight) +
					   (word1_lower_weight > word2_lower_weight) +
					   (word1_lower_weight > word2_upper_weight);
			word1_upper_sort = (word1_upper_weight > word1_lower_weight) +
					   (word1_upper_weight > word2_lower_weight) +
					   (word1_upper_weight > word2_upper_weight);
			word2_lower_sort = (word2_lower_weight > word1_lower_weight) +
					   (word2_lower_weight > word1_upper_weight) +
					   (word2_lower_weight > word2_upper_weight);
			word2_upper_sort = (word2_upper_weight > word1_lower_weight) +
					   (word2_upper_weight > word1_upper_weight) +
					   (word2_upper_weight > word2_lower_weight);

/* Call double-precision routine to compute first set of multipliers. */
/* We compute two-to-phi multiplied by the fudge factor so that normalize won't have to. */

			gwfft_weights_fudged (gwdata->dd_data, i + j, gwdata->b, &ttp, &ttmp, &ttp_over_b, &ttmp_times_b);

/* Set the LSW entries in an SSE2 pair */

			table[0] = ttmp;
			table[2] = word1_lower_sort < 3 ? ttmp : ttmp_times_b;
			table[4] = word1_lower_sort < 2 ? ttmp : ttmp_times_b;
			table[6] = word1_lower_sort < 1 ? ttmp : ttmp_times_b;
			table[8] = ttmp_times_b;
			tab10 = table + 10;
			tab10[0] = ttp;
			tab10[2] = word1_lower_sort < 3 ? ttp : ttp_over_b;
			tab10[4] = word1_lower_sort < 2 ? ttp : ttp_over_b;
			tab10[6] = word1_lower_sort < 1 ? ttp : ttp_over_b;
			tab10[8] = ttp_over_b;

/* Set the MSW entries in an SSE2 pair */

			gwfft_weights_fudged (gwdata->dd_data, i + j + upper_sse2_word, gwdata->b, &ttp, &ttmp, &ttp_over_b, &ttmp_times_b);

			table[1] = ttmp;
			table[3] = word1_upper_sort < 3 ? ttmp : ttmp_times_b;
			table[5] = word1_upper_sort < 2 ? ttmp : ttmp_times_b;
			table[7] = word1_upper_sort < 1 ? ttmp : ttmp_times_b;
			table[9] = ttmp_times_b;
			tab10[1] = ttp;
			tab10[3] = word1_upper_sort < 3 ? ttp : ttp_over_b;
			tab10[5] = word1_upper_sort < 2 ? ttp : ttp_over_b;
			tab10[7] = word1_upper_sort < 1 ? ttp : ttp_over_b;
			tab10[9] = ttp_over_b;

			table += 20;

/* Repeat for the second word in the cache line */

			gwfft_weights_fudged (gwdata->dd_data, i + j + second_word_in_cache_line, gwdata->b, &ttp, &ttmp, &ttp_over_b, &ttmp_times_b);

			table[0] = ttmp;
			table[2] = word2_lower_sort < 3 ? ttmp : ttmp_times_b;
			table[4] = word2_lower_sort < 2 ? ttmp : ttmp_times_b;
			table[6] = word2_lower_sort < 1 ? ttmp : ttmp_times_b;
			table[8] = ttmp_times_b;
			tab10 = table + 10;
			tab10[0] = ttp;
			tab10[2] = word2_lower_sort < 3 ? ttp : ttp_over_b;
			tab10[4] = word2_lower_sort < 2 ? ttp : ttp_over_b;
			tab10[6] = word2_lower_sort < 1 ? ttp : ttp_over_b;
			tab10[8] = ttp_over_b;

			gwfft_weights_fudged (gwdata->dd_data, i + j + second_word_in_cache_line + upper_sse2_word, gwdata->b, &ttp, &ttmp, &ttp_over_b, &ttmp_times_b);

			table[1] = ttmp;
			table[3] = word2_upper_sort < 3 ? ttmp : ttmp_times_b;
			table[5] = word2_upper_sort < 2 ? ttmp : ttmp_times_b;
			table[7] = word2_upper_sort < 1 ? ttmp : ttmp_times_b;
			table[9] = ttmp_times_b;
			tab10[1] = ttp;
			tab10[3] = word2_upper_sort < 3 ? ttp : ttp_over_b;
			tab10[5] = word2_upper_sort < 2 ? ttp : ttp_over_b;
			tab10[7] = word2_upper_sort < 1 ? ttp : ttp_over_b;
			tab10[9] = ttp_over_b;

			table += 20;
		}
	}
	return (table);
}

/* This routine builds a big/little flags table - used by radix-4 normalization */
/* routines.  It differs from the home-grown SSE2 routines in that the different */
/* memory layout for PFA makes this routine much simpler. */

double *r4_build_biglit_table (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	double	*table)		/* Pointer to the table to fill in */
{
	unsigned char *p;
	unsigned long h, i, j, k, m, u, gap;
	unsigned long hlimit, haddin, mmult, umult;

/* Handle one pass FFTs differently */

	if (gwdata->PASS2_SIZE == 0) {

/* Loop to build table */

		p = (unsigned char *) table;
		for (i = 0; i < gwdata->FFTLEN; i++) {
			unsigned long table_entry;

/* Find where this data appears in the FFT array and in the table we are building. */

			j = addr_offset (gwdata, i) / sizeof (double);
			table_entry = j >> 1;

/* Now set the biglit table entry for a LSW in an SSE2 pair */

			if ((j & 1) == 0) {
				p[table_entry] = is_big_word (gwdata, i) * 16;
			}

/* Otherwise, set the biglit table entry for a MSW in an SSE2 pair */

			else {
				if (is_big_word (gwdata, i)) p[table_entry] += 32;
			}
		}
		return ((double *) (p + gwdata->FFTLEN / 2));
	}

/* Determine the gap between XMM high and low words */

	gap = gwdata->PASS2_SIZE / 2;

/* Loop to build table in exactly the same order that it will be */
/* used by the assembly code. */

	p = (unsigned char *) table;
	umult = gwdata->FFTLEN / 2;
	hlimit = gwdata->FFTLEN / 4 / (2*gap);
	for (i = 0; i < gap; i += gwdata->PASS1_CACHE_LINES) {
		for (h = 0; h < hlimit; h++) {
			haddin = h * 2 * gap;
			mmult = gwdata->FFTLEN / 4;
			for (j = 0; j < gwdata->PASS1_CACHE_LINES; j++) {
				for (u = 0; u < 2; u++) {
					for (m = 0; m < 2; m++) {
						for (k = 0; k < 2 * gap; k += gap) {
							unsigned long word;

/* Now set the big/little flag for a LSW in an SSE2 pair */
/* Otherwise, set the big/little flag for a MSW in an SSE2 pair */

							word = haddin + i + j + u * umult + m * mmult + k;
							if (k == 0) *p = is_big_word (gwdata, word) * 16;
							else if (is_big_word (gwdata, word)) *p += 32;

/* Set the ttp and ttmp fudge flags for two pass FFTs.  The fudge flag is */
/* set if the col mult * the grp mult is twice the correct fft_weight, */
/* meaning a mul by 0.5 is required to generate the correct multiplier. */
/* Since we can't do equality compares on floats, this test is a little bit */
/* cryptic. */

							if (gwfft_weight_exponent (gwdata->dd_data, word) + 0.5 <
							    gwfft_weight_exponent (gwdata->dd_data, word % gap) +
							    gwfft_weight_exponent (gwdata->dd_data, word - word % gap)) {
								if (k == 0) *p += 64;
								else *p += 128;
							}

/* Set some offsets that help the assembly code step through the big/lit */
/* array in a non-traditional order.  Two pass-FFTs step through the array */
/* in chunks of PASS1_CACHE_LINES, but the add, sub, and carry propagation */
/* code need to access the big/lit array linearly.  Set two global variables */
/* that tell the assembly code the big/lit array distance between words */
/* 0 and 2, and words 0 and 4. */

							if (word == 2)
								((struct gwasm_data *) gwdata->asm_data)->BIGLIT_INCR2 =
									(uint32_t) ((char *) p - (char *) table);
							if (word == 4)
								((struct gwasm_data *) gwdata->asm_data)->BIGLIT_INCR4 =
									(uint32_t) ((char *) p - (char *) table);
						}
						p++;
					}
				}
			}
		}
	}
	return ((double *) p);
}

/* This routine builds the big/little flags table for a r4dwpn (radix-4 with partial normalization) FFT */

double *r4dwpn_build_biglit_table (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	double	*table)		/* Pointer to the table to fill in */
{
const	int	START_OF_CHAIN = 0x8000;
const	int	END_OF_CHAIN = 0x4000;
const	int	CHAIN_3_COMMON = 0x2000;
const	int	CHAIN_2_COMMON = 0x1000;
const	int	CHAIN_1_COMMON = 0x0800;
	int	combos[16], next_combo[16];
	unsigned char combo_recorded[256];
	unsigned char *p;
	unsigned long upper_sse2_word, fftlen_over_4, num_cols;
	unsigned long group, i, j, k, m, n, num_combos;

/* Big/lit flags form a very regular pattern.  For example, if there are 18.3 b's */
/* per FFT word then you get either a big word followed by two little words or a */
/* big  followed by three little words.  Here we determine which patterns of big/lit */
/* are possible in an xnorm_wpn macro which processes 4 SSE2 words.  There are */
/* 16 possible valid combinations which can be represented by indexing into an array */
/* of 32 SSE2 values. */

/* Generate all possible valid combinations of big/lit flags */

	memset (combo_recorded, 0, sizeof (combo_recorded));
	upper_sse2_word = gwdata->PASS2_SIZE / 2;
	fftlen_over_4 = gwdata->FFTLEN / 4;
	num_combos = 0;

	/* Loop over all cache lines */
	for (i = 0; i < upper_sse2_word; i++) {
		for (j = 0; j < fftlen_over_4; j += 2*upper_sse2_word) {
			int	combo;

			/* Generate combo for this cache line */
			combo = (is_big_word (gwdata, i + j + upper_sse2_word) << 7) +
				(is_big_word (gwdata, i + j) << 6) +
				(is_big_word (gwdata, i + j + fftlen_over_4 + upper_sse2_word) << 5) +
				(is_big_word (gwdata, i + j + fftlen_over_4) << 4) +
				(is_big_word (gwdata, i + j + 2*fftlen_over_4 + upper_sse2_word) << 3) +
				(is_big_word (gwdata, i + j + 2*fftlen_over_4) << 2) +
				(is_big_word (gwdata, i + j + 3*fftlen_over_4 + upper_sse2_word) << 1) +
				(is_big_word (gwdata, i + j + 3*fftlen_over_4));

			/* Ignore this combo if it is a duplicate.  Otherwise, add it to our combos collection. */
			if (! combo_recorded[combo]) {
				combo_recorded[combo] = 1;
				combos[num_combos++] = combo;
			}
		}
	}

/* Concatentate combos to save space.  Let's hope they fit in 48 entries. */

	/* Init the next-in-chain array */
	for (i = 0; i < num_combos; i++)
		next_combo[i] = START_OF_CHAIN + END_OF_CHAIN;

	/* Look for 2 chains where the end of one chain has elements in common with */
	/* the start of the other chain,  First look for 3 elements in common, then */
	/* two, then one. */
	for (n = 3; n != 0; n--) {

		/* Examine all chain starts */
		for (i = 0; i < num_combos; i++) {
			int	chain_end;

			/* Skip if not the start of a chain */
			if (! (next_combo[i] & START_OF_CHAIN)) continue;

			/* Find end of chain */
			for (chain_end = i; ! (next_combo[chain_end] & END_OF_CHAIN); chain_end = next_combo[chain_end] & 0xFF);

			/* Now look at all chain ends */
			for (j = 0; j < num_combos; j++) {

				/* Skip if not a chain end */
				if (! (next_combo[j] & END_OF_CHAIN)) continue;
				/* Can't chain to ourselves! */
				if (j == chain_end) continue;

				/* See if chain end has the proper number of common elements with the chain start */
				if (n == 3 && (combos[i] >> 2) == (combos[j] & 0x3F)) {
					next_combo[j] = (next_combo[j] & START_OF_CHAIN) + CHAIN_3_COMMON + i;
					next_combo[i] &= ~START_OF_CHAIN;
					break;
				}
				if (n == 2 && (combos[i] >> 4) == (combos[j] & 0x0F)) {
					next_combo[j] = (next_combo[j] & START_OF_CHAIN) + CHAIN_2_COMMON + i;
					next_combo[i] &= ~START_OF_CHAIN;
					break;
				}
				if (n == 1 && (combos[i] >> 6) == (combos[j] & 0x03)) {
					next_combo[j] = (next_combo[j] & START_OF_CHAIN) + CHAIN_1_COMMON + i;
					next_combo[i] &= ~START_OF_CHAIN;
					break;
				}
			}
		}
	}

/* HACK: Remember the chains in ASM_TIMERS so that we can properly build LIMIT_INVERSE, */
/* LIMIT_BIGMAX, and LIMIT_BIGMAX_NEG at a later time. */

	n = 0;
	for (i = 0; i < num_combos; i++) {
		if (! (next_combo[i] & START_OF_CHAIN)) continue;

		/* Output up to 4 bytes for each entry in the chain */
		for (j = i; ; j = next_combo[j] & 0xFF) {
			((char *)gwdata->ASM_TIMERS)[n++] = (combos[j] >> 6);
			if (! (next_combo[j] & CHAIN_3_COMMON))
				((char *)gwdata->ASM_TIMERS)[n++] = (combos[j] >> 4) & 3;
			if (! (next_combo[j] & (CHAIN_3_COMMON + CHAIN_2_COMMON)))
				((char *)gwdata->ASM_TIMERS)[n++] = (combos[j] >> 2) & 3;
			if (! (next_combo[j] & (CHAIN_3_COMMON + CHAIN_2_COMMON + CHAIN_1_COMMON)))
				((char *)gwdata->ASM_TIMERS)[n++] = combos[j] & 3;
			if (next_combo[j] & END_OF_CHAIN) break;
		}
	}
	GWASSERT (n <= 48);		// I've seen n get as high as 35
	if (n > 48) gwdata->GWERROR = GWERROR_INTERNAL + 6;

/* Determine the number of column two-to-phi multipliers */

	num_cols = gwdata->PASS2_SIZE * gwdata->wpn_count;

/* Loop to build table in exactly the same order that it will be */
/* used by the assembly code. */

	p = (unsigned char *) table;
	for (group = 0; group < upper_sse2_word; group += gwdata->PASS1_CACHE_LINES) {
	    for (n = 0; n < gwdata->FFTLEN / 4; n += num_cols) {
	        for (j = 0; j < num_cols; j += gwdata->PASS2_SIZE) {
		    for (i = 0; i < gwdata->PASS1_CACHE_LINES; i++) {
			for (k = 0; k < gwdata->FFTLEN; k += gwdata->FFTLEN / 4) {
			    for (m = 0; m < gwdata->PASS2_SIZE; m += upper_sse2_word) {
				unsigned long word;

/* Now set the big/little flag for a LSW in an SSE2 pair */
/* Otherwise, set the big/little flag for a MSW in an SSE2 pair */

				word = group + j + n + i + k + m;
				if (m == 0) *p = is_big_word (gwdata, word) * 16;
				else if (is_big_word (gwdata, word)) *p += 32;

/* Set the ttp and ttmp fudge flags for two pass FFTs.  The fudge flag is */
/* set if the col mult * the grp mult is b times the correct fft_weight, */
/* meaning a mul by 1/b is required to generate the correct multiplier. */
/* Since we can't do equality compares on floats, this test is a little bit */
/* cryptic. */

				if (gwfft_weight_exponent (gwdata->dd_data, word) + 0.5 <
				    gwfft_weight_exponent (gwdata->dd_data, n + k + m) +
				    gwfft_weight_exponent (gwdata->dd_data, group + j + i)) {
					if (m == 0) *p += 64;
					else *p += 128;
				}

/* Set some offsets that help the assembly code step through the big/lit */
/* array in a non-traditional order.  Two pass-FFTs step through the array */
/* in chunks of PASS1_CACHE_LINES, but the add, sub, and carry propagation */
/* code need to access the big/lit array linearly.  Set two global variables */
/* that tell the assembly code the big/lit array distance between words */
/* 0 and 2, and words 0 and 4. */

				if (word == 2)
					((struct gwasm_data *) gwdata->asm_data)->BIGLIT_INCR2 =
						(uint32_t) ((char *) p - (char *) table);
				if (word == 4)
					((struct gwasm_data *) gwdata->asm_data)->BIGLIT_INCR4 =
						(uint32_t) ((char *) p - (char *) table);
			    }

/* Apply our method for reducing fudge factor data by 25%.  We've observed that one */
/* of the 4 combinations never occurs. */

			    if (*p >= 128) *p -= 64;

			    p++;
			}

/* Combine last four big/lit 2-bit flag values into one 6-bit flags value. */

			p -= 4;
			for (m = 0; m <= 44; m++) {
				if (((char *)gwdata->ASM_TIMERS)[m]   == (p[0] & 48) >> 4 &&
				    ((char *)gwdata->ASM_TIMERS)[m+1] == (p[1] & 48) >> 4 &&
				    ((char *)gwdata->ASM_TIMERS)[m+2] == (p[2] & 48) >> 4 &&
				    ((char *)gwdata->ASM_TIMERS)[m+3] == (p[3] & 48) >> 4)
					break;
			}
			ASSERTG (m != 45);
			if (m == 45) gwdata->GWERROR = GWERROR_INTERNAL + 7;

/* Combine 1st and 2nd fudge factor flags.  Combine 3rd and 4th fudge factor flags. */
/* Output the 2 combined fudge factors and the 6-bit big/lit flags value using just 2 bytes. */
/* This is a more compact encoding than in previous versions. */

			p[0] = ((p[1] >> 6) + (p[0] >> 6)) * 16 + ((p[3] >> 6) + (p[2] >> 6)) * 2;
			p[1] = (unsigned char) (m << 2);
			p += 2;
		    }
		}
	    }
	}
	return ((double *) p);
}

/* This routine builds a sin/cos table - used by gwsetup */

double *hg_build_sin_cos_table (
	double	*table,		/* Pointer to the table to fill in */
	unsigned long N,	/* Number of DATA values processed by this */
				/* FFT level.  This explains the divide by 2 */
				/* for complex FFTs later in this routine */
	int	hermetian_skip,	/* True if some sin/cos values are skipped */
	int	type)		/* 0 = old style - a plain old array */
				/* 1 = SSE2 - data is duplicated */
				/* 2 = SSE2 - data is interleaved */
{
	unsigned long i;

/* Handle hermetian skip when interleaving.  First data slot is left */
/* undefined. */

	if (type == 2 && hermetian_skip) type = 3;

/* Special case the really small sin/cos tables.  If N is between 9 and 16 */
/* or between 33 and 64, then the assembly code is only doing one FFT level. */
/* In this case, the code just uses the middle sin/cos values of a 2N sized */
/* table.  We could optimize this inefficient memory usage at a later date. */

	if (N <= 8) return (table);
	if (N >= 9 && N <= 16) N = N * 2;
	if (N >= 33 && N <= 64 && type == 1 && hermetian_skip) N = N * 2;

/* In the all-complex case. build the same size table as the hermetian */
/* case which skips half the i values. */

	if (!hermetian_skip) N = N / 2;

/* Loop to build table. */

	for (i = hermetian_skip ? ((N & 4) ? 4 : 8) : 0; i < N; i += 4) {
		unsigned long shifted_i, shifted_N, flipped_i;
		double	sincos[6];

/* Flip the bits in i.  Our prime-factor-FFT makes this a little complex. */
/* The algorithm below works, but I've long since forgotten why. */

		shifted_i = i; shifted_N = N; flipped_i = 0;
		while ((shifted_N & 1) == 0) {
			flipped_i <<= 1;
			if (shifted_i & 1) flipped_i++;
			shifted_i >>= 1;
			shifted_N >>= 1;
		}
		flipped_i = (flipped_i * shifted_N) + shifted_i;

/* When the FFT is working on real data Hermetian symettry allows us to */
/* eliminate half of the FFT data and consequently half of the sin/cos data */
/* Case 1:  If shifted source is > shifted N/2, then we */
/* do not need these sin/cos values. */
/* Case 2:  If shifted source is zero, loop to find the top */
/* two bits.  Skip the number if the top two bits equal 3. */

		if (hermetian_skip) {
			if (shifted_i > shifted_N / 2) continue;
			if (shifted_i == 0) {
				unsigned long j;
				for (j = i; j > 3; j >>= 1);
				if (j == 3) continue;
			}
		}

/* Compute the 3 sin/cos values */

		gwsincos3 (flipped_i, N, (double *) &sincos);

/* Copy the sin/cos values in the appropriate way */

		if (type == 0) {
			memcpy (table, sincos, sizeof (sincos));
			table += 6;
		} else if (type == 1) {
			table[0] = table[1] = sincos[0];
			table[2] = table[3] = sincos[1];
			table[4] = table[5] = sincos[2];
			table[6] = table[7] = sincos[3];
			table[8] = table[9] = sincos[4];
			table[10] = table[11] = sincos[5];
			table += 12;
		} else if (type == 2) {
			table[0] = sincos[0];
			table[2] = sincos[1];
			table[4] = sincos[2];
			table[6] = sincos[3];
			table[8] = sincos[4];
			table[10] = sincos[5];
			type++;
		} else {
			table[1] = sincos[0];
			table[3] = sincos[1];
			table[5] = sincos[2];
			table[7] = sincos[3];
			table[9] = sincos[4];
			table[11] = sincos[5];
			type--;
			table += 12;
		}
	}
	return (table);
}

/* This routine builds a pass 2 premultiplier table - used by gwsetup */

double *hg_build_premult_table (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	double	*table)		/* Pointer to the table to fill in */
{
	unsigned long i, N, incr, type, m_fudge;

/* Build a premultiplier table for the second pass incrementing by */
/* the pre-calculated pass2_size. */

	N = gwdata->FFTLEN;
	incr = gwdata->PASS2_SIZE;
	if (gwdata->ALL_COMPLEX_FFT) N = N / 2;

/* Mod 2^N+1 arithmetic starts at first data set, */
/* mod 2^N-1 skips some data sets */

 	if (gwdata->ALL_COMPLEX_FFT) i = 0;
	else i = incr * 4;

/* To add in the flipped_m component, we want the sin/cos of flipped_m */
/* over pass2_size.  This fudge factor will convert flipped_m into something */
/* that can be divided by N. */

	m_fudge = N / gwdata->PASS2_SIZE;

/* Loop to build table. */

	type = 0;
	for ( ; i < N; i += incr) {
		unsigned long shifted_i, shifted_N, flipped_i, k, l, m;
		unsigned long grouping_size;
		double	*table_start;
		double	sincos[2];

/* Flip the bits in i.  Our prime-factor-FFT makes this a little complex. */
/* The algorithm below works, but I've long since forgotten why. */

		shifted_i = i; shifted_N = N; flipped_i = 0;
		while ((shifted_N & 1) == 0) {
			flipped_i <<= 1;
			if (shifted_i & 1) flipped_i++;
			shifted_i >>= 1;
			shifted_N >>= 1;
		}
		flipped_i = (flipped_i * shifted_N) + shifted_i;

/* When the FFT is working on real data Hermetian symettry allows us to */
/* eliminate half of the FFT data and consequently half of the sin/cos data */
/* Case 1:  If shifted source is > shifted N/2, then we */
/* do not need these sin/cos values. */
/* Case 2:  If shifted source is zero, loop to find the top */
/* two bits.  Skip the number if the top two bits equal 3. */

		if (!gwdata->ALL_COMPLEX_FFT) {
			if (shifted_i > shifted_N / 2) continue;
			if (shifted_i == 0) {
				unsigned long j;
				for (j = i; j > 3; j >>= 1);	
				if (j == 3) continue;
			}
		}

/* Generate the group multipliers.  We used to always create groups of 4, */
/* but to save memory we now group by different amounts based on pass 2 size */

		grouping_size = 4;
		if (gwdata->PASS2_SIZE == 1024) grouping_size = 8;
		if (gwdata->PASS2_SIZE == 2048) grouping_size = 8;
		if (gwdata->PASS2_SIZE == 4096) grouping_size = 16;
		if (gwdata->PASS2_SIZE == 8192) grouping_size = 16;
		table_start = table;
		for (k = 0; k < incr / 4; k += grouping_size) {

/* There are 4 multipliers in a XMM_PMD set */

			for (l = 0; l < 4; l++) {
				unsigned long real_k, pm;

/* Compute the sin/cos value (root of unity) */

				real_k = l * incr/4 + k;
				pm = real_k * flipped_i;
				if (!gwdata->ALL_COMPLEX_FFT) {
					gwsincos (pm % N, N, (double *) &sincos);
				}

/* If C > 0, then also multiply by the proper root of -1.  This is done */
/* by changing the value we are taking the sin/cos of */

				else {
					pm = pm * 4 + real_k;
					gwsincos (pm % (N*4), N*4, (double *) &sincos);
				}

/* Save the premultiplier value */

				table[l*4+type] = sincos[0];
				table[l*4+2+type] = sincos[1];
			}
			table += 16;
		}
	
/* Generate the 16 column multipliers * first 4 sin/cos values. */
/* Also multiply by the LAST 4 sin/cos values so that the xsincos_complex */
/* table can be 1/4 of it's usual size.  The extra room in the cache more */
/* than compensates for the 12 extra column multipliers. */

		for (m = 0; m < 4; m++) {
		    unsigned long flipped_m;
		    flipped_m = ((m & 1) << 1) + ((m & 2) >> 1);
		    for (k = 0; k < grouping_size; k++) {
			for (l = 0; l < 4; l++) {
				unsigned long pm;

/* Compute the sin/cos value (root of unity) */

				pm = k * flipped_i +
				     l * flipped_m * N/16 +
				     (k & 3) * flipped_m * m_fudge;
				if (!gwdata->ALL_COMPLEX_FFT) {
					gwsincos (pm % N, N, (double *) &sincos);
				}

/* If C > 0, then also multiply by the proper root of -1.  This is done */
/* by changing the value we are taking the sin/cos of */

				else {
					pm = pm * 4 + k;
					gwsincos (pm % (N*4), N*4, (double *) &sincos);
				}

/* Save the premultiplier value */

				table[l*4+type] = sincos[0];
				table[l*4+2+type] = sincos[1];
			}
			table += 16;
		    }
		}

		if (type == 0) table = table_start;
		type = 1 - type;
 	}

	return (table);
}

/* This routine builds a plus 1 premultiplier table - used by gwsetup */
/* when c is positive. */

double *hg_build_plus1_table (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	double	*table)		/* Pointer to the table to fill in */
{
	unsigned long i, j, k, l, N;
	int	pfa;

/* Set flag if this is a 3*2^n FFT */

	pfa = (gwdata->FFTLEN != pow_two_above_or_equal (gwdata->FFTLEN));

/* Adjust for two-pass FFTs */

	if (gwdata->PASS2_SIZE == 0) N = gwdata->FFTLEN;
	else N = gwdata->FFTLEN / (gwdata->PASS2_SIZE / 2);

/* Loop to build premultiplier table in the same order as the underlying */
/* assembly macro needs them.  The pfa macro operates on 3 cache lines */
/* while the power-of-two macro operates on 2 cache lines. */
/* A 64 length FFT needs 0,8,16,24 for the macro then 3 more iterations */
/* for the cache lines beginning with 2,4,6. */
/* A 48 length FFT needs 0,8,16 and 4,12,20 for the first macro then */
/* one more iteration for the cache lines beginning with 2. */

	for (i = 0; i < N / (pfa ? 24 : 32); i++) {
	for (l = 0; l < 2; l++) {
		double	sincos[2];

/* Generate the pre multipliers (roots of -1). */

		for (k = 0; k < (unsigned long) (pfa ? 3 : 4); k++) {
		for (j = 0; j < 2; j++) {
			long	temp;

/* Compute the sin/cos value */

			if (pfa)
				temp = (long) ((i * 2 + l * N/12 + j + k * N/6) % N);
			else
				temp = (long) ((i * 4 + l * 2 + j + k * N/8) % N);
			gwsincos (temp, N*2, (double *) &sincos);

/* Save the premultiplier value */

			table[0+j] = sincos[0];
			table[2+j] = sincos[1];

/* For two-pass FFTs we could apply the root of -1 for the upper SSE2 */
/* double here or in the pass 2 premultipliers.  We've arbitrarily chosen */
/* to do it in the pass 2 premults. */

			if (gwdata->PASS2_SIZE) {
				j = 1;
				table[0+j] = sincos[0];
				table[2+j] = sincos[1];
			}
		}
		table += 4;
		}
	}
 	}

	return (table);
}

/* This routine builds a normalization table - used by SSE2 normalization */
/* routines */

double *hg_build_norm_table (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	double	*table,		/* Pointer to the table to fill in */
	int	col)		/* TRUE if building column, not group, table */
{
	unsigned long i, k, num_cols;

/* Handle one-pass FFTs first, there are no group multipliers */

	if (gwdata->PASS2_SIZE == 0) {
		if (!col) return (table);

/* Loop to build table */

		for (i = 0; i < gwdata->FFTLEN; i++) {
			unsigned long j, table_entry;
			double	ttp, ttmp;

/* Call double-precision routine to compute the two multipliers */

			gwfft_weights3 (gwdata->dd_data, i, &ttp, NULL, &ttmp);

/* Find where this data appears in the FFT array and in the table we are building. */

			j = addr_offset (gwdata, i) / sizeof (double);
			table_entry = j >> 1;

/* Now set the entry for the MSW or LSW in an SSE2 pair */

			table[table_entry*4+(j&1)] = ttmp;
			table[table_entry*4+2+(j&1)] = ttp;
		}
		return (table + gwdata->FFTLEN + gwdata->FFTLEN);
	}

/* Two pass FFTs are handled here */

	num_cols = gwdata->PASS2_SIZE / 2;
	if (col) {

/* Loop to build table */

		for (i = 0; i < num_cols; i++) {
			double	ttp, ttmp;

/* Call double-precision routine to compute the two multipliers */

			gwfft_weights3 (gwdata->dd_data, i, &ttp, NULL, &ttmp);

/* Now set the entry for BOTH the MSW and LSW in an SSE2 pair */

			table[i*4] = ttmp;
			table[i*4+1] = ttmp;
			table[i*4+2] = ttp;
			table[i*4+3] = ttp;
		}
		return (table + num_cols * 4);
	}

/* Build the group multipliers table */

	else {
		unsigned long pfa, h, hlimit, haddin, m, mmult, u, umult;

/* Determine if this is a PFA 5, 6, 7, or 8 */

		for (pfa = gwdata->FFTLEN; pfa > 8; pfa >>= 1);

/* Loop to build table */

		umult = gwdata->FFTLEN / 2;
		hlimit = gwdata->FFTLEN / 4 / (2*num_cols);
		for (h = 0; h < hlimit; h++) {
			if (pfa == 5) {
				if (h < hlimit / 5) {
					haddin = h * 2 * num_cols;
					mmult = gwdata->FFTLEN / 20;
				} else {
					haddin = gwdata->FFTLEN/10 + (h - hlimit/5) * 2 * num_cols;
					mmult = gwdata->FFTLEN / 5;
				}
			} else if (pfa == 7) {
				if (h < hlimit / 7) {
					haddin = h * 2 * num_cols;
					mmult = gwdata->FFTLEN / 28;
				} else if (h < 3 * hlimit / 7) {
					haddin = gwdata->FFTLEN/14 + (h - hlimit/7) * 2 * num_cols;
					mmult = gwdata->FFTLEN / 14;
				} else {
					haddin = 3*gwdata->FFTLEN/14 + (h - 3*hlimit/7) * 2 * num_cols;
					mmult = gwdata->FFTLEN / 7;
				}
			} else {
				haddin = h * 2 * num_cols;
				mmult = gwdata->FFTLEN / 4;
			}
			for (u = 0; u < 2; u++) {
			for (m = 0; m < 2; m++) {
			for (k = 0; k < 2; k++) {
				double	ttp, ttmp;
				long	n;

/* Call double-precision routine to compute the two multipliers */

				n = haddin + u * umult + m * mmult + k * num_cols;
				gwfft_weights3 (gwdata->dd_data, n, &ttp, &ttmp, NULL);

/* Now set the entry for BOTH the MSW and LSW in an SSE2 pair */

				table[k] = ttmp;
				table[2+k] = ttp;
			}
			table += 4;
			}
			}
		}
		return (table);
	}
}

/* This routine builds a big/little flags table - used by SSE2 normalization */
/* routines */

double *hg_build_biglit_table (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	double	*table)		/* Pointer to the table to fill in */
{
	unsigned char *p;
	unsigned long h, i, j, k, m, u, gap;
	unsigned long pfa, hlimit, haddin, mmult, umult;

/* Handle one pass FFTs differently */

	if (gwdata->PASS2_SIZE == 0) {

/* Loop to build table */

		p = (unsigned char *) table;
		for (i = 0; i < gwdata->FFTLEN; i++) {
			unsigned long table_entry;

/* Find where this data appears in the FFT array and in the table we are building. */

			j = addr_offset (gwdata, i) / sizeof (double);
			table_entry = j >> 1;

/* Now set the biglit table entry for a LSW in an SSE2 pair */

			if ((j & 1) == 0) {
				p[table_entry] = is_big_word (gwdata, i) * 16;
			}

/* Otherwise, set the biglit table entry for a MSW in an SSE2 pair */

			else {
				if (is_big_word (gwdata, i)) p[table_entry] += 32;
			}
		}
		return ((double *) (p + gwdata->FFTLEN / 2));
	}

/* Determine if this is a PFA 5, 6, 7, or 8 */

	for (pfa = gwdata->FFTLEN; pfa > 8; pfa >>= 1);

/* Determine the gap between XMM high and low words */

	gap = gwdata->PASS2_SIZE / 2;

/* Loop to build table in exactly the same order that it will be */
/* used by the assembly code.  This is especially ugly in the PFA cases */

	p = (unsigned char *) table;
	umult = gwdata->FFTLEN / 2;
	hlimit = gwdata->FFTLEN / 4 / (2*gap);
	for (i = 0; i < gap; i += gwdata->PASS1_CACHE_LINES) {
	for (h = 0; h < hlimit; h++) {
		if (pfa == 5) {
			if (h < hlimit / 5) {
				haddin = h * 2 * gap;
				mmult = gwdata->FFTLEN / 20;
			} else {
				haddin = gwdata->FFTLEN/10 + (h - hlimit/5) * 2 * gap;
				mmult = gwdata->FFTLEN / 5;
			}
		} else if (pfa == 7) {
			if (h < hlimit / 7) {
				haddin = h * 2 * gap;
				mmult = gwdata->FFTLEN / 28;
			} else if (h < 3 * hlimit / 7) {
				haddin = gwdata->FFTLEN/14 + (h - hlimit/7) * 2 * gap;
				mmult = gwdata->FFTLEN / 14;
			} else {
				haddin = 3*gwdata->FFTLEN/14 + (h - 3*hlimit/7) * 2 * gap;
				mmult = gwdata->FFTLEN / 7;
			}
		} else {
			haddin = h * 2 * gap;
			mmult = gwdata->FFTLEN / 4;
		}
	for (j = 0; j < gwdata->PASS1_CACHE_LINES; j++) {
	for (u = 0; u < 2; u++) {
	for (m = 0; m < 2; m++) {
	for (k = 0; k < 2 * gap; k += gap) {
		unsigned long word;

/* Now set the big/little flag for a LSW in an SSE2 pair */
/* Otherwise, set the big/little flag for a MSW in an SSE2 pair */

		word = haddin + i + j + u * umult + m * mmult + k;
		if (k == 0) *p = is_big_word (gwdata, word) * 16;
		else if (is_big_word (gwdata, word)) *p += 32;

/* Set the ttp and ttmp fudge flags for two pass FFTs.  The fudge flag is */
/* set if the col mult * the grp mult is twice the correct fft_weight, */
/* meaning a mul by 0.5 is required to generate the correct multiplier. */
/* Since we can't do equality compares on floats, this test is a little bit */
/* cryptic. */

		if (gwfft_weight_exponent (gwdata->dd_data, word) + 0.5 <
		    gwfft_weight_exponent (gwdata->dd_data, word % gap) +
		    gwfft_weight_exponent (gwdata->dd_data, word - word % gap)) {
			if (k == 0) *p += 64;
			else *p += 128;
		}

/* Set some offsets that help the assembly code step through the big/lit */
/* array in a non-traditional order.  Two pass-FFTs step through the array */
/* in chunks of PASS1_CACHE_LINES, but the add, sub, and carry propagation */
/* code need to access the big/lit array linearly.  Set two global variables */
/* that tell the assembly code the big/lit array distance between words */
/* 0 and 2, and words 0 and 4. */

		if (word == 2)
			((struct gwasm_data *) gwdata->asm_data)->BIGLIT_INCR2 =
				(uint32_t) ((char *) p - (char *) table);
		if (word == 4)
			((struct gwasm_data *) gwdata->asm_data)->BIGLIT_INCR4 =
				(uint32_t) ((char *) p - (char *) table);
	}
	p++;
	}
	}
	}
	}
	}
	return ((double *) p);
}

/* Ancient x87 setup routines */

#ifndef X86_64

/* This routine builds an x87 sin/cos table - used by gwsetup */

double *x87_build_sin_cos_table (
	double	*table,		/* Pointer to the table to fill in */
	unsigned long N,
	int	hermetian_skip)	/* True if some sin/cos values are skipped */
{
	unsigned long i;

/* Special case the really small sin/cos tables.  If N is between 9 and 16 */
/* then the assembly code is only doing one FFT level. */
/* In this case, the code just uses the middle sin/cos values of a 2N sized */
/* table.  We could optimize this inefficient memory usage at a later date. */

	if (N <= 8) return (table);
	if (N >= 9 && N <= 16) N = N * 2;

/* The N value passed in represents the number of real numbers that are */
/* processed in a section.  If heremetian_skip is not set, then we are */
/* instead dealing with complex numbers and there are half as many complex */
/* numbers in a section.  For example, when doing 8 levels in pass 2, this */
/* routine is called with N=512.  The first real section has 512 values, */
/* while the remaining pass 2 sections have 256 complex values. */

	if (!hermetian_skip) N = N / 2;

/* Loop to build table */

	for (i = hermetian_skip ? ((N & 4) ? 4 : 8) : 0; i < N; i += 4) {
		unsigned long shifted_i, shifted_N, flipped_i;
		double	sincos[6];

/* Flip the bits in i.  Our prime-factor-FFT makes this a little complex. */
/* The algorithm below works, but I've long since forgotten why. */

		shifted_i = i; shifted_N = N; flipped_i = 0;
		while ((shifted_N & 1) == 0) {
			flipped_i <<= 1;
			if (shifted_i & 1) flipped_i++;
			shifted_i >>= 1;
			shifted_N >>= 1;
		}
		flipped_i = (flipped_i * shifted_N) + shifted_i;

/* When the FFT is working on real data Hermetian symettry allows us to */
/* eliminate half of the FFT data and consequently half of the sin/cos data */
/* Case 1:  If shifted source is > shifted N/2, then we */
/* do not need these sin/cos values. */
/* Case 2:  If shifted source is zero, loop to find the top */
/* two bits.  Skip the number if the top two bits equal 3. */

		if (hermetian_skip) {
			if (shifted_i > shifted_N / 2) continue;
			if (shifted_i == 0) {
				unsigned long j;
				for (j = i; j > 3; j >>= 1);
				if (j == 3) continue;
			}
		}

/* Compute the 3 sin/cos values */

		gwsincos3 (flipped_i, N, (double *) &sincos);

/* Copy the sin/cos values to the table */

		memcpy (table, sincos, sizeof (sincos));
		table += 6;
	}
	return (table);
}

/* This routine builds a pass 2 premultiplier table - used by gwsetup */

double *x87_build_premult_table (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	double	*table)		/* Pointer to the table to fill in */
{
	unsigned long i, N, incr;

/* Build a premultiplier table for the second pass incrementing by */
/* the pre-calculated pass2_size. */

	N = gwdata->FFTLEN;
	incr = gwdata->PASS2_SIZE;
	if (gwdata->ALL_COMPLEX_FFT) N = N / 2;

/* Mod 2^N+1 arithmetic starts at first data set, */
/* mod 2^N-1 skips some data sets */

	if (gwdata->ALL_COMPLEX_FFT) i = 0;
	else i = incr * 2;

/* Loop to build table */

	for ( ; i < N; i += incr) {
		unsigned long shifted_i, shifted_N, flipped_i, k, l;
		double	sincos[2];

/* Flip the bits in i.  Our prime-factor-FFT makes this a little complex. */
/* The algorithm below works, but I've long since forgotten why. */

		shifted_i = i; shifted_N = N; flipped_i = 0;
		while ((shifted_N & 1) == 0) {
			flipped_i <<= 1;
			if (shifted_i & 1) flipped_i++;
			shifted_i >>= 1;
			shifted_N >>= 1;
		}
		flipped_i = (flipped_i * shifted_N) + shifted_i;

/* When the FFT is working on real data Hermetian symettry allows us to */
/* eliminate half of the FFT data and consequently half of the sin/cos data */
/* Case 1:  If shifted source is > shifted N/2, then we */
/* do not need these sin/cos values. */
/* Case 2:  If shifted source is zero, loop to find the top */
/* two bits.  Skip the number if the top two bits equal 3. */

		if (!gwdata->ALL_COMPLEX_FFT) {
			if (shifted_i > shifted_N / 2) continue;
			if (shifted_i == 0) {
				unsigned long j;
				for (j = i; j > 3; j >>= 1);	
				if (j == 3) continue;
			}
		}

/* Generate the group multipliers */

		for (k = 0; k < incr / 4; k += 4) {

/* There are 4 multipliers in a PMD set */

			for (l = 0; l < 4; l++) {

/* Compute the sin/cos value (root of unity) */

				if (!gwdata->ALL_COMPLEX_FFT) {
					gwsincos (((l * incr/4 + k) * flipped_i) % N, N, (double *) &sincos);
				}

/* If C > 0, then also multiply by the proper root of -1.  This is done */
/* by changing the value we are taking the sin/cos of */

				else {
					gwsincos (((l * incr/4 + k) * flipped_i * 4 + l*incr/4+k) % (N*4), N*4, (double *) &sincos);
				}

/* Save the premultiplier values */

				table[l*2] = sincos[0];
				table[l*2+1] = sincos[1];
			}
			table += 8;
		}
	
/* Generate the 4 column multipliers * 4 sin/cos values */

		for (k = 0; k < 4; k++) {
			for (l = 0; l < 4; l++) {

/* Compute the sin/cos value (root of unity) */

				if (!gwdata->ALL_COMPLEX_FFT) {
					gwsincos ((k * flipped_i + l * N/16) % N, N, (double *) &sincos);
				}

/* If C > 0, then also multiply by the proper root of -1.  This is done */
/* by changing the value we are taking the sin/cos of */

				else {
					gwsincos (((k * flipped_i * 2 + l * N/8) *2 + k) % (N*4), N*4, (double *) &sincos);
				}

/* Save the premultiplier value */

				table[l*2] = sincos[0];
				table[l*2+1] = sincos[1];
			}
			table += 8;
		}
 	}

	return (table);
}

/* This routine builds a plus 1 premultiplier table - used by gwsetup */
/* when c is positive. */

double *x87_build_plus1_table (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	double	*table)		/* Pointer to the table to fill in */
{
	unsigned long i, k, N;
	int	pfa;

/* Set flag if this is a 3*2^n FFT */

	pfa = (gwdata->FFTLEN != pow_two_above_or_equal (gwdata->FFTLEN));

/* Adjust for two-pass FFTs */

	if (gwdata->PASS2_SIZE == 0) N = gwdata->FFTLEN;
	else N = gwdata->FFTLEN / gwdata->PASS2_SIZE;

/* Loop to build premultiplier table in the same order as the underlying */
/* assembly macro needs them. */

	for (i = 0; i < N / (pfa ? 6 : 8); i++) {
		double	sincos[2];

/* Generate the pre multipliers (roots of -1) used in one three_complex */
/* or four complex macro. */

		for (k = 0; k < (unsigned long) (pfa ? 3 : 4); k++) {
			long	temp;

/* Compute the sin/cos value */

			if (pfa)
				temp = (long) ((i + k * N/6) % N);
			else
				temp = (long) ((i + k * N/8) % N);
			gwsincos (temp, N*2, (double *) &sincos);

/* Save the premultiplier value */

			table[0] = sincos[0];
			table[1] = sincos[1];
			table += 2;
		}
	}

	return (table);
}

/* This routine builds a normalization table - used by x87 normalization */
/* routines */

double *x87_build_norm_table (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	double	*table,		/* Pointer to the table to fill in */
	int	col)		/* TRUE if building column, not group, table */
{
	unsigned long i, k, num_cols;

/* Handle one-pass FFTs first, there are no group multipliers */

	if (gwdata->PASS2_SIZE == 0) {
		if (!col) return (table);

/* Loop to build table */

		for (i = 0; i < gwdata->FFTLEN; i++) {
			unsigned long j;
			double	ttp, ttmp;

/* Call double-precision routine to compute the two multipliers */

			gwfft_weights3 (gwdata->dd_data, i, &ttp, NULL, &ttmp);

/* Find where this data appears in the FFT array and in the table we are building. */

			j = addr_offset (gwdata, i) / sizeof (double);

/* Now set the appropriate table entry.  These are put into the array */
/* in the same order that the normalization code needs them. */

			table[j*2] = ttmp;
			table[j*2+1] = ttp;
		}
		return (table + gwdata->FFTLEN + gwdata->FFTLEN);
	}

/* Two pass FFTs are handled here */

	num_cols = gwdata->PASS2_SIZE;
	if (col) {

/* Loop to build columns table */

		for (i = 0; i < num_cols; i++) {
			double	ttp, ttmp;

/* Call double-precision routine to compute the two multipliers */

			gwfft_weights3 (gwdata->dd_data, i, &ttp, NULL, &ttmp);

/* Now set the appropriate table entry.  These are put into the array */
/* in the same order that the normalization code needs them. */

			table[i+i] = ttmp;
			table[i+i+1] = ttp;
		}
		return (table + num_cols * 2);
	}

/* Build the group multipliers table */

	else {
		unsigned long num_grps;
		
/* Loop to build group table */

		num_grps = gwdata->FFTLEN / num_cols;
		for (i = 0; i < num_grps; i++) {
			double	ttp, ttmp;

/* Call double-precision routine to compute the two multipliers */

			gwfft_weights3 (gwdata->dd_data, i * num_cols, &ttp, &ttmp, NULL);

/* Now set the appropriate table entry.  These are put into the array */
/* in the same order that the normalization code needs them. */

			if (i < num_grps / 2) k = i * 2;
			else k = (i - num_grps / 2) * 2 + 1;
			table[k+k] = ttmp;
			table[k+k+1] = ttp;
		}
		return (table + num_grps * 2);
	}
}

/* This routine builds a big/little flags table - used by x87 normalization */
/* routines */

double *x87_build_biglit_table (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	double	*table)		/* Pointer to the table to fill in */
{
	unsigned char *p;
	unsigned long i, j, k, m;

/* Handle one pass FFTs differently */

	if (gwdata->PASS2_SIZE == 0) {

/* Loop to build table */

		p = (unsigned char *) table;
		for (i = 0; i < gwdata->FFTLEN; i++) {
			unsigned long table_entry;

/* Find where this data appears in the FFT array and in the table we are building. */

			j = addr_offset (gwdata, i) / sizeof (double);
			table_entry = j >> 1;

/* Now set the biglit table entry for a LSW in a pair */

			if ((j & 1) == 0) {
				p[table_entry] = is_big_word (gwdata, i) * 16;
			}

/* Otherwise, set the biglit table entry for a MSW in a pair */

			else {
				if (is_big_word (gwdata, i)) p[table_entry] += 32;
			}
		}
		return ((double *) (p + gwdata->FFTLEN / 2));
	}

/* Loop to build table in exactly the same order that it will be */
/* used by the assembly code. */

	p = (unsigned char *) table;
	for (i = 0; i < gwdata->PASS2_SIZE; i += gwdata->PASS1_CACHE_LINES * 2) {
	for (j = 0; j < gwdata->FFTLEN / 2; j += gwdata->PASS2_SIZE) {
	for (k = 0; k < gwdata->PASS1_CACHE_LINES * 2; k++) {
	for (m = 0; m < gwdata->FFTLEN; m += gwdata->FFTLEN / 2) {
		unsigned long word;

/* Now set the big/little flag for a LSW in a pair */
/* Otherwise, set the big/little flag for a MSW in a pair */

		word = i + j + k + m;
		if (m == 0) *p = is_big_word (gwdata, word) * 16;
		else if (is_big_word (gwdata, word)) *p += 32;

/* Set the ttp and ttmp fudge flags for two pass FFTs */
/* The fudge flag is set if col mult * grp mult will be greater than 2 */

		if (gwfft_weight_exponent (gwdata->dd_data, word) + 0.5 <
		    gwfft_weight_exponent (gwdata->dd_data, word % gwdata->PASS2_SIZE) +
		    gwfft_weight_exponent (gwdata->dd_data, word - word % gwdata->PASS2_SIZE)) {
			if (m == 0) *p += 64;
			else *p += 128;
		}

/* Set some offsets that help the assembly code step through the big/lit */
/* array in a non-traditional order.  Two pass-FFTs step through the array */
/* in chunks of PASS1_CACHE_LINES, but the add, sub, and carry propagation */
/* code need to access the big/lit array linearly.  Set two global variables */
/* that tell the assembly code the big/lit array distance between words */
/* 0 and 2, and words 0 and 4. */

		if (word == 2)
			((struct gwasm_data *) gwdata->asm_data)->BIGLIT_INCR2 =
				(uint32_t) ((char *) p - (char *) table);
		if (word == 4)
			((struct gwasm_data *) gwdata->asm_data)->BIGLIT_INCR4 =
				(uint32_t) ((char *) p - (char *) table);
	}
	p++;
	}
	}
	}
	return ((double *) p);
}

/* End ancient x87 setup routines */

#endif
