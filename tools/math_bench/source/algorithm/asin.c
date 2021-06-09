//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright(c) 2021 Intel Corporation. All rights reserved.
//
//  Author: Shriram Shastry <malladi.sastry@linux.intel.com>
//

/*  CORDIC-based approximation of inverse cosine */
/*  inverse cosine of x based on a CORDIC approximation */
/*  real scalar, vector, matrix, or N-dimensional array */
/*  containing real values from -1 to 1 input Q2.30 */
/*  returned cordicacos is in radians in Q3.29 */
/*  ARC SINE: 'Error (max = 0.000000026950137), THD+N  = -157.948952635422842' */

/* Include Files */
#include <stdint.h>
#include "trig.h"

/* Function Definitions */
/*
 * function th_rad_fxp = asin_fixpt(a)
 *
 * Arguments	: const int32_t a[21]
 *		  int32_t th_rad_fxp[21]
 * Return Type	: void
 */
void asin_fixpt(const int32_t a[21], int32_t th_rad_fxp[21])
{
	static const int32_t iv[30] = {
	    497837829, 263043836, 133525158, 67021686, 33543515, 16775850, 8388437, 4194282,
	    2097149,   1048575,	  524287,    262143,   131071,	 65535,	   32767,   16383,
	    8191,      4095,	  2047,	     1023,     511,	 255,	   127,	    63,
	    31,	       15,	  8,	     4,	       2,	 1};
	int32_t i;
	int32_t saturatedUnaryMinus;

	for (i = 0; i < 21; i++) {
		saturatedUnaryMinus = a[i];
		if (saturatedUnaryMinus >= 0) {
			th_rad_fxp[i] = iScalarCordicAsin(saturatedUnaryMinus, (int16_t)31, iv);
		} else {
			if (saturatedUnaryMinus <= INT32_MIN) {
				saturatedUnaryMinus = INT32_MAX;
			} else {
				saturatedUnaryMinus = -saturatedUnaryMinus;
			}
			th_rad_fxp[i] = -iScalarCordicAsin(saturatedUnaryMinus, (int16_t)31, iv);
		}
	}
}
