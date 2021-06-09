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
/*  'Error (max = 0.000000026950137), THD+N  = -157.948952635422842' */

/* Include Files */
#include <stdint.h>
#include "trig.h"
#include <string.h>

/* Function Definitions */
/*
 * function [a] = init_asin_fixpt()
 *
 * a = fi([-1:.1:1],1,32,30);
 *
 * Arguments	: int32_t a[21]Q2.30
 * Return Type	: void
 */
void init_asin_fixpt(int32_t a[21])
{
	static const int32_t iv[21] = {
	    -1073741824, -966367642, -858993459, -751619277, -644245094, -536870912, -429496730,
	    -322122547,	 -214748365, -107374182, 0,	     107374182,	 214748365,  322122547,
	    429496730,	 536870912,  644245094,	 751619277,  858993459,	 966367642,  1073741824};
	(void)memcpy(&a[0], &iv[0], 21U * (sizeof(int32_t)));
}

