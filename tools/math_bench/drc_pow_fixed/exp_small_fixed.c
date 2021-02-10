// SPDX - License - Identifier: BSD - 3 - Clause
//
//Copyright(c) 2020 Intel Corporation.All rights reserved.
//
//Author : Shriram Shastry <malladi.sastry@linux.intel.com>
#include <stdint.h>
#include "/../common/typdef.h"
/* Exponent function for small values of x. This function calculates
 * fairly accurately exponent for x in range -2.0 .. +2.0. The iteration
 * uses first 11 terms of Taylor series approximation for exponent
 * function. With the current scaling the numerator just remains under
 * 64 bits with the 11 terms.
 *
 * See https://en.wikipedia.org/wiki/Exponential_function#Computation
 *
 * The input is Q3.29
 * The output is Q9.23
 */

static int32_t exp_small_fixed(int32_t x)
{
	int64_t p;
	int64_t num = Q_SHIFT_RND(x, 29, 23);
	int32_t y0 = (int32_t)num;
	int32_t den = 1;
	int32_t inc;
	int k;

	/* Numerator is x^k, denominator is k! */
	for (k = 2; k < 12; k++) {
		p = num * x; /* Q9.23 x Q3.29 -> Q12.52 */
		num = Q_SHIFT_RND(p, 52, 23);
		den = den * k;
		inc = (int32_t)(num / den);
		y0 += inc;
	}

	return y0 + ONE_Q23;
}
