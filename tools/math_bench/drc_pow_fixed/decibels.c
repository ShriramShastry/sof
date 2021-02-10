// SPDX - License - Identifier: BSD - 3 - Clause
//
//Copyright(c) 2020 Intel Corporation.All rights reserved.
//
//Author : Shriram Shastry <malladi.sastry@linux.intel.com>
#include <stdint.h>
#include "/../common/typdef.h"
#include "exp_small_fixed.h"

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

/* Fixed point exponent function for approximate range -11.5 .. 7.6
 * that corresponds to decibels range -100 .. +66 dB.
 *
 * The functions uses rule exp(x) = exp(x/2) * exp(x/2) to reduce
 * the input argument for private small value exp() function that is
 * accurate with input range -2.0 .. +2.0. The number of possible
 * divisions by 2 is computed into variable n. The returned value is
 * exp()^(2^n).
 *
 * Input  is Q5.27, -16.0 .. +16.0, but note the input range limitation
 * Output is Q12.20, 0.0 .. +2048.0
 */

int32_t exp_fixed(int32_t x)
{
	int32_t xs;
	int32_t y;
	int32_t y0;
	int i;
	int n = 0;

	if (x < Q_CONVERT_FLOAT(-11.5, 27))
		return 0;

	if (x > Q_CONVERT_FLOAT(7.6245, 27))
		return INT32_MAX;

	/* x is Q5.27 */
	xs = x;
	while (xs >= TWO_Q27 || xs <= MINUS_TWO_Q27) {
		xs >>= 1;
		n++;
	}

	/* exp_small_fixed() input is Q3.29, while x1 is Q5.27
	 * exp_small_fixed() output is Q9.23, while y0 is Q12.20
	 */
	y0 = Q_SHIFT_RND(exp_small_fixed(Q_SHIFT_LEFT(xs, 27, 29)), 23, 20);
	y = ONE_Q20;
	for (i = 0; i < (1 << n); i++)
		y = (int32_t)Q_MULTSR_32X32((int64_t)y, y0, 20, 20, 20);

	return y;
}
