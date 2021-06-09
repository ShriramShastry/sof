//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright(c) 2021 Intel Corporation. All rights reserved.
//
//  Author: Shriram Shastry <malladi.sastry@linux.intel.com>
//

/* Include Files */
#include "trig.h"
#include <stdint.h>

/*
 * function [cdcCoSinTh] = cosine_fixpt(thRadFxp)
 *
 * Description
 *  The cordic cosine algorithm only converges when the angle is in the range
 *  [-pi/2, pi/2).If an angle is outside of this range, then a multiple of
 *  pi/2 is added or subtracted from the angle until it is within the range
 *  [-pi/2,pi/2).
 *
 * Arguments	: const int32_t thRadFxp[2047]
 *		  int32_t cdcCoSinTh[2047]
 * Return Type	: void
 */
void cosine_fixpt(const int32_t thRadFxp[2047], int32_t cdcCoSinTh[2047])
{
	static const int32_t inpLUT[31] = {
	    843314857, 497837829, 263043837, 133525159, 67021687, 33543516, 16775851, 8388437,
	    4194283,   2097149,	  1048576,   524288,	262144,	  131072,   65536,    32768,
	    16384,     8192,	  4096,	     2048,	1024,	  512,	    256,      128,
	    64,	       32,	  16,	     8,		4,	  2,	    1};
	int32_t b_idx;
	int32_t b_yn;
	int32_t idx;
	int32_t xn;
	int32_t xtmp;
	int32_t ytmp;
	int32_t z;
	int32_t sign;

	/*  +------------------------------------+------------------------------------+----------------+-----------+ */
	/*  |	     thRadFxp			 |	     cdcCoSinTh		      |	   thRadFxp    | cdcCoSinTh| */
	/*  +----------+---------------+---------+-----------+--------------+---------+----------------+-----------+ */
	/*  |WordLength| FractionLength| Signbit | WordLength|FractionLength|Signbit  |	    Qformat    |    Qformat| */
	/*  +----------+---------------+---------+-----------+--------------+---------+----------------+-----------+ */
	/*  | 32       |    28	       |    0	 |  32	     | 30	    | 1	      |	    4.28       |     2.30  | */
	/*  +----------+---------------+---------+-----------+--------------+---------+----------------+-----------+ */
	/*  Addition or subtraction by a multiple of pi/2 is done in the data type  */
	/*  of the input. When the fraction length is 29, then the quantization error */
	/*  introduced by the addition or subtraction of pi/2 is done with 29 bits of  */
	/*  precision */
	/*  Input range of cordiccos must be in the range [-2*pi, 2*pi), a signed type */
	/*  with fractionLength = wordLength-4 will fit this range without overflow */
	/*  So increase of fractionLength makes the addition or subtraction of a  */
	/*  multiple of pi/2  more precise */
	for (idx = 0; idx < 2047; idx++) {
		z = thRadFxp[idx];
		if (z > 421657428) {
			if ((z - 843314857) <= 421657428) {
				z -= 843314857;
				sign = 1;
			} else {
				z -= 1686629713;
				sign = -1;
			}
		} else if (z < -421657428) {
			if ((z + 843314857) >= -421657428) {
				z += 843314857;
				sign = 1;
			} else {
				z += 1686629713;
				sign = -1;
			}
		} else {
			sign = -1;
		}
		z *= 4;
		b_yn = 0;
		xn = 652032874;
		xtmp = 652032874;
		ytmp = 0;
		for (b_idx = 0; b_idx < 31; b_idx++) {
			if (z < 0) {
				z += inpLUT[b_idx];
				xn += ytmp;
				b_yn -= xtmp;
			} else {
				z -= inpLUT[b_idx];
				xn -= ytmp;
				b_yn += xtmp;
			}
			xtmp = xn >> (b_idx + 1);
			ytmp = b_yn >> (b_idx + 1);
		}
		if (sign) {
			cdcCoSinTh[idx] = -xn;
		} else {
			cdcCoSinTh[idx] = xn;
		}
	}
}

