//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright(c) 2021 Intel Corporation. All rights reserved.
//
//  Author: Shriram Shastry <malladi.sastry@linux.intel.com>
//

  /* Include Files */
#include "trig.h"
#include <stdint.h>

/* Function Definitions */
/*
 * function [y, x] = sincos_fixpt(theta)
 *
 * Arguments	: const int32_t theta[2047]
 *		  int32_t y[2047]
 *		  int32_t x[2047]
 * Return Type	: void
 */
void sincos_fixpt(const int32_t theta[2047], int32_t y[2047], int32_t x[2047])
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
	int32_t negate;

	for (idx = 0; idx < 2047; idx++) {
		z = theta[idx];
		if (z > 421657428) {
			if ((z - 843314857) <= 421657428) {
				z -= 843314857;
				negate = 1;
			} else {
				z -= 1686629713;
				negate = -1;
			}
		} else if (z < -421657428) {
			if ((z + 843314857) >= -421657428) {
				z += 843314857;
				negate = 1;
			} else {
				z += 1686629713;
				negate = -1;
			}
		} else {
			negate = -1;
		}
		z <<= 2ULL;
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
			xtmp = (xn >> (b_idx + 1));
			ytmp = (b_yn >> (b_idx + 1));
		}
		if (negate) {
			x[idx] = -xn;
			y[idx] = -b_yn;
		} else {
			x[idx] = xn;
			y[idx] = b_yn;
		}
	}


}

