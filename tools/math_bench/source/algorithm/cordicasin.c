//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright(c) 2021 Intel Corporation. All rights reserved.
//
//  Author: Shriram Shastry <malladi.sastry@linux.intel.com>
//

/* Include Files */
#include <stdint.h>
#include "trig.h"

/* Function Definitions */
/*
 * Arguments	: int32_t sinValue
 *		  int16_t nIters
 *		  const int32_t LUT[30]
 * Return Type	: int32_t
 */
int32_t iScalarCordicAsin(int32_t sinValue, int16_t nIters, const int32_t LUT[30])
{
	int32_t b_i;
	int32_t b_sinValue;
	int32_t i;
	int32_t x;
	int32_t xDShift;
	int32_t xShift;
	int32_t y;
	int32_t yDShift;
	int32_t yShift;
	int32_t z;
	int16_t j;
	int16_t k;
	b_sinValue = (sinValue >> 1U);
	if (b_sinValue > 379625062) {
		x = 0;
		y = 1073741824;
		z = 843314856;
	} else {
		x = 1073741824;
		y = 0;
		z = 0;
	}
	b_sinValue <<= 1U;
	i = (int32_t)((int16_t)(nIters - 1));
	for (b_i = 0; b_i < i; b_i++) {
		j = (int16_t)((b_i + 1) << 1U);
		if (j >= 31) {
			j = 31;
		}
		if (b_i < 31) {
			k = (int16_t)b_i;
		} else {
			k = 31;
		}
		xShift = (x >> ((uint32_t)k));
		xDShift = (x >> ((uint32_t)j));
		yShift = (y >> ((uint32_t)k));
		yDShift = (y >> ((uint32_t)j));
		if (y == b_sinValue) {
			x += xDShift;
			y += yDShift;
		} else if (((y >= b_sinValue) && (x >= 0)) || ((y < b_sinValue) && (x < 0))) {
			x = (x - xDShift) + yShift;
			y = (y - yDShift) - xShift;
			z -= LUT[b_i];
		} else {
			x = (x - xDShift) - yShift;
			y = (y - yDShift) + xShift;
			z += LUT[b_i];
		}
		b_sinValue += (b_sinValue >> ((uint32_t)j));
	}
	if (z < 0) {
		z = -z;
	}
	return z;
}

