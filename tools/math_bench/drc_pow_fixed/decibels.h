// SPDX - License - Identifier: BSD - 3 - Clause
//
//Copyright(c) 2021 Intel Corporation.All rights reserved.
//
//Author : Shriram Shastry <malladi.sastry@linux.intel.com>
#ifndef __SOF_MATH_DECIBELS_H__
#define __SOF_MATH_DECIBELS_H__

#include <stdint.h>
#include "/../common/typdef.h"

int32_t exp_fixed(int32_t x); /* Input is Q5.27, output is Q12.20 */

#endif /* __SOF_MATH_DECIBELS_H__ */