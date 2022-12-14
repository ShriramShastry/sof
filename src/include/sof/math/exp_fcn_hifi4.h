/* SPDX-License-Identifier: BSD-3-Clause
 *
 * Copyright(c) 2022 Intel Corporation. All rights reserved.
 *
 * Author: Shriram Shastry <malladi.sastry@linux.intel.com>
 *
 */
#ifndef __SOF_MATH_EXPONENTIAL_HIFI4_H__
#define __SOF_MATH_EXPONENTIAL_HIFI4_H__

#include <stdint.h>
#include <stddef.h>

#define EXPONENTIAL_GENERIC
#if defined(__XCC__)

#include <xtensa/config/core-isa.h>
#if XCHAL_HAVE_HIFI4
#undef EXPONENTIAL_GENERIC
#define EXPONENTIAL_HIFI4
#endif

#endif

/* define constants */
#define BIT_MASK_LOW_Q27P5 0x0000000008000000
#define BIT_MASK_Q62P2 0x4000000000000000
#define CONVERG_ERROR 28823037624320LL	// error smaller than 1e-4,1/2 ^ -44.7122876209085
#define QUOTIENT_SCALE 0x400000000
#define TERMS_Q23P9 0x800000
#define LSHIFT_BITS 0x2000
#define NUM_Q14P18 0x4000

int32_t exp_int32(int32_t x);

#endif /* __SOF_MATH_EXPONENTIAL_HIFI4_H__ */
