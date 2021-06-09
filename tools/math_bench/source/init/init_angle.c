//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright(c) 2021 Intel Corporation. All rights reserved.
//
//  Author: Shriram Shastry <malladi.sastry@linux.intel.com>
//


  /* Include Files */
#include "trig.h"
#include <string.h>
#include <stdint.h>

/* Function Definitions */
/*
 * function dblRandomVals = init_angle_fixpt()
 *
 * Arguments	: cint32_T dblRandomVals[20]
 * Return Type	: void
 */
void init_angle_fixpt(cint32_T dblRandomVals[20])
{
	static const cint32_T icv[20] = {{
					     547178834, /* re */
					     377742374	/* im */
					 },
					 {
					     481036337, /* re */
					     1550053697 /* im */
					 },
					 {
					     1434089580, /* re */
					     1016833508	 /* im */
					 },
					 {
					     1813335193, /* re */
					     327920753	 /* im */
					 },
					 {
					     739808117, /* re */
					     732506673	/* im */
					 },
					 {
					     1676110988, /* re */
					     1304381568	 /* im */
					 },
					 {
					     1450195708, /* re */
					     411672616	 /* im */
					 },
					 {
					     14388141,	/* re */
					     1585701926 /* im */
					 },
					 {
					     1293214653, /* re */
					     521409030	 /* im */
					 },
					 {
					     830646675, /* re */
					     1970101499 /* im */
					 },
					 {
					     1967095022, /* re */
					     577887850	 /* im */
					 },
					 {
					     2576981,	/* re */
					     1643898733 /* im */
					 },
					 {
					     992996439, /* re */
					     405230165	/* im */
					 },
					 {
					     911177312, /* re */
					     617401549	/* im */
					 },
					 {
					     989775214, /* re */
					     195635761	/* im */
					 },
					 {
					     1653991906, /* re */
					     1237380078	 /* im */
					 },
					 {
					     692563477, /* re */
					     1467590325 /* im */
					 },
					 {
					     1685130419, /* re */
					     1173814562	 /* im */
					 },
					 {
					     1012323792, /* re */
					     914183789	 /* im */
					 },
					 {
					     76879915,	/* re */
					     1383838463 /* im */
					 }};
	(void)memcpy(&dblRandomVals[0], &icv[0], 20U * (sizeof(cint32_T)));
}

