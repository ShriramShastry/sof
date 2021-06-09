//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright(c) 2021 Intel Corporation. All rights reserved.
//
//  Author: Shriram Shastry <malladi.sastry@linux.intel.com>
//

 /*
 *  CORDIC-based four quadrant inverse tangent
 *  Input Arguments
 *  computes the four quadrant arctangent of y and x
 *  which is in the range [-pi, pi] radians
 *  y and x are floating-point numbers,
 *  Output Arguments
 *  theta is the arctangent value,
 *  floating to fixed point error and signal to	 noise ratio
 *  Error = 0.000000001862645149230957
 *  THD+N = -156.296144004064558430400211 (dBc)
 */




/* Include Files */
#include "trig.h"
#include <string.h>
#include <stdint.h>

/* Function Declarations */
static void MultiWordAdd(const uint32_t u1[], const uint32_t u2[], uint32_t y[], int32_t n);

static void MultiWordNeg(const uint32_t u1[], uint32_t y[], int32_t n);

static void MultiWordSignedWrap(const uint32_t u1[], int32_t n1, uint32_t n2, uint32_t y[]);

static void MultiWordSub(const uint32_t u1[], const uint32_t u2[], uint32_t y[], int32_t n);

static void sLong2MultiWord(int32_t u, uint32_t y[], int32_t n);

static void sMultiWord2MultiWord(const uint32_t u1[], int32_t n1, uint32_t y[], int32_t n);

static int32_t sMultiWord2sLongSat(const uint32_t u1[], int32_t n1);

static void sMultiWord2sMultiWordSat(const uint32_t u1[], int32_t n1, uint32_t y[], int32_t n);

static int32_t sMultiWordCmp(const uint32_t u1[], const uint32_t u2[], int32_t n);

static int32_t sMultiWordLt(const uint32_t u1[], const uint32_t u2[], int32_t n);

static void sMultiWordShl(const uint32_t u1[], int32_t n1, uint32_t n2, uint32_t y[], int32_t n);

static void sMultiWordShr(const uint32_t u1[], int32_t n1, uint32_t n2, uint32_t y[], int32_t n);

/* Function Definitions */
/*
 * Arguments	: const uint32_t u1[]
 *		  const uint32_t u2[]
 *		  uint32_t y[]
 *		  int32_t n
 * Return Type	: void
 */
static void MultiWordAdd(const uint32_t u1[], const uint32_t u2[], uint32_t y[], int32_t n)
{
	int32_t carry = 0;
	int32_t i;
	uint32_t u1i;
	uint32_t yi;
	for (i = 0; i < n; i++) {
		u1i = u1[i];
		yi = (u1i + u2[i]) + ((uint32_t)carry);
		y[i] = yi;
		if (((uint32_t)carry) != 0U) {
			carry = (yi <= u1i) ? ((int32_t)1) : ((int32_t)0);
		} else {
			carry = (yi < u1i) ? ((int32_t)1) : ((int32_t)0);
		}
	}
}

/*
 * Arguments	: const uint32_t u1[]
 *		  uint32_t y[]
 *		  int32_t n
 * Return Type	: void
 */
static void MultiWordNeg(const uint32_t u1[], uint32_t y[], int32_t n)
{
	int32_t carry = 1;
	int32_t i;
	uint32_t yi;
	for (i = 0; i < n; i++) {
		yi = (~u1[i]) + ((uint32_t)carry);
		y[i] = yi;
		carry = (yi < ((uint32_t)carry)) ? ((int32_t)1) : ((int32_t)0);
	}
}

/*
 * Arguments	: const uint32_t u1[]
 *		  int32_t n1
 *		  uint32_t n2
 *		  uint32_t y[]
 * Return Type	: void
 */
static void MultiWordSignedWrap(const uint32_t u1[], int32_t n1, uint32_t n2, uint32_t y[])
{
	int32_t n1m1;
	uint32_t mask;
	uint32_t ys;
	n1m1 = n1 - 1;
	if (0 <= (n1m1 - 1)) {
		(void)memcpy(&y[0], &u1[0], ((uint32_t)n1m1) * (sizeof(uint32_t)));
	}
	mask = (((uint32_t)1U) << (31U - n2));
	if ((u1[n1 - 1] & mask) != 0U) {
		ys = INT32_MAX;
	} else {
		ys = 0U;
	}
	mask = (mask << 1U) - 1U;
	y[n1 - 1] = (u1[n1 - 1] & mask) | ((~mask) & ys);
}

/*
 * Arguments	: const uint32_t u1[]
 *		  const uint32_t u2[]
 *		  uint32_t y[]
 *		  int32_t n
 * Return Type	: void
 */
static void MultiWordSub(const uint32_t u1[], const uint32_t u2[], uint32_t y[], int32_t n)
{
	int32_t borrow = 0;
	int32_t i;
	uint32_t u1i;
	uint32_t yi;
	for (i = 0; i < n; i++) {
		u1i = u1[i];
		yi = (u1i - u2[i]) - ((uint32_t)borrow);
		y[i] = yi;
		if (((uint32_t)borrow) != 0U) {
			borrow = (yi >= u1i) ? ((int32_t)1) : ((int32_t)0);
		} else {
			borrow = (yi > u1i) ? ((int32_t)1) : ((int32_t)0);
		}
	}
}

/*
 * Arguments	: int32_t u
 *		  uint32_t y[]
 *		  int32_t n
 * Return Type	: void
 */
static void sLong2MultiWord(int32_t u, uint32_t y[], int32_t n)
{
	int32_t i;
	uint32_t yi;
	y[0] = (uint32_t)u;
	if (u < 0) {
		yi = INT32_MAX;
	} else {
		yi = 0U;
	}
	for (i = 1; i < n; i++) {
		y[i] = yi;
	}
}

/*
 * Arguments	: const uint32_t u1[]
 *		  int32_t n1
 *		  uint32_t y[]
 *		  int32_t n
 * Return Type	: void
 */
static void sMultiWord2MultiWord(const uint32_t u1[], int32_t n1, uint32_t y[], int32_t n)
{
	int32_t i;
	int32_t nm;
	uint32_t u1i;
	if (n1 < n) {
		nm = n1;
	} else {
		nm = n;
	}
	if (0 <= (nm - 1)) {
		(void)memcpy(&y[0], &u1[0], ((uint32_t)nm) * (sizeof(uint32_t)));
	}
	if (n > n1) {
		if ((u1[n1 - 1] & 2147483648U) != 0U) {
			u1i = INT32_MAX;
		} else {
			u1i = 0U;
		}
		for (i = nm; i < n; i++) {
			y[i] = u1i;
		}
	}
}

/*
 * Arguments	: const uint32_t u1[]
 *		  int32_t n1
 * Return Type	: int32_t
 */
static int32_t sMultiWord2sLongSat(const uint32_t u1[], int32_t n1)
{
	uint32_t y;
	sMultiWord2sMultiWordSat(u1, n1, (uint32_t *)(&y), 1);
	return (int32_t)y;
}

/*
 * Arguments	: const uint32_t u1[]
 *		  int32_t n1
 *		  uint32_t y[]
 *		  int32_t n
 * Return Type	: void
 */
static void sMultiWord2sMultiWordSat(const uint32_t u1[], int32_t n1, uint32_t y[], int32_t n)
{
	int32_t i;
	int32_t nm1;
	uint32_t ys;
	int32_t doSaturation = -1;
	nm1 = n - 1;
	if ((u1[n1 - 1] & 2147483648U) != 0U) {
		ys = INT32_MAX;
	} else {
		ys = 0U;
	}
	if (n1 > n) {
		doSaturation = (((u1[n1 - 1] ^ u1[n - 1]) & 2147483648U) != 0U);
		i = n1 - 1;
		while ((!doSaturation) && (i >= n)) {
			doSaturation = (u1[i] != ys);
			i--;
		}
	}
	if (doSaturation) {
		ys = ~ys;
		for (i = 0; i < nm1; i++) {
			y[i] = ys;
		}
		y[i] = ys ^ 2147483648U;
	} else {
		if (n1 < n) {
			nm1 = n1;
		} else {
			nm1 = n;
		}
		if (0 <= (nm1 - 1)) {
			(void)memcpy(&y[0], &u1[0], ((uint32_t)nm1) * (sizeof(uint32_t)));
		}
		for (i = 0; i < nm1; i++) {
		}
		while (i < n) {
			y[i] = ys;
			i++;
		}
	}
}

/*
 * Arguments	: const uint32_t u1[]
 *		  const uint32_t u2[]
 *		  int32_t n
 * Return Type	: int32_t
 */
static int32_t sMultiWordCmp(const uint32_t u1[], const uint32_t u2[], int32_t n)
{
	int32_t i;
	int32_t y;
	uint32_t su1;
	uint32_t u2i;
	su1 = u1[n - 1] & 2147483648U;
	if (su1 != (u2[n - 1] & 2147483648U)) {
		if (su1 != 0U) {
			y = -1;
		} else {
			y = 1;
		}
	} else {
		y = 0;
		i = n;
		while ((y == 0) && (i > 0)) {
			i--;
			su1 = u1[i];
			u2i = u2[i];
			if (su1 != u2i) {
				if (su1 > u2i) {
					y = 1;
				} else {
					y = -1;
				}
			}
		}
	}
	return y;
}

/*
 * Arguments	: const uint32_t u1[]
 *		  const uint32_t u2[]
 *		  int32_t n
 * Return Type	: int32_t
 */
static int32_t sMultiWordLt(const uint32_t u1[], const uint32_t u2[], int32_t n)
{
	return sMultiWordCmp(u1, u2, n) < 0;
}

/*
 * Arguments	: const uint32_t u1[]
 *		  int32_t n1
 *		  uint32_t n2
 *		  uint32_t y[]
 *		  int32_t n
 * Return Type	: void
 */
static void sMultiWordShl(const uint32_t u1[], int32_t n1, uint32_t n2, uint32_t y[], int32_t n)
{
	int32_t i;
	int32_t nb;
	int32_t nc;
	uint32_t nl;
	uint32_t u1i;
	uint32_t yi;
	uint32_t ys;
	nb = ((int32_t)n2) / 32;
	if ((u1[n1 - 1] & 2147483648U) != 0U) {
		ys = INT32_MAX;
	} else {
		ys = 0U;
	}
	if (nb > n) {
		nc = n;
	} else {
		nc = nb;
	}
	u1i = 0U;
	if (0 <= (nc - 1)) {
		(void)memset(&y[0], 0, ((uint32_t)nc) * (sizeof(uint32_t)));
	}
	for (i = 0; i < nc; i++) {
	}
	if (nb < n) {
		nl = n2 - (((uint32_t)nb) * 32U);
		nb += n1;
		if (nb > n) {
			nb = n;
		}
		nb -= i;
		if (nl > 0U) {
			for (nc = 0; nc < nb; nc++) {
				yi = (u1i >> (32U - nl));
				u1i = u1[nc];
				y[i] = yi | (u1i << nl);
				i++;
			}
			if (i < n) {
				y[i] = (u1i >> (32U - nl)) | (ys << nl);
				i++;
			}
		} else {
			for (nc = 0; nc < nb; nc++) {
				y[i] = u1[nc];
				i++;
			}
		}
	}
	while (i < n) {
		y[i] = ys;
		i++;
	}
}

/*
 * Arguments	: const uint32_t u1[]
 *		  int32_t n1
 *		  uint32_t n2
 *		  uint32_t y[]
 *		  int32_t n
 * Return Type	: void
 */
static void sMultiWordShr(const uint32_t u1[], int32_t n1, uint32_t n2, uint32_t y[], int32_t n)
{
	int32_t i;
	int32_t i1;
	int32_t nb;
	int32_t nc;
	uint32_t nr;
	uint32_t u1i;
	uint32_t yi;
	uint32_t ys;
	nb = ((int32_t)n2) / 32;
	i = 0;
	if ((u1[n1 - 1] & 2147483648U) != 0U) {
		ys = INT32_MAX;
	} else {
		ys = 0U;
	}
	if (nb < n1) {
		nc = n + nb;
		if (nc > n1) {
			nc = n1;
		}
		nr = n2 - (((uint32_t)nb) * 32U);
		if (nr > 0U) {
			u1i = u1[nb];
			for (i1 = nb + 1; i1 < nc; i1++) {
				yi = (u1i >> nr);
				u1i = u1[i1];
				y[i] = yi | (u1i << (32U - nr));
				i++;
			}
			if (nc < n1) {
				yi = u1[nc];
			} else {
				yi = ys;
			}
			y[i] = (u1i >> nr) | (yi << (32U - nr));
			i++;
		} else {
			for (i1 = nb; i1 < nc; i1++) {
				y[i] = u1[i1];
				i++;
			}
		}
	}
	while (i < n) {
		y[i] = ys;
		i++;
	}
}

/*
 * function [cdcatan2Th] = atan2_fixpt(thetaY,thetaX)
 *
 * Arguments	: const int32_t thetaY[359]Q.3.29
 *		  const int32_t thetaX[359]Q3.29
 *		  int32_t cdcatan2Th[359]Q3.29
 * Return Type	: void
 */
void atan2_fixpt(const int32_t thetaY[359], const int32_t thetaX[359], int32_t cdcatan2Th[359])
{
	static const int64m_t r3 = {
	    {0U, 0U} /* s_buff */
	};
	static const int64m_t r6 = {
	    {1686629713U, 0U} /* s_buff */
	};
	static const int32_t inpLUT[31] = {
	    421657428, 248918915, 131521918, 66762579, 33510843, 16771758, 8387925, 4194219,
	    2097141,   1048575,	  524288,    262144,   131072,	 65536,	   32768,   16384,
	    8192,      4096,	  2048,	     1024,     512,	 256,	   128,	    64,
	    32,	       16,	  8,	     4,	       2,	 1,	   1};
	int64m_t b_c;
	int64m_t c;
	int64m_t r;
	int64m_t r1;
	int64m_t r10;
	int64m_t r11;
	int64m_t r2;
	int64m_t r4;
	int64m_t r5;
	int64m_t r7;
	int64m_t xAcc;
	int64m_t yAcc;
	int96m_t r8;
	int96m_t r9;
	int32_t b_idx;
	int32_t idx;
	int32_t zN;
	int32_t x_quad_adjust;
	int32_t y_nonzero;
	int32_t y_quad_adjust;
	for (idx = 0; idx < 359; idx++) {
		b_idx = thetaY[idx];
		if (b_idx < 0) {
			sLong2MultiWord(b_idx, (uint32_t *)(&r1.s_buff[0U]), 2);
			MultiWordSignedWrap((uint32_t *)(&r1.s_buff[0U]), 2, 30U,
					    (uint32_t *)(&r2.s_buff[0U]));
			MultiWordNeg((uint32_t *)(&r2.s_buff[0U]), (uint32_t *)(&r.s_buff[0U]), 2);
			MultiWordSignedWrap((uint32_t *)(&r.s_buff[0U]), 2, 30U,
					    (uint32_t *)(&yAcc.s_buff[0U]));
			y_quad_adjust = 1;
			y_nonzero = 1;
		} else {
			sLong2MultiWord(b_idx, (uint32_t *)(&r.s_buff[0U]), 2);
			MultiWordSignedWrap((uint32_t *)(&r.s_buff[0U]), 2, 30U,
					    (uint32_t *)(&yAcc.s_buff[0U]));
			y_quad_adjust = -1;
			y_nonzero = (b_idx > 0);
		}
		b_idx = thetaX[idx];
		if (b_idx < 0) {
			sLong2MultiWord(b_idx, (uint32_t *)(&r1.s_buff[0U]), 2);
			MultiWordSignedWrap((uint32_t *)(&r1.s_buff[0U]), 2, 30U,
					    (uint32_t *)(&r2.s_buff[0U]));
			MultiWordNeg((uint32_t *)(&r2.s_buff[0U]), (uint32_t *)(&r.s_buff[0U]), 2);
			MultiWordSignedWrap((uint32_t *)(&r.s_buff[0U]), 2, 30U,
					    (uint32_t *)(&xAcc.s_buff[0U]));
			x_quad_adjust = 1;
		} else {
			sLong2MultiWord(b_idx, (uint32_t *)(&r.s_buff[0U]), 2);
			MultiWordSignedWrap((uint32_t *)(&r.s_buff[0U]), 2, 30U,
					    (uint32_t *)(&xAcc.s_buff[0U]));
			x_quad_adjust = -1;
		}
		zN = 0;
		c = xAcc;
		b_c = yAcc;
		r = r3;
		for (b_idx = 0; b_idx < 31; b_idx++) {
			if (sMultiWordLt((uint32_t *)(&yAcc.s_buff[0U]),
					 (uint32_t *)(&r3.s_buff[0U]), 2)) {
				MultiWordSub((uint32_t *)(&xAcc.s_buff[0U]),
					     (uint32_t *)(&b_c.s_buff[0U]),
					     (uint32_t *)(&r2.s_buff[0U]), 2);
				MultiWordSignedWrap((uint32_t *)(&r2.s_buff[0U]), 2, 30U,
						    (uint32_t *)(&xAcc.s_buff[0U]));
				MultiWordAdd((uint32_t *)(&yAcc.s_buff[0U]),
					     (uint32_t *)(&c.s_buff[0U]),
					     (uint32_t *)(&r2.s_buff[0U]), 2);
				MultiWordSignedWrap((uint32_t *)(&r2.s_buff[0U]), 2, 30U,
						    (uint32_t *)(&yAcc.s_buff[0U]));
				zN -= inpLUT[b_idx];
			} else {
				MultiWordAdd((uint32_t *)(&xAcc.s_buff[0U]),
					     (uint32_t *)(&b_c.s_buff[0U]),
					     (uint32_t *)(&r2.s_buff[0U]), 2);
				MultiWordSignedWrap((uint32_t *)(&r2.s_buff[0U]), 2, 30U,
						    (uint32_t *)(&xAcc.s_buff[0U]));
				MultiWordSub((uint32_t *)(&yAcc.s_buff[0U]),
					     (uint32_t *)(&c.s_buff[0U]),
					     (uint32_t *)(&r2.s_buff[0U]), 2);
				MultiWordSignedWrap((uint32_t *)(&r2.s_buff[0U]), 2, 30U,
						    (uint32_t *)(&yAcc.s_buff[0U]));
				zN += inpLUT[b_idx];
			}
			sMultiWordShr((uint32_t *)(&xAcc.s_buff[0U]), 2,
				      (uint32_t)((int32_t)(b_idx + 1)), (uint32_t *)(&c.s_buff[0U]),
				      2);
			sMultiWordShr((uint32_t *)(&yAcc.s_buff[0U]), 2,
				      (uint32_t)((int32_t)(b_idx + 1)),
				      (uint32_t *)(&b_c.s_buff[0U]), 2);
		}
		if (zN <= INT32_MIN) {
			b_idx = INT32_MAX;
		} else {
			b_idx = -zN;
		}
		sLong2MultiWord(zN, (uint32_t *)(&r4.s_buff[0U]), 2);
		MultiWordSignedWrap((uint32_t *)(&r4.s_buff[0U]), 2, 31U,
				    (uint32_t *)(&r5.s_buff[0U]));
		r4 = r6;
		MultiWordSub((uint32_t *)(&r5.s_buff[0U]), (uint32_t *)(&r6.s_buff[0U]),
			     (uint32_t *)(&r7.s_buff[0U]), 2);
		sMultiWord2MultiWord((uint32_t *)(&r7.s_buff[0U]), 2, (uint32_t *)(&r8.s_buff[0U]),
				     3);
		sMultiWordShl((uint32_t *)(&r8.s_buff[0U]), 3, 31U, (uint32_t *)(&r9.s_buff[0U]),
			      3);
		sMultiWord2sMultiWordSat((uint32_t *)(&r9.s_buff[0U]), 3,
					 (uint32_t *)(&r1.s_buff[0U]), 2);
		sMultiWordShr((uint32_t *)(&r1.s_buff[0U]), 2, 31U, (uint32_t *)(&r2.s_buff[0U]),
			      2);
		sLong2MultiWord(zN, (uint32_t *)(&r10.s_buff[0U]), 2);
		MultiWordSignedWrap((uint32_t *)(&r10.s_buff[0U]), 2, 31U,
				    (uint32_t *)(&r11.s_buff[0U]));
		MultiWordSub((uint32_t *)(&r6.s_buff[0U]), (uint32_t *)(&r11.s_buff[0U]),
			     (uint32_t *)(&r5.s_buff[0U]), 2);
		sMultiWord2MultiWord((uint32_t *)(&r5.s_buff[0U]), 2, (uint32_t *)(&r8.s_buff[0U]),
				     3);
		sMultiWordShl((uint32_t *)(&r8.s_buff[0U]), 3, 31U, (uint32_t *)(&r9.s_buff[0U]),
			      3);
		sMultiWord2sMultiWordSat((uint32_t *)(&r9.s_buff[0U]), 3,
					 (uint32_t *)(&r7.s_buff[0U]), 2);
		sMultiWordShr((uint32_t *)(&r7.s_buff[0U]), 2, 31U, (uint32_t *)(&r1.s_buff[0U]),
			      2);
		if (y_nonzero) {
			if (x_quad_adjust) {
				if (y_quad_adjust) {
					cdcatan2Th[idx] =
					    sMultiWord2sLongSat((uint32_t *)(&r2.s_buff[0U]), 2);
				} else {
					cdcatan2Th[idx] =
					    sMultiWord2sLongSat((uint32_t *)(&r1.s_buff[0U]), 2);
				}
			} else if (y_quad_adjust) {
				cdcatan2Th[idx] = b_idx;
			} else {
				cdcatan2Th[idx] = zN;
			}
		} else if (x_quad_adjust) {
			cdcatan2Th[idx] = 1686629713;
		} else {
			cdcatan2Th[idx] = 0;
		}
	}
}

