// SPDX-License-Identifier: BSD-3-Clause
//
// Copyright(c) 2022 Intel Corporation. All rights reserved.
//
// Author: Shriram Shastry <malladi.sastry@linux.intel.com>
//
//
#include <sof/math/exp_fcn_hifi4.h>
#include <sof/common.h>
#include <stdbool.h>
#include <stdint.h>
#include <stddef.h>

#ifdef EXPONENTIAL_HIFI4
#include <xtensa/tie/xt_hifi4.h>

/*
 * Arguments	: int64_t in_0
 *		  int64_t in_1
 *		  uint64_t *ptroutbitshi
 *		  uint64_t *ptroutbitslo
 * Return Type	: void
 * Description:Perform element-wise multiplication on in_0 and in_1
 * while keeping the required product word length and fractional
 * length in mind. mul_s64 function divide the 64-bit quantities
 * into two 32-bit words,multiply the low words to produce the
 * lowest and second-lowest words in the result, then both pairs
 * of low and high words from different numbers to produce the
 * second and third lowest words in the result, and finally both
 * high words to produce the two highest words in the outcome.
 * Add them all up, taking carry into consideration.
 *
 * The 64 x 64 bit multiplication of operands in_0 and in_1 is
 * shown in the image below. The 64-bit operand in_0,in_1 is
 * represented by the notation in0_H, in1_H for the top 32 bits
 * and in0_L, in1_L for the bottom 32 bits.
 *
 *				in0_H : in0_L
 *			x	in1_H : in1_L
 *			---------------------
 *	P0			in0_L x in1_L
 *	P1		in0_H x in1_L		64 bit inner multiplication
 *	P2		in0_L x in1_H		64 bit inner multiplication
 *	P3	in0_H x in1_H
 *			--------------------
 *			[64 x 64 bit multiplication] sum of inner products
 * All combinations are multiplied by one another and then added.
 * Each inner product is moved into its proper power location.given the names
 * of the inner products, redoing the addition where 000 represents 32 zero
 * bits.The inner products can be added together in 64 bit addition.The sum
 * of two 64-bit numbers yields a 65-bit output.
 *		   (P0H:P0L)
 *		P1H(P1L:000)
 *		P2H(P2L:000)
 *	P3H:P3L(000:000)
 *	.......(aaa:P0L)
 * By combining P0H:P0L and P1L:000. This can lead to a carry, denote as CRY0.
 * The partial result is then multiplied by P2L:000.
 * We call it CRY1 because it has the potential to carry again.
 *	(CRY0 + CRY1)P0H:P0L
 *	(	 P1H)P1L:000
 *	(	 P2H)P2L:000
 *	(P3H:	 P3L)000:000
 *	--------------------
 *	(ccc:bbb)aaa:P0L
 * P1H, P2H, and P3H:P3L are added to the carry CRY0 + CRY1.This increase will
 * not result in an overflow.
 *
 */
static void mul_s64(ae_int64 in_0, ae_int64 in_1, ae_int64 *ptroutbitshi, ae_int64 *ptroutbitslo)
{
	ae_int64 producthihi, producthilo, productlolo;
	ae_int64 producthi, carry, product_hl_lh_h, product_hl_lh_l;

	ae_int32x2 in0_32 = AE_MOVINT32X2_FROMINT64(in_0);
	ae_int32x2 in1_32 = AE_MOVINT32X2_FROMINT64(in_1);

	ae_ep ep_lolo = AE_MOVEA(0);
	ae_ep ep_hilo = AE_MOVEA(0);
	ae_ep ep_HL_LH = AE_MOVEA(0);

	producthihi = AE_MUL32_HH(in0_32, in1_32);
	AE_MULZAAD32USEP_HL_LH(ep_hilo, producthilo, in1_32, in0_32);
	productlolo = AE_MUL32U_LL(in0_32, in1_32);

	product_hl_lh_h = AE_SRAI72(ep_hilo, producthilo, 32);
	product_hl_lh_l = AE_SLAI64(producthilo, 32);

	AE_ADD72(ep_lolo, productlolo, ep_HL_LH, product_hl_lh_l);

	carry = AE_SRAI72(ep_lolo, productlolo, 32);

	carry = AE_SRLI64(carry, 32);
	producthi = AE_ADD64(producthihi, carry);

	producthi = AE_ADD64(producthi, product_hl_lh_h);

	*ptroutbitslo = productlolo;
	*ptroutbitshi = producthi;
}

/*
 * Arguments	: int64_t a
 *		  int64_t b
 * Return Type	: int64_t
 */
static int64_t mul_s64_sr_sat_near(int64_t a, int64_t b)
{
	ae_int64 result;
	ae_int64 u64_chi;
	ae_int64 u64_clo;
	ae_int64 temp;

	mul_s64(a, b, &u64_chi, &u64_clo);

	ae_int64 roundup = AE_AND64(u64_clo, BIT_MASK_LOW_Q27P5);

	roundup = AE_SRLI64(roundup, 27);
	temp = AE_OR64(AE_SLAI64(u64_chi, 36), AE_SRLI64(u64_clo, 28));

	return AE_ADD64(temp, roundup);
}

int64_t onebyfact_Q63[19] = {
		4611686018427387904LL,
		1537228672809129301LL,
		384307168202282325LL,
		76861433640456465LL,
		12810238940076077LL,
		1830034134296582LL,
		228754266787072LL,
		25417140754119LL,
		2541714075411LL,
		231064915946LL,
		19255409662LL,
		1481185358LL,
		105798954LL,
		7053264LL,
		440829LL,
		25931LL,
		1441LL,
		76LL,
		4LL
};

/* f(x) = a^x, x is variable and a is base
 *
 * Arguments	: int32_t x(Q4.28)
 * input range -5 to 5
 *
 * Return Type	: int32_t ts(Q9.23)
 * output range 0.0067465305 to 148.4131488800
 *+------------------+-----------------+--------+--------+
 *| x		     | ts (returntype) |   x	|  ts	 |
 *+----+-----+-------+----+----+-------+--------+--------+
 *|WLen| FLen|Signbit|WLen|FLen|Signbit| Qformat| Qformat|
 *+----+-----+-------+----+----+-------+--------+--------+
 *| 32 | 28  |	1    | 32 | 23 |   0   | 4.28	| 9.23	 |
 *+------------------+-----------------+--------+--------+
 */
int32_t exp_int32(int32_t x)
{
	ae_int64 mp;
	ae_int64 outhi;
	ae_int64 outlo;
	ae_int64 qt;
	ae_int64 onebyfact;
	ae_int64 temp;

	ae_int64 *ponebyfact_Q63 = (ae_int64 *)&onebyfact_Q63[0];
	ae_int64 ts = TERMS_Q23P9;
	xtbool flag;

	int64_t qt_temp;
	int64_t b_n;

	mp = (x + LSHIFT_BITS) >> 14; //x in Q50.14

	mul_s64(mp, BIT_MASK_Q62P2, &outhi, &outlo);
	qt = AE_OR64(AE_SLAI64(outhi, 46), AE_SRLI64(outlo, 18));

	temp = AE_SRAI64(AE_ADD64(qt, QUOTIENT_SCALE), 35);

	ts = AE_ADD64(ts, temp);

	mp = mul_s64_sr_sat_near(mp, (int64_t)x);

	for (b_n = 0; b_n < ARRAY_SIZE(onebyfact_Q63); b_n++) {
		ae_int64 sign;

		AE_L64_IP(onebyfact, ponebyfact_Q63, 8);

		mul_s64(mp, onebyfact, &outhi, &outlo);
		qt = AE_OR64(AE_SLAI64(outhi, 45), AE_SRLI64(outlo, 19));

		temp = AE_SRAI64(AE_ADD64(qt, QUOTIENT_SCALE), 35);
		ts = AE_ADD64(ts, temp);

		mp = mul_s64_sr_sat_near(mp, (int64_t)x);

		sign = AE_NEG64(qt);

		flag = AE_LT64(qt, 0);
		AE_MOVT64(qt, sign, flag);

		qt_temp = qt;
		if (qt_temp < CONVERG_ERROR)
			break;
	}

	return AE_MOVAD32_L(AE_MOVINT32X2_FROMINT64(ts));
}
#endif
