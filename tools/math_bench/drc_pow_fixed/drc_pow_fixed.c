#pragma warning (disable : 4013)
#pragma warning (disable : 4996)
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "init_data.h"
#include "decibels.h"
#include "/../Users/shastry/source/Work/Audio/SourceCode/a_v_03/Models/drc_Existing_CodeWrapper/common/typdef.h"
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include "decibels.h"
#include "exp_small_fixed.h"

#define  TEST_VECTOR 64
#define q_mult(a, b, qa, qb, qy) ((int32_t)Q_MULTSR_32X32((int64_t)(a), b, qa, qb, qy))
//int32_t exp_small_fixed(int32_t x);
/*
 * Input x is Q6.26: (1, 32.0) x = [>0:0.1:90] all positive
 *       y is Q2.30: (-2.0, 2.0)
 * Output is Q12.20: max 2048.0
 */
inline int32_t drc_pow_fixed(int32_t x, int32_t y, int16_t i, FILE *fd)
{
	/* x^y = expf(y * log(x)) */ /*[2^-2,2^2] */
	if (x <= 0)
		return 0;

	fprintf(fd, " %13li \n", exp_fixed(q_mult(y, drc_log_fixed(x), 30, 26, 27)));
	return exp_fixed(q_mult(y, drc_log_fixed(x), 30, 26, 27));
}

/*
 * Input is Q6.26: max 32.0
 * Output range ~ (-inf, 3.4657); regulated to Q6.26: (-32.0, 32.0)
 */
inline int32_t drc_log_fixed(int32_t x)
{
	const int32_t LOG10 = Q_CONVERT_FLOAT(2.3025850929940457f, 29);
	int32_t log10_x;

	if (x <= 0)
		return Q_CONVERT_FLOAT(-30.0f, 26);

	/* log(x) = log(10) * log10(x) */
	log10_x = log10_fixed(x); /* Q6.26 */
	return q_mult(LOG10, log10_x, 29, 26, 26);
}

/*
 * Input is Q6.26: max 32.0
 * Output range ~ (-inf, 1.505); regulated to Q6.26: (-32.0, 32.0)
 */
static inline int32_t log10_fixed(int32_t x)
{
#define QC 26
	/* Coefficients obtained from:
	 * fpminimax(log10(x), 5, [|SG...|], [1/2;sqrt(2)/2], absolute);
	 * max err ~= 6.088e-8
	 */
	const int32_t ONE_OVER_SQRT2 = Q_CONVERT_FLOAT(0.70710678118654752f, 30); /* 1/sqrt(2) */
	const int32_t A5 = Q_CONVERT_FLOAT(1.131880283355712890625f, QC);
	const int32_t A4 = Q_CONVERT_FLOAT(-4.258677959442138671875f, QC);
	const int32_t A3 = Q_CONVERT_FLOAT(6.81631565093994140625f, QC);
	const int32_t A2 = Q_CONVERT_FLOAT(-6.1185703277587890625f, QC);
	const int32_t A1 = Q_CONVERT_FLOAT(3.6505267620086669921875f, QC);
	const int32_t A0 = Q_CONVERT_FLOAT(-1.217894077301025390625f, QC);
	const int32_t LOG10_2 = Q_CONVERT_FLOAT(0.301029995663981195214f, QC);
	int32_t e;
	int32_t exp; /* Q31.1 */
	int32_t x2, x4; /* Q2.30 */
	int32_t A5Xx, A3Xx;

	x = rexp_fixed(x, 26, &e); /* Q2.30 */
	exp = (int32_t)e << 1; /* Q_CONVERT_FLOAT(e, 1) */

	if (x > ONE_OVER_SQRT2) {
		x = q_mult(x, ONE_OVER_SQRT2, 30, 30, 30);
		exp += 1; /* Q_CONVERT_FLOAT(0.5, 1) */
	}

	x2 = q_mult(x, x, 30, 30, 30);
	x4 = q_mult(x2, x2, 30, 30, 30);
	A5Xx = q_mult(A5, x, QC, 30, QC);
	A3Xx = q_mult(A3, x, QC, 30, QC);
	return q_mult((A5Xx + A4), x4, QC, 30, QC) + q_mult((A3Xx + A2), x2, QC, 30, QC)
		+ q_mult(A1, x, QC, 30, QC) + A0 + q_mult(exp, LOG10_2, 1, QC, QC);
#undef QC
}

/*
 * Input depends on precision_x
 * Output range [0.5, 1); regulated to Q2.30
 */
static inline int32_t rexp_fixed(int32_t x, int precision_x, int* e)
{
	int bit = 31 - norm_int32(x);

	*e = bit - precision_x;

	if (bit > 30)
		return Q_SHIFT_RND(x, bit, 30);
	if (bit < 30)
		return Q_SHIFT_LEFT(x, bit, 30);
	return x;
}

/* Count the left shift amount to normalize a 32 bit signed integer value
 * without causing overflow. Input value 0 will result to 31.
 */
int norm_int32(int32_t val)
{
	int s;
	int32_t n;

	if (!val)
		return 31;

	n = val << 1;
	s = 0;
	if (val > 0) {
		while (n > 0) {
			n = n << 1;
			s++;
		}
	}
	else { /* we don't goto else section */
		while (n < 0) {
			n = n << 1;
			s++;
		}
	}
	return s;
}

#undef q_mult
/*
 * Input x is Q6.26: (-32.0, 32.0)
 *       y is Q2.30: (-2.0, 2.0)
 * Output is Q12.20: max 2048.0
 *
 */
int main(void)
{
	/*int32_t x[TEST_VECTOR],y = 0;*/
	int32_t y[41];
	uint32_t x[320];
	int16_t i;
	int8_t j;
	mkdir("Results", 0777);
	FILE *fd = fopen("Results/drc_pow_fixed.txt", "w");
	if (!fd) {
		fprintf(fd, "error: unable to open file");
	}
	else {
		//*Input x is Q6.26: (1, 32.0) x = [> 0:0.1 : 90] all positive
		//	* y is Q2.30 : (-2.0, 2.0)
		//	* Output is Q12.20 : max 2048.0
		/*fprintf(fd, " %10s  %10s %11s %10s %14s \n-----------+----------+-----------+-----------+--------------+\n", "idx-i", "idx-j", "in-xQ6/26", "in-yQ2/30", "out-powQ12/20");*/
		fprintf(fd, " %10s  %10s %11s %10s %14s \n", "idx-i", "idx-j", "in-xQ6/26", "in-yQ2/30", "out-powQ12/20");
		init_data_fixpt(x, y);
		for (i = 0; i < 319; i++)
		{
			for (j = 0; j < 41; j++) {
				fprintf(fd, " %10d  %10d %11li %11li ", i, j, x[i], y[j]);
				drc_pow_fixed(x[i], y[j], i, fd);		/*return exp_fixed(q_mult(y, drc_log_fixed(x), 30, 26, 27));*/
			}
		}
	}
}
