// SPDX - License - Identifier: BSD - 3 - Clause
//
//Copyright(c) 2021 Intel Corporation.All rights reserved.
//
//Author : Shriram Shastry <malladi.sastry@linux.intel.com>
#include <stdint.h>
#include <stdint.h>
//#include "numbers.h"
#define SINE_C_Q20 341782638 /* 2*SINE_NQUART/pi in Q12.20 */
#define SINE_NQUART 512 /* Must be 2^N */
#define SINE_TABLE_SIZE (SINE_NQUART + 1)
#define TEST_VECTOR_SINE 21
#define ONE_Q20         Q_CONVERT_FLOAT(1.0, 20)	  /* Use Q12.20 */
#define ONE_Q23         Q_CONVERT_FLOAT(1.0, 23)	  /* Use Q9.23 */
#define TWO_Q27         Q_CONVERT_FLOAT(2.0, 27)	  /* Use Q5.27 */
#define MINUS_TWO_Q27   Q_CONVERT_FLOAT(-2.0, 27)	  /* Use Q5.27 */
#define Q_CONVERT_FLOAT(f, qy) \
	((int32_t)(((const double)f) * ((int64_t)1 << (const int)qy) + 0.5))
#define Q_SHIFT_RND(x, src_q, dst_q) \
	((((x) >> ((src_q) - (dst_q) - 1)) + 1) >> 1)

/* Fractional multiplication with shift and round
 * Note that the parameters px and py must be cast to (int64_t) if other type.
 */
#define Q_SHIFT_LEFT(x, src_q, dst_q) ((x) << ((dst_q) - (src_q)))
#define Q_MULTSR_32X32(px, py, qx, qy, qp) \
	((((px) * (py) >> ((qx) + (qy) - (qp) - 1)) + 1) >> 1)
//#define q_mult(a, b, qa, qb, qy) ((int32_t)Q_MULTSR_32X32((int64_t)(a), b, qa, qb, qy)

#define LOG10_DIV20_Q27 Q_CONVERT_FLOAT(0.1151292546, 27) /* Use Q5.27 */   
#define LOG10_DIV20_Q27 Q_CONVERT_FLOAT(0.1151292546, 27) /* Use Q5.27 */
#define Q_SHIFT_RND(x, src_q, dst_q) \
	((((x) >> ((src_q) - (dst_q) - 1)) + 1) >> 1)
#define q_mult(a, b, qa, qb, qy) ((int32_t)Q_MULTSR_32X32((int64_t)(a), b, qa, qb, qy))
#define Q_SHIFT_RND(x, src_q, dst_q) \
	((((x) >> ((src_q) - (dst_q) - 1)) + 1) >> 1)
/* Alternative version since compiler does not allow (x >> -1) */
#define Q_SHIFT_LEFT(x, src_q, dst_q) ((x) << ((dst_q) - (src_q)))
/* Compute the number of shifts
 * This will result in a compiler overflow error if shift bits are out of
 * range as INT64_MAX/MIN is greater than 32 bit Q shift parameter
 */
#define Q_SHIFT_BITS_64(qx, qy, qz) \
	((qx + qy - qz) <= 63 ? (((qx + qy - qz) >= 0) ? \
	 (qx + qy - qz) : INT64_MIN) : INT64_MAX)
#define TEST_VECTOR_SINE 21
#define TEST_VECTOR_db2lin 1551
//static inline int64_t q_mults_32x32(int32_t x, int32_t y, const int shift_bits)
//{
//	return ((int64_t)x * y) >> shift_bits;
//}
//#define TEST_VECTOR 1551
#define MAX_int32_t                    ((int32_t)(2147483647))
#define MIN_int32_t                    ((int32_t)(-2147483647-1))
#define ABS(a) ((a<0)?(-a):(a))
/* Type Definitions */
#ifndef typedef_struct0_T
#define typedef_struct0_T

typedef struct {
	uint32_t u1;
	int32_t u2[1551];
}struct0_T;
#endif  
extern struct0_T init_struc_fixpt(void);


