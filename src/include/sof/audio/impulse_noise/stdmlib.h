/*
 * SPDX-License-Identifier: BSD-3-Clause
 * File: impnse_fixpt.c
 *
 * Author :Sriram Shastry & Ingalsuo, Seppo
 * Coyright(c) 2017 Intel Corporation. All rights reserved.
 */

#ifndef STDLIBFCN_H
#define STDLIBFCN_H

/* Include Files */
#include <stddef.h>
#include <sof/audio/impulse_noise/stdmlib.h>
#include <sof/audio/impulse_noise/rtwtypes.h>



#include "stdint.h"


/* Function Definitions */
typedef union
{
  double value;
  struct
  {
    uint32_t lsw;
    uint32_t msw;
  } parts;
} ieee_double_shape_type;


/* Get two 32 bit ints from a double.  */

#define EXTRACT_WORDS(ix0,ix1,d)				\
do {								\
  ieee_double_shape_type ew_u;					\
  ew_u.value = (d);						\
  (ix0) = ew_u.parts.msw;					\
  (ix1) = ew_u.parts.lsw;					\
} while (/*CONSTCOND*/0)



 #define INSERT_WORDS(d,ix0,ix1)					\
 do {								\
   ieee_double_shape_type iw_u;					\
   iw_u.parts.msw = (ix0);					\
   iw_u.parts.lsw = (ix1);					\
   (d) = iw_u.value;						\
 } while (/*CONSTCOND*/0)

/* Get the more significant 32 bit int from a double.  */

#define GET_HIGH_WORD(i,d)					\
do {								\
  ieee_double_shape_type gh_u;					\
   gh_u.value = (d);						\
   (i) = gh_u.parts.msw;						\
 } while (/*CONSTCOND*/0)


/* Set the more significant 32 bits of a double from an int.  */

#define SET_HIGH_WORD(d,v)					\
do {								\
  ieee_double_shape_type sh_u;					\
  sh_u.value = (d);						\
  sh_u.parts.msw = (v);						\
  (d) = sh_u.value;						\
} while (/*CONSTCOND*/0)

//double ldexp(double, int);
#define	DBL_EXPBITS	11
#define	DBL_EXP_INFNAN	2047
#define	DBL_EXP_BIAS	1023
#define	DBL_FRACHBITS	20
#define	DBL_FRACLBITS	32
#define	DBL_FRACBITS	(DBL_FRACHBITS + DBL_FRACLBITS)


struct ieee_double {
#if _BYTE_ORDER == _BIG_ENDIAN
    uint32_T	dbl_sign : 1;
    uint32_T	dbl_exp : DBL_EXPBITS;
    uint32_T	dbl_frach : DBL_FRACHBITS;
    uint32_T	dbl_fracl : DBL_FRACLBITS;
#else
    uint32_T	dbl_fracl : DBL_FRACLBITS;
    uint32_T	dbl_frach : DBL_FRACHBITS;
    uint32_T	dbl_exp : DBL_EXPBITS;
    uint32_T	dbl_sign : 1;
#endif
};

union ieee_double_u {
    double			dblu_d;
    struct ieee_double	dblu_dbl;
};

/* Function Declarations */
//extern double  fabs(double x);
extern double  impnse_trunc(double x);
extern double  impnse_floor(double x);
extern double  impnse_mod(double s_num, double s_den);
#endif

/*
 * File trailer for stdmlib.h
 *
 * [EOF]
 */
