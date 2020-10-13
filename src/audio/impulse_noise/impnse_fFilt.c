/*
 * SPDX-License-Identifier: BSD-3-Clause
 * File: impnse.c
 *
 * Author :Sriram Shastry & Ingalsuo, Seppo
 * Coyright(c) 2017 Intel Corporation. All rights reserved.
 */

/* Include Files */
#include <sof/audio/impulse_noise/impnse.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <sof/audio/impulse_noise/rtwtypes.h>
#include <sof/audio/impulse_noise/stdmlib.h>

/* Type Definitions */
#ifndef struct_emxArray_int32_T_101
#define struct_emxArray_int32_T_101

struct emxArray_int32_T_101
{
  int data[101];
  int size[1];
};

#endif                                 /*struct_emxArray_int32_T_101*/

#ifndef typedef_emxArray_int32_T_101
#define typedef_emxArray_int32_T_101

typedef struct emxArray_int32_T_101 emxArray_int32_T_101;

#endif                                 /*typedef_emxArray_int32_T_101*/

/* Function Declarations */
static unsigned int MultiWord2uLong(const unsigned int u[]);
static void MultiWordAdd(const unsigned int u1[], const unsigned int u2[],
  unsigned int y[], int n);
static void MultiWordSignedWrap(const unsigned int u1[], int n1, unsigned int n2,
  unsigned int y[]);
static void MultiWordSub(const unsigned int u1[], const unsigned int u2[],
  unsigned int y[], int n);
static int asr_s32(int u, unsigned int n);
static void detectpeek(const short x_1[101], unsigned int minpeakh, unsigned
  char locs_data[], int locs_size[1], unsigned short pks_data[], int pks_size[1]);
static int div_s32_convergent(int numerator, int denominator);
static void merge(int idx_data[], short x_data[], int offset, int np, int nq,
                  int iwork_data[], short xwork_data[]);
static double rt_remd(double u0, double u1);
static void sLong2MultiWord(int u, unsigned int y[], int n);
static int sMultiWordCmp(const unsigned int u1[], const unsigned int u2[], int n);
static bool sMultiWordLe(const unsigned int u1[], const unsigned int u2[], int n);
static bool sMultiWordLt(const unsigned int u1[], const unsigned int u2[], int n);
static void sMultiWordMul(const unsigned int u1[], int n1, const unsigned int
  u2[], int n2, unsigned int y[], int n);
static void sMultiWordShl(const unsigned int u1[], int n1, unsigned int n2,
  unsigned int y[], int n);
static void sMultiWordShr(const unsigned int u1[], int n1, unsigned int n2,
  unsigned int y[], int n);
static void sMultiWordShrConv(const unsigned int u1[], int n1, unsigned int n2,
  unsigned int y[], int n);
static void sort(short x_data[], const int x_size[1]);
static void sortIdx(short x_data[], const int x_size[1], int idx_data[], int
                    idx_size[1]);
static void uLong2MultiWord(unsigned int u, unsigned int y[], int n);

/* Function Definitions */

/*
 * Arguments    : const unsigned int u[]
 * Return Type  : unsigned int
 */
static unsigned int MultiWord2uLong(const unsigned int u[])
{
  return u[0];
}

/*
 * Arguments    : const unsigned int u1[]
 *                const unsigned int u2[]
 *                unsigned int y[]
 *                int n
 * Return Type  : void
 */
static void MultiWordAdd(const unsigned int u1[], const unsigned int u2[],
  unsigned int y[], int n)
{
  int i;
  unsigned int u1i;
  unsigned int yi;
  int carry = 0;
  for (i = 0; i < n; i++)
  {
    u1i = u1[i];
    yi = (u1i + u2[i]) + ((unsigned int)carry);
    y[i] = yi;
    if (((unsigned int)carry) != 0U)
    {
      carry = (int)((yi <= u1i) ? 1 : 0);
    }
    else
    {
      carry = (int)((yi < u1i) ? 1 : 0);
    }
  }
}

/*
 * Arguments    : const unsigned int u1[]
 *                int n1
 *                unsigned int n2
 *                unsigned int y[]
 * Return Type  : void
 */
static void MultiWordSignedWrap(const unsigned int u1[], int n1, unsigned int n2,
  unsigned int y[])
{
  int n1m1;
  unsigned int mask;
  unsigned int ys;
  n1m1 = n1 - 1;
  if (0 <= (n1m1 - 1))
  {
    memcpy(&y[0], &u1[0], ((unsigned int)n1m1) * (sizeof(unsigned int)));
  }

  mask = (1U << (31U - n2));
  if ((u1[n1m1] & mask) != 0U)
  {
    ys = MAX_uint32_T;
  }
  else
  {
    ys = 0U;
  }

  mask = (mask << 1U) - 1U;
  y[n1m1] = (u1[n1m1] & mask) | ((~mask) & ys);
}

/*
 * Arguments    : const unsigned int u1[]
 *                const unsigned int u2[]
 *                unsigned int y[]
 *                int n
 * Return Type  : void
 */
static void MultiWordSub(const unsigned int u1[], const unsigned int u2[],
  unsigned int y[], int n)
{
  int i;
  unsigned int u1i;
  unsigned int yi;
  int borrow = 0;
  for (i = 0; i < n; i++)
  {
    u1i = u1[i];
    yi = (u1i - u2[i]) - ((unsigned int)borrow);
    y[i] = yi;
    if (((unsigned int)borrow) != 0U)
    {
      borrow = (int)((yi >= u1i) ? 1 : 0);
    }
    else
    {
      borrow = (int)((yi > u1i) ? 1 : 0);
    }
  }
}

/*
 * Arguments    : int u
 *                unsigned int n
 * Return Type  : int
 */
static int asr_s32(int u, unsigned int n)
{
  int y;
  if (u >= 0)
  {
    y = (int)((unsigned int)(((unsigned int)u) >> n));
  }
  else
  {
    y = (-((int)((unsigned int)(((unsigned int)((int)(-1 - u))) >> n)))) - 1;
  }

  return y;
}

/*
 * function [locs, pks]=detectpeek(x_1,minpeakdist,minpeakh)
 * Arguments    : const short x_1[101]
 *                unsigned int minpeakh
 *                unsigned char locs_data[]
 *                int locs_size[1]
 *                unsigned short pks_data[]
 *                int pks_size[1]
 * Return Type  : void
 */
static void detectpeek(const short x_1[101], unsigned int minpeakh, unsigned
  char locs_data[], int locs_size[1], unsigned short pks_data[], int pks_size[1])
{
  int i;
  unsigned short x[101];
  unsigned short u;
  int k0;
  bool bv[99];
  bool exitg1;
  bool bv1[99];
  int nxin;
  signed char ii_data[99];
  double b_u;
  bool c_data[99];
  int k;

  /*  detectpeek : detectpeek function */
  /*  args */
  /*  X           - Indata  raw data */
  /*  minpeakdist - minimum data sample for analysis */
  /*  minpeakh    - minimum heigth to declare a peek */
  for (i = 0; i < 101; i++)
  {
    x[i] = (unsigned short)(((unsigned short)x_1[i]) << ((unsigned int)1));
  }

  /*  Find all maxima and ties */
  for (i = 0; i < 99; i++)
  {
    u = x[i + 1];
    bv[i] = (u >= x[i]);
    bv1[i] = (u >= x[i + 2]);
  }

  i = 0;
  k0 = 0;
  exitg1 = false;
  while ((!exitg1) && (k0 < 99))
  {
    if ((bv[k0]) && (bv1[k0]))
    {
      i++;
      ii_data[i - 1] = (signed char)(k0 + 1);
      if (i >= 99)
      {
        exitg1 = true;
      }
      else
      {
        k0++;
      }
    }
    else
    {
      k0++;
    }
  }

  if (1 > i)
  {
    nxin = 0;
  }
  else
  {
    nxin = i;
  }

  for (k0 = 0; k0 < nxin; k0++)
  {
    b_u = ((double)ii_data[k0]) + 1.0;
    if (rt_remd(b_u, 2.0) != 0.5)
    {
      b_u += 0.5;
    }

    
    b_u = impnse_floor(b_u);
    
    locs_data[k0] = (unsigned char)b_u;
  }

  /*  If no minpeakdist specified, default to 1. */
  /*  If there's a minpeakheight */
  for (k0 = 0; k0 < nxin; k0++)
  {
    c_data[k0] = ((((unsigned int)x[((int)locs_data[k0]) - 1]) << ((unsigned int)
      1)) <= minpeakh);
  }

  i = 0;
  for (k = 0; k < nxin; k++)
  {
    i += c_data[k] ? 1 : 0;
  }

  i = nxin - i;
  k0 = -1;
  for (k = 0; k < nxin; k++)
  {
    if (((k + 1) > nxin) || (!c_data[k]))
    {
      k0++;
      locs_data[k0] = locs_data[k];
    }
  }

  if (1 > i)
  {
    locs_size[0] = 0;
  }
  else
  {
    locs_size[0] = i;
  }

  pks_size[0] = locs_size[0];
  i = locs_size[0];
  for (k0 = 0; k0 < i; k0++)
  {
    pks_data[k0] = x[((int)locs_data[k0]) - 1];
  }
}

/*
 * Arguments    : int numerator
 *                int denominator
 * Return Type  : int
 */
static int div_s32_convergent(int numerator, int denominator)
{
  int quotient;
  unsigned int absNumerator;
  unsigned int absDenominator;
  unsigned int tempAbsQuotient;
  if (denominator == 0)
  {
    if (numerator >= 0)
    {
      quotient = MAX_int32_T;
    }
    else
    {
      quotient = MIN_int32_T;
    }
  }
  else
  {
    if (numerator < 0)
    {
      absNumerator = (~((unsigned int)numerator)) + 1U;
    }
    else
    {
      absNumerator = (unsigned int)numerator;
    }

    if (denominator < 0)
    {
      absDenominator = (~((unsigned int)denominator)) + 1U;
    }
    else
    {
      absDenominator = (unsigned int)denominator;
    }

    tempAbsQuotient = absNumerator / absDenominator;
    absNumerator %= absDenominator;
    absNumerator <<= 1U;
    if ((absNumerator > absDenominator) || ((absNumerator == absDenominator) &&
         ((tempAbsQuotient & 1U) != 0U)))
    {
      tempAbsQuotient++;
    }

    if ((numerator < 0) != (denominator < 0))
    {
      quotient = -((int)tempAbsQuotient);
    }
    else
    {
      quotient = (int)tempAbsQuotient;
    }
  }

  return quotient;
}

/*
 * Arguments    : int idx_data[]
 *                short x_data[]
 *                int offset
 *                int np
 *                int nq
 *                int iwork_data[]
 *                short xwork_data[]
 * Return Type  : void
 */
static void merge(int idx_data[], short x_data[], int offset, int np, int nq,
                  int iwork_data[], short xwork_data[])
{
  int n_tmp;
  int iout;
  int p;
  int i;
  int q;
  int exitg1;
  if (nq != 0)
  {
    n_tmp = np + nq;
    for (iout = 0; iout < n_tmp; iout++)
    {
      i = offset + iout;
      iwork_data[iout] = idx_data[i];
      xwork_data[iout] = x_data[i];
    }

    p = 0;
    q = np;
    iout = offset - 1;
    do
    {
      exitg1 = 0;
      iout++;
      if (xwork_data[p] <= xwork_data[q])
      {
        idx_data[iout] = iwork_data[p];
        x_data[iout] = xwork_data[p];
        if ((p + 1) < np)
        {
          p++;
        }
        else
        {
          exitg1 = 1;
        }
      }
      else
      {
        idx_data[iout] = iwork_data[q];
        x_data[iout] = xwork_data[q];
        if ((q + 1) < n_tmp)
        {
          q++;
        }
        else
        {
          q = iout - p;
          for (iout = p + 1; iout <= np; iout++)
          {
            i = q + iout;
            idx_data[i] = iwork_data[iout - 1];
            x_data[i] = xwork_data[iout - 1];
          }

          exitg1 = 1;
        }
      }
    }
    while (exitg1 == 0);
  }
}

/*
 * Arguments    : double u0
 *                double u1
 * Return Type  : double
 */
static double rt_remd(double u0, double u1)
{
  double y;
  double q;
  
  if ((u1 != 0.0) && (u1 != impnse_trunc(u1)))
  {
    q = fabs(u0 / u1);
    // if (fabs(q - floor(q + 0.5)) <= (DBL_EPSILON * q))
    if (fabs(q - impnse_floor(q + 0.5)) <= (DBL_EPSILON * q))
    {
      y = 0.0;
    }
    else
    {
       y = mod16(u0, u1);
    }
  }
  else
  {
       y = mod16(u0, u1);
  }

  return y;
}

/*
 * Arguments    : int u
 *                unsigned int y[]
 *                int n
 * Return Type  : void
 */
static void sLong2MultiWord(int u, unsigned int y[], int n)
{
  unsigned int yi;
  int i;
  y[0] = (unsigned int)u;
  if (u < 0)
  {
    yi = MAX_uint32_T;
  }
  else
  {
    yi = 0U;
  }

  for (i = 1; i < n; i++)
  {
    y[i] = yi;
  }
}

/*
 * Arguments    : const unsigned int u1[]
 *                const unsigned int u2[]
 *                int n
 * Return Type  : int
 */
static int sMultiWordCmp(const unsigned int u1[], const unsigned int u2[], int n)
{
  int y;
  unsigned int su1;
  int i;
  unsigned int u2i;
  su1 = u1[n - 1] & 2147483648U;
  if ((su1 ^ (u2[n - 1] & 2147483648U)) != 0U)
  {
    if (su1 != 0U)
    {
      y = -1;
    }
    else
    {
      y = 1;
    }
  }
  else
  {
    y = 0;
    i = n;
    while ((y == 0) && (i > 0))
    {
      i--;
      su1 = u1[i];
      u2i = u2[i];
      if (su1 != u2i)
      {
        if (su1 > u2i)
        {
          y = 1;
        }
        else
        {
          y = -1;
        }
      }
    }
  }

  return y;
}

/*
 * Arguments    : const unsigned int u1[]
 *                const unsigned int u2[]
 *                int n
 * Return Type  : bool
 */
static bool sMultiWordLe(const unsigned int u1[], const unsigned int u2[], int n)
{
  return sMultiWordCmp(u1, u2, n) <= 0;
}

/*
 * Arguments    : const unsigned int u1[]
 *                const unsigned int u2[]
 *                int n
 * Return Type  : bool
 */
static bool sMultiWordLt(const unsigned int u1[], const unsigned int u2[], int n)
{
  return sMultiWordCmp(u1, u2, n) < 0;
}

/*
 * Arguments    : const unsigned int u1[]
 *                int n1
 *                const unsigned int u2[]
 *                int n2
 *                unsigned int y[]
 *                int n
 * Return Type  : void
 */
static void sMultiWordMul(const unsigned int u1[], int n1, const unsigned int
  u2[], int n2, unsigned int y[], int n)
{
  bool isNegative1;
  bool isNegative2;
  int cb1;
  int i;
  unsigned int cb;
  unsigned int u1i;
  int k;
  int a1;
  unsigned int yk;
  int a0;
  int cb2;
  int ni;
  int j;
  int b1;
  int b0;
  unsigned int w01;
  unsigned int t;
  isNegative1 = ((u1[n1 - 1] & 2147483648U) != 0U);
  isNegative2 = ((u2[n2 - 1] & 2147483648U) != 0U);
  cb1 = 1;

  /* Initialize output to zero */
  if (0 <= (n - 1))
  {
    memset(&y[0], 0, ((unsigned int)n) * (sizeof(unsigned int)));
  }

  for (i = 0; i < n1; i++)
  {
    cb = 0U;
    u1i = u1[i];
    if (isNegative1)
    {
      u1i = (~u1i) + ((unsigned int)cb1);
      cb1 = (int)((u1i < ((unsigned int)cb1)) ? 1 : 0);
    }

    a1 = (int)((unsigned int)(u1i >> 16U));
    a0 = (int)((unsigned int)(u1i & 65535U));
    cb2 = 1;
    ni = n - i;
    if (n2 <= ni)
    {
      ni = n2;
    }

    k = i;
    for (j = 0; j < ni; j++)
    {
      u1i = u2[j];
      if (isNegative2)
      {
        u1i = (~u1i) + ((unsigned int)cb2);
        cb2 = (int)((u1i < ((unsigned int)cb2)) ? 1 : 0);
      }

      b1 = (int)((unsigned int)(u1i >> 16U));
      b0 = (int)((unsigned int)(u1i & 65535U));
      u1i = ((unsigned int)a1) * ((unsigned int)b0);
      w01 = ((unsigned int)a0) * ((unsigned int)b1);
      yk = y[k] + cb;
      cb = (unsigned int)((yk < cb) ? 1 : 0);
      t = ((unsigned int)a0) * ((unsigned int)b0);
      yk += t;
      cb += (unsigned int)((yk < t) ? 1 : 0);
      t = (u1i << 16U);
      yk += t;
      cb += (unsigned int)((yk < t) ? 1 : 0);
      t = (w01 << 16U);
      yk += t;
      cb += (unsigned int)((yk < t) ? 1 : 0);
      y[k] = yk;
      cb += (u1i >> 16U);
      cb += (w01 >> 16U);
      cb += ((unsigned int)a1) * ((unsigned int)b1);
      k++;
    }

    if (k < n)
    {
      y[k] = cb;
    }
  }

  /* Apply sign */
  if (isNegative1 != isNegative2)
  {
    cb = 1U;
    for (k = 0; k < n; k++)
    {
      yk = (~y[k]) + cb;
      y[k] = yk;
      cb = (unsigned int)((yk < cb) ? 1 : 0);
    }
  }
}

/*
 * Arguments    : const unsigned int u1[]
 *                int n1
 *                unsigned int n2
 *                unsigned int y[]
 *                int n
 * Return Type  : void
 */
static void sMultiWordShl(const unsigned int u1[], int n1, unsigned int n2,
  unsigned int y[], int n)
{
  int nb;
  unsigned int ys;
  int nc;
  unsigned int u1i;
  int i;
  unsigned int nl;
  unsigned int nr;
  unsigned int yi;
  nb = ((int)n2) / 32;
  if ((u1[n1 - 1] & 2147483648U) != 0U)
  {
    ys = MAX_uint32_T;
  }
  else
  {
    ys = 0U;
  }

  if (nb > n)
  {
    nc = n;
  }
  else
  {
    nc = nb;
  }

  u1i = 0U;
  if (0 <= (nc - 1))
  {
    memset(&y[0], 0, ((unsigned int)nc) * (sizeof(unsigned int)));
  }

  for (i = 0; i < nc; i++)
  {
  }

  if (nb < n)
  {
    nl = n2 - (((unsigned int)nb) * 32U);
    nb += n1;
    if (nb > n)
    {
      nb = n;
    }

    nb -= i;
    if (nl > 0U)
    {
      nr = 32U - nl;
      for (nc = 0; nc < nb; nc++)
      {
        yi = (u1i >> nr);
        u1i = u1[nc];
        y[i] = yi | (u1i << nl);
        i++;
      }

      if (i < n)
      {
        y[i] = (u1i >> nr) | (ys << nl);
        i++;
      }
    }
    else
    {
      for (nc = 0; nc < nb; nc++)
      {
        y[i] = u1[nc];
        i++;
      }
    }
  }

  while (i < n)
  {
    y[i] = ys;
    i++;
  }
}

/*
 * Arguments    : const unsigned int u1[]
 *                int n1
 *                unsigned int n2
 *                unsigned int y[]
 *                int n
 * Return Type  : void
 */
static void sMultiWordShr(const unsigned int u1[], int n1, unsigned int n2,
  unsigned int y[], int n)
{
  int nb;
  int i;
  unsigned int ys;
  int nc;
  unsigned int nr;
  int i1;
  unsigned int nl;
  unsigned int u1i;
  unsigned int yi;
  nb = ((int)n2) / 32;
  i = 0;
  if ((u1[n1 - 1] & 2147483648U) != 0U)
  {
    ys = MAX_uint32_T;
  }
  else
  {
    ys = 0U;
  }

  if (nb < n1)
  {
    nc = n + nb;
    if (nc > n1)
    {
      nc = n1;
    }

    nr = n2 - (((unsigned int)nb) * 32U);
    if (nr > 0U)
    {
      nl = 32U - nr;
      u1i = u1[nb];
      for (i1 = nb + 1; i1 < nc; i1++)
      {
        yi = (u1i >> nr);
        u1i = u1[i1];
        y[i] = yi | (u1i << nl);
        i++;
      }

      if (nc < n1)
      {
        yi = u1[nc];
      }
      else
      {
        yi = ys;
      }

      y[i] = (u1i >> nr) | (yi << nl);
      i++;
    }
    else
    {
      for (i1 = nb; i1 < nc; i1++)
      {
        y[i] = u1[i1];
        i++;
      }
    }
  }

  while (i < n)
  {
    y[i] = ys;
    i++;
  }
}

/*
 * Arguments    : const unsigned int u1[]
 *                int n1
 *                unsigned int n2
 *                unsigned int y[]
 *                int n
 * Return Type  : void
 */
static void sMultiWordShrConv(const unsigned int u1[], int n1, unsigned int n2,
  unsigned int y[], int n)
{
  unsigned int n2m1;
  int nb;
  bool doRoundUp;
  unsigned int maskHalfLSB;
  unsigned int maskLSB;
  unsigned int mask;
  int i;
  unsigned int u1b;
  bool tieCond = false;
  n2m1 = n2 - 1U;
  nb = ((int)n2m1) / 32;
  if (nb < n1)
  {
    n2m1 -= (unsigned int)((int)(nb * 32));
    maskHalfLSB = (1U << n2m1);
    maskLSB = (2U << n2m1);
    mask = maskLSB - 1U;
    u1b = u1[nb];

    /* BitsBeingShiftedOff>LSB/2 || (Tie condition && Odd LSB condition) */
    doRoundUp = ((u1b & mask) > maskHalfLSB);
    tieCond = ((u1b & mask) == maskHalfLSB);
    if (tieCond && (nb > 0))
    {
      i = 0;
      while ((!doRoundUp) && (i < nb))
      {
        doRoundUp = (u1[i] != 0U);
        i++;
      }
    }

    if ((!doRoundUp) && tieCond)
    {
      if (nb > 0)
      {
        i = 0;
        while (tieCond && (i < nb))
        {
          tieCond = (u1[i] == 0U);
          i++;
        }
      }

      if (n2m1 < 31U)
      {
        doRoundUp = (tieCond && ((u1b & maskLSB) != 0U));
      }
      else
      {
        if ((nb < (n1 - 1)) && (n2m1 == 31U))
        {
          doRoundUp = (tieCond && ((u1[nb + 1] & 1U) != 0U));
        }
      }
    }
  }
  else
  {
    /* Roundup negative inputs when only sign bit is left after rshift */
    doRoundUp = ((u1[n1 - 1] & 2147483648U) != 0U);
  }

  sMultiWordShr(u1, n1, n2, y, n);
  i = 0;
  while (doRoundUp && (i < n))
  {
    y[i]++;
    doRoundUp = ((y[i] == 0U) && doRoundUp);
    i++;
  }
}

/*
 * Arguments    : short x_data[]
 *                const int x_size[1]
 * Return Type  : void
 */
static void sort(short x_data[], const int x_size[1])
{
  int dim;
  int j;
  int vlen;
  int vwork_size[1];
  int vstride;
  int k;
  short vwork_data[101];
  emxArray_int32_T_101 b_vwork_data;
  dim = 0;
  if (x_size[0] != 1)
  {
    dim = -1;
  }

  if ((dim + 2) <= 1)
  {
    j = x_size[0];
  }
  else
  {
    j = 1;
  }

  vlen = j - 1;
  vwork_size[0] = j;
  vstride = 1;
  for (k = 0; k <= dim; k++)
  {
    vstride *= x_size[0];
  }

  for (j = 0; j < vstride; j++)
  {
    for (k = 0; k <= vlen; k++)
    {
      vwork_data[k] = x_data[j + (k * vstride)];
    }

    sortIdx(vwork_data, vwork_size, b_vwork_data.data, b_vwork_data.size);
    for (k = 0; k <= vlen; k++)
    {
      x_data[j + (k * vstride)] = vwork_data[k];
    }
  }
}

/*
 * Arguments    : short x_data[]
 *                const int x_size[1]
 *                int idx_data[]
 *                int idx_size[1]
 * Return Type  : void
 */
static void sortIdx(short x_data[], const int x_size[1], int idx_data[], int
                    idx_size[1])
{
  signed char unnamed_idx_0;
  int i1;
  int n;
  short x4[4];
  unsigned char idx4[4];
  int iwork_data[101];
  short xwork_data[101];
  int nQuartets;
  int j;
  int i2;
  int i;
  int i3;
  int i4;
  signed char perm[4];
  short x4_tmp;
  short b_x4_tmp;
  short c_x4_tmp;
  int idx_tmp;
  unnamed_idx_0 = (signed char)x_size[0];
  idx_size[0] = (int)unnamed_idx_0;
  i1 = (int)unnamed_idx_0;
  if (0 <= (i1 - 1))
  {
    memset(&idx_data[0], 0, ((unsigned int)i1) * (sizeof(int)));
  }

  if (x_size[0] != 0)
  {
    n = x_size[0];
    x4[0] = 0;
    idx4[0] = 0U;
    x4[1] = 0;
    idx4[1] = 0U;
    x4[2] = 0;
    idx4[2] = 0U;
    x4[3] = 0;
    idx4[3] = 0U;
    i1 = (int)unnamed_idx_0;
    if (0 <= (i1 - 1))
    {
      memset(&iwork_data[0], 0, ((unsigned int)i1) * (sizeof(int)));
    }

    i1 = x_size[0];
    if (0 <= (i1 - 1))
    {
      memset(&xwork_data[0], 0, ((unsigned int)i1) * (sizeof(short)));
    }

    nQuartets = asr_s32(x_size[0], 2U);
    for (j = 0; j < nQuartets; j++)
    {
      i = j * 4;
      idx4[0] = (unsigned char)((int)(i + 1));
      idx4[1] = (unsigned char)((int)(i + 2));
      idx4[2] = (unsigned char)((int)(i + 3));
      idx4[3] = (unsigned char)((int)(i + 4));
      x4[0] = x_data[i];
      x4_tmp = x_data[i + 1];
      x4[1] = x4_tmp;
      b_x4_tmp = x_data[i + 2];
      x4[2] = b_x4_tmp;
      c_x4_tmp = x_data[i + 3];
      x4[3] = c_x4_tmp;
      if (x_data[i] <= x4_tmp)
      {
        i1 = 1;
        i2 = 2;
      }
      else
      {
        i1 = 2;
        i2 = 1;
      }

      if (b_x4_tmp <= c_x4_tmp)
      {
        i3 = 3;
        i4 = 4;
      }
      else
      {
        i3 = 4;
        i4 = 3;
      }

      x4_tmp = x4[i1 - 1];
      b_x4_tmp = x4[i3 - 1];
      if (x4_tmp <= b_x4_tmp)
      {
        x4_tmp = x4[i2 - 1];
        if (x4_tmp <= b_x4_tmp)
        {
          perm[0] = (signed char)i1;
          perm[1] = (signed char)i2;
          perm[2] = (signed char)i3;
          perm[3] = (signed char)i4;
        }
        else if (x4_tmp <= x4[i4 - 1])
        {
          perm[0] = (signed char)i1;
          perm[1] = (signed char)i3;
          perm[2] = (signed char)i2;
          perm[3] = (signed char)i4;
        }
        else
        {
          perm[0] = (signed char)i1;
          perm[1] = (signed char)i3;
          perm[2] = (signed char)i4;
          perm[3] = (signed char)i2;
        }
      }
      else
      {
        b_x4_tmp = x4[i4 - 1];
        if (x4_tmp <= b_x4_tmp)
        {
          if (x4[i2 - 1] <= b_x4_tmp)
          {
            perm[0] = (signed char)i3;
            perm[1] = (signed char)i1;
            perm[2] = (signed char)i2;
            perm[3] = (signed char)i4;
          }
          else
          {
            perm[0] = (signed char)i3;
            perm[1] = (signed char)i1;
            perm[2] = (signed char)i4;
            perm[3] = (signed char)i2;
          }
        }
        else
        {
          perm[0] = (signed char)i3;
          perm[1] = (signed char)i4;
          perm[2] = (signed char)i1;
          perm[3] = (signed char)i2;
        }
      }

      i4 = ((int)perm[0]) - 1;
      idx_data[i] = (int)idx4[i4];
      idx_tmp = ((int)perm[1]) - 1;
      idx_data[i + 1] = (int)idx4[idx_tmp];
      i1 = ((int)perm[2]) - 1;
      idx_data[i + 2] = (int)idx4[i1];
      i2 = ((int)perm[3]) - 1;
      idx_data[i + 3] = (int)idx4[i2];
      x_data[i] = x4[i4];
      x_data[i + 1] = x4[idx_tmp];
      x_data[i + 2] = x4[i1];
      x_data[i + 3] = x4[i2];
    }

    i2 = nQuartets * 4;
    i3 = (x_size[0] - i2) - 1;
    if ((i3 + 1) > 0)
    {
      for (nQuartets = 0; nQuartets <= i3; nQuartets++)
      {
        i1 = i2 + nQuartets;
        idx4[nQuartets] = (unsigned char)((int)(i1 + 1));
        x4[nQuartets] = x_data[i1];
      }

      perm[1] = 0;
      perm[2] = 0;
      perm[3] = 0;
      switch (i3 + 1)
      {
       case 1:
        perm[0] = 1;
        break;

       case 2:
        if (x4[0] <= x4[1])
        {
          perm[0] = 1;
          perm[1] = 2;
        }
        else
        {
          perm[0] = 2;
          perm[1] = 1;
        }
        break;

       default:
        if (x4[0] <= x4[1])
        {
          if (x4[1] <= x4[2])
          {
            perm[0] = 1;
            perm[1] = 2;
            perm[2] = 3;
          }
          else if (x4[0] <= x4[2])
          {
            perm[0] = 1;
            perm[1] = 3;
            perm[2] = 2;
          }
          else
          {
            perm[0] = 3;
            perm[1] = 1;
            perm[2] = 2;
          }
        }
        else if (x4[0] <= x4[2])
        {
          perm[0] = 2;
          perm[1] = 1;
          perm[2] = 3;
        }
        else if (x4[1] <= x4[2])
        {
          perm[0] = 2;
          perm[1] = 3;
          perm[2] = 1;
        }
        else
        {
          perm[0] = 3;
          perm[1] = 2;
          perm[2] = 1;
        }
        break;
      }

      for (nQuartets = 0; nQuartets <= i3; nQuartets++)
      {
        i4 = ((int)perm[nQuartets]) - 1;
        idx_tmp = i2 + nQuartets;
        idx_data[idx_tmp] = (int)idx4[i4];
        x_data[idx_tmp] = x4[i4];
      }
    }

    if (n > 1)
    {
      i4 = asr_s32(n, 2U);
      i3 = 4;
      while (i4 > 1)
      {
        if ((i4 & 1) != 0)
        {
          i4--;
          i1 = i3 * i4;
          i2 = n - i1;
          if (i2 > i3)
          {
            merge(idx_data, x_data, i1, i3, i2 - i3, iwork_data, xwork_data);
          }
        }

        i1 = i3 * 2;
        i4 = asr_s32(i4, 1U);
        for (nQuartets = 0; nQuartets < i4; nQuartets++)
        {
          merge(idx_data, x_data, nQuartets * i1, i3, i3, iwork_data, xwork_data);
        }

        i3 = i1;
      }

      if (n > i3)
      {
        merge(idx_data, x_data, 0, i3, n - i3, iwork_data, xwork_data);
      }
    }
  }
}

/*
 * Arguments    : unsigned int u
 *                unsigned int y[]
 *                int n
 * Return Type  : void
 */
static void uLong2MultiWord(unsigned int u, unsigned int y[], int n)
{
  y[0] = u;
  if (1 <= (n - 1))
  {
    memset(&y[1], 0, ((unsigned int)((int)(n + -1))) * (sizeof(unsigned int)));
  }
}

/*
 * function [AudioSteam]  = impnse_fixpt(ImpnseIsOut)
 * Arguments    : const struct0_T ImpnseIsOut
 * Return Type  : struct0_T
 */
struct0_T impnse_fixpt(const struct0_T ImpnseIsOut)
{
  struct0_T AudioSteam;
  int i;
  int idx;
  unsigned char b_idx;
  int i1;
  int i2;
  int sumx;
  int k;
  short b0;
  int i3;
  int b_i;
  int64m_T y[101];
  int64m_T a1;
  int64m_T r;
  int64m_T r1;
  static const int64m_T r2 =
  {
    {
      0U, 0U
    }                                  /* chunks */
  };

  int64m_T r3;
  static const int64m_T r4 =
  {
    {
      MAX_uint32_T, 1023U
    }                                  /* chunks */
  };

  static const int64m_T r5 =
  {
    {
      0U, 4294966272U
    }                                  /* chunks */
  };

  unsigned int u;
  unsigned int u1;
  int b_y;
  int64m_T r6;
  int64m_T r7;
  int64m_T r8;
  int64m_T r9;
  int64m_T r10;
  int64m_T r11;
  int64m_T r12;
  int64m_T r13;
  int64m_T r14;
  short c_y[101];
  unsigned char locs_data[99];
  int locs_size[1];
  unsigned short fmo_2_data[99];
  int fmo_2_size[1];
  int c_size_idx_0;
  bool x_data[9999];
  unsigned char c_data[9999];
  int nz[101];
  bool b;
  bool bv[101];
  signed char tmp_data[101];
  short b_x_data[101];
  signed char b_tmp_data[101];

  AudioSteam = ImpnseIsOut;

  /* 'impnse_fixpt:24' for idx = fi(1, 0, 1, 0, fm) :s_chan */
  i = (int)ImpnseIsOut.NumChan;
  for (idx = 0; idx < i; idx++)
  {
    b_idx = (unsigned char)((int)((idx + 1) & 3));

    i1 = ((int)b_idx) - 1;
    i2 = 101 * i1;
    if ((((int)AudioSteam.f[i2]) & 4194304) != 0)
    {
      sumx = ((int)AudioSteam.f[i2]) | -4194304;
    }
    else
    {
      sumx = ((int)AudioSteam.f[i2]) & 4194303;
    }

    for (k = 0; k < 100; k++)
    {
      i3 = sumx + ((int)AudioSteam.f[(k + i2) + 1]);
      if ((i3 & 4194304) != 0)
      {
        sumx = i3 | -4194304;
      }
      else
      {
        sumx = i3 & 4194303;
      }
    }

    b0 = (short)div_s32_convergent(sumx, 101);
    for (b_i = 0; b_i < 101; b_i++)
    {
      i3 = (int)AudioSteam.f[b_i + i2];
      if ((i3 & 65536) != 0)
      {
        i3 |= -65536;
      }
      else
      {
        i3 &= 65535;
      }

      if ((((int)b0) & 65536) != 0)
      {
        sumx = ((int)b0) | -65536;
      }
      else
      {
        sumx = ((int)b0) & 65535;
      }

      i3 -= sumx;
      if ((i3 & 65536) != 0)
      {
        i3 |= -65536;
      }
      else
      {
        i3 &= 65535;
      }

      u = (unsigned int)i3;
      u1 = (unsigned int)i3;
      sMultiWordMul((unsigned int *)(&u), 1, (unsigned int *)(&u1), 1, (unsigned
        int *)(&y[b_i].chunks[0U]), 2);
    }

    MultiWordSignedWrap((unsigned int *)(&y[0].chunks[0U]), 2, 23U, (unsigned
      int *)(&a1.chunks[0U]));
    for (k = 0; k < 100; k++)
    {
      MultiWordSignedWrap((unsigned int *)(&y[k + 1].chunks[0U]), 2, 23U,
                          (unsigned int *)(&r.chunks[0U]));
      MultiWordAdd((unsigned int *)(&a1.chunks[0U]), (unsigned int *)(&r.chunks
        [0U]), (unsigned int *)(&r1.chunks[0U]), 2);
      MultiWordSignedWrap((unsigned int *)(&r1.chunks[0U]), 2, 23U, (unsigned
        int *)(&a1.chunks[0U]));
    }

    r1 = r2;
    if (sMultiWordLt((unsigned int *)(&a1.chunks[0U]), (unsigned int *)
                     (&r2.chunks[0U]), 2))
    {
      r3 = r5;
    }
    else
    {
      r3 = r4;
    }

    b_y = 0;
    r = r2;
    if (!sMultiWordLe((unsigned int *)(&r3.chunks[0U]), (unsigned int *)
                      (&r2.chunks[0U]), 2))
    {
      for (b_i = 20; b_i >= 0; b_i--)
      {
        i3 = b_y | (1 << ((unsigned int)b_i));
        u = (unsigned int)i3;
        u1 = (unsigned int)i3;
        sMultiWordMul((unsigned int *)(&u), 1, (unsigned int *)(&u1), 1,
                      (unsigned int *)(&r6.chunks[0U]), 2);
        if (sMultiWordLe((unsigned int *)(&r6.chunks[0U]), (unsigned int *)
                         (&r3.chunks[0U]), 2))
        {
          b_y = i3;
        }
      }

      if (b_y < 2097151)
      {
        i3 = b_y + 1;
        u = (unsigned int)i3;
        u1 = (unsigned int)i3;
        sMultiWordMul((unsigned int *)(&u), 1, (unsigned int *)(&u1), 1,
                      (unsigned int *)(&r7.chunks[0U]), 2);
        MultiWordSignedWrap((unsigned int *)(&r3.chunks[0U]), 2, 20U, (unsigned
          int *)(&r8.chunks[0U]));
        MultiWordSub((unsigned int *)(&r7.chunks[0U]), (unsigned int *)
                     (&r8.chunks[0U]), (unsigned int *)(&r10.chunks[0U]), 2);
        MultiWordSignedWrap((unsigned int *)(&r10.chunks[0U]), 2, 20U, (unsigned
          int *)(&r12.chunks[0U]));
        MultiWordSignedWrap((unsigned int *)(&r3.chunks[0U]), 2, 20U, (unsigned
          int *)(&r8.chunks[0U]));
        u = (unsigned int)b_y;
        u1 = (unsigned int)b_y;
        sMultiWordMul((unsigned int *)(&u), 1, (unsigned int *)(&u1), 1,
                      (unsigned int *)(&r14.chunks[0U]), 2);
        MultiWordSub((unsigned int *)(&r8.chunks[0U]), (unsigned int *)
                     (&r14.chunks[0U]), (unsigned int *)(&r7.chunks[0U]), 2);
        MultiWordSignedWrap((unsigned int *)(&r7.chunks[0U]), 2, 20U, (unsigned
          int *)(&r10.chunks[0U]));
        if (sMultiWordLt((unsigned int *)(&r12.chunks[0U]), (unsigned int *)
                         (&r10.chunks[0U]), 2))
        {
          b_y = i3;
        }
      }
    }

    /*  Locate peaks */
    if ((((int)AudioSteam.f[i2]) & 4194304) != 0)
    {
      sumx = ((int)AudioSteam.f[i2]) | -4194304;
    }
    else
    {
      sumx = ((int)AudioSteam.f[i2]) & 4194303;
    }

    for (k = 0; k < 100; k++)
    {
      i3 = sumx + ((int)AudioSteam.f[(k + i2) + 1]);
      if ((i3 & 4194304) != 0)
      {
        sumx = i3 | -4194304;
      }
      else
      {
        sumx = i3 & 4194303;
      }
    }

    sLong2MultiWord((int)((short)div_s32_convergent(sumx, 101)), (unsigned int *)
                    (&r9.chunks[0U]), 2);
    sMultiWordShl((unsigned int *)(&r9.chunks[0U]), 2, 18U, (unsigned int *)
                  (&r11.chunks[0U]), 2);
    MultiWordSignedWrap((unsigned int *)(&r11.chunks[0U]), 2, 28U, (unsigned int
      *)(&r13.chunks[0U]));
    uLong2MultiWord(((unsigned int)((unsigned short)(((unsigned short)b_y) <<
      ((unsigned int)4)))) * 40960U, (unsigned int *)(&r9.chunks[0U]), 2);
    MultiWordSignedWrap((unsigned int *)(&r9.chunks[0U]), 2, 28U, (unsigned int *)
                        (&r11.chunks[0U]));
    MultiWordAdd((unsigned int *)(&r13.chunks[0U]), (unsigned int *)
                 (&r11.chunks[0U]), (unsigned int *)(&r14.chunks[0U]), 2);
    MultiWordSignedWrap((unsigned int *)(&r14.chunks[0U]), 2, 28U, (unsigned int
      *)(&r8.chunks[0U]));
    sMultiWordShrConv((unsigned int *)(&r8.chunks[0U]), 2, 16U, (unsigned int *)
                      (&r7.chunks[0U]), 2);
    AudioSteam.minHeight[i1] = MultiWord2uLong((unsigned int *)(&r7.chunks[0U]))
      & 262143U;

    /*  min height to be 'peak';  3 std above mean */
    for (k = 0; k < 101; k++)
    {
      i3 = k + i2;
      if (AudioSteam.f[i3] < 0)
      {
        c_y[k] = (short)(-AudioSteam.f[i3]);
      }
      else
      {
        c_y[k] = AudioSteam.f[i3];
      }
    }

    detectpeek(c_y, AudioSteam.minHeight[((int)b_idx) - 1], locs_data, locs_size,
               fmo_2_data, fmo_2_size);

    c_size_idx_0 = (int)((signed char)locs_size[0]);
    if (((signed char)locs_size[0]) != 0)
    {
      sumx = (int)((((signed char)locs_size[0]) != 1) ? 1 : 0);
      i3 = c_size_idx_0 - 1;
      for (k = 0; k < 101; k++)
      {
        for (b_y = 0; b_y <= i3; b_y++)
        {
          c_data[b_y + (c_size_idx_0 * k)] = (unsigned char)(((unsigned int)
            ((int)(k + 1))) - ((unsigned int)locs_data[sumx * b_y]));
        }
      }
    }

    sumx = c_size_idx_0 * 101;
    for (i3 = 0; i3 < sumx; i3++)
    {
      x_data[i3] = (((int)c_data[i3]) == 0);
    }

    if (c_size_idx_0 == 0)
    {
      memset(&nz[0], 0, 101U * (sizeof(int)));
    }
    else
    {
      for (b_i = 0; b_i < 101; b_i++)
      {
        sumx = b_i * c_size_idx_0;
        nz[b_i] = x_data[sumx] ? 1 : 0;
        for (k = 2; k <= c_size_idx_0; k++)
        {
          nz[b_i] += x_data[(sumx + k) - 1] ? 1 : 0;
        }
      }
    }

    for (i3 = 0; i3 < 101; i3++)
    {
      AudioSteam.isOut[i3 + i2] = (nz[i3] == 1);
    }

    /*  resize size mismatch */
    b_y = 0;
    for (b_i = 0; b_i < 101; b_i++)
    {
      b = !AudioSteam.isOut[b_i + i2];
      bv[b_i] = b;
      if (b)
      {
        b_y++;
      }
    }

    sumx = 0;
    for (b_i = 0; b_i < 101; b_i++)
    {
      if (bv[b_i])
      {
        tmp_data[sumx] = (signed char)(b_i + 1);
        sumx++;
      }
    }

    sumx = asr_s32(b_y + ((int)((b_y < 0) ? 1 : 0)), 1U);
    locs_size[0] = b_y;
    for (i3 = 0; i3 < b_y; i3++)
    {
      b_x_data[i3] = AudioSteam.f[(((int)tmp_data[i3]) + i2) - 1];
    }

    sort(b_x_data, locs_size);
    if (b_y == (sumx * 2))
    {
      i3 = ((int)b_x_data[sumx - 1]) + ((int)b_x_data[sumx]);
      if ((i3 & 65536) != 0)
      {
        i3 |= -65536;
      }
      else
      {
        i3 &= 65535;
      }

      b0 = (short)(asr_s32(i3, 1U) + ((int)((((i3 & 1) == 1) && ((i3 & 2) != 0))
        ? 1 : 0)));
    }
    else
    {
      b0 = b_x_data[sumx];
    }

    AudioSteam.repValue[i1] = (((unsigned int)b0) << ((unsigned int)7)) &
      8388607U;

    /*  you could also use NaN, 0, etc.... */
    /*  make a copy */
    b_y = 0;
    for (b_i = 0; b_i < 101; b_i++)
    {
      i3 = b_i + i2;
      if ((((int)AudioSteam.f[i3]) & 65536) != 0)
      {
        AudioSteam.fFilt[i3] = ((int)AudioSteam.f[i3]) | -65536;
      }
      else
      {
        AudioSteam.fFilt[i3] = ((int)AudioSteam.f[i3]) & 65535;
      }

      if (AudioSteam.isOut[i3])
      {
        b_y++;
      }
    }

    sumx = 0;
    for (b_i = 0; b_i < 101; b_i++)
    {
      if (AudioSteam.isOut[b_i + i2])
      {
        b_tmp_data[sumx] = (signed char)(b_i + 1);
        sumx++;
      }
    }

    i1 = (int)((unsigned int)(((AudioSteam.repValue[i1] >> ((unsigned int)7)) +
      ((unsigned int)((((AudioSteam.repValue[i1] & 64U) != 0U) &&
                       ((AudioSteam.repValue[i1] & 63U) != 0U)) ? 1 : 0))) +
                ((unsigned int)((((AudioSteam.repValue[i1] & 127U) == 64U) &&
      ((AudioSteam.repValue[i1] & 128U) != 0U)) ? 1 : 0))));
    if ((i1 & 65536) != 0)
    {
      sumx = i1 | -65536;
    }
    else
    {
      sumx = i1 & 65535;
    }

    for (i1 = 0; i1 < b_y; i1++)
    {
      AudioSteam.fFilt[(((int)b_tmp_data[i1]) + i2) - 1] = sumx;
    }
  }

  return AudioSteam;
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void impnse_initialize(void)
{
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void impnse_terminate(void)
{
  /* (no terminate code required) */
}

/*
 * function [ImpnseIsOut] = init_struc_fixpt
 * Arguments    : void
 * Return Type  : struct0_T
 */
struct0_T init_struc_fixpt(void)
{
  struct0_T ImpnseIsOut;
  int i;
  static const short iv[202] =
  {
    396, 465, 519, 398, 409, 404, 20349, 563, 497, 503, 477, 570, 785, 759, 576,
    446, 406, 397, 408, 558, 398, 417, 412, 413, 405, 465, 433, 440, 502, 441,
    406, 429, 583, 523, 396, 465, 519, 398, 409, 404, -20349, 563, 497, 503, 477,
    570, 785, 759, 576, 446, 406, 397, 408, 558, 398, 417, 412, 413, 405, 465,
    433, 440, 502, 441, 406, 429, 583, 523, 396, 465, 519, 398, 409, 404, 20349,
    563, 497, 503, 477, 570, 785, 759, 576, 446, 406, 397, 408, 558, 398, 417,
    412, 413, 405, 465, 433, 440, 502, 441, 406, 429, 583, 396, 465, 519, 398,
    409, 404, -20349, 563, 497, 503, 477, 570, 785, 759, 576, 446, 406, 397, 408,
    558, 398, 417, 412, 413, 405, 465, 433, 440, 502, 441, 406, 429, 583, 523,
    396, 465, 519, 398, 409, 404, 20349, 563, 497, 503, 477, 570, 785, 759, 576,
    446, 406, 397, 408, 558, 398, 417, 412, 413, 405, 465, 433, 440, 502, 441,
    406, 429, 583, 523, 396, 465, 519, 398, 409, 404, -20349, 563, 497, 503, 477,
    570, 785, 759, 576, 446, 406, 397, 408, 558, 398, 417, 412, 413, 405, 465,
    433, 440, 502, 441, 406, 429, 583
  };

  ImpnseIsOut.NumChan = 2U;

  ImpnseIsOut.minHeight[0] = 0U;
  ImpnseIsOut.repValue[0] = 0U;
  ImpnseIsOut.minHeight[1] = 0U;
  ImpnseIsOut.repValue[1] = 0U;

  for (i = 0; i < 202; i++)
  {
    ImpnseIsOut.f[i] = iv[i];
    ImpnseIsOut.fFilt[i] = 0;
    ImpnseIsOut.isOut[i] = false;
  }

  return ImpnseIsOut;
}



/*
 * File trailer for impnse.c
 *
 * [EOF]
 */
