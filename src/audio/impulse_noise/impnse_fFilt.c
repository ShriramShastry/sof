/*
 * SPDX-License-Identifier: BSD-3-Clause
 * File: impnse_fFilt.c
 *
 * Author :Sriram Shastry & Ingalsuo, Seppo
 * Coyright(c) 2017 Intel Corporation. All rights reserved.
 */
/* Include Files */
#include  <sof/audio/impulse_noise/impnse.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <sof/audio/impulse_noise/stdmlib.h>
/* Type Definitions */
#include "stdint.h"

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
static int MultiWord2sLong(const unsigned int u[]);
static unsigned int MultiWord2uLong(const unsigned int u[]);
static void MultiWordAdd(const unsigned int u1[], const unsigned int u2[],
  unsigned int y[], int n);
static void MultiWordIor(const unsigned int u1[], const unsigned int u2[],
  unsigned int y[], int n);
static void MultiWordNeg(const unsigned int u1[], unsigned int y[], int n);
static void MultiWordSetSignedMax(unsigned int y[], int n);
static void MultiWordSetSignedMin(unsigned int y[], int n);
static void MultiWordSignedWrap(const unsigned int u1[], int n1, unsigned int n2,
  unsigned int y[]);
static void MultiWordSub(const unsigned int u1[], const unsigned int u2[],
  unsigned int y[], int n);
static void MultiWordUnsignedWrap(const unsigned int u1[], int n1, unsigned int
  n2, unsigned int y[]);
static int asr_s32(int u, unsigned int n);
static void detectpeek(const int x_1[101], const uint64m_T minpeakh, unsigned
  char locs_data[], int locs_size[1], unsigned int pks_data[], int pks_size[1]);
static void merge(int idx_data[], int x_data[], int offset, int np, int nq, int
                  iwork_data[], int xwork_data[]);
static double rt_remd(double u0, double u1);
static void sLong2MultiWord(int u, unsigned int y[], int n);
static void sMultiWord2MultiWord(const unsigned int u1[], int n1, unsigned int
  y[], int n);
static int sMultiWordCmp(const unsigned int u1[], const unsigned int u2[], int n);
static void sMultiWordDivConv(const unsigned int u1[], int n1, const unsigned
  int u2[], int n2, unsigned int b_y1[], int m1, unsigned int y2[], int m2,
  unsigned int t1[], int l1, unsigned int t2[], int l2);
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
static void sort(int x_data[], const int x_size[1]);
static void sortIdx(int x_data[], const int x_size[1], int idx_data[], int
                    idx_size[1]);
static void uLong2MultiWord(unsigned int u, unsigned int y[], int n);
static void uMultiWord2MultiWord(const unsigned int u1[], int n1, unsigned int
  y[], int n);
static int uMultiWordCmp(const unsigned int u1[], const unsigned int u2[], int n);
static int uMultiWordCmpShr(const unsigned int u1[], const unsigned int u2[],
  int n);
static int uMultiWordDiv(unsigned int a[], int na, unsigned int b[], int nb,
  unsigned int q[], int nq, unsigned int r[], int nr);
static void uMultiWordInc(unsigned int y[], int n);
static bool uMultiWordLe(const unsigned int u1[], const unsigned int u2[], int n);
static void uMultiWordMul(const unsigned int u1[], int n1, const unsigned int
  u2[], int n2, unsigned int y[], int n);
static void uMultiWordShr(const unsigned int u1[], int n1, unsigned int n2,
  unsigned int y[], int n);
static void uMultiWordShrConv(const unsigned int u1[], int n1, unsigned int n2,
  unsigned int y[], int n);

/* Function Definitions */

/*
 * Arguments    : const unsigned int u[]
 * Return Type  : int
 */
static int MultiWord2sLong(const unsigned int u[])
{
  return (int)u[0];
}

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
 *                const unsigned int u2[]
 *                unsigned int y[]
 *                int n
 * Return Type  : void
 */
static void MultiWordIor(const unsigned int u1[], const unsigned int u2[],
  unsigned int y[], int n)
{
  int i;
  for (i = 0; i < n; i++)
  {
    y[i] = u1[i] | u2[i];
  }
}

/*
 * Arguments    : const unsigned int u1[]
 *                unsigned int y[]
 *                int n
 * Return Type  : void
 */
static void MultiWordNeg(const unsigned int u1[], unsigned int y[], int n)
{
  int i;
  unsigned int yi;
  int carry = 1;
  for (i = 0; i < n; i++)
  {
    yi = (~u1[i]) + ((unsigned int)carry);
    y[i] = yi;
    carry = (int)((yi < ((unsigned int)carry)) ? 1 : 0);
  }
}

/*
 * Arguments    : unsigned int y[]
 *                int n
 * Return Type  : void
 */
static void MultiWordSetSignedMax(unsigned int y[], int n)
{
  int n1;
  int i;
  n1 = n - 1;
  for (i = 0; i < n1; i++)
  {
    y[i] = MAX_uint32_T;
  }

  y[n1] = 2147483647U;
}

/*
 * Arguments    : unsigned int y[]
 *                int n
 * Return Type  : void
 */
static void MultiWordSetSignedMin(unsigned int y[], int n)
{
  int n1;
  n1 = n - 1;
  if (0 <= (n1 - 1))
  {
    memset(&y[0], 0, ((unsigned int)n1) * (sizeof(unsigned int)));
  }

  y[n1] = 2147483648U;
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
 * Arguments    : const unsigned int u1[]
 *                int n1
 *                unsigned int n2
 *                unsigned int y[]
 * Return Type  : void
 */
static void MultiWordUnsignedWrap(const unsigned int u1[], int n1, unsigned int
  n2, unsigned int y[])
{
  int n1m1;
  n1m1 = n1 - 1;
  if (0 <= (n1m1 - 1))
  {
    memcpy(&y[0], &u1[0], ((unsigned int)n1m1) * (sizeof(unsigned int)));
  }

  y[n1m1] = u1[n1m1] & ((1U << (32U - n2)) - 1U);
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
 * Arguments    : const int x_1[101]
 *                const uint64m_T minpeakh
 *                unsigned char locs_data[]
 *                int locs_size[1]
 *                unsigned int pks_data[]
 *                int pks_size[1]
 * Return Type  : void
 */
static void detectpeek(const int x_1[101], const uint64m_T minpeakh, unsigned
  char locs_data[], int locs_size[1], unsigned int pks_data[], int pks_size[1])
{
  int i;
  unsigned int x[101];
  unsigned int u;
  int k0;
  bool bv[99];
  bool exitg1;
  bool bv1[99];
  int nxin;
  signed char ii_data[99];
  double b_u;
  uint64m_T r;
  int k;
  uint64m_T r1;
  bool c_data[99];

  /*  function : detectpeek  */
  /*  args */
  /*  X           - Indata  raw data */
  /*  minpeakdist - minimum data sample for analysis */
  /*  minpeakh    - minimum heigth to declare a peek */
  /*  return */
  /*  peak value */
  /*  location for peak value */
  /* 'impnse_fixpt:63' fm = get_fimath(); */
  /* 'impnse_fixpt:64' x = fi(x_1, 0, 32, 32, fm); */
  /* 'impnse_fixpt:66' if size(x,fi(2, 0, 2, 0, fm))==1 */
  /* 'impnse_fixpt:66' x(:)=fi(x', 0, 32, 32, fm); */
  for (i = 0; i < 101; i++)
  {
    x[i] = (((unsigned int)x_1[i]) << ((unsigned int)1));
  }

  /*  Find all maxima and ties */
  /* 'impnse_fixpt:68' locs=fi(find(x(2:end-1)>=x(1:end-2) & x(2:end-1)>=x(3:end))+1, 0, 7, 0, fm); */
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
    uLong2MultiWord(x[((int)locs_data[k0]) - 1], (unsigned int *)(&r.chunks[0U]),
                    2);
    sMultiWordShl((unsigned int *)(&r.chunks[0U]), 2, 1U, (unsigned int *)
                  (&r1.chunks[0U]), 2);
    r = minpeakh;
    c_data[k0] = uMultiWordLe((unsigned int *)(&r1.chunks[0U]), (unsigned int *)
      (&minpeakh.chunks[0U]), 2);
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
 * Arguments    : int idx_data[]
 *                int x_data[]
 *                int offset
 *                int np
 *                int nq
 *                int iwork_data[]
 *                int xwork_data[]
 * Return Type  : void
 */
static void merge(int idx_data[], int x_data[], int offset, int np, int nq, int
                  iwork_data[], int xwork_data[])
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
    if (fabs(q - impnse_floor(q + 0.5)) <= (DBL_EPSILON * q))
    {
      y = 0.0;
    }
    else
    {
      y = impnse_mod(u0, u1);
    }
  }
  else
  {
    y = impnse_mod(u0, u1);
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
 *                int n1
 *                unsigned int y[]
 *                int n
 * Return Type  : void
 */
static void sMultiWord2MultiWord(const unsigned int u1[], int n1, unsigned int
  y[], int n)
{
  int nm;
  unsigned int u1i;
  int i;
  if (n1 < n)
  {
    nm = n1;
  }
  else
  {
    nm = n;
  }

  if (0 <= (nm - 1))
  {
    memcpy(&y[0], &u1[0], ((unsigned int)nm) * (sizeof(unsigned int)));
  }

  if (n > n1)
  {
    if ((u1[n1 - 1] & 2147483648U) != 0U)
    {
      u1i = MAX_uint32_T;
    }
    else
    {
      u1i = 0U;
    }

    for (i = nm; i < n; i++)
    {
      y[i] = u1i;
    }
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
 *                int n1
 *                const unsigned int u2[]
 *                int n2
 *                unsigned int b_y1[]
 *                int m1
 *                unsigned int y2[]
 *                int m2
 *                unsigned int t1[]
 *                int l1
 *                unsigned int t2[]
 *                int l2
 * Return Type  : void
 */
static void sMultiWordDivConv(const unsigned int u1[], int n1, const unsigned
  int u2[], int n2, unsigned int b_y1[], int m1, unsigned int y2[], int m2,
  unsigned int t1[], int l1, unsigned int t2[], int l2)
{
  bool numNeg;
  bool denNeg;
  int cmp;
  numNeg = ((u1[n1 - 1] & 2147483648U) != 0U);
  denNeg = ((u2[n2 - 1] & 2147483648U) != 0U);
  if (numNeg)
  {
    MultiWordNeg(u1, t1, n1);
  }
  else
  {
    sMultiWord2MultiWord(u1, n1, t1, l1);
  }

  if (denNeg)
  {
    MultiWordNeg(u2, t2, n2);
  }
  else
  {
    sMultiWord2MultiWord(u2, n2, t2, l2);
  }

  if (uMultiWordDiv(t1, l1, t2, l2, b_y1, m1, y2, m2) < 0)
  {
    if (numNeg)
    {
      MultiWordSetSignedMin(b_y1, m1);
    }
    else
    {
      MultiWordSetSignedMax(b_y1, m1);
    }
  }
  else
  {
    cmp = uMultiWordCmpShr(y2, t2, m2);
    if ((cmp > 0) || ((cmp == 0) && ((b_y1[0] & 1U) != 0U)))
    {
      uMultiWordInc(b_y1, m1);
    }

    if (numNeg != denNeg)
    {
      MultiWordNeg(b_y1, b_y1, m1);
    }
  }
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
 * Arguments    : int x_data[]
 *                const int x_size[1]
 * Return Type  : void
 */
static void sort(int x_data[], const int x_size[1])
{
  int dim;
  int j;
  int vlen;
  int vwork_size[1];
  int vstride;
  int k;
  int vwork_data[101];
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
 * Arguments    : int x_data[]
 *                const int x_size[1]
 *                int idx_data[]
 *                int idx_size[1]
 * Return Type  : void
 */
static void sortIdx(int x_data[], const int x_size[1], int idx_data[], int
                    idx_size[1])
{
  signed char unnamed_idx_0;
  int i3;
  int n;
  int x4[4];
  unsigned char idx4[4];
  int iwork_data[101];
  int xwork_data[101];
  int nQuartets;
  int j;
  int i4;
  int i;
  int nLeft;
  int i1;
  int nPairs;
  signed char perm[4];
  int i2;
  unnamed_idx_0 = (signed char)x_size[0];
  idx_size[0] = (int)unnamed_idx_0;
  i3 = (int)unnamed_idx_0;
  if (0 <= (i3 - 1))
  {
    memset(&idx_data[0], 0, ((unsigned int)i3) * (sizeof(int)));
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
    i3 = (int)unnamed_idx_0;
    if (0 <= (i3 - 1))
    {
      memset(&iwork_data[0], 0, ((unsigned int)i3) * (sizeof(int)));
    }

    i3 = x_size[0];
    if (0 <= (i3 - 1))
    {
      memset(&xwork_data[0], 0, ((unsigned int)i3) * (sizeof(int)));
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
      i3 = x_data[i + 1];
      x4[1] = i3;
      i4 = x_data[i + 2];
      x4[2] = i4;
      nLeft = x_data[i + 3];
      x4[3] = nLeft;
      if (x_data[i] <= i3)
      {
        i1 = 1;
        i2 = 2;
      }
      else
      {
        i1 = 2;
        i2 = 1;
      }

      if (i4 <= nLeft)
      {
        i3 = 3;
        i4 = 4;
      }
      else
      {
        i3 = 4;
        i4 = 3;
      }

      nLeft = x4[i1 - 1];
      nPairs = x4[i3 - 1];
      if (nLeft <= nPairs)
      {
        nLeft = x4[i2 - 1];
        if (nLeft <= nPairs)
        {
          perm[0] = (signed char)i1;
          perm[1] = (signed char)i2;
          perm[2] = (signed char)i3;
          perm[3] = (signed char)i4;
        }
        else if (nLeft <= x4[i4 - 1])
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
        nPairs = x4[i4 - 1];
        if (nLeft <= nPairs)
        {
          if (x4[i2 - 1] <= nPairs)
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

      nPairs = ((int)perm[0]) - 1;
      idx_data[i] = (int)idx4[nPairs];
      i2 = ((int)perm[1]) - 1;
      idx_data[i + 1] = (int)idx4[i2];
      i3 = ((int)perm[2]) - 1;
      idx_data[i + 2] = (int)idx4[i3];
      i4 = ((int)perm[3]) - 1;
      idx_data[i + 3] = (int)idx4[i4];
      x_data[i] = x4[nPairs];
      x_data[i + 1] = x4[i2];
      x_data[i + 2] = x4[i3];
      x_data[i + 3] = x4[i4];
    }

    i4 = nQuartets * 4;
    nLeft = (x_size[0] - i4) - 1;
    if ((nLeft + 1) > 0)
    {
      for (i1 = 0; i1 <= nLeft; i1++)
      {
        i3 = i4 + i1;
        idx4[i1] = (unsigned char)((int)(i3 + 1));
        x4[i1] = x_data[i3];
      }

      perm[1] = 0;
      perm[2] = 0;
      perm[3] = 0;
      switch (nLeft + 1)
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

      for (i1 = 0; i1 <= nLeft; i1++)
      {
        nPairs = ((int)perm[i1]) - 1;
        i2 = i4 + i1;
        idx_data[i2] = (int)idx4[nPairs];
        x_data[i2] = x4[nPairs];
      }
    }

    if (n > 1)
    {
      nPairs = asr_s32(n, 2U);
      nLeft = 4;
      while (nPairs > 1)
      {
        if ((nPairs & 1) != 0)
        {
          nPairs--;
          i3 = nLeft * nPairs;
          i4 = n - i3;
          if (i4 > nLeft)
          {
            merge(idx_data, x_data, i3, nLeft, i4 - nLeft, iwork_data,
                  xwork_data);
          }
        }

        i3 = nLeft * 2;
        nPairs = asr_s32(nPairs, 1U);
        for (i1 = 0; i1 < nPairs; i1++)
        {
          merge(idx_data, x_data, i1 * i3, nLeft, nLeft, iwork_data, xwork_data);
        }

        nLeft = i3;
      }

      if (n > nLeft)
      {
        merge(idx_data, x_data, 0, nLeft, n - nLeft, iwork_data, xwork_data);
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
 * Arguments    : const unsigned int u1[]
 *                int n1
 *                unsigned int y[]
 *                int n
 * Return Type  : void
 */
static void uMultiWord2MultiWord(const unsigned int u1[], int n1, unsigned int
  y[], int n)
{
  int nm;
  if (n1 < n)
  {
    nm = n1;
  }
  else
  {
    nm = n;
  }

  if (0 <= (nm - 1))
  {
    memcpy(&y[0], &u1[0], ((unsigned int)nm) * (sizeof(unsigned int)));
  }

  if ((n > n1) && (nm <= (n - 1)))
  {
    memset(&y[nm], 0, ((unsigned int)((int)(n - nm))) * (sizeof(unsigned int)));
  }
}

/*
 * Arguments    : const unsigned int u1[]
 *                const unsigned int u2[]
 *                int n
 * Return Type  : int
 */
static int uMultiWordCmp(const unsigned int u1[], const unsigned int u2[], int n)
{
  int y;
  int i;
  unsigned int u1i;
  unsigned int u2i;
  y = 0;
  i = n;
  while ((y == 0) && (i > 0))
  {
    i--;
    u1i = u1[i];
    u2i = u2[i];
    if (u1i != u2i)
    {
      if (u1i > u2i)
      {
        y = 1;
      }
      else
      {
        y = -1;
      }
    }
  }

  return y;
}

/*
 * Arguments    : const unsigned int u1[]
 *                const unsigned int u2[]
 *                int n
 * Return Type  : int
 */
static int uMultiWordCmpShr(const unsigned int u1[], const unsigned int u2[],
  int n)
{
  int y;
  unsigned int u2i;
  int i;
  unsigned int u1i;
  unsigned int y1i;
  u2i = 0U;
  y = 0;
  i = n;
  while ((y == 0) && (i > 0))
  {
    i--;
    u1i = u1[i];
    y1i = (u2i << 31U);
    u2i = u2[i];
    y1i |= (u2i >> 1U);
    if (u1i != y1i)
    {
      if (u1i > y1i)
      {
        y = 1;
      }
      else
      {
        y = -1;
      }
    }
  }

  if ((y == 0) && ((u2i & 1U) == 1U))
  {
    y = -1;
  }

  return y;
}

/*
 * Arguments    : unsigned int a[]
 *                int na
 *                unsigned int b[]
 *                int nb
 *                unsigned int q[]
 *                int nq
 *                unsigned int r[]
 *                int nr
 * Return Type  : int
 */
static int uMultiWordDiv(unsigned int a[], int na, unsigned int b[], int nb,
  unsigned int q[], int nq, unsigned int r[], int nr)
{
  int y;
  int nzb;
  int tpi;
  int nza;
  int nb1;
  int na1;
  unsigned int ak;
  unsigned int kbb;
  unsigned int bk;
  unsigned int t;
  unsigned int nbq;
  unsigned int kba;
  unsigned int nba;
  unsigned int nbb;
  unsigned int mask;
  unsigned int kbs;
  int kb;
  unsigned int tnb;
  int ka;
  nzb = nb;
  tpi = nb - 1;
  while ((nzb > 0) && (b[tpi] == 0U))
  {
    nzb--;
    tpi--;
  }

  if (nzb > 0)
  {
    nza = na;
    if (0 <= (nq - 1))
    {
      memset(&q[0], 0, ((unsigned int)nq) * (sizeof(unsigned int)));
    }

    tpi = na - 1;
    while ((nza > 0) && (a[tpi] == 0U))
    {
      nza--;
      tpi--;
    }

    if ((nza > 0) && (nza >= nzb))
    {
      nb1 = nzb - 1;
      na1 = nza - 1;
      if (0 <= (nr - 1))
      {
        memset(&r[0], 0, ((unsigned int)nr) * (sizeof(unsigned int)));
      }

      /* Quick return if dividend and divisor fit into single word. */
      if (nza == 1)
      {
        ak = a[0];
        bk = b[0];
        nbq = ak / bk;
        q[0] = nbq;
        r[0] = ak - (nbq * bk);
        y = 7;
      }
      else
      {
        /* Remove leading zeros from both, dividend and divisor. */
        kbb = 1U;
        t = (b[nb1] >> 1U);
        while (t != 0U)
        {
          kbb++;
          t >>= 1U;
        }

        kba = 1U;
        t = (a[na1] >> 1U);
        while (t != 0U)
        {
          kba++;
          t >>= 1U;
        }

        /* Quick return if quotient is zero. */
        if ((nza > nzb) || (kba >= kbb))
        {
          nba = (((unsigned int)na1) * 32U) + kba;
          nbb = (((unsigned int)nb1) * 32U) + kbb;

          /* Normalize b. */
          if (kbb != 32U)
          {
            bk = b[nb1];
            kbs = 32U - kbb;
            for (kb = nb1; kb > 0; kb--)
            {
              t = (bk << kbs);
              bk = b[kb - 1];
              t |= (bk >> kbb);
              b[kb] = t;
            }

            b[0] = (bk << kbs);
            mask = ~((1U << kbs) - 1U);
          }
          else
          {
            mask = MAX_uint32_T;
          }

          /* Initialize quotient to zero. */
          tnb = 0U;
          y = 0;

          /* Until exit conditions have been met, do */
          do
          {
            /* Normalize a */
            if (kba != 32U)
            {
              kbs = 32U - kba;
              tnb += kbs;
              ak = a[na1];
              for (ka = na1; ka > 0; ka--)
              {
                t = (ak << kbs);
                ak = a[ka - 1];
                t |= (ak >> kba);
                a[ka] = t;
              }

              a[0] = (ak << kbs);
            }

            /* Compare b against the a. */
            ak = a[na1];
            bk = b[nb1];
            if (nb1 == 0)
            {
              t = mask;
            }
            else
            {
              t = MAX_uint32_T;
            }

            if ((ak & t) == bk)
            {
              tpi = 0;
              ka = na1;
              kb = nb1;
              while ((tpi == 0) && (kb > 0))
              {
                ka--;
                ak = a[ka];
                kb--;
                bk = b[kb];
                if (kb == 0)
                {
                  t = mask;
                }
                else
                {
                  t = MAX_uint32_T;
                }

                if ((ak & t) != bk)
                {
                  if (ak > bk)
                  {
                    tpi = 1;
                  }
                  else
                  {
                    tpi = -1;
                  }
                }
              }
            }
            else if (ak > bk)
            {
              tpi = 1;
            }
            else
            {
              tpi = -1;
            }

            /* If the remainder in a is still greater or equal to b, subtract normalized divisor from a. */
            if ((tpi >= 0) || (nba > nbb))
            {
              nbq = nba - nbb;

              /* If the remainder and the divisor are equal, set remainder to zero. */
              if (tpi == 0)
              {
                ka = na1;
                for (kb = nb1; kb > 0; kb--)
                {
                  a[ka] = 0U;
                  ka--;
                }

                a[ka] -= b[0];
              }
              else
              {
                /* Otherwise, subtract the divisor from the remainder */
                if (tpi < 0)
                {
                  ak = a[na1];
                  kba = 31U;
                  for (ka = na1; ka > 0; ka--)
                  {
                    t = (ak << 1U);
                    ak = a[ka - 1];
                    t |= (ak >> 31U);
                    a[ka] = t;
                  }

                  a[0] = (ak << 1U);
                  tnb++;
                  nbq--;
                }

                tpi = 0;
                ka = na1 - nb1;
                for (kb = 0; kb < nzb; kb++)
                {
                  bk = b[kb];
                  t = a[ka];
                  ak = (t - bk) - ((unsigned int)tpi);
                  if (((unsigned int)tpi) != 0U)
                  {
                    tpi = (int)((ak >= t) ? 1 : 0);
                  }
                  else
                  {
                    tpi = (int)((ak > t) ? 1 : 0);
                  }

                  a[ka] = ak;
                  ka++;
                }
              }

              /* Update the quotient. */
              tpi = ((int)nbq) / 32;
              q[tpi] |= (1U << (nbq - (((unsigned int)tpi) * 32U)));

              /* Remove leading zeros from the remainder and check whether the exit conditions have been met. */
              tpi = na1;
              while ((nza > 0) && (a[tpi] == 0U))
              {
                nza--;
                tpi--;
              }

              if (nza >= nzb)
              {
                na1 = nza - 1;
                kba = 1U;
                t = (a[na1] >> 1U);
                while (t != 0U)
                {
                  kba++;
                  t >>= 1U;
                }

                nba = ((((unsigned int)na1) * 32U) + kba) - tnb;
                if (nba < nbb)
                {
                  y = 2;
                }
              }
              else if (nza == 0)
              {
                y = 1;
              }
              else
              {
                na1 = nza - 1;
                y = 4;
              }
            }
            else
            {
              y = 3;
            }
          }
          while (y == 0);

          /* Return the remainder. */
          if (y == 1)
          {
            r[0] = a[0];
          }
          else
          {
            kb = ((int)tnb) / 32;
            nbq = tnb - (((unsigned int)kb) * 32U);
            if (nbq == 0U)
            {
              ka = kb;
              for (tpi = 0; tpi <= nb1; tpi++)
              {
                r[tpi] = a[ka];
                ka++;
              }
            }
            else
            {
              kbs = 32U - nbq;
              ak = a[kb];
              tpi = 0;
              for (ka = kb + 1; ka <= na1; ka++)
              {
                t = (ak >> nbq);
                ak = a[ka];
                t |= (ak << kbs);
                r[tpi] = t;
                tpi++;
              }

              r[tpi] = (ak >> nbq);
            }
          }

          /* Restore b. */
          if (kbb != 32U)
          {
            bk = b[0];
            kbs = 32U - kbb;
            for (kb = 0; kb < nb1; kb++)
            {
              t = (bk >> kbs);
              bk = b[kb + 1];
              t |= (bk << kbb);
              b[kb] = t;
            }

            b[kb] = (bk >> kbs);
          }
        }
        else
        {
          for (tpi = 0; tpi < nr; tpi++)
          {
            r[tpi] = a[tpi];
          }

          y = 6;
        }
      }
    }
    else
    {
      for (tpi = 0; tpi < nr; tpi++)
      {
        r[tpi] = a[tpi];
      }

      y = 5;
    }
  }
  else
  {
    y = -1;
  }

  return y;
}

/*
 * Arguments    : unsigned int y[]
 *                int n
 * Return Type  : void
 */
static void uMultiWordInc(unsigned int y[], int n)
{
  int i;
  unsigned int yi;
  int carry = 1;
  for (i = 0; i < n; i++)
  {
    yi = y[i] + ((unsigned int)carry);
    y[i] = yi;
    carry = (int)((yi < ((unsigned int)carry)) ? 1 : 0);
  }
}

/*
 * Arguments    : const unsigned int u1[]
 *                const unsigned int u2[]
 *                int n
 * Return Type  : bool
 */
static bool uMultiWordLe(const unsigned int u1[], const unsigned int u2[], int n)
{
  return uMultiWordCmp(u1, u2, n) <= 0;
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
static void uMultiWordMul(const unsigned int u1[], int n1, const unsigned int
  u2[], int n2, unsigned int y[], int n)
{
  int i;
  unsigned int cb;
  unsigned int u1i;
  int a1;
  int a0;
  int ni;
  int k;
  int j;
  int b1;
  int b0;
  unsigned int w01;
  unsigned int yk;
  unsigned int t;

  /* Initialize output to zero */
  if (0 <= (n - 1))
  {
    memset(&y[0], 0, ((unsigned int)n) * (sizeof(unsigned int)));
  }

  for (i = 0; i < n1; i++)
  {
    cb = 0U;
    u1i = u1[i];
    a1 = (int)((unsigned int)(u1i >> 16U));
    a0 = (int)((unsigned int)(u1i & 65535U));
    ni = n - i;
    if (n2 <= ni)
    {
      ni = n2;
    }

    k = i;
    for (j = 0; j < ni; j++)
    {
      u1i = u2[j];
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
}

/*
 * Arguments    : const unsigned int u1[]
 *                int n1
 *                unsigned int n2
 *                unsigned int y[]
 *                int n
 * Return Type  : void
 */
static void uMultiWordShr(const unsigned int u1[], int n1, unsigned int n2,
  unsigned int y[], int n)
{
  int nb;
  int i;
  int nc;
  unsigned int nr;
  int i1;
  unsigned int nl;
  unsigned int u1i;
  unsigned int yi;
  nb = ((int)n2) / 32;
  i = 0;
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

      yi = (u1i >> nr);
      if (nc < n1)
      {
        yi |= (u1[nc] << nl);
      }

      y[i] = yi;
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
    y[i] = 0U;
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
static void uMultiWordShrConv(const unsigned int u1[], int n1, unsigned int n2,
  unsigned int y[], int n)
{
  unsigned int n2m1;
  int nb;
  unsigned int maskHalfLSB;
  int i;
  unsigned int maskLSB;
  unsigned int mask;
  bool doRoundUp = false;
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

  uMultiWordShr(u1, n1, n2, y, n);
  i = 0;
  while (doRoundUp && (i < n))
  {
    y[i]++;
    doRoundUp = ((y[i] == 0U) && doRoundUp);
    i++;
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
  int64m_T r;
  int64m_T sumx;
  int midm1;
  int64m_T r1;
  int96m_T r2;
  int64m_T r3;
  int64m_T r4;
  int b_i;
  int96m_T y[101];
  int96m_T a1;
  int64m_T r5;
  int96m_T r6;
  int96m_T r7;
  int64m_T r8;
  static const int96m_T r9 =
  {
    {
      0U, 0U, 0U
    }                                  /* chunks */
  };

  int96m_T r10;
  static const int96m_T r11 =
  {
    {
      MAX_uint32_T, MAX_uint32_T, 1023U
    }                                  /* chunks */
  };

  static const int96m_T r12 =
  {
    {
      0U, 0U, 4294966272U
    }                                  /* chunks */
  };

  int64m_T b_y;
  static const int64m_T r13 =
  {
    {
      0U, 0U
    }                                  /* chunks */
  };

  static const int64m_T r14 =
  {
    {
      1U, 0U
    }                                  /* chunks */
  };

  static const int64m_T r15 =
  {
    {
      MAX_uint32_T, 31U
    }                                  /* chunks */
  };

  int64m_T ytemp;
  int96m_T r16;
  int96m_T r17;
  int64m_T r18;
  int i3;
  int96m_T r19;
  int96m_T r20;
  int96m_T r21;
  int96m_T r22;
  int96m_T r23;
  int96m_T r24;
  int96m_T r25;
  unsigned int u;
  unsigned int u1;
  uint64m_T r26;
  int96m_T r27;
  uint64m_T r28;
  int c_y[101];
  unsigned char locs_data[99];
  int locs_size[1];
  unsigned int fmo_2_data[99];
  int fmo_2_size[1];
  bool s_data[101];
  int trueCount;
  bool b;
  signed char tmp_data[101];
  int x_data[101];
  int64m_T r29;
  int64m_T r30;
  int64m_T r31;
  int64m_T r32;
  uint64m_T r33;
  uint64m_T r34;
  signed char b_tmp_data[101];
  int64m_T hfi;

  AudioSteam = ImpnseIsOut;

  i = (int)ImpnseIsOut.NumChan;
  for (idx = 0; idx < i; idx++)
  {
    b_idx = (unsigned char)((int)((idx + 1) & 3));

    i1 = ((int)b_idx) - 1;
    i2 = 101 * i1;
    sLong2MultiWord(AudioSteam.f[i2], (unsigned int *)(&r.chunks[0U]), 2);
    MultiWordSignedWrap((unsigned int *)(&r.chunks[0U]), 2, 25U, (unsigned int *)
                        (&sumx.chunks[0U]));
    for (midm1 = 0; midm1 < 100; midm1++)
    {
      sLong2MultiWord(AudioSteam.f[(midm1 + i2) + 1], (unsigned int *)
                      (&r1.chunks[0U]), 2);
      MultiWordSignedWrap((unsigned int *)(&r1.chunks[0U]), 2, 25U, (unsigned
        int *)(&r3.chunks[0U]));
      MultiWordAdd((unsigned int *)(&sumx.chunks[0U]), (unsigned int *)
                   (&r3.chunks[0U]), (unsigned int *)(&r.chunks[0U]), 2);
      MultiWordSignedWrap((unsigned int *)(&r.chunks[0U]), 2, 25U, (unsigned int
        *)(&sumx.chunks[0U]));
    }

    uLong2MultiWord(101U, (unsigned int *)(&r.chunks[0U]), 2);
    sMultiWordDivConv((unsigned int *)(&sumx.chunks[0U]), 2, (unsigned int *)
                      (&r.chunks[0U]), 2, (unsigned int *)(&r2.chunks[0U]), 3,
                      (unsigned int *)(&r3.chunks[0U]), 2, (unsigned int *)
                      (&r1.chunks[0U]), 2, (unsigned int *)(&r4.chunks[0U]), 2);
    midm1 = MultiWord2sLong((unsigned int *)(&r2.chunks[0U]));
    for (b_i = 0; b_i < 101; b_i++)
    {
      sLong2MultiWord(AudioSteam.f[b_i + i2], (unsigned int *)(&r4.chunks[0U]),
                      2);
      MultiWordSignedWrap((unsigned int *)(&r4.chunks[0U]), 2, 31U, (unsigned
        int *)(&r1.chunks[0U]));
      sLong2MultiWord(midm1, (unsigned int *)(&r5.chunks[0U]), 2);
      MultiWordSignedWrap((unsigned int *)(&r5.chunks[0U]), 2, 31U, (unsigned
        int *)(&r4.chunks[0U]));
      MultiWordSub((unsigned int *)(&r1.chunks[0U]), (unsigned int *)
                   (&r4.chunks[0U]), (unsigned int *)(&r3.chunks[0U]), 2);
      MultiWordSignedWrap((unsigned int *)(&r3.chunks[0U]), 2, 31U, (unsigned
        int *)(&r8.chunks[0U]));
      sMultiWordMul((unsigned int *)(&r8.chunks[0U]), 2, (unsigned int *)
                    (&r8.chunks[0U]), 2, (unsigned int *)(&y[b_i].chunks[0U]), 3);
    }

    MultiWordSignedWrap((unsigned int *)(&y[0].chunks[0U]), 3, 23U, (unsigned
      int *)(&a1.chunks[0U]));
    for (midm1 = 0; midm1 < 100; midm1++)
    {
      MultiWordSignedWrap((unsigned int *)(&y[midm1 + 1].chunks[0U]), 3, 23U,
                          (unsigned int *)(&r6.chunks[0U]));
      MultiWordAdd((unsigned int *)(&a1.chunks[0U]), (unsigned int *)
                   (&r6.chunks[0U]), (unsigned int *)(&r7.chunks[0U]), 3);
      MultiWordSignedWrap((unsigned int *)(&r7.chunks[0U]), 3, 23U, (unsigned
        int *)(&a1.chunks[0U]));
    }

    r7 = r9;
    if (sMultiWordLt((unsigned int *)(&a1.chunks[0U]), (unsigned int *)
                     (&r9.chunks[0U]), 3))
    {
      r10 = r12;
    }
    else
    {
      r10 = r11;
    }

    b_y = r13;
    r6 = r9;
    if (!sMultiWordLe((unsigned int *)(&r10.chunks[0U]), (unsigned int *)
                      (&r9.chunks[0U]), 3))
    {
      r5 = r14;
      for (b_i = 36; b_i >= 0; b_i--)
      {
        sMultiWordShl((unsigned int *)(&r14.chunks[0U]), 2, (unsigned int)b_i,
                      (unsigned int *)(&r4.chunks[0U]), 2);
        MultiWordSignedWrap((unsigned int *)(&r4.chunks[0U]), 2, 26U, (unsigned
          int *)(&r1.chunks[0U]));
        MultiWordIor((unsigned int *)(&b_y.chunks[0U]), (unsigned int *)
                     (&r1.chunks[0U]), (unsigned int *)(&r3.chunks[0U]), 2);
        MultiWordSignedWrap((unsigned int *)(&r3.chunks[0U]), 2, 26U, (unsigned
          int *)(&ytemp.chunks[0U]));
        sMultiWordMul((unsigned int *)(&ytemp.chunks[0U]), 2, (unsigned int *)
                      (&ytemp.chunks[0U]), 2, (unsigned int *)(&r16.chunks[0U]),
                      3);
        if (sMultiWordLe((unsigned int *)(&r16.chunks[0U]), (unsigned int *)
                         (&r10.chunks[0U]), 3))
        {
          b_y = ytemp;
        }
      }

      r3 = r15;
      if (sMultiWordLt((unsigned int *)(&b_y.chunks[0U]), (unsigned int *)
                       (&r15.chunks[0U]), 2))
      {
        r4 = r14;
        MultiWordAdd((unsigned int *)(&b_y.chunks[0U]), (unsigned int *)
                     (&r14.chunks[0U]), (unsigned int *)(&r1.chunks[0U]), 2);
        MultiWordSignedWrap((unsigned int *)(&r1.chunks[0U]), 2, 26U, (unsigned
          int *)(&ytemp.chunks[0U]));
        sMultiWordMul((unsigned int *)(&ytemp.chunks[0U]), 2, (unsigned int *)
                      (&ytemp.chunks[0U]), 2, (unsigned int *)(&r17.chunks[0U]),
                      3);
        MultiWordSignedWrap((unsigned int *)(&r10.chunks[0U]), 3, 20U, (unsigned
          int *)(&r19.chunks[0U]));
        MultiWordSub((unsigned int *)(&r17.chunks[0U]), (unsigned int *)
                     (&r19.chunks[0U]), (unsigned int *)(&r21.chunks[0U]), 3);
        MultiWordSignedWrap((unsigned int *)(&r21.chunks[0U]), 3, 20U, (unsigned
          int *)(&r23.chunks[0U]));
        MultiWordSignedWrap((unsigned int *)(&r10.chunks[0U]), 3, 20U, (unsigned
          int *)(&r19.chunks[0U]));
        sMultiWordMul((unsigned int *)(&b_y.chunks[0U]), 2, (unsigned int *)
                      (&b_y.chunks[0U]), 2, (unsigned int *)(&r25.chunks[0U]), 3);
        MultiWordSub((unsigned int *)(&r19.chunks[0U]), (unsigned int *)
                     (&r25.chunks[0U]), (unsigned int *)(&r17.chunks[0U]), 3);
        MultiWordSignedWrap((unsigned int *)(&r17.chunks[0U]), 3, 20U, (unsigned
          int *)(&r21.chunks[0U]));
        if (sMultiWordLt((unsigned int *)(&r23.chunks[0U]), (unsigned int *)
                         (&r21.chunks[0U]), 3))
        {
          b_y = ytemp;
        }
      }
    }

    /*  Locate peaks */
    sLong2MultiWord(AudioSteam.f[i2], (unsigned int *)(&r1.chunks[0U]), 2);
    MultiWordSignedWrap((unsigned int *)(&r1.chunks[0U]), 2, 25U, (unsigned int *)
                        (&sumx.chunks[0U]));
    for (midm1 = 0; midm1 < 100; midm1++)
    {
      sLong2MultiWord(AudioSteam.f[(midm1 + i2) + 1], (unsigned int *)
                      (&r5.chunks[0U]), 2);
      MultiWordSignedWrap((unsigned int *)(&r5.chunks[0U]), 2, 25U, (unsigned
        int *)(&r4.chunks[0U]));
      MultiWordAdd((unsigned int *)(&sumx.chunks[0U]), (unsigned int *)
                   (&r4.chunks[0U]), (unsigned int *)(&r1.chunks[0U]), 2);
      MultiWordSignedWrap((unsigned int *)(&r1.chunks[0U]), 2, 25U, (unsigned
        int *)(&sumx.chunks[0U]));
    }

    uLong2MultiWord(101U, (unsigned int *)(&r1.chunks[0U]), 2);
    sMultiWordDivConv((unsigned int *)(&sumx.chunks[0U]), 2, (unsigned int *)
                      (&r1.chunks[0U]), 2, (unsigned int *)(&r17.chunks[0U]), 3,
                      (unsigned int *)(&r4.chunks[0U]), 2, (unsigned int *)
                      (&r5.chunks[0U]), 2, (unsigned int *)(&r18.chunks[0U]), 2);
    i3 = MultiWord2sLong((unsigned int *)(&r17.chunks[0U]));
    sLong2MultiWord(i3, (unsigned int *)(&r20.chunks[0U]), 3);
    sMultiWordShl((unsigned int *)(&r20.chunks[0U]), 3, 34U, (unsigned int *)
                  (&r22.chunks[0U]), 3);
    MultiWordSignedWrap((unsigned int *)(&r22.chunks[0U]), 3, 28U, (unsigned int
      *)(&r24.chunks[0U]));
    sMultiWordShl((unsigned int *)(&b_y.chunks[0U]), 2, 4U, (unsigned int *)
                  (&r4.chunks[0U]), 2);
    u = MultiWord2uLong((unsigned int *)(&r4.chunks[0U]));
    u1 = 2684354560U;
    uMultiWordMul((unsigned int *)(&u), 1, (unsigned int *)(&u1), 1, (unsigned
      int *)(&r26.chunks[0U]), 2);
    uMultiWord2MultiWord((unsigned int *)(&r26.chunks[0U]), 2, (unsigned int *)(
      &r20.chunks[0U]), 3);
    MultiWordSignedWrap((unsigned int *)(&r20.chunks[0U]), 3, 28U, (unsigned int
      *)(&r22.chunks[0U]));
    MultiWordAdd((unsigned int *)(&r24.chunks[0U]), (unsigned int *)
                 (&r22.chunks[0U]), (unsigned int *)(&r27.chunks[0U]), 3);
    MultiWordSignedWrap((unsigned int *)(&r27.chunks[0U]), 3, 28U, (unsigned int
      *)(&r25.chunks[0U]));
    sMultiWordShrConv((unsigned int *)(&r25.chunks[0U]), 3, 32U, (unsigned int *)
                      (&r19.chunks[0U]), 3);
    sMultiWord2MultiWord((unsigned int *)(&r19.chunks[0U]), 3, (unsigned int *)(
      &r28.chunks[0U]), 2);
    MultiWordUnsignedWrap((unsigned int *)(&r28.chunks[0U]), 2, 30U, (unsigned
      int *)(&r26.chunks[0U]));
    AudioSteam.minHeight[i1] = r26;

    /*  min height to be 'peak';  3 std above mean */

    for (midm1 = 0; midm1 < 101; midm1++)
    {
      i3 = midm1 + i2;
      if (AudioSteam.f[i3] < 0)
      {
        c_y[midm1] = -AudioSteam.f[i3];
      }
      else
      {
        c_y[midm1] = AudioSteam.f[i3];
      }

      s_data[midm1] = false;
    }

    detectpeek(c_y, AudioSteam.minHeight[((int)b_idx) - 1], locs_data, locs_size,
               fmo_2_data, fmo_2_size);

    midm1 = locs_size[0];
    for (i3 = 0; i3 < midm1; i3++)
    {
      s_data[((int)locs_data[i3]) - 1] = true;
    }

    /*      s_data = sum(bsxfun(@minus, reshape((1:size(AudioSteam.f(:,idx),1)),1,[]), reshape(locs,[],1)) == 0) == 1; */
    memcpy(&AudioSteam.isOut[i2], &s_data[0], 101U * (sizeof(bool)));

    /*  resize size mismatch */
    trueCount = 0;
    for (b_i = 0; b_i < 101; b_i++)
    {
      b = !AudioSteam.isOut[b_i + i2];
      s_data[b_i] = b;
      if (b)
      {
        trueCount++;
      }
    }

    midm1 = 0;
    for (b_i = 0; b_i < 101; b_i++)
    {
      if (s_data[b_i])
      {
        tmp_data[midm1] = (signed char)(b_i + 1);
        midm1++;
      }
    }

    midm1 = asr_s32(trueCount + ((int)((trueCount < 0) ? 1 : 0)), 1U);
    locs_size[0] = trueCount;
    for (i3 = 0; i3 < trueCount; i3++)
    {
      x_data[i3] = AudioSteam.f[(((int)tmp_data[i3]) + i2) - 1];
    }

    sort(x_data, locs_size);
    if (trueCount == (midm1 * 2))
    {
      sLong2MultiWord(x_data[midm1 - 1], (unsigned int *)(&r29.chunks[0U]), 2);
      MultiWordSignedWrap((unsigned int *)(&r29.chunks[0U]), 2, 31U, (unsigned
        int *)(&r30.chunks[0U]));
      sLong2MultiWord(x_data[midm1], (unsigned int *)(&r32.chunks[0U]), 2);
      MultiWordSignedWrap((unsigned int *)(&r32.chunks[0U]), 2, 31U, (unsigned
        int *)(&r29.chunks[0U]));
      MultiWordAdd((unsigned int *)(&r30.chunks[0U]), (unsigned int *)
                   (&r29.chunks[0U]), (unsigned int *)(&r31.chunks[0U]), 2);
      MultiWordSignedWrap((unsigned int *)(&r31.chunks[0U]), 2, 31U, (unsigned
        int *)(&r18.chunks[0U]));
      sMultiWordShrConv((unsigned int *)(&r18.chunks[0U]), 2, 1U, (unsigned int *)
                        (&r5.chunks[0U]), 2);
      midm1 = MultiWord2sLong((unsigned int *)(&r5.chunks[0U]));
    }
    else
    {
      midm1 = x_data[midm1];
    }

    sLong2MultiWord(midm1, (unsigned int *)(&r31.chunks[0U]), 2);
    sMultiWordShl((unsigned int *)(&r31.chunks[0U]), 2, 7U, (unsigned int *)
                  (&r18.chunks[0U]), 2);
    sMultiWord2MultiWord((unsigned int *)(&r18.chunks[0U]), 2, (unsigned int *)(
      &r28.chunks[0U]), 2);
    MultiWordUnsignedWrap((unsigned int *)(&r28.chunks[0U]), 2, 25U, (unsigned
      int *)(&r33.chunks[0U]));
    AudioSteam.repValue[i1] = r33;

    /*  you could also use NaN, 0, etc.... */
    /*  make a copy */
    trueCount = 0;
    for (b_i = 0; b_i < 101; b_i++)
    {
      i3 = b_i + i2;
      sLong2MultiWord(AudioSteam.f[i3], (unsigned int *)(&r18.chunks[0U]), 2);
      MultiWordSignedWrap((unsigned int *)(&r18.chunks[0U]), 2, 31U, (unsigned
        int *)(&r31.chunks[0U]));
      AudioSteam.fFilt[i3] = r31;
      if (AudioSteam.isOut[i3])
      {
        trueCount++;
      }
    }

    midm1 = 0;
    for (b_i = 0; b_i < 101; b_i++)
    {
      if (AudioSteam.isOut[b_i + i2])
      {
        b_tmp_data[midm1] = (signed char)(b_i + 1);
        midm1++;
      }
    }

    r34 = AudioSteam.repValue[i1];
    uMultiWordShrConv((unsigned int *)(&r34.chunks[0U]), 2, 7U, (unsigned int *)
                      (&r28.chunks[0U]), 2);
    uMultiWord2MultiWord((unsigned int *)(&r28.chunks[0U]), 2, (unsigned int *)(
      &r18.chunks[0U]), 2);
    MultiWordSignedWrap((unsigned int *)(&r18.chunks[0U]), 2, 31U, (unsigned int
      *)(&hfi.chunks[0U]));
    for (i1 = 0; i1 < trueCount; i1++)
    {
      AudioSteam.fFilt[(((int)b_tmp_data[i1]) + i2) - 1] = hfi;
    }

    /*  replace with median value */
    /*      IsOutStruct = AudioSteam; % copy struct params */
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
  static const uint64m_T r =
  {
    {
      0U, 0U
    }                                  /* chunks */
  };

  int i;
  static const int iv[202] =
  {
    25952256, 30474240, 34013184, 26083328, 26804224, 26476544, 1333592064,
    36896768, 32571392, 32964608, 31260672, 37355520, 51445760, 49741824,
    37748736, 29229056, 26607616, 26017792, 26738688, 36569088, 26083328,
    27328512, 27000832, 27066368, 26542080, 30474240, 28377088, 28835840,
    32899072, 28901376, 26607616, 28114944, 38207488, 34275328, 25952256,
    30474240, 34013184, 26083328, 26804224, 26476544, -1333592064, 36896768,
    32571392, 32964608, 31260672, 37355520, 51445760, 49741824, 37748736,
    29229056, 26607616, 26017792, 26738688, 36569088, 26083328, 27328512,
    27000832, 27066368, 26542080, 30474240, 28377088, 28835840, 32899072,
    28901376, 26607616, 28114944, 38207488, 34275328, 25952256, 30474240,
    34013184, 26083328, 26804224, 26476544, 1333592064, 36896768, 32571392,
    32964608, 31260672, 37355520, 51445760, 49741824, 37748736, 29229056,
    26607616, 26017792, 26738688, 36569088, 26083328, 27328512, 27000832,
    27066368, 26542080, 30474240, 28377088, 28835840, 32899072, 28901376,
    26607616, 28114944, 38207488, 25952256, 30474240, 34013184, 26083328,
    26804224, 26476544, -1333592064, 36896768, 32571392, 32964608, 31260672,
    37355520, 51445760, 49741824, 37748736, 29229056, 26607616, 26017792,
    26738688, 36569088, 26083328, 27328512, 27000832, 27066368, 26542080,
    30474240, 28377088, 28835840, 32899072, 28901376, 26607616, 28114944,
    38207488, 34275328, 25952256, 30474240, 34013184, 26083328, 26804224,
    26476544, 1333592064, 36896768, 32571392, 32964608, 31260672, 37355520,
    51445760, 49741824, 37748736, 29229056, 26607616, 26017792, 26738688,
    36569088, 26083328, 27328512, 27000832, 27066368, 26542080, 30474240,
    28377088, 28835840, 32899072, 28901376, 26607616, 28114944, 38207488,
    34275328, 25952256, 30474240, 34013184, 26083328, 26804224, 26476544,
    -1333592064, 36896768, 32571392, 32964608, 31260672, 37355520, 51445760,
    49741824, 37748736, 29229056, 26607616, 26017792, 26738688, 36569088,
    26083328, 27328512, 27000832, 27066368, 26542080, 30474240, 28377088,
    28835840, 32899072, 28901376, 26607616, 28114944, 38207488
  };

  static const int64m_T r1 =
  {
    {
      0U, 0U
    }                                  /* chunks */
  };

  ImpnseIsOut.NumChan = 2U;

  ImpnseIsOut.minHeight[0] = r;
  ImpnseIsOut.repValue[0] = r;
  ImpnseIsOut.minHeight[1] = r;
  ImpnseIsOut.repValue[1] = r;

  for (i = 0; i < 202; i++)
  {
    ImpnseIsOut.f[i] = iv[i];
    ImpnseIsOut.fFilt[i] = r1;
    ImpnseIsOut.isOut[i] = false;
  }

  return ImpnseIsOut;
}

/*
 * File trailer for impnse.c
 *
 * [EOF]
 */
