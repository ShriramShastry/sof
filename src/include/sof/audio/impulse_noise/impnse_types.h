 /*
 * SPDX - License - Identifier: BSD - 3 - Clause
 * File: impnse_types.h
 *
 * Author : Sriram Shastry & Ingalsuo, Seppo
 * Coyright(c) 2017 Intel Corporation.All rights reserved.
 */

#ifndef IMPNSE_TYPES_H
#define IMPNSE_TYPES_H

/* Include Files */
#include <sof/audio/impulse_noise/rtwtypes.h>

/* Type Definitions */
#ifndef typedef_int64m_T
#define typedef_int64m_T

typedef struct
{
  unsigned int chunks[2];
}
int64m_T;

#endif                                 /*typedef_int64m_T*/

#ifndef typedef_struct0_T
#define typedef_struct0_T

typedef struct
{
  unsigned char NumChan;
  short f[202];
  unsigned int minHeight[2];
  unsigned int repValue[2];
  int fFilt[202];
  bool isOut[202];
}
struct0_T;

#endif                                 /*typedef_struct0_T*/
#endif

/*
 * File trailer for impnse_types.h
 *
 * [EOF]
 */
