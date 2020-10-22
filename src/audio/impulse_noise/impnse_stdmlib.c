/*
 * SPDX-License-Identifier: BSD-3-Clause
 * File: impnse_stdmlib.c
 *
 * Author :Sriram Shastry & Ingalsuo, Seppo
 * Coyright(c) 2017 Intel Corporation. All rights reserved.
 */


/* Include Files */
#include <stddef.h>
#include <stdlib.h>
//#include <sof/audio/impulse_noise/rtwtypes.h>
//#include <sof/audio/impulse_noise/stdmlib.h>

#include <sof/audio/impulse_noise/rtwtypes.h>
#include <sof/audio/impulse_noise/stdmlib.h>



 static const double huge = 1.0e300;
 double  impnse_trunc(double x)
 {
 	int32_t i0,i1,jj0;
 	uint32_t i;
 	EXTRACT_WORDS(i0,i1,x);
 	jj0 = ((i0>>20)&0x7ff)-0x3ff;
 	if(jj0<20) {
 	    if(jj0<0) { 	/* raise inexact if x != 0 */
 		if(huge+x>0.0) {/* |x|<1, so return 0*sign(x) */
 		    i0 &= 0x80000000U;
 		    i1 = 0;
 		}
 	    } else {
 		i = (0x000fffff)>>jj0;
 		if(((i0&i)|i1)==0) return x; /* x is integral */
 		if(huge+x>0.0) {	/* raise inexact flag */
 		    i0 &= (~i); i1=0;
 		}
 	    }
 	} else if (jj0>51) {
 	    if(jj0==0x400) return x+x;	/* inf or NaN */
 	    else return x;		/* x is integral */
 	} else {
 	    i = ((uint32_t)(0xffffffff))>>(jj0-20);
 	    if((i1&i)==0) return x;	/* x is integral */
 	    if(huge+x>0.0)		/* raise inexact flag */
 		i1 &= (~i);
 	}
 	INSERT_WORDS(x,i0,i1);
 	return x;
 }


 /*
  * impnse_floor(x) := the largest integer no larger than x;
  *
  * Note: Inexact will be signaled if x is not an integer, as is
  *	customary for IEEE 754.  No other signal can be emitted.
  */

 double impnse_floor(double x)
 {
 	int32_t i0,i1,jj0;
 	uint32_t i,j;
 	EXTRACT_WORDS(i0,i1,x);
 	jj0 = ((i0>>20)&0x7ff)-0x3ff;

 	if(jj0<20) {
 	    if(jj0<0) { 	/* raise inexact if x != 0 */
 		if(huge+x>0.0) {/* return 0*sign(x) if |x|<1 */
 		    if(i0>=0) {i0=i1=0;}
 		    else if(((i0&0x7fffffff)|i1)!=0)
 			{ i0=0xbff00000;i1=0;}
 		}
 	    } else {
 		i = (0x000fffff)>>jj0;
 		if(((i0&i)|i1)==0) return x; /* x is integral */
 		if(huge+x>0.0) {	/* raise inexact flag */
 		    if(i0<0) i0 += (0x00100000)>>jj0;
 		    i0 &= (~i); i1=0;
 		}
 	    }
 	} else if (jj0>51) {
 	    if(jj0==0x400) return x+x;	/* inf or NaN */
 	    else return x;		/* x is integral */
 	} else {
 	    i = ((uint32_t)(0xffffffff))>>(jj0-20);
 	    if((i1&i)==0) return x;	/* x is integral */
 	    if(huge+x>0.0) { 		/* raise inexact flag */
 		if(i0<0) {
 		    if(jj0==20) i0+=1;
 		    else {
 			j = i1+(1<<(52-jj0));
 			if(j<(uint32_t)i1) i0 +=1 ; 	/* got a carry */
 			i1=j;
 		    }
 		}
 		i1 &= (~i);
 	    }
 	}
   
 	INSERT_WORDS(x,i0,i1);
 	return x;
 }

 
 
 /*
*-------------------------------------------------------------------------------
*
*  Prototype: double impnse_mod(double s_num, double s_den)
*
*  This function returns s_num mod s_den. Inputs are assumed to be positive,
*  and s_den should not be zero.
*
*  Input Arguments:
*     s_num: 16-bit number
*     s_den: 16-bit number
*
*  Output Arguments:
*
*  Returns:
*     s_num mod s_den
*
*  Global Variables:
*
*-------------------------------------------------------------------------------
*/

 double impnse_mod(double s_num, double s_den)
 {
     if (s_den == 1)
         return (int32_T)0;

     if (s_den == 0) // error
         return (int32_T)(0x7FFFFFFF);

     while (s_num >= s_den)
         s_num -= s_den;

     return (int32_T)s_num;
 }

 /*
 * File trailer for impnse_stdmlib.c
 *
 * [EOF]
 */
