//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright(c) 2021 Intel Corporation. All rights reserved.
//
//  Author: Shriram Shastry <malladi.sastry@linux.intel.com>
//
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "trig.h"

/* function call to arc cosine handler and return the next arc cosine state */
cordicstate cordic_acosine_handler(void)
{
	int32_t a[21];
	int32_t th_rad_fxp[21];
	init_acos_fixpt(a);
	acos_fixpt(a, th_rad_fxp);

	return cordic_cosine_state;
}
/* function call to cosine handler and return the next cosine state */
cordicstate cordic_cosine_handler(void)
{
	int32_t cdcCoSinTh[2047];
	int32_t thRadFxp[2047];
	init_cosine_fixpt(thRadFxp);
	cosine_fixpt(thRadFxp, cdcCoSinTh);

	return cordic_sincos_state;
}
/* function call to sine cosine handler and return the next sine cosine state */
cordicstate cordic_sincos_handler(void)
{
	int32_t theta[2047];
	int32_t x[2047];
	int32_t y[2047];
	init_sincos_fixpt(theta);
	sincos_fixpt(theta, y, x);

	return cordic_angle_state;
}
/* function call to angle handler and return the next angle state */
cordicstate cordic_angle_handler(void)
{
	cint32_T dblRandomVals[20];
	int32_t theta_dbl_cdc[20];
	init_angle_fixpt(dblRandomVals);
	angle_fixpt(dblRandomVals, theta_dbl_cdc);

	return cordic_asine_state;
}
/* function call to arc sine handler and return the next arc sine state */
cordicstate cordic_asine_handler(void)
{
	int32_t a[21];
	int32_t th_rad_fxp[21];
	init_asin_fixpt(a);
	asin_fixpt(a, th_rad_fxp);

	return cordic_atan_state;
}
/* function call to arc tangent handler and return the next arc tangent state */
cordicstate cordic_atan_handler(void)
{
	int32_t thetax[359];
	int32_t thetay[359];
	uint16_t z[359];
	int32_t cdcatan2Th[359];
	init_atan2_fixpt(thetax, thetay, z);
	atan2_fixpt(thetay, thetax, cdcatan2Th);

	return cordic_cexp_state;
}
/* function call to cmpx exponential handler and return the next cmpx exponential state */
cordicstate cordic_cexp_handler(void)
{
	cint32_T cdcCoSinTh[2047];
	int32_t thRadFxp[2047];
	init_cmpx_exp_fixpt(thRadFxp);
	cmpx_exp_fixpt(thRadFxp, cdcCoSinTh);

	return cordic_sine_state;
}
/* function call to sine handler and return the next sine state */
cordicstate cordic_sine_handler(void)
{
	int32_t cdcCoSinTh[2047];
	int32_t thRadFxp[2047];
	init_sine_fixpt(thRadFxp);
	sine_fixpt(thRadFxp, cdcCoSinTh);

	return cordic_acosine_state;
}

int main(int argc, char *argv[])
{
	cordicstate cnextstate = cordic_acosine_state;
	/*cordic_event cnewevent;*/
	int32_t idx = 0;
	/* Lookup table to define valid states and event of finite state machine */
	static tb_cdc_eventhandler finite_state_machine =
	{
	/**
	* finite_state_machine	   |	     event	     |		function	 |
	* transition from current state to next state
	*/
	[cordic_acosine_state]	= {[cordic_acosine_event]    = cordic_acosine_handler },
	    [cordic_cosine_state]   = {[cordic_cosine_event]	 = cordic_cosine_handler },
	[cordic_sincos_state]	= {[cordic_sincos_event]     = cordic_sincos_handler },
	[cordic_angle_state]	= {[cordic_angle_event]	     = cordic_angle_handler },
	[cordic_asine_state]	= {[cordic_asine_event]	     = cordic_asine_handler },
	[cordic_atan_state]	= {[cordic_atan_event]	     = cordic_atan_handler },
	[cordic_cexp_state]	= {[cordic_cexp_event]	     = cordic_cexp_handler },
	[cordic_sine_state]	= {[cordic_sine_event]	     = cordic_sine_handler },
	};
	while(1)
	{
	/*  assuming cmocka api to read the next event */
	cordic_event cnewevent = read_fsm_handler(idx); 
	/*  To be done error */
	/*Description- exceptionalhandler or enum type*/
	cordic_event cnewerr = fsmreaderr(idx);	 
	
	if( ( cnextstate < cordic_last_state) && (cnewevent < cordic_last_event) && finite_state_machine[cnextstate][cnewevent]!= NULL)
	{
	    
	    cnextstate = (*finite_state_machine[cnextstate][cnewevent])();
	}
	else if ((cnextstate >= cordic_last_state) || (cnewevent >= cordic_last_event))
	{
		break; /* Invalid, we are done!  */
	}
		idx++; /* increment state for testing purpose */
	}
	return 0;
}
/* testing purpose incrementing the counter */
int32_t read_fsm_handler(int32_t idx)
{
    return (idx);
}
/*  To be implemented - hardcoded for now */
int32_t fsmreaderr(int32_t idx)
{
    enum  ErrNo;
    return CORDIC_ACOSINE;
}

