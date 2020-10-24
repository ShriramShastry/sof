 /*
 * SPDX - License - Identifier: BSD - 3 - Clause
  * File: main.c
 *
 * Author : Sriram Shastry & Ingalsuo, Seppo
 * Coyright(c) 2017 Intel Corporation.All rights reserved.
 */


/* Include Files */
#include <sof/audio/impulse_noise/main.h>
#include <sof/audio/impulse_noise/impnse.h>
#if 0
int impnse_main(int argc, const char * const argv[])
{
  (void)argc;
  (void)argv;

  /* The initialize function is being called automatically from your entry-point function. So, a call to initialize is not included here. */
  /* Invoke the entry-point functions.
     You can call entry-point functions multiple times. */
  #if 0
  struct0_T ImpnseIsOut;
  struct0_T AudioSteam;
  ImpnseIsOut = init_struc_fixpt();         // function initialization - this is  required for wrapper unit testing
  AudioSteam = impnse_fixpt(ImpnseIsOut);   // function initialization - this is  required for impulse noise detection and cancellation
  #endif


  /* Terminate the application.
     You do not need to do this more than one time. */
  impnse_terminate();
  return 0;
}
#endif
/*
 * File trailer for main.c
 *
 * [EOF]
 */
