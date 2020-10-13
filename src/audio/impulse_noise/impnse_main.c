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

/* Function Declarations */
static void argInit_101x2_boolean_T(bool result[202]);
static void argInit_101x2_sfix16_En15(short result[202]);
static void argInit_101x2_sfix17_En15(int result[202]);
static void argInit_1x2_ufix18_En17(unsigned int result[2]);
static bool argInit_boolean_T(void);
static short argInit_sfix16_En15(void);
static int argInit_sfix17_En15(void);
static void argInit_struct0_T(struct0_T *result);
static unsigned int argInit_ufix18_En17(void);
static unsigned char argInit_ufix2(void);
static void main_impnse_fixpt(void);
static void main_init_struc_fixpt(void);

/* Function Definitions */

/*
 * Arguments    : bool result[202]
 * Return Type  : void
 */
static void argInit_101x2_boolean_T(bool result[202])
{
  int idx0;
  bool result_tmp;

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 101; idx0++)
  {
    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result_tmp = argInit_boolean_T();
    result[idx0] = result_tmp;

    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result[idx0 + 101] = result_tmp;
  }
}

/*
 * Arguments    : short result[202]
 * Return Type  : void
 */
static void argInit_101x2_sfix16_En15(short result[202])
{
  int idx0;
  short result_tmp;

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 101; idx0++)
  {
    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result_tmp = argInit_sfix16_En15();
    result[idx0] = result_tmp;

    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result[idx0 + 101] = result_tmp;
  }
}

/*
 * Arguments    : int result[202]
 * Return Type  : void
 */
static void argInit_101x2_sfix17_En15(int result[202])
{
  int idx0;
  int result_tmp;

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 101; idx0++)
  {
    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result_tmp = argInit_sfix17_En15();
    result[idx0] = result_tmp;

    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result[idx0 + 101] = result_tmp;
  }
}

/*
 * Arguments    : unsigned int result[2]
 * Return Type  : void
 */
static void argInit_1x2_ufix18_En17(unsigned int result[2])
{
  unsigned int result_tmp;

  /* Loop over the array to initialize each element. */
  /* Set the value of the array element.
     Change this value to the value that the application requires. */
  result_tmp = argInit_ufix18_En17();
  result[0] = result_tmp;

  /* Set the value of the array element.
     Change this value to the value that the application requires. */
  result[1] = result_tmp;
}

/*
 * Arguments    : void
 * Return Type  : bool
 */
static bool argInit_boolean_T(void)
{
  return false;
}

/*
 * Arguments    : void
 * Return Type  : short
 */
static short argInit_sfix16_En15(void)
{
  return 0;
}

/*
 * Arguments    : void
 * Return Type  : int
 */
static int argInit_sfix17_En15(void)
{
  return 0;
}

/*
 * Arguments    : struct0_T *result
 * Return Type  : void
 */
static void argInit_struct0_T(struct0_T *result)
{
  /* Set the value of each structure field.
     Change this value to the value that the application requires. */
  result->NumChan = argInit_ufix2();
  argInit_101x2_sfix16_En15(result->f);
  argInit_1x2_ufix18_En17(result->minHeight);
  argInit_1x2_ufix18_En17(result->repValue);
  argInit_101x2_sfix17_En15(result->fFilt);
  argInit_101x2_boolean_T(result->isOut);
}

/*
 * Arguments    : void
 * Return Type  : unsigned int
 */
static unsigned int argInit_ufix18_En17(void)
{
  return 0U;
}

/*
 * Arguments    : void
 * Return Type  : unsigned char
 */
static unsigned char argInit_ufix2(void)
{
  return 0U;
}

/*
 * Arguments    : void
 * Return Type  : void
 */
static void main_impnse_fixpt(void)
{
  struct0_T r;
  struct0_T AudioSteam;

  /* Initialize function 'impnse_fixpt' input arguments. */
  /* Initialize function input argument 'ImpnseIsOut'. */
  /* Call the entry-point 'impnse_fixpt'. */
  argInit_struct0_T(&r);
  AudioSteam = impnse_fixpt(r);
}

/*
 * Arguments    : void
 * Return Type  : void
 */
static void main_init_struc_fixpt(void)
{
  struct0_T ImpnseIsOut;

  /* Call the entry-point 'init_struc_fixpt'. */
  ImpnseIsOut = init_struc_fixpt();
}

/*
 * Arguments    : int argc
 *                const char * const argv[]
 * Return Type  : int
 */
int impnse_main(int argc, const char * const argv[])
{
  (void)argc;
  (void)argv;

  /* The initialize function is being called automatically from your entry-point function. So, a call to initialize is not included here. */
  /* Invoke the entry-point functions.
     You can call entry-point functions multiple times. */
  main_impnse_fixpt();                      // function initialization - this is not required
  main_init_struc_fixpt();                  // function initialization - this is not required
  struct0_T ImpnseIsOut;
  struct0_T AudioSteam;
  ImpnseIsOut = init_struc_fixpt();         // function initialization - this is  required
  AudioSteam = impnse_fixpt(ImpnseIsOut);   // function initialization - this is  required



  /* Terminate the application.
     You do not need to do this more than one time. */
  impnse_terminate();
  return 0;
}

/*
 * File trailer for main.c
 *
 * [EOF]
 */
