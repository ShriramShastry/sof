#include "../test_runner.h"
#include <string.h>
#include <stdint.h>
#include <math.h>

// can add a method to init data if we want to store in a separate file
#define START 0
#define STOP 10
#define NUM_VAL 20

static  float data[] = {
	(float)0.0,
	(float)0.1
};

static float reference_data[] = {
	(float)0.43,
	(float)0.32
};

static int32_t data_32[] = {
	(int32_t)2147483648,
	(int32_t)2147483647
};

static int32_t reference_data_32[] = {
	(int32_t)2147483648,
	(int32_t)2147483647
};

static int16_t data_16[] = {
	(int16_t)65536,
	(int16_t)65535
};

static int16_t reference_data_16[] = {
	(int16_t)65536,
	(int16_t)65535
};


static int8_t data_8[] = {
	(int8_t)256,
	(int8_t)255
};

static int8_t reference_data_8[] = {
	(int8_t)256,
	(int8_t)255
};

float test_func(size_t index) {
	return func_under_test(data[index]);
}

float reference_func(size_t index) {
	return reference_data[index];
}

struct test_case test_1 = {
	.type = TEST_DATA_FLOAT,
	.is_hardcoded_data = true,
	.name = "test func 1",
	.num_val = NUM_VAL, 
	.funcs = {
		.float_type = {
			.test_func = test_func,
			.reference_func = reference_func,
		}
	}
};


/*****************************************************************************
*
* FUNCTION     :   UnwrapPhase
*
* PROTOTYPE    :   void UnwrapPhase(Float32 *pfa_Array, int16 s_pts);
*
* DESCRIPTION  :   Unwraps an array of phases (in radians) assumed to be between 0 and 2*Pi
*                  so that relative jumps between values are no bigger than Pi
*
*
* PARAMETERS   :   r0 = pfa_Array , a0 = s_pts
*
*
* RETURN-VALUE :
*
*
* NOTES        :   Similar to Matlab's unwrap function
*****************************************************************************/
// might be manual register initially
REGISTER_TEST(test_1);