#include "../test_runner.h"

#define START 0
#define STOP 10
#define NUM_VAL 20

int8_t test_func(size_t index) {
	// custom generators can be used to produce various test ranges, linear, log, exp, etc.
	return func_under_test(log_generator(START, STOP, NUM_VAL, index));
}

int8_t reference_func(size_t index) {
	// e.g. func from libc
	return reference_function(log_generator(START, STOP, NUM_VAL, index));
}

struct test_case test_2 = {
	.type = TEST_DATA_S8;
	.is_hardcoded_data = true;
	.name = "test func 2";
	.num_val = NUM_VAL;
	.funcs = {
		.s8_type = {
			.test_func = test_func
			.reference_func = reference_func;
		};
	};
}

// might be manual register initially
REGISTER_TEST(test_2);
