#include "../test_runner.h"

// can add a method to init data if we want to store in a separate file

static float data[] = {
	0.0,
	0.1
};

static float reference_data[] = {
	0.43,
	0.32
};

float test_func(size_t index) {
	return func_under_test(data[index]);
}

float reference_func(size_t index) {
	return reference_data[index];
}

struct test_case test_1 = {
	.type = TEST_DATA_FLOAT;
	.is_hardcoded_data = true;
	.name = "test func 1";
	.num_val = ARRAY_SIZE(data);
	.funcs = {
		.float_type = {
			.test_func = test_func
			.reference_func = reference_func;
		};
	};
}

// might be manual register initially
REGISTER_TEST(test_1);
