#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#ifndef __cplusplus
#include <stdbool.h>
#endif

enum test_data_type {
	TEST_DATA_FLOAT,
	TEST_DATA_S8,
	TEST_DATA_S16,
	TEST_DATA_S24,
	TEST_DATA_S32,
	// doubles, unsigned, etc.
};

struct float_functions {
	float(*test_func)(size_t index);
	float(*reference_func)(size_t index);
};

struct s32_functions {
	int32_t(*test_func)(size_t index);
	int32_t(*reference_func)(size_t index);
};


struct s16_functions {
	int16_t(*test_func)(size_t index);
	int16_t(*reference_func)(size_t index);
};

struct s8_functions {
	int8_t(*test_func)(size_t index);
	int8_t(*reference_func)(size_t index);
};
struct test_case {
	enum test_data_type type;
	bool  is_hardcoded_data;
	const char* name;
	size_t num_val;

	union {
		// pass in void so we can pass in multdimentational data if we want
		struct float_functions float_type;
		struct s32_functions s32_type;
		struct s16_functions s16_type;
		struct s8_functions s8_type;
	} funcs;
};
