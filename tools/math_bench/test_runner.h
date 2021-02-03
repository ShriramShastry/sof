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
}

struct s16_functions {
	int16_t(*test_func)(size_t index);
	int16_t(*reference_func)(size_t index);
}

struct test_case {
	enum test_data_type type;
	bool is_hardcoded_data;
	const char *name;
	size_t num_val;
	union {
		// pass in void so we can pass in multdimentational data if we want
		struct float_functions float_type;
		struct s16_functions s16_type;
	} funcs;
}
