void float_runner(struct test_case *test) {
	struct float_functions *f_func = &test->funcs.float_type;
	for (size_t i = 0; i < test->num_val; i++) {
		// dimmy diff check
		if (abs(f_func->test_func(i) - f_func->reference_func(i)) > 0.5) {
			// print error
		}
		// log to file?
	}
}

void s8_runner(struct test_case *test) {
	// similar to above but types are adjusted
}

// common utility functions here like writing to csv file

int main(int argc, char *argv[]) {
	// parse flags
	// e.g. flags
	// -n test name
	// register tests
	// lookup test based on name
	switch (type_of_test) {
	case TEST_TYPE_FLOAT:
		float_runner(test);
		break;
	case TEST_TYPE_S8:
		break;
	case TEST_TYPE_S16:
		break;
	case TEST_TYPE_S24:
		break;
	case TEST_TYPE_S32:
		break;
	}
	// differ to test runner type
}
