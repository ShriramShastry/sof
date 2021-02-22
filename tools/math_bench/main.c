#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "test_runner.h"
#include "common/const.h"
#include <errno.h>

static void float_runner(struct test_case *test)
{
	struct float_functions *f_func = &test->funcs.float_type;
	for (size_t i = 0; i < test->num_val; i++) {
		// dimmy diff check
		if (abs(f_func->test_func(i) - f_func->reference_func(i)) > 0.5) {
			// print error
		}
		// log to file?
	}
}

static void s32_runner(struct test_case *test)
{
	struct s32_functions* f_func = &test->funcs.s32_type;
	for (size_t i = 0; i < test->num_val; i++) {
		// dimmy diff check
		if (abs(f_func->test_func(i) - f_func->reference_func(i)) > 0.5) {
			// print error
		}
		// log to file?
	}
}

static void s16_runner(struct test_case *test)
{
	struct s16_functions* f_func = &test->funcs.s16_type;
	for (size_t i = 0; i < test->num_val; i++) {
		// dimmy diff check
		if (abs(f_func->test_func(i) - f_func->reference_func(i)) > 0.5) {
			// print error
		}
		// log to file?
	}
}

static void s8_runner(struct test_case *test)
{
	struct s8_functions* f_func = &test->funcs.s8_type;
	for (size_t i = 0; i < test->num_val; i++) {
		// dimmy diff check
		if (abs(f_func->test_func(i) - f_func->reference_func(i)) > 0.5) {
			// print error
		}
		// log to file?
	}
}

void (*TB_Handler_states[]) (void) = { &float_runner,
							  &s8_runner,
							  &s16_runner,
							  &s32_runner
};

// common utility functions here like writing to csv file

int main(int argc, char *argv[])
{
	// parse flags
	// e.g. flags
	// -n test name
	// register tests
	// lookup test based on name
	int type_of_test = 1, test = 0; /* This should be global*/

	(*TB_Handler_states[type_of_test - TEST_TYPE_FLOAT])();

	switch (type_of_test) {
	case TEST_TYPE_FLOAT:
		float_runner(test);
		break;
	case TEST_TYPE_S8:
		s8_runner(test);
		break;
	case TEST_TYPE_S16:
		s16_runner(test);
		break;
	case TEST_TYPE_S24:
		s24_runner(test);
		break;
	case TEST_TYPE_S32:
		s32_runner(test);
		break;
	}
	// differ to test runner type
}
