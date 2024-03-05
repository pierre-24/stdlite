#include "tests_suite.h"
#include "user_input.h"

void setUp() {
    stdl_set_debug_level(3);
}

void test_user_input_ok() {
    stdl_user_input* inp = NULL;
    stdl_user_input_new(&inp);
    TEST_ASSERT_NOT_NULL(inp);

    ASSERT_STDL_OK(stdl_user_input_fill_from_toml(inp, "../app/tests/test_files/test.toml"));

    TEST_ASSERT_EQUAL_STRING("test calculation", inp->title);
    TEST_ASSERT_EQUAL_STRING("water_sto3g.fchk", inp->ctx_source_path);
    TEST_ASSERT_EQUAL_STRING("context.h5", inp->ctx_output);

    ASSERT_STDL_OK(stdl_user_input_delete(inp));
}
