#include <stdlite/helpers.h>

#include "tests_suite.h"
#include "user_input.h"

void setUp() {
    stdl_set_debug_level(-1);
}

void test_user_input_ok() {
    stdl_user_input* inp = NULL;
    stdl_user_input_new(&inp);
    TEST_ASSERT_NOT_NULL(inp);

    ASSERT_STDL_OK(stdl_user_input_fill_from_toml(inp, "../app/tests/test_files/test.toml"));

    TEST_ASSERT_EQUAL_STRING("test calculation", inp->title);
    TEST_ASSERT_EQUAL_STRING("water_sto3g.fchk", inp->ctx_source_path);
    TEST_ASSERT_EQUAL_STRING("context.h5", inp->ctx_output);
    TEST_ASSERT_FLOAT_WITHIN(1e-4, 12.f / STDL_CONST_AU_TO_EV, inp->ctx_ethr);

    ASSERT_STDL_OK(stdl_user_input_delete(inp));
}


void test_parse_frequency_ok() {
    double value;

    ASSERT_STDL_OK(stdl_user_input_parse_frequency("1.25au", &value));
    TEST_ASSERT_DOUBLE_WITHIN(1e-6, 1.25, value);

    ASSERT_STDL_OK(stdl_user_input_parse_frequency("1.25eV", &value));
    TEST_ASSERT_DOUBLE_WITHIN(1e-6, 1.25 / STDL_CONST_AU_TO_EV, value);

    ASSERT_STDL_OK(stdl_user_input_parse_frequency("125nm", &value));
    TEST_ASSERT_DOUBLE_WITHIN(1e-6, STDL_CONST_HC / 125, value);

    ASSERT_STDL_KO(stdl_user_input_parse_frequency("1.25 au", &value));
    ASSERT_STDL_KO(stdl_user_input_parse_frequency("au", &value));
    ASSERT_STDL_KO(stdl_user_input_parse_frequency("1.25aux", &value));
}