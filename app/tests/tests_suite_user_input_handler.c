#include <stdlite/helpers.h>
#include <unistd.h>

#include "tests_suite.h"
#include "user_input_handler.h"

FILE* stream;

void setUp(void) {
    stdl_set_debug_level(-1);

    stream = tmpfile();
}

void tearDown(void) {
    fclose(stream);
}


void test_user_input_context_fill_from_toml_ok() {
    stdl_user_input_handler* inp = NULL;
    stdl_user_input_handler_new(&inp);
    TEST_ASSERT_NOT_NULL(inp);

    fputs("title = \"test calculation\"\n"
          "data_output = \"water_631g.h5\"\n"
          "[context]\n"
          "source = \"../tests/test_files/water_631g.fchk\"\n"
          "source_type = \"FCHK\"\n"
          "tda = 0\n"
          "ethr = '12eV'\n"
          "e2thr=1e-3\n"
          "ax = 1.0\n"
          "gammaJ=1.0\n"
          "gammaK=0.5",
          stream);
    rewind(stream);

    ASSERT_STDL_OK(stdl_user_input_handler_fill_from_toml(inp, stream));
    ASSERT_STDL_OK(stdl_user_input_handler_check(inp));

    TEST_ASSERT_EQUAL_STRING("test calculation", inp->title);
    TEST_ASSERT_EQUAL_STRING("water_631g.h5", inp->data_output);
    TEST_ASSERT_EQUAL_STRING("../tests/test_files/water_631g.fchk", inp->ctx_source);
    TEST_ASSERT_EQUAL(STDL_SRC_FCHK, inp->ctx_source_type);
    TEST_ASSERT_EQUAL(STDL_METHOD_MONOPOLE, inp->ctx_method);
    TEST_ASSERT_EQUAL(0, inp->ctx_tda);
    TEST_ASSERT_FLOAT_WITHIN(1e-4, 12.f / STDL_CONST_AU_TO_EV, inp->ctx_ethr);
    TEST_ASSERT_FLOAT_WITHIN(1e-4, 1e-3, inp->ctx_e2thr);
    TEST_ASSERT_FLOAT_WITHIN(1e-4, 1.0, inp->ctx_ax);
    TEST_ASSERT_FLOAT_WITHIN(1e-4, 1.0, inp->ctx_gammaJ);
    TEST_ASSERT_FLOAT_WITHIN(1e-4, 0.5, inp->ctx_gammaK);

    ASSERT_STDL_OK(stdl_user_input_handler_delete(inp));
}

void test_user_input_context_fill_from_args_ok() {
    stdl_user_input_handler* inp = NULL;
    stdl_user_input_handler_new(&inp);
    TEST_ASSERT_NOT_NULL(inp);

    char* args[] =  {
            "self",
            "--data_output=test.h5",
            "--ctx_source=../tests/test_files/water_631g.molden",
            "--ctx_source_type=MOLDEN",
            "--ctx_gammaJ=0.5",
            "--ctx_gammaK=1.0",
            "--ctx_ethr=12eV",
            "--ctx_e2thr=1e-3",
            "--ctx_ax=1.0",
            "--ctx_tda=0",
    };

    ASSERT_STDL_OK(stdl_user_input_handler_fill_from_args(inp, sizeof(args) / sizeof(char *), args));
    ASSERT_STDL_OK(stdl_user_input_handler_check(inp));

    TEST_ASSERT_EQUAL_STRING("test.h5", inp->data_output);
    TEST_ASSERT_EQUAL_STRING("../tests/test_files/water_631g.molden", inp->ctx_source);
    TEST_ASSERT_EQUAL(STDL_SRC_MOLDEN, inp->ctx_source_type);
    TEST_ASSERT_FLOAT_WITHIN(1e-4, 0.5, inp->ctx_gammaJ);
    TEST_ASSERT_FLOAT_WITHIN(1e-4, 1.0, inp->ctx_gammaK);
    TEST_ASSERT_FLOAT_WITHIN(1e-4, 12.f / STDL_CONST_AU_TO_EV, inp->ctx_ethr);
    TEST_ASSERT_FLOAT_WITHIN(1e-4, 1e-3, inp->ctx_e2thr);
    TEST_ASSERT_FLOAT_WITHIN(1e-4, 1.0, inp->ctx_ax);
    TEST_ASSERT_EQUAL(0, inp->ctx_tda);

    ASSERT_STDL_OK(stdl_user_input_handler_delete(inp));
}

void test_user_input_context_fill_from_args_and_file_ok() {
    stdl_user_input_handler* inp = NULL;
    stdl_user_input_handler_new(&inp);
    TEST_ASSERT_NOT_NULL(inp);

    char* args[] =  {
            "self",
            "../app/tests/test_files/test.toml",
            "--ctx_e2thr=1e-3",
    };

    ASSERT_STDL_OK(stdl_user_input_handler_fill_from_args(inp, sizeof(args) / sizeof(char *), args));
    ASSERT_STDL_OK(stdl_user_input_handler_check(inp));

    TEST_ASSERT_EQUAL_STRING("test calculation", inp->title);
    TEST_ASSERT_FLOAT_WITHIN(1e-4, 1e-3, inp->ctx_e2thr);

    ASSERT_STDL_OK(stdl_user_input_handler_delete(inp));
}

void test_user_input_incorrect_toml_ko() {
    stdl_user_input_handler* inp = NULL;
    stdl_user_input_handler_new(&inp);
    TEST_ASSERT_NOT_NULL(inp);

    fputs("[context]\n"
          "source = \"whatever\"\n"
          "source_type = \"whatever\"", // incorrect source type
          stream);
    rewind(stream);

    ASSERT_STDL_KO(stdl_user_input_handler_fill_from_toml(inp, stream));
    ASSERT_STDL_OK(stdl_user_input_handler_delete(inp));
}

void test_user_input_incorrect_check_args_ko() {
    stdl_user_input_handler* inp = NULL;
    stdl_user_input_handler_new(&inp);
    TEST_ASSERT_NOT_NULL(inp);

    fputs("[context]\n"
          "source = \"test.fchk\"\n"
          "ax = -1",  // ax should be in [0,1]
          stream);
    rewind(stream);

    ASSERT_STDL_OK(stdl_user_input_handler_fill_from_toml(inp, stream));
    ASSERT_STDL_KO(stdl_user_input_handler_check(inp));

    ASSERT_STDL_OK(stdl_user_input_handler_delete(inp));
}


void test_parse_frequency_ok() {
    double value;

    ASSERT_STDL_OK(stdl_user_input_handler_parse_frequency("1.25au", &value));
    TEST_ASSERT_DOUBLE_WITHIN(1e-6, 1.25, value);

    ASSERT_STDL_OK(stdl_user_input_handler_parse_frequency("1.25", &value));
    TEST_ASSERT_DOUBLE_WITHIN(1e-6, 1.25, value);

    ASSERT_STDL_OK(stdl_user_input_handler_parse_frequency("1.25eV", &value));
    TEST_ASSERT_DOUBLE_WITHIN(1e-6, 1.25 / STDL_CONST_AU_TO_EV, value);

    ASSERT_STDL_OK(stdl_user_input_handler_parse_frequency("125nm", &value));
    TEST_ASSERT_DOUBLE_WITHIN(1e-6, STDL_CONST_HC / 125, value);

    ASSERT_STDL_KO(stdl_user_input_handler_parse_frequency("1.25 au", &value));
    ASSERT_STDL_KO(stdl_user_input_handler_parse_frequency("au", &value));
    ASSERT_STDL_KO(stdl_user_input_handler_parse_frequency("1.25aux", &value));
}

void test_user_input_make_context() {
    stdl_user_input_handler* inp = NULL;

    char* args[] =  {
            "self",
            "--data_output=test_make_context.h5",
            "--ctx_source=../tests/test_files/water_631g.molden",
            "--ctx_source_type=MOLDEN",
            "--ctx_gammaJ=0.5",
            "--ctx_gammaK=1.0",
            "--ctx_ethr=20eV",
            "--ctx_ax=1.0",
    };

    ASSERT_STDL_OK(stdl_user_input_handler_new_from_args(sizeof(args) / sizeof(char *), args, &inp));

    TEST_ASSERT_EQUAL_STRING("../tests/test_files/water_631g.molden", inp->ctx_source);
    TEST_ASSERT_EQUAL(STDL_SRC_MOLDEN, inp->ctx_source_type);
    TEST_ASSERT_EQUAL_STRING("test_make_context.h5", inp->data_output);
    TEST_ASSERT_FLOAT_WITHIN(1e-4, 0.5, inp->ctx_gammaJ);
    TEST_ASSERT_FLOAT_WITHIN(1e-4, 1.0, inp->ctx_gammaK);
    TEST_ASSERT_FLOAT_WITHIN(1e-4, 20.f / STDL_CONST_AU_TO_EV, inp->ctx_ethr);
    TEST_ASSERT_FLOAT_WITHIN(1e-4, 1.0, inp->ctx_ax);
    TEST_ASSERT_EQUAL(1, inp->ctx_tda);

    // create context
    stdl_context* ctx = NULL;
    ASSERT_STDL_OK(stdl_user_input_handler_make_context(inp, &ctx));
    TEST_ASSERT_NOT_NULL(ctx);

    TEST_ASSERT_FLOAT_WITHIN(1e-4, 0.5, ctx->gammaJ);
    TEST_ASSERT_FLOAT_WITHIN(1e-4, 1.0, ctx->gammaK);
    TEST_ASSERT_FLOAT_WITHIN(1e-4, 20.f / STDL_CONST_AU_TO_EV, ctx->ethr);
    TEST_ASSERT_FLOAT_WITHIN(1e-4, inp->ctx_e2thr, ctx->e2thr);
    TEST_ASSERT_FLOAT_WITHIN(1e-4, 1.0, ctx->ax);
    TEST_ASSERT_NULL(ctx->B); // because ctx_tda=0 by default

    // just check context
    hid_t file_id = H5Fopen(inp->data_output, H5F_ACC_RDONLY, H5P_DEFAULT);
    TEST_ASSERT_NOT_EQUAL(file_id, H5I_INVALID_HID);
    stdl_context* ctx2 = NULL;
    ASSERT_STDL_OK(stdl_context_load_h5(file_id, &ctx2));
    TEST_ASSERT_EQUAL(ctx2->nmo, ctx->nmo);
    ASSERT_STDL_OK(stdl_context_delete(ctx2));
    H5Fclose(file_id);

    ASSERT_STDL_OK(stdl_context_delete(ctx));
    ASSERT_STDL_OK(stdl_user_input_handler_delete(inp));
}
