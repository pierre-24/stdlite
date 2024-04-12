#include <unistd.h>

#include <stdlite/helpers.h>

#include "user_input_handler.h"
#include "responses_handler.h"
#include "tests_suite.h"

FILE* stream;

void setUp(void) {
    stdl_set_debug_level(0);
    stdl_set_log_level(0);

    stream = tmpfile();
}

void tearDown(void) {
    fclose(stream);
}

void test_user_input_compute_responses() {
    stdl_user_input_handler *inp = NULL;
    stdl_user_input_handler_new(&inp);
    TEST_ASSERT_NOT_NULL(inp);

    fputs("title = \"test calculation\"\n"
          "data_output=\"test_compute_response.h5\"\n"
          "[context]\n"
          "source = \"../tests/test_files/chiral_sto3g.molden\"\n"
          "source_type = \"MOLDEN\"\n"
          "ethr = '40eV'\n"
          "gammaJ = 2.0\n"
          "gammaK = 4.0\n"
          "ax = 1.0\n"
          "tda=0\n"
          "[responses]\n"
          "linear = [{opA = 'dipl', opB = 'dipl', wB = '532nm'}, {opA = 'angm', opB = 'dipl', wB = '532nm'}]\n"
          "#quadratic = [{opA = 'dipl', opB = 'dipl', opC = 'dipl', wB = '1064nm', wC = '1064nm'}]\n"
          "linear_sr = [{opA = 'dipl', nroots = -1}, {opA = 'angm', nroots = -1}]",
          stream);
    rewind(stream);

    ASSERT_STDL_OK(stdl_user_input_handler_fill_from_toml(inp, stream));
    ASSERT_STDL_OK(stdl_user_input_handler_check(inp));

    // create context
    stdl_context *ctx = NULL;
    ASSERT_STDL_OK(stdl_user_input_handler_make_context(inp, &ctx));
    TEST_ASSERT_NOT_NULL(ctx);

    // prepare responses
    stdl_responses_handler *rh = NULL;
    ASSERT_STDL_OK(stdl_responses_handler_new_from_input(inp, ctx, &rh));

    stdl_op_data* op_data_dipl = rh->lrvs_data[STDL_OP_DIPL];
    TEST_ASSERT_NOT_NULL(op_data_dipl);

    TEST_ASSERT_EQUAL(STDL_OP_DIPL, op_data_dipl->op);
    TEST_ASSERT_EQUAL(1, op_data_dipl->nlrvs);
    TEST_ASSERT_FLOAT_ARRAY_WITHIN(1e-6, ((float []) {STDL_CONST_HC / 532.f}), op_data_dipl->w, 1);

    stdl_op_data* op_data_angm = rh->lrvs_data[STDL_OP_ANGM];
    TEST_ASSERT_NOT_NULL(op_data_angm);

    TEST_ASSERT_EQUAL(STDL_OP_ANGM, op_data_angm->op);
    TEST_ASSERT_EQUAL(1, op_data_angm->nlrvs);
    TEST_ASSERT_FLOAT_ARRAY_WITHIN(1e-6, ((float []) {STDL_CONST_HC / 532.f}), op_data_angm->w, 1);

    // compute responses
    ASSERT_STDL_OK(stdl_responses_handler_compute(rh, inp, ctx));

    TEST_ASSERT_EQUAL(STDL_OP_DIPL, op_data_dipl->lrvs[0].op);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, STDL_CONST_HC / 532.f, op_data_dipl->lrvs[0].w);

    TEST_ASSERT_EQUAL(STDL_OP_ANGM, op_data_angm->lrvs[0].op);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, STDL_CONST_HC / 532.f, op_data_angm->lrvs[0].w);

    // compute properties
    ASSERT_STDL_OK(stdl_response_handler_compute_properties(rh, inp, ctx));

    // delete op_data_dipl output
    unlink(inp->data_output);

    ASSERT_STDL_OK(stdl_responses_handler_delete(rh));
    ASSERT_STDL_OK(stdl_context_delete(ctx));
    ASSERT_STDL_OK(stdl_user_input_handler_delete(inp));
}
