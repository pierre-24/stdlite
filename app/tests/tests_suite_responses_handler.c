#include <unistd.h>

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
          "source = \"../tests/test_files/water_631g.fchk\"\n"
          "source_type = \"FCHK\"\n"
          "ethr = '12eV'\n"
          "gammaJ = 2.0\n"
          "gammaK = 4.0\n"
          "ax = 1.0\n"
          "[responses]\n"
          "linear = [{opA = 'dipl', opB = 'dipl', wB = '1064nm'}, {opA = 'dipl', opB = 'dipl', wB = '532nm'}]\n"
          "quadratic = [{opA = 'dipl', opB = 'dipl', opC = 'dipl', wB = '1064nm', wC = '1064nm'}]\n"
          "linear_sr = [{opA = 'dipl', nroots = -1}]",
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

    // compute responses
    ASSERT_STDL_OK(stdl_responses_handler_compute(rh, inp, ctx));

    /*
    // compute properties
    ASSERT_STDL_OK(stdl_user_input_handler_compute_properties(inp, ctx, rh));*/

    // delete data output
    unlink(inp->data_output);

    ASSERT_STDL_OK(stdl_responses_handler_delete(rh));
    ASSERT_STDL_OK(stdl_context_delete(ctx));
    ASSERT_STDL_OK(stdl_user_input_handler_delete(inp));
}
