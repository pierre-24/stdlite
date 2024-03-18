#include <stdlite/helpers.h>

#include "user_input_handler.h"
#include "responses_handler.h"
#include "tests_suite.h"

FILE* stream;

void setUp(void) {
    stdl_set_debug_level(-1);

    stream = tmpfile();
}

void tearDown(void) {
    fclose(stream);
}


void test_user_input_prepare_responses() {
    stdl_user_input_handler* inp = NULL;
    stdl_user_input_handler_new(&inp);
    TEST_ASSERT_NOT_NULL(inp);

    fputs("title = \"test calculation\"\n"
          "[context]\n"
          "source = \"../tests/test_files/water_631g.fchk\"\n"
          "source_type = \"FCHK\"\n"
          "ethr = '12eV'\n"
          "[responses]\n"
          "linear = [{opA = 'dipl', opB = 'dipl', wB = '1064nm'}, {opA = 'dipl', opB = 'dipl', wB = '532nm'}]\n"
          "quadratic = [{opA = 'dipl', opB = 'dipl', opC = 'dipl', wB = '1064nm', wC = '1064nm'}]\n"
          "linear_sr = [{opA = 'dipl', nroots = -1}]",
          stream);
    rewind(stream);

    ASSERT_STDL_OK(stdl_user_input_handler_fill_from_toml(inp, stream));
    ASSERT_STDL_OK(stdl_user_input_handler_check(inp));

    // create context
    stdl_context* ctx = NULL;
    ASSERT_STDL_OK(stdl_user_input_handler_make_context(inp, &ctx));
    TEST_ASSERT_NOT_NULL(ctx);

    // prepare responses
    stdl_responses_handler* rh = NULL;
    ASSERT_STDL_OK(stdl_user_input_handler_prepare_responses(inp, ctx, &rh));

    TEST_ASSERT_EQUAL(1, rh->nops);
    TEST_ASSERT_EQUAL(STDL_OP_DIPL, rh->ops[0]);
    TEST_ASSERT_EQUAL(1, rh->nlrvreqs);
    TEST_ASSERT_EQUAL(ctx->ncsfs, rh->nexci);

    stdl_lrv_request * lrv_req = rh->lrvreqs[0];
    TEST_ASSERT_NOT_NULL(lrv_req);
    TEST_ASSERT_EQUAL(STDL_OP_DIPL, lrv_req->op);
    TEST_ASSERT_EQUAL(2, lrv_req->nw);
    TEST_ASSERT_FLOAT_WITHIN(1e-4, 1064.f, STDL_CONST_HC / lrv_req->w[0]);
    TEST_ASSERT_FLOAT_WITHIN(1e-4, 532.f, STDL_CONST_HC / lrv_req->w[1]);

    stdl_response_request* req = inp->res_resreqs;
    TEST_ASSERT_NOT_NULL(req);

    TEST_ASSERT_EQUAL(1, req->resp_order);
    TEST_ASSERT_EQUAL(0, req->res_order);
    TEST_ASSERT_FLOAT_WITHIN(1e-4, STDL_CONST_HC / req->w[0], 1064.f);
    TEST_ASSERT_EQUAL(req->ops[0], STDL_OP_DIPL);
    TEST_ASSERT_EQUAL(0, req->nroots);
    TEST_ASSERT_EQUAL(lrv_req, req->lrvreqs[0]);
    TEST_ASSERT_EQUAL(0, req->wpos[0]);

    TEST_ASSERT_NOT_NULL(req->next);
    req = req->next;
    TEST_ASSERT_EQUAL(1, req->resp_order);
    TEST_ASSERT_EQUAL(0, req->res_order);
    TEST_ASSERT_FLOAT_WITHIN(1e-4, STDL_CONST_HC / req->w[0], 532.f);
    TEST_ASSERT_EQUAL(req->ops[0], STDL_OP_DIPL);
    TEST_ASSERT_EQUAL(0, req->nroots);
    TEST_ASSERT_EQUAL(lrv_req, req->lrvreqs[0]);
    TEST_ASSERT_EQUAL(1, req->wpos[0]);

    TEST_ASSERT_NOT_NULL(req->next);
    req = req->next;
    TEST_ASSERT_EQUAL(2, req->resp_order);
    TEST_ASSERT_EQUAL(0, req->res_order);
    TEST_ASSERT_FLOAT_WITHIN(1e-4, STDL_CONST_HC / req->w[0], 532.f);
    TEST_ASSERT_FLOAT_WITHIN(1e-4, STDL_CONST_HC / req->w[1], 1064.f);
    TEST_ASSERT_FLOAT_WITHIN(1e-4, STDL_CONST_HC / req->w[2], 1064.f);
    TEST_ASSERT_EQUAL(req->ops[0], STDL_OP_DIPL);
    TEST_ASSERT_EQUAL(0, req->nroots);
    TEST_ASSERT_EQUAL(lrv_req, req->lrvreqs[0]);
    TEST_ASSERT_EQUAL(1, req->wpos[0]);
    TEST_ASSERT_EQUAL(0, req->wpos[1]);
    TEST_ASSERT_EQUAL(0, req->wpos[2]);

    TEST_ASSERT_NOT_NULL(req->next);
    req = req->next;
    TEST_ASSERT_EQUAL(1, req->resp_order);
    TEST_ASSERT_EQUAL(1, req->res_order);
    TEST_ASSERT_EQUAL(req->ops[0], STDL_OP_DIPL);
    TEST_ASSERT_EQUAL(-1, req->nroots);
    TEST_ASSERT_NULL(req->w);
    TEST_ASSERT_NULL(req->lrvreqs);
    TEST_ASSERT_NULL(req->wpos);

    ASSERT_STDL_OK(stdl_responses_handler_delete(rh));
    ASSERT_STDL_OK(stdl_context_delete(ctx));
    ASSERT_STDL_OK(stdl_user_input_handler_delete(inp));
}

void test_user_input_compute_responses() {
    stdl_user_input_handler *inp = NULL;
    stdl_user_input_handler_new(&inp);
    TEST_ASSERT_NOT_NULL(inp);

    fputs("title = \"test calculation\"\n"
          "[context]\n"
          "source = \"../tests/test_files/water_631g.fchk\"\n"
          "source_type = \"FCHK\"\n"
          "ethr = '12eV'\n"
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
    ASSERT_STDL_OK(stdl_user_input_handler_prepare_responses(inp, ctx, &rh));

    // compute responses
    ASSERT_STDL_OK(stdl_responses_handler_compute(rh, ctx));

    // compute properties
    ASSERT_STDL_OK(stdl_user_input_handler_compute_properties(inp, ctx, rh));

    ASSERT_STDL_OK(stdl_responses_handler_delete(rh));
    ASSERT_STDL_OK(stdl_context_delete(ctx));
    ASSERT_STDL_OK(stdl_user_input_handler_delete(inp));
}
