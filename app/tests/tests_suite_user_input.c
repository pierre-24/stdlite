#include <stdlite/helpers.h>
#include <unistd.h>

#include "tests_suite.h"
#include "user_input.h"

FILE* stream;

void setUp(void) {
    stdl_set_debug_level(-1);

    stream = tmpfile();
}

void tearDown(void) {
    fclose(stream);
}


void test_user_input_context_fill_from_toml_ok() {
    stdl_user_input* inp = NULL;
    stdl_user_input_new(&inp);
    TEST_ASSERT_NOT_NULL(inp);

    fputs("title = \"test calculation\"\n"
          "[context]\n"
          "source = \"../tests/test_files/water_631g.fchk\"\n"
          "source_type = \"FCHK\"\n"
          "output = \"context_water_631g.h5\"\n"
          "method = \"monopole_direct\"\n"
          "tda = 0\n"
          "ethr = '12eV'\n"
          "e2thr=1e-3\n"
          "ax = 1.0\n"
          "gammaJ=1.0\n"
          "gammaK=0.5",
          stream);
    rewind(stream);

    ASSERT_STDL_OK(stdl_user_input_fill_from_toml(inp, stream));
    ASSERT_STDL_OK(stdl_user_input_check(inp));

    TEST_ASSERT_EQUAL_STRING("test calculation", inp->title);
    TEST_ASSERT_EQUAL_STRING("../tests/test_files/water_631g.fchk", inp->ctx_source);
    TEST_ASSERT_EQUAL(STDL_SRC_FCHK, inp->ctx_source_type);
    TEST_ASSERT_EQUAL_STRING("context_water_631g.h5", inp->ctx_output);
    TEST_ASSERT_EQUAL(STDL_METHOD_MONOPOLE_DIRECT, inp->ctx_method);
    TEST_ASSERT_EQUAL(0, inp->ctx_tda);
    TEST_ASSERT_FLOAT_WITHIN(1e-4, 12.f / STDL_CONST_AU_TO_EV, inp->ctx_ethr);
    TEST_ASSERT_FLOAT_WITHIN(1e-4, 1e-3, inp->ctx_e2thr);
    TEST_ASSERT_FLOAT_WITHIN(1e-4, 1.0, inp->ctx_ax);
    TEST_ASSERT_FLOAT_WITHIN(1e-4, 1.0, inp->ctx_gammaJ);
    TEST_ASSERT_FLOAT_WITHIN(1e-4, 0.5, inp->ctx_gammaK);

    ASSERT_STDL_OK(stdl_user_input_delete(inp));
}

void test_user_input_context_fill_from_args_ok() {
    stdl_user_input* inp = NULL;
    stdl_user_input_new(&inp);
    TEST_ASSERT_NOT_NULL(inp);

    char* args[] =  {
            "self",
            "--ctx_source=../tests/test_files/water_631g.molden",
            "--ctx_source_type=MOLDEN",
            "--ctx_output=test.h5",
            "--ctx_gammaJ=0.5",
            "--ctx_gammaK=1.0",
            "--ctx_ethr=12eV",
            "--ctx_e2thr=1e-3",
            "--ctx_ax=1.0",
            "--ctx_tda=0",
    };

    ASSERT_STDL_OK(stdl_user_input_fill_from_args(inp, sizeof(args) / sizeof(char*), args));
    ASSERT_STDL_OK(stdl_user_input_check(inp));

    TEST_ASSERT_EQUAL_STRING("../tests/test_files/water_631g.molden", inp->ctx_source);
    TEST_ASSERT_EQUAL(STDL_SRC_MOLDEN, inp->ctx_source_type);
    TEST_ASSERT_EQUAL_STRING("test.h5", inp->ctx_output);
    TEST_ASSERT_FLOAT_WITHIN(1e-4, 0.5, inp->ctx_gammaJ);
    TEST_ASSERT_FLOAT_WITHIN(1e-4, 1.0, inp->ctx_gammaK);
    TEST_ASSERT_FLOAT_WITHIN(1e-4, 12.f / STDL_CONST_AU_TO_EV, inp->ctx_ethr);
    TEST_ASSERT_FLOAT_WITHIN(1e-4, 1e-3, inp->ctx_e2thr);
    TEST_ASSERT_FLOAT_WITHIN(1e-4, 1.0, inp->ctx_ax);
    TEST_ASSERT_EQUAL(0, inp->ctx_tda);

    ASSERT_STDL_OK(stdl_user_input_delete(inp));
}

void test_user_input_context_fill_from_args_and_file_ok() {
    stdl_user_input* inp = NULL;
    stdl_user_input_new(&inp);
    TEST_ASSERT_NOT_NULL(inp);

    char* args[] =  {
            "self",
            "../app/tests/test_files/test.toml",
            "--ctx_e2thr=1e-3",
    };

    ASSERT_STDL_OK(stdl_user_input_fill_from_args(inp, sizeof(args) / sizeof(char*), args));
    ASSERT_STDL_OK(stdl_user_input_check(inp));

    TEST_ASSERT_EQUAL_STRING("test calculation", inp->title);
    TEST_ASSERT_FLOAT_WITHIN(1e-4, 1e-3, inp->ctx_e2thr);

    ASSERT_STDL_OK(stdl_user_input_delete(inp));
}

void test_user_input_incorrect_toml_ko() {
    stdl_user_input* inp = NULL;
    stdl_user_input_new(&inp);
    TEST_ASSERT_NOT_NULL(inp);

    fputs("[context]\n"
          "source = \"whatever\"\n"
          "source_type = \"whatever\"", // incorrect source type
          stream);
    rewind(stream);

    ASSERT_STDL_KO(stdl_user_input_fill_from_toml(inp, stream));
    ASSERT_STDL_OK(stdl_user_input_delete(inp));
}

void test_user_input_incorrect_check_args_ko() {
    stdl_user_input* inp = NULL;
    stdl_user_input_new(&inp);
    TEST_ASSERT_NOT_NULL(inp);

    fputs("[context]\n"
          "source = \"test.fchk\"\n"
          "ax = -1",  // ax should be in [0,1]
          stream);
    rewind(stream);

    ASSERT_STDL_OK(stdl_user_input_fill_from_toml(inp, stream));
    ASSERT_STDL_KO(stdl_user_input_check(inp));

    ASSERT_STDL_OK(stdl_user_input_delete(inp));
}


void test_parse_frequency_ok() {
    double value;

    ASSERT_STDL_OK(stdl_user_input_parse_frequency("1.25au", &value));
    TEST_ASSERT_DOUBLE_WITHIN(1e-6, 1.25, value);

    ASSERT_STDL_OK(stdl_user_input_parse_frequency("1.25", &value));
    TEST_ASSERT_DOUBLE_WITHIN(1e-6, 1.25, value);

    ASSERT_STDL_OK(stdl_user_input_parse_frequency("1.25eV", &value));
    TEST_ASSERT_DOUBLE_WITHIN(1e-6, 1.25 / STDL_CONST_AU_TO_EV, value);

    ASSERT_STDL_OK(stdl_user_input_parse_frequency("125nm", &value));
    TEST_ASSERT_DOUBLE_WITHIN(1e-6, STDL_CONST_HC / 125, value);

    ASSERT_STDL_KO(stdl_user_input_parse_frequency("1.25 au", &value));
    ASSERT_STDL_KO(stdl_user_input_parse_frequency("au", &value));
    ASSERT_STDL_KO(stdl_user_input_parse_frequency("1.25aux", &value));
}

void test_user_input_make_context() {
    stdl_user_input* inp = NULL;

    char* args[] =  {
            "self",
            "--ctx_source=../tests/test_files/water_631g.molden",
            "--ctx_source_type=MOLDEN",
            "--ctx_output=test.h5",
            "--ctx_gammaJ=0.5",
            "--ctx_gammaK=1.0",
            "--ctx_ethr=20eV",
            "--ctx_ax=1.0",
    };

    ASSERT_STDL_OK(stdl_user_input_new_from_args(sizeof(args) / sizeof(char*), args, &inp));

    TEST_ASSERT_EQUAL_STRING("../tests/test_files/water_631g.molden", inp->ctx_source);
    TEST_ASSERT_EQUAL(STDL_SRC_MOLDEN, inp->ctx_source_type);
    TEST_ASSERT_EQUAL_STRING("test.h5", inp->ctx_output);
    TEST_ASSERT_FLOAT_WITHIN(1e-4, 0.5, inp->ctx_gammaJ);
    TEST_ASSERT_FLOAT_WITHIN(1e-4, 1.0, inp->ctx_gammaK);
    TEST_ASSERT_FLOAT_WITHIN(1e-4, 20.f / STDL_CONST_AU_TO_EV, inp->ctx_ethr);
    TEST_ASSERT_FLOAT_WITHIN(1e-4, 1.0, inp->ctx_ax);
    TEST_ASSERT_EQUAL(1, inp->ctx_tda);

    // create context
    stdl_context* ctx = NULL;
    ASSERT_STDL_OK(stdl_user_input_make_context(inp, &ctx));
    TEST_ASSERT_NOT_NULL(ctx);

    TEST_ASSERT_FLOAT_WITHIN(1e-4, 0.5, ctx->gammaJ);
    TEST_ASSERT_FLOAT_WITHIN(1e-4, 1.0, ctx->gammaK);
    TEST_ASSERT_FLOAT_WITHIN(1e-4, 20.f / STDL_CONST_AU_TO_EV, ctx->ethr);
    TEST_ASSERT_FLOAT_WITHIN(1e-4, inp->ctx_e2thr, ctx->e2thr);
    TEST_ASSERT_FLOAT_WITHIN(1e-4, 1.0, ctx->ax);
    TEST_ASSERT_NULL(ctx->B); // because ctx_tda=0 by default

    // just check context
    hid_t file_id = H5Fopen(inp->ctx_output, H5F_ACC_RDONLY, H5P_DEFAULT);
    TEST_ASSERT_NOT_EQUAL(file_id, H5I_INVALID_HID);
    stdl_context* ctx2 = NULL;
    ASSERT_STDL_OK(stdl_context_load_h5(file_id, &ctx2));
    TEST_ASSERT_EQUAL(ctx2->nmo, ctx->nmo);
    ASSERT_STDL_OK(stdl_context_delete(ctx2));
    H5Fclose(file_id);

    ASSERT_STDL_OK(stdl_context_delete(ctx));
    ASSERT_STDL_OK(stdl_user_input_delete(inp));
}


void test_user_input_lresp_ok() {
    stdl_user_input* inp = NULL;
    stdl_user_input_new(&inp);
    TEST_ASSERT_NOT_NULL(inp);

    fputs("title = \"test calculation\"\n"
          "[context]\n"
          "source = \"../tests/test_files/water_631g.fchk\"\n"
          "source_type = \"FCHK\"\n"
          "[responses]\n"
          "linear = [{opA = 'dipl', opB = 'dipl', wB = '1064nm'}, {opA = 'dipl', opB = 'dipl', wB = '532nm'}]\n"
          "quadratic = [{opA = 'dipl', opB = 'dipl', opC = 'dipl', wB = '1064nm', wC = '1064nm'}]\n"
          "linear_sr = [{opA = 'dipl', root = -1}]",
          stream);
    rewind(stream);

    ASSERT_STDL_OK(stdl_user_input_fill_from_toml(inp, stream));
    ASSERT_STDL_OK(stdl_user_input_check(inp));

    stdl_response_request* req = inp->res_resreqs;
    TEST_ASSERT_NOT_NULL(req);

    TEST_ASSERT_EQUAL(1, req->resp_order);
    TEST_ASSERT_EQUAL(0, req->res_order);
    TEST_ASSERT_FLOAT_WITHIN(1e-4, STDL_CONST_HC / req->w[0], 1064.f);
    TEST_ASSERT_FLOAT_WITHIN(1e-4, STDL_CONST_HC / req->w[1], 1064.f);
    TEST_ASSERT_EQUAL(req->ops[0], STDL_OP_DIPL);
    TEST_ASSERT_EQUAL(0, req->nroot);

    TEST_ASSERT_NOT_NULL(req->next);
    req = req->next;
    TEST_ASSERT_EQUAL(1, req->resp_order);
    TEST_ASSERT_EQUAL(0, req->res_order);
    TEST_ASSERT_FLOAT_WITHIN(1e-4, STDL_CONST_HC / req->w[0], 532.f);
    TEST_ASSERT_FLOAT_WITHIN(1e-4, STDL_CONST_HC / req->w[1], 532.f);
    TEST_ASSERT_EQUAL(req->ops[0], STDL_OP_DIPL);
    TEST_ASSERT_EQUAL(0, req->nroot);

    TEST_ASSERT_NOT_NULL(req->next);
    req = req->next;
    TEST_ASSERT_EQUAL(2, req->resp_order);
    TEST_ASSERT_EQUAL(0, req->res_order);
    TEST_ASSERT_FLOAT_WITHIN(1e-4, STDL_CONST_HC / req->w[0], 532.f);
    TEST_ASSERT_FLOAT_WITHIN(1e-4, STDL_CONST_HC / req->w[1], 1064.f);
    TEST_ASSERT_FLOAT_WITHIN(1e-4, STDL_CONST_HC / req->w[2], 1064.f);
    TEST_ASSERT_EQUAL(req->ops[0], STDL_OP_DIPL);
    TEST_ASSERT_EQUAL(0, req->nroot);

    TEST_ASSERT_NOT_NULL(req->next);
    req = req->next;
    TEST_ASSERT_EQUAL(1, req->resp_order);
    TEST_ASSERT_EQUAL(1, req->res_order);
    TEST_ASSERT_EQUAL(req->ops[0], STDL_OP_DIPL);
    TEST_ASSERT_EQUAL(-1, req->nroot);

    ASSERT_STDL_OK(stdl_user_input_delete(inp));
}

void test_user_input_prepare_responses() {
    stdl_user_input* inp = NULL;
    stdl_user_input_new(&inp);
    TEST_ASSERT_NOT_NULL(inp);

    fputs("title = \"test calculation\"\n"
          "[context]\n"
          "source = \"../tests/test_files/water_631g.fchk\"\n"
          "source_type = \"FCHK\"\n"
          "ethr = '12eV'\n"
          "[responses]\n"
          "linear = [{opA = 'dipl', opB = 'dipl', wB = '1064nm'}, {opA = 'dipl', opB = 'dipl', wB = '532nm'}]\n"
          "quadratic = [{opA = 'dipl', opB = 'dipl', opC = 'dipl', wB = '1064nm', wC = '1064nm'}]\n",
          stream);
    rewind(stream);

    ASSERT_STDL_OK(stdl_user_input_fill_from_toml(inp, stream));
    ASSERT_STDL_OK(stdl_user_input_check(inp));

    // create context
    stdl_context* ctx = NULL;
    ASSERT_STDL_OK(stdl_user_input_make_context(inp, &ctx));
    TEST_ASSERT_NOT_NULL(ctx);

    // prepare responses
    ASSERT_STDL_OK(stdl_user_input_prepare_responses(inp, ctx));

    TEST_ASSERT_EQUAL(1, inp->res_nops);
    TEST_ASSERT_EQUAL(STDL_OP_DIPL, inp->res_ops[0]);
    TEST_ASSERT_EQUAL(1, inp->res_nlrvreq);
    TEST_ASSERT_EQUAL(0, inp->res_nexci);

    stdl_lrv_request * lrv_req = inp->res_lrvreqs[0];
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
    TEST_ASSERT_EQUAL(0, req->nroot);
    TEST_ASSERT_EQUAL(lrv_req, req->requests[0]);
    TEST_ASSERT_EQUAL(0, req->wpos[0]);

    TEST_ASSERT_NOT_NULL(req->next);
    req = req->next;
    TEST_ASSERT_EQUAL(1, req->resp_order);
    TEST_ASSERT_EQUAL(0, req->res_order);
    TEST_ASSERT_FLOAT_WITHIN(1e-4, STDL_CONST_HC / req->w[0], 532.f);
    TEST_ASSERT_EQUAL(req->ops[0], STDL_OP_DIPL);
    TEST_ASSERT_EQUAL(0, req->nroot);
    TEST_ASSERT_EQUAL(lrv_req, req->requests[0]);
    TEST_ASSERT_EQUAL(1, req->wpos[0]);

    TEST_ASSERT_NOT_NULL(req->next);
    req = req->next;
    TEST_ASSERT_EQUAL(2, req->resp_order);
    TEST_ASSERT_EQUAL(0, req->res_order);
    TEST_ASSERT_FLOAT_WITHIN(1e-4, STDL_CONST_HC / req->w[0], 532.f);
    TEST_ASSERT_FLOAT_WITHIN(1e-4, STDL_CONST_HC / req->w[1], 1064.f);
    TEST_ASSERT_FLOAT_WITHIN(1e-4, STDL_CONST_HC / req->w[2], 1064.f);
    TEST_ASSERT_EQUAL(req->ops[0], STDL_OP_DIPL);
    TEST_ASSERT_EQUAL(0, req->nroot);
    TEST_ASSERT_EQUAL(lrv_req, req->requests[0]);
    TEST_ASSERT_EQUAL(1, req->wpos[0]);
    TEST_ASSERT_EQUAL(0, req->wpos[1]);
    TEST_ASSERT_EQUAL(0, req->wpos[2]);

    ASSERT_STDL_OK(stdl_context_delete(ctx));
    ASSERT_STDL_OK(stdl_user_input_delete(inp));
}
