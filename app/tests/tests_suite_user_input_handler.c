#include <stdlite/helpers.h>
#include <unistd.h>

#include "tests_suite.h"
#include "user_input_handler.h"
#include <stdlite/utils/matrix.h>

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
          "gammaK=0.5 \n"
          "gauge_origin=[1.0, 1.0, 1.0]",
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
    TEST_ASSERT_EQUAL(inp->ctx_gauge_origin, STDL_GAUGE_ORIGIN_CUSTOM);
    TEST_ASSERT_DOUBLE_ARRAY_WITHIN(1e-6, inp->ctx_gauge_origin_custom, ((double []) {1. / STDL_CONST_AU_TO_ANG, 1. / STDL_CONST_AU_TO_ANG, 1. / STDL_CONST_AU_TO_ANG}), 3);

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
    TEST_ASSERT_EQUAL(inp->ctx_gauge_origin, STDL_GAUGE_ORIGIN_CM);

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
    TEST_ASSERT_NULL(ctx->AmB); // because ctx_tda=0 by default

    // just check in H5
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

void test_user_input_reuse_context() {
    stdl_user_input_handler* inp = NULL;

    char* args[] =  {
            "self",
            "--data_output=test_remake_context.h5",
            "--ctx_source=../app/tests/test_files/context_water_sto3g.h5",
            "--ctx_source_type=STDL_CTX",
    };

    ASSERT_STDL_OK(stdl_user_input_handler_new_from_args(sizeof(args) / sizeof(char *), args, &inp));

    TEST_ASSERT_EQUAL_STRING("../app/tests/test_files/context_water_sto3g.h5", inp->ctx_source);
    TEST_ASSERT_EQUAL(STDL_SRC_CTX, inp->ctx_source_type);
    TEST_ASSERT_EQUAL_STRING("test_remake_context.h5", inp->data_output);

    // create context
    stdl_context* ctx = NULL;
    ASSERT_STDL_OK(stdl_user_input_handler_make_context(inp, &ctx));
    TEST_ASSERT_NOT_NULL(ctx);

    TEST_ASSERT_FLOAT_WITHIN(1e-4, 4.0, ctx->gammaJ);
    TEST_ASSERT_FLOAT_WITHIN(1e-4, 2.0, ctx->gammaK);
    TEST_ASSERT_FLOAT_WITHIN(1e-4, 50.f / STDL_CONST_AU_TO_EV, ctx->ethr);
    TEST_ASSERT_FLOAT_WITHIN(1e-4, 0.5, ctx->ax);
    TEST_ASSERT_NOT_NULL(ctx->AmB);

    // just check
    hid_t file_id = H5Fopen(inp->ctx_source, H5F_ACC_RDONLY, H5P_DEFAULT);
    TEST_ASSERT_NOT_EQUAL(file_id, H5I_INVALID_HID);

    hid_t ctx_group_id = H5Gopen1(file_id, "context");

    float* AmB = malloc(STDL_MATRIX_SP_SIZE(ctx->ncsfs) * sizeof(float ));
    TEST_ASSERT_NOT_NULL(AmB);

    herr_t status = H5LTread_dataset(ctx_group_id, "A-B", H5T_NATIVE_FLOAT, AmB);
    TEST_ASSERT_TRUE(status == 0);

    TEST_ASSERT_FLOAT_ARRAY_WITHIN(1e-6, ctx->AmB, AmB, STDL_MATRIX_SP_SIZE(ctx->ncsfs));

    H5Gclose(ctx_group_id);
    H5Fclose(file_id);

    STDL_FREE_IF_USED(AmB);

    ASSERT_STDL_OK(stdl_context_delete(ctx));
    ASSERT_STDL_OK(stdl_user_input_handler_delete(inp));
}

void test_user_input_response_requests_ok() {
    stdl_user_input_handler *inp = NULL;
    stdl_user_input_handler_new(&inp);
    TEST_ASSERT_NOT_NULL(inp);

    fputs("title = \"test calculation\"\n"
          "data_output=\"test_prepare_response.h5\"\n"
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

    TEST_ASSERT_EQUAL(inp->res_nw, 2);
    TEST_ASSERT_FLOAT_ARRAY_WITHIN(1e-6, ((float[]) {STDL_CONST_HC / 1064.f, STDL_CONST_HC / 532.f}), inp->res_w, 2);

    stdl_response_request *req = inp->res_resreqs;
    TEST_ASSERT_NOT_NULL(req);

    TEST_ASSERT_EQUAL(1, req->resp_order);
    TEST_ASSERT_EQUAL(0, req->resi_order);
    TEST_ASSERT_EQUAL(2, req->nlrvs);
    TEST_ASSERT_EQUAL_UINT64_ARRAY(((size_t[]) {0, 0}), req->iw, 2);
    TEST_ASSERT_EQUAL(2, req->nops);
    TEST_ASSERT_EQUAL_INT_ARRAY(((int[]) {STDL_OP_DIPL, STDL_OP_DIPL}), req->ops, 2);
    TEST_ASSERT_EQUAL(0, req->nroots);

    TEST_ASSERT_NOT_NULL(req->next);
    req = req->next;
    TEST_ASSERT_EQUAL(1, req->resp_order);
    TEST_ASSERT_EQUAL(0, req->resi_order);
    TEST_ASSERT_EQUAL(2, req->nlrvs);
    TEST_ASSERT_EQUAL_UINT64_ARRAY(((size_t[]){1, 1}), req->iw, 2);
    TEST_ASSERT_EQUAL(2, req->nops);
    TEST_ASSERT_EQUAL_INT_ARRAY(((int[]) {STDL_OP_DIPL, STDL_OP_DIPL}), req->ops, 2);
    TEST_ASSERT_EQUAL(0, req->nroots);

    TEST_ASSERT_NOT_NULL(req->next);
    req = req->next;
    TEST_ASSERT_EQUAL(2, req->resp_order);
    TEST_ASSERT_EQUAL(0, req->resi_order);
    TEST_ASSERT_EQUAL(3, req->nlrvs);
    TEST_ASSERT_EQUAL_UINT64_ARRAY(((size_t[]){1, 0, 0}), req->iw, 3);
    TEST_ASSERT_EQUAL(3, req->nops);
    TEST_ASSERT_EQUAL_INT_ARRAY(((int[]) {STDL_OP_DIPL, STDL_OP_DIPL, STDL_OP_DIPL}), req->ops, 3);
    TEST_ASSERT_EQUAL(0, req->nroots);

    TEST_ASSERT_NOT_NULL(req->next);
    req = req->next;
    TEST_ASSERT_EQUAL(1, req->resp_order);
    TEST_ASSERT_EQUAL(1, req->resi_order);
    TEST_ASSERT_EQUAL(0, req->nlrvs);
    TEST_ASSERT_EQUAL(2, req->nops);
    TEST_ASSERT_EQUAL(req->ops[0], STDL_OP_DIPL);
    TEST_ASSERT_EQUAL(req->ops[1], STDL_OP_DIPL);
    TEST_ASSERT_EQUAL(-1, req->nroots);
    TEST_ASSERT_NULL(req->iw);

    // delete data output
    unlink(inp->data_output);

    ASSERT_STDL_OK(stdl_user_input_handler_delete(inp));
}
