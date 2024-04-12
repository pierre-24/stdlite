#include <string.h>

#include <stdlite/utils/fchk_parser.h>
#include <stdlite/utils/matrix.h>
#include <stdlite/integrals.h>
#include <stdlite/helpers.h>

#include "tests_suite.h"

stdl_wavefunction * wf = NULL;
stdl_basis* bs = NULL;

void setUp(void) {
    stdl_set_debug_level(0);
    stdl_set_log_level(2);

    read_molden("../tests/test_files/water_sto3g_dalton.molden", &wf, &bs);
}

void tearDown(void) {
    ASSERT_STDL_OK(stdl_wavefunction_delete(wf));
    ASSERT_STDL_OK(stdl_basis_delete(bs));
}

void test_ovlp_ok() {

    double* S = malloc(STDL_MATRIX_SP_SIZE(wf->nao) * sizeof(double));
    TEST_ASSERT_NOT_NULL(S);

    ASSERT_STDL_OK(stdl_operator_int1e_dsp(bs, STDL_OP_OVLP, 1., S));

    stdl_matrix_dsp_print(2, wf->nao, S, "S");

    // check that <i|i> = 1 (not so obvious in libcint with cartesian functions)
    for (size_t i = 0; i < wf->nao; ++i)
        TEST_ASSERT_DOUBLE_WITHIN(1e-7, 1.0, S[STDL_MATRIX_SP_IDX(i, i)]);

    STDL_FREE_ALL(S);
}

void test_diplens_ok() {

    double* diplens = malloc(3 * STDL_MATRIX_SP_SIZE(wf->nao) * sizeof(double));
    TEST_ASSERT_NOT_NULL(diplens);

    ASSERT_STDL_OK(stdl_operator_int1e_dsp(bs, STDL_OP_DIPL, -1., diplens));

    stdl_matrix_dsp_print(2, wf->nao, diplens + 0 * STDL_MATRIX_SP_SIZE(wf->nao), "xdipl");
    stdl_matrix_dsp_print(2, wf->nao, diplens + 1 * STDL_MATRIX_SP_SIZE(wf->nao), "ydipl");
    stdl_matrix_dsp_print(2, wf->nao, diplens + 2 * STDL_MATRIX_SP_SIZE(wf->nao), "zdipl");

    free(diplens);
}

void test_dipvels_ok() {
    double* dipvels = malloc(3 * STDL_MATRIX_SP_SIZE(wf->nao) * sizeof(double));
    TEST_ASSERT_NOT_NULL(dipvels);

    ASSERT_STDL_OK(stdl_operator_int1e_dsp(bs, STDL_OP_DIPV, 1., dipvels));

    stdl_matrix_dsp_print(2, wf->nao, dipvels + 0 * STDL_MATRIX_SP_SIZE(wf->nao), "xdipv");
    stdl_matrix_dsp_print(2, wf->nao, dipvels + 1 * STDL_MATRIX_SP_SIZE(wf->nao), "ydipv");
    stdl_matrix_dsp_print(2, wf->nao, dipvels + 2 * STDL_MATRIX_SP_SIZE(wf->nao), "zdipv");

    free(dipvels);
}

void test_angmoms_ok() {

    double* angmoms = malloc(3 * STDL_MATRIX_SP_SIZE(wf->nao) * sizeof(double));
    TEST_ASSERT_NOT_NULL(angmoms);

    ASSERT_STDL_OK(stdl_operator_int1e_dsp(bs, STDL_OP_ANGM, 1., angmoms));

    stdl_matrix_dsp_print(2, wf->nao, angmoms + 0 * STDL_MATRIX_SP_SIZE(wf->nao), "xangm");
    stdl_matrix_dsp_print(2, wf->nao, angmoms + 1 * STDL_MATRIX_SP_SIZE(wf->nao), "yangm");
    stdl_matrix_dsp_print(2, wf->nao, angmoms + 2 * STDL_MATRIX_SP_SIZE(wf->nao), "zangm");

    free(angmoms);
}
