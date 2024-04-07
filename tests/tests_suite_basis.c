#include <string.h>

#include <stdlite/utils/fchk_parser.h>
#include <stdlite/utils/matrix.h>

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

    stdl_basis_dsp_ovlp(bs, S);

    stdl_matrix_dsp_print(wf->nao, S, "S");

    // check that <i|i> = 1 (not so obvious in libcint with cartesian functions)
    for (size_t i = 0; i < wf->nao; ++i)
        TEST_ASSERT_DOUBLE_WITHIN(1e-7, 1.0, S[STDL_MATRIX_SP_IDX(i, i)]);

    free(S);
}

void test_dipoles_ok() {

    double* dipoles = malloc(3 * STDL_MATRIX_SP_SIZE(wf->nao) * sizeof(double));
    TEST_ASSERT_NOT_NULL(dipoles);

    stdl_basis_dsp_diplen(bs, dipoles);

    stdl_matrix_dsp_print(wf->nao, dipoles + 0 * STDL_MATRIX_SP_SIZE(wf->nao), "xdipl");
    stdl_matrix_dsp_print(wf->nao, dipoles + 1 * STDL_MATRIX_SP_SIZE(wf->nao), "ydipl");
    stdl_matrix_dsp_print(wf->nao, dipoles + 2 * STDL_MATRIX_SP_SIZE(wf->nao), "zdipl");

    free(dipoles);
}

void test_angmoms_ok() {

    double* angmoms = malloc(3 * STDL_MATRIX_SP_SIZE(wf->nao) * sizeof(double));
    TEST_ASSERT_NOT_NULL(angmoms);

    stdl_basis_dsp_angmom(bs, angmoms);


    stdl_matrix_dsp_print(wf->nao, angmoms + 0 * STDL_MATRIX_SP_SIZE(wf->nao), "xmagom");
    stdl_matrix_dsp_print(wf->nao, angmoms + 1 * STDL_MATRIX_SP_SIZE(wf->nao), "ymagmom");
    stdl_matrix_dsp_print(wf->nao, angmoms + 2 * STDL_MATRIX_SP_SIZE(wf->nao), "zmagmom");

    free(angmoms);
}

void test_dipvels_ok() {
    double* dipvels = malloc(3 * STDL_MATRIX_SP_SIZE(wf->nao) * sizeof(double));
    TEST_ASSERT_NOT_NULL(dipvels);

    stdl_basis_dsp_dipvel(bs, dipvels);

    stdl_matrix_dsp_print(wf->nao, dipvels + 0 * STDL_MATRIX_SP_SIZE(wf->nao), "xdipvel");
    stdl_matrix_dsp_print(wf->nao, dipvels + 1 * STDL_MATRIX_SP_SIZE(wf->nao), "ydipvel");
    stdl_matrix_dsp_print(wf->nao, dipvels + 2 * STDL_MATRIX_SP_SIZE(wf->nao), "zdipvel");

    free(dipvels);
}
