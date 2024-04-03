#include <string.h>

#include <stdlite/utils/fchk_parser.h>
#include <stdlite/utils/matrix.h>

#include "tests_suite.h"

void setUp(void) {
    stdl_set_debug_level(0);
    stdl_set_log_level(2);
}


void _check_basis(stdl_basis* bs) {

    // check that <i|i> = 1.
    double* buf = malloc(3 * 3 * sizeof(double));

    for(int i=0; i < bs->nbas; i++) {
        int1e_ovlp_cart(buf, NULL, (int[]) {i, i}, bs->atm, bs->natm, bs->bas, bs->nbas, bs->env, NULL, NULL);
        int n = CINTcgtos_cart(i, bs->bas);
        for(int j=0; j < n; j++) {
            TEST_ASSERT_DOUBLE_WITHIN(1e-8, 1.f, buf[j * n + j]);
        }
    }

    free(buf);
}

void test_basis_functions_ok() {
    stdl_wavefunction * wf = NULL;
    stdl_basis* bs = NULL;
    read_fchk("../tests/test_files/water_sto3g.fchk", &wf, &bs);
    ASSERT_STDL_OK(stdl_wavefunction_delete(wf));

    _check_basis(bs);

    ASSERT_STDL_OK(stdl_basis_delete(bs));
}

void test_ovlp_ok() {
    stdl_wavefunction * wf = NULL;
    stdl_basis* bs = NULL;
    read_fchk("../tests/test_files/water_631g.fchk", &wf, &bs);

    double* S = malloc(STDL_MATRIX_SP_SIZE(wf->nao) * sizeof(double));
    TEST_ASSERT_NOT_NULL(S);

    stdl_basis_dsp_ovlp(bs, S);

    // check that <i|i> = 1
    for (size_t i = 0; i < wf->nao; ++i)
        TEST_ASSERT_DOUBLE_WITHIN(1e-8, 1.0, S[STDL_MATRIX_SP_IDX(i, i)]);

    free(S);

    ASSERT_STDL_OK(stdl_wavefunction_delete(wf));
    ASSERT_STDL_OK(stdl_basis_delete(bs));
}

void test_dipoles_ok() {
    stdl_wavefunction * wf = NULL;
    stdl_basis* bs = NULL;
    read_molden("../tests/test_files/He.molden", &wf, &bs);

    double* dipoles = malloc(3 * STDL_MATRIX_SP_SIZE(wf->nao) * sizeof(double));
    TEST_ASSERT_NOT_NULL(dipoles);

    stdl_basis_dsp_dipole(bs, dipoles);

    stdl_matrix_dsp_print(wf->nao, dipoles + 0 * STDL_MATRIX_SP_SIZE(wf->nao), "xdipl");
    stdl_matrix_dsp_print(wf->nao, dipoles + 1 * STDL_MATRIX_SP_SIZE(wf->nao), "ydipl");
    stdl_matrix_dsp_print(wf->nao, dipoles + 2 * STDL_MATRIX_SP_SIZE(wf->nao), "zdipl");

    free(dipoles);

    ASSERT_STDL_OK(stdl_wavefunction_delete(wf));
    ASSERT_STDL_OK(stdl_basis_delete(bs));
}

void test_angmoms_ok() {
    stdl_wavefunction * wf = NULL;
    stdl_basis* bs = NULL;
    read_molden("../tests/test_files/He.molden", &wf, &bs);

    double* angmoms = malloc(3 * STDL_MATRIX_SP_SIZE(wf->nao) * sizeof(double));
    TEST_ASSERT_NOT_NULL(angmoms);

    stdl_basis_dsp_angmom(bs, angmoms);


    stdl_matrix_dsp_print(wf->nao, angmoms + 0 * STDL_MATRIX_SP_SIZE(wf->nao), "xmagom");
    stdl_matrix_dsp_print(wf->nao, angmoms + 1 * STDL_MATRIX_SP_SIZE(wf->nao), "ymagmom");
    stdl_matrix_dsp_print(wf->nao, angmoms + 2 * STDL_MATRIX_SP_SIZE(wf->nao), "zmagmom");

    free(angmoms);

    ASSERT_STDL_OK(stdl_wavefunction_delete(wf));
    ASSERT_STDL_OK(stdl_basis_delete(bs));
}
