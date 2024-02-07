#include <string.h>
#include <unistd.h>

#include <stdlite/utils/fchk_parser.h>
#include <stdlite/utils/matrix.h>

#include "tests_suite.h"


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
    char* fchk_path = "../tests/test_files/water_sto3g.fchk";

    FILE* f = fopen(fchk_path, "r");
    TEST_ASSERT_NOT_NULL(f);

    stdl_lexer* lx = NULL;
    STDL_OK(stdl_lexer_new(f, &lx));
    TEST_ASSERT_NOT_NULL(lx);

    STDL_OK(stdl_fchk_parser_skip_intro(lx));

    stdl_wavefunction * wf = NULL;
    stdl_basis* bs = NULL;
    STDL_OK(stdl_fchk_parser_extract(lx, &wf, &bs));
    STDL_OK(stdl_wavefunction_delete(wf));

    _check_basis(bs);

    STDL_OK(stdl_basis_delete(bs));
    STDL_OK(stdl_lexer_delete(lx));

    fclose(f);
}

void test_ovlp_ok() {
    char* fchk_path = "../tests/test_files/water_631g.fchk";

    FILE* f = fopen(fchk_path, "r");
    TEST_ASSERT_NOT_NULL(f);

    stdl_lexer* lx = NULL;
    STDL_OK(stdl_lexer_new(f, &lx));
    TEST_ASSERT_NOT_NULL(lx);

    STDL_OK(stdl_fchk_parser_skip_intro(lx));

    stdl_wavefunction * wf = NULL;
    stdl_basis* bs = NULL;
    STDL_OK(stdl_fchk_parser_extract(lx, &wf, &bs));

    double* S = malloc(STDL_MATRIX_SP_SIZE(wf->nao) * sizeof(double));
    TEST_ASSERT_NOT_NULL(S);

    stdl_basis_compute_dsp_ovlp(bs, S);

    // check that <i|i> = 1
    for (size_t i = 0; i < wf->nao; ++i)
        TEST_ASSERT_DOUBLE_WITHIN(1e-8, 1.0, S[STDL_MATRIX_SP_IDX(i, i)]);

    free(S);

    STDL_OK(stdl_wavefunction_delete(wf));
    STDL_OK(stdl_basis_delete(bs));
    STDL_OK(stdl_lexer_delete(lx));

    fclose(f);
}

void test_dipoles_ok() {
    char* fchk_path = "../tests/test_files/water_sto3g.fchk";

    FILE* f = fopen(fchk_path, "r");
    TEST_ASSERT_NOT_NULL(f);

    stdl_lexer* lx = NULL;
    STDL_OK(stdl_lexer_new(f, &lx));
    TEST_ASSERT_NOT_NULL(lx);

    STDL_OK(stdl_fchk_parser_skip_intro(lx));

    stdl_wavefunction * wf = NULL;
    stdl_basis* bs = NULL;
    STDL_OK(stdl_fchk_parser_extract(lx, &wf, &bs));

    float* dipoles = malloc(3 * STDL_MATRIX_SP_SIZE(wf->nao) * sizeof(float ));
    TEST_ASSERT_NOT_NULL(dipoles);

    stdl_basis_compute_ssp_dipole(bs, dipoles);

    stdl_matrix_ssp_print(wf->nao, dipoles + 0 * STDL_MATRIX_SP_SIZE(wf->nao), "xlint");
    stdl_matrix_ssp_print(wf->nao, dipoles + 1 * STDL_MATRIX_SP_SIZE(wf->nao), "ylint");
    stdl_matrix_ssp_print(wf->nao, dipoles + 2 * STDL_MATRIX_SP_SIZE(wf->nao), "zlint");

    free(dipoles);

    STDL_OK(stdl_wavefunction_delete(wf));
    STDL_OK(stdl_basis_delete(bs));
    STDL_OK(stdl_lexer_delete(lx));

    fclose(f);
}
