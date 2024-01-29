#include <string.h>
#include <unistd.h>

#include <stdlite/utils/fchk_parser.h>

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
    STDL_OK(stdl_lexer_new(&lx, f));
    TEST_ASSERT_NOT_NULL(lx);

    STDL_OK(stdl_fchk_parser_skip_intro(lx));

    stdl_wavefunction * wf = NULL;
    stdl_basis* bs = NULL;
    STDL_OK(stdl_fchk_parser_extract(&wf, &bs, lx));
    STDL_OK(stdl_wavefunction_delete(wf));

    _check_basis(bs);

    STDL_OK(stdl_basis_delete(bs));
    STDL_OK(stdl_lexer_delete(lx));

    fclose(f);
}
