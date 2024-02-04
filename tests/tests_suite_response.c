#include <stdlite/response.h>
#include <stdlite/utils/fchk_parser.h>

#include "tests_suite.h"

void test_response_TDA_full_ok() {
    char* fchk_path = "../tests/test_files/water_631g.fchk";

    FILE* f = fopen(fchk_path, "r");
    TEST_ASSERT_NOT_NULL(f);

    stdl_lexer* lx = NULL;
    STDL_OK(stdl_lexer_new(f, &lx));
    STDL_OK(stdl_fchk_parser_skip_intro(lx));

    stdl_wavefunction * wf = NULL;
    stdl_basis * bs = NULL;

    STDL_OK(stdl_fchk_parser_extract(lx, &wf, &bs));

    stdl_context* ctx = NULL;
    STDL_OK(stdl_context_new(wf, bs, 2.0, 4.0, 12. / 27.212, 1e-4, 1.0, &ctx));

    TEST_ASSERT_EQUAL_INT(ctx->nmo,7);
    TEST_ASSERT_EQUAL_INT(ctx->nocc, 4);

    size_t nselected = 0;
    size_t* csfs = NULL;
    float * A = NULL;
    STDL_OK(stdl_context_select_csfs_monopole(ctx, &nselected, &csfs, &A, NULL));

    float* energies = NULL;
    float* amplitudes = NULL;
    size_t nexci = 0;

    stdl_response_casida_TDA_full(ctx, nselected, A, &energies, &amplitudes);

    // in this case, the eigenvalues are more or less the diagonal elements of A.
    for (size_t kia = 0; kia < nselected; ++kia) {
        TEST_ASSERT_FLOAT_WITHIN(1e-2, A[kia * nselected + kia], energies[kia]);
    }

    free(csfs);
    free(A);
    free(energies);
    free(amplitudes);

    STDL_OK(stdl_context_delete(ctx));
    STDL_OK(stdl_lexer_delete(lx));

    fclose(f);
}
