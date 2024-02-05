#include <stdlite/response.h>
#include <stdlite/utils/fchk_parser.h>
#include <stdlite/utils/matrix.h>
#include <stdlite/helpers.h>
#include <string.h>

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

    // save diagonal elements (i.e., the energies of the CSFs)
    float* A_diag = malloc(nselected * sizeof(float ));
    TEST_ASSERT_NOT_NULL(A_diag);

    for (size_t kia = 0; kia < nselected; ++kia)
        A_diag[kia] = A[STDL_MATRIX_SP_IDX(kia, kia)];

    float* energies = NULL;
    float* amplitudes = NULL;

    STDL_OK(stdl_response_casida_TDA_full(ctx, nselected, A, &energies, &amplitudes));

    for (size_t kia = 0; kia < nselected; ++kia) {
        // in this case, the eigenvalues are more or less the diagonal elements of A.
        TEST_ASSERT_FLOAT_WITHIN(1e-2, A_diag[kia], energies[kia]);

        // check that it is normed
        float sum = .0f;
        for (size_t kjb = 0; kjb < nselected; kjb++) {
            sum += powf(amplitudes[kia * nselected + kjb], 2);
        }

        TEST_ASSERT_FLOAT_WITHIN(1e-6, 1.0f, sum);

    }

    STDL_FREE_ALL(csfs, A, energies, amplitudes, A_diag);

    STDL_OK(stdl_context_delete(ctx));
    STDL_OK(stdl_lexer_delete(lx));

    fclose(f);
}

void test_response_TDA_ok() {
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

    size_t nselected = 0;
    size_t* csfs = NULL;
    float * A = NULL;
    STDL_OK(stdl_context_select_csfs_monopole(ctx, &nselected, &csfs, &A, NULL));

    float* Ap = malloc(STDL_MATRIX_SP_SIZE(nselected) * sizeof(float ));
    TEST_ASSERT_NOT_NULL(Ap);
    memcpy(Ap, A, STDL_MATRIX_SP_SIZE(nselected) * sizeof(float ));

    float* energies = NULL;
    float* amplitudes = NULL;

    STDL_OK(stdl_response_casida_TDA_full(ctx, nselected, A, &energies, &amplitudes));

    size_t nrequested = 5;
    float* first_energies = NULL;
    float* first_amplitudes = NULL;

    STDL_OK(stdl_response_casida_TDA(ctx, nselected, Ap, nrequested, &first_energies, &first_amplitudes));

    for (size_t kia = 0; kia < nrequested; ++kia) {
        // the same eigenvalues should have been obtained
        TEST_ASSERT_FLOAT_WITHIN(1e-5, energies[kia], first_energies[kia]);
    }

    // stdl_matrix_sge_print(nrequested, nselected, first_amplitudes, "CSFS");

    STDL_FREE_ALL(csfs, A, Ap, energies, amplitudes, first_energies, first_amplitudes);

    STDL_OK(stdl_context_delete(ctx));
    STDL_OK(stdl_lexer_delete(lx));

    fclose(f);
}
