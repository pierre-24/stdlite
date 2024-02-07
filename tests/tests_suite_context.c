#include <string.h>

#include <stdlite/utils/fchk_parser.h>
#include <stdlite/utils/matrix.h>
#include <stdlite/context.h>

#include "tests_suite.h"

void setUp() {
    stdl_set_debug_level(-1);
}

void test_context_select_MO_ok() {
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
    STDL_OK(stdl_context_new(wf, bs, 2.0, 4.0, 10. / 27.212, 1e-4, 1.0, &ctx));

    TEST_ASSERT_EQUAL_INT(5, ctx->nmo);
    TEST_ASSERT_EQUAL_INT(3, ctx->nocc);

    // check that the MO are normalized
    for (size_t i = 0; i < ctx->nmo; ++i) {
        double sum = .0;

        for (size_t j = 0; j < wf->nao; ++j) {
            sum += pow(ctx->C[i * wf->nao + j], 2.0);
        }

        TEST_ASSERT_DOUBLE_WITHIN(1e-8, 1., sum);
    }

    // compute the density matrix
    double* P = malloc(wf->nao * wf->nao * sizeof(double));
    TEST_ASSERT_NOT_NULL(P);

    STDL_OK(stdl_wavefunction_compute_density_dsy(ctx->nocc, ctx->nmo, wf->nao, ctx->C, P));

    // check that density is symmetric
    for (size_t i = 0; i < wf->nao; ++i) {
        for (size_t j = 0; j <=i ; ++j) {
            TEST_ASSERT_EQUAL_DOUBLE(P[i * wf->nao + j], P[j * wf->nao + i]);
        }
    }

    // count electrons
    double total = .0;
    for(size_t i=0; i < wf->nao; i++)
        total += P[i * wf->nao + i];

    TEST_ASSERT_DOUBLE_WITHIN(1e-6, 2.f * ctx->nocc, total);

    // stdl_matrix_dge_print(wf->nao, 0, P, "D");

    free(P);

    STDL_OK(stdl_context_delete(ctx));
    STDL_OK(stdl_lexer_delete(lx));

    fclose(f);
}


void test_context_select_csfs_ok() {
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

    TEST_ASSERT_EQUAL_INT(0, ctx->ncsfs);

    STDL_OK(stdl_context_select_csfs_monopole(ctx, 0));

    TEST_ASSERT_EQUAL_INT(10, ctx->ncsfs);
    TEST_ASSERT_NULL(ctx->B);

    // check that energies are in increasing order
    for (size_t kia = 0; kia < ctx->ncsfs; ++kia) {
        TEST_ASSERT_EQUAL_FLOAT(ctx->A[STDL_MATRIX_SP_IDX(kia, kia)], ctx->ecsfs[kia]);

        if(kia > 0)
            TEST_ASSERT_TRUE(ctx->ecsfs[kia-1] <= ctx->ecsfs[kia]);
    }

    // stdl_matrix_ssp_print(ctx->ncsfs, ctx->A, "A");

    STDL_OK(stdl_context_delete(ctx));
    STDL_OK(stdl_lexer_delete(lx));

    fclose(f);
}
