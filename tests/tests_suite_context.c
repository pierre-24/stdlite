#include <string.h>

#include <stdlite/utils/fchk_parser.h>
#include <stdlite/matrix.h>
#include <stdlite/context.h>

#include "tests_suite.h"

void setUp() {
    stdl_set_debug_level(3);
}


void test_context_select_MO_ok() {
    char* fchk_path = "../tests/test_files/water_631g.fchk";

    FILE* f = fopen(fchk_path, "r");
    TEST_ASSERT_NOT_NULL(f);

    stdl_lexer* lx = NULL;
    STDL_OK(stdl_lexer_new(&lx, f));
    STDL_OK(stdl_fchk_parser_skip_intro(lx));

    stdl_wavefunction * wf = NULL;
    stdl_basis * bs = NULL;

    STDL_OK(stdl_fchk_parser_extract(&wf, &bs, lx));

    stdl_context* ctx = NULL;
    STDL_OK(stdlite_context_new(&ctx, wf, bs, 2.0, 4.0, 7. / 27.212, 1.0));

    TEST_ASSERT_EQUAL_INT(ctx->nmo, 5);
    TEST_ASSERT_EQUAL_INT(ctx->nocc, 3);
    TEST_ASSERT_EQUAL_INT(ctx->nvirt, 2);

    // check that the MO are normalized
    for (size_t i = 0; i < ctx->nmo; ++i) {
        double sum = .0;

        for (size_t j = 0; j < wf->nao; ++j) {
            sum += pow(ctx->C[i * wf->nao + j], 2.0);
        }

        TEST_ASSERT_DOUBLE_WITHIN(1e-8, 1., sum);
    }

    // compute the density matrix
    double* density_mat = NULL;
    STDL_OK(stdl_wavefunction_compute_density(&density_mat, ctx->C, ctx->nocc, ctx->nmo, wf->nao));

    // check that density is symmetric
    for (size_t i = 0; i < wf->nao; ++i) {
        for (size_t j = 0; j <=i ; ++j) {
            TEST_ASSERT_EQUAL_DOUBLE(density_mat[i * wf->nao + j], density_mat[j * wf->nao + i]);
        }
    }

    // count electrons
    double total = .0;
    for(size_t i=0; i < wf->nao; i++)
        total += density_mat[i * wf->nao + i];

    TEST_ASSERT_DOUBLE_WITHIN(1e-8, 2.f * ctx->nocc, total);

    // stdl_matrix_dge_print(wf->nao, 0, density_mat, "D");

    free(density_mat);

    STDL_OK(stdlite_context_delete(ctx));
    STDL_OK(stdl_lexer_delete(lx));

    fclose(f);
}
