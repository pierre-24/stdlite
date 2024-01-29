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

    stdl_matrix_dge_print(wf->nao, ctx->nmo, ctx->C, "CA");

    // check that the MO are normalized
    for (size_t i = 0; i < ctx->nmo; ++i) {
        double sum = .0;

        for (size_t j = 0; j < wf->nao; ++j) {
            sum += pow(ctx->C[i * wf->nao + j], 2.0);
        }

        TEST_ASSERT_DOUBLE_WITHIN(1e-8, 1., sum);
    }

    STDL_OK(stdlite_context_delete(ctx));

    STDL_OK(stdl_lexer_delete(lx));

    fclose(f);
}
