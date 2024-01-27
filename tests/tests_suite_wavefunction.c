#include <string.h>
#include <unistd.h>

#include <stdlite/utils/fchk_parser.h>
#include <stdlite/matrix.h>

#include "tests_suite.h"

#include <cblas.h>

void _check_wavefunction(stdl_wavefunction* wf) {

    // compute the density matrix, assuming a closed-shell WF.
    // D is `double[nao*nao]`
    // Elements: `D_ij = sum_k n_i * C_ik * C_jk ~ 2 * sum_k^occ * C_ik (C^T)_kj`.
    double* density_mat = malloc(wf->nao * wf->nao * sizeof(double));
    TEST_ASSERT_NOT_NULL(density_mat);

    cblas_dgemm(
            CblasRowMajor, CblasTrans, CblasNoTrans,
            (int) wf->nao, (int) wf->nao, (int) wf->nelec / 2,
            2.f, wf->C, (int) wf->nmo,
            wf->C, (int) wf->nmo,
            .0, density_mat, (int) wf->nao
    );

    // check that density is symmetric
    for (size_t i = 0; i < wf->nao; ++i) {
        for (size_t j = 0; j <=i ; ++j) {
            TEST_ASSERT_EQUAL_DOUBLE(density_mat[i * wf->nao + j], density_mat[j * wf->nao + i]);
        }
    }

    // compute population
    // `N = tr(DS)`.
    // For full Mulliken population  (i.e., `PS-SP`), see, e.g., https://doi.org/10.26434/chemrxiv.12722072.v1
    double* mulliken_pop = malloc(wf->nao * wf->nao * sizeof(double));
    TEST_ASSERT_NOT_NULL(mulliken_pop);

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                (int) wf->nao, (int) wf->nao, (int) wf->nao,
                1.f, density_mat, (int) wf->nmo,
                wf->S, (int) wf->nao,
                .0, mulliken_pop, (int) wf->nao
    );

    double total = .0;
    for(size_t i=0; i < wf->nao; i++)
        total += mulliken_pop[i * wf->nao + i];

    TEST_ASSERT_DOUBLE_WITHIN(1e-8, (double) wf->nelec, total);

    free(mulliken_pop);
    free(density_mat);
}

void test_content_ok() {
    char cwd[512], fchk_path[1024];
    TEST_ASSERT_NOT_NULL(getcwd(cwd, 512));

    sprintf(fchk_path, "%s/../tests/test_files/water_sto3g.fchk", cwd);

    FILE* f = fopen(fchk_path, "r");
    TEST_ASSERT_NOT_NULL(f);

    stdl_lexer* lx = NULL;
    STDL_OK(stdl_lexer_new(&lx, f));
    STDL_OK(stdl_fchk_parser_skip_intro(lx));

    stdl_wavefunction * wf = NULL;
    stdl_basis * bs = NULL;
    STDL_OK(stdl_fchk_parser_extract(&wf, &bs, lx));
    STDL_OK(stdl_basis_delete(bs));

    _check_wavefunction(wf);

    STDL_OK(stdl_wavefunction_delete(wf));
    STDL_OK(stdl_lexer_delete(lx));

    fclose(f);
}
