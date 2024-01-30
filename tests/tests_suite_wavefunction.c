#include <string.h>
#include <unistd.h>

#include <stdlite/utils/fchk_parser.h>
#include <stdlite/matrix.h>

#include "tests_suite.h"

#include <cblas.h>

// check that the data are correct by computing the Mulliken population.
// It should sum up to the number of electrons
void compute_population_and_check(stdl_wavefunction* wf, int sym) {

    // compute the density matrix
    double* density_mat = NULL;
    STDL_OK(stdl_wavefunction_compute_density(&density_mat, wf->C, wf->nocc, wf->nmo, wf->nao));

    // stdl_matrix_dge_print(original_wf->nao, 0, density_mat, "D");

    // check that density is symmetric
    for (size_t i = 0; i < wf->nao; ++i) {
        for (size_t j = 0; j <=i ; ++j) {
            TEST_ASSERT_EQUAL_DOUBLE(density_mat[i * wf->nao + j], density_mat[j * wf->nao + i]);
        }
    }

    // compute Mulliken population
    // See, e.g., https://doi.org/10.26434/chemrxiv.12722072.v1
    double* mulliken_pop = malloc(wf->nao * wf->nao * sizeof(double));
    TEST_ASSERT_NOT_NULL(mulliken_pop);

    if(!sym) {
        double* tmp = malloc(wf->nao * wf->nao * sizeof(double));
        TEST_ASSERT_NOT_NULL(tmp);

        cblas_dsymm(CblasRowMajor, CblasRight, CblasLower,
                    (int) wf->nao, (int) wf->nao,
                    1.f, wf->S, (int) wf->nao,
                    density_mat, (int) wf->nao,
                    .0, tmp, (int) wf->nao
        );

        // `D*S+S*D = D*S+S^T*D^T = D*S+(D*S)^T`
        for (size_t i = 0; i < wf->nao; ++i) {
            for (size_t j = 0; j < wf->nao; ++j) {
                mulliken_pop[i * wf->nao + j] = .5 * (tmp[i * wf->nao + j] + tmp[j * wf->nao + i]);
            }
        }

        free(tmp);
    } else // S=1, so 1/2*(S*D+D*S) = D
        memcpy(mulliken_pop, density_mat, wf->nao * wf->nao * sizeof(double));

    // stdl_matrix_dge_print(original_wf->nao, 0, mulliken_pop, "1/2*(DS+SD)");

    double total = .0;
    for(size_t i=0; i < wf->nao; i++)
        total += mulliken_pop[i * wf->nao + i];

    TEST_ASSERT_DOUBLE_WITHIN(1e-8, (double) wf->nocc * 2, total);

    free(mulliken_pop);
    free(density_mat);
}

void test_content_ok() {

    char* fchk_path = "../tests/test_files/water_sto3g.fchk";

    FILE* f = fopen(fchk_path, "r");
    TEST_ASSERT_NOT_NULL(f);

    stdl_lexer* lx = NULL;
    STDL_OK(stdl_lexer_new(&lx, f));
    STDL_OK(stdl_fchk_parser_skip_intro(lx));

    stdl_wavefunction * wf = NULL;
    stdl_basis * bs = NULL;
    STDL_OK(stdl_fchk_parser_extract(&wf, &bs, lx));
    STDL_OK(stdl_basis_delete(bs));

    compute_population_and_check(wf, 0);

    STDL_OK(stdl_wavefunction_delete(wf));
    STDL_OK(stdl_lexer_delete(lx));

    fclose(f);
}

void test_orthogonalize_ok() {

    char* fchk_path = "../tests/test_files/water_sto3g.fchk";

    FILE* f = fopen(fchk_path, "r");
    TEST_ASSERT_NOT_NULL(f);

    stdl_lexer* lx = NULL;
    STDL_OK(stdl_lexer_new(&lx, f));
    STDL_OK(stdl_fchk_parser_skip_intro(lx));

    stdl_wavefunction * wf = NULL;
    stdl_basis * bs = NULL;
    STDL_OK(stdl_fchk_parser_extract(&wf, &bs, lx));
    STDL_OK(stdl_basis_delete(bs));

    STDL_OK(stdl_wavefunction_orthogonalize_C(wf->C, wf->S, wf->nmo, wf->nao));

    // check that the MO are normalized
    for (size_t i = 0; i < wf->nmo; ++i) {
        double sum = .0;

        for (size_t j = 0; j < wf->nao; ++j) {
            sum += pow(wf->C[i * wf->nao + j], 2.0);
        }

        TEST_ASSERT_DOUBLE_WITHIN(1e-8, 1., sum);
    }

    // set S to identity
    for (size_t i = 0; i < wf->nao; ++i) {
        for (size_t j = 0; j < wf->nao; ++j)
            wf->S[i * wf->nao + j] = i == j ? 1.:.0;
    }

    compute_population_and_check(wf, 1);

    STDL_OK(stdl_wavefunction_delete(wf));
    STDL_OK(stdl_lexer_delete(lx));

    fclose(f);
}


void test_remove_mo_ok() {

    char* fchk_path = "../tests/test_files/water_sto3g.fchk";

    FILE* f = fopen(fchk_path, "r");
    TEST_ASSERT_NOT_NULL(f);

    stdl_lexer* lx = NULL;
    STDL_OK(stdl_lexer_new(&lx, f));
    STDL_OK(stdl_fchk_parser_skip_intro(lx));

    stdl_wavefunction * wf = NULL;
    stdl_basis * bs = NULL;
    STDL_OK(stdl_fchk_parser_extract(&wf, &bs, lx));
    STDL_OK(stdl_basis_delete(bs));

    // "remove" the two first MO
    double* tmp = wf->C;

    wf->C = wf->C + 2 * wf->nao;
    wf->nmo = 5;
    wf->nocc = 3;

    compute_population_and_check(wf, 0);

    // put back C
    wf->C = tmp;

    STDL_OK(stdl_wavefunction_delete(wf));
    STDL_OK(stdl_lexer_delete(lx));

    fclose(f);
}
