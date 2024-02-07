#include <string.h>
#include <unistd.h>
#include <lapacke.h>

#include <stdlite/utils/fchk_parser.h>
#include <stdlite/utils/matrix.h>
#include <stdlite/helpers.h>

#include "tests_suite.h"

#include <cblas.h>

// check that the data are correct by computing the Mulliken population.
// It should sum up to the number of electrons
void compute_population_and_check(stdl_wavefunction* wf, int sym) {

    // compute the density matrix
    double* P = malloc(wf->nao * wf->nao * sizeof(double));
    TEST_ASSERT_NOT_NULL(P);

    STDL_OK(stdl_wavefunction_compute_density_dsy(wf->nocc, wf->nmo, wf->nao, wf->C, P));

    // stdl_matrix_dge_print(wf->nao, wf->nao, P, "P");

    // check that density is symmetric
    for (size_t i = 0; i < wf->nao; ++i) {
        for (size_t j = 0; j <=i ; ++j) {
            TEST_ASSERT_EQUAL_DOUBLE(P[i * wf->nao + j], P[j * wf->nao + i]);
        }
    }

    // compute Mulliken population as `M = 1/2*(P*S+S*P)`
    // See, e.g., https://doi.org/10.26434/chemrxiv.12722072.v1
    double* mulliken_pop = malloc(wf->nao * wf->nao * sizeof(double));
    TEST_ASSERT_NOT_NULL(mulliken_pop);

    if(!sym) {
        double* Ssy = malloc(wf->nao * wf->nao * sizeof(double));
        TEST_ASSERT_NOT_NULL(Ssy);

        LAPACKE_dtpttr(LAPACK_ROW_MAJOR, 'L', (int) wf->nao, wf->S, Ssy, (int) wf->nao);

        double* tmp = malloc(wf->nao * wf->nao * sizeof(double));
        TEST_ASSERT_NOT_NULL(tmp);

        cblas_dsymm(CblasRowMajor, CblasRight, CblasLower,
                    (int) wf->nao, (int) wf->nao,
                    1.f, Ssy, (int) wf->nao,
                    P, (int) wf->nao,
                    .0, tmp, (int) wf->nao
        );

        // `1/2*(P*S+S*P) = 1/2*(P*S+S^T*P^T) = 1/2*(P*S+(P*S)^T)`
        for (size_t i = 0; i < wf->nao; ++i) {
            for (size_t j = 0; j < wf->nao; ++j) {
                mulliken_pop[i * wf->nao + j] = .5 * (tmp[i * wf->nao + j] + tmp[j * wf->nao + i]);
            }
        }

        STDL_FREE_ALL(tmp, Ssy);
    } else // S=1, so 1/2*(S*D+D*S) = D
        memcpy(mulliken_pop, P, wf->nao * wf->nao * sizeof(double));

    // stdl_matrix_dge_print(wf->nao, 0, mulliken_pop, "1/2*(PS+SP)");

    double total = .0;
    for(size_t i=0; i < wf->nao; i++)
        total += mulliken_pop[i * wf->nao + i];

    TEST_ASSERT_DOUBLE_WITHIN(1e-8, (double) wf->nocc * 2, total);

    free(mulliken_pop);
    free(P);
}

void test_content_ok() {

    char* fchk_path = "../tests/test_files/water_sto3g.fchk";

    FILE* f = fopen(fchk_path, "r");
    TEST_ASSERT_NOT_NULL(f);

    stdl_lexer* lx = NULL;
    STDL_OK(stdl_lexer_new(f, &lx));
    STDL_OK(stdl_fchk_parser_skip_intro(lx));

    stdl_wavefunction * wf = NULL;
    stdl_basis * bs = NULL;
    STDL_OK(stdl_fchk_parser_extract(lx, &wf, &bs));
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
    STDL_OK(stdl_lexer_new(f, &lx));
    STDL_OK(stdl_fchk_parser_skip_intro(lx));

    stdl_wavefunction * wf = NULL;
    stdl_basis * bs = NULL;
    STDL_OK(stdl_fchk_parser_extract(lx, &wf, &bs));
    STDL_OK(stdl_basis_delete(bs));

    STDL_OK(stdl_wavefunction_orthogonalize_C_dge(wf->nmo, wf->nao, wf->S, wf->C));

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
        for (size_t j = 0; j <= i; ++j)
            wf->S[STDL_MATRIX_SP_IDX(i, j)] = i == j ? 1.:.0;
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
    STDL_OK(stdl_lexer_new(f, &lx));
    STDL_OK(stdl_fchk_parser_skip_intro(lx));

    stdl_wavefunction * wf = NULL;
    stdl_basis * bs = NULL;
    STDL_OK(stdl_fchk_parser_extract(lx, &wf, &bs));
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
