#include <string.h>
#include <unistd.h>
#include <lapacke.h>

#include <stdlite/utils/fchk_parser.h>
#include <stdlite/utils/matrix.h>
#include <stdlite/helpers.h>

#include "tests_suite.h"

#include <cblas.h>

// check that the data are correct by computing the Mulliken population.
// It should sum up to the number of electrons (i.e., 2 * wf->nocc).
void compute_population_and_check(stdl_wavefunction* wf, int sym) {

    // compute the density matrix
    double* P = malloc(wf->nao * wf->nao * sizeof(double));
    TEST_ASSERT_NOT_NULL(P);

    ASSERT_STDL_OK(stdl_wavefunction_compute_density_dsy(wf->nocc, wf->nmo, wf->nao, wf->C, P));

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

// Read a FCHK file and extract both basis set and wavefunction
void read_fchk(char* fchk_path, stdl_wavefunction** wf, stdl_basis** bs) {

    FILE* f = fopen(fchk_path, "r");
    TEST_ASSERT_NOT_NULL(f);

    stdl_lexer* lx = NULL;
    ASSERT_STDL_OK(stdl_lexer_new(f, &lx));
    ASSERT_STDL_OK(stdl_fchk_parser_skip_intro(lx));

    ASSERT_STDL_OK(stdl_fchk_parser_extract(lx, wf, bs));

    ASSERT_STDL_OK(stdl_lexer_delete(lx));

    fclose(f);
}

void test_content_ok() {
    stdl_wavefunction * wf = NULL;
    stdl_basis * bs = NULL;
    read_fchk("../tests/test_files/water_sto3g.fchk", &wf, &bs);
    ASSERT_STDL_OK(stdl_basis_delete(bs));

    compute_population_and_check(wf, 0);

    ASSERT_STDL_OK(stdl_wavefunction_delete(wf));
}

void test_sqrtS_ok() {
    stdl_wavefunction * wf = NULL;
    stdl_basis * bs = NULL;
    read_fchk("../tests/test_files/water_sto3g.fchk", &wf, &bs);

    ASSERT_STDL_OK(stdl_basis_delete(bs));

    // sy
    double* sqrtS = malloc(wf->nao * wf->nao * sizeof(double));
    TEST_ASSERT_NOT_NULL(sqrtS);

    stdl_matrix_dsp_sqrt_sy(wf->nao, wf->S, sqrtS);

    // Compare in-place procedure and sy
    double* Sp = malloc(STDL_MATRIX_SP_SIZE(wf->nao) * sizeof(double));
    memcpy(Sp, wf->S, STDL_MATRIX_SP_SIZE(wf->nao) * sizeof(double));
    TEST_ASSERT_NOT_NULL(Sp);

    stdl_matrix_dsp_sqrt(wf->nao, Sp);

    for (size_t mu = 0; mu < wf->nao; ++mu) {
        for (size_t nu = 0; nu <= mu; ++nu) {
            TEST_ASSERT_DOUBLE_WITHIN(1e-8, Sp[STDL_MATRIX_SP_IDX(mu, nu)], sqrtS[mu * wf->nao + nu]);
        }
    }

    // compute S' = S^1/2 * S^1/2 and compare to S
    double* reS = malloc(wf->nao * wf->nao * sizeof(double));
    TEST_ASSERT_NOT_NULL(reS);

    cblas_dsymm(CblasRowMajor, CblasLeft, CblasLower,
                (int) wf->nao, (int) wf->nao, 1.0, sqrtS, (int) wf->nao,
                sqrtS, (int) wf->nao,
                .0, reS, (int) wf->nao
                );

    for(size_t i = 0; i < wf->nao; i++) {
        for (size_t j = 0; j <= i; ++j) {
            TEST_ASSERT_DOUBLE_WITHIN(1e-8, reS[i * wf->nao + j], wf->S[STDL_MATRIX_SP_IDX(i, j)]);
        }
    }

    STDL_FREE_ALL(sqrtS, reS, Sp);
    ASSERT_STDL_OK(stdl_wavefunction_delete(wf));
}

void test_orthogonalize_ok() {
    stdl_wavefunction * wf = NULL;
    stdl_basis * bs = NULL;
    read_fchk("../tests/test_files/water_sto3g.fchk", &wf, &bs);

    ASSERT_STDL_OK(stdl_basis_delete(bs));

    ASSERT_STDL_OK(stdl_wavefunction_orthogonalize_C_dge(wf->nmo, wf->nao, wf->S, wf->C));

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

    ASSERT_STDL_OK(stdl_wavefunction_delete(wf));
}


void test_remove_mo_ok() {
    stdl_wavefunction * wf = NULL;
    stdl_basis * bs = NULL;
    read_fchk("../tests/test_files/water_sto3g.fchk", &wf, &bs);

    ASSERT_STDL_OK(stdl_basis_delete(bs));

    // "remove" the two first MO
    double* tmp = wf->C;

    wf->C = wf->C + 2 * wf->nao;
    wf->nmo = 5;
    wf->nocc = 3;

    compute_population_and_check(wf, 0);

    // put back C
    wf->C = tmp;

    ASSERT_STDL_OK(stdl_wavefunction_delete(wf));
}

void test_dipole() {
    stdl_wavefunction * wf = NULL;
    stdl_basis * bs = NULL;
    read_fchk("../tests/test_files/water_sto3g.fchk", &wf, &bs);

    double* dipoles_sp = malloc(3 * STDL_MATRIX_SP_SIZE(wf->nao) * sizeof(double));
    TEST_ASSERT_NOT_NULL(dipoles_sp);

    // compute dipole integrals
    ASSERT_STDL_OK(stdl_basis_compute_dsp_dipole(bs, dipoles_sp));

    // compute explicitly the electronic dipole moment along z
    double dipz1 = .0;
    for (size_t p = 0; p < wf->nocc; ++p) {
        for (size_t mu = 0; mu < wf->nao; ++mu) {
            for (size_t nu = 0; nu < wf->nao; ++nu) {
                dipz1 += 2 * wf->C[p * wf->nao + mu] * wf->C[p * wf->nao + nu] * dipoles_sp[2 * STDL_MATRIX_SP_SIZE(wf->nao) + STDL_MATRIX_SP_IDX(mu, nu)];
            }
        }
    }

    TEST_ASSERT_DOUBLE_WITHIN(1e-6, -0.675072, dipz1);

    // blow the matrix
    double* dipole_z_sy = malloc(wf->nao * wf->nao * sizeof(double));
    TEST_ASSERT_NOT_NULL(dipole_z_sy);

    stdl_matrix_dsp_blowsy(wf->nao, 'L', dipoles_sp + 2 * STDL_MATRIX_SP_SIZE(wf->nao), dipole_z_sy);

    // compute through density
    // dipz = sum_mu (P * D)_mu,mu
    double* P = malloc(wf->nao * wf->nao * sizeof(double ));
    TEST_ASSERT_NOT_NULL(P);

    double* result = malloc(wf->nao * wf->nao * sizeof(double ));
    TEST_ASSERT_NOT_NULL(result);

    ASSERT_STDL_OK(stdl_wavefunction_compute_density_dsy(wf->nocc, wf->nmo, wf->nao, wf->C, P));

    cblas_dsymm(CblasRowMajor, CblasRight, CblasLower,
                (int) wf->nao, (int) wf->nao,
                1.f, dipole_z_sy, (int) wf->nao,
                P, (int) wf->nao,
                .0, result, (int) wf->nao
    );

    double dipz2 = .0;
    for (size_t mu = 0; mu < wf->nao; ++mu) {
        dipz2 += result[mu * wf->nao + mu];
    }

    TEST_ASSERT_DOUBLE_WITHIN(1e-6, dipz1, dipz2);

    // compute through MO basis
    // dipz = sum_p n_p * D_pp (where n_p is the occupancy of MO p)
    double* dipole_z_mo_sy = malloc(wf->nmo * wf->nmo * sizeof(double));
    TEST_ASSERT_NOT_NULL(dipole_z_mo_sy);

    ASSERT_STDL_OK(stdl_wavefunction_dsy_ao_to_mo(wf->nao, wf->nmo, wf->C, dipole_z_sy, dipole_z_mo_sy));

    double dipz3 = .0;
    for (size_t p = 0; p < wf->nocc; ++p) {
        dipz3 += 2 * dipole_z_mo_sy[p * wf->nmo + p];
    }

    TEST_ASSERT_DOUBLE_WITHIN(1e-6, dipz1, dipz3);

    STDL_FREE_ALL(dipoles_sp, dipole_z_sy, dipole_z_mo_sy, P, result);

    ASSERT_STDL_OK(stdl_basis_delete(bs));
    ASSERT_STDL_OK(stdl_wavefunction_delete(wf));
}
