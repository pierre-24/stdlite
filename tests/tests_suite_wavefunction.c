#include <string.h>
#include <cblas.h>

#include <stdlite/utils/fchk_parser.h>
#include <stdlite/utils/matrix.h>
#include <stdlite/helpers.h>
#include <stdlite/utils/molden_parser.h>

#include "tests_suite.h"

void setUp(void) {
    stdl_set_debug_level(-1);
}

// check that the data are correct by computing the Mulliken population.
// It should sum up to the number of electrons (i.e., 2 * wf->nocc), within `prec`.
void compute_population_and_check(stdl_wavefunction *wf, int sym, double prec) {

    // compute the density matrix
    double* P = malloc(STDL_MATRIX_SP_SIZE(wf->nao) * sizeof(double));
    TEST_ASSERT_NOT_NULL(P);

    ASSERT_STDL_OK(stdl_wavefunction_compute_density_dsp(wf->nocc, wf->nmo, wf->nao, wf->C, P));

    // stdl_matrix_dsp_print(wf->nao, P, "P");

    // compute Mulliken population as `M = 1/2*(P*S+S*P)`
    // See, e.g., https://doi.org/10.26434/chemrxiv.12722072.v1
    double* mulliken_pop = malloc(STDL_MATRIX_SP_SIZE(wf->nao) * sizeof(double));
    TEST_ASSERT_NOT_NULL(mulliken_pop);

    if(!sym) {
        double* Ssy = malloc(wf->nao * wf->nao * sizeof(double));
        TEST_ASSERT_NOT_NULL(Ssy);

        // stdl_matrix_dsp_print(wf->nao, wf->S, "S");

        stdl_matrix_dsp_blowsy(wf->nao, 'L', wf->S, Ssy);

        double* Pge = malloc(wf->nao * wf->nao * sizeof(double));
        TEST_ASSERT_NOT_NULL(Pge);

        stdl_matrix_dsp_blowge(1, wf->nao, P, Pge);

        double* tmp = malloc(wf->nao * wf->nao * sizeof(double));
        TEST_ASSERT_NOT_NULL(tmp);

        cblas_dsymm(CblasRowMajor, CblasRight, CblasLower,
                    (int) wf->nao, (int) wf->nao,
                    1.f, Ssy, (int) wf->nao,
                    Pge, (int) wf->nao,
                    .0, tmp, (int) wf->nao
        );

        // `1/2*(P*S+S*P) = 1/2*(P*S+S^T*P^T) = 1/2*(P*S+(P*S)^T)`
        for (size_t mu = 0; mu < wf->nao; ++mu) {
            for (size_t nu = 0; nu < wf->nao; ++nu) {
                mulliken_pop[STDL_MATRIX_SP_IDX(mu, nu)] = .5 * (tmp[mu * wf->nao + nu] + tmp[nu * wf->nao + mu]);
            }
        }

        STDL_FREE_ALL(tmp, Ssy, Pge);
    } else // S=1, so 1/2*(S*D+D*S) = D
        memcpy(mulliken_pop, P, STDL_MATRIX_SP_SIZE(wf->nao) * sizeof(double));

    // stdl_matrix_dsp_print(wf->nao, mulliken_pop, "1/2*(PS+SP)");

    double total = .0;
    for(size_t mu=0; mu < wf->nao; mu++)
        total += mulliken_pop[STDL_MATRIX_SP_IDX(mu, mu)];

    TEST_ASSERT_DOUBLE_WITHIN(prec, (double) wf->nocc * 2, total);

    free(mulliken_pop);
    free(P);
}

void test_content_sto3g_ok() {

    stdl_wavefunction * wf = NULL;
    stdl_basis * bs = NULL;
    read_fchk("../tests/test_files/water_sto3g.fchk", &wf, &bs);
    ASSERT_STDL_OK(stdl_basis_delete(bs));

    compute_population_and_check(wf, 0, 1e-6);

    ASSERT_STDL_OK(stdl_wavefunction_delete(wf));
}

void test_content_631g_ok() {
    stdl_wavefunction * wf = NULL;
    stdl_basis * bs = NULL;
    read_fchk("../tests/test_files/water_631g.fchk", &wf, &bs);
    ASSERT_STDL_OK(stdl_basis_delete(bs));

    compute_population_and_check(wf, 0, 1e-6);

    ASSERT_STDL_OK(stdl_wavefunction_delete(wf));
}

void test_content_cart_6d10f_ok() {
    stdl_wavefunction * wf = NULL;
    stdl_basis * bs = NULL;
    read_fchk("../tests/test_files/water_631gdf.fchk", &wf, &bs);
    ASSERT_STDL_OK(stdl_basis_delete(bs));

    compute_population_and_check(wf, 0, 1e-6);

    ASSERT_STDL_OK(stdl_wavefunction_delete(wf));
}

void test_content_cart_5d7f_ok() {
    stdl_wavefunction * wf = NULL;
    stdl_basis * bs = NULL;
    read_fchk("../tests/test_files/water_631gdf_sph.fchk", &wf, &bs);
    ASSERT_STDL_OK(stdl_basis_delete(bs));

    compute_population_and_check(wf, 0, 1e-6);

    ASSERT_STDL_OK(stdl_wavefunction_delete(wf));
}

void test_molden_cart() {
    char* molden_path = "../tests/test_files/water_631gdf.molden";

    FILE* f = fopen(molden_path, "r");
    TEST_ASSERT_NOT_NULL(f);

    stdl_lexer* lx = NULL;
    ASSERT_STDL_OK(stdl_lexer_new(f, &lx));

    stdl_wavefunction* wf = NULL;
    stdl_basis* bs = NULL;

    ASSERT_STDL_OK(stdl_molden_parser_extract(lx, &wf, &bs));

    fclose(f);
    ASSERT_STDL_OK(stdl_lexer_delete(lx));

    ASSERT_STDL_OK(stdl_basis_delete(bs));

    compute_population_and_check(wf, 0, 1e-4); /* LCAO coefficients are not very precise in this MOLDEN */

    ASSERT_STDL_OK(stdl_wavefunction_delete(wf));
}

void test_molden_sph() {
    char* molden_path = "../tests/test_files/water_631gdf_sph.molden";

    FILE* f = fopen(molden_path, "r");
    TEST_ASSERT_NOT_NULL(f);

    stdl_lexer* lx = NULL;
    ASSERT_STDL_OK(stdl_lexer_new(f, &lx));

    stdl_wavefunction* wf = NULL;
    stdl_basis* bs = NULL;

    ASSERT_STDL_OK(stdl_molden_parser_extract(lx, &wf, &bs));

    fclose(f);
    ASSERT_STDL_OK(stdl_lexer_delete(lx));

    ASSERT_STDL_OK(stdl_basis_delete(bs));

    compute_population_and_check(wf, 0, 1e-4); /* LCAO coefficients are not very precise in this MOLDEN */

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
    read_fchk("../tests/test_files/water_631gdf_sph.fchk", &wf, &bs);

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

    compute_population_and_check(wf, 1, 1e-6);

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

    compute_population_and_check(wf, 0, 1e-6);

    // put back C
    wf->C = tmp;

    ASSERT_STDL_OK(stdl_wavefunction_delete(wf));
}

void test_ao_to_mo() {
    stdl_wavefunction * wf = NULL;
    stdl_basis * bs = NULL;
    read_fchk("../tests/test_files/water_sto3g.fchk", &wf, &bs);

    double* prop_ao_sp = malloc(STDL_MATRIX_SP_SIZE(wf->nao) * sizeof(double));
    TEST_ASSERT_NOT_NULL(prop_ao_sp);

    for (size_t mu = 0; mu < wf->nao; ++mu) {
        for (size_t nu = 0; nu < wf->nao; ++nu) {
            prop_ao_sp[STDL_MATRIX_SP_IDX(mu, nu)] = (double) mu;
        }
    }

    double* prop_mo_sp = malloc(STDL_MATRIX_SP_SIZE(wf->nmo) * sizeof(double));
    TEST_ASSERT_NOT_NULL(prop_mo_sp);

    ASSERT_STDL_OK(stdl_wavefunction_dsp_ao_to_dsp_mo(wf->nao, wf->nmo, wf->C, prop_ao_sp, prop_mo_sp));

    double* prop_ao_ge = malloc(wf->nao * wf->nao * sizeof(double ));
    TEST_ASSERT_NOT_NULL(prop_ao_ge);

    stdl_matrix_dsp_blowge(1, wf->nao, prop_ao_sp, prop_ao_ge);

    // test that prop in AO basis are the same
    for (size_t mu = 0; mu < wf->nao; ++mu) {
        for (size_t nu = 0; nu < wf->nao; ++nu) {
            TEST_ASSERT_EQUAL_DOUBLE(prop_ao_sp[STDL_MATRIX_SP_IDX(mu, nu)], prop_ao_ge[mu * wf->nao + nu]);
            TEST_ASSERT_EQUAL_DOUBLE(prop_ao_sp[STDL_MATRIX_SP_IDX(mu, nu)], prop_ao_ge[nu * wf->nao + mu]);
        }
    }

    double* prop_mo_ge = malloc(wf->nmo * wf->nmo * sizeof(double ));
    TEST_ASSERT_NOT_NULL(prop_mo_ge);

    ASSERT_STDL_OK(stdl_wavefunction_dge_ao_to_dge_mo(wf->nao, wf->nmo, wf->C, prop_ao_ge, prop_mo_ge));

    // stdl_matrix_dsp_print(wf->nmo, prop_mo_sp, "P");
    // stdl_matrix_dge_print(wf->nmo, wf->nmo, prop_mo_ge, "P'");

    // test that prop in MO basis are the same
    for (size_t p = 0; p < wf->nmo; ++p) {
        for (size_t q = 0; q < wf->nmo; ++q) {
            TEST_ASSERT_EQUAL_DOUBLE(prop_mo_sp[STDL_MATRIX_SP_IDX(p, q)], prop_mo_ge[p * wf->nmo + q]);
            TEST_ASSERT_EQUAL_DOUBLE(prop_mo_sp[STDL_MATRIX_SP_IDX(p, q)], prop_mo_ge[q * wf->nmo + p]);
        }
    }

    STDL_FREE_ALL(prop_ao_sp, prop_mo_sp, prop_ao_ge, prop_mo_ge);

    ASSERT_STDL_OK(stdl_basis_delete(bs));
    ASSERT_STDL_OK(stdl_wavefunction_delete(wf));
}

void test_dipole() {
    stdl_wavefunction * wf = NULL;
    stdl_basis * bs = NULL;
    read_fchk("../tests/test_files/water_sto3g.fchk", &wf, &bs);

    double* dipoles_sp = malloc(3 * STDL_MATRIX_SP_SIZE(wf->nao) * sizeof(double));
    TEST_ASSERT_NOT_NULL(dipoles_sp);

    // compute dipole integrals
    ASSERT_STDL_OK(stdl_basis_dsp_dipole(bs, dipoles_sp));

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
    double* dipole_z_ge = malloc(wf->nao * wf->nao * sizeof(double));
    TEST_ASSERT_NOT_NULL(dipole_z_ge);

    stdl_matrix_dsp_blowge(1, wf->nao, dipoles_sp + 2 * STDL_MATRIX_SP_SIZE(wf->nao), dipole_z_ge);

    // compute through density
    // dipz = sum_mu (P * D)_mu,mu
    double* P = malloc(STDL_MATRIX_SP_SIZE(wf->nao) * sizeof(double ));
    TEST_ASSERT_NOT_NULL(P);

    ASSERT_STDL_OK(stdl_wavefunction_compute_density_dsp(wf->nocc, wf->nmo, wf->nao, wf->C, P));

    double* Pge = malloc(wf->nao * wf->nao * sizeof(double ));
    TEST_ASSERT_NOT_NULL(Pge);

    stdl_matrix_dsp_blowge(1, wf->nao, P, Pge);

    double* result = malloc(wf->nao * wf->nao * sizeof(double ));
    TEST_ASSERT_NOT_NULL(result);

    cblas_dsymm(CblasRowMajor, CblasRight, CblasLower,
                (int) wf->nao, (int) wf->nao,
                1.f, dipole_z_ge, (int) wf->nao,
                Pge, (int) wf->nao,
                .0, result, (int) wf->nao
    );

    double dipz2 = .0;
    for (size_t mu = 0; mu < wf->nao; ++mu) {
        dipz2 += result[mu * wf->nao + mu];
    }

    TEST_ASSERT_DOUBLE_WITHIN(1e-6, dipz1, dipz2);

    // compute through MO basis
    // dipz = sum_p n_p * D_pp (where n_p is the occupancy of MO p)
    double* dipole_z_mo_sp = malloc(STDL_MATRIX_SP_SIZE(wf->nmo) * sizeof(double));
    TEST_ASSERT_NOT_NULL(dipole_z_mo_sp);

    ASSERT_STDL_OK(stdl_wavefunction_dsp_ao_to_dsp_mo(wf->nao, wf->nmo, wf->C, dipoles_sp + 2 * STDL_MATRIX_SP_SIZE(wf->nao), dipole_z_mo_sp));

    double dipz3 = .0;
    for (size_t p = 0; p < wf->nocc; ++p) {
        dipz3 += 2 * dipole_z_mo_sp[STDL_MATRIX_SP_IDX(p, p)];
    }

    TEST_ASSERT_DOUBLE_WITHIN(1e-6, dipz1, dipz3);

    STDL_FREE_ALL(dipoles_sp, dipole_z_ge, dipole_z_mo_sp, P, Pge, result);

    ASSERT_STDL_OK(stdl_basis_delete(bs));
    ASSERT_STDL_OK(stdl_wavefunction_delete(wf));
}
