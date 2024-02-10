#include <string.h>

#include <stdlite/utils/fchk_parser.h>
#include <stdlite/utils/matrix.h>
#include <stdlite/context.h>
#include <stdlite/helpers.h>
#include <cblas.h>

#include "tests_suite.h"

void setUp() {
    stdl_set_debug_level(-1);
}

void read_fchk(char* fchk_path,stdl_wavefunction** wf, stdl_basis** bs) {

    FILE* f = fopen(fchk_path, "r");
    TEST_ASSERT_NOT_NULL(f);

    stdl_lexer* lx = NULL;
    ASSERT_STDL_OK(stdl_lexer_new(f, &lx));
    ASSERT_STDL_OK(stdl_fchk_parser_skip_intro(lx));

    ASSERT_STDL_OK(stdl_fchk_parser_extract(lx, wf, bs));

    ASSERT_STDL_OK(stdl_lexer_delete(lx));

    fclose(f);
}

void test_context_select_MO_ok() {
    stdl_wavefunction * wf = NULL;
    stdl_basis * bs = NULL;
    read_fchk("../tests/test_files/water_631g.fchk", &wf, &bs);

    stdl_context* ctx = NULL;
    ASSERT_STDL_OK(stdl_context_new(wf, bs, 2.0, 4.0, 10. / 27.212, 1e-4, 1.0, &ctx));

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
    double* P = malloc(STDL_MATRIX_SP_SIZE(wf->nao) * sizeof(double));
    TEST_ASSERT_NOT_NULL(P);

    ASSERT_STDL_OK(stdl_wavefunction_compute_density_dsp(ctx->nocc, ctx->nmo, wf->nao, ctx->C, P));

    // count electrons
    double total = .0;
    for(size_t mu=0; mu < wf->nao; mu++)
        total += P[STDL_MATRIX_SP_IDX(mu, mu)];

    TEST_ASSERT_DOUBLE_WITHIN(1e-6, 2.f * ctx->nocc, total);

    // stdl_matrix_dge_print(wf->nao, 0, P, "D");

    free(P);

    ASSERT_STDL_OK(stdl_context_delete(ctx));
}

void test_dipole() {
    stdl_wavefunction * wf = NULL;
    stdl_basis * bs = NULL;
    read_fchk("../tests/test_files/water_631g.fchk", &wf, &bs);

    stdl_context* ctx = NULL;
    ASSERT_STDL_OK(stdl_context_new(wf, bs, 2.0, 4.0, 12. / 27.212, 1e-4, 1.0, &ctx));

    double* dipoles_sp = malloc(3 * STDL_MATRIX_SP_SIZE(wf->nao) * sizeof(double));
    TEST_ASSERT_NOT_NULL(dipoles_sp);

    // compute dipole integrals
    ASSERT_STDL_OK(stdl_basis_dsp_dipole(bs, dipoles_sp));

    // compute explicitly the electronic dipole moment along z
    double dipz1 = .0;
    for (size_t p = 0; p < ctx->nocc; ++p) {
        for (size_t mu = 0; mu < wf->nao; ++mu) {
            for (size_t nu = 0; nu < wf->nao; ++nu) {
                dipz1 += 2 * ctx->C_orig[p * wf->nao + mu] * ctx->C_orig[p * wf->nao + nu] * dipoles_sp[2 * STDL_MATRIX_SP_SIZE(wf->nao) + STDL_MATRIX_SP_IDX(mu, nu)];
            }
        }
    }

    // blow the matrix
    double* dipole_z_sy = malloc(wf->nao * wf->nao * sizeof(double));
    TEST_ASSERT_NOT_NULL(dipole_z_sy);

    stdl_matrix_dsp_blowsy(wf->nao, 'L', dipoles_sp + 2 * STDL_MATRIX_SP_SIZE(wf->nao), dipole_z_sy);

    // compute through density
    // dipz = sum_mu (P * D)_mu,mu
    double* P = malloc(STDL_MATRIX_SP_SIZE(wf->nao) * sizeof(double ));
    TEST_ASSERT_NOT_NULL(P);

    ASSERT_STDL_OK(stdl_wavefunction_compute_density_dsp(ctx->nocc, ctx->nmo, wf->nao, ctx->C_orig, P));

    double* Pge = malloc(wf->nao * wf->nao * sizeof(double ));
    TEST_ASSERT_NOT_NULL(Pge);

    stdl_matrix_dsp_blowge(1, wf->nao, P, Pge);

    double* result = malloc(wf->nao * wf->nao * sizeof(double ));
    TEST_ASSERT_NOT_NULL(result);

    cblas_dsymm(CblasRowMajor, CblasRight, CblasLower,
                (int) wf->nao, (int) wf->nao,
                1.f, dipole_z_sy, (int) wf->nao,
                Pge, (int) wf->nao,
                .0, result, (int) wf->nao
    );

    double dipz2 = .0;
    for (size_t mu = 0; mu < wf->nao; ++mu) {
        dipz2 += result[mu * wf->nao + mu];
    }

    TEST_ASSERT_DOUBLE_WITHIN(1e-6, dipz1, dipz2);

    // compute through MO basis
    // dipz = sum_p n_p * D_pp = 2 * sum^occ_p D_pp (where n_p is the occupancy of MO p)
    double* dipole_z_mo_sy = malloc(ctx->nocc * ctx->nocc * sizeof(double));
    TEST_ASSERT_NOT_NULL(dipole_z_mo_sy);

    // only use occupied MOs
    ASSERT_STDL_OK(stdl_wavefunction_dsy_ao_to_mo(wf->nao, ctx->nocc, ctx->C_orig, dipole_z_sy, dipole_z_mo_sy));

    double dipz3 = .0;
    for (size_t p = 0; p < ctx->nocc; ++p) {
        dipz3 += 2 * dipole_z_mo_sy[p * ctx->nocc + p];
    }

    TEST_ASSERT_DOUBLE_WITHIN(1e-6, dipz1, dipz3);

    STDL_FREE_ALL(dipoles_sp, dipole_z_sy, dipole_z_mo_sy, P, Pge, result);

    ASSERT_STDL_OK(stdl_context_delete(ctx));
}


void test_context_select_csfs_ok() {
    stdl_wavefunction * wf = NULL;
    stdl_basis * bs = NULL;
    read_fchk("../tests/test_files/water_631g.fchk", &wf, &bs);

    stdl_context* ctx = NULL;
    ASSERT_STDL_OK(stdl_context_new(wf, bs, 2.0, 4.0, 12. / 27.212, 1e-4, 1.0, &ctx));

    TEST_ASSERT_EQUAL_INT(ctx->nmo,7);
    TEST_ASSERT_EQUAL_INT(ctx->nocc, 4);

    TEST_ASSERT_EQUAL_INT(0, ctx->ncsfs);

    ASSERT_STDL_OK(stdl_context_select_csfs_monopole(ctx, 0));

    TEST_ASSERT_EQUAL_INT(10, ctx->ncsfs);
    TEST_ASSERT_NULL(ctx->B);

    // check that energies are in increasing order
    for (size_t kia = 0; kia < ctx->ncsfs; ++kia) {
        TEST_ASSERT_EQUAL_FLOAT(ctx->A[STDL_MATRIX_SP_IDX(kia, kia)], ctx->ecsfs[kia]);

        if(kia > 0)
            TEST_ASSERT_TRUE(ctx->ecsfs[kia-1] <= ctx->ecsfs[kia]);
    }

    // stdl_matrix_ssp_print(ctx->ncsfs, ctx->A, "A");

    ASSERT_STDL_OK(stdl_context_delete(ctx));
}
