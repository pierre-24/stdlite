#include <string.h>

#include <stdlite/utils/fchk_parser.h>
#include <stdlite/utils/matrix.h>
#include <stdlite/context.h>
#include <stdlite/helpers.h>
#include <cblas.h>
#include <unistd.h>

#include "tests_suite.h"

void setUp() {
    stdl_set_debug_level(-1);
}

void test_context_select_MO_ok() {
    stdl_wavefunction * wf = NULL;
    stdl_basis * bs = NULL;
    read_fchk("../tests/test_files/water_631g.fchk", &wf, &bs);

    stdl_context* ctx = NULL;
    ASSERT_STDL_OK(stdl_context_new(wf, bs, 2.0, 4.0, 10. / STDL_CONST_AU_TO_EV, 1e-4, 1.0, &ctx));

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

void test_context_select_MO_631gdf_ok() {
    stdl_wavefunction * wf = NULL;
    stdl_basis * bs = NULL;
    read_fchk("../tests/test_files/water_631gdf.fchk", &wf, &bs);

    stdl_context* ctx = NULL;
    ASSERT_STDL_OK(stdl_context_new(wf, bs, 2.0, 4.0, 12. / STDL_CONST_AU_TO_EV, 1e-4, 1.0, &ctx));

    TEST_ASSERT_EQUAL_INT(9, ctx->nmo);
    TEST_ASSERT_EQUAL_INT(4, ctx->nocc);

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
    ASSERT_STDL_OK(stdl_context_new(wf, bs, 2.0, 4.0, 12. / STDL_CONST_AU_TO_EV, 1e-4, 1.0, &ctx));

    double* dipoles_sp = malloc(3 * STDL_MATRIX_SP_SIZE(wf->nao) * sizeof(double));
    TEST_ASSERT_NOT_NULL(dipoles_sp);

    // compute dipole integrals
    ASSERT_STDL_OK(stdl_basis_dsp_dipole(bs, dipoles_sp));

    // compute explicitly the electronic dipole moment along z
    double dipz1 = .0;
    for (size_t p = 0; p < ctx->nocc; ++p) {
        for (size_t mu = 0; mu < wf->nao; ++mu) {
            for (size_t nu = 0; nu < wf->nao; ++nu) {
                dipz1 += 2 * ctx->C_ptr[p * wf->nao + mu] * ctx->C_ptr[p * wf->nao + nu] * dipoles_sp[2 * STDL_MATRIX_SP_SIZE(wf->nao) + STDL_MATRIX_SP_IDX(mu, nu)];
            }
        }
    }

    // blow the matrix
    double* dipole_z_ge = malloc(wf->nao * wf->nao * sizeof(double));
    TEST_ASSERT_NOT_NULL(dipole_z_ge);

    stdl_matrix_dsp_blowge(1, wf->nao, dipoles_sp + 2 * STDL_MATRIX_SP_SIZE(wf->nao), dipole_z_ge);

    // compute through density
    // dipz = tr(P*D)
    double* P = malloc(STDL_MATRIX_SP_SIZE(wf->nao) * sizeof(double ));
    TEST_ASSERT_NOT_NULL(P);

    ASSERT_STDL_OK(stdl_wavefunction_compute_density_dsp(ctx->nocc, ctx->nmo, wf->nao, ctx->C_ptr, P));

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
    // dipz = sum_p n_p * D_pp = 2 * sum^occ_p D_pp (where n_p is the occupancy of MO p)
    double* dipole_z_mo = malloc(STDL_MATRIX_SP_SIZE(ctx->nocc) * sizeof(double));
    TEST_ASSERT_NOT_NULL(dipole_z_mo);

    // only use occupied MOs!
    ASSERT_STDL_OK(stdl_wavefunction_dsp_ao_to_dsp_mo(wf->nao, ctx->nocc, ctx->C_ptr, dipoles_sp + 2 * STDL_MATRIX_SP_SIZE(wf->nao), dipole_z_mo));

    double dipz3 = .0;
    for (size_t p = 0; p < ctx->nocc; ++p) {
        dipz3 += 2 * dipole_z_mo[STDL_MATRIX_SP_IDX(p, p)];
    }

    TEST_ASSERT_DOUBLE_WITHIN(1e-6, dipz1, dipz3);

    STDL_FREE_ALL(dipoles_sp, dipole_z_ge, dipole_z_mo, P, Pge, result);

    ASSERT_STDL_OK(stdl_context_delete(ctx));
}


void test_context_select_csfs_ok() {
    stdl_wavefunction * wf = NULL;
    stdl_basis * bs = NULL;
    read_fchk("../tests/test_files/water_631g.fchk", &wf, &bs);

    stdl_context* ctx = NULL;
    ASSERT_STDL_OK(stdl_context_new(wf, bs, 2.0, 4.0, 12. / STDL_CONST_AU_TO_EV, 1e-4, 1.0, &ctx));

    TEST_ASSERT_EQUAL_INT(ctx->nmo,7);
    TEST_ASSERT_EQUAL_INT(ctx->nocc, 4);

    TEST_ASSERT_EQUAL_INT(0, ctx->ncsfs);

    ASSERT_STDL_OK(stdl_context_select_csfs_monopole(ctx, 1));

    TEST_ASSERT_EQUAL_INT(10, ctx->ncsfs);
    TEST_ASSERT_NULL(ctx->B);

    stdl_matrix_ssp_print(ctx->ncsfs, ctx->A, "A");

    // check that energies are in increasing order
    for (size_t lia = 1; lia < ctx->ncsfs; ++lia) {
        TEST_ASSERT_TRUE(ctx->A[STDL_MATRIX_SP_IDX(lia - 1, lia - 1)] <= ctx->A[STDL_MATRIX_SP_IDX(lia, lia)]);
    }

    // stdl_matrix_ssp_print(ctx->ncsfs, ctx->A, "A");

    ASSERT_STDL_OK(stdl_context_delete(ctx));
}


void test_context_select_csfs_direct_ok() {
    stdl_wavefunction * wf1 = NULL;
    stdl_basis * bs1 = NULL;
    read_fchk("../tests/test_files/water_631g.fchk", &wf1, &bs1);

    stdl_context* ctx1 = NULL;
    ASSERT_STDL_OK(stdl_context_new(wf1, bs1, 2.0, 4.0, 12. / STDL_CONST_AU_TO_EV, 1e-4, 1.0, &ctx1));
    ASSERT_STDL_OK(stdl_context_select_csfs_monopole(ctx1, 1));

    stdl_wavefunction * wf2 = NULL;
    stdl_basis * bs2 = NULL;
    read_fchk("../tests/test_files/water_631g.fchk", &wf2, &bs2);

    stdl_context* ctx2 = NULL;
    ASSERT_STDL_OK(stdl_context_new(wf2, bs2, 2.0, 4.0, 12. / STDL_CONST_AU_TO_EV, 1e-4, 1.0, &ctx2));
    ASSERT_STDL_OK(stdl_context_select_csfs_monopole_direct(ctx2, 1));

    // check that both matrix are identical
    for (size_t kia = 0; kia < ctx1->ncsfs; ++kia) {
        for (size_t kjb = 0; kjb < ctx1->ncsfs; ++kjb) {
            TEST_ASSERT_FLOAT_WITHIN(1e-6, ctx1->A[STDL_MATRIX_SP_IDX(kia, kjb)], ctx2->A[STDL_MATRIX_SP_IDX(kia, kjb)]);
            TEST_ASSERT_FLOAT_WITHIN(1e-6, ctx1->B[STDL_MATRIX_SP_IDX(kia, kjb)], ctx2->B[STDL_MATRIX_SP_IDX(kia, kjb)]);
        }
    }

    ASSERT_STDL_OK(stdl_context_delete(ctx1));
    ASSERT_STDL_OK(stdl_context_delete(ctx2));
}


void test_context_dump_load_h5_ok() {
    stdl_wavefunction * wf1 = NULL;
    stdl_basis * bs1 = NULL;
    read_fchk("../tests/test_files/water_631g.fchk", &wf1, &bs1);

    stdl_context* ctx1 = NULL;
    ASSERT_STDL_OK(stdl_context_new(wf1, bs1, 2.0, 4.0, 12. / STDL_CONST_AU_TO_EV, 1e-4, 1.0, &ctx1));
    ASSERT_STDL_OK(stdl_context_select_csfs_monopole(ctx1, 1));

    // dump & load
    char tmp_path[512];
    sprintf(tmp_path, "%s/tmp.h5", P_tmpdir);

    ASSERT_STDL_OK(stdl_context_dump_h5(ctx1, tmp_path));

    stdl_context* ctx2 = NULL;
    ASSERT_STDL_OK(stdl_context_load_h5(tmp_path, &ctx2));
    unlink(tmp_path);

    // test a few equalities
    TEST_ASSERT_EQUAL(ctx1->original_wf->nao, ctx2->original_wf->nao);
    TEST_ASSERT_EQUAL(ctx1->original_wf->nmo, ctx2->original_wf->nmo);

    TEST_ASSERT_DOUBLE_ARRAY_WITHIN(1e-12, ctx1->original_wf->C, ctx2->original_wf->C, ctx1->original_wf->nmo * ctx1->original_wf->nao);
    TEST_ASSERT_DOUBLE_ARRAY_WITHIN(1e-12, ctx1->original_wf->e, ctx2->original_wf->e, ctx1->original_wf->nao);
    TEST_ASSERT_DOUBLE_ARRAY_WITHIN(1e-12, ctx1->original_wf->S, ctx2->original_wf->S, STDL_MATRIX_SP_SIZE(ctx1->original_wf->nao));

    TEST_ASSERT_EQUAL(ctx1->bs->natm, ctx2->bs->natm);
    TEST_ASSERT_EQUAL(ctx1->bs->nbas, ctx2->bs->nbas);
    TEST_ASSERT_EQUAL(ctx1->bs->env_size, ctx2->bs->env_size);
    TEST_ASSERT_EQUAL_INT_ARRAY(ctx1->bs->atm, ctx2->bs->atm, 6 * ctx1->bs->natm);
    TEST_ASSERT_EQUAL_INT_ARRAY(ctx1->bs->bas, ctx2->bs->bas, 8 * ctx1->bs->nbas);
    TEST_ASSERT_DOUBLE_ARRAY_WITHIN(1e-12, ctx1->bs->env, ctx2->bs->env, ctx1->bs->env_size);

    TEST_ASSERT_EQUAL(ctx1->nmo, ctx2->nmo);
    TEST_ASSERT_EQUAL(ctx1->nocc, ctx2->nocc);
    TEST_ASSERT_EQUAL(ctx1->ncsfs, ctx2->ncsfs);

    TEST_ASSERT_EQUAL_INT_ARRAY(ctx1->csfs, ctx2->csfs, ctx1->ncsfs);
    TEST_ASSERT_FLOAT_ARRAY_WITHIN(1e-6, ctx1->A, ctx2->A, STDL_MATRIX_SP_SIZE(ctx1->ncsfs));
    TEST_ASSERT_FLOAT_ARRAY_WITHIN(1e-6, ctx1->B, ctx2->B, STDL_MATRIX_SP_SIZE(ctx1->ncsfs));
    TEST_ASSERT_DOUBLE_ARRAY_WITHIN(1e-12, ctx1->C, ctx2->C, ctx1->nmo * ctx1->original_wf->nao);

    ASSERT_STDL_OK(stdl_context_delete(ctx1));
    ASSERT_STDL_OK(stdl_context_delete(ctx2));
}
