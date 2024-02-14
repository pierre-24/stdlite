#include <stdlite/response.h>
#include <stdlite/utils/matrix.h>
#include <stdlite/helpers.h>
#include <stdlite/property.h>
#include <string.h>

#include "tests_suite.h"



void test_property_print_excitations_ok() {
    TEST_IGNORE_MESSAGE("for manual check only");

    /* MUST MATCH:
     * state    eV      nm       fL        Rv(corr)
     *    1    7.354   168.6     0.0105     0.0000    -1.00(   4->   5)  0.04(   4->   9) -0.00(   3->   5)
     *    2    9.190   134.9     0.1303    -0.0000     1.00(   3->   5) -0.03(   3->   9)  0.02(   1->   5)
     *    3   10.022   123.7     0.0000    -0.0000    -1.00(   4->   6) -0.04(   4->  11) -0.03(   4->   7)
     *    4   11.859   104.6     0.2028     0.0000    -1.00(   3->   6) -0.07(   2->   5) -0.03(   3->  11)
     *    5   13.791    89.9     0.3154     0.0000    -1.00(   2->   5)  0.07(   3->   6)  0.02(   3->   7)
     *    6   16.578    74.8     0.2052     0.0000     1.00(   2->   6)  0.02(   2->  11)  0.02(   3->   9)
     */

    stdl_wavefunction * wf = NULL;
    stdl_basis * bs = NULL;
    read_fchk("../tests/test_files/water_631g.fchk", &wf, &bs);

    stdl_context* ctx = NULL;
    ASSERT_STDL_OK(stdl_context_new(wf, bs, 2.0, 4.0, 25. / 27.212, 1e-4, 1.0, &ctx));

    ASSERT_STDL_OK(stdl_context_select_csfs_monopole(ctx, 1));

    // copy A for latter
    float* Ap = malloc(STDL_MATRIX_SP_SIZE(ctx->ncsfs) * sizeof(float));
    TEST_ASSERT_NOT_NULL(Ap);
    memcpy(Ap, ctx->A, STDL_MATRIX_SP_SIZE(ctx->ncsfs) * sizeof(float));

    // compute the dipole moments
    double* dipoles_sp_AO = malloc(3 * STDL_MATRIX_SP_SIZE(wf->nao) * sizeof(double));
    TEST_ASSERT_NOT_NULL(dipoles_sp_AO);

    ASSERT_STDL_OK(stdl_basis_dsp_dipole(bs, dipoles_sp_AO));

    double* dipoles_sp_MO = malloc(3 * STDL_MATRIX_SP_SIZE(wf->nmo) * sizeof(double));
    TEST_ASSERT_NOT_NULL(dipoles_sp_MO);

    for (int cpt = 0; cpt < 3; ++cpt)
        ASSERT_STDL_OK(stdl_wavefunction_dsp_ao_to_dsp_mo(
                wf->nao,
                ctx->nmo,
                ctx->C_orig,
                dipoles_sp_AO + cpt * STDL_MATRIX_SP_SIZE(wf->nao),
                dipoles_sp_MO + cpt * STDL_MATRIX_SP_SIZE(ctx->nmo))
        );

    // request all excitations
    size_t nrequested =  ctx->ncsfs;
    float* etda = malloc(nrequested * sizeof(float));
    float* Xtda = malloc(nrequested * ctx->ncsfs * sizeof(float));

    ASSERT_STDL_OK(stdl_response_TDA_casida(ctx, nrequested, etda, Xtda));

    float* tdipstda = malloc(nrequested * 3 * sizeof(float ));
    stdl_property_transition_dipoles(ctx, nrequested, dipoles_sp_MO, Xtda, NULL, tdipstda);

    ASSERT_STDL_OK(stdl_property_print_excitations(ctx, nrequested, etda, tdipstda));
    ASSERT_STDL_OK(stdl_property_print_excitations_contribs(ctx, nrequested, etda, Xtda, NULL, .0001f));

    // replace A
    free(ctx->A);
    ctx->A = Ap;

    float* etd = malloc(ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(etd);

    float* Xtd = malloc(ctx->ncsfs * ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(Xtd);

    float* Ytd = malloc(ctx->ncsfs * ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(Ytd);

    ASSERT_STDL_OK(stdl_response_TD_casida(ctx, ctx->ncsfs, etd, Xtd, Ytd));

    float* tdipstd = malloc(nrequested * 3 * sizeof(float ));
    stdl_property_transition_dipoles(ctx, nrequested, dipoles_sp_MO, Xtd, Ytd, tdipstd);

    ASSERT_STDL_OK(stdl_property_print_excitations(ctx, nrequested, etd, tdipstd));

    STDL_FREE_ALL(dipoles_sp_MO, dipoles_sp_AO, etda, Xtda, tdipstda, etd, Xtd, Ytd, tdipstd);

    ASSERT_STDL_OK(stdl_context_delete(ctx));
}

void test_polarizability_ok() {
    stdl_wavefunction * wf = NULL;
    stdl_basis * bs = NULL;
    read_fchk("../tests/test_files/water_631g.fchk", &wf, &bs);

    stdl_context* ctx = NULL;
    ASSERT_STDL_OK(stdl_context_new(wf, bs, 2.0, 4.0, 12. / 27.212, 1e-4, 1.0, &ctx));

    TEST_ASSERT_EQUAL_INT(ctx->nmo,7);
    TEST_ASSERT_EQUAL_INT(ctx->nocc, 4);

    ASSERT_STDL_OK(stdl_context_select_csfs_monopole(ctx, 1));

    TEST_ASSERT_EQUAL_INT(ctx->ncsfs, 10);

    // compute dipole integrals and convert to MO
    double* dipoles_sp_AO = malloc(3 * STDL_MATRIX_SP_SIZE(wf->nao) * sizeof(double));
    TEST_ASSERT_NOT_NULL(dipoles_sp_AO);

    double* dipoles_sp_MO = malloc(3 * STDL_MATRIX_SP_SIZE(wf->nmo) * sizeof(double));
    TEST_ASSERT_NOT_NULL(dipoles_sp_MO);

    ASSERT_STDL_OK(stdl_basis_dsp_dipole(bs, dipoles_sp_AO));

    for (int cpt = 0; cpt < 3; ++cpt)
        ASSERT_STDL_OK(stdl_wavefunction_dsp_ao_to_dsp_mo(
                wf->nao,
                ctx->nmo,
                ctx->C_orig,
                dipoles_sp_AO + cpt * STDL_MATRIX_SP_SIZE(wf->nao),
                dipoles_sp_MO + cpt * STDL_MATRIX_SP_SIZE(ctx->nmo))
        );

    // build egrad
    float* egrad = malloc(3 * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(egrad);

    stdl_response_perturbed_gradient(ctx, 3, dipoles_sp_MO, egrad);

    // solve response
    float* X = malloc(2 * 3 * ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(X);

    float* Y = malloc(2 * 3 * ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(Y);

    ASSERT_STDL_OK(stdl_response_TD_linear(ctx, 2, (float[]) {0, 4.2822696E-02}, 3, egrad, X, Y));

    // compute polarizabilities
    float alpha_static[6], alpha_dynamic[6], static_mean, dynamic_mean;
    stdl_property_polarizability(ctx, dipoles_sp_MO, X, Y, alpha_static);
    stdl_property_mean_polarizability(alpha_static, &static_mean);
    TEST_ASSERT_DOUBLE_WITHIN(1e-2, 4.06, static_mean);

    stdl_property_polarizability(ctx, dipoles_sp_MO, X + 3 * ctx->ncsfs, Y + 3 * ctx->ncsfs, alpha_dynamic);
    stdl_property_mean_polarizability(alpha_dynamic, &dynamic_mean);
    TEST_ASSERT_DOUBLE_WITHIN(1e-2, 4.102, dynamic_mean);

    // stdl_matrix_ssp_print(3, alpha_dynamic, "a(-w;w)");

    STDL_FREE_ALL(dipoles_sp_AO, dipoles_sp_MO, egrad, X, Y);

    ASSERT_STDL_OK(stdl_context_delete(ctx));
}


void test_polarizability_TDA_ok() {
    stdl_wavefunction * wf = NULL;
    stdl_basis * bs = NULL;
    read_fchk("../tests/test_files/water_631g.fchk", &wf, &bs);

    stdl_context* ctx = NULL;
    ASSERT_STDL_OK(stdl_context_new(wf, bs, 2.0, 4.0, 12. / 27.212, 1e-4, 1.0, &ctx));

    TEST_ASSERT_EQUAL_INT(ctx->nmo,7);
    TEST_ASSERT_EQUAL_INT(ctx->nocc, 4);

    ASSERT_STDL_OK(stdl_context_select_csfs_monopole(ctx, 1));

    TEST_ASSERT_EQUAL_INT(ctx->ncsfs, 10);

    // compute dipole integrals and convert to MO
    double* dipoles_sp_AO = malloc(3 * STDL_MATRIX_SP_SIZE(wf->nao) * sizeof(double));
    TEST_ASSERT_NOT_NULL(dipoles_sp_AO);

    double* dipoles_sp_MO = malloc(3 * STDL_MATRIX_SP_SIZE(wf->nmo) * sizeof(double));
    TEST_ASSERT_NOT_NULL(dipoles_sp_MO);

    ASSERT_STDL_OK(stdl_basis_dsp_dipole(bs, dipoles_sp_AO));

    for (int cpt = 0; cpt < 3; ++cpt)
        ASSERT_STDL_OK(stdl_wavefunction_dsp_ao_to_dsp_mo(
                wf->nao,
                ctx->nmo,
                ctx->C_orig,
                dipoles_sp_AO + cpt * STDL_MATRIX_SP_SIZE(wf->nao),
                dipoles_sp_MO + cpt * STDL_MATRIX_SP_SIZE(ctx->nmo))
        );

    // build egrad
    float* egrad = malloc(3 * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(egrad);

    stdl_response_perturbed_gradient(ctx, 3, dipoles_sp_MO, egrad);

    // copy A for latter
    float* Ap = malloc(STDL_MATRIX_SP_SIZE(ctx->ncsfs) * sizeof(float));
    TEST_ASSERT_NOT_NULL(Ap);
    memcpy(Ap, ctx->A, STDL_MATRIX_SP_SIZE(ctx->ncsfs) * sizeof(float));

    // solve response
    float* X = malloc(2 * 3 * ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(X);

    float* Y = malloc(2 * 3 * ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(Y);

    ASSERT_STDL_OK(stdl_response_TD_linear(ctx, 2, (float[]) {0, 4.2822696E-02 * 2}, 3, egrad, X, Y));

    // compute polarizability
    float alpha_static[6], alpha_dynamic[6], static_mean, dynamic_mean;

    stdl_property_polarizability(ctx, dipoles_sp_MO, X, Y, alpha_static);
    stdl_property_mean_polarizability(alpha_static, &static_mean);
    stdl_property_polarizability(ctx, dipoles_sp_MO, X + 3 * ctx->ncsfs, Y + 3 * ctx->ncsfs, alpha_dynamic);
    stdl_property_mean_polarizability(alpha_dynamic, &dynamic_mean);

    // stdl_matrix_ssp_print(3, alpha_dynamic, "a(-w;w)");

    // replace A
    free(ctx->A);
    ctx->A = Ap;

    // solve
    float* Xtda = malloc(2 * 3 * ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(Xtda);

    ASSERT_STDL_OK(stdl_response_TDA_linear(ctx, 2, (float[]) {0, 4.2822696E-02 * 2}, 3, egrad, Xtda));

    // compute polarizability
    float alpha_static_tda[6], alpha_dynamic_tda[6], static_mean_tda, dynamic_mean_tda;
    stdl_property_polarizability(ctx, dipoles_sp_MO, Xtda, NULL, alpha_static_tda);
    stdl_property_mean_polarizability(alpha_static_tda, &static_mean_tda);
    stdl_property_polarizability(ctx, dipoles_sp_MO, Xtda + 3 * ctx->ncsfs, NULL, alpha_dynamic_tda);
    stdl_property_mean_polarizability(alpha_dynamic_tda, &dynamic_mean_tda);

    //printf("static: %f → %f\n", static_mean, static_mean_tda);
    //printf("dynamic: %f → %f\n", dynamic_mean, dynamic_mean_tda);

    // stdl_matrix_ssp_print(3, alpha_dynamic_tda, "a(-w;w)");

    // TEST_ASSERT_FLOAT_WITHIN(1e-1, dynamic_mean / static_mean, dynamic_mean_tda / static_mean_tda);

    STDL_FREE_ALL(dipoles_sp_AO, dipoles_sp_MO, egrad, X, Y, Xtda);

    ASSERT_STDL_OK(stdl_context_delete(ctx));
}
