#include <stdlite/response.h>
#include <stdlite/utils/matrix.h>
#include <stdlite/helpers.h>
#include <stdlite/property.h>
#include <stdlite/utils/experimental_quantity.h>
#include <stdlite/utils/permutations.h>

#include <string.h>

#include "tests_suite.h"

void setUp() {
    stdl_set_debug_level(-1);
    stdl_set_log_level(0);
}

void test_property_polarizability_TD_ok() {
    stdl_wavefunction * wf = NULL;
    stdl_basis * bs = NULL;
    read_fchk("../tests/test_files/water_631g.fchk", &wf, &bs);

    stdl_context* ctx = NULL;
    ASSERT_STDL_OK(stdl_context_new(wf, bs, 2.0, 4.0, 12. / STDL_CONST_AU_TO_EV, 1e-4, 1.0, &ctx));
    ASSERT_STDL_OK(stdl_context_select_csfs_monopole(ctx, 1));

    // compute dipole integrals and convert to MO
    double* dipoles_mat = malloc(3 * STDL_MATRIX_SP_SIZE(ctx->nmo) * sizeof(double));
    TEST_ASSERT_NOT_NULL(dipoles_mat);

    make_dipoles_MO(wf, bs, ctx, dipoles_mat);

    // build egrad
    float* egrad = malloc(3 * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(egrad);

    stdl_response_perturbed_gradient(ctx, 3, dipoles_mat, egrad);

    // solve response
    size_t nw = 3;
    float w[] = {0, 4.282270E-2f, 4.282270E-2f * 2};

    float* X = malloc(nw * 3 * ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(X);

    float* Y = malloc(nw * 3 * ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(Y);

    ASSERT_STDL_OK(stdl_response_TD_linear(ctx, nw, w, 3, egrad, X, Y));

    // compute polarizabilities
    float alpha[6], alpha_iso, alpha_aniso, result[] = {4.061f, 4.103f, 4.236f};

    for (size_t iw = 0; iw < nw; ++iw) {
        stdl_property_polarizability(ctx, dipoles_mat, X + iw * 3 * ctx->ncsfs, Y + iw * 3 * ctx->ncsfs, alpha);
        stdl_qexp_polarizability(alpha, &alpha_iso, &alpha_aniso);
        TEST_ASSERT_FLOAT_WITHIN(1e-2, result[iw], alpha_iso);
    }

    STDL_FREE_ALL(dipoles_mat, egrad, X, Y);

    ASSERT_STDL_OK(stdl_context_delete(ctx));
}

void test_property_polarizability_TD_fchk_vs_molden_ok() {
    // FCHK
    stdl_wavefunction * wf_fchk = NULL;
    stdl_basis * bs_fchk = NULL;
    read_fchk("../tests/test_files/water_631gdf.fchk", &wf_fchk, &bs_fchk);

    stdl_context* ctx_fchk = NULL;
    ASSERT_STDL_OK(stdl_context_new(wf_fchk, bs_fchk, 2.0, 4.0, 12. / STDL_CONST_AU_TO_EV, 1e-4, 1.0, &ctx_fchk));
    ASSERT_STDL_OK(stdl_context_select_csfs_monopole(ctx_fchk, 1));

    double* dipoles_mat_fchk = malloc(3 * STDL_MATRIX_SP_SIZE(ctx_fchk->nmo) * sizeof(double));
    TEST_ASSERT_NOT_NULL(dipoles_mat_fchk);

    make_dipoles_MO(wf_fchk, bs_fchk, ctx_fchk, dipoles_mat_fchk);

    float* egrad_fchk = malloc(3 * ctx_fchk->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(egrad_fchk);

    stdl_response_perturbed_gradient(ctx_fchk, 3, dipoles_mat_fchk, egrad_fchk);

    size_t nw = 3;
    float w[] = {0, 4.282270E-2f, 2 * 4.282270E-2f};

    float* X_fchk = malloc(nw * 3 * ctx_fchk->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(X_fchk);

    float* Y_fchk = malloc(nw * 3 * ctx_fchk->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(Y_fchk);

    ASSERT_STDL_OK(stdl_response_TD_linear(ctx_fchk, nw, w, 3, egrad_fchk, X_fchk, Y_fchk));

    // molden
    stdl_wavefunction * wf_molden = NULL;
    stdl_basis * bs_molden = NULL;
    read_molden("../tests/test_files/water_631gdf.molden", &wf_molden, &bs_molden);

    stdl_context* ctx_molden = NULL;
    ASSERT_STDL_OK(stdl_context_new(wf_molden, bs_molden, 2.0, 4.0, 12. / STDL_CONST_AU_TO_EV, 1e-4, 1.0, &ctx_molden));
    ASSERT_STDL_OK(stdl_context_select_csfs_monopole(ctx_molden, 1));

    double* dipoles_mat_molden = malloc(3 * STDL_MATRIX_SP_SIZE(ctx_molden->nmo) * sizeof(double));
    TEST_ASSERT_NOT_NULL(dipoles_mat_molden);

    make_dipoles_MO(wf_molden, bs_molden, ctx_molden, dipoles_mat_molden);

    float* egrad_molden = malloc(3 * ctx_molden->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(egrad_molden);

    stdl_response_perturbed_gradient(ctx_molden, 3, dipoles_mat_molden, egrad_molden);

    float* X_molden = malloc(nw * 3 * ctx_molden->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(X_molden);

    float* Y_molden = malloc(nw * 3 * ctx_molden->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(Y_molden);

    ASSERT_STDL_OK(stdl_response_TD_linear(ctx_molden, nw, w, 3, egrad_molden, X_molden, Y_molden));

    // compute polarizabilities
    float alpha[6], alpha_iso_fchk, alpha_aniso_fchk, alpha_iso_molden, alpha_aniso_molden;

    for (size_t iw = 0; iw < nw; ++iw) {
        stdl_property_polarizability(ctx_fchk, dipoles_mat_fchk, X_fchk + iw * 3 * ctx_fchk->ncsfs, Y_fchk + iw * 3 * ctx_fchk->ncsfs, alpha);
        stdl_qexp_polarizability(alpha, &alpha_iso_fchk, &alpha_aniso_fchk);

        stdl_property_polarizability(ctx_molden, dipoles_mat_molden, X_molden + iw * 3 * ctx_molden->ncsfs, Y_molden + iw * 3 * ctx_molden->ncsfs, alpha);
        stdl_qexp_polarizability(alpha, &alpha_iso_molden, &alpha_aniso_molden);

        TEST_ASSERT_FLOAT_WITHIN(1e-2, alpha_iso_molden, alpha_iso_fchk);
        TEST_ASSERT_FLOAT_WITHIN(1e-2, alpha_aniso_molden, alpha_aniso_fchk);
    }

    STDL_FREE_ALL(dipoles_mat_fchk, egrad_fchk, X_fchk, Y_fchk);
    ASSERT_STDL_OK(stdl_context_delete(ctx_fchk));

    STDL_FREE_ALL(dipoles_mat_molden, egrad_molden, X_molden, Y_molden);
    ASSERT_STDL_OK(stdl_context_delete(ctx_molden));
}

void test_property_polarizability_TD_SOS_ok() {
    stdl_wavefunction * wf = NULL;
    stdl_basis * bs = NULL;
    read_molden("../tests/test_files/water_631gdf_sph.molden", &wf, &bs);

    stdl_context* ctx = NULL;
    ASSERT_STDL_OK(stdl_context_new(wf, bs, 2.0, 4.0, 25. / STDL_CONST_AU_TO_EV, 1e-4, 1.0, &ctx));
    ASSERT_STDL_OK(stdl_context_select_csfs_monopole(ctx, 1));

    // fetch excitations
    float* etd = malloc(ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(etd);

    float* Xamptd = malloc(ctx->ncsfs * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(Xamptd);

    float* Yamptd = malloc(ctx->ncsfs * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(Yamptd);

    ASSERT_STDL_OK(stdl_response_TD_casida(ctx, ctx->ncsfs, etd, Xamptd, Yamptd));

    // compute dipole integrals and convert to MO
    double* dipoles_mat = malloc(3 * STDL_MATRIX_SP_SIZE(wf->nmo) * sizeof(double));
    TEST_ASSERT_NOT_NULL(dipoles_mat);

    make_dipoles_MO(wf, bs, ctx, dipoles_mat);

    // get transition dipoles
    float* tdipstd = malloc(ctx->ncsfs * 3 * sizeof(float ));
    stdl_property_transition_dipoles(ctx, ctx->ncsfs, dipoles_mat, Xamptd, Yamptd, tdipstd);

    // build egrad
    float* egrad = malloc(3 * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(egrad);

    stdl_response_perturbed_gradient(ctx, 3, dipoles_mat, egrad);

    // solve linear response
    size_t nw = 3;
    float w[] = {0, STDL_CONST_HC / 1064.f, STDL_CONST_HC / 532.f};

    float* Xtd = malloc(nw * 3 * ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(Xtd);

    float* Ytd = malloc(nw * 3 * ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(Ytd);

    ASSERT_STDL_OK(stdl_response_TD_linear(ctx, nw, w, 3, egrad, Xtd, Ytd));

    // compute polarizabilities
    float alpha[6], alpha_zz;

    for (size_t iw = 0; iw < nw; ++iw) {
        stdl_property_polarizability(ctx, dipoles_mat, Xtd + iw * 3 * ctx->ncsfs, Ytd + iw * 3 * ctx->ncsfs, alpha);

        // stdl_matrix_ssp_print(3, alpha, "alpha");

        // compute alpha_zz from SOS
        alpha_zz = .0f;
        for (size_t iexci = 0; iexci < ctx->ncsfs; ++iexci)
            alpha_zz += powf(tdipstd[2 * ctx->ncsfs + iexci], 2) / (etd[iexci] - w[iw]) + powf(tdipstd[2 * ctx->ncsfs + iexci], 2) / (etd[iexci] + w[iw]);

        TEST_ASSERT_FLOAT_WITHIN(1e-2, alpha[STDL_MATRIX_SP_IDX(2, 2)], alpha_zz);
    }

    STDL_FREE_ALL(etd, Xamptd, Yamptd, dipoles_mat, tdipstd, egrad, Xtd, Ytd);
    ASSERT_STDL_OK(stdl_context_delete(ctx));
}


void test_property_first_hyperpolarizability_TD_ok() {
    stdl_wavefunction * wf = NULL;
    stdl_basis * bs = NULL;
    read_fchk("../tests/test_files/water_631g.fchk", &wf, &bs);

    stdl_context* ctx = NULL;
    ASSERT_STDL_OK(stdl_context_new(wf, bs, 2.0, 4.0, 12. / STDL_CONST_AU_TO_EV, 1e-4, 1.0, &ctx));
    ASSERT_STDL_OK(stdl_context_select_csfs_monopole(ctx, 1));

    // compute dipole integrals and convert to MO
    double* dipoles_mat = malloc(3 * STDL_MATRIX_SP_SIZE(ctx->nmo) * sizeof(double));
    TEST_ASSERT_NOT_NULL(dipoles_mat);

    make_dipoles_MO(wf, bs, ctx, dipoles_mat);

    // build egrad
    float* egrad = malloc(3 * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(egrad);

    stdl_response_perturbed_gradient(ctx, 3, dipoles_mat, egrad);

    // solve response
    size_t nw = 4;
    float w[] = {0, STDL_CONST_HC / 1064.f, -STDL_CONST_HC / 532.f, STDL_CONST_HC / 532.f};

    float* X = malloc(nw * 3 * ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(X);

    float* Y = malloc(nw * 3 * ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(Y);

    ASSERT_STDL_OK(stdl_response_TD_linear(ctx, nw, w, 3, egrad, X, Y));

    // compute beta
    float beta[3][3][3], beta2_ZZZ, beta2_ZXX;
    ASSERT_STDL_OK(stdl_property_first_hyperpolarizability(
            ctx,
            dipoles_mat,
            (float* []) {X, X, X},
            (float* []) {Y, Y, Y},
            (float *) beta
            ));

    // stdl_matrix_sge_print(9, 3, beta, "beta");
    stdl_qexp_first_hyperpolarizability_hrs(beta, &beta2_ZZZ, &beta2_ZXX);

    TEST_ASSERT_FLOAT_WITHIN(1e-2, 234.07, beta2_ZZZ);
    TEST_ASSERT_FLOAT_WITHIN(1e-2, 38.16, beta2_ZXX);

    ASSERT_STDL_OK(stdl_property_first_hyperpolarizability(
            ctx,
            dipoles_mat,
            (float* []) {X + 2 * 3 * ctx->ncsfs, X + 1 * 3 * ctx->ncsfs, X + 1 * 3 * ctx->ncsfs},
            (float* []) {Y + 2 * 3 * ctx->ncsfs, Y + 1 * 3 * ctx->ncsfs, Y + 1 * 3 * ctx->ncsfs},
            (float*) beta
            ));

    // stdl_matrix_sge_print(9, 3, beta, "beta");
    stdl_qexp_first_hyperpolarizability_hrs(beta, &beta2_ZZZ, &beta2_ZXX);

    TEST_ASSERT_FLOAT_WITHIN(1e-2, 272.19, beta2_ZZZ);
    TEST_ASSERT_FLOAT_WITHIN(1e-2, 43.85, beta2_ZXX);

    // check that one gets the same result by using X(-w) = Y(w)
    ASSERT_STDL_OK(stdl_property_first_hyperpolarizability(
            ctx,
            dipoles_mat,
            (float* []) {Y + 3 * 3 * ctx->ncsfs, X + 1 * 3 * ctx->ncsfs, X + 1 * 3 * ctx->ncsfs},
            (float* []) {X + 3 * 3 * ctx->ncsfs, Y + 1 * 3 * ctx->ncsfs, Y + 1 * 3 * ctx->ncsfs},
            (float*) beta
    ));

    stdl_qexp_first_hyperpolarizability_hrs(beta, &beta2_ZZZ, &beta2_ZXX);

    TEST_ASSERT_FLOAT_WITHIN(1e-2, 272.19, beta2_ZZZ);
    TEST_ASSERT_FLOAT_WITHIN(1e-2, 43.85, beta2_ZXX);

    STDL_FREE_ALL(dipoles_mat, egrad, X, Y);

    ASSERT_STDL_OK(stdl_context_delete(ctx));
}

void test_property_polarizability_TDA_ok() {
    stdl_wavefunction * wf = NULL;
    stdl_basis * bs = NULL;
    read_fchk("../tests/test_files/water_631g.fchk", &wf, &bs);

    stdl_context* ctx = NULL;
    ASSERT_STDL_OK(stdl_context_new(wf, bs, 2.0, 4.0, 12. / STDL_CONST_AU_TO_EV, 1e-4, 1.0, &ctx));
    ASSERT_STDL_OK(stdl_context_select_csfs_monopole(ctx, 1));

    // compute dipole integrals and convert to MO
    double* dipoles_mat = malloc(3 * STDL_MATRIX_SP_SIZE(ctx->nmo) * sizeof(double));
    TEST_ASSERT_NOT_NULL(dipoles_mat);

    make_dipoles_MO(wf, bs, ctx, dipoles_mat);

    // build egrad
    float* egrad = malloc(3 * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(egrad);

    stdl_response_perturbed_gradient(ctx, 3, dipoles_mat, egrad);

    // solve response
    size_t nw = 3;
    float w[] = {0, STDL_CONST_HC / 1064.f, STDL_CONST_HC / 532.f};

    // solve
    float* Xtda = malloc(nw * 3 * ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(Xtda);

    float* Ytda = malloc(nw * 3 * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(Ytda);

    // compute polarizabilities
    float alpha[6], alpha_iso, alpha_aniso, result[] = {4.128f, 4.170f, 4.303f};
    ASSERT_STDL_OK(stdl_response_TDA_linear(ctx, nw, w, 3, egrad, Xtda, Ytda));

    for (size_t iw = 0; iw < nw; ++iw) {
        stdl_property_polarizability(ctx, dipoles_mat, Xtda + iw * 3 * ctx->ncsfs, Ytda + iw * 3 * ctx->ncsfs, alpha);
        stdl_qexp_polarizability(alpha, &alpha_iso, &alpha_aniso);
        TEST_ASSERT_FLOAT_WITHIN(1e-2, result[iw], alpha_iso);
    }

    STDL_FREE_ALL(dipoles_mat, egrad, Xtda, Ytda);

    ASSERT_STDL_OK(stdl_context_delete(ctx));
}

void test_property_polarizability_TDA_SOS_ok() {
    stdl_wavefunction * wf = NULL;
    stdl_basis * bs = NULL;
    read_fchk("../tests/test_files/water_631gdf.fchk", &wf, &bs);

    stdl_context* ctx = NULL;
    ASSERT_STDL_OK(stdl_context_new(wf, bs, 2.0, 4.0, 25. / STDL_CONST_AU_TO_EV, 1e-4, 1.0, &ctx));

    ASSERT_STDL_OK(stdl_context_select_csfs_monopole(ctx, 0));

    // compute dipole integrals and convert to MO
    double* dipoles_mat = malloc(3 * STDL_MATRIX_SP_SIZE(wf->nmo) * sizeof(double));
    TEST_ASSERT_NOT_NULL(dipoles_mat);

    make_dipoles_MO(wf, bs, ctx, dipoles_mat);

    // build egrad
    float* egrad = malloc(3 * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(egrad);

    stdl_response_perturbed_gradient(ctx, 3, dipoles_mat, egrad);

    // fetch all excitations
    float* etda = malloc(ctx->ncsfs * sizeof(float ));
    float* Xamptda = malloc(ctx->ncsfs * ctx->ncsfs * sizeof(float ));
    ASSERT_STDL_OK(stdl_response_TDA_casida(ctx, ctx->ncsfs, etda, Xamptda));

    // get transition dipoles
    float* tdipstda = malloc(ctx->ncsfs * 3 * sizeof(float ));
    ASSERT_STDL_OK(stdl_property_transition_dipoles(ctx, ctx->ncsfs, dipoles_mat, Xamptda, NULL, tdipstda));

    // solve linear response
    size_t nw = 3;
    float w[] = {0, STDL_CONST_HC / 1064.f, STDL_CONST_HC / 532.f};

    float* Xtda = malloc(nw * 3 * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(Xtda);

    float* Ytda = malloc(nw * 3 * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(Ytda);

    ASSERT_STDL_OK(stdl_response_TDA_linear(ctx, nw, w, 3, egrad, Xtda, Ytda));

    // compute polarizabilities
    float alpha[6], alpha_zz;

    for (size_t iw = 0; iw < nw; ++iw) {
        stdl_property_polarizability(ctx, dipoles_mat, Xtda + iw * 3 * ctx->ncsfs, Ytda + iw * 3 * ctx->ncsfs, alpha);

        // stdl_matrix_ssp_print(3, alpha, "alpha");

        // compute alpha_zz from SOS
        alpha_zz = .0f;
        for (size_t iexci = 0; iexci < ctx->ncsfs; ++iexci)
            alpha_zz += powf(tdipstda[2 * ctx->ncsfs + iexci], 2) / (etda[iexci] - w[iw]) + powf(tdipstda[2 * ctx->ncsfs + iexci], 2) / (w[iw] + etda[iexci]);

        TEST_ASSERT_FLOAT_WITHIN(1e-2, alpha[STDL_MATRIX_SP_IDX(2, 2)], alpha_zz);
    }

    STDL_FREE_ALL(dipoles_mat, etda, Xamptda, egrad, Xtda, Ytda, tdipstda);

    ASSERT_STDL_OK(stdl_context_delete(ctx));
}

float oscillator_strength(size_t i, size_t j, size_t nexci, float* e, float * e2etdips) {
    float x = .0f;
    for (int cpt = 0; cpt < 3; ++cpt) {
        x += powf(e2etdips[cpt * STDL_MATRIX_SP_SIZE(nexci) + STDL_MATRIX_SP_IDX(i, j)], 2);
    }

    return 2.f/3 * (e[j] - e[i]) * x;
}

void test_property_e2e_transition_dipoles_ok() {
    stdl_wavefunction * wf = NULL;
    stdl_basis * bs = NULL;
    read_fchk("../tests/test_files/water_sto3g.fchk", &wf, &bs);

    stdl_context* ctx = NULL;
    ASSERT_STDL_OK(stdl_context_new(wf, bs, 2.0, 4.0, 20. / STDL_CONST_AU_TO_EV, 1e-4, 1.0, &ctx));

    ASSERT_STDL_OK(stdl_context_select_csfs_monopole(ctx, 1));

    // compute dipole integrals and convert to MO
    double* dipoles_sp_MO = malloc(3 * STDL_MATRIX_SP_SIZE(ctx->nmo) * sizeof(double));
    TEST_ASSERT_NOT_NULL(dipoles_sp_MO);

    make_dipoles_MO(wf, bs, ctx, dipoles_sp_MO);

    // request all excitations
    size_t nrequested = 4;

    float* etd = malloc(ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(etd);

    float* Xtd = malloc(ctx->ncsfs * ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(Xtd);

    float* Ytd = malloc(ctx->ncsfs * ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(Ytd);

    ASSERT_STDL_OK(stdl_response_TD_casida(ctx, ctx->ncsfs, etd, Xtd, Ytd));

    float* e2etdips = malloc(3 * STDL_MATRIX_SP_SIZE(nrequested) * sizeof(float ));
    TEST_ASSERT_NOT_NULL(e2etdips);

    stdl_property_e2e_transition_dipoles(ctx, nrequested, dipoles_sp_MO, Xtd, Ytd, e2etdips);

    // checked against stda
    TEST_ASSERT_FLOAT_WITHIN(1e-4, .0005f, oscillator_strength(0, 1, nrequested, etd, e2etdips));
    TEST_ASSERT_FLOAT_WITHIN(1e-4, 0.0235f, oscillator_strength(0, 2, nrequested, etd, e2etdips));
    TEST_ASSERT_FLOAT_WITHIN(1e-4, .0f, oscillator_strength(0, 3, nrequested, etd, e2etdips));

    STDL_FREE_ALL(dipoles_sp_MO, etd, Xtd, Ytd, e2etdips);

    ASSERT_STDL_OK(stdl_context_delete(ctx));
}

void test_property_first_hyperpolarizability_TD_SOS_ok() {
    stdl_wavefunction * wf = NULL;
    stdl_basis * bs = NULL;
    read_fchk("../tests/test_files/water_sto3g.fchk", &wf, &bs);

    stdl_context* ctx = NULL;
    ASSERT_STDL_OK(stdl_context_new(wf, bs, 2.0, 4.0, 20. / STDL_CONST_AU_TO_EV, 1e-4, 1.0, &ctx));
    ASSERT_STDL_OK(stdl_context_select_csfs_monopole(ctx, 1));

    // fetch excitations
    float* etd = malloc(ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(etd);

    float* Xamptd = malloc(ctx->ncsfs * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(Xamptd);

    float* Yamptd = malloc(ctx->ncsfs * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(Yamptd);

    ASSERT_STDL_OK(stdl_response_TD_casida(ctx, ctx->ncsfs, etd, Xamptd, Yamptd));

    // compute dipole integrals and convert to MO
    double* dipoles_mat = malloc(3 * STDL_MATRIX_SP_SIZE(wf->nmo) * sizeof(double));
    TEST_ASSERT_NOT_NULL(dipoles_mat);

    make_dipoles_MO(wf, bs, ctx, dipoles_mat);

    // get 0→m transition dipoles
    float* t0mdipstd = malloc(ctx->ncsfs * 3 * sizeof(float ));
    stdl_property_transition_dipoles(ctx, ctx->ncsfs, dipoles_mat, Xamptd, Yamptd, t0mdipstd);

    // get m→n transition dipoles
    float* tmndipstd = malloc(3 * STDL_MATRIX_SP_SIZE(ctx->ncsfs) * sizeof(float ));
    TEST_ASSERT_NOT_NULL(tmndipstd);

    stdl_property_e2e_transition_dipoles(ctx, ctx->ncsfs , dipoles_mat, Xamptd, Yamptd, tmndipstd);

    // build egrad
    float* egrad = malloc(3 * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(egrad);

    stdl_response_perturbed_gradient(ctx, 3, dipoles_mat, egrad);

    // solve linear response
    size_t nw = 3;
    float w[] = {0, STDL_CONST_HC / 1064.f, -STDL_CONST_HC / 532.f};

    float* Xtd = malloc(nw * 3 * ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(Xtd);

    float* Ytd = malloc(nw * 3 * ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(Ytd);

    ASSERT_STDL_OK(stdl_response_TD_linear(ctx, nw, w, 3, egrad, Xtd, Ytd));

    // compute hyperpolarizabilities
    float beta[3][3][3], beta_component;
    int to_permute[6];

    to_permute[1] = to_permute[3] =  to_permute[5] = 2; // zzz

    for (size_t i = 0; i < 2; ++i) {
        if(i == 0) {
            to_permute[0] = to_permute[2] = to_permute[4] = 0; // static
        } else {
            to_permute[0] = 2; // -2w
            to_permute[2] = to_permute[4] = 1; // w
        }

        ASSERT_STDL_OK(stdl_property_first_hyperpolarizability(
                ctx,
                dipoles_mat,
                (float *[]) {Xtd + to_permute[0] * 3 * ctx->ncsfs, Xtd + to_permute[2] * 3 * ctx->ncsfs, Xtd + to_permute[4] * 3 * ctx->ncsfs},
                (float *[]) {Ytd + to_permute[0] * 3 * ctx->ncsfs, Ytd + to_permute[2] * 3 * ctx->ncsfs, Ytd + to_permute[4] * 3 * ctx->ncsfs},
                beta
        ));

        // stdl_matrix_sge_print(9, 3, beta, "beta");

        // compare with SOS
        beta_component = 0;
        stdl_permutations* perms = NULL;
        stdl_permutations_new(to_permute, 3, 2 * sizeof(int ), &perms);
        stdl_permutations_remove_duplicates(perms, 3, 2 * sizeof(int));

        stdl_permutations* current = perms;
        size_t nperms = 0;
        while (current != NULL) {
            int* r = (int*) current->perm;
            int zeta = r[1], sigma = r[3], tau = r[5], iw = r[0], kw = r[4];

            for (size_t m = 0; m < ctx->ncsfs; ++m) {
                for (size_t n = 0; n < ctx->ncsfs; ++n) {
                    beta_component += t0mdipstd[zeta * ctx->ncsfs + m] * tmndipstd[sigma * STDL_MATRIX_SP_SIZE(ctx->ncsfs) + STDL_MATRIX_SP_IDX(m, n)] * t0mdipstd[tau * ctx->ncsfs + n] / ((etd[m] + w[iw]) * (etd[n] - w[kw]));
                }
            }

            nperms++;
            current = current->next;
        }

        // printf("beta_component = %f\n", 6 / (float) nperms * beta_component);
        TEST_ASSERT_FLOAT_WITHIN(1e-1, beta[to_permute[1]][to_permute[3]][to_permute[5]], 6 / (float) nperms * beta_component);

        stdl_permutations_delete(perms);
    }

    STDL_FREE_ALL(etd, Xamptd, Yamptd, dipoles_mat, t0mdipstd, tmndipstd, egrad, Xtd, Ytd);
    ASSERT_STDL_OK(stdl_context_delete(ctx));
}

void test_property_first_hyperpolarizability_TDA_SOS_ok() {
    stdl_wavefunction * wf = NULL;
    stdl_basis * bs = NULL;
    read_fchk("../tests/test_files/water_631gdf_sph.fchk", &wf, &bs);

    stdl_context* ctx = NULL;
    ASSERT_STDL_OK(stdl_context_new(wf, bs, 2.0, 4.0, 20. / STDL_CONST_AU_TO_EV, 1e-4, 1.0, &ctx));
    ASSERT_STDL_OK(stdl_context_select_csfs_monopole(ctx, 1));

    // fetch excitations
    float* etda = malloc(ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(etda);

    float* Xamptda = malloc(ctx->ncsfs * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(Xamptda);

    ASSERT_STDL_OK(stdl_response_TDA_casida(ctx, ctx->ncsfs, etda, Xamptda));

    // compute dipole integrals and convert to MO
    double* dipoles_mat = malloc(3 * STDL_MATRIX_SP_SIZE(wf->nmo) * sizeof(double));
    TEST_ASSERT_NOT_NULL(dipoles_mat);

    make_dipoles_MO(wf, bs, ctx, dipoles_mat);

    // get 0→m transition dipoles
    float* t0mdipstda = malloc(ctx->ncsfs * 3 * sizeof(float ));
    stdl_property_transition_dipoles(ctx, ctx->ncsfs, dipoles_mat, Xamptda, NULL, t0mdipstda);

    // get m→n transition dipoles
    float* tmndipstda = malloc(3 * STDL_MATRIX_SP_SIZE(ctx->ncsfs) * sizeof(float ));
    TEST_ASSERT_NOT_NULL(tmndipstda);

    stdl_property_e2e_transition_dipoles(ctx, ctx->ncsfs , dipoles_mat, Xamptda, NULL, tmndipstda);

    // build egrad
    float* egrad = malloc(3 * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(egrad);

    stdl_response_perturbed_gradient(ctx, 3, dipoles_mat, egrad);

    // solve linear response
    size_t nw = 3;
    float w[] = {0, STDL_CONST_HC / 1064.f, - STDL_CONST_HC / 532.f};

    float* Xtda = malloc(nw * 3 * ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(Xtda);

    float* Ytda = malloc(nw * 3 * ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(Ytda);

    ASSERT_STDL_OK(stdl_response_TDA_linear(ctx, nw, w, 3, egrad, Xtda, Ytda));

    // compute hyperpolarizabilities
    float beta[3][3][3], beta_component;
    int to_permute[6];

    to_permute[1] = to_permute[3] =  to_permute[5] = 2; // zzz

    for (size_t i = 0; i < 2; ++i) {
        if(i == 0) {
            to_permute[0] = to_permute[2] = to_permute[4] = 0; // static
        } else {
            to_permute[0] = 2; // -2w
            to_permute[2] = to_permute[4] = 1; // w
        }

        ASSERT_STDL_OK(stdl_property_first_hyperpolarizability(
                ctx,
                dipoles_mat,
                (float *[]) {Xtda + to_permute[0] * 3 * ctx->ncsfs, Xtda + to_permute[2] * 3 * ctx->ncsfs, Xtda + to_permute[4] * 3 * ctx->ncsfs},
                (float *[]) {Ytda + to_permute[0] * 3 * ctx->ncsfs, Ytda + to_permute[2] * 3 * ctx->ncsfs, Ytda + to_permute[4] * 3 * ctx->ncsfs},
                beta
        ));

        // stdl_matrix_sge_print(9, 3, beta, "beta");

        // compare with SOS
        beta_component = 0;
        stdl_permutations* perms = NULL;
        stdl_permutations_new(to_permute, 3, 2 * sizeof(int ), &perms);
        stdl_permutations_remove_duplicates(perms, 3, 2 * sizeof(int));

        stdl_permutations* current = perms;
        size_t nperms = 0;
        while (current != NULL) {
            int* r = (int*) current->perm;
            int zeta = r[1], sigma = r[3], tau = r[5], iw = r[0], kw = r[4];

            for (size_t m = 0; m < ctx->ncsfs; ++m) {
                for (size_t n = 0; n < ctx->ncsfs; ++n) {
                    beta_component += t0mdipstda[zeta * ctx->ncsfs + m] * tmndipstda[sigma * STDL_MATRIX_SP_SIZE(ctx->ncsfs) + STDL_MATRIX_SP_IDX(m, n)] * t0mdipstda[tau * ctx->ncsfs + n] / ((etda[m] + w[iw]) * (etda[n] - w[kw]));
                }
            }

            nperms++;
            current = current->next;
        }

        // printf("beta_component = %f\n", 6 / (float) nperms * beta_component);
        TEST_ASSERT_FLOAT_WITHIN(1e-1, beta[to_permute[1]][to_permute[3]][to_permute[5]], 6 / (float) nperms * beta_component);

        stdl_permutations_delete(perms);
    }

    STDL_FREE_ALL(etda, Xamptda, dipoles_mat, t0mdipstda, tmndipstda, egrad, Xtda, Ytda);
    ASSERT_STDL_OK(stdl_context_delete(ctx));
}