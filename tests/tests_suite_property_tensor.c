#include <stdlite/response.h>
#include <stdlite/utils/matrix.h>
#include <stdlite/helpers.h>
#include <stdlite/utils/experimental_quantity.h>
#include <stdlite/property.h>
#include <stdlite/property_tensor.h>
#include <stdlite/integrals.h>

#include <string.h>

#include "tests_suite.h"

void setUp() {
    stdl_set_debug_level(-1);
    stdl_set_log_level(0);
}

void test_response_polarizability_TD_ok() {

    stdl_wavefunction * wf = NULL;
    stdl_basis * bs = NULL;
    read_fchk("../tests/test_files/water_631g.fchk", &wf, &bs);

    stdl_context* ctx = NULL;
    ASSERT_STDL_OK(stdl_context_new(wf, bs, 2.0, 4.0, 12. / STDL_CONST_AU_TO_EV, 1e-4, 1.0, &ctx));
    ASSERT_STDL_OK(stdl_context_select_csfs_monopole(ctx, 1));

    // compute dipole integrals and convert to MO
    double* dipoles_mat = malloc(3 * STDL_MATRIX_SP_SIZE(ctx->nmo) * sizeof(double));
    TEST_ASSERT_NOT_NULL(dipoles_mat);

    make_int1e_MO(wf, bs, STDL_OP_DIPL, -1., ctx, dipoles_mat);

    // build egrad
    float* egrad = malloc(3 * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(egrad);

    stdl_response_perturbed_gradient(ctx, 3, 1, dipoles_mat, egrad);

    // solve response
    size_t nw = 3;
    float w[] = {0, 4.282270E-2f, 4.282270E-2f * 2};

    float* X = malloc(nw * 3 * ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(X);

    float* Y = malloc(nw * 3 * ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(Y);

    ASSERT_STDL_OK(stdl_response_TD_linear(ctx, nw, w, 3, 1, egrad, X, Y));

    // compute polarizabilities
    float alpha[9], alpha_iso, alpha_aniso, result[] = {4.061f, 4.103f, 4.236f};

    for (size_t iw = 0; iw < nw; ++iw) {
        ASSERT_STDL_OK(stdl_property_tensor_linear(
                ctx,
                (stdl_lrv[]) {
                        {STDL_OP_DIPL, dipoles_mat, -w[iw], Y + iw * 3 * ctx->ncsfs, X + iw * 3 * ctx->ncsfs},
                        {STDL_OP_DIPL, dipoles_mat, w[iw], X + iw * 3 * ctx->ncsfs, Y + iw * 3 * ctx->ncsfs}
                },
                alpha
        ));
        stdl_qexp_polarizability(alpha, &alpha_iso, &alpha_aniso);
        TEST_ASSERT_FLOAT_WITHIN(1e-2, result[iw], alpha_iso);
    }

    STDL_FREE_ALL(dipoles_mat, egrad, X, Y);

    ASSERT_STDL_OK(stdl_context_delete(ctx));
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

    make_int1e_MO(wf, bs, STDL_OP_DIPL, -1., ctx, dipoles_mat);

    // get transition dipoles
    float* tdipstd = malloc(ctx->ncsfs * 3 * sizeof(float ));
    stdl_property_transition_dipoles(ctx, ctx->ncsfs, dipoles_mat, Xamptd, Yamptd, tdipstd);

    // build egrad
    float* egrad = malloc(3 * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(egrad);

    stdl_response_perturbed_gradient(ctx, 3, 1, dipoles_mat, egrad);

    // solve linear response
    size_t nw = 3;
    float w[] = {0, STDL_CONST_HC / 1064.f, STDL_CONST_HC / 532.f};

    float* Xtd = malloc(nw * 3 * ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(Xtd);

    float* Ytd = malloc(nw * 3 * ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(Ytd);

    ASSERT_STDL_OK(stdl_response_TD_linear(ctx, nw, w, 3, 1, egrad, Xtd, Ytd));

    // compute polarizabilities
    float alpha[9], alpha_zz;

    for (size_t iw = 0; iw < nw; ++iw) {
        ASSERT_STDL_OK(stdl_property_tensor_linear(
                ctx,
                (stdl_lrv[]) {
                        {STDL_OP_DIPL, dipoles_mat, -w[iw], Ytd + iw * 3 * ctx->ncsfs, Xtd + iw * 3 * ctx->ncsfs},
                        {STDL_OP_DIPL, dipoles_mat, w[iw], Xtd + iw * 3 * ctx->ncsfs, Ytd + iw * 3 * ctx->ncsfs}
                },
                alpha
        ));

        // stdl_matrix_ssp_print(3, alpha, "alpha");

        // compute alpha_zz from SOS
        alpha_zz = .0f;
        for (size_t iexci = 0; iexci < ctx->ncsfs; ++iexci)
            alpha_zz += powf(tdipstd[2 * ctx->ncsfs + iexci], 2) / (etd[iexci] - w[iw]) + powf(tdipstd[2 * ctx->ncsfs + iexci], 2) / (etd[iexci] + w[iw]);

        TEST_ASSERT_FLOAT_WITHIN(1e-2, alpha[8], alpha_zz);
    }

    STDL_FREE_ALL(etd, Xamptd, Yamptd, dipoles_mat, tdipstd, egrad, Xtd, Ytd);
    ASSERT_STDL_OK(stdl_context_delete(ctx));
}


void test_property_tensor_polarizability_TD_fchk_vs_molden_ok() {
    // FCHK
    stdl_wavefunction * wf_fchk = NULL;
    stdl_basis * bs_fchk = NULL;
    read_fchk("../tests/test_files/water_631gdf.fchk", &wf_fchk, &bs_fchk);

    stdl_context* ctx_fchk = NULL;
    ASSERT_STDL_OK(stdl_context_new(wf_fchk, bs_fchk, 2.0, 4.0, 12. / STDL_CONST_AU_TO_EV, 1e-4, 1.0, &ctx_fchk));
    ASSERT_STDL_OK(stdl_context_select_csfs_monopole(ctx_fchk, 1));

    double* dipoles_mat_fchk = malloc(3 * STDL_MATRIX_SP_SIZE(ctx_fchk->nmo) * sizeof(double));
    TEST_ASSERT_NOT_NULL(dipoles_mat_fchk);

    make_int1e_MO(wf_fchk, bs_fchk, STDL_OP_DIPL, -1., ctx_fchk, dipoles_mat_fchk);

    float* egrad_fchk = malloc(3 * ctx_fchk->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(egrad_fchk);

    stdl_response_perturbed_gradient(ctx_fchk, 3, 1, dipoles_mat_fchk, egrad_fchk);

    size_t nw = 3;
    float w[] = {0, 4.282270E-2f, 2 * 4.282270E-2f};

    float* X_fchk = malloc(nw * 3 * ctx_fchk->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(X_fchk);

    float* Y_fchk = malloc(nw * 3 * ctx_fchk->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(Y_fchk);

    ASSERT_STDL_OK(stdl_response_TD_linear(ctx_fchk, nw, w, 3, 1, egrad_fchk, X_fchk, Y_fchk));

    // molden
    stdl_wavefunction * wf_molden = NULL;
    stdl_basis * bs_molden = NULL;
    read_molden("../tests/test_files/water_631gdf.molden", &wf_molden, &bs_molden);

    stdl_context* ctx_molden = NULL;
    ASSERT_STDL_OK(stdl_context_new(wf_molden, bs_molden, 2.0, 4.0, 12. / STDL_CONST_AU_TO_EV, 1e-4, 1.0, &ctx_molden));
    ASSERT_STDL_OK(stdl_context_select_csfs_monopole(ctx_molden, 1));

    double* dipoles_mat_molden = malloc(3 * STDL_MATRIX_SP_SIZE(ctx_molden->nmo) * sizeof(double));
    TEST_ASSERT_NOT_NULL(dipoles_mat_molden);

    make_int1e_MO(wf_molden, bs_molden, STDL_OP_DIPL, -1., ctx_molden, dipoles_mat_molden);

    float* egrad_molden = malloc(3 * ctx_molden->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(egrad_molden);

    stdl_response_perturbed_gradient(ctx_molden, 3, 1, dipoles_mat_molden, egrad_molden);

    float* X_molden = malloc(nw * 3 * ctx_molden->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(X_molden);

    float* Y_molden = malloc(nw * 3 * ctx_molden->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(Y_molden);

    ASSERT_STDL_OK(stdl_response_TD_linear(ctx_molden, nw, w, 3, 1, egrad_molden, X_molden, Y_molden));

    // compute polarizabilities
    float alpha[9], alpha_iso_fchk, alpha_aniso_fchk, alpha_iso_molden, alpha_aniso_molden;

    for (size_t iw = 0; iw < nw; ++iw) {

        ASSERT_STDL_OK(stdl_property_tensor_linear(
                ctx_fchk,
                (stdl_lrv[]) {
                        {STDL_OP_DIPL, dipoles_mat_fchk, -w[iw], Y_fchk + iw * 3 * ctx_fchk->ncsfs, X_fchk+ iw * 3 * ctx_fchk->ncsfs},
                        {STDL_OP_DIPL, dipoles_mat_fchk, w[iw], X_fchk + iw * 3 * ctx_fchk->ncsfs, Y_fchk + iw * 3 * ctx_fchk->ncsfs}
                },
                alpha
        ));
        stdl_qexp_polarizability(alpha, &alpha_iso_fchk, &alpha_aniso_fchk);

        ASSERT_STDL_OK(stdl_property_tensor_linear(
                ctx_molden,
                (stdl_lrv[]) {
                        {STDL_OP_DIPL, dipoles_mat_molden, -w[iw], Y_molden + iw * 3 * ctx_molden->ncsfs, X_molden + iw * 3 * ctx_molden->ncsfs},
                        {STDL_OP_DIPL, dipoles_mat_molden, w[iw], X_molden + iw * 3 * ctx_molden->ncsfs, Y_molden + iw * 3 * ctx_molden->ncsfs}
                },
                alpha
        ));

        stdl_qexp_polarizability(alpha, &alpha_iso_molden, &alpha_aniso_molden);

        TEST_ASSERT_FLOAT_WITHIN(1e-2, alpha_iso_molden, alpha_iso_fchk);
        TEST_ASSERT_FLOAT_WITHIN(1e-2, alpha_aniso_molden, alpha_aniso_fchk);
    }

    STDL_FREE_ALL(dipoles_mat_fchk, egrad_fchk, X_fchk, Y_fchk);
    ASSERT_STDL_OK(stdl_context_delete(ctx_fchk));

    STDL_FREE_ALL(dipoles_mat_molden, egrad_molden, X_molden, Y_molden);
    ASSERT_STDL_OK(stdl_context_delete(ctx_molden));
}


void test_response_polarizability_TDA_ok() {
    stdl_wavefunction * wf = NULL;
    stdl_basis * bs = NULL;
    read_fchk("../tests/test_files/water_631g.fchk", &wf, &bs);

    stdl_context* ctx = NULL;
    ASSERT_STDL_OK(stdl_context_new(wf, bs, 2.0, 4.0, 12. / STDL_CONST_AU_TO_EV, 1e-4, 1.0, &ctx));
    ASSERT_STDL_OK(stdl_context_select_csfs_monopole(ctx, 1));

    // compute dipole integrals and convert to MO
    double* dipoles_mat = malloc(3 * STDL_MATRIX_SP_SIZE(ctx->nmo) * sizeof(double));
    TEST_ASSERT_NOT_NULL(dipoles_mat);

    make_int1e_MO(wf, bs, STDL_OP_DIPL, -1., ctx, dipoles_mat);

    // build egrad
    float* egrad = malloc(3 * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(egrad);

    stdl_response_perturbed_gradient(ctx, 3, 1, dipoles_mat, egrad);

    // solve response
    size_t nw = 3;
    float w[] = {0, STDL_CONST_HC / 1064.f, STDL_CONST_HC / 532.f};

    // solve
    float* Xtda = malloc(nw * 3 * ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(Xtda);

    float* Ytda = malloc(nw * 3 * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(Ytda);

    // compute polarizabilities
    float alpha[9], alpha_iso, alpha_aniso, result[] = {4.128f, 4.170f, 4.303f};
    ASSERT_STDL_OK(stdl_response_TDA_linear(ctx, nw, w, 3, 1, egrad, Xtda, Ytda));

    for (size_t iw = 0; iw < nw; ++iw) {
        ASSERT_STDL_OK(stdl_property_tensor_linear(
                ctx,
                (stdl_lrv[]) {
                        {STDL_OP_DIPL, dipoles_mat, -w[iw], Ytda + iw * 3 * ctx->ncsfs, Xtda + iw * 3 * ctx->ncsfs},
                        {STDL_OP_DIPL, dipoles_mat, w[iw], Xtda + iw * 3 * ctx->ncsfs, Ytda + iw * 3 * ctx->ncsfs}
                },
                alpha
        ));

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

    make_int1e_MO(wf, bs, STDL_OP_DIPL, -1., ctx, dipoles_mat);

    // build egrad
    float* egrad = malloc(3 * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(egrad);

    stdl_response_perturbed_gradient(ctx, 3, 1, dipoles_mat, egrad);

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

    ASSERT_STDL_OK(stdl_response_TDA_linear(ctx, nw, w, 3, 1, egrad, Xtda, Ytda));

    // compute polarizabilities
    float alpha[9], alpha_zz;

    for (size_t iw = 0; iw < nw; ++iw) {
        ASSERT_STDL_OK(stdl_property_tensor_linear(
                ctx,
                (stdl_lrv[]) {
                        {STDL_OP_DIPL, dipoles_mat, -w[iw], Ytda + iw * 3 * ctx->ncsfs, Xtda + iw * 3 * ctx->ncsfs},
                        {STDL_OP_DIPL, dipoles_mat, w[iw], Xtda + iw * 3 * ctx->ncsfs, Ytda + iw * 3 * ctx->ncsfs}
                },
                alpha
        ));

        // stdl_matrix_ssp_print(3, alpha, "alpha");

        // compute alpha_zz from SOS
        alpha_zz = .0f;
        for (size_t iexci = 0; iexci < ctx->ncsfs; ++iexci)
            alpha_zz += powf(tdipstda[2 * ctx->ncsfs + iexci], 2) / (etda[iexci] - w[iw]) + powf(tdipstda[2 * ctx->ncsfs + iexci], 2) / (w[iw] + etda[iexci]);

        TEST_ASSERT_FLOAT_WITHIN(1e-2, alpha[8], alpha_zz);
    }

    STDL_FREE_ALL(dipoles_mat, etda, Xamptda, egrad, Xtda, Ytda, tdipstda);

    ASSERT_STDL_OK(stdl_context_delete(ctx));
}
