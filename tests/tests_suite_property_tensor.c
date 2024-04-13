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

        stdl_lrv lrv = {STDL_OP_DIPL, dipoles_mat, w[iw], X + iw * 3 * ctx->ncsfs, Y + iw * 3 * ctx->ncsfs};

        ASSERT_STDL_OK(stdl_property_tensor_linear(
                ctx,
                (stdl_lrv*[]) {&lrv, &lrv},
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
    read_molden("../tests/test_files/chiral_sto3g.molden", &wf, &bs);

    stdl_context* ctx = NULL;
    ASSERT_STDL_OK(stdl_context_new(wf, bs, 2.0, 4.0, 12. / STDL_CONST_AU_TO_EV, 1e-4, 1.0, &ctx));
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
    ASSERT_STDL_OK(stdl_property_tensor_g2e_moments(ctx, STDL_OP_DIPL, dipoles_mat, ctx->ncsfs, Xamptd, Yamptd, tdipstd));

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
    float alpha[9], alpha_sos[9];

    for (size_t iw = 0; iw < nw; ++iw) {
        stdl_lrv lrv = {STDL_OP_DIPL, dipoles_mat, w[iw], Xtd + iw * 3 * ctx->ncsfs, Ytd + iw * 3 * ctx->ncsfs};

        ASSERT_STDL_OK(stdl_property_tensor_linear(
                ctx,
                (stdl_lrv*[]) {&lrv, &lrv},
                alpha
        ));

        stdl_matrix_sge_print(2, 3, 3, alpha, "alpha");

        // compute property_tensor from SOS
        for (int cpt_i = 0; cpt_i < 3; ++cpt_i) {
            for (int cpt_j = 0; cpt_j < 3; ++cpt_j) {
                alpha_sos[cpt_i * 3 + cpt_j] = .0f;
                for (size_t iexci = 0; iexci < ctx->ncsfs; ++iexci)
                    alpha_sos[cpt_i * 3 + cpt_j] += (tdipstd[cpt_i * ctx->ncsfs + iexci] * tdipstd[cpt_j * ctx->ncsfs + iexci]) / (etd[iexci] - w[iw]) + (tdipstd[cpt_i * ctx->ncsfs + iexci] * tdipstd[cpt_j * ctx->ncsfs + iexci]) / (etd[iexci] + w[iw]);
            }
        }

        stdl_matrix_sge_print(2, 3, 3, alpha_sos, "alpha (SOS)");

        TEST_ASSERT_FLOAT_ARRAY_WITHIN(5e-1, alpha, alpha_sos, 9);
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
    ASSERT_STDL_OK(stdl_context_new(wf_fchk, bs_fchk, 2.0, 4.0, 10. / STDL_CONST_AU_TO_EV, 1e-4, 1.0, &ctx_fchk));
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
        stdl_lrv lrv_fchk = {STDL_OP_DIPL, dipoles_mat_fchk, w[iw], X_fchk + iw * 3 * ctx_fchk->ncsfs, Y_fchk + iw * 3 * ctx_fchk->ncsfs};
        ASSERT_STDL_OK(stdl_property_tensor_linear(
                ctx_fchk,
                (stdl_lrv*[]) {&lrv_fchk, &lrv_fchk},
                alpha
        ));
        stdl_qexp_polarizability(alpha, &alpha_iso_fchk, &alpha_aniso_fchk);

        stdl_lrv lrv_molden = {STDL_OP_DIPL, dipoles_mat_molden, w[iw], X_molden + iw * 3 * ctx_molden->ncsfs, Y_molden + iw * 3 * ctx_molden->ncsfs};
        ASSERT_STDL_OK(stdl_property_tensor_linear(
                ctx_molden,
                (stdl_lrv*[]) {&lrv_molden, &lrv_molden},
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
        stdl_lrv lrv = {STDL_OP_DIPL, dipoles_mat, w[iw], Xtda + iw * 3 * ctx->ncsfs, Ytda + iw * 3 * ctx->ncsfs};
        ASSERT_STDL_OK(stdl_property_tensor_linear(
                ctx,
                (stdl_lrv*[]) {&lrv, &lrv},
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
    read_molden("../tests/test_files/chiral_sto3g.molden", &wf, &bs);

    stdl_context* ctx = NULL;
    ASSERT_STDL_OK(stdl_context_new(wf, bs, 2.0, 4.0, 10. / STDL_CONST_AU_TO_EV, 1e-4, 1.0, &ctx));

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
    ASSERT_STDL_OK(stdl_property_tensor_g2e_moments(ctx, STDL_OP_DIPL, dipoles_mat, ctx->ncsfs, Xamptda, NULL, tdipstda));

    // solve linear response
    size_t nw = 3;
    float w[] = {0, STDL_CONST_HC / 1064.f, STDL_CONST_HC / 532.f};

    float* Xtda = malloc(nw * 3 * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(Xtda);

    float* Ytda = malloc(nw * 3 * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(Ytda);

    ASSERT_STDL_OK(stdl_response_TDA_linear(ctx, nw, w, 3, 1, egrad, Xtda, Ytda));

    // compute polarizabilities
    float alpha[9], alpha_sos[9];

    for (size_t iw = 0; iw < nw; ++iw) {
        stdl_lrv lrv = {STDL_OP_DIPL, dipoles_mat, w[iw], Xtda + iw * 3 * ctx->ncsfs, Ytda + iw * 3 * ctx->ncsfs};
        ASSERT_STDL_OK(stdl_property_tensor_linear(
                ctx,
                (stdl_lrv*[]) {&lrv, &lrv},
                alpha
        ));

        stdl_matrix_sge_print(2, 3, 3, alpha, "alpha");

        // compute property_tensor from SOS
        for (int cpt_i = 0; cpt_i < 3; ++cpt_i) {
            for (int cpt_j = 0; cpt_j < 3; ++cpt_j) {
                alpha_sos[cpt_i * 3 + cpt_j] = .0f;
                for (size_t iexci = 0; iexci < ctx->ncsfs; ++iexci)
                    alpha_sos[cpt_i * 3 + cpt_j] += (tdipstda[cpt_i * ctx->ncsfs + iexci] * tdipstda[cpt_j * ctx->ncsfs + iexci]) / (etda[iexci] - w[iw]) + (tdipstda[cpt_i * ctx->ncsfs + iexci] * tdipstda[cpt_j * ctx->ncsfs + iexci]) / (etda[iexci] + w[iw]);
            }
        }

        stdl_matrix_sge_print(2, 3, 3, alpha_sos, "alpha (SOS)");

        TEST_ASSERT_FLOAT_ARRAY_WITHIN(1e-1, alpha, alpha_sos, 9);
    }

    STDL_FREE_ALL(dipoles_mat, etda, Xamptda, egrad, Xtda, Ytda, tdipstda);

    ASSERT_STDL_OK(stdl_context_delete(ctx));
}

void test_property_tensor_angmom_dipl_TD_SOS_ok() {
    stdl_wavefunction * wf = NULL;
    stdl_basis * bs = NULL;
    read_molden("../tests/test_files/chiral_sto3g.molden", &wf, &bs);

    stdl_context* ctx = NULL;
    ASSERT_STDL_OK(stdl_context_new(wf, bs, 2.0, 4.0, 10. / STDL_CONST_AU_TO_EV, 1e-4, 1.0, &ctx));
    ASSERT_STDL_OK(stdl_context_select_csfs_monopole(ctx, 1));

    // fetch excitations
    float* etd = malloc(ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(etd);

    float* Xamptd = malloc(ctx->ncsfs * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(Xamptd);

    float* Yamptd = malloc(ctx->ncsfs * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(Yamptd);

    ASSERT_STDL_OK(stdl_response_TD_casida(ctx, ctx->ncsfs, etd, Xamptd, Yamptd));

    // compute dipl integrals and convert to MO
    double* dipl_mat = malloc(3 * STDL_MATRIX_SP_SIZE(wf->nmo) * sizeof(double));
    TEST_ASSERT_NOT_NULL(dipl_mat);

    make_int1e_MO(wf, bs, STDL_OP_DIPL, -1., ctx, dipl_mat);

    // get transition dipoles
    float* tg2edipl = malloc(ctx->ncsfs * 3 * sizeof(float ));
    ASSERT_STDL_OK(stdl_property_tensor_g2e_moments(ctx, STDL_OP_DIPL, dipl_mat, ctx->ncsfs, Xamptd, Yamptd, tg2edipl));

    // build egrad
    float* egrad_dipl = malloc(3 * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(egrad_dipl);

    stdl_response_perturbed_gradient(ctx, 3, 1, dipl_mat, egrad_dipl);

    // solve linear response
    size_t nw = 3;
    float w[] = {0, 0.025f, STDL_CONST_HC / 532.f};

    float* X_dipl = malloc(nw * 3 * ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(X_dipl);

    float* Y_dipl = malloc(nw * 3 * ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(Y_dipl);

    ASSERT_STDL_OK(stdl_response_TD_linear(ctx, nw, w, 3, 1, egrad_dipl, X_dipl, Y_dipl));

    // compute angm integrals and convert to MO
    double* angm_mat = malloc(3 * STDL_MATRIX_SP_SIZE(wf->nmo) * sizeof(double));
    TEST_ASSERT_NOT_NULL(angm_mat);

    make_int1e_MO(wf, bs, STDL_OP_ANGM, 1., ctx, angm_mat);

    // get transition dipoles
    float* tg2eangm = malloc(ctx->ncsfs * 3 * sizeof(float ));
    ASSERT_STDL_OK(stdl_property_tensor_g2e_moments(ctx, STDL_OP_ANGM, angm_mat, ctx->ncsfs, Xamptd, Yamptd, tg2eangm));

    // build egrad
    float* egrad_angm = malloc(3 * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(egrad_angm);

    stdl_response_perturbed_gradient(ctx, 3, 0, angm_mat, egrad_angm);

    // solve linear response
    float* X_angm = malloc(nw * 3 * ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(X_angm);

    float* Y_angm = malloc(nw * 3 * ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(Y_angm);

    ASSERT_STDL_OK(stdl_response_TD_linear(ctx, nw, w, 3, 1, egrad_angm, X_angm, Y_angm));

    // compute tensor
    float tensor[9], tensor_sos[9];

    for (size_t iw = 0; iw < nw; ++iw) {
        stdl_lrv lrv0 = {STDL_OP_ANGM, angm_mat, w[iw], X_angm + iw * 3 * ctx->ncsfs, Y_angm + iw * 3 * ctx->ncsfs},
            lrv1 = {STDL_OP_DIPL, dipl_mat, w[iw], X_dipl + iw * 3 * ctx->ncsfs, Y_dipl + iw * 3 * ctx->ncsfs};

        ASSERT_STDL_OK(stdl_property_tensor_linear(
                ctx,
                (stdl_lrv*[]) {&lrv0, &lrv1},
                tensor
        ));

        stdl_matrix_sge_print(2, 3, 3, tensor, "tensor");

        // compute property_tensor from SOS
        for (int cpt_i = 0; cpt_i < 3; ++cpt_i) {
            for (int cpt_j = 0; cpt_j < 3; ++cpt_j) {
                tensor_sos[cpt_i * 3 + cpt_j] = .0f;
                for (size_t iexci = 0; iexci < ctx->ncsfs; ++iexci) {
                    float v = tg2eangm[cpt_i * ctx->ncsfs + iexci] * tg2edipl[cpt_j * ctx->ncsfs + iexci];
                    tensor_sos[cpt_i * 3 + cpt_j] += v * (1 / (etd[iexci] - w[iw]) - 1 / (etd[iexci] + w[iw]));
                }
            }
        }

        stdl_matrix_sge_print(2, 3, 3, tensor_sos, "tensor (SOS)");

        TEST_ASSERT_FLOAT_ARRAY_WITHIN(5e-1, tensor, tensor_sos, 9);
    }

    STDL_FREE_ALL(etd, Xamptd, Yamptd, dipl_mat, tg2edipl, angm_mat, tg2eangm, egrad_dipl, X_dipl, Y_dipl, egrad_angm, X_angm, Y_angm);
    ASSERT_STDL_OK(stdl_context_delete(ctx));
}
