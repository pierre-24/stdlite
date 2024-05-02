#include <stdlite/response.h>
#include <stdlite/utils/matrix.h>
#include <stdlite/helpers.h>
#include <stdlite/utils/experimental_quantity.h>
#include <stdlite/property_tensor.h>
#include <stdlite/integrals.h>

#include <string.h>

#include "tests_suite.h"

void setUp() {
    stdl_set_debug_level(-1);
    stdl_set_log_level(0);
}

void _tensor_sos2(stdl_operator ops[2], float w, size_t nexci, float* e, float* tg2e, float* tensor_sos) {
    size_t dim0 = STDL_OPERATOR_DIM[ops[0]], dim1 = STDL_OPERATOR_DIM[ops[1]];
    float trs = STDL_OPERATOR_TRS[ops[0]] * STDL_OPERATOR_TRS[ops[1]];

    for (size_t cpt_i = 0; cpt_i < dim0; ++cpt_i) {
        for (size_t cpt_j = 0; cpt_j < dim1; ++cpt_j) {
            tensor_sos[cpt_i * dim0 + cpt_j] = .0f;
            for (size_t iexci = 0; iexci < nexci; ++iexci) {
                float v = tg2e[cpt_i * nexci + iexci] * tg2e[(dim0 + cpt_j) * nexci + iexci];
                tensor_sos[cpt_i * dim0 + cpt_j] += -v * (1 / (w - e[iexci]) - trs / (w + e[iexci] ));
            }
        }
    }
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

    float* XpY = malloc(nw * 3 * ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(XpY);

    float* XmY = malloc(nw * 3 * ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(XmY);

    ASSERT_STDL_OK(stdl_response_TD_linear(ctx, nw, w, 3, 1, egrad, XpY, XmY));

    // compute polarizabilities
    float alpha[9], alpha_iso, alpha_aniso, result[] = {4.061f, 4.103f, 4.236f};

    for (size_t iw = 0; iw < nw; ++iw) {

        stdl_lrv lrv = {STDL_OP_DIPL, dipoles_mat, w[iw], XpY + iw * 3 * ctx->ncsfs, XmY + iw * 3 * ctx->ncsfs};

        ASSERT_STDL_OK(stdl_property_tensor_linear(
                ctx,
                (stdl_lrv*[]) {&lrv, &lrv},
                alpha
        ));
        stdl_qexp_polarizability(alpha, &alpha_iso, &alpha_aniso);
        TEST_ASSERT_FLOAT_WITHIN(1e-2, result[iw], alpha_iso);
    }

    STDL_FREE_ALL(dipoles_mat, egrad, XpY, XmY);

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

    float* XpYamptd = malloc(ctx->ncsfs * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(XpYamptd);

    float* XmYamptd = malloc(ctx->ncsfs * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(XmYamptd);

    ASSERT_STDL_OK(stdl_response_TD_casida(ctx, ctx->ncsfs, etd, XpYamptd, XmYamptd));

    // compute dipole integrals and convert to MO
    double* dipoles_mat = malloc(3 * STDL_MATRIX_SP_SIZE(wf->nmo) * sizeof(double));
    TEST_ASSERT_NOT_NULL(dipoles_mat);

    make_int1e_MO(wf, bs, STDL_OP_DIPL, -1., ctx, dipoles_mat);

    // get transition dipoles
    float* tg2etd = malloc(ctx->ncsfs * 6 * sizeof(float ));
    ASSERT_STDL_OK(stdl_property_tensor_g2e_moments(ctx, (stdl_operator []) {STDL_OP_DIPL, STDL_OP_DIPL}, (double*[]) {dipoles_mat, dipoles_mat}, ctx->ncsfs, XpYamptd, XmYamptd, tg2etd));

    // build egrad
    float* egrad = malloc(3 * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(egrad);

    stdl_response_perturbed_gradient(ctx, 3, 1, dipoles_mat, egrad);

    // solve linear response
    size_t nw = 3;
    float w[] = {0, STDL_CONST_HC / 1064.f, STDL_CONST_HC / 532.f};

    float* XpYtd = malloc(nw * 3 * ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(XpYtd);

    float* XmYtd = malloc(nw * 3 * ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(XmYtd);

    ASSERT_STDL_OK(stdl_response_TD_linear(ctx, nw, w, 3, 1, egrad, XpYtd, XmYtd));

    // compute polarizabilities
    float alpha[9], alpha_sos[9];

    for (size_t iw = 0; iw < nw; ++iw) {
        stdl_lrv lrv = {STDL_OP_DIPL, dipoles_mat, w[iw], XpYtd + iw * 3 * ctx->ncsfs, XmYtd + iw * 3 * ctx->ncsfs};

        ASSERT_STDL_OK(stdl_property_tensor_linear(
                ctx,
                (stdl_lrv*[]) {&lrv, &lrv},
                alpha
        ));

        stdl_matrix_sge_print(2, 3, 3, alpha, "alpha");

        // compute tensor from SOS and compare
        _tensor_sos2((stdl_operator[]) {STDL_OP_DIPL, STDL_OP_DIPL}, w[iw], ctx->ncsfs, etd, tg2etd, alpha_sos);

        stdl_matrix_sge_print(2, 3, 3, alpha_sos, "alpha (SOS)");

        TEST_ASSERT_FLOAT_ARRAY_WITHIN(5e-1, alpha, alpha_sos, 9);
    }

    STDL_FREE_ALL(etd, XpYamptd, XmYamptd, dipoles_mat, tg2etd, egrad, XpYtd, XmYtd);
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
    float* tg2etda = malloc(ctx->ncsfs * 6 * sizeof(float ));
    ASSERT_STDL_OK(stdl_property_tensor_g2e_moments(ctx, (stdl_operator []) {STDL_OP_DIPL, STDL_OP_DIPL}, (double*[]) {dipoles_mat, dipoles_mat}, ctx->ncsfs, Xamptda, NULL, tg2etda));

    // solve linear response
    size_t nw = 3;
    float w[] = {0, STDL_CONST_HC / 1064.f, STDL_CONST_HC / 532.f};

    float* XpYtda = malloc(nw * 3 * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(XpYtda);

    float* XmYtda = malloc(nw * 3 * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(XmYtda);

    ASSERT_STDL_OK(stdl_response_TDA_linear(ctx, nw, w, 3, 1, egrad, XpYtda, XmYtda));

    // compute polarizabilities
    float alpha[9], alpha_sos[9];

    for (size_t iw = 0; iw < nw; ++iw) {
        stdl_lrv lrv = {STDL_OP_DIPL, dipoles_mat, w[iw], XpYtda + iw * 3 * ctx->ncsfs, XmYtda + iw * 3 * ctx->ncsfs};
        ASSERT_STDL_OK(stdl_property_tensor_linear(
                ctx,
                (stdl_lrv*[]) {&lrv, &lrv},
                alpha
        ));

        stdl_matrix_sge_print(2, 3, 3, alpha, "alpha");

        // compute tensor from SOS and compare
        _tensor_sos2((stdl_operator[]) {STDL_OP_DIPL, STDL_OP_DIPL}, w[iw], ctx->ncsfs, etda, tg2etda, alpha_sos);

        stdl_matrix_sge_print(2, 3, 3, alpha_sos, "alpha (SOS)");

        TEST_ASSERT_FLOAT_ARRAY_WITHIN(1e-1, alpha, alpha_sos, 9);
    }

    STDL_FREE_ALL(dipoles_mat, etda, Xamptda, egrad, XpYtda, XmYtda, tg2etda);

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
    float* tg2e = malloc(ctx->ncsfs * 6 * sizeof(float ));
    ASSERT_STDL_OK(stdl_property_tensor_g2e_moments(ctx, (stdl_operator []) {STDL_OP_ANGM, STDL_OP_DIPL}, (double*[]) {angm_mat, dipl_mat}, ctx->ncsfs, Xamptd, Yamptd, tg2e));

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

        // compute tensor from SOS
        _tensor_sos2((stdl_operator[]) {STDL_OP_ANGM, STDL_OP_DIPL}, w[iw], ctx->ncsfs, etd, tg2e, tensor_sos);

        stdl_matrix_sge_print(2, 3, 3, tensor_sos, "tensor (SOS)");

        TEST_ASSERT_FLOAT_ARRAY_WITHIN(5e-1, tensor, tensor_sos, 9);
    }

    STDL_FREE_ALL(etd, Xamptd, Yamptd, dipl_mat, angm_mat, tg2e, egrad_dipl, X_dipl, Y_dipl, egrad_angm, X_angm, Y_angm);
    ASSERT_STDL_OK(stdl_context_delete(ctx));
}

void test_response_first_polarizability_TD_ok() {
    stdl_wavefunction * wf = NULL;
    stdl_basis * bs = NULL;
    read_molden("../tests/test_files/chiral_sto3g.molden", &wf, &bs);

    stdl_context* ctx = NULL;
    ASSERT_STDL_OK(stdl_context_new(wf, bs, 4.0, 2.0, 25. / STDL_CONST_AU_TO_EV, 1e-4, 1.0, &ctx));
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
    float w[] = {0, STDL_CONST_HC / 1064.f, STDL_CONST_HC / 1064.f * 2};

    float* XpY = malloc(nw * 3 * ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(XpY);

    float* XmY = malloc(nw * 3 * ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(XmY);

    ASSERT_STDL_OK(stdl_response_TD_linear(ctx, nw, w, 3, 1, egrad, XpY, XmY));

    // compute first hyperpolarizabilities
    float beta[27], beta2_ZZZ, beta2_ZXX, beta2_J1, beta2_J3;

    stdl_lrv lrv0 = {STDL_OP_DIPL, dipoles_mat, w[0], XpY, XmY};
    stdl_lrv lrv1 = {STDL_OP_DIPL, dipoles_mat, w[1], XpY + 1 * 3 * ctx->ncsfs, XmY + 1 * 3 * ctx->ncsfs};
    stdl_lrv lrv2 = {STDL_OP_DIPL, dipoles_mat, w[2], XpY + 2 * 3 * ctx->ncsfs, XmY + 2 * 3 * ctx->ncsfs};

    ASSERT_STDL_OK(stdl_property_tensor_quadratic(
            ctx,
            (stdl_lrv*[]) {&lrv0, &lrv0, &lrv0},
            beta
    ));

    stdl_matrix_sge_print(2, 9, 3, beta, "beta(0)");
    stdl_qexp_first_hyperpolarizability_hrs(beta, &beta2_ZZZ, &beta2_ZXX, &beta2_J1, &beta2_J3);

    TEST_ASSERT_FLOAT_WITHIN(1e-2, 68.606, beta2_ZZZ);
    TEST_ASSERT_FLOAT_WITHIN(1e-2, 25.955, beta2_ZXX);
    TEST_ASSERT_FLOAT_WITHIN(1e-2, 178.042, beta2_J1);
    TEST_ASSERT_FLOAT_WITHIN(1e-2, 577.453, beta2_J3);

    ASSERT_STDL_OK(stdl_property_tensor_quadratic(
            ctx,
            (stdl_lrv*[]) {&lrv2, &lrv1, &lrv1},
            beta
    ));

    stdl_matrix_sge_print(2, 9, 3, beta, "beta(-2w;w,w)");
    stdl_qexp_first_hyperpolarizability_hrs(beta, &beta2_ZZZ, &beta2_ZXX, &beta2_J1, &beta2_J3);

    TEST_ASSERT_FLOAT_WITHIN(1e-2, 77.871, beta2_ZZZ);
    TEST_ASSERT_FLOAT_WITHIN(1e-2, 29.350, beta2_ZXX);
    TEST_ASSERT_FLOAT_WITHIN(1e-2, 199.166, beta2_J1);
    TEST_ASSERT_FLOAT_WITHIN(1e-2, 659.151, beta2_J3);

    STDL_FREE_ALL(dipoles_mat, egrad, XpY, XmY);

    ASSERT_STDL_OK(stdl_context_delete(ctx));
}

void test_response_antisym_TD_ok() {
    stdl_wavefunction * wf = NULL;
    stdl_basis * bs = NULL;
    read_molden("../tests/test_files/chiral_sto3g.molden", &wf, &bs);

    stdl_context* ctx = NULL;
    ASSERT_STDL_OK(stdl_context_new(wf, bs, 4.0, 2.0, 40. / STDL_CONST_AU_TO_EV, 1e-4, 1.0, &ctx));
    ASSERT_STDL_OK(stdl_context_select_csfs_monopole(ctx, 1));

    size_t nw = 3;
    float w[] = {0, STDL_CONST_HC / 1064.f, STDL_CONST_HC / 532.f};

    // compute dipole integrals and convert to MO, then solve response
    double* dipl_mat = malloc(3 * STDL_MATRIX_SP_SIZE(ctx->nmo) * sizeof(double));
    TEST_ASSERT_NOT_NULL(dipl_mat);

    make_int1e_MO(wf, bs, STDL_OP_DIPL, -1., ctx, dipl_mat);

    float* egrad = malloc(3 * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(egrad);

    stdl_response_perturbed_gradient(ctx, 3, 1, dipl_mat, egrad);

    float* XpY_dipl = malloc(nw * 3 * ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(XpY_dipl);

    float* XmY_dipl = malloc(nw * 3 * ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(XmY_dipl);

    ASSERT_STDL_OK(stdl_response_TD_linear(ctx, nw, w, 3, 1, egrad, XpY_dipl, XmY_dipl));

    // compute angm integrals and convert to MO, then solve response
    double* angm_mat = malloc(3 * STDL_MATRIX_SP_SIZE(ctx->nmo) * sizeof(double));
    TEST_ASSERT_NOT_NULL(angm_mat );

    make_int1e_MO(wf, bs, STDL_OP_ANGM, 1., ctx, angm_mat );

    stdl_response_perturbed_gradient(ctx, 3, 1, angm_mat , egrad);

    float* XpY_angm = malloc(nw * 3 * ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(XpY_angm);

    float* XmY_angm = malloc(nw * 3 * ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(XmY_angm);

    ASSERT_STDL_OK(stdl_response_TD_linear(ctx, nw, w, 3, 0, egrad, XpY_angm, XmY_angm));

    // compute tensor
    float tensor[27];

    stdl_lrv lrv0_dipl = {STDL_OP_DIPL, dipl_mat, w[0], XpY_dipl, XmY_dipl};
    stdl_lrv lrv1_dipl = {STDL_OP_DIPL, dipl_mat, w[1], XpY_dipl + 1 * 3 * ctx->ncsfs, XmY_dipl + 1 * 3 * ctx->ncsfs};
    stdl_lrv lrv2_dipl = {STDL_OP_DIPL, dipl_mat, w[2], XpY_dipl + 2 * 3 * ctx->ncsfs, XmY_dipl + 2 * 3 * ctx->ncsfs};
    stdl_lrv lrv0_angm = {STDL_OP_ANGM, angm_mat, w[0], XpY_angm, XmY_angm};

    ASSERT_STDL_OK(stdl_property_tensor_quadratic(
            ctx,
            (stdl_lrv*[]) {&lrv0_dipl, &lrv0_dipl, &lrv0_angm},
            tensor
    ));

    stdl_matrix_sge_print(2, 9, 3, tensor, "<<dipl;dipl,angm>>(0)");

    ASSERT_STDL_OK(stdl_property_tensor_quadratic(
            ctx,
            (stdl_lrv*[]) {&lrv1_dipl, &lrv1_dipl, &lrv0_angm},
            tensor
    ));

    stdl_matrix_sge_print(2, 9, 3, tensor, "<<dipl;dipl,angm>>(w)");

    ASSERT_STDL_OK(stdl_property_tensor_quadratic(
            ctx,
            (stdl_lrv*[]) {&lrv2_dipl, &lrv2_dipl, &lrv0_angm},
            tensor
    ));

    stdl_matrix_sge_print(2, 9, 3, tensor, "<<dipl;dipl,angm>>(2w)");

    STDL_FREE_ALL(dipl_mat, angm_mat, egrad, XpY_dipl, XmY_dipl, XmY_angm, XpY_angm);

    ASSERT_STDL_OK(stdl_context_delete(ctx));
}

float oscillator_strength(size_t shift, size_t m, size_t n, size_t nexci, float* e, float * e2etdips) {
    float x = .0f;
    for (int cpt = 0; cpt < 3; ++cpt) {
        // printf("%ld→%ld[%d] = %f\n", m, n, cpt, e2etdips[cpt * STDL_MATRIX_SP_SIZE(nexci) + STDL_MATRIX_SP_IDX(m, n)]);
        x += powf(e2etdips[(shift+cpt) * STDL_MATRIX_SP_SIZE(nexci) + STDL_MATRIX_SP_IDX(m, n)], 2);
    }

    return 2.f / 3 * (e[n] - e[m]) * x;
}

void test_property_e2e_ok() {
    stdl_wavefunction * wf = NULL;
    stdl_basis * bs = NULL;
    read_molden("../tests/test_files/chiral_sto3g.molden", &wf, &bs);

    stdl_context* ctx = NULL;
    ASSERT_STDL_OK(stdl_context_new(wf, bs, 2.0, 4.0, 20. / STDL_CONST_AU_TO_EV, 1e-4, 1.0, &ctx));
    ASSERT_STDL_OK(stdl_context_select_csfs_monopole(ctx, 1));

    // fetch excitations
    float* etd = malloc(ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(etd);

    float* XpYamptd = malloc(ctx->ncsfs * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(XpYamptd);

    float* XmYamptd = malloc(ctx->ncsfs * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(XmYamptd);

    ASSERT_STDL_OK(stdl_response_TD_casida(ctx, ctx->ncsfs, etd, XpYamptd, XmYamptd));

    // compute dipole integrals and convert to MO
    double* dipoles_mat = malloc(3 * STDL_MATRIX_SP_SIZE(wf->nmo) * sizeof(double));
    TEST_ASSERT_NOT_NULL(dipoles_mat);

    make_int1e_MO(wf, bs, STDL_OP_DIPL, -1., ctx, dipoles_mat);

    // get g2e dipoles
    float* tg2etd = malloc(ctx->ncsfs * 6 * sizeof(float ));
    TEST_ASSERT_NOT_NULL(tg2etd);
    // ASSERT_STDL_OK(stdl_property_tensor_g2e_moments(ctx, (stdl_operator []) {STDL_OP_DIPL, STDL_OP_DIPL}, (double*[]) {dipoles_mat, dipoles_mat}, ctx->ncsfs, XpYamptd, XmYamptd, tg2etd));

    // get ge2e dipoles
    size_t nexci = 23;
    float* eg2etd = malloc(9 * STDL_MATRIX_SP_SIZE(nexci) * sizeof(float ));
    TEST_ASSERT_NOT_NULL(eg2etd);
    ASSERT_STDL_OK(stdl_property_tensor_e2e_moments(ctx, (stdl_operator []) {STDL_OP_DIPL, STDL_OP_DIPL, STDL_OP_DIPL}, (double*[]) {dipoles_mat, dipoles_mat, dipoles_mat}, nexci, XpYamptd, XmYamptd, eg2etd));

    // oscillator strength
    for(size_t n = 0; n < nexci; n++) {
        printf("1→%ld E=%.3f, f=%f\n",
               n + 1,
               (etd[n] - etd[0]) * STDL_CONST_AU_TO_EV,
               oscillator_strength(6, 0, n, nexci, etd, eg2etd)
               );
    }

    STDL_FREE_ALL(etd, XpYamptd, XmYamptd, dipoles_mat, tg2etd, eg2etd);
    ASSERT_STDL_OK(stdl_context_delete(ctx));
}
