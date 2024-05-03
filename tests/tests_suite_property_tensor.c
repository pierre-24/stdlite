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
                tensor_sos[cpt_i * dim0 + cpt_j] += -v * (1 / (w - e[iexci]) - trs / (w + e[iexci]));
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

size_t test_linear_n = 4;
stdl_operator test_linear_ops[] = {
        STDL_OP_DIPL, STDL_OP_DIPL,
        STDL_OP_DIPL, STDL_OP_ANGM,
        STDL_OP_DIPV, STDL_OP_ANGM,
        STDL_OP_ANGM, STDL_OP_ANGM,
};

size_t test_linear_nw = 3;
float test_linear_w[] = {
        0.f, .05f, 0.1f,
};

void test_property_linear_tensor_TD_SOS_ok() {
    stdl_wavefunction * wf = NULL;
    stdl_basis * bs = NULL;
    read_molden("../tests/test_files/chiral_sto3g.molden", &wf, &bs);

    stdl_context* ctx = NULL;
    ASSERT_STDL_OK(stdl_context_new(wf, bs, 2.0, 4.0, 25. / STDL_CONST_AU_TO_EV, 1e-4, 1.0, &ctx));
    ASSERT_STDL_OK(stdl_context_select_csfs_monopole(ctx, 1));

    // fetch excitations
    float* eexci = malloc(ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(eexci);

    float* XpYamp = malloc(ctx->ncsfs * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(XpYamp);

    float* XmYamp = malloc(ctx->ncsfs * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(XmYamp);

    ASSERT_STDL_OK(stdl_response_TD_casida(ctx, ctx->ncsfs, eexci, XpYamp, XmYamp));

    for (size_t itest = 0; itest < test_linear_n; ++itest) {
        stdl_operator opA = test_linear_ops[itest * 2], opB = test_linear_ops[itest * 2 + 1];

        // compute integrals in MO basis
        double* opA_ints = malloc(STDL_OPERATOR_DIM[opA] * STDL_MATRIX_SP_SIZE(wf->nmo) * sizeof(double));
        TEST_ASSERT_NOT_NULL(opA_ints);

        make_int1e_MO(wf, bs, opA, opA == STDL_OP_DIPL ? -1 : 1, ctx, opA_ints);

        double* opB_ints = malloc(STDL_OPERATOR_DIM[opB] * STDL_MATRIX_SP_SIZE(wf->nmo) * sizeof(double));
        TEST_ASSERT_NOT_NULL(opB_ints);

        make_int1e_MO(wf, bs, opB, opB == STDL_OP_DIPL ? -1 : 1, ctx, opB_ints);

        // get transition moments
        float* tg2e = malloc(ctx->ncsfs * (STDL_OPERATOR_DIM[opA] + STDL_OPERATOR_DIM[opB]) * sizeof(float ));
        ASSERT_STDL_OK(stdl_property_tensor_g2e_moments(ctx, (stdl_operator []) {opA, opB}, (double*[]) {opA_ints, opB_ints}, ctx->ncsfs, XpYamp, XmYamp, tg2e));

        // compute LRV for opA
        float* egrad = malloc(3 * ctx->ncsfs * sizeof(float));
        TEST_ASSERT_NOT_NULL(egrad);

        ASSERT_STDL_OK(stdl_response_perturbed_gradient(ctx, STDL_OPERATOR_DIM[opA], STDL_OPERATOR_ISSYM[opA], opA_ints, egrad));

        float* XpY_opA = malloc(test_linear_nw * STDL_OPERATOR_DIM[opA] * ctx->ncsfs * sizeof(float ));
        TEST_ASSERT_NOT_NULL(XpY_opA);

        float* XmY_opA = malloc(test_linear_nw * STDL_OPERATOR_DIM[opA] * ctx->ncsfs * sizeof(float ));
        TEST_ASSERT_NOT_NULL(XmY_opA);

        ASSERT_STDL_OK(stdl_response_TD_linear(ctx, test_linear_nw, test_linear_w, STDL_OPERATOR_DIM[opA], STDL_OPERATOR_HERMITIAN[opA], egrad, XpY_opA, XmY_opA));

        // compute LRV for opB
        ASSERT_STDL_OK(stdl_response_perturbed_gradient(ctx, STDL_OPERATOR_DIM[opB], STDL_OPERATOR_ISSYM[opB], opB_ints, egrad));

        float* XpY_opB = malloc(test_linear_nw * STDL_OPERATOR_DIM[opB] * ctx->ncsfs * sizeof(float ));
        TEST_ASSERT_NOT_NULL(XpY_opB);

        float* XmY_opB = malloc(test_linear_nw * STDL_OPERATOR_DIM[opB] * ctx->ncsfs * sizeof(float ));
        TEST_ASSERT_NOT_NULL(XmY_opB);

        ASSERT_STDL_OK(stdl_response_TD_linear(ctx, test_linear_nw, test_linear_w, STDL_OPERATOR_DIM[opB], STDL_OPERATOR_HERMITIAN[opB], egrad, XpY_opB, XmY_opB));

        // compute tensor and compare
        float* tensor_resp = malloc(STDL_OPERATOR_DIM[opA] * STDL_OPERATOR_DIM[opB] * sizeof(float ));
        float* tensor_sos = malloc(STDL_OPERATOR_DIM[opA] * STDL_OPERATOR_DIM[opB] * sizeof(float ));

        for (size_t iw = 0; iw < test_linear_nw; ++iw) {
            stdl_lrv lrvA = {opA, opA_ints, test_linear_w[iw], XpY_opA + iw * 3 * ctx->ncsfs, XmY_opA + iw * 3 * ctx->ncsfs};
            stdl_lrv lrvB = {opB, opB_ints, test_linear_w[iw], XpY_opB + iw * 3 * ctx->ncsfs, XmY_opB + iw * 3 * ctx->ncsfs};

            ASSERT_STDL_OK(stdl_property_tensor_linear(
                    ctx,
                    (stdl_lrv*[]) {&lrvA, &lrvB},
                    tensor_resp
            ));

            stdl_matrix_sge_print(2, 3, 3, tensor_resp, "tensor (resp)");

            _tensor_sos2((stdl_operator[]) {opA, opB}, test_linear_w[iw], ctx->ncsfs, eexci, tg2e, tensor_sos);

            stdl_matrix_sge_print(2, 3, 3, tensor_sos, "tensor (SOS)");

            for (int ielm = 0; ielm < 9; ++ielm) {
                if(fabsf(tensor_resp[ielm]) > 5e-1)
                    TEST_ASSERT_FLOAT_WITHIN(3e-1, 1., tensor_resp[ielm] / tensor_sos[ielm]);
            }
        }

        STDL_FREE_ALL(opA_ints, opB_ints, tg2e, egrad, XpY_opA, XmY_opA, XpY_opB, XmY_opB, tensor_resp, tensor_sos);

    }

    STDL_FREE_ALL(eexci, XpYamp, XmYamp);
    ASSERT_STDL_OK(stdl_context_delete(ctx));
}

void test_property_linear_tensor_TDA_SOS_ok() {
    stdl_wavefunction * wf = NULL;
    stdl_basis * bs = NULL;
    read_molden("../tests/test_files/chiral_sto3g.molden", &wf, &bs);

    stdl_context* ctx = NULL;
    ASSERT_STDL_OK(stdl_context_new(wf, bs, 2.0, 4.0, 25. / STDL_CONST_AU_TO_EV, 1e-4, 1.0, &ctx));
    ASSERT_STDL_OK(stdl_context_select_csfs_monopole(ctx, 0));

    // fetch excitations
    float* eexci = malloc(ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(eexci);

    float* XpYamp = malloc(ctx->ncsfs * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(XpYamp);

    ASSERT_STDL_OK(stdl_response_TDA_casida(ctx, ctx->ncsfs, eexci, XpYamp));

    for (size_t itest = 0; itest < test_linear_n; ++itest) {
        stdl_operator opA = test_linear_ops[itest * 2], opB = test_linear_ops[itest * 2 + 1];

        // compute integrals in MO basis
        double* opA_ints = malloc(STDL_OPERATOR_DIM[opA] * STDL_MATRIX_SP_SIZE(wf->nmo) * sizeof(double));
        TEST_ASSERT_NOT_NULL(opA_ints);

        make_int1e_MO(wf, bs, opA, opA == STDL_OP_DIPL ? -1 : 1, ctx, opA_ints);

        double* opB_ints = malloc(STDL_OPERATOR_DIM[opB] * STDL_MATRIX_SP_SIZE(wf->nmo) * sizeof(double));
        TEST_ASSERT_NOT_NULL(opB_ints);

        make_int1e_MO(wf, bs, opB, opB == STDL_OP_DIPL ? -1 : 1, ctx, opB_ints);

        // get transition moments
        float* tg2e = malloc(ctx->ncsfs * (STDL_OPERATOR_DIM[opA] + STDL_OPERATOR_DIM[opB]) * sizeof(float ));
        ASSERT_STDL_OK(stdl_property_tensor_g2e_moments(ctx, (stdl_operator []) {opA, opB}, (double*[]) {opA_ints, opB_ints}, ctx->ncsfs, XpYamp, NULL, tg2e));

        // compute LRV for opA
        float* egrad = malloc(3 * ctx->ncsfs * sizeof(float));
        TEST_ASSERT_NOT_NULL(egrad);

        ASSERT_STDL_OK(stdl_response_perturbed_gradient(ctx, STDL_OPERATOR_DIM[opA], STDL_OPERATOR_ISSYM[opA], opA_ints, egrad));

        float* XpY_opA = malloc(test_linear_nw * STDL_OPERATOR_DIM[opA] * ctx->ncsfs * sizeof(float ));
        TEST_ASSERT_NOT_NULL(XpY_opA);

        float* XmY_opA = malloc(test_linear_nw * STDL_OPERATOR_DIM[opA] * ctx->ncsfs * sizeof(float ));
        TEST_ASSERT_NOT_NULL(XmY_opA);

        ASSERT_STDL_OK(stdl_response_TDA_linear(ctx, test_linear_nw, test_linear_w, STDL_OPERATOR_DIM[opA], STDL_OPERATOR_HERMITIAN[opA], egrad, XpY_opA, XmY_opA));

        // compute LRV for opB
        ASSERT_STDL_OK(stdl_response_perturbed_gradient(ctx, STDL_OPERATOR_DIM[opB], STDL_OPERATOR_ISSYM[opB], opB_ints, egrad));

        float* XpY_opB = malloc(test_linear_nw * STDL_OPERATOR_DIM[opB] * ctx->ncsfs * sizeof(float ));
        TEST_ASSERT_NOT_NULL(XpY_opB);

        float* XmY_opB = malloc(test_linear_nw * STDL_OPERATOR_DIM[opB] * ctx->ncsfs * sizeof(float ));
        TEST_ASSERT_NOT_NULL(XmY_opB);

        ASSERT_STDL_OK(stdl_response_TDA_linear(ctx, test_linear_nw, test_linear_w, STDL_OPERATOR_DIM[opB], STDL_OPERATOR_HERMITIAN[opB], egrad, XpY_opB, XmY_opB));

        // compute tensor and compare
        float* tensor_resp = malloc(STDL_OPERATOR_DIM[opA] * STDL_OPERATOR_DIM[opB] * sizeof(float ));
        float* tensor_sos = malloc(STDL_OPERATOR_DIM[opA] * STDL_OPERATOR_DIM[opB] * sizeof(float ));

        for (size_t iw = 0; iw < test_linear_nw; ++iw) {
            stdl_lrv lrvA = {opA, opA_ints, test_linear_w[iw], XpY_opA + iw * 3 * ctx->ncsfs, XmY_opA + iw * 3 * ctx->ncsfs};
            stdl_lrv lrvB = {opB, opB_ints, test_linear_w[iw], XpY_opB + iw * 3 * ctx->ncsfs, XmY_opB + iw * 3 * ctx->ncsfs};

            ASSERT_STDL_OK(stdl_property_tensor_linear(
                    ctx,
                    (stdl_lrv*[]) {&lrvA, &lrvB},
                    tensor_resp
            ));

            stdl_matrix_sge_print(2, 3, 3, tensor_resp, "tensor (resp)");

            _tensor_sos2((stdl_operator[]) {opA, opB}, test_linear_w[iw], ctx->ncsfs, eexci, tg2e, tensor_sos);

            stdl_matrix_sge_print(2, 3, 3, tensor_sos, "tensor (SOS)");

            for (int ielm = 0; ielm < 9; ++ielm) {
                if(fabsf(tensor_resp[ielm]) > 5e-1)
                    TEST_ASSERT_FLOAT_WITHIN(3e-1, 1., tensor_resp[ielm] / tensor_sos[ielm]);
            }
        }

        STDL_FREE_ALL(opA_ints, opB_ints, tg2e, egrad, XpY_opA, XmY_opA, XpY_opB, XmY_opB, tensor_resp, tensor_sos);

    }

    STDL_FREE_ALL(eexci, XpYamp);
    ASSERT_STDL_OK(stdl_context_delete(ctx));
}


void test_response_first_polarizability_TD_ok() {
    TEST_IGNORE_MESSAGE("tmp");

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
    TEST_IGNORE_MESSAGE("tmp");

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
    TEST_IGNORE_MESSAGE("tmp");

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
    /*for(size_t n = 0; n < nexci; n++) {
        printf("1→%ld E=%.3f, f=%f\n",
               n + 1,
               (etd[n] - etd[0]) * STDL_CONST_AU_TO_EV,
               oscillator_strength(6, 0, n, nexci, etd, eg2etd)
               );
    }*/

    STDL_FREE_ALL(etd, XpYamptd, XmYamptd, dipoles_mat, tg2etd, eg2etd);
    ASSERT_STDL_OK(stdl_context_delete(ctx));
}
