#include <stdlite/response.h>
#include <stdlite/utils/matrix.h>
#include <stdlite/helpers.h>
#include <stdlite/property.h>
#include <stdlite/utils/experimental_quantity.h>
#include <stdlite/utils/permutations.h>

#include <string.h>

#include "tests_suite.h"
#include "stdlite/integrals.h"

void setUp() {
    stdl_set_debug_level(-1);
    stdl_set_log_level(0);
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
        stdl_response_lr_tensor(
                ctx,
                (size_t[]) {3, 3}, (int[]) {1, 1},
                dipoles_mat,
                Xtd + iw * 3 * ctx->ncsfs,
                Ytd + iw * 3 * ctx->ncsfs,
                0, alpha);

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

    make_int1e_MO(wf, bs, STDL_OP_DIPL, -1., ctx, dipoles_mat);

    // build egrad
    float* egrad = malloc(3 * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(egrad);

    stdl_response_perturbed_gradient(ctx, 3, 1, dipoles_mat, egrad);

    // solve response
    size_t nw = 4;
    float w[] = {0, STDL_CONST_HC / 1064.f, -STDL_CONST_HC / 532.f, STDL_CONST_HC / 532.f};

    float* X = malloc(nw * 3 * ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(X);

    float* Y = malloc(nw * 3 * ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(Y);

    ASSERT_STDL_OK(stdl_response_TD_linear(ctx, nw, w, 3, 1, egrad, X, Y));

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
        stdl_response_lr_tensor(
                ctx,
                (size_t[]) {3, 3}, (int[]) {1, 1},
                dipoles_mat,
                Xtda + iw * 3 * ctx->ncsfs,
                Ytda + iw * 3 * ctx->ncsfs,
                0, alpha);

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

    make_int1e_MO(wf, bs, STDL_OP_DIPL, -1., ctx, dipoles_sp_MO);

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

    make_int1e_MO(wf, bs, STDL_OP_DIPL, -1., ctx, dipoles_mat);

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

    stdl_response_perturbed_gradient(ctx, 3, 1, dipoles_mat, egrad);

    // solve linear response
    size_t nw = 3;
    float w[] = {0, STDL_CONST_HC / 1064.f, -STDL_CONST_HC / 532.f};

    float* Xtd = malloc(nw * 3 * ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(Xtd);

    float* Ytd = malloc(nw * 3 * ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(Ytd);

    ASSERT_STDL_OK(stdl_response_TD_linear(ctx, nw, w, 3, 1, egrad, Xtd, Ytd));

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

    make_int1e_MO(wf, bs, STDL_OP_DIPL, -1., ctx, dipoles_mat);

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

    stdl_response_perturbed_gradient(ctx, 3, 1, dipoles_mat, egrad);

    // solve linear response
    size_t nw = 3;
    float w[] = {0, STDL_CONST_HC / 1064.f, - STDL_CONST_HC / 532.f};

    float* Xtda = malloc(nw * 3 * ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(Xtda);

    float* Ytda = malloc(nw * 3 * ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(Ytda);

    ASSERT_STDL_OK(stdl_response_TDA_linear(ctx, nw, w, 3, 1, egrad, Xtda, Ytda));

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