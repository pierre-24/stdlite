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

    // compute dipole integrals and convert to MO
    double* dipoles_sp_MO = malloc(3 * STDL_MATRIX_SP_SIZE(ctx->nmo) * sizeof(double));
    TEST_ASSERT_NOT_NULL(dipoles_sp_MO);

    make_dipoles_MO(wf, bs, ctx, dipoles_sp_MO);

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

    STDL_FREE_ALL(dipoles_sp_MO, etda, Xtda, tdipstda, etd, Xtd, Ytd, tdipstd);

    ASSERT_STDL_OK(stdl_context_delete(ctx));
}

void test_polarizability_TD_ok() {
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
    float alpha[6], alpha_iso, result[] = {4.061f, 4.103f, 4.236f};

    for (size_t iw = 0; iw < nw; ++iw) {
        stdl_property_polarizability(ctx, dipoles_mat, X + iw * 3 * ctx->ncsfs, Y + iw * 3 * ctx->ncsfs, alpha);
        stdl_property_mean_polarizability(alpha, &alpha_iso);
        TEST_ASSERT_FLOAT_WITHIN(1e-2, result[iw], alpha_iso);
    }

    STDL_FREE_ALL(dipoles_mat, egrad, X, Y);

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

    // solve
    float* Xtda = malloc(nw * 3 * ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(Xtda);

    // compute polarizabilities
    float alpha[6], alpha_iso, result[] = {2.064f, 2.283f, 2.566f};
    ASSERT_STDL_OK(stdl_response_TDA_linear(ctx, nw, w, 3, egrad, Xtda));

    for (size_t iw = 0; iw < nw; ++iw) {
        stdl_property_polarizability(ctx, dipoles_mat, Xtda + iw * 3 * ctx->ncsfs, NULL, alpha);
        stdl_property_mean_polarizability(alpha, &alpha_iso);
        TEST_ASSERT_FLOAT_WITHIN(1e-2, result[iw], alpha_iso);
    }

    STDL_FREE_ALL(dipoles_mat, egrad, Xtda);

    ASSERT_STDL_OK(stdl_context_delete(ctx));
}
