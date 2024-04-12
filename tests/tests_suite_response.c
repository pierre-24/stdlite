#include <stdlite/response.h>
#include <stdlite/utils/matrix.h>
#include <stdlite/helpers.h>
#include <stdlite/integrals.h>

#include <string.h>

#include "tests_suite.h"

void setUp() {
    stdl_set_debug_level(-1);
    stdl_set_log_level(0);
}


// Use the spectral representation of a linear response vector to re-create `Wlin`, which is either X(w) or Y(w) depending on `getX`
void _make_response_vector(stdl_context *ctx, float w, int is_hermitian, double *ints_MO, float *e, float *Xamp, float *Yamp, int getX, float *Wlin) {
    size_t nvirt = ctx->nmo - ctx->nocc;

    float mu_ia, x_ia, y_ia;

    for (size_t lia = 0; lia < ctx->ncsfs; ++lia) {
        for (int zeta = 0; zeta < 3; ++zeta)
            Wlin[lia * 3 + zeta] = .0f;
    }

    for (int zeta = 0; zeta < 3; ++zeta) {
        for (size_t iexci = 0; iexci < ctx->ncsfs; ++iexci) {
            for (size_t lia = 0; lia < ctx->ncsfs; ++lia) {
                size_t i = ctx->csfs[lia] / nvirt, a = ctx->csfs[lia] % nvirt + ctx->nocc;

                mu_ia = (!is_hermitian ? -1.f : 1.f) * (float) ints_MO[zeta * STDL_MATRIX_SP_SIZE(ctx->nmo) + STDL_MATRIX_SP_IDX(i, a)];

                x_ia = Xamp[iexci * ctx->ncsfs + lia];
                if (Yamp != NULL)
                    y_ia = Yamp[iexci * ctx->ncsfs + lia];
                else
                    y_ia = .0f;

                if(is_hermitian)
                    Wlin[lia * 3 + zeta] += mu_ia * (x_ia + y_ia) * (((getX ? x_ia : y_ia) / (w - e[iexci])) - ((getX ? y_ia : x_ia) / (w + e[iexci])));
                else
                    Wlin[lia * 3 + zeta] += mu_ia * (x_ia - y_ia) * (((getX ? x_ia : y_ia) / (w - e[iexci])) + ((getX ? y_ia : x_ia) / (w + e[iexci])));

            }
        }
    }
}

void test_response_egrad_non_hermitian_ok() {
    stdl_wavefunction *wf = NULL;
    stdl_basis *bs = NULL;
    read_fchk("../tests/test_files/water_sto3g.fchk", &wf, &bs);

    stdl_context *ctx = NULL;
    ASSERT_STDL_OK(stdl_context_new(wf, bs, 2.0, 4.0, 20. / STDL_CONST_AU_TO_EV, 1e-4, 1.0, &ctx));

    ASSERT_STDL_OK(stdl_context_select_csfs_monopole(ctx, 0));

    // compute integrals and convert to MO
    double *op_mat_sp = malloc(3 * STDL_MATRIX_SP_SIZE(wf->nmo) * sizeof(double));
    TEST_ASSERT_NOT_NULL(op_mat_sp);

    make_int1e_MO(wf, bs, STDL_OP_DIPV, 1., ctx, op_mat_sp);

    // build egrad
    float *egrad = malloc(3 * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(egrad);

    stdl_response_perturbed_gradient(ctx, 3, 0, op_mat_sp, egrad);

    // blow and re-build egrad
    double *op_mat_ge = malloc(3 * ctx->nmo * ctx->nmo * sizeof(double));
    TEST_ASSERT_NOT_NULL(op_mat_ge);

    for (int cpt = 0; cpt < 3; ++cpt)
        ASSERT_STDL_OK(stdl_matrix_dsp_blowge(ctx->nmo, 0, op_mat_sp + cpt * STDL_MATRIX_SP_SIZE(ctx->nmo), op_mat_ge + cpt * ctx->nmo * ctx->nmo));

    float *egrad_from_ge = malloc(3 * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(egrad_from_ge);

    size_t nvirt = ctx->nmo - ctx->nocc;
    for (size_t lia = 0; lia < ctx->ncsfs; ++lia) {
        size_t i = ctx->csfs[lia] / nvirt, a = ctx->csfs[lia] % nvirt + ctx->nocc;
        for (size_t cpt = 0; cpt < 3; ++cpt) {
            egrad_from_ge[lia * 3 + cpt] = -2.f * (float) op_mat_ge[cpt * ctx->nmo * ctx->nmo + i * ctx->nmo + a];
        }
    }

    stdl_matrix_sge_print(2, ctx->nmo, 3, egrad, "egrad");
    stdl_matrix_sge_print(2, ctx->nmo, 3, egrad_from_ge, "egrad'");

    // check if equivalent
    TEST_ASSERT_FLOAT_ARRAY_WITHIN(1e-6, egrad, egrad_from_ge, 3 * ctx->ncsfs);

    STDL_FREE_ALL(op_mat_sp, op_mat_ge, egrad, egrad_from_ge);
    ASSERT_STDL_OK(stdl_context_delete(ctx));
}

void test_response_TDA_ok() {
    stdl_wavefunction * wf = NULL;
    stdl_basis * bs = NULL;
    read_fchk("../tests/test_files/water_631g.fchk", &wf, &bs);

    stdl_context* ctx = NULL;
    ASSERT_STDL_OK(stdl_context_new(wf, bs, 2.0, 4.0, 12. / STDL_CONST_AU_TO_EV, 1e-4, 1.0, &ctx));

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

    for (size_t lia = 0; lia < ctx->ncsfs; ++lia) {
        // The eigenvalues are more or less the diagonal elements of A.
        TEST_ASSERT_FLOAT_WITHIN(1e-2, ctx->A[STDL_MATRIX_SP_IDX(lia, lia)], etda[lia]);

        // check that X^2 is normed
        float sum = .0f;
        for (size_t kjb = 0; kjb < ctx->ncsfs; kjb++) {
            sum += powf(Xamptda[lia * ctx->ncsfs + kjb], 2);
        }

        TEST_ASSERT_FLOAT_WITHIN(1e-4, 1.0f, sum);
    }

    // solve linear response
    size_t nw = 4;
    float w[] = {0, STDL_CONST_HC / 1064.f, STDL_CONST_HC / 532.f, -STDL_CONST_HC / 532.f};

    float* Xtda = malloc(nw * 3 * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(Xtda);

    float* Ytda = malloc(nw * 3 * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(Ytda);

    ASSERT_STDL_OK(stdl_response_TDA_linear(ctx, nw, w, 3, 1, egrad, Xtda, Ytda));

    // check that static X and Y are equals
    TEST_ASSERT_FLOAT_ARRAY_WITHIN(1e-6, Xtda, Ytda,  3 * ctx->ncsfs);

    // check that X(-w) = Y(w) (and conversely)
    TEST_ASSERT_FLOAT_ARRAY_WITHIN(1e-6, Xtda + 2 * 3 * ctx->ncsfs, Ytda + 3 * 3 * ctx->ncsfs,  3 * ctx->ncsfs);
    TEST_ASSERT_FLOAT_ARRAY_WITHIN(1e-6, Ytda + 2 * 3 * ctx->ncsfs, Xtda + 3 * 3 * ctx->ncsfs,  3 * ctx->ncsfs);

    float* reXtda = malloc(nw * 3 * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(reXtda);

    float* reYtda = malloc(nw * 3 * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(reYtda);

    for (size_t iexci = 0; iexci < nw; ++iexci) {
        _make_response_vector(ctx, w[iexci], 1, dipoles_mat, etda, Xamptda, NULL, 1, reXtda);
        _make_response_vector(ctx, w[iexci], 1, dipoles_mat, etda, Xamptda, NULL, 0, reYtda);

        TEST_ASSERT_FLOAT_ARRAY_WITHIN(1e-1, Xtda + iexci * 3 * ctx->ncsfs, reXtda,  3 * ctx->ncsfs);
        TEST_ASSERT_FLOAT_ARRAY_WITHIN(1e-1, Ytda + iexci * 3 * ctx->ncsfs, reYtda,  3 * ctx->ncsfs);
    }

    STDL_FREE_ALL(dipoles_mat, etda, Xamptda, egrad, Xtda, Ytda, reXtda, reYtda);

    ASSERT_STDL_OK(stdl_context_delete(ctx));
}

void test_response_TDA_nonhermitian_ok() {
    stdl_wavefunction * wf = NULL;
    stdl_basis * bs = NULL;
    read_fchk("../tests/test_files/water_sto3g.fchk", &wf, &bs);

    stdl_context* ctx = NULL;
    ASSERT_STDL_OK(stdl_context_new(wf, bs, 2.0, 4.0, 20. / STDL_CONST_AU_TO_EV, 1e-4, 1.0, &ctx));

    ASSERT_STDL_OK(stdl_context_select_csfs_monopole(ctx, 0));

    // compute dipole integrals and convert to MO
    double* op_ints_MO = malloc(3 * STDL_MATRIX_SP_SIZE(wf->nmo) * sizeof(double));
    TEST_ASSERT_NOT_NULL(op_ints_MO);

    make_int1e_MO(wf, bs, STDL_OP_ANGM, 1., ctx, op_ints_MO);

    // build egrad
    float* egrad = malloc(3 * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(egrad);

    stdl_response_perturbed_gradient(ctx, 3, 0, op_ints_MO, egrad);

    // fetch all excitations
    float* etda = malloc(ctx->ncsfs * sizeof(float ));
    float* Xamptda = malloc(ctx->ncsfs * ctx->ncsfs * sizeof(float ));
    ASSERT_STDL_OK(stdl_response_TDA_casida(ctx, ctx->ncsfs, etda, Xamptda));

    // solve linear response
    size_t nw = 4;
    float w[] = {0, STDL_CONST_HC / 1064.f, STDL_CONST_HC / 532.f, -STDL_CONST_HC / 532.f};

    float* Xtda = malloc(nw * 3 * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(Xtda);

    float* Ytda = malloc(nw * 3 * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(Ytda);

    ASSERT_STDL_OK(stdl_response_TDA_linear(ctx, nw, w, 3, 0, egrad, Xtda, Ytda));

    // compare to spectral representation
    float* reXtda = malloc(nw * 3 * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(reXtda);

    float* reYtda = malloc(nw * 3 * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(reYtda);

    for (size_t iexci = 0; iexci < nw; ++iexci) {
        _make_response_vector(ctx, w[iexci], 0, op_ints_MO, etda, Xamptda, NULL, 1, reXtda);
        _make_response_vector(ctx, w[iexci], 0, op_ints_MO, etda, Xamptda, NULL, 0, reYtda);

        TEST_ASSERT_FLOAT_ARRAY_WITHIN(1e-1, Xtda + iexci * 3 * ctx->ncsfs, reXtda, 3 * ctx->ncsfs);
        TEST_ASSERT_FLOAT_ARRAY_WITHIN(1e-1, Ytda + iexci * 3 * ctx->ncsfs, reYtda, 3 * ctx->ncsfs);
    }

    STDL_FREE_ALL(op_ints_MO, etda, Xamptda, egrad, Xtda, Ytda, reXtda, reYtda);

    ASSERT_STDL_OK(stdl_context_delete(ctx));
}

void test_response_TD_ok() {
    stdl_wavefunction * wf = NULL;
    stdl_basis * bs = NULL;
    read_fchk("../tests/test_files/water_631g.fchk", &wf, &bs);

    stdl_context* ctx = NULL;
    ASSERT_STDL_OK(stdl_context_new(wf, bs, 2.0, 4.0, 12. / STDL_CONST_AU_TO_EV, 1e-4, 1.0, &ctx));

    TEST_ASSERT_EQUAL_INT(ctx->nmo,7);
    TEST_ASSERT_EQUAL_INT(ctx->nocc, 4);

    ASSERT_STDL_OK(stdl_context_select_csfs_monopole(ctx, 1));
    TEST_ASSERT_EQUAL_INT(ctx->ncsfs, 10);

    // fetch excitations
    float* etd = malloc(ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(etd);

    float* Xamptd = malloc(ctx->ncsfs * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(Xamptd);

    float* Yamptd = malloc(ctx->ncsfs * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(Yamptd);

    ASSERT_STDL_OK(stdl_response_TD_casida(ctx, ctx->ncsfs, etd, Xamptd, Yamptd));

    for (size_t lia = 0; lia < ctx->ncsfs; ++lia) {
        // The eigenvalues are more or less the diagonal elements of A.
        TEST_ASSERT_FLOAT_WITHIN(1e-2, ctx->A[STDL_MATRIX_SP_IDX(lia, lia)], etd[lia]);

        // check that X^2-Y^2 is normed
        float sum = .0f;
        for (size_t kjb = 0; kjb < ctx->ncsfs; kjb++) {
            sum += powf(Xamptd[lia * ctx->ncsfs + kjb], 2) - powf(Yamptd[lia * ctx->ncsfs + kjb], 2);
        }

        TEST_ASSERT_FLOAT_WITHIN(1e-4, 1.0f, sum);
    }

    // compute dipole integrals and convert to MO
    double* dipoles_mat = malloc(3 * STDL_MATRIX_SP_SIZE(wf->nmo) * sizeof(double));
    TEST_ASSERT_NOT_NULL(dipoles_mat);

    make_int1e_MO(wf, bs, STDL_OP_DIPL, -1., ctx, dipoles_mat);

    // build egrad
    float* egrad = malloc(3 * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(egrad);

    stdl_response_perturbed_gradient(ctx, 3, 1, dipoles_mat, egrad);

    // solve linear response
    size_t nw = 4;
    float w[] = {0, STDL_CONST_HC / 1064.f, STDL_CONST_HC / 532.f, - STDL_CONST_HC / 532.f};

    float* Xtd = malloc(nw * 3 * ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(Xtd);

    float* Ytd = malloc(nw * 3 * ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(Ytd);

    ASSERT_STDL_OK(stdl_response_TD_linear(ctx, nw, w, 3, 1, egrad, Xtd, Ytd));

    // check that static X and Y are equals
    TEST_ASSERT_FLOAT_ARRAY_WITHIN(1e-6, Xtd, Ytd,  3 * ctx->ncsfs);

    // check that X(-w) = Y(w) (and conversely)
    TEST_ASSERT_FLOAT_ARRAY_WITHIN(1e-6, Xtd + 2 * 3 * ctx->ncsfs, Ytd + 3 * 3 * ctx->ncsfs,  3 * ctx->ncsfs);
    TEST_ASSERT_FLOAT_ARRAY_WITHIN(1e-6, Ytd + 2 * 3 * ctx->ncsfs, Xtd + 3 * 3 * ctx->ncsfs,  3 * ctx->ncsfs);

    // compare to spectral representation
    float* reXtd = malloc(nw * 3 * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(reXtd);

    float* reYtd = malloc(nw * 3 * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(reYtd);

    for (size_t iexci = 0; iexci < nw; ++iexci) {
        _make_response_vector(ctx, w[iexci], 1, dipoles_mat, etd, Xamptd, Yamptd, 1, reXtd);
        _make_response_vector(ctx, w[iexci], 1, dipoles_mat, etd, Xamptd, Yamptd, 0, reYtd);

        TEST_ASSERT_FLOAT_ARRAY_WITHIN(1e-1, Xtd + iexci * 3 * ctx->ncsfs, reXtd,  3 * ctx->ncsfs);
        TEST_ASSERT_FLOAT_ARRAY_WITHIN(1e-1, Ytd + iexci * 3 * ctx->ncsfs, reYtd,  3 * ctx->ncsfs);
    }

    STDL_FREE_ALL(etd, Xamptd, Yamptd, dipoles_mat, egrad, Xtd, Ytd, reXtd, reYtd);
    ASSERT_STDL_OK(stdl_context_delete(ctx));
}

void test_response_TD_nonhermitian_ok() {
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
    double* op_ints_MO = malloc(3 * STDL_MATRIX_SP_SIZE(wf->nmo) * sizeof(double));
    TEST_ASSERT_NOT_NULL(op_ints_MO);

    make_int1e_MO(wf, bs, STDL_OP_ANGM, 1., ctx, op_ints_MO);

    // build egrad
    float* egrad = malloc(3 * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(egrad);

    stdl_response_perturbed_gradient(ctx, 3, 0, op_ints_MO, egrad);

    // solve linear response
    size_t nw = 4;
    float w[] = {0, STDL_CONST_HC / 1064.f, STDL_CONST_HC / 532.f, - STDL_CONST_HC / 532.f};

    float* Xtd = malloc(nw * 3 * ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(Xtd);

    float* Ytd = malloc(nw * 3 * ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(Ytd);

    ASSERT_STDL_OK(stdl_response_TD_linear(ctx, nw, w, 3, 0, egrad, Xtd, Ytd));

    // compare to spectral representation
    float* reXtd = malloc(nw * 3 * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(reXtd);

    float* reYtd = malloc(nw * 3 * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(reYtd);
    for (size_t iexci = 0; iexci < nw; ++iexci) {
        _make_response_vector(ctx, w[iexci], 0, op_ints_MO, etd, Xamptd, Yamptd, 1, reXtd);
        _make_response_vector(ctx, w[iexci], 0, op_ints_MO, etd, Xamptd, Yamptd, 0, reYtd);

        TEST_ASSERT_FLOAT_ARRAY_WITHIN(1e-1, Xtd + iexci * 3 * ctx->ncsfs, reXtd, 3 * ctx->ncsfs);
        TEST_ASSERT_FLOAT_ARRAY_WITHIN(1e-1, Ytd + iexci * 3 * ctx->ncsfs, reYtd, 3 * ctx->ncsfs);
    }

    // check if x(w) = -y(w) (and conversely)
    for(size_t ielm=0; ielm < 3 * ctx->ncsfs; ielm++) {
        TEST_ASSERT_FLOAT_WITHIN(1e-6, Xtd[ielm], -Ytd[ielm]); // static
        TEST_ASSERT_FLOAT_WITHIN(1e-6, Xtd[2 * 3 * ctx->ncsfs + ielm], -Ytd[3 * 3 * ctx->ncsfs +  ielm]); // dynamic
    }

    STDL_FREE_ALL(etd, Xamptd, Yamptd, op_ints_MO, egrad, Xtd, Ytd, reXtd, reYtd);
    ASSERT_STDL_OK(stdl_context_delete(ctx));
}
