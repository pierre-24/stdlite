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


// Use the spectral representation of a linear response vector to re-create `Wlin`, which is either (X+Y)(w) or (X-Y)(w) depending on `getXpY`
void _make_response_vector(stdl_context *ctx, float w, int issym, int isherm, double *ints_MO, float *e, float *XpYamp, float *XmYamp, int getXpY, float *Wlin) {
    size_t nvirt = ctx->nmo - ctx->nocc;

    assert(isherm); // not implemented for anti-hermitian operators yet

    float mu_ia, t_ia, u_ia;

    for (size_t lia = 0; lia < ctx->ncsfs; ++lia) {
        for (int zeta = 0; zeta < 3; ++zeta)
            Wlin[lia * 3 + zeta] = .0f;
    }

    for (int zeta = 0; zeta < 3; ++zeta) {
        for (size_t iexci = 0; iexci < ctx->ncsfs; ++iexci) {
            for (size_t lia = 0; lia < ctx->ncsfs; ++lia) {
                size_t i = ctx->csfs[lia] / nvirt, a = ctx->csfs[lia] % nvirt + ctx->nocc;

                mu_ia = (issym ? 1.f : -1.f) * (float) ints_MO[zeta * STDL_MATRIX_SP_SIZE(ctx->nmo) + STDL_MATRIX_SP_IDX(i, a)];

                t_ia = XpYamp[iexci * ctx->ncsfs + lia];
                if (XmYamp != NULL)
                    u_ia = XmYamp[iexci * ctx->ncsfs + lia];
                else
                    u_ia = XpYamp[iexci * ctx->ncsfs + lia];

                Wlin[lia * 3 + zeta] += mu_ia * t_ia * (getXpY? t_ia: u_ia) * (1 / (w - e[iexci]) + (getXpY ? -1.f : 1.f) / (w + e[iexci]));

            }
        }
    }
}

void test_response_egrad_antisym_ok() {
    stdl_wavefunction *wf = NULL;
    stdl_basis *bs = NULL;
    read_molden("../tests/test_files/chiral_sto3g.molden", &wf, &bs);

    stdl_context *ctx = NULL;
    ASSERT_STDL_OK(stdl_context_new(wf, bs, 2.0, 4.0, 20. / STDL_CONST_AU_TO_EV, 1e-4, 1.0, &ctx));

    ASSERT_STDL_OK(stdl_context_select_csfs_monopole(ctx, 0));

    // compute integrals and convert to MO
    double *op_mat_sp = malloc(3 * STDL_MATRIX_SP_SIZE(wf->nmo) * sizeof(double));
    TEST_ASSERT_NOT_NULL(op_mat_sp);

    make_int1e_MO(wf, bs, STDL_OP_ANGM, ctx, op_mat_sp);

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

    make_int1e_MO(wf, bs, STDL_OP_DIPL, ctx, dipoles_mat);

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
        TEST_ASSERT_FLOAT_WITHIN(1e-2, ctx->ApB[STDL_MATRIX_SP_IDX(lia, lia)], etda[lia]);

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

    float* XpYtda = malloc(nw * 3 * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(XpYtda);

    float* XmYtda = malloc(nw * 3 * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(XmYtda);

    ASSERT_STDL_OK(stdl_response_TDA_linear(ctx, nw, w, 3, 1, egrad, XpYtda, XmYtda));

    // check some equalities
    for (size_t ielm = 0; ielm < 3 * ctx->ncsfs; ++ielm) {
        TEST_ASSERT_FLOAT_WITHIN(1e-6, 0, XmYtda[ielm]); // check that (X-Y)(0) = 0
        TEST_ASSERT_FLOAT_WITHIN(1e-6, XpYtda[2 * 3 * ctx->ncsfs + ielm], XpYtda[3 * 3 * ctx->ncsfs + ielm]); // check that (X+Y)(-w) = (X+Y)(w)
        TEST_ASSERT_FLOAT_WITHIN(1e-6, XmYtda[2 * 3 * ctx->ncsfs + ielm], -XmYtda[3 * 3 * ctx->ncsfs + ielm]); // check that (X-Y)(-w) = -(X-Y)(w)
    }

    // compare to spectral representation
    float* reXpYtda = malloc(nw * 3 * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(reXpYtda);

    float* reYmYtda = malloc(nw * 3 * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(reYmYtda);

    for (size_t iexci = 0; iexci < nw; ++iexci) {
        _make_response_vector(ctx, w[iexci], 1, 1, dipoles_mat, etda, Xamptda, NULL, 1, reXpYtda);
        _make_response_vector(ctx, w[iexci], 1, 1, dipoles_mat, etda, Xamptda, NULL, 0, reYmYtda);

        TEST_ASSERT_FLOAT_ARRAY_WITHIN(1e-1, XpYtda + iexci * 3 * ctx->ncsfs, reXpYtda, 3 * ctx->ncsfs);
        TEST_ASSERT_FLOAT_ARRAY_WITHIN(1e-1, XmYtda + iexci * 3 * ctx->ncsfs, reYmYtda, 3 * ctx->ncsfs);
    }

    STDL_FREE_ALL(dipoles_mat, etda, Xamptda, egrad, XpYtda, XmYtda, reXpYtda, reYmYtda);

    ASSERT_STDL_OK(stdl_context_delete(ctx));
}

void test_response_TDA_antisym_ok() {
    stdl_wavefunction * wf = NULL;
    stdl_basis * bs = NULL;
    read_fchk("../tests/test_files/water_sto3g.fchk", &wf, &bs);

    stdl_context* ctx = NULL;
    ASSERT_STDL_OK(stdl_context_new(wf, bs, 2.0, 4.0, 20. / STDL_CONST_AU_TO_EV, 1e-4, 1.0, &ctx));

    ASSERT_STDL_OK(stdl_context_select_csfs_monopole(ctx, 0));

    // compute dipole integrals and convert to MO
    double* op_ints_MO = malloc(3 * STDL_MATRIX_SP_SIZE(wf->nmo) * sizeof(double));
    TEST_ASSERT_NOT_NULL(op_ints_MO);

    make_int1e_MO(wf, bs, STDL_OP_ANGM, ctx, op_ints_MO);

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

    ASSERT_STDL_OK(stdl_response_TDA_linear(ctx, nw, w, 3, 1, egrad, Xtda, Ytda));

    // compare to spectral representation
    float* reXtda = malloc(nw * 3 * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(reXtda);

    float* reYtda = malloc(nw * 3 * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(reYtda);

    for (size_t iexci = 0; iexci < nw; ++iexci) {
        _make_response_vector(ctx, w[iexci], 0, 1, op_ints_MO, etda, Xamptda, NULL, 1, reXtda);
        _make_response_vector(ctx, w[iexci], 0, 1, op_ints_MO, etda, Xamptda, NULL, 0, reYtda);

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

    float* XpYamptd = malloc(ctx->ncsfs * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(XpYamptd);

    float* XmYamptd = malloc(ctx->ncsfs * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(XmYamptd);

    ASSERT_STDL_OK(stdl_response_TD_casida(ctx, ctx->ncsfs, etd, XpYamptd, XmYamptd));

    // compute dipole integrals and convert to MO
    double* dipoles_mat = malloc(3 * STDL_MATRIX_SP_SIZE(wf->nmo) * sizeof(double));
    TEST_ASSERT_NOT_NULL(dipoles_mat);

    make_int1e_MO(wf, bs, STDL_OP_DIPL, ctx, dipoles_mat);

    // build egrad
    float* egrad = malloc(3 * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(egrad);

    stdl_response_perturbed_gradient(ctx, 3, 1, dipoles_mat, egrad);

    // solve linear response
    size_t nw = 4;
    float w[] = {0, STDL_CONST_HC / 1064.f, STDL_CONST_HC / 532.f, - STDL_CONST_HC / 532.f};

    float* XpYtd = malloc(nw * 3 * ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(XpYtd);

    float* XmYtd = malloc(nw * 3 * ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(XmYtd);

    ASSERT_STDL_OK(stdl_response_TD_linear(ctx, nw, w, 3, 1, egrad, XpYtd, XmYtd));

    // check some equalities
    for (size_t ielm = 0; ielm < 3 * ctx->ncsfs; ++ielm) {
        TEST_ASSERT_FLOAT_WITHIN(1e-6, 0, XmYtd[ielm]); // check that (X-Y)(0) = 0
        TEST_ASSERT_FLOAT_WITHIN(1e-6, XpYtd[2 * 3 * ctx->ncsfs + ielm], XpYtd[3 * 3 * ctx->ncsfs + ielm]); // check that (X+Y)(-w) = (X+Y)(w)
        TEST_ASSERT_FLOAT_WITHIN(1e-6, XmYtd[2 * 3 * ctx->ncsfs + ielm], -XmYtd[3 * 3 * ctx->ncsfs + ielm]); // check that (X-Y)(-w) = -(X-Y)(w)
    }

    // compare to spectral representation
    float* reXtd = malloc(nw * 3 * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(reXtd);

    float* reYtd = malloc(nw * 3 * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(reYtd);

   for (size_t iexci = 0; iexci < nw; ++iexci) {
        _make_response_vector(ctx, w[iexci], 1, 1, dipoles_mat, etd, XpYamptd, XmYamptd, 1, reXtd);
        _make_response_vector(ctx, w[iexci], 1, 1, dipoles_mat, etd, XpYamptd, XmYamptd, 0, reYtd);

        TEST_ASSERT_FLOAT_ARRAY_WITHIN(1e-1, XpYtd + iexci * 3 * ctx->ncsfs, reXtd, 3 * ctx->ncsfs);
        TEST_ASSERT_FLOAT_ARRAY_WITHIN(1e-1, XmYtd + iexci * 3 * ctx->ncsfs, reYtd, 3 * ctx->ncsfs);
    }

    STDL_FREE_ALL(etd, XpYamptd, XmYamptd, dipoles_mat, egrad, XpYtd, XmYtd, reXtd, reYtd);
    ASSERT_STDL_OK(stdl_context_delete(ctx));
}

void test_response_TD_antisym_ok() {

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

    make_int1e_MO(wf, bs, STDL_OP_ANGM, ctx, op_ints_MO);

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

    ASSERT_STDL_OK(stdl_response_TD_linear(ctx, nw, w, 3, 1, egrad, Xtd, Ytd));

    // compare to spectral representation
    float* reXtd = malloc(nw * 3 * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(reXtd);

    float* reYtd = malloc(nw * 3 * ctx->ncsfs * sizeof(float));
    TEST_ASSERT_NOT_NULL(reYtd);
    for (size_t iexci = 0; iexci < nw; ++iexci) {
        _make_response_vector(ctx, w[iexci], 0, 1, op_ints_MO, etd, Xamptd, Yamptd, 1, reXtd);
        _make_response_vector(ctx, w[iexci], 0, 1, op_ints_MO, etd, Xamptd, Yamptd, 0, reYtd);

        TEST_ASSERT_FLOAT_ARRAY_WITHIN(1e-1, Xtd + iexci * 3 * ctx->ncsfs, reXtd, 3 * ctx->ncsfs);
        TEST_ASSERT_FLOAT_ARRAY_WITHIN(1e-1, Ytd + iexci * 3 * ctx->ncsfs, reYtd, 3 * ctx->ncsfs);
    }

    STDL_FREE_ALL(etd, Xamptd, Yamptd, op_ints_MO, egrad, Xtd, Ytd, reXtd, reYtd);
    ASSERT_STDL_OK(stdl_context_delete(ctx));
}
