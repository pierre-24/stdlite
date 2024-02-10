#include <stdlite/response.h>
#include <stdlite/utils/fchk_parser.h>
#include <stdlite/utils/matrix.h>
#include <stdlite/helpers.h>
#include <stdlite/property.h>
#include <string.h>

#include "tests_suite.h"


void read_fchk(char* fchk_path, stdl_wavefunction** wf, stdl_basis** bs) {

    FILE* f = fopen(fchk_path, "r");
    TEST_ASSERT_NOT_NULL(f);

    stdl_lexer* lx = NULL;
    ASSERT_STDL_OK(stdl_lexer_new(f, &lx));
    ASSERT_STDL_OK(stdl_fchk_parser_skip_intro(lx));

    ASSERT_STDL_OK(stdl_fchk_parser_extract(lx, wf, bs));

    ASSERT_STDL_OK(stdl_lexer_delete(lx));

    fclose(f);
}

void test_response_TDA_casida_ok() {
    stdl_wavefunction * wf = NULL;
    stdl_basis * bs = NULL;
    read_fchk("../tests/test_files/water_631g.fchk", &wf, &bs);

    stdl_context* ctx = NULL;
    ASSERT_STDL_OK(stdl_context_new(wf, bs, 2.0, 4.0, 12. / 27.212, 1e-4, 1.0, &ctx));

    ASSERT_STDL_OK(stdl_context_select_csfs_monopole(ctx, 0));

    // copy A for latter
    float* Ap = malloc(STDL_MATRIX_SP_SIZE(ctx->ncsfs) * sizeof(float));
    TEST_ASSERT_NOT_NULL(Ap);
    memcpy(Ap, ctx->A, STDL_MATRIX_SP_SIZE(ctx->ncsfs) * sizeof(float));

    // fetch all excitations
    float* energies = malloc(ctx->ncsfs * sizeof(float ));
    float* amplitudes = malloc(ctx->ncsfs * ctx->ncsfs * sizeof(float ));
    ASSERT_STDL_OK(stdl_response_TDA_casida(ctx, ctx->ncsfs, energies, amplitudes));

    // replace A, which has now been severely damaged.
    STDL_FREE_ALL(ctx->A);
    ctx->A = Ap;

    // request the 5 first excitations
    size_t nrequested = 5;
    float* first_energies = malloc(nrequested * sizeof(float));
    float* first_amplitudes = malloc(nrequested * ctx->ncsfs * sizeof(float));

    ASSERT_STDL_OK(stdl_response_TDA_casida(ctx, nrequested, first_energies, first_amplitudes));

    for (size_t kia = 0; kia < nrequested; ++kia) {
        // the same eigenvalues should have been obtained
        TEST_ASSERT_FLOAT_WITHIN(1e-5, energies[kia], first_energies[kia]);
    }

    STDL_FREE_ALL(energies, amplitudes, first_energies, first_amplitudes);

    ASSERT_STDL_OK(stdl_context_delete(ctx));
}

void test_response_TDA_print_ok() {
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

    ASSERT_STDL_OK(stdl_context_select_csfs_monopole(ctx, 0));

    // request the first excitations
    size_t nrequested = 10;
    float* e = malloc(nrequested * sizeof(float));
    float* X = malloc(nrequested * ctx->ncsfs * sizeof(float));

    ASSERT_STDL_OK(stdl_response_TDA_casida(ctx, nrequested, e, X));

    // compute the dipole moments
    double* dipoles_sp_AO = malloc(3 * STDL_MATRIX_SP_SIZE(wf->nao) * sizeof(double));
    TEST_ASSERT_NOT_NULL(dipoles_sp_AO);

    double* dipoles_sp_MO = malloc(3 * STDL_MATRIX_SP_SIZE(wf->nmo) * sizeof(double));
    TEST_ASSERT_NOT_NULL(dipoles_sp_MO);

    // compute dipole integrals and convert to MO
    ASSERT_STDL_OK(stdl_basis_dsp_dipole(bs, dipoles_sp_AO));

    for (int cpt = 0; cpt < 3; ++cpt)
        ASSERT_STDL_OK(stdl_wavefunction_dsp_ao_to_dsp_mo(
                wf->nao,
                ctx->nmo,
                ctx->C_orig,
                dipoles_sp_AO + cpt * STDL_MATRIX_SP_SIZE(wf->nao),
                dipoles_sp_MO + cpt * STDL_MATRIX_SP_SIZE(ctx->nmo))
                );

    float* tdips = malloc(nrequested * 3 * sizeof(float ));
    stdl_property_transition_dipoles(ctx, nrequested, dipoles_sp_MO, X, NULL, tdips);

    ASSERT_STDL_OK(stdl_property_print_excitations(ctx, nrequested, e, tdips));
    ASSERT_STDL_OK(stdl_property_print_excitations_contribs(ctx, nrequested, e, X, NULL,.001f));

    STDL_FREE_ALL(e, X, dipoles_sp_AO, dipoles_sp_MO, tdips);

    ASSERT_STDL_OK(stdl_context_delete(ctx));
}

void test_response_RPA_casida_ok() {
    stdl_wavefunction * wf = NULL;
    stdl_basis * bs = NULL;
    read_fchk("../tests/test_files/water_631g.fchk", &wf, &bs);

    stdl_context* ctx = NULL;
    ASSERT_STDL_OK(stdl_context_new(wf, bs, 2.0, 4.0, 12. / 27.212, 1e-4, 1.0, &ctx));

    TEST_ASSERT_EQUAL_INT(ctx->nmo,7);
    TEST_ASSERT_EQUAL_INT(ctx->nocc, 4);

    ASSERT_STDL_OK(stdl_context_select_csfs_monopole(ctx, 1));

    TEST_ASSERT_EQUAL_INT(ctx->ncsfs, 10);

    float* energies = malloc(ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(energies);

    float* X = malloc(ctx->ncsfs * ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(X);

    float* Y = malloc(ctx->ncsfs * ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(Y);

    ASSERT_STDL_OK(stdl_response_RPA_casida(ctx, ctx->ncsfs, energies, X, Y));

    for (size_t kia = 0; kia < ctx->ncsfs; ++kia) {
        // in this case, the eigenvalues are more or less the diagonal elements of A.
        TEST_ASSERT_FLOAT_WITHIN(1e-2, ctx->ecsfs[kia], energies[kia]);

        // check that X^2-Y^2 is normed
        float sum = .0f;
        for (size_t kjb = 0; kjb < ctx->ncsfs; kjb++) {
            sum += powf(X[kia * ctx->ncsfs + kjb], 2) - powf(Y[kia * ctx->ncsfs + kjb], 2);
        }

        TEST_ASSERT_FLOAT_WITHIN(1e-4, 1.0f, sum);
    }

    STDL_FREE_ALL(energies, X, Y);

    ASSERT_STDL_OK(stdl_context_delete(ctx));
}



void test_response_RPA_linear_ok() {
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

    float* X = malloc(3 * ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(X);

    float* Y = malloc(3 * ctx->ncsfs * sizeof(float ));
    TEST_ASSERT_NOT_NULL(Y);

    ASSERT_STDL_OK(stdl_response_RPA_linear(ctx, 4.2822696E-02, 3, egrad, X, Y));

    STDL_FREE_ALL(dipoles_sp_AO, dipoles_sp_MO, egrad, X, Y);

    ASSERT_STDL_OK(stdl_context_delete(ctx));
}
