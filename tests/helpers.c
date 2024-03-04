#include <stdlite/utils/fchk_parser.h>
#include <stdlite/utils/molden_parser.h>
#include <stdlite/utils/matrix.h>

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


void read_molden(char* molden_path, stdl_wavefunction** wf, stdl_basis** bs) {
    FILE* f = fopen(molden_path, "r");
    TEST_ASSERT_NOT_NULL(f);

    stdl_lexer* lx = NULL;
    ASSERT_STDL_OK(stdl_lexer_new(f, &lx));

    ASSERT_STDL_OK(stdl_molden_parser_extract(lx, wf, bs));

    ASSERT_STDL_OK(stdl_lexer_delete(lx));

    fclose(f);
}

void make_dipoles_MO(stdl_wavefunction* wf, stdl_basis* bs, stdl_context* ctx, double* dipoles_sp_MO) {
    // compute dipole integrals and convert to MO
    double* dipoles_sp_AO = malloc(3 * STDL_MATRIX_SP_SIZE(wf->nao) * sizeof(double));
    TEST_ASSERT_NOT_NULL(dipoles_sp_AO);

    ASSERT_STDL_OK(stdl_basis_dsp_dipole(bs, dipoles_sp_AO));

    for (int cpt = 0; cpt < 3; ++cpt)
        ASSERT_STDL_OK(stdl_wavefunction_dsp_ao_to_dsp_mo(
                wf->nao,
                ctx->nmo,
                ctx->C_ptr,
                dipoles_sp_AO + cpt * STDL_MATRIX_SP_SIZE(wf->nao),
                dipoles_sp_MO + cpt * STDL_MATRIX_SP_SIZE(ctx->nmo))
        );

    free(dipoles_sp_AO);
}
