#include <assert.h>
#include <lapacke.h>
#include <string.h>
#include <stdio.h>

#include "stdlite/response.h"
#include "stdlite/logging.h"
#include "stdlite/helpers.h"
#include "stdlite/utils/matrix.h"


int stdl_response_casida_TDA_full(stdl_context *ctx, float *energies, float *amplitudes) {
    assert(ctx != NULL && ctx->ncsfs > 0 && energies != NULL && amplitudes != NULL);

    int err = LAPACKE_sspev(
            LAPACK_ROW_MAJOR, 'V', 'L',
            (int) ctx->ncsfs, ctx->A,
            energies, amplitudes, (int) ctx->ncsfs
            );

    STDL_ERROR_HANDLE_AND_REPORT(err != 0, return STDL_ERR_RESPONSE, "error while sspev(): %d", err);

    stdl_matrix_sge_transpose(ctx->ncsfs, ctx->ncsfs, amplitudes);

    return STDL_ERR_OK;
}

int stdl_response_casida_TDA(stdl_context *ctx, size_t nexci, float *energies, float *amplitudes) {
    assert(ctx != NULL && ctx->ncsfs > 0 && nexci > 0 && energies != NULL && amplitudes != NULL);

    int* ifail = malloc(ctx->ncsfs * sizeof(float ));

    STDL_ERROR_HANDLE_AND_REPORT(ifail == NULL, return STDL_ERR_MALLOC, "malloc");

    int found = 0;

    int err = LAPACKE_sspevx(
            LAPACK_ROW_MAJOR, 'V', 'I', 'L',
            (int) ctx->ncsfs, ctx->A,
            .0f, .0f,
            1 /* even though we are in C, it starts at 1 */, (int) nexci, STDL_RESPONSE_EIGV_ABSTOL,
            &found, energies, amplitudes, (int) nexci, ifail
            );

    STDL_ERROR_HANDLE_AND_REPORT(err != 0, STDL_FREE_ALL(ifail); return STDL_ERR_RESPONSE, "error while sspevx(): %d", err);

    stdl_matrix_sge_transpose(ctx->ncsfs, nexci, amplitudes);

    STDL_FREE_ALL(ifail);

    return STDL_ERR_OK;
}

int stdl_response_print_excitations(stdl_context *ctx, size_t nexci, float *energies, float *amplitudes, double* dipoles_MO) {
    size_t nvirt = ctx->nmo - ctx->nocc;

    printf("---- ------- Energy ------- ------ Transition dipole ---------\n");
    printf("       (Eh)   (eV)    (nm)      X        Y        Z       fL  \n");
    printf("---- ------- ------ ------- -------- -------- -------- -------\n");
    for (size_t iexci = 0; iexci < nexci; ++iexci) {
        // print energies
        printf("%4ld %7.5f %6.3f %7.2f", iexci + 1, energies[iexci], energies[iexci] * STDL_CONST_AU_TO_EV, STDL_CONST_HCL / energies[iexci]);

        // compute transition dipole
        float dip[3] = {.0f, .0f, .0f};
        for (size_t kia = 0; kia < ctx->ncsfs; ++kia) {
            size_t i = ctx->csfs[kia] / nvirt, a = ctx->csfs[kia] % nvirt + ctx->nocc;

            for (int cpt = 0; cpt < 3; ++cpt)
                dip[cpt] += amplitudes[iexci * ctx->ncsfs + kia] * ((float) dipoles_MO[cpt * STDL_MATRIX_SP_SIZE(ctx->nmo) + STDL_MATRIX_SP_IDX(i, a)]);
        }

        // print dipole
        printf(" % 8.5f % 8.5f % 8.5f %7.5f",
               dip[0],
               dip[1],
               dip[2],
               energies[iexci] * (powf(dip[0], 2) + powf(dip[1], 2) + powf(dip[2], 2)) * 4.f / 3 /* TODO: 4/3?!? */
               );

        printf("\n");
    }

    printf("--------------------------------------------------------------\n");

    return STDL_ERR_OK;
}
