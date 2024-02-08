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
    float s2 = sqrtf(2);

    printf("---- -------- Energy -------- ------ Transition dipole ---------\n");
    printf("       (Eh)     (eV)    (nm)      X        Y        Z      fL   \n");
    printf("---- -------- ------- ------- -------- -------- -------- -------\n");
    for (size_t iexci = 0; iexci < nexci; ++iexci) {
        // print energies
        printf("%4ld %8.5f %7.3f %7.2f", iexci + 1, energies[iexci], energies[iexci] * STDL_CONST_AU_TO_EV, STDL_CONST_HCL / energies[iexci]);

        // compute transition dipole
        float dip[3] = {.0f, .0f, .0f};
        for (size_t lia = 0; lia < ctx->ncsfs; ++lia) {
            size_t i = ctx->csfs[lia] / nvirt, a = ctx->csfs[lia] % nvirt + ctx->nocc;

            for (int cpt = 0; cpt < 3; ++cpt)
                dip[cpt] += s2 * amplitudes[iexci * ctx->ncsfs + lia] * ((float) dipoles_MO[cpt * STDL_MATRIX_SP_SIZE(ctx->nmo) + STDL_MATRIX_SP_IDX(i, a)]);
        }

        // print dipole
        printf(" % 8.5f % 8.5f % 8.5f %7.5f",
               dip[0],
               dip[1],
               dip[2],
               energies[iexci] * (powf(dip[0], 2) + powf(dip[1], 2) + powf(dip[2], 2)) * 2.f / 3 /* TODO: 4/3?!? */
               );

        printf("\n");
    }

    printf("----------------------------------------------------------------\n");

    return STDL_ERR_OK;
}

int stdl_response_print_excitations_contribs(stdl_context *ctx, size_t nexci, float* energies, float *amplitudes, float thresh) {
    size_t nvirt = ctx->nmo - ctx->nocc;
    float s2o2 = sqrtf(2) / 2;

    printf("---- -- E --- ------- Contributions -------\n");
    for (size_t iexci = 0; iexci < nexci; ++iexci) {
        // print energies
        printf("%4ld %8.5f", iexci + 1, energies[iexci]);

        // print contributions
        size_t printed = 0;
        for (size_t kia = 0; kia < ctx->ncsfs; ++kia) {
            float c = powf(amplitudes[iexci * ctx->ncsfs + kia], 2);
            if(c >= thresh) {
                size_t i = ctx->csfs[kia] / nvirt, a = ctx->csfs[kia] % nvirt, hi = ctx->nocc - i - 1;

                if(printed > 0)
                    printf("             ");

                printf(" %5.1f%% (% 6.4f) H", c * 100, amplitudes[iexci * ctx->ncsfs + kia] * s2o2);

                if(hi > 0)
                    printf("-%ld", hi);
                printf("→L");
                if (a > 0)
                    printf("+%ld", a);

                printf("\n");

                printed += 1;
            }
        }

        if(iexci < nexci -1)
            printf("---- -------- -----------------------------\n");
    }

    printf("-------------------------------------------\n");

    return STDL_ERR_OK;
}


// turn ctx->A and ctx-B to A+B and A-B, respectively
void _make_apb_amb(stdl_context* ctx) {
    for (size_t kia = 0; kia < ctx->ncsfs; ++kia) {
        for (size_t kjb = 0; kjb <= kia; ++kjb) {
            float a = ctx->A[STDL_MATRIX_SP_IDX(kia, kjb)], b = ctx->B[STDL_MATRIX_SP_IDX(kia, kjb)];
            ctx->A[STDL_MATRIX_SP_IDX(kia, kjb)] = a + b;
            ctx->B[STDL_MATRIX_SP_IDX(kia, kjb)] = a - b;
        }
    }
}

int stdl_response_casida_TD_full(stdl_context *ctx, float *energies, float *amplitudes) {
    assert(ctx != NULL && ctx->ncsfs > 0 && ctx->B != NULL && energies != NULL && amplitudes != NULL);

    // compute A+B and A-B
    _make_apb_amb(ctx);

    // get (A-B)^1/2
    stdl_matrix_ssp_sqrt(ctx->ncsfs, ctx->B);

    return STDL_ERR_OK;
}
