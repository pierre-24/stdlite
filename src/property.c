#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "stdlite/response.h"
#include "stdlite/logging.h"
#include "stdlite/helpers.h"
#include "stdlite/utils/matrix.h"

int stdl_property_transition_dipoles(stdl_context *ctx, size_t nexci, double* dips_MO, float* X, float* Y, float * tdips) {
    assert(ctx != NULL && nexci > 0 && dips_MO != NULL && X != NULL && tdips != NULL);

    size_t nvirt = ctx->nmo - ctx->nocc;
    float s2 = sqrtf(2);

    for (size_t iexci = 0; iexci < nexci; ++iexci) {
        tdips[0 * nexci + iexci] =  tdips[1 * nexci + iexci] =  tdips[2 * nexci + iexci] = .0f;

        for (size_t lia = 0; lia < ctx->ncsfs; ++lia) {
            size_t i = ctx->csfs[lia] / nvirt, a = ctx->csfs[lia] % nvirt + ctx->nocc;
            float amplitude = X[iexci * ctx->ncsfs + lia];
            if(Y != NULL)
                amplitude += Y[iexci * ctx->ncsfs + lia];

            for (size_t cpt = 0; cpt < 3; ++cpt)
                tdips[cpt * nexci + iexci] += s2 * amplitude * ((float) dips_MO[cpt * STDL_MATRIX_SP_SIZE(ctx->nmo) + STDL_MATRIX_SP_IDX(i, a)]);
        }
    }

    return STDL_ERR_OK;
}

int stdl_property_print_excitations(stdl_context *ctx, size_t nexci, float *energies, float *tdips) {
    assert(ctx != NULL && nexci > 0 && energies != NULL && tdips != NULL);

    printf("---- -------- Energy -------- ------ Transition dipole ---------\n");
    printf("       (Eh)     (eV)    (nm)      X        Y        Z      fL   \n");
    printf("---- -------- ------- ------- -------- -------- -------- -------\n");
    for (size_t iexci = 0; iexci < nexci; ++iexci) {
        // print energies
        printf("%4ld %8.5f %7.3f %7.2f", iexci + 1, energies[iexci], energies[iexci] * STDL_CONST_AU_TO_EV, STDL_CONST_HCL / energies[iexci]);

        // print transition dipole & oscillator strength
        printf(" % 8.5f % 8.5f % 8.5f %7.5f",
               tdips[0 * nexci + iexci],
               tdips[1 * nexci + iexci],
               tdips[2 * nexci + iexci],
               energies[iexci] * (powf(tdips[0 * nexci + iexci], 2) + powf(tdips[1 * nexci + iexci], 2) + powf(tdips[2 * nexci + iexci], 2)) * 2.f / 3
        );

        printf("\n");
    }

    printf("----------------------------------------------------------------\n");

    return STDL_ERR_OK;
}

int stdl_property_print_excitations_contribs(stdl_context *ctx, size_t nexci, float *energies, float *X, float *Y, float thresh) {
    size_t nvirt = ctx->nmo - ctx->nocc;
    float s2o2 = sqrtf(2) / 2;

    printf("---- -- E --- ------- Contributions -------\n");
    for (size_t iexci = 0; iexci < nexci; ++iexci) {
        // print energies
        printf("%4ld %8.5f", iexci + 1, energies[iexci]);

        // print contributions
        size_t printed = 0;
        for (size_t kia = 0; kia < ctx->ncsfs; ++kia) {
            float c = powf(X[iexci * ctx->ncsfs + kia], 2);

            if(Y != NULL)
                c -= powf(Y[iexci * ctx->ncsfs + kia], 2);

            if(c >= thresh) {
                size_t i = ctx->csfs[kia] / nvirt, a = ctx->csfs[kia] % nvirt, hi = ctx->nocc - i - 1;

                if(printed > 0)
                    printf("             ");

                printf(" %5.1f%% (% 6.4f) H", c * 100, X[iexci * ctx->ncsfs + kia] * s2o2);

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