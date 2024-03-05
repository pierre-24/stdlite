#include <assert.h>
#include <stdio.h>
#include <math.h>

#include "stdlite/utils/experimental_quantity.h"
#include "stdlite/utils/matrix.h"
#include "stdlite/logging.h"
#include "stdlite/helpers.h"


int stdl_qexp_polarizability(float * alpha, float* iso, float* aniso)  {
    assert(alpha != NULL && iso != NULL && aniso != NULL);

    *iso = .0f;
    *aniso = .0f;

    for (int i = 0; i < 3; ++i) {
        *iso += alpha[STDL_MATRIX_SP_IDX(i, i)];

        for (int j = 0; j < 3; ++j) {
            *aniso += 3 * powf(alpha[STDL_MATRIX_SP_IDX(i, j)], 2) - alpha[STDL_MATRIX_SP_IDX(i, i)] * alpha[STDL_MATRIX_SP_IDX(j, j)];
        }
    }

    *iso /= 3;
    *aniso = powf(.5f * (*aniso), .5f);

    return STDL_ERR_OK;
}

int stdl_qexp_first_hyperpolarizability_hrs(float beta[3][3][3], float* beta2_ZZZ, float* beta2_ZXX)  {
    assert(beta != NULL && beta2_ZZZ != NULL && beta2_ZXX);

    *beta2_ZZZ = .0f;
    *beta2_ZXX = .0f;

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                *beta2_ZZZ += 2 * powf(beta[i][j][k], 2) + 1 * beta[i][j][j] * beta[i][k][k] + 4 * (beta[i][i][j] * beta[j][k][k] + beta[i][i][j] * beta[k][j][k] + beta[i][j][k] * beta[j][i][k]);
                *beta2_ZXX += 6 * powf(beta[i][j][k], 2) + 3 * beta[i][j][j] * beta[i][k][k] - 2 * (beta[i][i][j] * beta[j][k][k] + beta[i][i][j] * beta[k][j][k] + beta[i][j][k] * beta[j][i][k]);
            }
        }
    }

    *beta2_ZZZ /= 105;
    *beta2_ZXX /= 105;

    return STDL_ERR_OK;
}

int stdl_qexp_excitations_print(stdl_context *ctx, size_t nexci, float *energies, float *tdips) {
    assert(ctx != NULL && nexci > 0 && energies != NULL && tdips != NULL);

    printf("---- -------- Energy -------- ------ Transition dipole ---------\n");
    printf("       (Eh)     (eV)    (nm)      X        Y        Z      fL   \n");
    printf("---- -------- ------- ------- -------- -------- -------- -------\n");
    for (size_t iexci = 0; iexci < nexci; ++iexci) {
        // print energies
        printf("%4ld %8.5f %7.3f %7.2f", iexci + 1, energies[iexci], energies[iexci] * STDL_CONST_AU_TO_EV, STDL_CONST_HC / energies[iexci]);

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

int stdl_qexp_excitations_contribs_print(stdl_context *ctx, size_t nexci, float *energies, float *X, float *Y, float thresh) {
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
