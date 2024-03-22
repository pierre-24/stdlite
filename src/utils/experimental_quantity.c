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
