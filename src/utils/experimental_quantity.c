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
        *iso += alpha[i * 3 + i];

        for (int j = 0; j < 3; ++j) {
            *aniso += 3 * powf(alpha[i * 3 + j], 2) - alpha[i * 3 + i] * alpha[j * 3 + j];
        }
    }

    *iso /= 3;
    *aniso = powf(.5f * (*aniso), .5f);

    return STDL_ERR_OK;
}

int stdl_qexp_first_hyperpolarizability(float beta[27], float beta_vec[3]) {
    assert(beta != NULL && beta_vec != NULL);

    for (int i = 0; i < 3; ++i) {
        beta_vec[i] = 0;
        for (int j = 0; j < 3; ++j) {
            beta_vec[i] += beta[i*9 + j*3 + j] + beta[j*9 + i*3 + j] + beta[j*9 + j*3 + i];
        }
        beta_vec[i] *= 1.f/3;
    }

    return STDL_ERR_OK;
}


int stdl_qexp_first_hyperpolarizability_hrs(float beta[27], float *beta2_ZZZ, float *beta2_ZXX, float *beta2_J1, float *beta2_J3) {
    assert(beta != NULL && beta2_ZZZ != NULL && beta2_ZXX && beta2_J1 != NULL && beta2_J3 != NULL);

    *beta2_ZZZ = .0f;
    *beta2_ZXX = .0f;
    *beta2_J1 = .0f;
    *beta2_J3 = .0f;

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                *beta2_ZZZ += 2 * powf(beta[i*9 + j*3 + k], 2) + 1 * beta[i*9 + j*3 + j] * beta[i*9 + k*3 + k] + 4 * (beta[i*9 + i*3 + j] * beta[j*9 + k*3 + k] + beta[i*9 + i*3 + j] * beta[k*9 + j*3 + k] + beta[i*9 + j*3 + k] * beta[j*9 + i*3 + k]);
                *beta2_ZXX += 6 * powf(beta[i*9 + j*3 + k], 2) + 3 * beta[i*9 + j*3 + j] * beta[i*9 + k*3 + k] - 2 * (beta[i*9 + i*3 + j] * beta[j*9 + k*3 + k] + beta[i*9 + i*3 + j] * beta[k*9 + j*3 + k] + beta[i*9 + j*3 + k] * beta[j*9 + i*3 + k]);
                *beta2_J1 += .6f * beta[i*9 + j*3 + j] * beta[i*9 + k*3 + k];
                *beta2_J3 += powf(beta[i*9 + j*3 + k], 2) - .6f * beta[i*9 + j*3 + j] * beta[i*9 + k*3 + k];
            }
        }
    }

    *beta2_ZZZ /= 105;
    *beta2_ZXX /= 105;

    return STDL_ERR_OK;
}
