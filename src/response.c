#include <assert.h>
#include <lapacke.h>
#include <string.h>

#include "stdlite/response.h"
#include "stdlite/logging.h"
#include "stdlite/helpers.h"
#include "stdlite/utils/matrix.h"


int stdl_response_casida_TDA_full(stdl_context *ctx, size_t ncsfs, float *A, float **energies, float **amplitudes) {
    assert(ctx != NULL && ncsfs > 0 && A != NULL && energies != NULL && amplitudes != NULL);

    float* tmp = malloc(STDL_MATRIX_SP_SIZE(ncsfs) * sizeof(float ));
    *energies = malloc(ncsfs * sizeof(float ));
    *amplitudes = malloc(ncsfs * ncsfs * sizeof(float ));

    STDL_ERROR_HANDLE_AND_REPORT(*energies == NULL || *amplitudes == NULL || tmp == NULL, STDL_FREE_ALL(*energies, *amplitudes); return STDL_ERR_MALLOC, "malloc");

    memcpy(tmp, A, STDL_MATRIX_SP_SIZE(ncsfs) * sizeof(float ));

    int err = LAPACKE_sspev(
            LAPACK_ROW_MAJOR, 'V', 'L',
            (int) ncsfs, tmp,
            *energies, *amplitudes, (int) ncsfs
            );
    STDL_ERROR_HANDLE_AND_REPORT(err != 0, STDL_FREE_ALL(*energies, *amplitudes, tmp); return STDL_ERR_MALLOC, "error while sspev(): %d", err);

    STDL_FREE_ALL(tmp);

    return STDL_ERR_OK;
}

int stdl_response_casida_TDA(stdl_context* ctx, size_t ncsfs, float *A, size_t nexci, float** energies, float** amplitudes){
    assert(ctx != NULL && ncsfs > 0 && A != NULL && energies != NULL && amplitudes != NULL);

    float* tmp = malloc(STDL_MATRIX_SP_SIZE(ncsfs) * sizeof(float ));
    *energies = malloc(nexci * sizeof(float ));
    *amplitudes = malloc(nexci * ncsfs * sizeof(float ));

    int* ifail = malloc(ncsfs * sizeof(float ));

    STDL_ERROR_HANDLE_AND_REPORT(*energies == NULL || *amplitudes == NULL || tmp == NULL || ifail == NULL, STDL_FREE_ALL(*energies, *amplitudes, tmp, ifail); return STDL_ERR_MALLOC, "malloc");

    memcpy(tmp, A, STDL_MATRIX_SP_SIZE(ncsfs) * sizeof(float ));

    int found = 0;

    int err = LAPACKE_sspevx(
            LAPACK_ROW_MAJOR, 'V', 'I', 'L',
            (int) ncsfs, tmp,
            .0f, .0f,
            1 /* even though we are in C, it starts at 1 */, (int) nexci, STDL_RESPONSE_EIGV_ABSTOL,
            &found,*energies, *amplitudes, (int) nexci, ifail
            );

    STDL_ERROR_HANDLE_AND_REPORT(err != 0, STDL_FREE_ALL(*energies, *amplitudes, tmp, ifail); return STDL_ERR_MALLOC, "error while sspevx(): %d", err);

    STDL_FREE_ALL(tmp, ifail);

    return STDL_ERR_OK;
}
