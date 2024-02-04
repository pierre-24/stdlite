#include <assert.h>
#include <lapacke.h>
#include <string.h>

#include "stdlite/response.h"
#include "stdlite/logging.h"
#include "stdlite/helpers.h"


int stdl_response_casida_TDA_full(stdl_context *ctx, size_t ncsfs, float *A, float **energies, float **amplitudes) {
    assert(ctx != NULL && ncsfs > 0 && A != NULL && energies != NULL && amplitudes != NULL);

    *energies = malloc(ncsfs * sizeof(float ));
    *amplitudes = malloc(ncsfs * ncsfs * sizeof(float ));

    memcpy(*amplitudes, A, ncsfs * ncsfs * sizeof(float ));

    STDL_ERROR_HANDLE_AND_REPORT(*energies == NULL || *amplitudes == NULL, STDL_FREE_ALL(*energies, *amplitudes); return STDL_ERR_MALLOC, "malloc");

    int err = LAPACKE_ssyev(LAPACK_ROW_MAJOR, 'V', 'L', (int) ncsfs, *amplitudes, (int) ncsfs, *energies);
    STDL_ERROR_HANDLE_AND_REPORT(err != 0, STDL_FREE_ALL(*energies, *amplitudes); return STDL_ERR_MALLOC, "error while ssyev: %d", err);

    return STDL_ERR_OK;
}
