#include <assert.h>
#include <lapacke.h>
#include <string.h>

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
