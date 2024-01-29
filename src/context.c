#include <assert.h>
#include <string.h>
#include <cblas.h>

#include "stdlite/context.h"
#include "stdlite/logging.h"
#include "stdlite/helpers.h"
#include "stdlite/matrix.h"

int stdlite_context_new(stdl_context ** ctx, stdl_wavefunction* wf, stdl_basis* bs, float gammaJ, float  gammaK, float ethr, float ax) {
    assert(ctx != NULL && wf != NULL && bs != NULL && gammaJ > 0 && gammaK > 0 && ethr > 0 && ax >= 0 && ax <= 1);

    *ctx = malloc(sizeof(stdl_context));
    STDL_ERROR_HANDLE_AND_REPORT(*ctx == NULL, return STDL_ERR_MALLOC, "malloc");

    (*ctx)->wf = wf;
    (*ctx)->bs = bs;
    (*ctx)->gammaJ = gammaJ;
    (*ctx)->gammaK = gammaK;
    (*ctx)->ethr = ethr;
    (*ctx)->ax = ax;

    // select MO to include
    STDL_DEBUG("range: %f a.u. (%.3f eV)", ethr, ethr * 27.212);

    size_t ohomo = (int) wf->nelec / 2 - 1, omin = 0, omax = 0;
    double ehomo = wf->e[ohomo], elumo = wf->e[ohomo + 1], emin = elumo - 2 * (1 + .8 * ax) * ethr, emax = ehomo+ 2 * (1 + .8 * ax) * ethr;

    STDL_DEBUG("occ MO cutoff: %f a.u. (%.3feV)", emin, emin * 27.212);
    STDL_DEBUG("virt MO cutoff: %f a.u. (%.3feV)", emax, emax * 27.212);

    for(size_t i=0; i < wf->nmo; i++) {
        if(wf->e[i] >= emin && omin == 0)
            omin = (int) i;

        if(wf->e[i] <= emax)
            omax = (int) i;
        else
            break;
    }

    (*ctx)->nmo = omax - omin + 1;
    (*ctx)->nocc = ohomo - omin;
    (*ctx)->nvirt = omax - ohomo + 1;

    STDL_DEBUG("Resulting partition: [%d || %d | %d || %d] (occ + virt = %d)", omin, (*ctx)->nocc, (*ctx)->nvirt, wf->nmo - omax - 1, (*ctx)->nmo);

    (*ctx)->e = malloc((*ctx)->nmo * sizeof(double ));
    (*ctx)->C = malloc((*ctx)->nmo * wf->nao * sizeof(double));

    STDL_ERROR_HANDLE_AND_REPORT((*ctx)->e == NULL || (*ctx)->C == NULL, stdlite_context_delete(*ctx); return STDL_ERR_MALLOC, "malloc");

    // copy coefficients
    for(size_t i=0; i < (*ctx)->nmo; i++) {
        (*ctx)->e[i] = wf->e[omin + i];
        memcpy(&((*ctx)->C[i * wf->nao]), &(wf->C[(i + omin) * wf->nao]), wf->nao * sizeof(double));
    }

    STDL_DEBUG("Compute S^(1/2)");

    // compute S^(1/2)
    double* sqrtS = malloc(wf->nao * wf->nao * sizeof(double));
    STDL_ERROR_HANDLE_AND_REPORT(sqrtS == NULL, return STDL_ERR_MALLOC, "malloc");
    memcpy(sqrtS, wf->S, wf->nao * wf->nao * sizeof(double));

    int error = stdl_matrix_dge_sqrt(&sqrtS, wf->nao);
    STDL_ERROR_CODE_HANDLE(error, return  error);

    STDL_DEBUG("Orthogonalize MOs");

    // C' = C * S^1/2 (the side depends on if first index is MO or AO)
    double* tmp = malloc((*ctx)->nmo * wf->nao * sizeof(double));
    STDL_ERROR_HANDLE_AND_REPORT(tmp == NULL, return STDL_ERR_MALLOC, "malloc");

    cblas_dsymm(CblasRowMajor, CblasRight, CblasLower,
                (int) (*ctx)->nmo, (int) wf->nao,
                1.f, sqrtS, (int) wf->nao,
                (*ctx)->C, (int) wf->nao,
                .0, tmp, (int) wf->nao
    );

    memcpy((*ctx)->C, tmp, (*ctx)->nmo * wf->nao * sizeof(double));

    STDL_FREE_ALL(sqrtS, tmp);

    return STDL_ERR_OK;
}

int stdlite_context_delete(stdl_context* ctx) {
    assert(ctx != NULL);

    if(ctx->wf != NULL)
        stdl_wavefunction_delete(ctx->wf);

    if(ctx->bs != NULL)
        stdl_basis_delete(ctx->bs);

    STDL_FREE_ALL(ctx->e, ctx->C, ctx);

    return STDL_ERR_OK;
}
