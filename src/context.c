#include <assert.h>
#include <string.h>

#include "stdlite/context.h"
#include "stdlite/logging.h"
#include "stdlite/helpers.h"
#include "stdlite/matrix.h"

int stdlite_context_new(stdl_context ** ctx, stdl_wavefunction* wf, stdl_basis* bs, float gammaJ, float  gammaK, float ethr, float ax) {
    assert(ctx != NULL && wf != NULL && bs != NULL && gammaJ > 0 && gammaK > 0 && ethr > 0 && ax >= 0 && ax <= 1);

    *ctx = malloc(sizeof(stdl_context));
    STDL_ERROR_HANDLE_AND_REPORT(*ctx == NULL, return STDL_ERR_MALLOC, "malloc");

    (*ctx)->original_wf = wf;
    (*ctx)->bs = bs;
    (*ctx)->gammaJ = gammaJ;
    (*ctx)->gammaK = gammaK;
    (*ctx)->ethr = ethr;
    (*ctx)->ax = ax;

    // select MO to include
    STDL_DEBUG("range: %f a.u. (%.3f eV)", ethr, ethr * 27.212);

    size_t ohomo = (int) wf->nocc - 1, omin = 0, omax = 0;
    double ehomo = wf->e[ohomo], elumo = wf->e[ohomo + 1], ewin = 2 * (1 + .8 * ax)  * ethr, emin = elumo -ewin, emax = ehomo+ ewin;

    STDL_DEBUG("window: %f a.u. (%.3feV)", ewin, ewin * 27.212);
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
    (*ctx)->nocc = ohomo - omin + 1;
    size_t nvirt = omax - ohomo;

    STDL_DEBUG("Resulting partition: [%d || %d | %d || %d] (occ + virt = %d, %.2f%% of %d MOs)", omin, (*ctx)->nocc, nvirt, wf->nmo - omax - 1, (*ctx)->nmo, (double) (*ctx)->nmo / wf->nmo * 100, wf->nmo);

    (*ctx)->e = malloc((*ctx)->nmo * sizeof(double ));
    (*ctx)->C = malloc((*ctx)->nmo * wf->nao * sizeof(double));

    STDL_ERROR_HANDLE_AND_REPORT((*ctx)->e == NULL || (*ctx)->C == NULL, stdlite_context_delete(*ctx); return STDL_ERR_MALLOC, "malloc");

    // copy coefficients and energies
    for(size_t i=0; i < (*ctx)->nmo; i++) {
        (*ctx)->e[i] = wf->e[omin + i];
        memcpy(&((*ctx)->C[i * wf->nao]), &(wf->C[(i + omin) * wf->nao]), wf->nao * sizeof(double));
    }

    STDL_DEBUG("Orthogonalize MOs");

    int error = stdl_wavefunction_orthogonalize_C((*ctx)->C, wf->S, (*ctx)->nmo, wf->nao);
    STDL_ERROR_CODE_HANDLE(error, stdlite_context_delete(*ctx); return error);

    return STDL_ERR_OK;
}

int stdlite_context_delete(stdl_context* ctx) {
    assert(ctx != NULL);

    if(ctx->original_wf != NULL)
        stdl_wavefunction_delete(ctx->original_wf);

    if(ctx->bs != NULL)
        stdl_basis_delete(ctx->bs);

    STDL_FREE_ALL(ctx->e, ctx->C, ctx);

    return STDL_ERR_OK;
}
