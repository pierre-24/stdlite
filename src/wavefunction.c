#include <assert.h>
#include <math.h>
#include <lapacke.h>
#include <cblas.h>
#include <string.h>

#include "stdlite/wavefunction.h"
#include "stdlite/logging.h"
#include "stdlite/helpers.h"
#include "stdlite/matrix.h"

int stdl_wavefunction_new(stdl_wavefunction **wf_ptr, size_t natm, size_t nelec, size_t nao, size_t nmo) {
    assert(wf_ptr != NULL && natm > 0 && nao > 0 && nmo > 0 && nmo <= nao && 2 * nmo >= nelec && nelec % 2 == 0);

    *wf_ptr = malloc(sizeof(stdl_wavefunction));
    STDL_ERROR_HANDLE_AND_REPORT(*wf_ptr == NULL, return STDL_ERR_MALLOC, "malloc");

    (*wf_ptr)->isortho = 0;
    (*wf_ptr)->natm = natm;
    (*wf_ptr)->nao = nao;
    (*wf_ptr)->nmo = nmo;
    (*wf_ptr)->nelec = nelec;

    (*wf_ptr)->atm = (*wf_ptr)->S = (*wf_ptr)->C = (*wf_ptr)->e = NULL;
    (*wf_ptr)->aotoatm = NULL;

    (*wf_ptr)->atm = malloc(4 * natm * sizeof(double));
    STDL_ERROR_HANDLE_AND_REPORT((*wf_ptr)->atm == NULL, stdl_wavefunction_delete(*wf_ptr); return STDL_ERR_MALLOC, "malloc");

    (*wf_ptr)->e = malloc(nao * sizeof(double));
    STDL_ERROR_HANDLE_AND_REPORT((*wf_ptr)->e == NULL, stdl_wavefunction_delete(*wf_ptr); return STDL_ERR_MALLOC, "malloc");

    (*wf_ptr)->aotoatm = malloc(nao * sizeof(size_t));
    STDL_ERROR_HANDLE_AND_REPORT((*wf_ptr)->aotoatm == NULL, stdl_wavefunction_delete(*wf_ptr); return STDL_ERR_MALLOC, "malloc");

    (*wf_ptr)->S = malloc(nao * nao * sizeof(double));
    STDL_ERROR_HANDLE_AND_REPORT((*wf_ptr)->S == NULL, stdl_wavefunction_delete(*wf_ptr); return STDL_ERR_MALLOC, "malloc");

    (*wf_ptr)->C = malloc(nao * nmo * sizeof(double));
    STDL_ERROR_HANDLE_AND_REPORT((*wf_ptr)->C == NULL, stdl_wavefunction_delete(*wf_ptr); return STDL_ERR_MALLOC, "malloc");

    return STDL_ERR_OK;
}

int stdl_wavefunction_delete(stdl_wavefunction *wf) {
    assert(wf != NULL);

    STDL_FREE_ALL(wf->atm, wf->e, wf->aotoatm, wf->S, wf->C, wf);

    return STDL_ERR_OK;
}

int stdl_wavefunction_orthogonalize(stdl_wavefunction *wf) {
    assert(wf != NULL);

    if(wf->isortho) {
        STDL_DEBUG("wavefunction is already orthogonal");
        return STDL_ERR_OK;
    }

    int error = stdl_matrix_ge_sqrt(&(wf->S), wf->nao);
    STDL_ERROR_CODE_HANDLE(error, return  error);

    // C' = C * S^1/2 (the side depends on if first index is MO or AO)
    double* tmp = malloc(wf->nmo * wf->nao * sizeof(double ));
    STDL_ERROR_HANDLE_AND_REPORT(tmp == NULL, return STDL_ERR_MALLOC, "malloc");

    cblas_dsymm(CblasRowMajor, CblasRight, CblasLower,
                (int) wf->nmo, (int) wf->nao,
                1.f, wf->S, (int) wf->nao,
                wf->C, (int) wf->nao,
                .0, tmp, (int) wf->nao
    );

    /* Or:
     * cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                (int) wf->nmo, (int) wf->nao, (int) wf->nao,
                1.f, wf->C, (int) wf->nao,
                wf->S, (int) wf->nao,
                .0f, w, (int) wf->nao
    );*/

    memcpy(wf->C, tmp, wf->nmo * wf->nao * sizeof(double));

    STDL_FREE_ALL(tmp);

    // identity S!
    for(size_t i = 0; i < wf->nao; i++) {
        for(size_t j=0; j < wf->nao; j++)
            wf->S[i * wf->nao + j] = i == j ? 1. : 0.;
    }

    wf->isortho = 1;

    return STDL_ERR_OK;
}

int stdl_wavefunction_compute_density(stdl_wavefunction *wf, double **D) {
    assert(wf != NULL && D != NULL);

    *D = malloc(wf->nao * wf->nao * sizeof(double));
    STDL_ERROR_HANDLE_AND_REPORT(*D == NULL, return STDL_ERR_MALLOC, "malloc");

    // X_ik = n_k*C_ik
    double* X = malloc(wf->nao * wf->nmo * sizeof(double));
    STDL_ERROR_HANDLE_AND_REPORT(X == NULL, return STDL_ERR_MALLOC, "malloc");

    for (size_t i = 0; i < wf->nmo; ++i) {
        for (size_t j = 0; j < wf->nao; ++j) {
            X[i * wf->nao + j] = wf->C[i * wf->nao + j] * ((i < wf->nelec / 2) ? 2 : 0);
        }
    }

    // D = X^T * C
    cblas_dgemm(
            CblasRowMajor, CblasTrans, CblasNoTrans,
            (int) wf->nao, (int) wf->nao, (int) wf->nmo,
            1.f, X, (int) wf->nao,
            wf->C, (int) wf->nao,
            .0, *D, (int) wf->nao
    );

    free(X);

    return STDL_ERR_OK;
}
