#include <assert.h>
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

int stdl_wavefunction_orthogonalize_C(double* C, double* S, size_t nmo, size_t nao) {
    assert(C != NULL && S != NULL && nmo > 0 && nao > 0);

    // compute S^(1/2)
    double* sqrtS = malloc(nao * nao * sizeof(double));
    STDL_ERROR_HANDLE_AND_REPORT(sqrtS == NULL, return STDL_ERR_MALLOC, "malloc");
    memcpy(sqrtS, S, nao * nao * sizeof(double));

    int error = stdl_matrix_dge_sqrt(&sqrtS, nao);
    STDL_ERROR_CODE_HANDLE(error, return  error);

    // C' = C * S^1/2
    double* tmp = malloc(nmo * nao * sizeof(double ));
    STDL_ERROR_HANDLE_AND_REPORT(tmp == NULL, return STDL_ERR_MALLOC, "malloc");

    cblas_dsymm(CblasRowMajor, CblasRight, CblasLower,
                (int) nmo, (int) nao,
                1.f, sqrtS, (int) nao,
                C, (int) nao,
                .0, tmp, (int) nao
    );

    memcpy(C, tmp, nmo * nao * sizeof(double));

    STDL_FREE_ALL(sqrtS, tmp);

    return STDL_ERR_OK;
}

int stdl_wavefunction_compute_density(double **D, double* C, size_t nelec, size_t nmo, size_t nao) {
    assert(C != NULL && D != NULL && nelec > 0 && nmo > 0 && nao > 0);

    *D = malloc(nao * nao * sizeof(double));
    STDL_ERROR_HANDLE_AND_REPORT(*D == NULL, return STDL_ERR_MALLOC, "malloc");

    // X_ik = n_k*C_ik
    double* X = malloc(nmo * nao * sizeof(double));
    STDL_ERROR_HANDLE_AND_REPORT(X == NULL, return STDL_ERR_MALLOC, "malloc");

    for (size_t i = 0; i < nmo; ++i) {
        for (size_t j = 0; j < nao; ++j) {
            X[i * nao + j] = C[i * nao + j] * ((i < nelec / 2) ? 2 : 0);
        }
    }

    // D = X^T * C
    cblas_dgemm(
            CblasRowMajor, CblasTrans, CblasNoTrans,
            (int) nao, (int) nao, (int) nmo,
            1.f, X, (int) nao,
            C, (int) nao,
            .0, *D, (int) nao
    );

    free(X);

    return STDL_ERR_OK;
}
