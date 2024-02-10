#include <assert.h>
#include <cblas.h>
#include <string.h>

#include "stdlite/wavefunction.h"
#include "stdlite/logging.h"
#include "stdlite/helpers.h"
#include "stdlite/utils/matrix.h"

int stdl_wavefunction_new(size_t natm, size_t nocc, size_t nao, size_t nmo, stdl_wavefunction **wf_ptr) {
    assert(wf_ptr != NULL && natm > 0 && nao > 0 && nmo > 0 && nmo <= nao && nmo >= nocc);

    *wf_ptr = malloc(sizeof(stdl_wavefunction));
    STDL_ERROR_HANDLE_AND_REPORT(*wf_ptr == NULL, return STDL_ERR_MALLOC, "malloc");

    (*wf_ptr)->natm = natm;
    (*wf_ptr)->nao = nao;
    (*wf_ptr)->nmo = nmo;
    (*wf_ptr)->nocc = nocc;

    (*wf_ptr)->atm = (*wf_ptr)->S = (*wf_ptr)->C = (*wf_ptr)->e = NULL;
    (*wf_ptr)->aotoatm = NULL;

    (*wf_ptr)->atm = malloc(4 * natm * sizeof(double));
    STDL_ERROR_HANDLE_AND_REPORT((*wf_ptr)->atm == NULL, stdl_wavefunction_delete(*wf_ptr); return STDL_ERR_MALLOC, "malloc");

    (*wf_ptr)->e = malloc(nao * sizeof(double));
    STDL_ERROR_HANDLE_AND_REPORT((*wf_ptr)->e == NULL, stdl_wavefunction_delete(*wf_ptr); return STDL_ERR_MALLOC, "malloc");

    (*wf_ptr)->aotoatm = malloc(nao * sizeof(size_t));
    STDL_ERROR_HANDLE_AND_REPORT((*wf_ptr)->aotoatm == NULL, stdl_wavefunction_delete(*wf_ptr); return STDL_ERR_MALLOC, "malloc");

    (*wf_ptr)->S = malloc(STDL_MATRIX_SP_SIZE(nao) * sizeof(double));
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

int stdl_wavefunction_orthogonalize_C_dge(size_t nmo, size_t nao, double *S, double *C) {
    assert(C != NULL && S != NULL && nmo > 0 && nao > 0);

    // compute S^1/2
    double* sqrtS = malloc(nao * nao * sizeof(double));
    STDL_ERROR_HANDLE_AND_REPORT(sqrtS == NULL, return STDL_ERR_MALLOC, "malloc");

    int error = stdl_matrix_dsp_sqrt_sy(nao, S, sqrtS);
    STDL_ERROR_CODE_HANDLE(error, free(sqrtS); return error);

    // C' = C * S^1/2
    double* tmp = malloc(nmo * nao * sizeof(double ));
    STDL_ERROR_HANDLE_AND_REPORT(tmp == NULL, STDL_FREE_ALL(sqrtS); return STDL_ERR_MALLOC, "malloc");

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

int stdl_wavefunction_compute_density_dsp(size_t nocc, size_t nmo, size_t nao, double *C, double *D) {
    assert(C != NULL && D != NULL && nocc > 0 && nmo > 0 && nao > 0 && D != NULL);

    for (size_t mu = 0; mu < nao; ++mu) {
        for (size_t nu = 0; nu <= mu; ++nu) {
            double sum = .0f;
            for (size_t p = 0; p < nmo; ++p) {
                sum += C[p * nao + mu] * ((p < nocc) ? 2 : 0) *  C[p * nao + nu] ;
            }
            D[STDL_MATRIX_SP_IDX(mu, nu)] = sum;
        }
    }

    return STDL_ERR_OK;
}

int stdl_wavefunction_dsy_ao_to_mo(size_t nao, size_t nmo, double* C, double* X_AO, double* X_MO) {
    double* X_tmp = malloc(nao * nmo * sizeof(double));
    STDL_ERROR_HANDLE_AND_REPORT(X_tmp == NULL, return STDL_ERR_MALLOC, "malloc");

    // X_tmp = C * X_AO
    cblas_dsymm(CblasRowMajor, CblasRight, CblasLower,
                (int) nmo, (int) nao,
                1.f, X_AO, (int) nao,
                C, (int) nao,
                .0, X_tmp, (int) nao
    );

    // X_MO = X * C^T
    cblas_dgemm(
            CblasRowMajor, CblasNoTrans, CblasTrans,
            (int) nmo, (int) nmo, (int) nao,
            1.f, X_tmp, (int) nao,
            C, (int) nao,
            .0, X_MO, (int) nmo
    );

    STDL_FREE_ALL(X_tmp);

    return STDL_ERR_OK;
}


