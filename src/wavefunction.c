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

    STDL_DEBUG("create wavefunction %p", *wf_ptr);

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

    STDL_DEBUG("delete wavefunction %p", wf);

    STDL_FREE_ALL(wf->atm, wf->e, wf->aotoatm, wf->S, wf->C, wf);

    return STDL_ERR_OK;
}

int stdl_wavefunction_orthogonalize_C_dge(size_t nmo, size_t nao, double *S, double *C) {
    assert(C != NULL && S != NULL && nmo > 0 && nao > 0);

    STDL_DEBUG("Orthogonalize LCAO coefs");

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

    STDL_DEBUG("Compute density");

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

int stdl_wavefunction_dsp_ao_to_dsp_mo(size_t nao, size_t nmo, double* C, double* X_AO, double* X_MO) {
    assert(C != NULL && X_AO != NULL && nmo > 0 && nao > 0 && X_MO != NULL);

    for (size_t p = 0; p < nmo; ++p) {
        for (size_t q = 0; q <= p; ++q) {
            double sum = 0;
            for (size_t mu = 0; mu < nao; ++mu) {
                for (size_t nu = 0; nu < nao; ++nu) {
                    sum += C[p * nao + mu] * X_AO[STDL_MATRIX_SP_IDX(mu, nu)] * C[q * nao + nu];
                }
            }

            X_MO[STDL_MATRIX_SP_IDX(p, q)] = sum;
        }
    }

    return STDL_ERR_OK;
}

int stdl_wavefunction_dge_ao_to_dge_mo(size_t nao, size_t nmo, double *C, double *X_AO, double *X_MO) {
    assert(C != NULL && X_AO != NULL && nmo > 0 && nao > 0 && X_MO != NULL);

    for (size_t p = 0; p < nmo; ++p) {
        for (size_t q = 0; q < nmo; ++q) {
            double sum = 0;
            for (size_t mu = 0; mu < nao; ++mu) {
                for (size_t nu = 0; nu < nao; ++nu) {
                    sum += C[p * nao + mu] * X_AO[mu * nao + nu] * C[q * nao + nu];
                }
            }

            X_MO[p * nmo + q] = sum;
        }
    }

    return STDL_ERR_OK;
}


