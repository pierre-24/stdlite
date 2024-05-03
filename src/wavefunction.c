#include <assert.h>
#include <string.h>

#include "stdlite/wavefunction.h"
#include "stdlite/logging.h"
#include "stdlite/helpers.h"
#include "stdlite/utils/matrix.h"
#include "stdlite/linear_algebra.h"

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
                (STDL_LA_INT) nmo, (STDL_LA_INT) nao,
                1.f, sqrtS, (STDL_LA_INT) nao,
                C, (STDL_LA_INT) nao,
                .0, tmp, (STDL_LA_INT) nao
    );

    memcpy(C, tmp, nmo * nao * sizeof(double));

    STDL_FREE_ALL(sqrtS, tmp);

    return STDL_ERR_OK;
}

int stdl_wavefunction_compute_density_dsp(size_t nocc, size_t nmo, size_t nao, double *C, double *D) {
    assert(C != NULL && D != NULL && nocc > 0 && nmo > 0 && nao > 0 && D != NULL);

    STDL_DEBUG("Compute density");

    #pragma  omp parallel for schedule(guided)
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

int stdl_wavefunction_dsp_ao_to_dsp_mo(size_t nao, size_t nmo, int issym, double *C, double *X_AO, double *X_MO) {
    assert(C != NULL && X_AO != NULL && nmo > 0 && nao > 0 && X_MO != NULL);

    STDL_DEBUG("AO to MO");

    #pragma omp parallel for schedule(guided)
    for (size_t p = 0; p < nmo; ++p) {
        for (size_t q = 0; q <= p; ++q) {
            double sum = 0;
            for (size_t mu = 0; mu < nao; ++mu) {
                for (size_t nu = 0; nu < nao; ++nu) {
                    sum += (!issym && mu < nu ? -1:1) * C[p * nao + mu] * X_AO[STDL_MATRIX_SP_IDX(mu, nu)] * C[q * nao + nu];
                }
            }

            X_MO[STDL_MATRIX_SP_IDX(p, q)] = sum;
        }
    }

    return STDL_ERR_OK;
}

int stdl_wavefunction_dge_ao_to_dge_mo(size_t nao, size_t nmo, double *C, double *X_AO, double *X_MO) {
    assert(C != NULL && X_AO != NULL && nmo > 0 && nao > 0 && X_MO != NULL);

    STDL_DEBUG("AO to MO");

    #pragma omp parallel for
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


int stdl_wavefunction_dump_h5(stdl_wavefunction* wf, hid_t file_id) {
    assert(wf != NULL && file_id != H5I_INVALID_HID);

    hid_t wf_group_id;
    herr_t status;

    wf_group_id = H5Gcreate(file_id, "wavefunction", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    STDL_ERROR_HANDLE_AND_REPORT(wf_group_id == H5I_INVALID_HID, return STDL_ERR_WRITE, "cannot create group");

    // info
    status = H5LTmake_dataset(wf_group_id, "info", 1, (hsize_t[]) {4}, H5T_NATIVE_ULONG, (size_t[]) {wf->natm, wf->nocc, wf->nao, wf->nmo});
    STDL_ERROR_HANDLE_AND_REPORT(status < 0, return STDL_ERR_WRITE, "cannot create dataset in group %d", wf_group_id);

    // atm
    status = H5LTmake_dataset(wf_group_id, "atm", 2, (hsize_t[]) {wf->natm, 4}, H5T_NATIVE_DOUBLE, wf->atm);
    STDL_ERROR_HANDLE_AND_REPORT(status < 0, return STDL_ERR_WRITE, "cannot create dataset in group %d", wf_group_id);

    // ao2atm
    status = H5LTmake_dataset(wf_group_id, "aotoatm", 1, (hsize_t[]) {wf->nao}, H5T_NATIVE_ULONG, wf->aotoatm);
    STDL_ERROR_HANDLE_AND_REPORT(status < 0, return STDL_ERR_WRITE, "cannot create dataset in group %d", wf_group_id);

    // S
    status = H5LTmake_dataset(wf_group_id, "S", 1, (hsize_t[]) {STDL_MATRIX_SP_SIZE(wf->nao)}, H5T_NATIVE_DOUBLE, wf->S);
    STDL_ERROR_HANDLE_AND_REPORT(status < 0, return STDL_ERR_WRITE, "cannot create dataset in group %d", wf_group_id);

    // C
    status = H5LTmake_dataset(wf_group_id, "C", 2, (hsize_t[]) {wf->nmo, wf->nao}, H5T_NATIVE_DOUBLE, wf->C);
    STDL_ERROR_HANDLE_AND_REPORT(status < 0, return STDL_ERR_WRITE, "cannot create dataset in group %d", wf_group_id);

    // e
    status = H5LTmake_dataset(wf_group_id, "e", 1, (hsize_t[]) {wf->nmo}, H5T_NATIVE_DOUBLE, wf->e);
    STDL_ERROR_HANDLE_AND_REPORT(status < 0, return STDL_ERR_WRITE, "cannot create dataset in group %d", wf_group_id);

    // ... and close
    status = H5Gclose(wf_group_id);
    STDL_ERROR_HANDLE_AND_REPORT(status < 0, return STDL_ERR_WRITE, "cannot close group %d", wf_group_id);

    // Set some attributes
    status = H5LTset_attribute_string(file_id, "wavefunction", "type", "stdl_wavefunction");
    STDL_ERROR_HANDLE_AND_REPORT(status < 0, return STDL_ERR_WRITE, "cannot write attribute");

    status = H5LTset_attribute_uint(file_id, "wavefunction", "version", (unsigned int[]) {1, 0}, 2);
    STDL_ERROR_HANDLE_AND_REPORT(status < 0, return STDL_ERR_WRITE, "cannot write attribute");

    return STDL_ERR_OK;
}

int stdl_wavefunction_load_h5(hid_t file_id, stdl_wavefunction **wf_ptr) {
    assert(wf_ptr != NULL && file_id != H5I_INVALID_HID);

    hid_t wf_group_id;
    herr_t status;
    char strbuff[128];
    int version[2];
    size_t ulongbuff[32];
    int err;

    // check attributes
    status = H5LTget_attribute_string(file_id, "wavefunction", "type", strbuff);
    STDL_ERROR_HANDLE_AND_REPORT(status < 0 || strcmp(strbuff, "stdl_wavefunction") != 0, err = STDL_ERR_READ; goto _end, "missing or incorrect attribute for `wavefunction`");

    status = H5LTget_attribute_int(file_id, "wavefunction", "version", version);
    STDL_ERROR_HANDLE_AND_REPORT(status < 0 || version[0] != 1 || version[1] != 0, err = STDL_ERR_READ; goto _end, "missing or incorrect attribute for `wavefunction`");

    // read wavefunction
    wf_group_id = H5Gopen1(file_id, "wavefunction");
    STDL_ERROR_HANDLE_AND_REPORT(wf_group_id == H5I_INVALID_HID, err = STDL_ERR_READ; goto _end, "unable to open group");

    status = H5LTread_dataset(wf_group_id, "info", H5T_NATIVE_ULONG, ulongbuff);
    STDL_ERROR_HANDLE_AND_REPORT(status < 0, err = STDL_ERR_READ; goto _end, "cannot read dataset");

    err = stdl_wavefunction_new(ulongbuff[0], ulongbuff[1], ulongbuff[2], ulongbuff[3], wf_ptr);
    STDL_ERROR_CODE_HANDLE(err, goto _end);

    status = H5LTread_dataset(wf_group_id, "atm", H5T_NATIVE_DOUBLE, (*wf_ptr)->atm);
    STDL_ERROR_HANDLE_AND_REPORT(status < 0, err = STDL_ERR_READ; goto _end, "cannot read dataset");

    status = H5LTread_dataset(wf_group_id, "aotoatm", H5T_NATIVE_ULLONG, (*wf_ptr)->aotoatm);
    STDL_ERROR_HANDLE_AND_REPORT(status < 0, err = STDL_ERR_READ; goto _end, "cannot read dataset");

    status = H5LTread_dataset(wf_group_id, "S", H5T_NATIVE_DOUBLE, (*wf_ptr)->S);
    STDL_ERROR_HANDLE_AND_REPORT(status < 0, err = STDL_ERR_READ; goto _end, "cannot read dataset");

    status = H5LTread_dataset(wf_group_id, "C", H5T_NATIVE_DOUBLE, (*wf_ptr)->C);
    STDL_ERROR_HANDLE_AND_REPORT(status < 0, err = STDL_ERR_READ; goto _end, "cannot read dataset");

    status = H5LTread_dataset(wf_group_id, "e", H5T_NATIVE_DOUBLE, (*wf_ptr)->e);
    STDL_ERROR_HANDLE_AND_REPORT(status < 0, err = STDL_ERR_READ; goto _end, "cannot read dataset");

    status = H5Gclose(wf_group_id);
    STDL_ERROR_HANDLE_AND_REPORT(status < 0, return STDL_ERR_WRITE, "cannot close group %d", wf_group_id);

    _end:
    if(err != STDL_ERR_OK && *wf_ptr != NULL)
        stdl_wavefunction_delete(*wf_ptr);

    return err;
}

int stdl_wavefunction_approximate_size(stdl_wavefunction *wf, size_t *sz) {
    assert(wf != NULL && sz != NULL);

    *sz = sizeof(stdl_wavefunction)
            + (wf->natm * 4 + wf->nao + STDL_MATRIX_SP_SIZE(wf->nao)
            + wf->nao * wf->nmo + wf->nmo) * sizeof(double );

    return STDL_ERR_OK;
}
