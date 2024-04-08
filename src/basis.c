#include <assert.h>
#include <stdio.h>
#include <string.h>
#include "stdlite/basis.h"
#include "stdlite/logging.h"
#include "stdlite.h"


int stdl_basis_new(int natm, int nbas, size_t env_size, int use_spherical, stdl_basis **bs_ptr) {
    assert(bs_ptr != NULL && natm > 0 && nbas > 0 && env_size > (3 * (size_t) natm + PTR_ENV_START));

    *bs_ptr = malloc(sizeof(stdl_basis));
    STDL_ERROR_HANDLE_AND_REPORT(*bs_ptr == NULL, return STDL_ERR_MALLOC, "malloc");

    STDL_DEBUG("create basis %p", *bs_ptr);

    (*bs_ptr)->natm = natm;
    (*bs_ptr)->nbas = nbas;
    (*bs_ptr)->use_spherical = use_spherical;
    (*bs_ptr)->env_size = env_size;

    (*bs_ptr)->atm = (*bs_ptr)->bas = NULL;
    (*bs_ptr)->env = NULL;

    (*bs_ptr)->atm = calloc(6 * natm, sizeof(int));
    (*bs_ptr)->bas = calloc(8 * nbas, sizeof(int));
    STDL_ERROR_HANDLE_AND_REPORT((*bs_ptr)->atm == NULL || (*bs_ptr)->bas == NULL, stdl_basis_delete(*bs_ptr); return STDL_ERR_MALLOC, "malloc");

    (*bs_ptr)->env = malloc(env_size * sizeof(double));
    STDL_ERROR_HANDLE_AND_REPORT((*bs_ptr)->env == NULL, stdl_basis_delete(*bs_ptr); return STDL_ERR_MALLOC, "malloc");

    // fill the first 20 elements of env with zeroes
    for (int i = 0; i < PTR_ENV_START; ++i) {
        (*bs_ptr)->env[i] = .0;
    }

    return STDL_ERR_OK;
}

int stdl_basis_delete(stdl_basis *bs) {
    assert(bs != NULL);

    STDL_DEBUG("delete basis %p", bs);

    STDL_FREE_ALL(bs->atm, bs->bas, bs->env, bs);

    return STDL_ERR_OK;
}

char _toc(int i) {
    switch (i) {
        case 0:
            return 's';
        case 1:
            return 'p';
        case 2:
            return 'd';
        case 3:
            return 'f';
        case 4:
            return 'g';
        case 5:
            return 'h';
        case 6:
            return 'i';
        default:
            return '?';
    }
}

int stdl_basis_reorder_C(size_t nmo, size_t nao, double *C, stdl_basis *bs, size_t maxshell, int **transpose) {
    assert(nao > 0 && nmo <= nao && C != NULL && bs != NULL && maxshell > 0 && transpose != NULL);

    STDL_DEBUG("reorder orbitals");

    double buff[CART_MAX] = {0}; // maximum that can be handled at the moment

    int sj, joffset;

    for (size_t i = 0; i < nmo; ++i) {

        joffset = 0;

        for(int j=0; j < bs->nbas; j++) {
            int angmom = bs->bas[j * 8 + 1];
            STDL_ERROR_HANDLE_AND_REPORT((size_t) angmom > maxshell, return STDL_ERR_BASIS, "Angular momentum %d is larger than possible (%ld)", angmom, maxshell);

            if (bs->use_spherical)
                sj = CINTcgto_spheric(j, bs->bas);
            else
                sj = CINTcgtos_cart(j, bs->bas);

            if(angmom > 1) { // s & p functions are ok
                STDL_DEBUG("reordering basis function #%d", j);
                // reorder
                for(int mu = 0; mu < sj; mu++) {
                    buff[transpose[angmom][mu]] = C[i * nao + joffset + mu];
                }

                // and copy back
                memcpy(C + i * nao + joffset, buff, sj * sizeof(double ));
            }

            joffset += sj;
        }
    }

    return STDL_ERR_OK;
}

int stdl_basis_dump_h5(stdl_basis *bs, hid_t file_id) {
    assert(file_id != H5I_INVALID_HID && bs != NULL);

    hid_t bs_group_id;
    herr_t status;

    bs_group_id = H5Gcreate(file_id, "basis", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    STDL_ERROR_HANDLE_AND_REPORT(bs_group_id == H5I_INVALID_HID, return STDL_ERR_WRITE, "cannot create group");

    // info
    status = H5LTmake_dataset(bs_group_id, "info", 1, (hsize_t[]) {4}, H5T_NATIVE_ULONG, (size_t[]) {(size_t) bs->natm, (size_t) bs->nbas, bs->env_size,  (size_t) bs->use_spherical != 0});
    STDL_ERROR_HANDLE_AND_REPORT(status < 0, return STDL_ERR_WRITE, "cannot create dataset in group %d", bs_group_id);

    // atm
    status = H5LTmake_dataset(bs_group_id, "atm", 2, (hsize_t[]) {bs->natm, 6}, H5T_NATIVE_INT, bs->atm);
    STDL_ERROR_HANDLE_AND_REPORT(status < 0, return STDL_ERR_WRITE, "cannot create dataset in group %d", bs_group_id);

    // bas
    status = H5LTmake_dataset(bs_group_id, "bas", 2, (hsize_t[]) {bs->nbas, 8}, H5T_NATIVE_INT, bs->bas);
    STDL_ERROR_HANDLE_AND_REPORT(status < 0, return STDL_ERR_WRITE, "cannot create dataset in group %d", bs_group_id);

    // env
    status = H5LTmake_dataset(bs_group_id, "env", 1, (hsize_t[]) {bs->env_size}, H5T_NATIVE_DOUBLE, bs->env);
    STDL_ERROR_HANDLE_AND_REPORT(status < 0, return STDL_ERR_WRITE, "cannot create dataset in group %d", bs_group_id);

    // ... and close
    status = H5Gclose(bs_group_id);
    STDL_ERROR_HANDLE_AND_REPORT(status < 0, return STDL_ERR_WRITE, "cannot close group %d", bs_group_id);

    // set attributes
    status = H5LTset_attribute_string(file_id, "basis", "type", "stdl_basis");
    STDL_ERROR_HANDLE_AND_REPORT(status < 0, return STDL_ERR_WRITE, "cannot write attribute");

    status = H5LTset_attribute_uint(file_id, "basis", "version", (unsigned int[]) {1, 0}, 2);
    STDL_ERROR_HANDLE_AND_REPORT(status < 0, return STDL_ERR_WRITE, "cannot write attribute");

    return STDL_ERR_OK;
}

int stdl_basis_load_h5(hid_t file_id, stdl_basis **bs_ptr) {
    assert(file_id != H5I_INVALID_HID && bs_ptr != NULL);

    hid_t bs_group_id;
    herr_t status;
    int err = STDL_ERR_OK;
    long ulongbuff[32];
    char strbuff[128];
    int version[2];

    // check attributes
    status = H5LTget_attribute_string(file_id, "basis", "type", strbuff);
    STDL_ERROR_HANDLE_AND_REPORT(status < 0 || strcmp(strbuff, "stdl_basis") != 0, err = STDL_ERR_READ; goto _end, "missing or incorrect attribute for `basis`");

    status = H5LTget_attribute_int(file_id, "basis", "version", version);
    STDL_ERROR_HANDLE_AND_REPORT(status < 0 || version[0] != 1 || version[1] != 0, err = STDL_ERR_READ; goto _end, "missing or incorrect attribute for `basis`");

    // read basis
    bs_group_id = H5Gopen1(file_id, "basis");
    STDL_ERROR_HANDLE_AND_REPORT(bs_group_id == H5I_INVALID_HID, err = STDL_ERR_READ; goto _end, "unable to open group");

    status = H5LTread_dataset(bs_group_id, "info", H5T_NATIVE_ULONG, ulongbuff);
    STDL_ERROR_HANDLE_AND_REPORT(status < 0, err = STDL_ERR_READ; goto _end, "cannot read dataset");

    err = stdl_basis_new((int) ulongbuff[0], (int) ulongbuff[1], ulongbuff[2], ulongbuff[3] != 0, bs_ptr);
    STDL_ERROR_CODE_HANDLE(err, goto _end);

    status = H5LTread_dataset(bs_group_id, "atm", H5T_NATIVE_INT, (*bs_ptr)->atm);
    STDL_ERROR_HANDLE_AND_REPORT(status < 0, err = STDL_ERR_READ; goto _end, "cannot read dataset");

    status = H5LTread_dataset(bs_group_id, "bas", H5T_NATIVE_INT, (*bs_ptr)->bas);
    STDL_ERROR_HANDLE_AND_REPORT(status < 0, err = STDL_ERR_READ; goto _end, "cannot read dataset");

    status = H5LTread_dataset(bs_group_id, "env", H5T_NATIVE_DOUBLE, (*bs_ptr)->env);
    STDL_ERROR_HANDLE_AND_REPORT(status < 0, err = STDL_ERR_READ; goto _end, "cannot read dataset");

    status = H5Gclose(bs_group_id);
    STDL_ERROR_HANDLE_AND_REPORT(status < 0, return STDL_ERR_WRITE, "cannot close group %d",bs_group_id);

    _end:
    if(err != STDL_ERR_OK && *bs_ptr != NULL)
        stdl_basis_delete(*bs_ptr);

    return err;
}

int stdl_basis_approximate_size(stdl_basis *bs, size_t *sz) {
    assert(bs != NULL && sz != NULL);

    *sz = sizeof(stdl_basis)
            + (bs->natm * 6 + bs->natm * 8) * sizeof(int)
            + bs->env_size * sizeof(double);

    return STDL_ERR_OK;
}

