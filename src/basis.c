#include <assert.h>
#include <stdio.h>
#include "stdlite/basis.h"
#include "stdlite/logging.h"
#include "stdlite.h"


int stdl_basis_new(stdl_basis **bs_ptr, int natm, int nbas, size_t env_size, int use_spherical) {
    assert(bs_ptr != NULL && natm > 0 && nbas > 0 && env_size > 0);

    *bs_ptr = malloc(sizeof(stdl_basis));
    STDL_ERROR_HANDLE_AND_REPORT(*bs_ptr == NULL, return STDL_ERR_MALLOC, "malloc");

    (*bs_ptr)->natm = natm;
    (*bs_ptr)->nbas = nbas;
    (*bs_ptr)->use_spherical = use_spherical;

    (*bs_ptr)->atm = (*bs_ptr)->bas = NULL;
    (*bs_ptr)->env = NULL;

    (*bs_ptr)->atm = malloc(6 * natm * sizeof(int));
    (*bs_ptr)->bas = malloc(8 * nbas * sizeof(int));
    STDL_ERROR_HANDLE_AND_REPORT((*bs_ptr)->atm == NULL || (*bs_ptr)->bas == NULL, stdl_basis_delete(*bs_ptr); return STDL_ERR_MALLOC, "malloc");

    (*bs_ptr)->env = malloc(env_size * sizeof(double));
    STDL_ERROR_HANDLE_AND_REPORT((*bs_ptr)->env == NULL, stdl_basis_delete(*bs_ptr); return STDL_ERR_MALLOC, "malloc");

    return STDL_ERR_OK;
}

int stdl_basis_delete(stdl_basis *bs) {
    assert(bs != NULL);

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

int stdl_basis_print(stdl_basis *bs, int denormalize) {
    assert(bs != NULL);

    printf("-- basis set containing %d (%s) basis functions, on %d centers --\n", bs->nbas, bs->use_spherical? "spherical": "cartesian", bs->natm);
    if(denormalize)
        printf("-- note: coefficients that are printed are denormalized --\n");

    for(int i=0; i < bs->nbas; i++) {
        int nprim = bs->bas[i*8+2], ncont = bs->bas[i*8+3];
        printf("%d -- %c-type cGTOs (%d*%d AOs), centered on atom #%d.\n", i + 1, _toc(bs->bas[i*8+1]), ncont, CINTcgtos_cart(i, bs->bas), bs->bas[i*8+0] + 1);
        printf("  Defines %d cGTOs composed of %d GTOs:\n", ncont, nprim);
        printf("      (exp)          (conts)\n");
        for(int j=0; j < bs->bas[i*8+2]; j++){
            printf("  %.8e", bs->env[bs->bas[i*8+5]+j]);
            for(int h=0; h < ncont; h++) {
                printf(" % .8e", bs->env[bs->bas[i*8+6]+j * ncont + h] * (denormalize ? 1. / CINTgto_norm(bs->bas[i*8+1], bs->env[bs->bas[i*8+5]+j]) : 1));
            }
            printf("\n");
        }
    }
    printf("-- end --\n");

    return STDL_ERR_OK;
}

int stdl_basis_compute_dsy_ovlp(stdl_basis* bs, size_t nao, double* S) {
    assert(bs != NULL && S != NULL && nao >= (size_t) bs->nbas);

    STDL_DEBUG("computing <i|j> to create the S matrix");

    int si, sj, ioffset=0, joffset;

    double* buff= malloc(28 * 28 * sizeof(double)); // the maximum libcint can handle
    STDL_ERROR_HANDLE_AND_REPORT(buff == NULL, return STDL_ERR_MALLOC, "malloc");

    for(int i=0; i < bs->nbas; i++) {
        if(bs->use_spherical)
            si = CINTcgto_spheric(i, bs->bas);
        else
            si = CINTcgtos_cart(i, bs->bas);

        joffset = 0;

        for(int j=0; j <= i; j++) {
            if(bs->use_spherical) {
                sj = CINTcgto_spheric(j, bs->bas);
                int1e_ovlp_sph(buff, NULL, (int[]) {i, j}, bs->atm, bs->natm, bs->bas, bs->nbas, bs->env, NULL, NULL);
            }
            else {
                sj = CINTcgtos_cart(j, bs->bas);
                int1e_ovlp_cart(buff, NULL, (int[]) {i, j}, bs->atm, bs->natm, bs->bas, bs->nbas, bs->env, NULL, NULL);
            }

            for(int iprim=0; iprim < si; iprim++) {
                for(int jprim=0; jprim < sj && joffset + jprim <= ioffset + iprim; jprim++) {
                    S[(ioffset + iprim) * nao + joffset + jprim] = S[(joffset + jprim) * nao + ioffset + iprim] = buff[iprim * sj + jprim];
                }
            }

            joffset += sj;
        }

        ioffset += si;
    }

    free(buff);

    return STDL_ERR_OK;
}
