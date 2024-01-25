#include <assert.h>
#include <stdio.h>
#include "stdlite/basis.h"
#include "stdlite/errors.h"
#include "stdlite.h"


int stdl_basis_new(stdl_basis **bs_ptr, int natm, int nbas, size_t env_size, int use_spherical) {
    assert(bs_ptr != NULL && natm > 0 && nbas > 0 && env_size > 0);

    *bs_ptr = malloc(sizeof(stdl_basis));

    if(*bs_ptr != NULL) {
        (*bs_ptr)->natm = natm;
        (*bs_ptr)->nbas = nbas;
        (*bs_ptr)->use_spherical = use_spherical;

        (*bs_ptr)->atm = (*bs_ptr)->bas = NULL;
        (*bs_ptr)->env = NULL;

        (*bs_ptr)->atm = malloc(6 * natm * sizeof(int));
        (*bs_ptr)->bas = malloc(8 * nbas * sizeof(int));
        if((*bs_ptr)->atm == NULL || (*bs_ptr)->bas == NULL) {
            stdl_basis_delete(*bs_ptr);
            return STDL_ERR_MALLOC;
        }

        (*bs_ptr)->env = malloc(env_size * sizeof(double));
        if((*bs_ptr)->env == NULL) {
            stdl_basis_delete(*bs_ptr);
            return STDL_ERR_MALLOC;
        }

        return STDL_ERR_OK;
    } else
        return STDL_ERR_MALLOC;
}

int stdl_basis_delete(stdl_basis *bs) {
    assert(bs != NULL);

    STDL_FREE_IF_USED(bs->atm);
    STDL_FREE_IF_USED(bs->bas);
    STDL_FREE_IF_USED(bs->env);

    free(bs);

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
