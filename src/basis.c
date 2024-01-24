#include <assert.h>
#include <stdio.h>
#include "stdlite/basis.h"
#include "stdlite/errors.h"
#include "stdlite.h"


stdl_basis *stdl_basis_new(int natm, int nbas, size_t env_size) {
    assert(natm > 0 && nbas > 0 && env_size > 0);

    stdl_basis* bs = malloc(sizeof(stdl_basis));

    if(bs != NULL) {
        bs->natm = natm;
        bs->nbas = nbas;

        bs->atm = bs->bas = NULL;
        bs->env = NULL;

        bs->atm = malloc(6 * natm * sizeof(int));
        bs->bas = malloc(8 * nbas * sizeof(int));
        if(bs->atm == NULL || bs->bas == NULL) {
            stdl_basis_delete(bs);
            return NULL;
        }

        bs->env = malloc(env_size * sizeof(double));
        if(bs->env == NULL) {
            stdl_basis_delete(bs);
            return NULL;
        }
    }

    return bs;
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
        default:
            return '?';
    }
}

int stdl_basis_print(stdl_basis *bs) {
    assert(bs != NULL);

    printf("-- basis set containing %d basis functions on %d atoms --\n", bs->nbas, bs->natm);
    for(int i=0; i < bs->nbas; i++) {
        int nprim = bs->bas[i*8+2], ncont = bs->bas[i*8+3];
        printf("%d -- %c-like function, centered on atom #%d.\n", i + 1, _toc(bs->bas[i*8+1]), bs->bas[i*8+0] + 1);
        printf("  Defined by %d GTOs and %d cGTOs:\n", nprim, ncont);
        printf("      (exp)          (conts)\n");
        for(int j=0; j < bs->bas[i*8+2]; j++){
            printf("  %.8e", bs->env[bs->bas[i*8+5]+j]);
            for(int h=0; h < ncont; h++) {
                printf(" % .8e", bs->env[bs->bas[i*8+6]+j * ncont + h]);
            }
            printf("\n");
        }
    }
    printf("-- end --\n");

    return STDL_ERR_OK;
}
