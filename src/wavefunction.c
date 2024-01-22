#include <assert.h>

#include "stdlite/wavefunction.h"
#include "stdlite/errors.h"
#include "stdlite/utils/utils.h"

stdl_wavefunction *stdl_wavefunction_new(size_t natm, size_t nelec, size_t nao, size_t nmo) {
    assert(natm > 0 && nao > 0 && nmo > 0 && nmo <= nao && nmo >= 2 * nelec);

    stdl_wavefunction* wf = malloc(sizeof(stdl_wavefunction));

    if(wf != NULL) {
        wf->natm = natm;
        wf->nao = nao;
        wf->nmo = nmo;
        wf->nelec = nelec;

        wf->atm = wf->S = wf->C = wf->e = NULL;
        wf->aotoatm = NULL;

        wf->atm = malloc(3 * sizeof(double));
        if(wf->atm == NULL) {
            stdl_wavefunction_delete(wf);
            return NULL;
        }

        wf->e = malloc(nao * sizeof(double));
        if(wf->e == NULL){
            stdl_wavefunction_delete(wf);
            return NULL;
        }

        wf->aotoatm = malloc(nao * sizeof(size_t));
        if(wf->aotoatm == NULL){
            stdl_wavefunction_delete(wf);
            return NULL;
        }

        wf->S = malloc(nao * nao * sizeof(double));
        if(wf->S == NULL) {
            stdl_wavefunction_delete(wf);
            return NULL;
        }

        wf->C = malloc(nao * nmo * sizeof(double));
        if(wf->C == NULL) {
            stdl_wavefunction_delete(wf);
            return NULL;
        }
    }

    return wf;
}

int stdl_wavefunction_delete(stdl_wavefunction *wf) {
    assert(wf != NULL);

    STDL_FREE_IFUSED(wf->atm);
    STDL_FREE_IFUSED(wf->e);
    STDL_FREE_IFUSED(wf->aotoatm);
    STDL_FREE_IFUSED(wf->S);
    STDL_FREE_IFUSED(wf->C);

    free(wf);

    return STDL_ERR_OK;
}
