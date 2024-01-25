#include <assert.h>

#include "stdlite/wavefunction.h"
#include "stdlite/errors.h"
#include "stdlite.h"

int stdl_wavefunction_new(stdl_wavefunction **wf, size_t natm, size_t nelec, size_t nao, size_t nmo) {
    assert(wf != NULL && natm > 0 && nao > 0 && nmo > 0 && nmo <= nao && 2*nmo >= nelec);

    *wf = malloc(sizeof(stdl_wavefunction));

    if(*wf != NULL) {
        (*wf)->natm = natm;
        (*wf)->nao = nao;
        (*wf)->nmo = nmo;
        (*wf)->nelec = nelec;

        (*wf)->atm = (*wf)->S = (*wf)->C = (*wf)->e = NULL;
        (*wf)->aotoatm = NULL;

        (*wf)->atm = malloc(4 * natm * sizeof(double));
        if((*wf)->atm == NULL) {
            stdl_wavefunction_delete(*wf);
            return STDL_ERR_MALLOC;
        }

        (*wf)->e = malloc(nao * sizeof(double));
        if((*wf)->e == NULL){
            stdl_wavefunction_delete(*wf);
            return STDL_ERR_MALLOC;
        }

        (*wf)->aotoatm = malloc(nao * sizeof(size_t));
        if((*wf)->aotoatm == NULL){
            stdl_wavefunction_delete(*wf);
            return STDL_ERR_MALLOC;
        }

        (*wf)->S = malloc(nao * nao * sizeof(double));
        if((*wf)->S == NULL) {
            stdl_wavefunction_delete(*wf);
            return STDL_ERR_MALLOC;
        }

        (*wf)->C = malloc(nao * nmo * sizeof(double));
        if((*wf)->C == NULL) {
            stdl_wavefunction_delete(*wf);
            return STDL_ERR_MALLOC;
        }

        return STDL_ERR_OK;

    } else
        return STDL_ERR_MALLOC;
}

int stdl_wavefunction_delete(stdl_wavefunction *wf) {
    assert(wf != NULL);

    STDL_FREE_IF_USED(wf->atm);
    STDL_FREE_IF_USED(wf->e);
    STDL_FREE_IF_USED(wf->aotoatm);
    STDL_FREE_IF_USED(wf->S);
    STDL_FREE_IF_USED(wf->C);

    free(wf);

    return STDL_ERR_OK;
}
