#include <assert.h>

#include "stdlite/wavefunction.h"
#include "stdlite/errors.h"
#include "stdlite.h"

int stdl_wavefunction_new(stdl_wavefunction **wf_ptr, size_t natm, size_t nelec, size_t nao, size_t nmo) {
    assert(wf_ptr != NULL && natm > 0 && nao > 0 && nmo > 0 && nmo <= nao && 2 * nmo >= nelec);

    *wf_ptr = malloc(sizeof(stdl_wavefunction));

    if(*wf_ptr != NULL) {
        (*wf_ptr)->natm = natm;
        (*wf_ptr)->nao = nao;
        (*wf_ptr)->nmo = nmo;
        (*wf_ptr)->nelec = nelec;

        (*wf_ptr)->atm = (*wf_ptr)->S = (*wf_ptr)->C = (*wf_ptr)->e = NULL;
        (*wf_ptr)->aotoatm = NULL;

        (*wf_ptr)->atm = malloc(4 * natm * sizeof(double));
        if((*wf_ptr)->atm == NULL) {
            stdl_wavefunction_delete(*wf_ptr);
            return STDL_ERR_MALLOC;
        }

        (*wf_ptr)->e = malloc(nao * sizeof(double));
        if((*wf_ptr)->e == NULL){
            stdl_wavefunction_delete(*wf_ptr);
            return STDL_ERR_MALLOC;
        }

        (*wf_ptr)->aotoatm = malloc(nao * sizeof(size_t));
        if((*wf_ptr)->aotoatm == NULL){
            stdl_wavefunction_delete(*wf_ptr);
            return STDL_ERR_MALLOC;
        }

        (*wf_ptr)->S = malloc(nao * nao * sizeof(double));
        if((*wf_ptr)->S == NULL) {
            stdl_wavefunction_delete(*wf_ptr);
            return STDL_ERR_MALLOC;
        }

        (*wf_ptr)->C = malloc(nao * nmo * sizeof(double));
        if((*wf_ptr)->C == NULL) {
            stdl_wavefunction_delete(*wf_ptr);
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
