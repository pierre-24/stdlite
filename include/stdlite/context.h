#ifndef STDLITE_CONTEXT_H
#define STDLITE_CONTEXT_H

#include <stdlite/wavefunction.h>
#include <stdlite/basis.h>

/**
 * Object for the configuration and storage of intermediates of a sTD-DFT calculation.
 * Contains a subset of MOs from `original_wf`.
 * @ingroup context
 */
struct stdlite_context_ {
    /// Input wavefunction
    stdl_wavefunction* original_wf;

    /// Input basis set
    stdl_basis* bs;

    /// Parameter for the method: $\gamma_J$.
    float gammaJ;

    /// Parameter for the method: $\gamma_K$.
    float gammaK;

    /// Parameter for the method: $\varepsilon_{thr}$, the threshold for selecting MO.
    float ethr;

    /// Parameter for the method: $a_x$, the amount of HF exchange.
    float ax;

    /// Number of MO considered in the calculation, so that `nmo <= original_wf->nmo && nocc + nvirt == nmo`.
    size_t nmo;

    /// Number of occupied MO considered in the calculation, with `nocc < nmo && nocc + nvirt == nmo`
    size_t nocc;

    /// `double[nmo]` Energy of MOs
    double* e;

    /// `double[nmo*original_wf->nao]` orthogonal MO coefficients for the selected MO.
    double* C;
};

typedef struct stdlite_context_ stdl_context;

/**
 * Create a new context.
 * Based on the information provided as input, a subset of MO from the wavefunction is selected, following the selection rules of sTD-DFT.
 * @param ctx The context to be created.
 * @param wf A valid wavefunction.
 * @param gammaJ parameter for Coulomb integrals approximation
 * @param gammaK parameter for exchange integrals approximation
 * @param ethr energy threshold
 * @param ax amount of HF exchange
 * @return error code
 * @ingroup context
 */
int stdlite_context_new(stdl_context ** ctx, stdl_wavefunction* wf, stdl_basis* basis, float gammaJ, float  gammaK, float ethr, float ax);


/**
 * Delete a context. Also delete the `stdl_wavefunction` and `stdl_basis` that it contains.
 * @param ctx the context to be deleted
 * @return error code
 * @ingroup context
 */
int stdlite_context_delete(stdl_context* ctx);


#endif //STDLITE_CONTEXT_H
