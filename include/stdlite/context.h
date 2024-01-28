#ifndef STDLITE_CONTEXT_H
#define STDLITE_CONTEXT_H

#include <stdlite/wavefunction.h>

/**
 * Object for the configuration and storage of intermediates of a sTD-DFT calculation.
 * All internal variables uses `float` instead of `double` for speed and memory (to the price of accuracy).
 * @ingroup context
 */
struct stdlite_context_ {
    /// Input wavefunction
    stdl_wavefunction* wf;

    /// Parameter for the method: $\gamma_J$.
    float gamma_J;

    /// Parameter for the method: $\gamma_K$.
    float gamma_K;

    /// Parameter for the method: $\varepsilon_{thr}$, the threshold for selecting MO.
    float e_thr;

    /// Parameter for the method: $a_x$, the amount of HF exchange.
    float ax;

    /// Number of MO considered in the calculation, so that `nmo <= wf->nmo && nocc + nvirt == nmo`.
    size_t nmo;

    /// Number of occupied MO considered in the calculation, with `nocc < nmo && nocc + nvirt == nmo`
    size_t nocc;

    /// Number of virtual (non-occupied) MO, with `nvirt < nmo && nocc + nvirt == nmo`
    size_t nvirt;

    /// `float[nmo]` Energy of MOs
    float* e;

    /// `float[nmo*wf->nao]` orthogonal MO coefficients for the selected MO.
    float* C;
};

typedef struct stdlite_context_ stdlite_context;

/**
 * Create a new context.
 * Based on the information provided as input, a subset of MO from the wavefunction is selected, following the selection rules of sTD-DFT.
 * @param ctx The context to be created.
 * @param wf A valid wavefunction
 * @param gamma_J parameter for Coulomb integrals approximation
 * @param gamma_K parameter for exchange integrals approximation
 * @param e_thr energy threshold
 * @param ax amount of HF exchange
 * @return error code
 */
int stdlite_context_new(stdlite_context ** ctx, stdl_wavefunction* wf, float gamma_J, float  gamma_K, float e_thr, float ax);


/**
 * Delete a context. Also delete the `stdl_wavefunction` that it contains.
 * @param ctx the context to be deleted
 * @return error code
 */
int stdlite_context_delete(stdlite_context* ctx);


#endif //STDLITE_CONTEXT_H
