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

    /// Parameter for the method: $E_{thr}$, the threshold for selecting MOs and primary CSFs.
    float ethr;

    /// Parameter for the method: $E^{(2)}_{thr}$, the threshold for (perturbatively) selecting secondary CSFs.
    float e2thr;

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
 * @param ethr energy threshold for primary CSFs, must be >0
 * @param e2thr energy threshold for secondary CSFs (selected pertubatively), must be >0
 * @param ax amount of HF exchange
 * @return error code
 * @ingroup context
 */
int stdl_context_new(stdl_context **ctx, stdl_wavefunction *wf, stdl_basis *bs, float gammaJ, float gammaK, float ethr,
                     float e2thr, float ax);


/**
 * Delete a context. Also delete the `stdl_wavefunction` and `stdl_basis` that it contains.
 * @param ctx the context to be deleted
 * @return error code
 * @ingroup context
 */
int stdl_context_delete(stdl_context* ctx);

/**
 * Select CSFs and build the $\mathbf A$ and $\mathbf B$ matrices, using the monopole approximation (original sTD-DFT).
 * If `B` is set to `NULL`, then only $\mathbf A$ is filled (Tamm-Dancoff approximation).
 *
 * @param ctx a valid context
 * @param[out] nselected number of CSFs that were selected. If equals to 0, then `csfs` and `A` are not initialized.
 * @param[out] `size_t[nselected]`, the indices (`i*ctx->nvirt + a`) of each selected CSF `i→a`, as `i = csfs[k] / ctx->nvirt; a = csfs[k] % ctx->nvirt`. They are sorted in increasing energy order, the energy being available at `A[k * nselected + k]`.
 * @return error code
 * @ingroup context
 */
int stdl_context_select_csfs_monopole(stdl_context *ctx, size_t *nselected, size_t **csfs, float **A, float **B);


#endif //STDLITE_CONTEXT_H
