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

    /// pointer to the energies for the selected MOs
    double* e_ptr;

    /// pointer to the (non orthogonal) MO coefficient for the selected MO
    double* C_ptr;

    /// `double[nmo*original_wf->nao]` **orthogonal** MO coefficients for the selected MO.
    double* C;

    /// number of CSFs selected by a given scheme (monopole, ...). Zero as long as no CSFs has been selected.
    size_t ncsfs;

    /// `size_t[ncsfs]`, the indices (`kia = i * nvirt + a - nocc`) of each selected CSF `i→a` (as `i = csfs[kia] / nvirt; a = csfs[kia] % nvirt + nocc`).
    /// They are given in increasing energy order, the energy being available at `ecsfs[kia]`.
    /// `NULL` as long as no CSFs has been selected.
    size_t* csfs;

    /// `float[STDL_MATRIX_SP_SIZE(ncfs)]`, part of the electronic Hessian matrix. `NULL` as long as no CSFs has been selected.
    float* A;

    /// `float[STDL_MATRIX_SP_SIZE(ncfs)]`, part of the electronic Hessian matrix. Might be `NULL` if only `A` is required.
    float* B;
};

typedef struct stdlite_context_ stdl_context;

/**
 * Create a new context.
 * Based on the information provided as input, a subset of MO from the wavefunction is selected, following the selection rules of sTD-DFT.
 * @param wf A valid wavefunction.
 * @param gammaJ parameter for Coulomb ops_integrals approximation
 * @param gammaK parameter for exchange ops_integrals approximation
 * @param ethr energy threshold for primary CSFs, must be >0
 * @param e2thr energy threshold for secondary CSFs (selected pertubatively), must be >0
 * @param ax amount of HF exchange
 * @param[out] ctx_ptr The context to be created.
 * @return error code
 * @ingroup context
 */
int stdl_context_new(stdl_wavefunction *wf, stdl_basis *bs, float gammaJ, float gammaK, float ethr, float e2thr, float ax, stdl_context **ctx_ptr);

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
 * @param[out] csfs `size_t[nselected]`, the indices (`i*ctx->nvirt + a`) of each selected CSF `i→a`, as `i = csfs[k] / ctx->nvirt; a = csfs[k] % ctx->nvirt`.
 *             They are sorted in increasing energy order, the energy being available at `A[k * nselected + k]`.
 * @param[out] A `float[STDL_MATRIX_SP_SIZE(nslected)]`, part of the electronic Hessian matrix, in `sp` format, to be created. Caller is responsible for free'ing it.
 * @param[out] B `float[STDL_MATRIX_SP_SIZE(nslected)]`, part of the electronic Hessian matrix, in `sp` format, to be created. Might be `NULL` if only `A` is required. If not, caller is responsible for free'ing it.
 * @return error code
 * @ingroup context
 */
int stdl_context_select_csfs_monopole(stdl_context *ctx, int compute_B);


/**
 * Dump a context in a H5 file
 *
 * @param ctx a valid context
 * @param file_id a valid H5 file id, opened in write mode
 * @return error code
 * @ingroup context
 */
int stdl_context_dump_h5(stdl_context* ctx, hid_t file_id);

/**
 * Load a context from a H5 file
 *
 * @param file_id a valid H5 file id
 * @param[out] ctx_ptr resulting context
 * @return error code
 * @ingroup context
 */
int stdl_context_load_h5(hid_t file_id, stdl_context** ctx_ptr);

/**
 * Get the approximate space in memory
 * @param ctx a valid context
 * @param[out] sz the total size (including basis set and wavefunction)
 * @param[out] bs_sz the size of the contained basis set
 * @param[out] wf_sz the size of the contained wavefunction
 * @return error code
 * @ingroup context
 */
int stdl_context_approximate_size(stdl_context* ctx, size_t* sz, size_t* bs_sz, size_t* wf_sz);


#endif //STDLITE_CONTEXT_H
