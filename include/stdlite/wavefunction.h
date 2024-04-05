#ifndef STDLITE_WAVEFUNCTION_H
#define STDLITE_WAVEFUNCTION_H

#include <stdlib.h>

#ifdef USE_HDF5_SERIAL
#include <hdf5/serial/hdf5.h>
#include <hdf5/serial/hdf5_hl.h>
# else
#include <hdf5.h>
#include <hdf5_hl.h>
#endif


/**
 * Structure that represent a (closed-shell!) wavefunction.
 * So, in practice, it stores the results of a QC calculation for subsequent (s)TD-DFT calculation.
 * Thus, it contains geometrical information, as well as the `S` (overlap), `C` (LCAO coefficients), and `e` (MO energies) matrices.
 * It is assumed that the MO are ordered in increasing energy values, and thus that electrons sits on the lower levels.
 *
 * @ingroup wavefunction
 */
struct stdl_wavefunction_ {
    /// Number of atoms
    size_t natm;

    /// `double[natm*4]`, list of atoms with nuclear charge (0) and coordinates (1:3)
    double* atm;

    /// Number of occupied MO. Should be >0.
    size_t nocc;

    /// Number of AO.
    size_t nao;

    /// Number of MO. Should fulfill `nocc <= nmo <= nao`.
    size_t nmo;

    /// `size_t[nao]`, 0-based list of corresponding atom for each AO.
    size_t* aotoatm;

    /// `double[STDL_MATRIX_SP_SIZE(nao)]`, the overlap matrix $S_{\mu\nu}$.
    double* S;

    /// `double[nmo*nao]`, the LCAO coefficients matrix $C_{pq}$.
    /// Stored (from Gaussian) as `[c_mo0ao0, c_mo0ao1, ..., c_mo0aoN, c_mo1ao0, ..., c_moNaoN]`: the first index refers to the MO, the second to the AO.
    double* C;

    /// `double[nmo]`, the MO energy vector $\varepsilon_p$ for each MO.
    double* e;
};

typedef struct stdl_wavefunction_ stdl_wavefunction;

/**
 * Create a new wavefunction.
 * @param natm number of atom (should be >0)
 * @param nocc number of occupied orbitals, must be >0.
 * @param nao number of atomic orbitals (should be >0)
 * @param nmo number of molecular orbitals (should fulfill `nelec <= 2*nmo <= 2*nao`)
 * @param[out] wf_ptr wavefunction object to be initialized
 * @return `SDTL_ERR_OK` if everything was ok.
 * @ingroup wavefunction
 */
int stdl_wavefunction_new(size_t natm, size_t nocc, size_t nao, size_t nmo, stdl_wavefunction **wf_ptr);

/**
 * Free the wavefunction.
 * @param wf The wavefunction object
 * @return `SDTL_ERR_OK`
 * @ingroup wavefunction
 */
int stdl_wavefunction_delete(stdl_wavefunction* wf);

/**
 * Symmetrize the LCAO coefficients in place using a Löwdin orthogonalization.
 *
 * In other words:
 *
 * $$C^{\perp} = S^{1/2}\,C.$$
 *
 * See, *e.g.*, [there](https://booksite.elsevier.com/9780444594365/downloads/16755_10030.pdf).
 *
 * @param nmo number of MO, must be >0.
 * @param nao number of AO, must be >0.
 * @param S `double[STLD_MATRIX_SP_SIZE(nao)]`, the overlap matrix
 * @param[in,out] C `double[nmo*nao]`, the coefficients to be orthogonalized
 * @return the error code
 * @ingroup wavefunction
 */
int stdl_wavefunction_orthogonalize_C_dge(size_t nmo, size_t nao, double *S, double *C);

/**
 * Compute the density matrix.
 *
 * $$P_{\mu\nu} = \sum_r^{MO} n_r\,C_{r\mu}\,C_{r\nu},$$
 *
 * where $n_r$ is the occupation number of MO $r$.
 *
 * @param C `double[nmo*nao]`, the coefficients
 * @param nocc number of occupied orbitals, must be >0.
 * @param nmo number of MO, must be >0.
 * @param nao number of AO, must be >0.
 * @param[out] D `double[STDL_MATRIX_SP_SIZE(nao)]` the density matrix.
 * @return error code.
 * @ingroup wavefunction
 */
int stdl_wavefunction_compute_density_dsp(size_t nocc, size_t nmo, size_t nao, double *C, double *D);


/**
 * Convert $X$ expressed in AO basis to MO basis. Assume a **symmetric** property (i.e., `X_AO[i,j] = X_AO[j,i]`).
 *
 * @param nao number of AO, must be >0.
 * @param nmo number of MO, must be `0 < nmo <= nao`.
 * @param C the LCAO coefficients
 * @param X_AO `double[STDL_MATRIX_SP_SIZE(nao)]`, the matrix in AO basis
 * @param[out] X_MO `double[STDL_MATRIX_SP_SIZE(nmo)]`, the matrix in MO basis
 * @return error code
 * @ingroup wavefunction
 */
int stdl_wavefunction_dsp_ao_to_dsp_mo(size_t nao, size_t nmo, double* C, double* X_AO, double* X_MO);


/**
 * Convert $X$ expressed from AO to MO basis.
 *
 * @param nao number of AO, must be >0.
 * @param nmo number of MO, must be `0 < nmo <= nao`.
 * @param C the LCAO coefficients
 * @param X_AO `double[nao,nao]`, the matrix in AO basis
 * @param[out] X_MO `double[nmo,nmo]`, the matrix in MO basis
 * @return error code
 * @ingroup wavefunction
 */
int stdl_wavefunction_dge_ao_to_dge_mo(size_t nao, size_t nmo, double *C, double *X_AO, double *X_MO);

/**
 * Dump a wavefunction in a H5 file_id
 *
 * @param wf the wavefunction
 * @param file_id a valid H5 file_id identifier
 * @return error code
 * @ingroup wavefunction
 */
int stdl_wavefunction_dump_h5(stdl_wavefunction* wf, hid_t file_id);

/**
 * Load a wavefunction from a H5 file_id
 * @param file_id a valid H5 file_id identifier
 * @param[out] wf_ptr the resulting wavefunction
 * @return error code
 * @ingroup wavefunction
 */
int stdl_wavefunction_load_h5(hid_t file_id, stdl_wavefunction **wf_ptr);

/**
 * Get the approximate space in memory
 *
 * @param wf the wavefunction
 * @param[out] sz an approximate size (in bytes)
 * @return error code
 * @ingroup wavefunction
 */
int stdl_wavefunction_approximate_size(stdl_wavefunction* wf, size_t* sz);


#endif //STDLITE_WAVEFUNCTION_H
