#ifndef STDLITE_WAVEFUNCTION_H
#define STDLITE_WAVEFUNCTION_H

#include <stdlib.h>

/**
 * Structure that represent a (closed-shell!) wavefunction.
 * So, in practice, it stores the results of a QC calculation for subsequent (s)TD-DFT calculation.
 * Thus, it contains geometrical information, as well as the `S` (overlap), `C` (LCAO coefficients), and `e` (MO energies) matrices.
 * It is assumed that the MO are ordered in increasing energy values, and thus that electrons sits on the lower levels.
 *
 * Note: this structure takes about `(3 * natm + nmo + nao * (nao + nmo)) * sizeof(double) + nao * sizeof(size_t)` bytes of space when correctly allocated.
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
 * Compute the density matrix (`sy` format).
 *
 * $$P_{\mu\nu} = \sum_r^{MO} n_r\,C_{r\mu}\,C_{r\nu},$$
 *
 * where $n_r$ is the occupation number of MO $r$.
 *
 * @param C `double[nmo*nao]`, the coefficients
 * @param nocc number of occupied orbitals, must be >0.
 * @param nmo number of MO, must be >0.
 * @param nao number of AO, must be >0.
 * @param[out] D `double[nao*nao]` the density matrix.
 * @return error code.
 * @ingroup wavefunction
 */
int stdl_wavefunction_compute_density_dsy(size_t nocc, size_t nmo, size_t nao, double *C, double *D);


/**
 * Convert $X$ expressed from AO to MO basis.
 *
 * $$X^{MO}_{pq} = \sum_{\mu\nu} C_{p\mu}\,X^{AO}_{\mu\nu}\,C_{q\nu}.$$
 *
 * @param nao number of AO, must be >0.
 * @param nmo number of MO, must be `0 < nmo <= nao`.
 * @param C the LCAO coefficients
 * @param X_AO `double[nao,nao]`, the matrix in AO basis
 * @param[out] X_MO `double[nmo,nmo]`, the matrix in MO basis
 * @return error code
 * @ingroup wavefunction
 */
int stdl_wavefunction_dsy_ao_to_mo(size_t nao, size_t nmo, double* C, double* X_AO, double* X_MO);


#endif //STDLITE_WAVEFUNCTION_H
