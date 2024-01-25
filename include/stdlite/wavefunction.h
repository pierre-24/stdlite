#ifndef STDLITE_WAVEFUNCTION_H
#define STDLITE_WAVEFUNCTION_H

#include <stdlib.h>

/**
 * Structure that represent a wavefunction.
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

    /// Number of electrons in the system.
    /// Note that those electron should populate MO that are considered, so the number might be different from the one of the actual system (if, say, core AO are not considered, then `nelect` should be equal the number of valence electron, not the total).
    size_t nelec;

    /// Number of AO.
    size_t nao;

    /// Number of MO. Should fulfill `nelec <= 2*nmo <= 2*nao`.
    size_t nmo;

    /// `size_t[nao]`, 0-based list of corresponding atom for each AO.
    size_t* aotoatm;

    /// `double[nao*nao]`, the (symmetric) overlap matrix (`S`)
    double* S;

    /// `double[nao*nmo]`, the LCAO coefficients matrix (`C`)
    double* C;

    /// `double[nmo]`, the MO energy vector (`e`) for each MO.
    double* e;
};

typedef struct stdl_wavefunction_ stdl_wavefunction;

/**
 * Create a new wavefunction.
 * @param wf wavefunction to be initialized
 * @param natm number of atom (should be >0)
 * @param nelec number of electrons (should be >0)
 * @param nao number of atomic orbitals (should be >0)
 * @param nmo number of molecular orbitals (should fulfill `nelec <= 2*nmo <= 2*nao`)
 * @return `SDTL_ERR_OK` if everything was ok.
 * @ingroup wavefunction
 */
int stdl_wavefunction_new(stdl_wavefunction **wf, size_t natm, size_t nelec, size_t nao, size_t nmo);

/**
 * Free the wavefunction.
 * @param wf The wavefunction object
 * @return `SDTL_ERR_OK`
 * @ingroup wavefunction
 */
int stdl_wavefunction_delete(stdl_wavefunction* wf);

#endif //STDLITE_WAVEFUNCTION_H
