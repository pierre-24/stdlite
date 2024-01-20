#ifndef STDLITE_WAVEFUNCTION_H
#define STDLITE_WAVEFUNCTION_H

#include <stdlib.h>

/**
 * Structure that represent a wavefunction.
 * So, in practice, it stores the results of a QC calculation for subsequent (s)TD-DFT calculation.
 * Thus, it contains geometrical information, as well as the `S` (overlap), `C` (LCAO coefficients), and `e` (MO energies) matrices.
 *
 * @ingroup wavefunction
 */
struct stdl_wavefunction_ {
    /// Number of atoms
    size_t natm;

    /// `double[natm*4]`, list of atoms with nuclear charge (0) and coordinates (1:3)
    double* atm;

    /// Number of AO.
    size_t nao;

    /// Number of MO
    size_t nmo;

    /// `int[nao]`, 0-based list of corresponding atom for each AO.
    int* aotoatm;

    /// `double[nao*nao]`, the overlap matrix (`S`)
    double* S;

    /// `double[nao*nmo]`, the LCAO coefficients matrix (`C`)
    double* C;

    /// `double[nmo]`, the MO energy vector (`e`)
    double* e;
};

typedef struct stdl_wavefunction_ stdl_wavefunction;

/**
 * Create a new wavefunction.
 * @param natm number of atom
 * @param nao number of atomic orbitals
 * @param nmo number of molecular orbitals
 * @return The initialized wavefunction object, if any, with all array initialized to their respective size.
 * @ingroup wavefunction
 */
stdl_wavefunction* stdl_wavefunction_new(size_t natm, size_t nao, size_t nmo);

/**
 * Free the wavefunction.
 * @param wf The wavefunction object
 * @return `SDTL_ERR_OK`
 * @ingroup wavefunction
 */
int stdl_wavefunction_delete(stdl_wavefunction* wf);

#endif //STDLITE_WAVEFUNCTION_H
