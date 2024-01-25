#ifndef STDLITE_BASIS_H
#define STDLITE_BASIS_H

#include <stdlib.h>

/**
 * Describe a set of basis functions, centered on atoms.
 * Follows the structure of [`libcint`](https://github.com/sunqm/libcint/blob/master/doc/program_ref.txt), so that it can be used to compute extra integrals.
 *
 * @ingroup basis
 */
struct stdl_basis_ {
    /// Number of atoms
    int natm;

    /**
     * `int[natm*6]` list of atom.
     *
     * For each (ith) atom,
     *
     * - `atm[i*6+0]`: nuclear charge,
     * - `atm[i*6+1]`: offset at which the coordinates are given in `env`,
     * - slots 2-5 are irrelevant
     */
    int* atm;

    /// Number of basis functions
    int nbas;

    /// Indicates that spherical functions should be assumed.
    int use_spherical;

    /**
     * `int[nbas*8]`, list of basis functions.
     *
     * For each (ith) basis function,
     *
     * - `bas[i*8+0]`: corresponding atom (0-based index)
     * - `bas[i*8+1]`: angular momentum
     * - `bas[i*8+2]`: number of primitive GTO, `nprim`,
     * - `bas[i*8+3]`: number of contracted GTO, `ncont`,
     * - `bas[i*8+4]`: not relevant here
     * - `bas[i*8+5]`: offset at which the `nprim` exponents of GTO are found in `env`,
     * - `bas[i*8+6]`: offset at which the `nprim*ncont` contraction coefficients of GTO are found in `env`,
     * - slot 7 is irrelevant.
     */
    int* bas;

    /**
     * `double[natm*3+??]`, array which stores the coordinates, GTO exponents, and contraction coefficients.
     * In practice, the `natm*3` coordinates are stored first.
     * Then for each basis function, the `nprim` exponents followed by the `nprim*ncont` coefficients are stored.
     */
    double* env;
};

typedef struct stdl_basis_ stdl_basis;

/**
 * Create a new basis set.
 * Initialize the arrays
 * @param natm Number of atoms
 * @param nbas Number of basis functions
 * @param env_size size of the `env` array, with `env_size > 3*natm`.
 * @return The initialized basis set object, if any, with all array initialized to their respective size.
 * @ingroup basis
 */
stdl_basis *stdl_basis_new(int natm, int nbas, size_t env_size, int use_spherical);

/**
 * Delete the basis set.
 * @param bs a valid basis set
 * @return `STDL_ERR_OK`
 * @ingroup basis
 */
int stdl_basis_delete(stdl_basis* bs);

/**
 * Print the content of the basis set.
 * @param bs a valid basis set
 * @return `STDL_ERR_OK`
 * @ingroup basis
 */
int stdl_basis_print(stdl_basis* bs);

#endif //STDLITE_BASIS_H
