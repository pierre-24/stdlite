#ifndef STDLITE_BASIS_H
#define STDLITE_BASIS_H

#include <stdlib.h>

#include <cint_funcs.h>

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
     * - `bas[i*8+4]`: not relevant here,
     * - `bas[i*8+5]`: offset at which the `nprim` exponents of GTO are found in `env`,
     * - `bas[i*8+6]`: offset at which the `nprim*ncont` contraction coefficients of GTO are found in `env`,
     * - slot 7 is irrelevant.
     */
    int* bas;

    /**
     * `double[natm*3+??]`, array which stores the coordinates, GTO exponents, and contraction coefficients.
     * In practice, the `natm*3` coordinates are stored first.
     * Then for each basis function, the `nprim` exponents followed by the `nprim*ncont` coefficients are stored.
     * For simplicity, each primitive should be normalized.
     */
    double* env;
};

typedef struct stdl_basis_ stdl_basis;

/**
 * Create a new basis set.
 * Initialize all the arrays.
 *
 * @warning According to [this](https://github.com/sunqm/libcint/issues/76), the first 20 values of `env` are reserved.
 *          One should thus account for that when computing `env_size` (20 should be added to what is required) and when filling it (first 20th values should be zero, and the first non-zero value should be at `env[20]`).
 *
 * @param natm Number of atoms
 * @param nbas Number of basis functions
 * @param env_size size of the `env` array, with `env_size > 3*natm + 20`.
 * @param[out] bs_ptr Basis set object to be created.
 * @return `STDL_ERR_OK` if everything went well.
 * @ingroup basis
 */
int stdl_basis_new(int natm, int nbas, size_t env_size, int use_spherical, stdl_basis **bs_ptr);

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
 * @param denormalize If set tpo 0, prints the actual coefficients stored in `env`. If not, prints without normalization.
 * @return `STDL_ERR_OK`
 * @ingroup basis
 */
int stdl_basis_print(stdl_basis *bs, int denormalize);

/**
 * Reorder the coefficients so that they match the one of libcint, using `transpose`.
 * @param nmo Number of MO, must be `nmo <= nao`
 * @param nao Number of AO must be > 0.
 * @param C `double[nmo*nao]`, the LCAO coefficients
 * @param bs a valid basis
 * @param maxshell the maximum shell that the transposition array contains
 * @param transpose `int*[maxshells]`, an array of transpositions
 * @ingroup basis
 */
int stdl_basis_reorder_C(size_t nmo, size_t nao, double *C, stdl_basis *bs, size_t maxshell, int **transpose);

/**
 * Compute the overlap matrix, $S_{\mu\nu} = \braket{\mu|\nu}$ (double precision).
 * @param bs a valid basis set
 * @param[out] S `double[STDL_MATRIX_SP_SIZE(nao)]` the resulting overlap matrix.
 * @return error code.
 * @ingroup basis
 */
int stdl_basis_dsp_ovlp(stdl_basis *bs, double *S);


/**
 * Compute the electronic dipole matrix in AO basis, $D_{\mu\nu} = \braket{\mu|e\,(\hat r - R_0)|\nu}$.
 * @param bs a valid basis set
 * @param[out] dipoles `float[3, STDL_MATRIX_SP_SIZE(nao)]` the resulting dipole matrix. The component of the dipole is thus the slowest varying index.
 * @return error code.
 * @ingroup basis
 */
int stdl_basis_dsp_dipole(stdl_basis *bs, double *dipoles);

#endif //STDLITE_BASIS_H
