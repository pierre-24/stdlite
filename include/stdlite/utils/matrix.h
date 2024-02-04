#ifndef STDLITE_MATRIX_H
#define STDLITE_MATRIX_H

#include <stdlib.h>

/**
 * Maximum number of columns per batch
 * @ingroup matrix
 */
#define STDL_MATRIX_MAX_COLS 5

/**
 * Get the index corresponding to an element of a symmetric packed (SP) matrix.
 * `i` is the row, while `j` is the column.
 * @ingroup matrix
 */
static inline size_t STDL_MATRIX_SP_IDX(size_t i, size_t j) {
    return (i >= j) ? i*(i+1) / 2 + j : j*(j+1) / 2 + i;
}

/**
 * Get the size of a symmetric packed (SP) matrix of size `n`.
 * @ingroup matrix
 */
static inline size_t STDL_MATRIX_SP_SIZE(size_t nx) {
    return nx*(nx+1)/2;
}

/**
 * Print a double precision matrix.
 * @param rows number of rows, must be >0.
 * @param columns number of columns. If 0, assume that the matrix is symmetric.
 * @param matrix the matrix
 * @return `STDL_ERR_OK`
 * @ingroup matrix
 */
int stdl_matrix_dge_print(size_t rows, size_t columns, double *matrix, char *title);

/**
 * Print a single precision matrix.
 * @param rows number of rows, must be >0.
 * @param columns number of columns. If 0, assume that the matrix is symmetric.
 * @param matrix the matrix
 * @return `STDL_ERR_OK`
 * @ingroup matrix
 */
int stdl_matrix_sge_print(size_t rows, size_t columns, float *matrix, char *title);

/**
 * Print a symmetric (thus square) packed matrix.
 * Storage is in the [lower triangle form](https://netlib.org/lapack/lug/node123.html) (`L`), so: `[a_00, a_10, a_11, a_20, a_21, ..., a_NN]`.
 * @param n side length of the matrix
 * @param matrix the matrix
 * @return `STDL_ERR_OK`
 * @ingroup matrix
 */
int stdl_matrix_dsp_print(size_t n, double *matrix);

/**
 * Print a symmetric (thus square) packed matrix.
 * Storage is in the [lower triangle form](https://netlib.org/lapack/lug/node123.html) (`L`), so: `[a_00, a_10, a_11, a_20, a_21, ..., a_NN]`.
 * @param n side length of the matrix
 * @param matrix the matrix
 * @return `STDL_ERR_OK`
 * @ingroup matrix
 */
int stdl_matrix_ssp_print(size_t n, float *matrix);


/**
 * Compute the square root of a square unitary matrix in double precision, `*mat`, in place.
 *
 * $$X^{1/2} = U\,\varepsilon^{1/2}\,U^T,$$
 *
 * where $\varepsilon$ are the eigenvalues, and $U$ are the eigenvectors.
 *
 * @param[in,out] mat Unitary matrix to be modified.
 * @param n side length of the matrix
 * @return error code
 * @ingroup matrix
 */
int stdl_matrix_dge_sqrt(size_t n, double **mat);

#endif //STDLITE_MATRIX_H
