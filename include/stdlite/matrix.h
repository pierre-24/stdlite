#ifndef STDLITE_MATRIX_H
#define STDLITE_MATRIX_H

#include <stdlib.h>

/// Maximum number of columns per batch
/// @ingroup matrix
#define STDL_MATRIX_MAX_COLS 8

/**
 * Print a matrix.
 * @param rows number of rows
 * @param columns number of columns
 * @param matrix the matrix
 * @param is_symmetric only prints the lower half.
 * @return `STDL_ERR_OK`
 * @ingroup matrix
 */
int stdl_matrix_ge_print(size_t rows, size_t columns, double *matrix, int is_symmetric);

/**
 * Print a symmetric (thus square) packed matrix.
 * Storage is in the [lower triangle form](https://netlib.org/lapack/lug/node123.html) (`L`), so: `[a_00, a_10, a_11, a_20, a_21, ..., a_NN]`.
 * @param n size of the matrix
 * @param matrix the matrix
 * @return `STDL_ERR_OK`
 * @ingroup matrix
 */
int stdl_matrix_sp_print(size_t n, double *matrix);


/**
 * Get the size of a symmetric packed (SP) matrix.
 * @ingroup matrix
 */
# define STDL_MATRIX_SP_SIZE(n) (n)*((n)+1)/2

/**
 * Get the index corresponding to an element of a symmetric packed (SP) matrix.
 * `i` is the row, while `j` is the column.
 * Assume `i < j`.
 * @ingroup matrix
 */
#define STDL_MATRIX_SP_IDX(i, j) (i)*(i+1)/2+j

#endif //STDLITE_MATRIX_H
