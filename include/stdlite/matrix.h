#ifndef STDLITE_MATRIX_H
#define STDLITE_MATRIX_H

#include <stdlib.h>

/// Maximum number of columns per batch
/// @ingroup matrix
#define STDL_MATRIX_MAX_COLS 5

/**
 * Print a matrix.
 * @param rows number of rows
 * @param columns number of columns
 * @param matrix the matrix
 * @param is_symmetric only prints the lower half.
 * @return `STDL_ERR_OK`
 * @ingroup matrix
 */
int stdl_matrix_numbers_print(size_t rows, size_t columns, double *matrix, int is_symmetric);

#endif //STDLITE_MATRIX_H
