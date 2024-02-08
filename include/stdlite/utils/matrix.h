#ifndef STDLITE_MATRIX_H
#define STDLITE_MATRIX_H

#include <stdlib.h>

/**
 * Maximum number of columns per batch
 * @ingroup matrix
 */
#ifndef STDL_MATRIX_MAX_COLS
#define STDL_MATRIX_MAX_COLS 5
#endif

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
 * @param matrix `double[rows*columns]` the matrix
 * @param title title to be printed if not `NULL`.
 * @return `STDL_ERR_OK`
 * @ingroup matrix
 */
int stdl_matrix_dge_print(size_t rows, size_t columns, double *matrix, char *title);

/**
 * Print a single precision matrix.
 * @param rows number of rows, must be >0.
 * @param columns number of columns. If 0, assume that the matrix is symmetric.
 * @param matrix `float[rows*columns]`  the matrix
 * @param title title to be printed if not `NULL`.
 * @return `STDL_ERR_OK`
 * @ingroup matrix
 */
int stdl_matrix_sge_print(size_t rows, size_t columns, float *matrix, char *title);

/**
 * Print a symmetric (thus square) packed matrix.
 *
 * @param n side length of the matrix
 * @param `double[STDL_MATRIX_SP_SIZE(n)]`  matrix the matrix
 * @param title title to be printed if not `NULL`.
 * @return `STDL_ERR_OK`
 * @ingroup matrix
 */
int stdl_matrix_dsp_print(size_t n, double *matrix, char *title);

/**
 * Print a symmetric (thus square) packed matrix.
 *
 * @param n side length of the matrix
 * @param `float[STDL_MATRIX_SP_SIZE(n)]` matrix the matrix
 * @param title title to be printed if not `NULL`.
 * @return `STDL_ERR_OK`
 * @ingroup matrix
 */
int stdl_matrix_ssp_print(size_t n, float *matrix, char *title);


/**
 * Compute the square root of a `sp` matrix (in double precision), `*mat`, in place.
 *
 * $$X^{1/2} = U\,\varepsilon^{1/2}\,U^T,$$
 *
 * where $\varepsilon$ are its eigenvalues, and $U$ are its eigenvectors.
 *
 * @param[in,out] mat `double[STDL_MATRIX_SP_SIZE(n)]` Unitary matrix to be modified.
 * @param n side length of the matrix
 * @return error code
 * @ingroup matrix
 */
int stdl_matrix_dsp_sqrt(size_t n, double *mat);


/**
 * Compute the square root of a `sp` matrix (in double precision), and output in `sy` format.
 *
 * $$X^{1/2} = U\,\varepsilon^{1/2}\,U^T,$$
 *
 * where $\varepsilon$ are the eigenvalues, and $U$ are the eigenvectors.
 *
 * @param mat `double[STDL_MATRIX_SP_SIZE(n)]` a matrix
 * @param n side length of the matrix
 * @param matsy `double[n,n]` the square root of `mat`
 * @return error code
 * @ingroup matrix
 */
int stdl_matrix_dsp_sqrt_sy(size_t n, double *mat, double* matsy);


/**
 * Compute the square root of a `sp` matrix (in single precision), `*mat`, in place.
 *
 * $$X^{1/2} = U\,\varepsilon^{1/2}\,U^T,$$
 *
 * where $\varepsilon$ are its eigenvalues, and $U$ are its eigenvectors.
 *
 * @param[in,out] mat `float[STDL_MATRIX_SP_SIZE(n)]` Unitary matrix to be modified.
 * @param n side length of the matrix
 * @return error code
 * @ingroup matrix
 */
int stdl_matrix_ssp_sqrt(size_t n, float *mat);

/**
 * Compute the square root of a `sp` matrix (in single precision), and output in `sy` format.
 *
 * $$X^{1/2} = U\,\varepsilon^{1/2}\,U^T,$$
 *
 * where $\varepsilon$ are the eigenvalues, and $U$ are the eigenvectors.
 *
 * @param mat `floay[STDL_MATRIX_SP_SIZE(n)]` a matrix
 * @param n side length of the matrix
 * @param matsy `float[n,n]` the square root of `mat`
 * @return error code
 * @ingroup matrix
 */
int stdl_matrix_ssp_sqrt_sy(size_t n, float *mat, float * matsy);


/**
 * Transpose a (single precision) rectangular matrix in place.
 * From <https://rosettacode.org/wiki/Matrix_transposition#C> and <https://en.wikipedia.org/wiki/In-place_matrix_transposition>.
 *
 * @param nrows number of rows
 * @param ncols number of columns
 * @param mat `float[nrows * ncols]`, matrix to be transposed
 * @return error code
 * @ingroup matrix
 */
int stdl_matrix_sge_transpose(size_t nrows, size_t ncols, float* mat);

/**
 * Blow a (double-precision) symmetry packed matrix (`in`) into a full storage (symmetric, `sy`) matrix (`out`).
 *
 * @param n order of `in` and `out`
 * @param uplo whether `in` and `out` are upper (`U`) or lower (`L`) triangular.
 * @param in `double[STDL_SP_SIZE(n)]` a matrix in symmetry packed storage
 * @param[in,out] out `double[n,n]` the resulting matrix.
 * @return error code
 * @ingroup matrix
 */
int stdl_matrix_dsp_blowsy(size_t n, char uplo, double *in, double *out);

/**
 * Blow a (single-precision) symmetry packed matrix (`in`) into a full storage (symmetric, `sy`) matrix (`out`).
 *
 * @param n order of `in` and `out`
 * @param uplo whether `in` and `out` are upper (`U`) or lower (`L`) triangular.
 * @param in `float[STDL_SP_SIZE(n)]` a matrix in symmetry packed storage
 * @param[in,out] out `float[n,n]` the resulting matrix.
 * @return error code
 * @ingroup matrix
 */
int stdl_matrix_ssp_blowsy(size_t n, char uplo, float *in, float *out);

/**
 * Shrink a (double-precision) full storage (symmetric, `sy`) matrix (`in`) into a symmetry packed matrix (`out`).
 *
 * @param n order of `in` and `out`
 * @param uplo whether `in` and `out` are upper (`U`) or lower (`L`) triangular.
 * @param in `double[n,n]` a matrix in full symmetric (`sp`) storage
 * @param[in,out] out `double[STDL_SP_SIZE(n)]` the resulting matrix.
 * @return error code
 * @ingroup matrix
 */
int stdl_matrix_dsy_shrinksp(size_t n, char uplo, double *in, double *out);

/**
 * Shrink a (single-precision) full storage (symmetric, `sy`) matrix (`in`) into a symmetry packed matrix (`out`).
 *
 * @param n order of `in` and `out`
 * @param uplo whether `in` and `out` are upper (`U`) or lower (`L`) triangular.
 * @param in `float[n,n]` a matrix in full symmetric (`sp`) storage
 * @param[in,out] out `float[STDL_SP_SIZE(n)]` the resulting matrix.
 * @return error code
 * @ingroup matrix
 */
int stdl_matrix_ssy_shrinksp(size_t n, char uplo, float *in, float *out);

#endif //STDLITE_MATRIX_H
