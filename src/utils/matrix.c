#include <stdio.h>
#include <string.h>
#include <lapacke.h>
#include <math.h>
#include <cblas.h>
#include <assert.h>

#include "include/stdlite/utils/matrix.h"
#include "stdlite/logging.h"
#include "stdlite/helpers.h"

int stdl_matrix_dge_print(size_t rows, size_t columns, double *matrix, char *title) {
    assert(rows > 0 && matrix != NULL);

    size_t i = 0, j, c;
    int is_symmetric = columns == 0;

    if(is_symmetric)
        columns = rows;

    if(title != NULL)
        printf(" %s\n\n", title);

    while (i < columns) {
        // header
        j = 0;
        printf("    ");
        while (j < STDL_MATRIX_MAX_COLS && (i + j) < columns){
            c = i + j;
            printf("    %6ld    ", c);
            j++;
        }
        printf("\n");

        // content
        for(size_t r = (is_symmetric? i:0); r < rows; r++) {
            printf("%4ld", r);
            j = 0;
            while (j < STDL_MATRIX_MAX_COLS && i + j < columns && ((is_symmetric && i + j <= r) || !is_symmetric)){
                c = i + j;
                printf(" % .6e", matrix[r * columns + c]);
                j++;
            }
            printf("\n");
        }


        i += STDL_MATRIX_MAX_COLS;
    }

    return STDL_ERR_OK;
}

int stdl_matrix_sge_print(size_t rows, size_t columns, float *matrix, char *title) {
    assert(rows > 0 && matrix != NULL);

    size_t i = 0, j, c;
    int is_symmetric = columns == 0;

    if(is_symmetric)
        columns = rows;

    if(title != NULL)
        printf(" %s\n\n", title);

    while (i < columns) {
        // header
        j = 0;
        printf("    ");
        while (j < STDL_MATRIX_MAX_COLS && (i + j) < columns){
            c = i + j;
            printf("    %6ld    ", c);
            j++;
        }
        printf("\n");

        // content
        for(size_t r = (is_symmetric? i:0); r < rows; r++) {
            printf("%4ld", r);
            j = 0;
            while (j < STDL_MATRIX_MAX_COLS && i + j < columns && ((is_symmetric && i + j <= r) || !is_symmetric)){
                c = i + j;
                printf(" % .6e", matrix[r * columns + c]);
                j++;
            }
            printf("\n");
        }


        i += STDL_MATRIX_MAX_COLS;
    }

    return STDL_ERR_OK;
}

int stdl_matrix_dsp_print(size_t n, double *matrix, char *title) {
    assert(n > 0 && matrix != NULL && matrix != NULL);

    if(title != NULL)
        printf(" %s\n\n", title);

    size_t i = 0, j, c;
    while (i < n) {
        // header
        j = 0;
        printf("    ");
        while (j < STDL_MATRIX_MAX_COLS && (i + j) < n){
            c = i + j;
            printf("    %6ld    ", c);
            j++;
        }
        printf("\n");

        // content
        for(size_t r = i; r < n; r++) {
            printf("%4ld", r);
            j = 0;
            while (j < STDL_MATRIX_MAX_COLS && i + j < n && i + j <= r){
                c = i + j;
                printf(" % .6e", matrix[STDL_MATRIX_SP_IDX(r, c)]);
                j++;
            }
            printf("\n");
        }

        i += STDL_MATRIX_MAX_COLS;
    }

    return STDL_ERR_OK;
}


int stdl_matrix_ssp_print(size_t n, float *matrix, char *title) {
    assert(n > 0 && matrix != NULL && matrix != NULL);

    if(title != NULL)
        printf(" %s\n\n", title);

    size_t i = 0, j, c;
    while (i < n) {
        // header
        j = 0;
        printf("    ");
        while (j < STDL_MATRIX_MAX_COLS && (i + j) < n){
            c = i + j;
            printf("    %6ld    ", c);
            j++;
        }
        printf("\n");

        // content
        for(size_t r = i; r < n; r++) {
            printf("%4ld", r);
            j = 0;
            while (j < STDL_MATRIX_MAX_COLS && i + j < n && i + j <= r){
                c = i + j;
                printf(" % .6e", matrix[STDL_MATRIX_SP_IDX(r, c)]);
                j++;
            }
            printf("\n");
        }

        i += STDL_MATRIX_MAX_COLS;
    }

    return STDL_ERR_OK;
}


int stdl_matrix_dsp_sqrt(size_t n, double *mat) {
    assert(mat != NULL && mat != NULL && n > 0);

    double* matsy = malloc(n * n * sizeof(double));
    STDL_ERROR_HANDLE_AND_REPORT(matsy == NULL, return STDL_ERR_MALLOC, "malloc");

    stdl_matrix_dsp_sqrt_sy(n, mat, matsy);

    LAPACKE_dtrttp(LAPACK_ROW_MAJOR, 'L', (int) n, matsy, (int) n, mat);

    STDL_FREE_ALL(matsy);

    return STDL_ERR_OK;
}


int stdl_matrix_dsp_sqrt_sy(size_t n, double *mat, double* matsy) {
    assert(mat != NULL && mat != NULL && n > 0 && matsy != NULL);

    size_t sz = n * n * sizeof(double);

    double* e = malloc(n * sizeof(double));
    double* w = malloc(sz);
    double* wcc = malloc(sz);

    STDL_ERROR_HANDLE_AND_REPORT(e == NULL || w == NULL || wcc == NULL, STDL_FREE_ALL(e, w, wcc); return STDL_ERR_MALLOC, "malloc");

    // eig
    int info = LAPACKE_dspev(LAPACK_ROW_MAJOR, 'V', 'L', (int) n, mat, e, w, (int) n);

    STDL_ERROR_HANDLE_AND_REPORT(info != 0, STDL_FREE_ALL(e, w, wcc); return STDL_ERR_MALLOC, "dsyev() returned %d", info);

    // compute the square root of eigenvalues
    for(size_t i = 0; i < n; i++) {
        if(e[i] < .0) {
            STDL_WARN("eigenvalue of S #%d is < .0, will be set to 0", i);
            e[i] = 0;
        } else
            e[i] = sqrt(e[i]);
    }

    // wcc = e * w
    memcpy(wcc, w, sz);
    for(size_t i = 0; i < n; i++) {
        for(size_t j=0; j < n; j++)
            wcc[i * n + j]  *= e[j];
    }

    // (*mat)^1/2 = w * wcc^T
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
                (int) n, (int) n, (int) n,
                1.f, w, (int) n,
                wcc, (int) n,
                .0f, matsy, (int) n
    );

    STDL_FREE_ALL(e, w, wcc);

    return STDL_ERR_OK;
}

int stdl_matrix_sge_transpose(size_t nrows, size_t ncols, float* mat) {
    assert(nrows > 0 && ncols > 0 && mat != NULL);

    size_t start, next, i;
    float tmp;

    for (start = 0; start <= nrows * ncols - 1; start++) {
        next = start;
        i = 0;
        do {
            i++;
            next = (next % nrows) * ncols + next / nrows;
        } while (next > start);
        if (next < start || i == 1) continue;

        tmp = mat[next = start];
        do {
            i = (next % nrows) * ncols + next / nrows;
            mat[next] = (i == start) ? tmp : mat[i];
            next = i;
        } while (next > start);
    }

    return STDL_ERR_OK;
}

int stdl_matrix_dsp_blowsy(size_t n, double* in, double* out) {
    assert(n > 0 && in != NULL && out != NULL);

    LAPACKE_dtpttr(LAPACK_ROW_MAJOR, 'L', (int) n, in, out, (int) n);

    for(size_t i = 0; i < n; i++) {
        for (size_t j = i; j < n; ++j) {
            out[i * n + j] = out[j * n + i];
        }
    }

    return STDL_ERR_OK;
}

int stdl_matrix_ssp_blowsy(size_t n, float * in, float * out) {
    assert(n > 0 && in != NULL && out != NULL);

    LAPACKE_stpttr(LAPACK_ROW_MAJOR, 'L', (int) n, in, out, (int) n);

    for(size_t i = 0; i < n; i++) {
        for (size_t j = i; j < n; ++j) {
            out[i * n + j] = out[j * n + i];
        }
    }

    return STDL_ERR_OK;
}
