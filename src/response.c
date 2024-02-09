#include <assert.h>
#include <lapacke.h>
#include <string.h>
#include <stdio.h>
#include <cblas.h>

#include "stdlite/response.h"
#include "stdlite/logging.h"
#include "stdlite/helpers.h"
#include "stdlite/utils/matrix.h"

int stdl_response_TDA_casida(stdl_context *ctx, size_t nexci, float *e, float *X) {
    assert(ctx != NULL && ctx->ncsfs > 0 && nexci > 0 && nexci <= ctx->ncsfs && e != NULL && X != NULL);

    int err;

    if (nexci < ctx->ncsfs) {
        int found = 0;

        int* ifail = malloc(ctx->ncsfs * sizeof(float ));
        STDL_ERROR_HANDLE_AND_REPORT(ifail == NULL, return STDL_ERR_MALLOC, "malloc");

        err = LAPACKE_sspevx(
                LAPACK_ROW_MAJOR, 'V', 'I', 'L',
                (int) ctx->ncsfs, ctx->A,
                .0f, .0f,
                1 /* even though we are in C, it starts at 1 */, (int) nexci, STDL_RESPONSE_EIGV_ABSTOL,
                &found, e, X, (int) nexci, ifail
        );

        STDL_FREE_ALL(ifail);
    } else {
        err = LAPACKE_sspev(
                LAPACK_ROW_MAJOR, 'V', 'L',
                (int) ctx->ncsfs, ctx->A,
                e, X, (int) ctx->ncsfs
        );
    }

    STDL_ERROR_HANDLE_AND_REPORT(err != 0, return STDL_ERR_RESPONSE, "error while sspevx(): %d", err);

    err = stdl_matrix_sge_transpose(ctx->ncsfs, nexci, X);
    STDL_ERROR_CODE_HANDLE(err, return err);

    return STDL_ERR_OK;
}


// turn ctx->A and ctx-B to A+B and A-B, respectively
void _make_apb_amb(stdl_context* ctx) {
    for (size_t kia = 0; kia < ctx->ncsfs; ++kia) {
        for (size_t kjb = 0; kjb <= kia; ++kjb) {
            float a = ctx->A[STDL_MATRIX_SP_IDX(kia, kjb)], b = ctx->B[STDL_MATRIX_SP_IDX(kia, kjb)];
            ctx->A[STDL_MATRIX_SP_IDX(kia, kjb)] = a + b;
            ctx->B[STDL_MATRIX_SP_IDX(kia, kjb)] = a - b;
        }
    }
}

int stdl_response_RPA_casida(stdl_context *ctx, size_t nexci, float *e, float *X, float* Y) {
    assert(ctx != NULL && ctx->ncsfs > 0 && nexci <= ctx->ncsfs && ctx->B != NULL && e != NULL && X != NULL && Y != NULL);

    size_t sz = ctx->ncsfs * ctx->ncsfs;
    float * wrk = malloc(3 * sz * sizeof(float));
    STDL_ERROR_HANDLE_AND_REPORT(wrk == NULL, return STDL_ERR_MALLOC, "malloc");

    float *U = wrk, *V = wrk + sz, *W = wrk + 2 * sz;

    // compute A+B and A-B
    _make_apb_amb(ctx);

    // U=A+B
    stdl_matrix_ssp_blowsy(ctx->ncsfs, 'L', ctx->A, U);

    // get V=(A-B)^1/2
    stdl_matrix_ssp_sqrt_sy(ctx->ncsfs, ctx->B, V);

    // W = (A+B)*(A-B)^1/2 = U*V
    cblas_ssymm(CblasRowMajor, CblasRight, CblasLower,
                (int) ctx->ncsfs, (int)  ctx->ncsfs,
                1.f, U, (int)  ctx->ncsfs,
                V, (int) ctx->ncsfs,
                .0f, W, (int) ctx->ncsfs
    );

    // U = (A-B)^1/2*(A+B)*(A-B)^1/2 = V*W
    cblas_ssymm(CblasRowMajor, CblasRight, CblasLower,
                (int) ctx->ncsfs, (int)  ctx->ncsfs,
                1.f, V, (int)  ctx->ncsfs,
                W, (int) ctx->ncsfs,
                .0f, U, (int) ctx->ncsfs
    );

    int err;

    if (nexci < ctx->ncsfs) {
        int found = 0;

        int* ifail = malloc(ctx->ncsfs * sizeof(float ));
        STDL_ERROR_HANDLE_AND_REPORT(ifail == NULL, return STDL_ERR_MALLOC, "malloc");

        err = LAPACKE_ssyevx(
                LAPACK_ROW_MAJOR, 'V', 'I', 'L',
                (int) ctx->ncsfs, U, (int) ctx->ncsfs,
                .0f, .0f,
                1 /* even though we are in C, it starts at 1 */, (int) nexci, STDL_RESPONSE_EIGV_ABSTOL,
                &found, e, X, (int) nexci, ifail
        );

        STDL_FREE_ALL(ifail);
    } else {
        err = LAPACKE_ssyev(
                LAPACK_ROW_MAJOR, 'V', 'L',
                (int) ctx->ncsfs, U, (int) ctx->ncsfs,
                e
        );

        memcpy(X, U, sz * sizeof(float));
    }

    STDL_ERROR_HANDLE_AND_REPORT(err != 0, return STDL_ERR_RESPONSE, "error while ssyevx(): %d", err);

    err = stdl_matrix_sge_transpose(ctx->ncsfs, nexci, X); // Now, X' = (A-B)^(-1/2)*(X+Y).

    // stdl_matrix_sge_print(nexci, ctx->ncsfs, X, "Z");

    STDL_ERROR_CODE_HANDLE(err, free(wrk); return err);

    // Get energies and compute X and Y
    for (size_t iexci = 0; iexci < nexci; ++iexci) {
        e[iexci]  = sqrtf(e[iexci]);

        // X'_i = (A-B)^(1/2) * X'_i / sqrt(w)
        for(size_t kia = 0; kia < ctx->ncsfs; kia++) {
            float sum = .0f;
            for(size_t kjb = 0; kjb < ctx->ncsfs; kjb++) {
                sum += V[kjb * ctx->ncsfs + kia] * X[iexci * ctx->ncsfs + kjb] / sqrtf(e[iexci]);
            }

            X[iexci * ctx->ncsfs + kia] = sum;
        }

        // Y'_i = 1/w*(A+B)*X'_i
        for(size_t kia = 0; kia < ctx->ncsfs; kia++) {
            float sum = .0f;
            for(size_t kjb = 0; kjb < ctx->ncsfs; kjb++) {
                sum += ctx->A[STDL_MATRIX_SP_IDX(kia, kjb)] * X[iexci * ctx->ncsfs + kjb] / e[iexci];
            }

            Y[iexci * ctx->ncsfs + kia] = sum;
        }

        // X_i = 1/2*(X'_i + Y'_i) && Y_i = 1/2*(X'_i - Y'_i)
        for(size_t kia = 0; kia < ctx->ncsfs; kia++) {
            float u = X[iexci * ctx->ncsfs + kia], v = Y[iexci * ctx->ncsfs + kia];
            X[iexci * ctx->ncsfs + kia] = .5f * (u + v);
            Y[iexci * ctx->ncsfs + kia] = .5f * (u - v);
        }
    }

    STDL_FREE_ALL(wrk);

    return STDL_ERR_OK;
}

int stdl_response_RPA_linear(stdl_context* ctx, float w, size_t dim, float* egrad, float* X, float* Y) {
    // compute A+B and A-B
    _make_apb_amb(ctx);

    // allocate space for the L size of the linear response equation
    float* L = malloc(STDL_MATRIX_SP_SIZE(ctx->ncsfs) * sizeof(float));
    STDL_ERROR_HANDLE_AND_REPORT(L == NULL, return STDL_ERR_MALLOC, "malloc");

    // invert A-B, taking advantage of its `sp` storage
    int* ipiv = malloc(ctx->ncsfs * sizeof(int));
    STDL_ERROR_HANDLE_AND_REPORT(ipiv == NULL, STDL_FREE_ALL(L); return STDL_ERR_MALLOC, "malloc");

    int err = LAPACKE_ssptrf(LAPACK_ROW_MAJOR, 'L', (int) ctx->ncsfs, ctx->B, ipiv);
    STDL_ERROR_HANDLE_AND_REPORT(err != 0,  STDL_FREE_ALL(L, ipiv); return STDL_ERR_RESPONSE, "error while ssptrf(): %d", err);

    err = LAPACKE_ssptri(LAPACK_ROW_MAJOR, 'L', (int) ctx->ncsfs, ctx->B, ipiv);
    STDL_ERROR_HANDLE_AND_REPORT(err != 0, STDL_FREE_ALL(L, ipiv); return STDL_ERR_RESPONSE, "error while ssptri(): %d", err);

    // now, ctx->B contains (A-B)^(-1)
    // make left side: L = (A+B)-w^2*(A-B)^(-1)
    for (size_t kia = 0; kia < ctx->ncsfs; ++kia) {
        for(size_t kjb = 0; kjb <= kia; ++kjb)
            L[STDL_MATRIX_SP_IDX(kia, kjb)] = ctx->A[STDL_MATRIX_SP_IDX(kia, kjb)] - powf(w, 2) * ctx->B[STDL_MATRIX_SP_IDX(kia, kjb)];
    }

    // copy egrad in X, to keep it for latter
    memcpy(X, egrad, ctx->ncsfs * dim * sizeof(float ));

    // solve the problem
    err = LAPACKE_ssptrf(LAPACK_ROW_MAJOR, 'L', (int) ctx->ncsfs, L, ipiv);
    STDL_ERROR_HANDLE_AND_REPORT(err != 0,  STDL_FREE_ALL(L, ipiv); return STDL_ERR_RESPONSE, "error while ssptrf(): %d", err);

    err = LAPACKE_ssptrs(LAPACK_ROW_MAJOR, 'L', (int) ctx->ncsfs, (int) dim, L, ipiv, X, (int) dim);
    STDL_ERROR_HANDLE_AND_REPORT(err != 0, STDL_FREE_ALL(L, ipiv); return STDL_ERR_RESPONSE, "error while ssptrs(): %d", err);

    stdl_matrix_sge_print(ctx->ncsfs, dim, X, "X'");

    // separate X and Y
    // Y' = w*(A-B)^(-1)*X'
    for(size_t kia = 0; kia < ctx->ncsfs; kia++) {
        for(size_t cpt = 0; cpt < dim; cpt++) {
            float sum = .0f;
            for (size_t kjb = 0; kjb < ctx->ncsfs; kjb++)
                sum += ctx->B[STDL_MATRIX_SP_IDX(kia, kjb)] * w * X[kia * dim + cpt];

            Y[kia * dim + cpt] = sum;
        }
    }

    // X = 1/2*(X' + Y') && Y = 1/2*(X' - Y')
    for(size_t kia = 0; kia < ctx->ncsfs; kia++) {
        for(size_t cpt = 0; cpt < dim; cpt++) {
             float u = X[kia * dim + cpt], v = Y[kia * dim + cpt];
             X[kia * dim + cpt] = .5f * (u + v);
             Y[kia * dim + cpt] = .5f * (u - v);
        }
    }

    stdl_matrix_sge_print(ctx->ncsfs, dim, X, "X");
    stdl_matrix_sge_print(ctx->ncsfs, dim, Y, "Y");

    STDL_FREE_ALL(L, ipiv);

    return STDL_ERR_OK;
}
