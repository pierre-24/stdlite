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

int stdl_response_TD_casida(stdl_context *ctx, size_t nexci, float *e, float *X, float* Y) {
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

int stdl_response_perturbed_gradient(stdl_context* ctx, size_t dim, double* eta_MO, float *egrad) {
    size_t nvirt = ctx->nmo - ctx->nocc;
    for (size_t lia = 0; lia < ctx->ncsfs; ++lia) {
        size_t i = ctx->csfs[lia] / nvirt, a = ctx->csfs[lia] % nvirt + ctx->nocc;
        for (size_t cpt = 0; cpt < dim; ++cpt) {
            egrad[lia * dim + cpt] = -2.f * (float) eta_MO[cpt * STDL_MATRIX_SP_SIZE(ctx->nmo) + STDL_MATRIX_SP_IDX(i, a)];
        }
    }

    return STDL_ERR_OK;
}


int stdl_response_TD_linear(stdl_context *ctx, size_t nw, float *w, size_t ndim, float *egrad, float *X, float *Y) {
    assert(ctx != NULL && ctx->ncsfs > 0 && nw > 0 && ndim > 0 && ctx->B != NULL && egrad != NULL && X != NULL && Y != NULL);

    size_t szXY = ctx->ncsfs * ndim;

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

    for (size_t iw = 0; iw < nw; ++iw) {
        // make left side: L = (A+B)-w^2*(A-B)^(-1)
        for (size_t kia = 0; kia < ctx->ncsfs; ++kia) {
            for(size_t kjb = 0; kjb <= kia; ++kjb)
                L[STDL_MATRIX_SP_IDX(kia, kjb)] = ctx->A[STDL_MATRIX_SP_IDX(kia, kjb)] - powf(w[iw], 2) * ctx->B[STDL_MATRIX_SP_IDX(kia, kjb)];
        }

        float *Xi = X + iw * szXY, *Yi = Y + iw * szXY;

        // copy egrad in X, to keep it for latter
        memcpy(Xi, egrad, ctx->ncsfs * ndim * sizeof(float ));

        // solve the problem
        err = LAPACKE_ssptrf(LAPACK_ROW_MAJOR, 'L', (int) ctx->ncsfs, L, ipiv);
        STDL_ERROR_HANDLE_AND_REPORT(err != 0,  STDL_FREE_ALL(L, ipiv); return STDL_ERR_RESPONSE, "error while ssptrf(): %d", err);

        err = LAPACKE_ssptrs(LAPACK_ROW_MAJOR, 'L', (int) ctx->ncsfs, (int) ndim, L, ipiv, Xi, (int) ndim);
        STDL_ERROR_HANDLE_AND_REPORT(err != 0, STDL_FREE_ALL(L, ipiv); return STDL_ERR_RESPONSE, "error while ssptrs(): %d", err);

        // stdl_matrix_sge_print(ctx->ncsfs, ndim, Xi, "X'");

        // separate X and Y
        // Yi' = w*(A-B)^(-1)*Xi'
        for(size_t kia = 0; kia < ctx->ncsfs; kia++) {
            for(size_t cpt = 0; cpt < ndim; cpt++) {
                float sum = .0f;
                for (size_t kjb = 0; kjb < ctx->ncsfs; kjb++)
                    sum += ctx->B[STDL_MATRIX_SP_IDX(kia, kjb)] * w[iw] * Xi[kia * ndim + cpt];

                Yi[kia * ndim + cpt] = sum;
            }
        }

        // Xi = 1/2*(Xi' + Yi') && Yi = 1/2*(Xi' - Yi')
        for(size_t kia = 0; kia < ctx->ncsfs; kia++) {
            for(size_t cpt = 0; cpt < ndim; cpt++) {
                float u = Xi[kia * ndim + cpt], v = Yi[kia * ndim + cpt];
                Xi[kia * ndim + cpt] = .5f * (u + v);
                Yi[kia * ndim + cpt] = .5f * (u - v);
            }
        }
    }

    // stdl_matrix_sge_print(nw * ctx->ncsfs, ndim, X, "X");
    // stdl_matrix_sge_print(nw * ctx->ncsfs, ndim, Y, "Y");

    STDL_FREE_ALL(L, ipiv);

    return STDL_ERR_OK;
}


int stdl_response_TDA_linear(stdl_context *ctx, size_t nw, float *w, size_t ndim, float *egrad, float *X) {
    assert(ctx != NULL && ctx->ncsfs > 0 && nw > 0 && ndim > 0 && egrad != NULL && X != NULL);

    size_t szXY = ctx->ncsfs * ndim;
    int err;

    // allocate space for the left size of the linear response equation
    float* L = malloc(STDL_MATRIX_SP_SIZE(ctx->ncsfs) * sizeof(float));
    STDL_ERROR_HANDLE_AND_REPORT(L == NULL, return STDL_ERR_MALLOC, "malloc");

    // tmp
    int* ipiv = malloc(ctx->ncsfs * sizeof(int));
    STDL_ERROR_HANDLE_AND_REPORT(ipiv == NULL, STDL_FREE_ALL(L); return STDL_ERR_MALLOC, "malloc");

    for (size_t iw = 0; iw < nw; ++iw) {
        // make left side: L = 2*(A-w*1)
        for (size_t kia = 0; kia < ctx->ncsfs; ++kia) {
            for(size_t kjb = 0; kjb <= kia; ++kjb)
                L[STDL_MATRIX_SP_IDX(kia, kjb)] = 2 * ctx->A[STDL_MATRIX_SP_IDX(kia, kjb)] - (kia == kjb ? 2 * w[iw] : 0);
        }

        float *Xi = X + iw * szXY;

        // copy egrad in X, to keep it for latter
        memcpy(Xi, egrad, ctx->ncsfs * ndim * sizeof(float ));

        // solve the problem
        err = LAPACKE_ssptrf(LAPACK_ROW_MAJOR, 'L', (int) ctx->ncsfs, L, ipiv);
        STDL_ERROR_HANDLE_AND_REPORT(err != 0,  STDL_FREE_ALL(L, ipiv); return STDL_ERR_RESPONSE, "error while ssptrf(): %d", err);

        err = LAPACKE_ssptrs(LAPACK_ROW_MAJOR, 'L', (int) ctx->ncsfs, (int) ndim, L, ipiv, Xi, (int) ndim);
        STDL_ERROR_HANDLE_AND_REPORT(err != 0, STDL_FREE_ALL(L, ipiv); return STDL_ERR_RESPONSE, "error while ssptrs(): %d", err);
    }

    STDL_FREE_ALL(L, ipiv);

    return STDL_ERR_OK;
}
