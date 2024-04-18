#include <assert.h>
#include <string.h>

#include "stdlite/response.h"
#include "stdlite/logging.h"
#include "stdlite/helpers.h"
#include "stdlite/utils/matrix.h"
#include "stdlite/linear_algebra.h"

void _log_memory(size_t sz) {
    double wrk_asz;
    char* wrk_usz;
    stdl_convert_size(sz, &wrk_asz, &wrk_usz);
    stdl_log_msg(0, "Extra memory required (temporary matrices): %.1f%s\n", wrk_asz, wrk_usz);
}

int stdl_response_TDA_casida(stdl_context *ctx, size_t nexci, float *e, float *X) {
    assert(ctx != NULL && ctx->ncsfs > 0 && nexci > 0 && nexci <= ctx->ncsfs && e != NULL && X != NULL);

    size_t wrk_sz = (STDL_MATRIX_SP_SIZE(ctx->ncsfs) + 8 * ctx->ncsfs) * sizeof(float ), iwrk_sz = 6 * ctx->ncsfs * sizeof(STDL_LA_INT);
    _log_memory(wrk_sz + (nexci < ctx->ncsfs? iwrk_sz : 0));

    stdl_log_msg(1, "+ ");
    stdl_log_msg(0, "Compute %ld excitation amplitude vectors (TDA-DFT) >", nexci);
    stdl_log_msg(1, "\n  | ");

    STDL_LA_INT lapack_err;

    float * wrk = malloc(wrk_sz);
    STDL_ERROR_HANDLE_AND_REPORT(wrk == NULL, return STDL_ERR_MALLOC, "malloc");

    memcpy(wrk, ctx->ApB, STDL_MATRIX_SP_SIZE(ctx->ncsfs) * sizeof(float ));
    float* lapack_wrk = wrk + STDL_MATRIX_SP_SIZE(ctx->ncsfs);

    if (nexci < ctx->ncsfs) {
        stdl_log_msg(1, "use sspevx ");

        STDL_LA_INT found = 0;
        STDL_LA_INT* iwrk = malloc(iwrk_sz);
        STDL_ERROR_HANDLE_AND_REPORT(iwrk == NULL, free(wrk); return STDL_ERR_MALLOC, "malloc");

        STDL_LA_INT* ifail = iwrk, *lapack_iwrk = iwrk + ctx->ncsfs;

        lapack_err = LAPACKE_sspevx_work(
                LAPACK_ROW_MAJOR, 'V', 'I', 'L',
                (STDL_LA_INT) ctx->ncsfs, wrk,
                .0f, .0f,
                1 /* even though we are in C, it starts at 1 */, (STDL_LA_INT) nexci, STDL_RESPONSE_EIGV_ABSTOL,
                &found, e, X, (STDL_LA_INT) nexci, lapack_wrk, lapack_iwrk, ifail
        );

        STDL_FREE_ALL(iwrk);
    } else {
        stdl_log_msg(1, "use sspev ");
        lapack_err = LAPACKE_sspev_work(
                LAPACK_ROW_MAJOR, 'V', 'L',
                (STDL_LA_INT) ctx->ncsfs, wrk,
                e, X, (STDL_LA_INT) ctx->ncsfs,
                lapack_wrk
        );
    }

    STDL_ERROR_HANDLE_AND_REPORT(lapack_err != 0, free(wrk); return STDL_ERR_RESPONSE, "error while sspevx(): %d", lapack_err);

    int err = stdl_matrix_sge_transpose(ctx->ncsfs, nexci, X);
    STDL_ERROR_CODE_HANDLE(err, free(wrk); return lapack_err);

    STDL_FREE_ALL(wrk);

    stdl_log_msg(0, "< done\n");

    stdl_matrix_sge_print(2, nexci, ctx->ncsfs, X, "Amplitude X");

    return STDL_ERR_OK;
}

int stdl_response_TD_casida(stdl_context *ctx, size_t nexci, float *e, float *X, float* Y) {
    assert(ctx != NULL && ctx->ncsfs > 0 && nexci <= ctx->ncsfs && ctx->AmB != NULL && e != NULL && X != NULL && Y != NULL);

    size_t sz = ctx->ncsfs * ctx->ncsfs;
    size_t wrk_sz = (3 * sz + 8 * ctx->ncsfs) * sizeof(float), iwrk_sz = 6 * ctx->ncsfs * sizeof(STDL_LA_INT);
    _log_memory(wrk_sz + (nexci < ctx->ncsfs? iwrk_sz : 0) + (ctx->ncsfs * (ctx->ncsfs + 4)) * sizeof(float ) /* <- stdl_matrix_ssp_sqrt_sy */);

    stdl_log_msg(1, "+ ");
    stdl_log_msg(0, "Compute %ld excitation amplitude vectors (TD-DFT) >", nexci);
    stdl_log_msg(1, "\n  | Make (A-B)^½ ");

    float * wrk = malloc(wrk_sz);
    STDL_ERROR_HANDLE_AND_REPORT(wrk == NULL, return STDL_ERR_MALLOC, "malloc");

    float *U = wrk, *V = wrk + sz, *W = wrk + 2 * sz, *lapack_wrk = wrk + 3 * sz;

    // U=A+B
    stdl_matrix_ssp_blowsy(ctx->ncsfs, 'L', ctx->ApB, U);

    // V=(A-B)^1/2
    stdl_matrix_ssp_sqrt_sy(ctx->ncsfs, ctx->AmB, V);

    stdl_log_msg(0, "-");
    stdl_log_msg(1, "\n  | Make (A-B)^½*(A+B)*(A-B)^½ ");

    // W = (A+B)*(A-B)^1/2 = U*V
    cblas_ssymm(CblasRowMajor, CblasRight, CblasLower,
                (STDL_LA_INT) ctx->ncsfs, (STDL_LA_INT)  ctx->ncsfs,
                1.f, U, (STDL_LA_INT)  ctx->ncsfs,
                V, (STDL_LA_INT) ctx->ncsfs,
                .0f, W, (STDL_LA_INT) ctx->ncsfs
    );

    // U = (A-B)^1/2*(A+B)*(A-B)^1/2 = V*W
    cblas_ssymm(CblasRowMajor, CblasRight, CblasLower,
                (STDL_LA_INT) ctx->ncsfs, (STDL_LA_INT)  ctx->ncsfs,
                1.f, V, (STDL_LA_INT)  ctx->ncsfs,
                W, (STDL_LA_INT) ctx->ncsfs,
                .0f, U, (STDL_LA_INT) ctx->ncsfs
    );

    stdl_log_msg(0, "-");
    stdl_log_msg(1, "\n  | ");

    STDL_LA_INT lapack_err;

    if (nexci < ctx->ncsfs) {
        stdl_log_msg(1, "use ssyevx ");

        STDL_LA_INT found = 0;
        STDL_LA_INT* iwrk = malloc(iwrk_sz);
        STDL_ERROR_HANDLE_AND_REPORT(iwrk == NULL, return STDL_ERR_MALLOC, "malloc");

        STDL_LA_INT* ifail = iwrk, *lapack_iwrk = iwrk + ctx->ncsfs;

        lapack_err = LAPACKE_ssyevx_work(
                LAPACK_ROW_MAJOR, 'V', 'I', 'L',
                (STDL_LA_INT) ctx->ncsfs, U, (STDL_LA_INT) ctx->ncsfs,
                .0f, .0f,
                1 /* even though we are in C, it starts at 1 */, (STDL_LA_INT) nexci, STDL_RESPONSE_EIGV_ABSTOL,
                &found, e, X, (STDL_LA_INT) nexci, lapack_wrk, 8 * (STDL_LA_INT) ctx->ncsfs, lapack_iwrk, ifail
        );

        STDL_FREE_ALL(ifail);
    } else {
        stdl_log_msg(1, "use ssyev ");
        lapack_err = LAPACKE_ssyev_work(
                LAPACK_ROW_MAJOR, 'V', 'L',
                (STDL_LA_INT) ctx->ncsfs, U, (STDL_LA_INT) ctx->ncsfs,
                e, lapack_wrk,  8 * (STDL_LA_INT) ctx->ncsfs
        );

        memcpy(X, U, sz * sizeof(float));
    }

    STDL_ERROR_HANDLE_AND_REPORT(lapack_err != 0, free(wrk); return STDL_ERR_RESPONSE, "error while ssyev[x](): %d", lapack_err);

    stdl_log_msg(0, "-");
    stdl_log_msg(1, "\n  | Transpose ");

    int err = stdl_matrix_sge_transpose(ctx->ncsfs, nexci, X); // Now, X' = (A-B)^(-1/2)*(X+Y).
    STDL_ERROR_CODE_HANDLE(err, free(wrk); return lapack_err);

    // stdl_matrix_sge_print(nexci, ctx->ncsfs, X, "Z");

    stdl_log_msg(0, "-");
    stdl_log_msg(1, "\n  | Compute X and Y ");

    // Get energies and compute X and Y
    #pragma omp parallel for
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
                sum += ctx->ApB[STDL_MATRIX_SP_IDX(kia, kjb)] * X[iexci * ctx->ncsfs + kjb] / e[iexci];
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

    stdl_log_msg(0, "< done\n");

    stdl_matrix_sge_print(2, nexci, ctx->ncsfs, X, "Amplitude X");
    stdl_matrix_sge_print(2, nexci, ctx->ncsfs, Y, "Amplitude Y");

    STDL_FREE_ALL(wrk);

    return STDL_ERR_OK;
}

int stdl_response_perturbed_gradient(stdl_context *ctx, size_t dim, int issym, double *op_ints_MO, float *egrad) {
    assert(ctx != NULL && dim > 0 && op_ints_MO != NULL && egrad != NULL);

    stdl_log_msg(1, "+ ");
    stdl_log_msg(0, "Compute perturbed gradient >");
    stdl_log_msg(1, "\n  | Looping through CSFs ");

    size_t nvirt = ctx->nmo - ctx->nocc;

    #pragma omp parallel for
    for (size_t lia = 0; lia < ctx->ncsfs; ++lia) {
        size_t i = ctx->csfs[lia] / nvirt, a = ctx->csfs[lia] % nvirt + ctx->nocc;
        for (size_t cpt = 0; cpt < dim; ++cpt) {
            egrad[lia * dim + cpt] = (!issym ? 2.f : -2.f) * (float) op_ints_MO[cpt * STDL_MATRIX_SP_SIZE(ctx->nmo) + STDL_MATRIX_SP_IDX(i, a)];
        }
    }

    stdl_log_msg(0, "< done\n");

    return STDL_ERR_OK;
}


int stdl_response_TD_linear(stdl_context *ctx, size_t nw, float *w, size_t ndim, int isherm, float *egrad, float *X, float *Y) {
    assert(ctx != NULL && ctx->ncsfs > 0 && nw > 0 && ndim > 0 && ctx->AmB != NULL && egrad != NULL && X != NULL && Y != NULL);

    size_t wrk_sz = (2 * STDL_MATRIX_SP_SIZE(ctx->ncsfs) + ctx->ncsfs) * sizeof(float), iwrk_sz = ctx->ncsfs * sizeof(STDL_LA_INT);
    _log_memory(wrk_sz + iwrk_sz);

    stdl_log_msg(1, "+ ");
    stdl_log_msg(0, "Compute %ld linear response vectors (TD-DFT) >", nw);
    stdl_log_msg(1, "\n  | Invert %s ", isherm ? "A-B" : "A+B");

    size_t szXY = ctx->ncsfs * ndim;

    float * wrk = malloc(wrk_sz);
    STDL_ERROR_HANDLE_AND_REPORT(wrk == NULL, return STDL_ERR_MALLOC, "malloc");

    float *AxB = wrk, *L = wrk + STDL_MATRIX_SP_SIZE(ctx->ncsfs), *lapack_wrk = wrk + 2 * STDL_MATRIX_SP_SIZE(ctx->ncsfs);

    // invert A-B [hermitian] or A+B [non-hermitian], taking advantage of its `sp` storage
    memcpy(AxB, isherm? ctx->AmB : ctx->ApB, STDL_MATRIX_SP_SIZE(ctx->ncsfs) * sizeof(float ));

    STDL_LA_INT* ipiv = malloc(iwrk_sz);
    STDL_ERROR_HANDLE_AND_REPORT(ipiv == NULL, STDL_FREE_ALL(wrk); return STDL_ERR_MALLOC, "malloc");

    STDL_LA_INT lapack_err = LAPACKE_ssptrf_work(LAPACK_ROW_MAJOR, 'L', (STDL_LA_INT) ctx->ncsfs, AxB, ipiv);
    STDL_ERROR_HANDLE_AND_REPORT(lapack_err != 0, STDL_FREE_ALL(wrk, ipiv); return STDL_ERR_RESPONSE, "error while ssptrf(): %d", lapack_err);

    lapack_err = LAPACKE_ssptri_work(LAPACK_ROW_MAJOR, 'L', (STDL_LA_INT) ctx->ncsfs, AxB, ipiv, lapack_wrk);
    STDL_ERROR_HANDLE_AND_REPORT(lapack_err != 0, STDL_FREE_ALL(wrk, ipiv); return STDL_ERR_RESPONSE, "error while ssptri(): %d", lapack_err);
    // now, AxB contains its inverse

    for (size_t iw = 0; iw < nw; ++iw) {
        stdl_log_msg(0, "-");
        stdl_log_msg(1, "\n  | Compute linear response vector for w=%f ", w[iw]);

        // make left side: L = (A+B)-w^2*(A-B)^(-1) [hermitian] (A-B)-w^2*(A+B)^(-1) [non-hermitian]
        #pragma omp parallel for schedule(guided)
        for (size_t kia = 0; kia < ctx->ncsfs; ++kia) {
            for(size_t kjb = 0; kjb <= kia; ++kjb)
                L[STDL_MATRIX_SP_IDX(kia, kjb)] = (isherm ? ctx->ApB : ctx->AmB)[STDL_MATRIX_SP_IDX(kia, kjb)] - powf(w[iw], 2) * AxB[STDL_MATRIX_SP_IDX(kia, kjb)];
        }

        float *Xi = X + iw * szXY, *Yi = Y + iw * szXY;

        // copy egrad in X, to keep it for latter
        memcpy(Xi, egrad, ctx->ncsfs * ndim * sizeof(float ));

        // solve the problem
        lapack_err = LAPACKE_ssptrf_work(LAPACK_ROW_MAJOR, 'L', (STDL_LA_INT) ctx->ncsfs, L, ipiv);
        STDL_ERROR_HANDLE_AND_REPORT(lapack_err != 0, STDL_FREE_ALL(wrk, ipiv); return STDL_ERR_RESPONSE, "error while ssptrf(): %d", lapack_err);

        lapack_err = LAPACKE_ssptrs_work(LAPACK_ROW_MAJOR, 'L', (STDL_LA_INT) ctx->ncsfs, (STDL_LA_INT) ndim, L, ipiv, Xi, (STDL_LA_INT) ndim);
        STDL_ERROR_HANDLE_AND_REPORT(lapack_err != 0, STDL_FREE_ALL(wrk, ipiv); return STDL_ERR_RESPONSE, "error while ssptrs(): %d", lapack_err);

        // separate X and Y
        // Xi' = Xi + Yi [hermitian] or Xi - Yi [non-hermitian]
        // Yi' = w*(A-B)^(-1)*Xi' = Xi - Yi [hermitian] or w*(A+B)^(-1)*Xi' = Xi + Yi [non-hermitian]
        #pragma omp parallel for
        for(size_t kia = 0; kia < ctx->ncsfs; kia++) {
            for(size_t cpt = 0; cpt < ndim; cpt++) {
                float sum = .0f;
                for (size_t kjb = 0; kjb < ctx->ncsfs; kjb++)
                    sum += AxB[STDL_MATRIX_SP_IDX(kia, kjb)] * w[iw] * Xi[kia * ndim + cpt];

                Yi[kia * ndim + cpt] = sum;
            }
        }

        // Xi = 1/2*(Xi' + Yi') && Yi = 1/2*(Xi' - Yi') [hermitian] or 1/2*(Yi-Xi) [non-hermitian]
        #pragma omp parallel for
        for(size_t kia = 0; kia < ctx->ncsfs; kia++) {
            for(size_t cpt = 0; cpt < ndim; cpt++) {
                float u = Xi[kia * ndim + cpt], v = Yi[kia * ndim + cpt];
                Xi[kia * ndim + cpt] = .5f * (u + v);
                Yi[kia * ndim + cpt] = .5f * (isherm ? (u - v) : (v - u));
            }
        }

        stdl_log_msg(2, "\n");
        stdl_matrix_sge_print(2, ctx->ncsfs, ndim, Xi, "x(w)");
        stdl_matrix_sge_print(2, ctx->ncsfs, ndim, Yi, "y(w)");
    }

    stdl_log_msg(0, "< done\n");

    STDL_FREE_ALL(wrk, ipiv);

    return STDL_ERR_OK;
}


int stdl_response_TDA_linear(stdl_context *ctx, size_t nw, float *w, size_t ndim, int isherm, float *egrad, float *X, float *Y) {
    assert(ctx != NULL && ctx->ncsfs > 0 && nw > 0 && ndim > 0 && egrad != NULL && X != NULL && Y != NULL);

    size_t wrk_sz = (2 * STDL_MATRIX_SP_SIZE(ctx->ncsfs) + ctx->ncsfs)* sizeof(float), iwrk_sz = ctx->ncsfs * sizeof(STDL_LA_INT);
    _log_memory(wrk_sz + iwrk_sz);

    stdl_log_msg(1, "+ ");
    stdl_log_msg(0, "Compute %ld linear response vectors (TDA-DFT) >", nw);
    stdl_log_msg(1, "\n  | Invert A ");

    size_t szXY = ctx->ncsfs * ndim;
    STDL_LA_INT lapack_err;

    // allocate space
    float* wrk = malloc(wrk_sz);
    STDL_ERROR_HANDLE_AND_REPORT(wrk == NULL, return STDL_ERR_MALLOC, "malloc");

    float* L = wrk, *Ai = wrk + STDL_MATRIX_SP_SIZE(ctx->ncsfs), *lapack_wrk = wrk + 2 *STDL_MATRIX_SP_SIZE(ctx->ncsfs);

    // invert A
    memcpy(Ai, ctx->ApB, STDL_MATRIX_SP_SIZE(ctx->ncsfs) * sizeof(float));

    STDL_LA_INT* ipiv = malloc(iwrk_sz);
    STDL_ERROR_HANDLE_AND_REPORT(ipiv == NULL, STDL_FREE_ALL(L); return STDL_ERR_MALLOC, "malloc");

    lapack_err = LAPACKE_ssptrf_work(LAPACK_ROW_MAJOR, 'L', (STDL_LA_INT) ctx->ncsfs, Ai, ipiv);
    STDL_ERROR_HANDLE_AND_REPORT(lapack_err != 0, STDL_FREE_ALL(L, ipiv, Ai); return STDL_ERR_RESPONSE, "error while ssptrf(): %d", lapack_err);

    lapack_err = LAPACKE_ssptri_work(LAPACK_ROW_MAJOR, 'L', (STDL_LA_INT) ctx->ncsfs, Ai, ipiv, lapack_wrk);
    STDL_ERROR_HANDLE_AND_REPORT(lapack_err != 0, STDL_FREE_ALL(L, ipiv, Ai); return STDL_ERR_RESPONSE, "error while ssptri(): %d", lapack_err);

    for (size_t iw = 0; iw < nw; ++iw) {
        stdl_log_msg(0, "-");
        stdl_log_msg(1, "\n  | Compute linear response vector for w=%f ", w[iw]);

        // make left side: L = A-w^2*A^(-1)
        #pragma omp parallel for schedule(guided)
        for (size_t kia = 0; kia < ctx->ncsfs; ++kia) {
            for(size_t kjb = 0; kjb <= kia; ++kjb)
                L[STDL_MATRIX_SP_IDX(kia, kjb)] = ctx->ApB[STDL_MATRIX_SP_IDX(kia, kjb)] - powf(w[iw], 2) * Ai[STDL_MATRIX_SP_IDX(kia, kjb)];
        }

        float *Xi = X + iw * szXY, *Yi = Y + iw * szXY;

        // copy egrad in Xi, to keep it for latter
        memcpy(Xi, egrad, ctx->ncsfs * ndim * sizeof(float ));

        // solve the problem
        lapack_err = LAPACKE_ssptrf_work(LAPACK_ROW_MAJOR, 'L', (STDL_LA_INT) ctx->ncsfs, L, ipiv);
        STDL_ERROR_HANDLE_AND_REPORT(lapack_err != 0, STDL_FREE_ALL(L, ipiv); return STDL_ERR_RESPONSE, "error while ssptrf(): %d", lapack_err);

        lapack_err = LAPACKE_ssptrs_work(LAPACK_ROW_MAJOR, 'L', (STDL_LA_INT) ctx->ncsfs, (STDL_LA_INT) ndim, L, ipiv, Xi, (STDL_LA_INT) ndim);
        STDL_ERROR_HANDLE_AND_REPORT(lapack_err != 0, STDL_FREE_ALL(L, ipiv); return STDL_ERR_RESPONSE, "error while ssptrs(): %d", lapack_err);

        // separate X and Y
        // Xi' = Xi + Yi [hermitian] or Xi - Yi [non-hermitian]
        // Yi' = w*A^(-1)*Xi' = Xi - Yi [hermitian] or Xi + Yi [non-hermitian]
        #pragma omp parallel for
        for(size_t kia = 0; kia < ctx->ncsfs; kia++) {
            for(size_t cpt = 0; cpt < ndim; cpt++) {
                float sum = .0f;
                for (size_t kjb = 0; kjb < ctx->ncsfs; kjb++)
                    sum += Ai[STDL_MATRIX_SP_IDX(kia, kjb)] * w[iw] * Xi[kia * ndim + cpt];

                Yi[kia * ndim + cpt] = sum;
            }
        }

        // Xi = 1/2*(Xi' + Yi') && Yi = 1/2*(Xi' - Yi') [hermitian] or 1/2*(Yi-Xi) [non-hermitian]
        #pragma omp parallel for
        for(size_t kia = 0; kia < ctx->ncsfs; kia++) {
            for(size_t cpt = 0; cpt < ndim; cpt++) {
                float u = Xi[kia * ndim + cpt], v = Yi[kia * ndim + cpt];
                Xi[kia * ndim + cpt] = .5f * (u + v);
                Yi[kia * ndim + cpt] = .5f * (isherm ? (u - v) : (v - u));
            }
        }

        stdl_log_msg(2, "\n");
        stdl_matrix_sge_print(2, ctx->ncsfs, ndim, Xi, "x(w)");
        stdl_matrix_sge_print(2, ctx->ncsfs, ndim, Xi, "y(w)");
    }

    stdl_log_msg(0, "< done\n");

    STDL_FREE_ALL(wrk, ipiv);

    return STDL_ERR_OK;
}
