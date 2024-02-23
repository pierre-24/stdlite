#include <stdio.h>
#include <assert.h>
#include <math.h>

#include "stdlite/response.h"
#include "stdlite/logging.h"
#include "stdlite/utils/matrix.h"
#include "stdlite/utils/permutations.h"


int stdl_property_polarizability(stdl_context* ctx, double* dips_MO, float* X, float* Y, float* alpha) {
    assert(ctx != NULL && dips_MO != NULL && X != NULL && Y != NULL && alpha != NULL);

    stdl_log_msg(0, "Compute polarizability tensor >");

    size_t nvirt = ctx->nmo - ctx->nocc;

    float s, d;

    for (int zeta = 0; zeta < 3; ++zeta) {
        for (int sigma = 0; sigma <= zeta; ++sigma) {
            alpha[STDL_MATRIX_SP_IDX(zeta, sigma)] = .0f;

            stdl_log_msg(0, "-");

            for (size_t lia = 0; lia < ctx->ncsfs; ++lia) {
                size_t i = ctx->csfs[lia] / nvirt, a = ctx->csfs[lia] % nvirt + ctx->nocc;

                d = (float) dips_MO[zeta * STDL_MATRIX_SP_SIZE(ctx->nmo) + STDL_MATRIX_SP_IDX(i, a)];
                s = X[lia * 3 + sigma];
                if(Y != NULL)
                    s += Y[lia * 3 + sigma];

                alpha[STDL_MATRIX_SP_IDX(zeta, sigma)] += -2.f * d * s;
            }
        }
    }

    stdl_log_msg(0, "< done\n");

    return STDL_ERR_OK;
}

int stdl_property_transition_dipoles(stdl_context *ctx, size_t nexci, double* dips_MO, float* X, float* Y, float * tdips) {
    assert(ctx != NULL && nexci > 0 && dips_MO != NULL && X != NULL && tdips != NULL);

    stdl_log_msg(0, "Compute ground to excited transition dipole moments >");

    size_t nvirt = ctx->nmo - ctx->nocc;
    float s2 = sqrtf(2);

    for (size_t iexci = 0; iexci < nexci; ++iexci) {
        tdips[0 * nexci + iexci] =  tdips[1 * nexci + iexci] =  tdips[2 * nexci + iexci] = .0f;

        for (size_t lia = 0; lia < ctx->ncsfs; ++lia) {
            size_t i = ctx->csfs[lia] / nvirt, a = ctx->csfs[lia] % nvirt + ctx->nocc;
            float amplitude = X[iexci * ctx->ncsfs + lia];
            if(Y != NULL)
                amplitude += Y[iexci * ctx->ncsfs + lia];

            for (size_t cpt = 0; cpt < 3; ++cpt)
                tdips[cpt * nexci + iexci] += s2 * amplitude * ((float) dips_MO[cpt * STDL_MATRIX_SP_SIZE(ctx->nmo) + STDL_MATRIX_SP_IDX(i, a)]);
        }
    }

    stdl_log_msg(0, "< done\n");

    return STDL_ERR_OK;
}

// a small structure to hold permutations in beta tensor calculations
typedef struct _bperm_ {
    float* X;
    float* Y;
    size_t cpt; // size_t is chosen for padding
} _bperm;

// Compute an element of the first hyperpolarizability tensor
int _first_hyperpolarizability_component(stdl_context* ctx, int component[3], double* dips_MO, float * X[3], float * Y[3], float* val) {
    assert(ctx != NULL && component != NULL && dips_MO != NULL && X[0] != NULL && X[1] != NULL && X[2] != NULL && Y[0] != NULL && Y[1] != NULL && Y[2] != NULL && val != NULL);
    size_t nvirt = ctx->nmo - ctx->nocc;

    stdl_permutations* set = NULL;

    int err = stdl_permutations_new(
            (_bperm []) {
                    {X[0], Y[0], component[0]},
                    {X[1], Y[1], component[1]},
                    {X[2], Y[2], component[2]}
            },
            3, sizeof(_bperm),
            &set
    );
    STDL_ERROR_CODE_HANDLE(err, return err);

    err = stdl_permutations_remove_duplicates(set, 3, sizeof(_bperm));
    STDL_ERROR_CODE_HANDLE(err, return err);

    stdl_permutations* current = set;
    size_t nperm = 0;

    *val = .0f;
    float A = .0f, B = .0f;

    while(current != NULL) {
        _bperm* e0 = (_bperm*) current->perm, *e1 = e0 + 1, *e2 = e0 + 2;

        float* cX = e0->X, *cY = e2->Y;
        size_t zeta = e0->cpt, sigma = e1->cpt, tau = e2->cpt;

        for (size_t lia = 0; lia < ctx->ncsfs; ++lia) {
            size_t i = ctx->csfs[lia] / nvirt, a = ctx->csfs[lia] % nvirt;
            float x = cX[lia * 3 + zeta];

            for (size_t ljb = 0; ljb < ctx->ncsfs; ++ljb) {
                size_t j = ctx->csfs[ljb] / nvirt, b = ctx->csfs[ljb] % nvirt;
                float y = cY[ljb * 3 + tau];

                if(b == a) { // jb == ja, so A'
                    float d = (float) -dips_MO[sigma * STDL_MATRIX_SP_SIZE(ctx->nmo) + STDL_MATRIX_SP_IDX(i, j)];
                    A += x * d * y;
                }

                if(j == i) {// jb == ib so B'
                    float d = (float) -dips_MO[sigma * STDL_MATRIX_SP_SIZE(ctx->nmo) + STDL_MATRIX_SP_IDX(ctx->nocc + a, ctx->nocc + b)];
                    B += x * d * y;
                }
            }
        }

        current = current->next;
        nperm++;
    }

    *val = (A - B) * (6 / (float) nperm);

    STDL_DEBUG("Computed element (%ld,%ld,%ld) = %f (nperm=%ld)", component[0], component[1], component[2], *val, nperm);

    err = stdl_permutations_delete(set);
    STDL_ERROR_CODE_HANDLE(err, return err);

    return STDL_ERR_OK;
}

int stdl_property_first_hyperpolarizability(stdl_context* ctx, double* dips_MO, float * Xs[3], float * Ys[3], float* beta) {
    assert(ctx != NULL && dips_MO != NULL && Xs[0] != NULL && Xs[1] != NULL && Xs[2] != NULL && Ys[0] != NULL && Ys[1] != NULL && Ys[2] != NULL && beta != NULL);

    stdl_log_msg(0, "Compute first hyperpolarizability tensor >");

    int isset[3][3][3] = {0};
    stdl_permutations* set = NULL;
    int err;

    for (int zeta = 0; zeta < 3; ++zeta) {
        for (int sigma = 0; sigma < 3; ++sigma) {
            for (int tau = 0; tau < 3; ++tau) {
                if(!isset[zeta][sigma][tau]) {
                    stdl_log_msg(0, "-");

                    float value = .0f;
                    _first_hyperpolarizability_component(ctx, (int[]) {zeta, sigma, tau}, dips_MO, Xs, Ys, &value);

                    // permute over the set of coordinates and use intrinsic permutations if any
                    err = stdl_permutations_new(
                            (_bperm []) {
                                    {Xs[0], Ys[0], zeta},
                                    {Xs[1], Ys[1], sigma},
                                    {Xs[2], Ys[2], tau}
                            },
                            3, sizeof(_bperm),
                            &set
                    );
                    STDL_ERROR_CODE_HANDLE(err, return err);

                    err = stdl_permutations_remove_duplicates(set, 3, sizeof(_bperm));
                    STDL_ERROR_CODE_HANDLE(err, return err);
                    stdl_permutations* current = set;
                    while(current != NULL) {
                        _bperm* e0 = (_bperm*) current->perm, *e1 = e0 + 1, *e2 = e0 + 2;
                        size_t rzeta = e0->cpt, rsigma = e1->cpt, rtau = e2->cpt;

                        if(e0->X == Xs[0] && e1->X == Xs[1] && e2->X == Xs[2]) {
                            beta[rzeta * 9 + rsigma * 3 + rtau] = value;
                            isset[rzeta][rsigma][rtau] = 1;

                            STDL_DEBUG("Set element (%ld,%ld,%ld) of beta via permutation", rzeta, rsigma, rtau);
                        }

                        current = current->next;
                    }

                    err = stdl_permutations_delete(set);
                    STDL_ERROR_CODE_HANDLE(err, return err);
                }
            }
        }
    }

    stdl_log_msg(0, "< done\n");

    return STDL_ERR_OK;
}
