#include <assert.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <cblas.h>

#include "stdlite/context.h"
#include "stdlite/logging.h"
#include "stdlite/helpers.h"
#include "stdlite/matrix.h"

// Chemical hardness (in a.u.), from 10.1002/qua.22202
float eta[] = {
        .0f, // no Z=0
        0.472592880f,
        0.922033910f,
        0.174528880f,
        0.257007330f,
        0.339490860f,
        0.421954120f,
        0.504381930f,
        0.586918630f,
        0.669313510f,
        0.751916070f,
        0.179641050f,
        0.221572760f,
        0.263485780f,
        0.305396450f,
        0.347340140f,
        0.389247250f,
        0.431156700f,
        0.473082690f,
        0.171054690f,
        0.202762440f,
        0.210073220f,
        0.217396470f,
        0.224710390f,
        0.232015010f,
        0.239339690f,
        0.246656380f,
        0.253982550f,
        0.261288630f,
        0.268594760f,
        0.275925650f,
        0.307629990f,
        0.339315800f,
        0.372359850f,
        0.402735490f,
        0.434457760f,
        0.466117080f,
        0.155850790f,
        0.186493240f,
        0.193562100f,
        0.200633110f,
        0.207705220f,
        0.214772540f,
        0.221846140f,
        0.228918720f,
        0.235986210f,
        0.243056120f,
        0.250130180f,
        0.257199370f,
        0.287847800f,
        0.318486730f,
        0.349124310f,
        0.379765930f,
        0.410408080f,
        0.441057770f,
        0.050193320f,
        0.067625700f,
        0.085044450f,
        0.102477360f,
        0.119911050f,
        0.137327720f,
        0.154762970f,
        0.172182650f,
        0.189612880f,
        0.207047600f,
        0.224467520f,
        0.241896450f,
        0.259325030f,
        0.276760940f,
        0.294182310f,
        0.311595870f,
        0.329022740f,
        0.345922980f,
        0.363880480f,
        0.381305860f,
        0.398774760f,
        0.416142980f,
        0.433645100f,
        0.451040140f,
        0.468489860f,
        0.485845500f,
        0.125267300f,
        0.142686770f,
        0.160116150f,
        0.177558890f,
        0.194975570f,
        0.212407780f,
        0.072635250f,
        0.094221580f,
        0.099202950f,
        0.104186210f,
        0.142356330f,
        0.163942940f,
        0.185519410f,
        0.223701390f,
};

int stdl_context_new(stdl_context ** ctx, stdl_wavefunction* wf, stdl_basis* bs, float gammaJ, float  gammaK, float ethr, float ax) {
    assert(ctx != NULL && wf != NULL && bs != NULL && gammaJ > 0 && gammaK > 0 && ethr > 0 && ax >= 0 && ax <= 1);

    *ctx = malloc(sizeof(stdl_context));
    STDL_ERROR_HANDLE_AND_REPORT(*ctx == NULL, return STDL_ERR_MALLOC, "malloc");

    (*ctx)->original_wf = wf;
    (*ctx)->bs = bs;
    (*ctx)->gammaJ = gammaJ;
    (*ctx)->gammaK = gammaK;
    (*ctx)->ethr = ethr;
    (*ctx)->ax = ax;

    // select MO to include
    STDL_DEBUG("range: %f a.u. (%.3f eV)", ethr, ethr * 27.212);

    size_t ohomo = (int) wf->nocc - 1, omin = 0, omax = 0;
    double ehomo = wf->e[ohomo], elumo = wf->e[ohomo + 1], ewin = 2 * (1 + .8 * ax)  * ethr, emin = elumo -ewin, emax = ehomo+ ewin;

    STDL_DEBUG("window: %f a.u. (%.3f eV)", ewin, ewin * 27.212);
    STDL_DEBUG("occ MO cutoff: %f a.u. (%.3f eV)", emin, emin * 27.212);
    STDL_DEBUG("virt MO cutoff: %f a.u. (%.3f eV)", emax, emax * 27.212);

    for(size_t i=0; i < wf->nmo; i++) {
        if(wf->e[i] >= emin && omin == 0)
            omin = (int) i;

        if(wf->e[i] <= emax)
            omax = (int) i;
        else
            break;
    }

    (*ctx)->nmo = omax - omin + 1;
    (*ctx)->nocc = ohomo - omin + 1;
    size_t nvirt = omax - ohomo;

    STDL_DEBUG("Resulting partition: [%d || %d | %d || %d] (occ + virt = %d, %.2f%% of %d MOs)", omin, (*ctx)->nocc, nvirt, wf->nmo - omax - 1, (*ctx)->nmo, (double) (*ctx)->nmo / wf->nmo * 100, wf->nmo);

    (*ctx)->e = malloc((*ctx)->nmo * sizeof(double ));
    (*ctx)->C = malloc((*ctx)->nmo * wf->nao * sizeof(double));

    STDL_ERROR_HANDLE_AND_REPORT((*ctx)->e == NULL || (*ctx)->C == NULL, stdl_context_delete(*ctx); return STDL_ERR_MALLOC, "malloc");

    // copy coefficients and energies
    for(size_t i=0; i < (*ctx)->nmo; i++) {
        (*ctx)->e[i] = wf->e[omin + i];
        memcpy(&((*ctx)->C[i * wf->nao]), &(wf->C[(i + omin) * wf->nao]), wf->nao * sizeof(double));
    }

    STDL_DEBUG("Orthogonalize MOs");

    int error = stdl_wavefunction_orthogonalize_C((*ctx)->C, wf->S, (*ctx)->nmo, wf->nao);
    STDL_ERROR_CODE_HANDLE(error, stdl_context_delete(*ctx); return error);

    return STDL_ERR_OK;
}

int stdl_context_delete(stdl_context* ctx) {
    assert(ctx != NULL);

    if(ctx->original_wf != NULL)
        stdl_wavefunction_delete(ctx->original_wf);

    if(ctx->bs != NULL)
        stdl_basis_delete(ctx->bs);

    STDL_FREE_ALL(ctx->e, ctx->C, ctx);

    return STDL_ERR_OK;
}

int stdl_context_select_csf(stdl_context *ctx) {

    /**
     * 1) To select primary CSFs, one needs to evaluate A_iaia, and so: (ii|aa) and (ia|ia).
     */

    size_t natm = ctx->original_wf->natm, nvirt = ctx->nmo - ctx->nocc, nexci = ctx->nocc * nvirt;
    double* atm = ctx->original_wf->atm;

    // Coulomb and exchange-like integrals (AA|BB)
    float * AABB_J = malloc(natm * natm * sizeof(float));
    float * AABB_K = malloc(natm * natm * sizeof(float));
    STDL_ERROR_HANDLE_AND_REPORT(AABB_J == NULL || AABB_K == NULL, return STDL_ERR_MALLOC, "malloc");

    for(size_t A=0; A < natm; A++) {
        for(size_t B=0; B <= A; B++) {
            float r_AB = 0;
            if(A != B) {
                r_AB = (float) sqrt(pow(atm[A * 4 + 1] - atm[B * 4 + 1], 2) + pow(atm[A * 4 + 2] - atm[B * 4 + 2], 2) + pow(atm[A * 4 + 3] - atm[B * 4 + 3], 2));
            }

            float etaAB = .5f * (eta[(int) atm[A * 4 + 0]] + eta[(int) atm[B * 4 + 0]]);

            AABB_J[A * natm + B] = 1.f / powf(powf(r_AB, ctx->gammaJ) + powf(ctx->ax * etaAB, -ctx->gammaJ), 1.f / ctx->gammaJ);
            AABB_J[B * natm + A] = AABB_J[A * natm + B];

            AABB_K[A * natm + B] = 1.f / powf(powf(r_AB, ctx->gammaK) + powf(etaAB, -ctx->gammaK), 1.f / ctx->gammaK);
            AABB_K[B * natm + A] = AABB_K[A * natm + B];
        }
    }

    // stdl_matrix_sge_print(natm, 0, AABB_J, "(AA|BB)_J");
    // stdl_matrix_sge_print(natm, 0, AABB_K, "(AA|BB)_K");

    // density charges for Coulomb terms: Q_ii and Q_aa
    float* qApp = malloc(ctx->nmo * natm * sizeof(float));
    STDL_ERROR_HANDLE_AND_REPORT(qApp == NULL, return STDL_ERR_MALLOC, "malloc");

    for(size_t p=0; p < ctx->nmo; p++) {
        for(size_t A=0; A < natm; A++)
            qApp[p * natm + A] = .0f;

        for(size_t mu=0; mu < ctx->original_wf->nao; mu++) {
            qApp[p * natm + ctx->original_wf->aotoatm[mu]] += powf((float) ctx->C[p * ctx->original_wf->nao + mu], 2.f);
        }
    }

    // use pointers, so that qAii is [nocc * natm] and qAaa is [nvirt * natm]
    float* qAii = qApp, *qAaa = qApp + ctx->nocc * natm;

    // stdl_matrix_sge_print(ctx->nocc, natm, qAii, "Q^A_ii");

    // stdl_matrix_sge_print(nvirt, natm, qAaa, "Q^A_aa");

    // tmpAii = qAii * AABB_J
    float* tmpAii = malloc(ctx->nocc * natm * sizeof(float));
    STDL_ERROR_HANDLE_AND_REPORT(tmpAii == NULL, return STDL_ERR_MALLOC, "malloc");

    cblas_ssymm(
            CblasRowMajor, CblasRight, CblasLower,
            (int) ctx->nocc, (int) natm,
            1.0f, AABB_J, (int) natm,
            qAii, (int) natm,
            .0f, tmpAii, (int) natm
            );

    // stdl_matrix_sge_print(ctx->nocc, natm, tmpAii, "tmp^A_ii");

    // (ii|aa) = tmpAii * qAaa^T
    float* iiaa = malloc(ctx->nocc * nvirt * sizeof(float));
    STDL_ERROR_HANDLE_AND_REPORT(iiaa == NULL, return STDL_ERR_MALLOC, "malloc");

    cblas_sgemm(
            CblasRowMajor, CblasNoTrans, CblasTrans,
            (int) ctx->nocc, (int) nvirt, (int) natm,
            1.0f, tmpAii, (int) natm,
            qAaa, (int) natm,
            .0f, iiaa, (int) nvirt
            );

    // stdl_matrix_sge_print(ctx->nocc, nvirt, iiaa, "(ii|aa)");

    STDL_FREE_ALL(tmpAii, qApp);

    // compute density charges for exchange term: Q_ia

    float* qAia = malloc( nexci * natm * sizeof(float));
    STDL_ERROR_HANDLE_AND_REPORT(qAia == NULL, return STDL_ERR_MALLOC, "malloc");

    for(size_t i=0; i < ctx->nocc; i++) {
        for (size_t a = 0; a < nvirt; ++a) {
            size_t k = i * nvirt + a;

            for(size_t A=0; A < natm; A++)
                qAia[k * natm + A] = .0f;

            for(size_t mu=0; mu < ctx->original_wf->nao; mu++) {
                qAia[k * natm + ctx->original_wf->aotoatm[mu]] += (float) (ctx->C[i * ctx->original_wf->nao + mu] * ctx->C[(ctx->nocc + a) * ctx->original_wf->nao + mu]);
            }
        }
    }

    // stdl_matrix_sge_print(nexci, natm, qAia, "Q^A_ia");

    // tmpAia = qAia * AABB_K
    float* tmpAia = malloc( nexci * natm * sizeof(float));
    STDL_ERROR_HANDLE_AND_REPORT(tmpAia == NULL, return STDL_ERR_MALLOC, "malloc");

    cblas_ssymm(
            CblasRowMajor, CblasRight, CblasLower,
            (int) nexci, (int) natm,
            1.0f, AABB_K, (int) natm,
            qAia, (int) natm,
            .0f, tmpAia, (int) natm
    );

    // stdl_matrix_sge_print(nexci, natm, tmpAia, "tmp^A_ia");

    // (ia|ia) = tmpAia * qAia^T
    float* iaia = malloc(nexci * nexci * sizeof(float));
    STDL_ERROR_HANDLE_AND_REPORT(iaia == NULL, return STDL_ERR_MALLOC, "malloc");

    cblas_sgemm(
            CblasRowMajor, CblasNoTrans, CblasTrans,
            (int) nexci, (int) nexci, (int) natm,
            1.0f, tmpAia, (int) natm,
            qAia, (int) natm,
            .0f, iaia, (int)nexci
    );

    // stdl_matrix_sge_print(nexci, nexci, iaia, "(ia|ia)");

    STDL_FREE_ALL(tmpAia, qAia);

    STDL_FREE_ALL(AABB_J, AABB_K, iiaa, iaia);

    return STDL_ERR_OK;
}
