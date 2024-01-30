#include <assert.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

#include "stdlite/context.h"
#include "stdlite/logging.h"
#include "stdlite/helpers.h"
#include "stdlite/matrix.h"

// Chemical hardness, from 10.1002/qua.22202
float eta[] = {
        .0, // no Z=0
        0.472592880,
        0.922033910,
        0.174528880,
        0.257007330,
        0.339490860,
        0.421954120,
        0.504381930,
        0.586918630,
        0.669313510,
        0.751916070,
        0.179641050,
        0.221572760,
        0.263485780,
        0.305396450,
        0.347340140,
        0.389247250,
        0.431156700,
        0.473082690,
        0.171054690,
        0.202762440,
        0.210073220,
        0.217396470,
        0.224710390,
        0.232015010,
        0.239339690,
        0.246656380,
        0.253982550,
        0.261288630,
        0.268594760,
        0.275925650,
        0.307629990,
        0.339315800,
        0.372359850,
        0.402735490,
        0.434457760,
        0.466117080,
        0.155850790,
        0.186493240,
        0.193562100,
        0.200633110,
        0.207705220,
        0.214772540,
        0.221846140,
        0.228918720,
        0.235986210,
        0.243056120,
        0.250130180,
        0.257199370,
        0.287847800,
        0.318486730,
        0.349124310,
        0.379765930,
        0.410408080,
        0.441057770,
        0.050193320,
        0.067625700,
        0.085044450,
        0.102477360,
        0.119911050,
        0.137327720,
        0.154762970,
        0.172182650,
        0.189612880,
        0.207047600,
        0.224467520,
        0.241896450,
        0.259325030,
        0.276760940,
        0.294182310,
        0.311595870,
        0.329022740,
        0.345922980,
        0.363880480,
        0.381305860,
        0.398774760,
        0.416142980,
        0.433645100,
        0.451040140,
        0.468489860,
        0.485845500,
        0.125267300,
        0.142686770,
        0.160116150,
        0.177558890,
        0.194975570,
        0.212407780,
        0.072635250,
        0.094221580,
        0.099202950,
        0.104186210,
        0.142356330,
        0.163942940,
        0.185519410,
        0.223701390,
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

    // Coulomb and exchange-like integrals (AA|BB)
    size_t natm = ctx->original_wf->natm;
    double* atm = ctx->original_wf->atm;
    float * AABB_J = malloc(natm * natm * sizeof(float));
    float * AABB_K = malloc(natm * natm * sizeof(float));
    STDL_ERROR_HANDLE_AND_REPORT(AABB_J == NULL || AABB_K == NULL, return STDL_ERR_MALLOC, "malloc");

    for(size_t A=0; A < natm; A++) {
        for(size_t B=0; B < natm; B++) {
            float r_AB = 0;
            if(A != B) {
                r_AB = (float) sqrt(pow(atm[A * 4 + 1] - atm[B * 4 + 1], 2) + pow(atm[A * 4 + 2] - atm[B * 4 + 2], 2) + pow(atm[A * 4 + 3] - atm[B * 4 + 3], 2));
            }

            float etaAB = .5f * (eta[(int) atm[A * 4 + 0]] + eta[(int) atm[B * 4 + 0]]);

            AABB_J[A * natm + B] = powf(1.f / ((float) powf(r_AB, ctx->gammaJ) + powf(ctx->ax * etaAB, -ctx->gammaJ)), 1.f / ctx->gammaJ);
            AABB_J[B * natm + A] = AABB_J[A * natm + B];

            AABB_K[A * natm + B] = powf(1.f / ((float) powf(r_AB, ctx->gammaK) + powf(etaAB, -ctx->gammaK)), 1.f / ctx->gammaK);
            AABB_K[B * natm + A] = AABB_K[A * natm + B];
        }
    }

    // stdl_matrix_sge_print(natm, 0, AABB_J, "(AA|BB)_J");
    // stdl_matrix_sge_print(natm, 0, AABB_K, "(AA|BB)_K");

    // density charges for Coulomb terms Q_ii and Q_aa
    size_t nvirt = ctx->nmo - ctx->nocc;
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

    // compute exchange term Q_ia
    size_t nexci = ctx->nocc * nvirt;

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

    STDL_FREE_ALL(AABB_J, AABB_K, qApp, qAia);

    return STDL_ERR_OK;
}
