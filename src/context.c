#include <assert.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <cblas.h>

#include "stdlite/context.h"
#include "stdlite/logging.h"
#include "stdlite/helpers.h"
#include "stdlite/matrix.h"

// Chemical hardness (in Eh), from 10.1002/qua.22202
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

int stdl_context_new(stdl_context **ctx, stdl_wavefunction *wf, stdl_basis *bs, float gammaJ, float gammaK, float ethr,
                     float e2thr, float ax) {
    assert(ctx != NULL && wf != NULL && bs != NULL && gammaJ > 0 && gammaK > 0 && ethr > 0 && e2thr > 0 && ax >= 0 && ax <= 1);

    *ctx = malloc(sizeof(stdl_context));
    STDL_ERROR_HANDLE_AND_REPORT(*ctx == NULL, return STDL_ERR_MALLOC, "malloc");

    (*ctx)->original_wf = wf;
    (*ctx)->bs = bs;
    (*ctx)->gammaJ = gammaJ;
    (*ctx)->gammaK = gammaK;
    (*ctx)->ethr = ethr;
    (*ctx)->e2thr = e2thr;
    (*ctx)->ax = ax;

    // select MO to include
    STDL_DEBUG("range: %f Eh (%.3f eV)", ethr, ethr * 27.212);

    size_t ohomo = (int) wf->nocc - 1, omin = 0, omax = 0;
    double ehomo = wf->e[ohomo], elumo = wf->e[ohomo + 1], ewin = 2 * (1 + .8 * ax)  * ethr, emin = elumo -ewin, emax = ehomo+ ewin;

    STDL_DEBUG("window: %f Eh (%.3f eV)", ewin, ewin * 27.212);
    STDL_DEBUG("occ MO cutoff: %f Eh (%.3f eV)", emin, emin * 27.212);
    STDL_DEBUG("virt MO cutoff: %f Eh (%.3f eV)", emax, emax * 27.212);

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

size_t _lin(size_t i, size_t j) {
    return (i >= j) ? i*(i+1) / 2 + j : j*(j+1) / 2 + i;
}

int stdl_context_select_csf(stdl_context *ctx) {
    assert(ctx != NULL);


    size_t natm = ctx->original_wf->natm,
            nvirt = ctx->nmo - ctx->nocc,
            nexci_ij = STDL_MATRIX_SP_SIZE(ctx->nocc),
            nexci_ab = STDL_MATRIX_SP_SIZE(nvirt),
            nexci_ia = ctx->nocc * nvirt;

    /*
     * 1) Prepare charges and intermediates, as one big block of (continuous) memory.
     *
     *  <--- natm --------->
     * +--------------------+ 0
     * |  AABB_J            |
     * +--------------------+ natm
     * |  AABB_K            |
     * +--------------------+ natm
     * |  qAij [packed]     |
     * +--------------------+ nexci_ij
     * |  qAab [packed]     |
     * +--------------------+ nexci_ab
     * |  qAia              |
     * +--------------------+ nexci_ia
     * |  ijBB_J [packed]   |
     * +--------------------+ nexci_ij
     * |  iaBB_K            |
     * +--------------------+ nexci_ia
     *
     * TOTAL = natm * (2 * natm + 2 * nexci_ij + 2 * nexci_ia + nexci_ab)
     */

    double* atm = ctx->original_wf->atm;

    float* env = malloc(natm * (2 * natm + 2 * nexci_ij + 2 * nexci_ia + nexci_ab) * sizeof(float));
    STDL_ERROR_HANDLE_AND_REPORT(env == NULL, return STDL_ERR_MALLOC, "malloc");

    // Coulomb and exchange-like integrals (AA|BB), `float[natm * natm]`
    float * AABB_J = env;
    float * AABB_K = env + natm * natm;

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

    // 1. density charges for Coulomb terms: Q_A^ij [in packed form, float[STDL_SP_SIZE(nocc)]], Q_A^ab [in packed form, float[STDL_SP_SIZE(nocc)]], and Q_A^ia [float[nocc * nvirt]].
    float* qAij = env + (2 * natm) * natm;
    float* qAab = env + (2 * natm + nexci_ij) * natm;
    float* qAia = env + (2 * natm + nexci_ij + nexci_ab) * natm;

    for(size_t i=0; i < ctx->nocc; i++) {
        for(size_t j=0; j <= i; j++) {
            size_t k = i * (i + 1) / 2 + j;

            for(size_t A=0; A < natm; A++)
                qAij[k * natm + A] = .0f;

            for(size_t mu=0; mu < ctx->original_wf->nao; mu++) {
                qAij[k * natm + ctx->original_wf->aotoatm[mu]] += (float) (ctx->C[i * ctx->original_wf->nao + mu] * ctx->C[j * ctx->original_wf->nao + mu]);
            }
        }
    }

    // stdl_matrix_sge_print(nexci_ij, natm, qAij, "q_A^ij");

    for(size_t a=0; a < nvirt; a++) {
        for(size_t b=0; b <= a; b++) {
            size_t k = a * (a + 1) / 2 + b;

            for(size_t A=0; A < natm; A++)
                qAab[k * natm + A] = .0f;

            for(size_t mu=0; mu < ctx->original_wf->nao; mu++) {
                qAab[k * natm + ctx->original_wf->aotoatm[mu]] += (float) (ctx->C[(a + ctx->nocc) * ctx->original_wf->nao + mu] * ctx->C[(b + ctx->nocc) * ctx->original_wf->nao + mu]);
            }
        }
    }

    // stdl_matrix_sge_print(nexci_ab, natm, qAab, "q^A_ab");

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

    // stdl_matrix_sge_print(nexci_ia, natm, qAia, "Q_A^ia");

    // 2. Intermediates: (ij|BB)_J [in packed form], and (ia|BB)_K:
    float* ijBB_J = env + (2 * natm + nexci_ij + nexci_ab + nexci_ia) * natm;
    float* iaBB_K = env + (2 * natm + 2 * nexci_ij + nexci_ab + nexci_ia) * natm;

    cblas_ssymm(
            CblasRowMajor, CblasRight, CblasLower,
            (int) nexci_ij, (int) natm,
            1.0f, AABB_J, (int) natm,
            qAij, (int) natm,
            .0f, ijBB_J, (int) natm
    );

    // stdl_matrix_sge_print(nexci_ij, natm, ijBB_J, "(ij|BB)_J");

    cblas_ssymm(
            CblasRowMajor, CblasRight, CblasLower,
            (int) nexci_ia, (int) natm,
            1.0f, AABB_K, (int) natm,
            qAia, (int) natm,
            .0f, iaBB_K, (int) natm
    );

    // stdl_matrix_sge_print(nexci_ia, natm, iaBB_K, "(ia|BB)_K");

    /*
     *  2) To select primary CSFs i→a, one needs to evaluate A'_ia,ia = (e_a - e_i) + 2*(ia|ia)' - (ii|aa)'.
     *     Then, CSFs are selected if A'_ia,ia <= E_thr.
     */

    // marks csfs as not-included (0), primary (1), or secondary (2).
    char* csfs = malloc(nexci_ia * sizeof(short));
    STDL_ERROR_HANDLE_AND_REPORT(csfs == NULL, free(env); return STDL_ERR_MALLOC, "malloc");

    // store diagonal components
    float* A_diag = malloc(nexci_ia * sizeof(float ));
    STDL_ERROR_HANDLE_AND_REPORT(csfs == NULL, STDL_FREE_ALL(env, csfs); return STDL_ERR_MALLOC, "malloc");

    size_t ncsfs = 0;

    for(size_t i=0; i < ctx->nocc; i++) {
        size_t kii = STDL_MATRIX_SP_IDX(i, i);

        for (size_t a = 0; a < nvirt; ++a) {
            size_t kaa = STDL_MATRIX_SP_IDX(a, a);
            size_t kia = i * nvirt + a;
            float iaia = .0f;
            float iiaa = .0f;

            for(size_t A=0; A < natm; A++) { // scalar products to compute (ia|ia)' and (ii|aa)'.
                iaia += iaBB_K[kia * natm + A] * qAia[kia * natm + A];
                iiaa += ijBB_J[kii * natm + A] * qAab[kaa * natm + A];
            }

            A_diag[kia] = (float) (ctx->e[ctx->nocc + a] - ctx->e[i]) + 2 * iaia - iiaa;

            if(A_diag[kia] <= ctx->ethr) {
                csfs[kia] = 1;
                ncsfs++;
                STDL_DEBUG("selected primary:: %ld→%ld (E=%f Eh)", i, ctx->nocc + a, A_diag[kia]);
            } else {// ... the rest is selected to be considered in perturbation.
                csfs[kia] = 2; // mark as potentially secondary for the moment
            }
        }
    }

    // while we're at it, sort the CSFs
    size_t* csfs_sorted_indices = malloc(nexci_ia * sizeof(size_t));
    STDL_ERROR_HANDLE_AND_REPORT(csfs_sorted_indices == NULL, STDL_FREE_ALL(env, csfs, A_diag); return STDL_ERR_MALLOC, "malloc");

    for(size_t kia=0; kia < nexci_ia; kia++)
        csfs_sorted_indices[kia] = kia;

    // heap sort, https://en.wikipedia.org/wiki/Heapsort#Standard_implementation
    size_t start = nexci_ia / 2, end = nexci_ia, tmp_swap, root, child;
    while(end > 1) {
        if(start) {
            start -= 1;
        } else {
            end -= 1;
            tmp_swap = csfs_sorted_indices[end];
            csfs_sorted_indices[end] = csfs_sorted_indices[0];
            csfs_sorted_indices[0] = tmp_swap;
        }

        root = start;
        while((2 * root + 1) < end) {
            child = 2 * root + 1;
            if(child + 1 < end && A_diag[csfs_sorted_indices[child]] < A_diag[csfs_sorted_indices[child + 1]])
                child += 1;

            if(A_diag[csfs_sorted_indices[root]] < A_diag[csfs_sorted_indices[child]]) {
                tmp_swap = csfs_sorted_indices[child];
                csfs_sorted_indices[child] = csfs_sorted_indices[root];
                csfs_sorted_indices[root] = tmp_swap;
                root = child;
            } else {
                break;
            }
        }
    }

    if(ncsfs > 0) {
        /*
         * 3) Now, select S-CSFs j→b so that E^(2)_jb > E^(2)_thr.
         */

        for(size_t kjb=0; kjb < nexci_ia; kjb++) { // loop over possible S-CSFs
            if(csfs[kjb] == 2) {
                size_t b = kjb % nvirt, j = kjb / nvirt;
                float e2 = .0f; // perturbation energy

                for(size_t kia=0; kia < nexci_ia; kia++) { // loop over P-CSFs
                    if(csfs[kia] == 1) {
                        float iajb = .0f;
                        float ijab = .0f;

                        size_t a = kia % nvirt, i = kia / nvirt;
                        size_t kij = _lin(i, j);
                        size_t kab = _lin(a, b);

                        for(size_t A=0; A < natm; A++) { // scalar products to compute (ia|jb)' and (ij|ab)'.
                            iajb += iaBB_K[kia * natm + A] * qAia[kjb * natm + A];
                            ijab += ijBB_J[kij * natm + A] * qAab[kab * natm + A];
                        }

                        float A_iajb = 2 * iajb - ijab;

                        e2 += powf(A_iajb, 2) / (A_diag[kjb] - A_diag[kia]);
                    }
                }

                if(e2 < ctx->e2thr) {
                    csfs[kjb] = 0; // discarded
                } else {
                    STDL_DEBUG("selected secondary:: %ld→%ld (E=%f Eh)", j, ctx->nocc + b, A_diag[kjb]);
                    ncsfs++;
                }
            }
        }
    } else {
        STDL_WARN("no CSFs selected. `E_thr` should be at least %f Eh!", A_diag[csfs_sorted_indices[0]]);
    }

    STDL_DEBUG("selected %d CSFs (%.2f%% of %d CSFs)", ncsfs, (float) ncsfs / (float) nexci_ia * 100, nexci_ia);

    STDL_FREE_ALL(env, csfs, A_diag, csfs_sorted_indices);

    return STDL_ERR_OK;
}
