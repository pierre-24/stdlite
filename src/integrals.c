#include "stdlite/integrals.h"
#include "stdlite/utils/matrix.h"

#include <math.h>

// compute renormalization factors for cartesian functions with l > 1.
void _compute_renormalization(stdl_basis* bs, double* renorm, double* buff) {
    int si, ioffset = 0;

    for(int ibas=0; ibas < bs->nbas; ibas++) {
        int angmom = bs->bas[ibas * 8 + 1];

        if(bs->use_spherical)
            si = CINTcgto_spheric(ibas, bs->bas);
        else
            si = CINTcgtos_cart(ibas, bs->bas);

        if(bs->use_spherical || angmom < 2) {
            for (int mu = 0; mu < si; ++mu)
                renorm[ioffset + mu] = 1.;
        } else if (angmom == 2) { // 6d
            double vxx = sqrt(1.25 / M_PI), vxy = sqrt(3.75 / M_PI);
            renorm[ioffset + 0] = vxx;
            renorm[ioffset + 1] = vxy;
            renorm[ioffset + 2] = vxy;
            renorm[ioffset + 3] = vxx;
            renorm[ioffset + 4] = vxy;
            renorm[ioffset + 5] = vxx;
        } else if (angmom == 3) { // 10f
            double vxxx = sqrt(1.75 / M_PI), vxxy = sqrt(8.75 / M_PI), vxyz = sqrt(26.25 / M_PI);
            renorm[ioffset + 0] = vxxx;
            renorm[ioffset + 1] = vxxy;
            renorm[ioffset + 2] = vxxy;
            renorm[ioffset + 3] = vxxy;
            renorm[ioffset + 4] = vxyz;
            renorm[ioffset + 5] = vxxy;
            renorm[ioffset + 6] = vxxx;
            renorm[ioffset + 7] = vxxy;
            renorm[ioffset + 8] = vxxy;
            renorm[ioffset + 9] = vxxx;
        } else {
            int1e_ovlp_cart(buff, NULL, (int[]) {ibas, ibas}, bs->atm, bs->natm, bs->bas, bs->nbas, bs->env, NULL, NULL);
            for (int mu = 0; mu < si; ++mu)
                renorm[ioffset + mu] = 1 / sqrt(buff[mu * si + mu]);
        }

        ioffset += si;
    }
}

int stdl_operator_int1e_dsp(stdl_basis *bs, stdl_operator op, double fac, double *values) {
    assert(bs != NULL && values != NULL && op < STDL_OP_COUNT);

    stdl_log_msg(0, "Computing 1e integral elements for `%s` (dim=%d) >", STDL_OPERATOR_NAME[op], STDL_OPERATOR_DIM[op]);

    size_t nao = 0;
    for(int ibas=0; ibas < bs->nbas; ibas++) {
        if (bs->use_spherical)
            nao += CINTcgto_spheric(ibas, bs->bas);
        else
            nao += CINTcgtos_cart(ibas, bs->bas);
    }

    int si, sj, ioffset=0, joffset;

    double buff[STDL_OPERATOR_DIM[op] * CART_MAX * CART_MAX];
    double* renorm = malloc(nao * sizeof(double));
    STDL_ERROR_HANDLE_AND_REPORT(renorm == NULL, return STDL_ERR_MALLOC, "malloc");

    _compute_renormalization(bs, renorm, buff);

    size_t ndips = STDL_MATRIX_SP_SIZE(nao);

    for(int ibas=0; ibas < bs->nbas; ibas++) {
        if(bs->use_spherical)
            si = CINTcgto_spheric(ibas, bs->bas);
        else
            si = CINTcgtos_cart(ibas, bs->bas);

        joffset = 0;

        for(int jbas=0; jbas <= ibas; jbas++) {
            if(bs->use_spherical) {
                sj = CINTcgto_spheric(jbas, bs->bas);
                STDL_OPERATOR_TO_CINT_SPH[op](buff, NULL, (int[]) {ibas, jbas}, bs->atm, bs->natm, bs->bas, bs->nbas, bs->env, NULL, NULL);
            }
            else {
                sj = CINTcgtos_cart(jbas, bs->bas);
                STDL_OPERATOR_TO_CINT_CART[op](buff, NULL, (int[]) {ibas, jbas}, bs->atm, bs->natm, bs->bas, bs->nbas, bs->env, NULL, NULL);
            }

            for(int iprim=0; iprim < si; iprim++) {
                for(int jprim=0; jprim < sj && joffset + jprim <= ioffset + iprim; jprim++) {
                    for (size_t idim = 0; idim < STDL_OPERATOR_DIM[op]; ++idim) {
                        values[idim * ndips + STDL_MATRIX_SP_IDX(ioffset + iprim, joffset + jprim)] = fac * buff[idim * si * sj + jprim * si + iprim] * renorm[ioffset + iprim] * renorm[joffset + jprim];
                    }
                }
            }

            joffset += sj;
        }

        ioffset += si;
    }

    free(renorm);

    stdl_log_msg(0, "< done\n");

    return STDL_ERR_OK;
}
