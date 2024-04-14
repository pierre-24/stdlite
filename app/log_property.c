#include <stdlite/logging.h>
#include <stdlite/helpers.h>
#include <stdlite/utils/experimental_quantity.h>
#include <assert.h>

#include "log_property.h"


int stdl_log_property_polarizability(stdl_response_request* req, float* alpha, float w) {
    assert(req != NULL && alpha != NULL);

    stdl_log_msg(0, "** alpha(-w;w) @ w=%f (%.2f nm)\n", w, STDL_CONST_HC / w);
    stdl_log_msg(0, "         x            y            z\n");

    for (size_t zeta = 0; zeta < 3; ++zeta) {
        switch (zeta) {
            case 0:
                stdl_log_msg(0, "x");
                break;
            case 1:
                stdl_log_msg(0, "y");
                break;
            case 2:
                stdl_log_msg(0, "z");
                break;
        }
        for (size_t sigma = 0; sigma < 3; ++sigma) {
            stdl_log_msg(0, " % 12.5f", alpha[zeta * 3 + sigma]);
        }
        stdl_log_msg(0, "\n");
    }

    float iso, aniso;
    stdl_qexp_polarizability(alpha, &iso, &aniso);
    stdl_log_msg(0, "iso =   % 12.5f\naniso = % 12.5f\n", iso, aniso);

    return STDL_ERR_OK;
}

int stdl_log_property_linear_tensor(stdl_response_request* req, float* tensor, float w) {
    assert(req != NULL && tensor != NULL);

    size_t dim0 = STDL_OPERATOR_DIM[req->ops[0]], dim1 = STDL_OPERATOR_DIM[req->ops[1]];
    float trs = STDL_OPERATOR_TRS[req->ops[0]] * STDL_OPERATOR_TRS[req->ops[1]];

    stdl_log_msg(0, "** -%s<<%s;%s>>_w @ w=%f Eh (%.2f nm)\n", trs < 0? "Im": "", STDL_OPERATOR_NAME[req->ops[0]], STDL_OPERATOR_NAME[req->ops[1]], w, STDL_CONST_HC / w);

    stdl_log_msg(0, "    ");
    for (size_t sigma = 0; sigma < dim1; ++sigma)
        stdl_log_msg(0, "     %3d     ", sigma);

    stdl_log_msg(0, "\n");

    for (size_t zeta = 0; zeta < dim0; ++zeta) {
        stdl_log_msg(0, "%-3d ", zeta);
        for (size_t sigma = 0; sigma < dim1; ++sigma) {
            stdl_log_msg(0, " % 12.5f", tensor[zeta * dim0 + sigma]);
        }
        stdl_log_msg(0, "\n");
    }

    return STDL_ERR_OK;
}


/*int stdl_log_property_first_hyperpolarizability(stdl_response_request* req, float beta[3][3][3]) {
    assert(req != NULL && beta != NULL);

    stdl_log_msg(0,
                 "** beta(-w1-w2;w1,w2), w1=%f (%.2f nm), w2=%f (%.2f nm)\n",
                 req->lrvreqs[1]->w[req->wpos[1]], STDL_CONST_HC / req->lrvreqs[1]->w[req->wpos[1]],
                 req->lrvreqs[2]->w[req->wpos[2]], STDL_CONST_HC / req->lrvreqs[2]->w[req->wpos[2]]
    );

    stdl_log_msg(0, "          x            y            z\n");
    for (size_t zeta = 0; zeta < 3; ++zeta) {
        for (size_t sigma = 0; sigma < 3; ++sigma) {
            switch (zeta) {
                case 0:
                    stdl_log_msg(0, "x");
                    break;
                case 1:
                    stdl_log_msg(0, "y");
                    break;
                case 2:
                    stdl_log_msg(0, "z");
                    break;
            }
            switch (sigma) {
                case 0:
                    stdl_log_msg(0, "x");
                    break;
                case 1:
                    stdl_log_msg(0, "y");
                    break;
                case 2:
                    stdl_log_msg(0, "z");
                    break;
            }

            for (size_t tau = 0; tau < 3; ++tau) {
                stdl_log_msg(0, " % 12.5f", beta[zeta][sigma][tau]);
            }

            stdl_log_msg(0, "\n");
        }
    }

    if(stdl_float_equals(req->lrvreqs[1]->w[req->wpos[1]], req->lrvreqs[2]->w[req->wpos[2]], 1e-6)) {
        float b2ZZZ, b2ZXX;
        stdl_qexp_first_hyperpolarizability_hrs(beta, &b2ZZZ, &b2ZXX);
        stdl_log_msg(0, "<B2ZZZ> = % 12.5f\n<B2ZXX> = % 12.5f\nBHRS    = % 12.5f\nDR      = % 12.5f\n", b2ZZZ, b2ZXX, sqrtf(b2ZZZ +  b2ZXX), b2ZZZ / b2ZXX);

    }

    return STDL_ERR_OK;
}*/

void stdl_log_property_amplitude_contributions(stdl_responses_handler *rh, stdl_context *ctx, float thresh) {
    assert(rh != NULL && ctx != NULL && thresh > 0);

    float s2o2 = sqrtf(2) / 2;
    size_t nvirt = ctx->nmo - ctx->nocc;

    stdl_log_msg(1, "**   -- E --- ---- Contributions (H=%4ld) ---\n", ctx->nocc);
    for (size_t iexci = 0; iexci < rh->nexci; ++iexci) {
        // print energies
        stdl_log_msg(1, "%4ld %8.5f", iexci + 1, rh->eexci[iexci]);

        // print contributions
        size_t printed = 0;
        for (size_t kia = 0; kia < ctx->ncsfs; ++kia) {
            float c = powf(rh->Xamp[iexci * ctx->ncsfs + kia], 2);

            if(rh->Yamp != NULL)
                c -= powf(rh->Yamp[iexci * ctx->ncsfs + kia], 2);

            if(c >= thresh) {
                size_t i = ctx->csfs[kia] / nvirt, a = ctx->csfs[kia] % nvirt, hi = ctx->nocc - i - 1;

                if(printed > 0)
                    stdl_log_msg(1, "             ");

                stdl_log_msg(1, " %5.1f%% (% 6.4f) H", c * 100, rh->Xamp[iexci * ctx->ncsfs + kia] * s2o2);

                if(hi > 0)
                    stdl_log_msg(1, "-%ld", hi);
                stdl_log_msg(1, "→L");
                if (a > 0)
                    stdl_log_msg(1, "+%ld", a);

                stdl_log_msg(1, " (%ld→%ld)", i + 1, ctx->nocc + a + 1);

                stdl_log_msg(1, "\n");

                printed += 1;
            }
        }

        if(iexci < rh->nexci -1)
            stdl_log_msg(1, "     -------- -------------------------------\n");
    }
}

int stdl_log_property_g2e_moments(stdl_responses_handler *rh, stdl_context *ctx, stdl_operator ops[2], float *tg2e) {
    assert(rh != NULL && ctx != NULL && tg2e != NULL);

    size_t dim0 = STDL_OPERATOR_DIM[ops[0]], dim1 = STDL_OPERATOR_DIM[ops[1]];

    // headers
    stdl_log_msg(0, "**    -------- Energy ------- ");
    for (size_t cpt = 0; cpt < dim0 * 2; ++cpt)
        stdl_log_msg(0, "-");
    stdl_log_msg(0, " <0|%4s|m> ", STDL_OPERATOR_NAME[ops[0]]);
    for (size_t cpt = 0; cpt < dim0 * 7 - 13; ++cpt)
        stdl_log_msg(0, "-");

    if(ops[0] != ops[1]) {
        stdl_log_msg(0, " ");
        for (size_t cpt = 0; cpt < dim0 * 2; ++cpt)
            stdl_log_msg(0, "-");
        stdl_log_msg(0, " <0|%4s|m> ", STDL_OPERATOR_NAME[ops[1]]);
        for (size_t cpt = 0; cpt < dim0 * 7 - 13; ++cpt)
            stdl_log_msg(0, "-");
    }

    if(dim0 == dim1) {
        stdl_log_msg(0, " --------");
        if ((ops[0] == ops[1] && ops[0] == STDL_OP_DIPL) // fL
            || (ops[0] == ops[1] && ops[0] == STDL_OP_DIPV) // fV
            || ((ops[0] == STDL_OP_DIPL && ops[1] == STDL_OP_ANGM) || (ops[1] == STDL_OP_DIPL && ops[0] == STDL_OP_ANGM)) // RL
            || ((ops[0] == STDL_OP_DIPV && ops[1] == STDL_OP_ANGM) || (ops[1] == STDL_OP_DIPV && ops[0] == STDL_OP_ANGM)) // RV
            )
            stdl_log_msg(0, "---------");
    }

    stdl_log_msg(0, "\n");

    stdl_log_msg(0, "       (Eh)     (eV)    (nm) ");
    for (size_t cpt = 0; cpt < dim0; ++cpt)
        stdl_log_msg(0, "   %-3d   ", cpt);

    if (ops[0] != ops[1]) {
        for (size_t cpt = 0; cpt < dim1; ++cpt)
            stdl_log_msg(0, "   %-3d   ", cpt);
    }

    if(dim0 == dim1) {
        stdl_log_msg(0, "    A·B  ");

        if (ops[0] == ops[1] && ops[0] == STDL_OP_DIPL)
            stdl_log_msg(0, "    fL   ");
        else if (ops[0] == ops[1] && ops[0] == STDL_OP_DIPV)
            stdl_log_msg(0, "    fV   ");
        else if ((ops[0] == STDL_OP_DIPL && ops[1] == STDL_OP_ANGM) || (ops[1] == STDL_OP_DIPL && ops[0] == STDL_OP_ANGM))
            stdl_log_msg(0, "    RL   ");
        else if ((ops[0] == STDL_OP_DIPV && ops[1] == STDL_OP_ANGM) || (ops[1] == STDL_OP_DIPV && ops[0] == STDL_OP_ANGM))
            stdl_log_msg(0, "    RV   ");
    }

    stdl_log_msg(0, "\n");

    for (size_t iexci = 0; iexci < rh->nexci; ++iexci) {
        // print energies
        stdl_log_msg(0, "%4ld %8.5f %7.3f %7.2f", iexci + 1, rh->eexci[iexci], rh->eexci[iexci] * STDL_CONST_AU_TO_EV, STDL_CONST_HC / rh->eexci[iexci]);

        for (size_t cpt = 0; cpt < dim0; ++cpt) {
            stdl_log_msg(0, " % 8.5f",tg2e[cpt * rh->nexci + iexci]);
        }

        if (ops[0] != ops[1]) {
            for (size_t cpt = 0; cpt < dim1; ++cpt) {
                stdl_log_msg(0, " % 8.5f",tg2e[rh->nexci * dim0 + cpt * rh->nexci + iexci]);
            }

        }

        if(dim0 == dim1) {
            float dotp = .0f;
            for (size_t cpt = 0; cpt < dim1; ++cpt) {
                dotp += tg2e[cpt * rh->nexci + iexci] * tg2e[rh->nexci * dim0 + cpt * rh->nexci + iexci];
            }

            stdl_log_msg(0, " % 8.5f", dotp);

            if(ops[0] == ops[1] && ops[0] == STDL_OP_DIPL)
                stdl_log_msg(0, " % 8.5f", 2.f/3 * rh->eexci[iexci] * dotp);
            else if(ops[0] == ops[1] && ops[0] == STDL_OP_DIPV)
                stdl_log_msg(0, " % 8.5f", 2.f/3 / rh->eexci[iexci] * dotp);
            else if ((ops[0] == STDL_OP_DIPL && ops[1] == STDL_OP_ANGM) || (ops[1] == STDL_OP_DIPL && ops[0] == STDL_OP_ANGM))
                stdl_log_msg(0, " % 8.5f", -.5* dotp);
            else if ((ops[0] == STDL_OP_DIPV && ops[1] == STDL_OP_ANGM) || (ops[1] == STDL_OP_DIPV && ops[0] == STDL_OP_ANGM))
                stdl_log_msg(0, " % 8.5f", -.5 * 1 / rh->eexci[iexci] * dotp);

        }

        stdl_log_msg(0, "\n");
    }

    return STDL_ERR_OK;
}

