#include <stdlite/logging.h>
#include <stdlite/helpers.h>
#include <stdlite/utils/experimental_quantity.h>
#include <stdlite/utils/matrix.h>
#include <assert.h>

#include "log_property.h"


int stdl_log_property_polarizability(stdl_response_request* req, float* alpha, float w) {
    assert(req != NULL && alpha != NULL);

    stdl_log_msg(0, "** alpha(-wB;wB) @ wB=%f Eh (%.2f nm)\n", w, STDL_CONST_HC / w);
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
    stdl_log_msg(0, "# iso   = % 12.5f\n# aniso = % 12.5f\n", iso, aniso);

    return STDL_ERR_OK;
}

int stdl_log_property_linear_tensor(stdl_response_request* req, float* tensor, float wB) {
    assert(req != NULL && tensor != NULL);

    size_t dim0 = STDL_OPERATOR_DIM[req->ops[0]], dim1 = STDL_OPERATOR_DIM[req->ops[1]];
    float trs = STDL_OPERATOR_TRS[req->ops[0]] * STDL_OPERATOR_TRS[req->ops[1]];

    stdl_log_msg(0, "** -%s<<%s;%s>>_wB @ wB=%f Eh (%.2f nm)\n", trs < 0? "Im": "", STDL_OPERATOR_NAME[req->ops[0]], STDL_OPERATOR_NAME[req->ops[1]], wB, STDL_CONST_HC / wB);

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

int stdl_log_property_first_hyperpolarizability(stdl_response_request* req, float* beta, float wB, float wC) {
    assert(req != NULL && beta != NULL);

    size_t dim0 = STDL_OPERATOR_DIM[req->ops[0]], dim1 = STDL_OPERATOR_DIM[req->ops[1]], dim2 = STDL_OPERATOR_DIM[req->ops[2]];
    float trs = STDL_OPERATOR_TRS[req->ops[0]] * STDL_OPERATOR_TRS[req->ops[1]] * STDL_OPERATOR_TRS[req->ops[2]];

    stdl_log_msg(0, "** beta(-wB-wC;wB,wC) @ wB=%f Eh (%.2f nm), wC=%f Eh (%.2f nm)\n", trs < 0? "Im": "", STDL_OPERATOR_NAME[req->ops[0]], STDL_OPERATOR_NAME[req->ops[1]], STDL_OPERATOR_NAME[req->ops[2]], wB, STDL_CONST_HC / wB, wC, STDL_CONST_HC / wC);

    stdl_log_msg(0, "          x            y            z\n");

    for (size_t zeta = 0; zeta < dim0; ++zeta) {
        for (size_t sigma = 0; sigma < dim1; ++sigma) {
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
            for (size_t tau = 0; tau < dim2; ++tau) {
                stdl_log_msg(0, " % 12.5f", beta[zeta * dim0 * dim1 + sigma * dim1 + tau]);
            }
            stdl_log_msg(0, "\n");
        }
    }

    float beta_vec[3] = {.0f};
    stdl_qexp_first_hyperpolarizability(beta, beta_vec);
    stdl_log_msg(0, "# Bvec    = (%.4f, %.4f, %.4f)\n# |Bvec|  = %12.5f\n",
                 beta_vec[0], beta_vec[1], beta_vec[2],
                 sqrtf(powf(beta_vec[0], 2) + powf(beta_vec[1], 2) + powf(beta_vec[2], 2))
                 );

    if(stdl_float_equals(wB, wC, 1e-6f)) { // SHG hyperpolarizability
        float b2ZZZ, b2ZXX, b2J1, b2J3;
        stdl_qexp_first_hyperpolarizability_hrs(beta, &b2ZZZ, &b2ZXX, &b2J1, &b2J3);
        stdl_log_msg(0,
                     "# <B2ZZZ> = % 12.5f\n"
                     "# <B2ZXX> = % 12.5f\n"
                     "# BHRS    = % 12.5f\n"
                     "# DR      = % 12.5f\n"
                     "# |Bj=1|  = % 12.5f\n"
                     "# |Bj=3|  = % 12.5f\n"
                     "# rho     = % 12.5f\n",
                     b2ZZZ,
                     b2ZXX,
                     sqrtf(b2ZZZ +  b2ZXX),
                     b2ZZZ / b2ZXX,
                     sqrtf(b2J1),
                     sqrtf(b2J3),
                     sqrtf(b2J3 / b2J1)
                     );

    }

    return STDL_ERR_OK;
}

int stdl_log_property_quadratic_tensor(stdl_response_request* req, float* tensor, float wB, float wC) {
    assert(req != NULL && tensor != NULL);

    size_t dim0 = STDL_OPERATOR_DIM[req->ops[0]], dim1 = STDL_OPERATOR_DIM[req->ops[1]], dim2 = STDL_OPERATOR_DIM[req->ops[2]];
    float trs = STDL_OPERATOR_TRS[req->ops[0]] * STDL_OPERATOR_TRS[req->ops[1]] * STDL_OPERATOR_TRS[req->ops[2]];

    stdl_log_msg(0, "** -%s<<%s;%s,%s>>_wB,wC @ wB=%f Eh (%.2f nm), wC=%f Eh (%.2f nm)\n", trs < 0? "Im": "", STDL_OPERATOR_NAME[req->ops[0]], STDL_OPERATOR_NAME[req->ops[1]], STDL_OPERATOR_NAME[req->ops[2]], wB, STDL_CONST_HC / wB, wC, STDL_CONST_HC / wC);

    stdl_log_msg(0, "        ");
    for (size_t tau = 0; tau < dim2; ++tau)
        stdl_log_msg(0, "     %3d     ", tau);

    stdl_log_msg(0, "\n");

    for (size_t zeta = 0; zeta < dim0; ++zeta) {
        for (size_t sigma = 0; sigma < dim1; ++sigma) {
            stdl_log_msg(0, "%3d,%-3d ", zeta, sigma);
            for (size_t tau = 0; tau < dim2; ++tau) {
                stdl_log_msg(0, " % 12.5f", tensor[zeta * dim0 * dim1 + sigma *dim1 + tau]);
            }
            stdl_log_msg(0, "\n");
        }
    }

    return STDL_ERR_OK;
}

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
            float c;

            if(rh->XmYamp != NULL)
                c = rh->XpYamp[iexci * ctx->ncsfs + kia] * rh->XmYamp[iexci * ctx->ncsfs + kia];
            else
                c = powf(rh->XpYamp[iexci * ctx->ncsfs + kia], 2);

            if(c >= thresh) {
                size_t i = ctx->csfs[kia] / nvirt, a = ctx->csfs[kia] % nvirt, hi = ctx->nocc - i - 1;

                if(printed > 0)
                    stdl_log_msg(1, "             ");

                stdl_log_msg(1, " %5.1f%% (% 6.4f) H", c * 100, rh->XpYamp[iexci * ctx->ncsfs + kia] * s2o2);

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

static inline size_t _biop(size_t op0, size_t op1) {
    return (op0 < op1) ? op0 * STDL_OP_COUNT + op1 : op1 * STDL_OP_COUNT + op0;
}

int stdl_log_property_g2e_moments(stdl_responses_handler *rh, stdl_context *ctx, stdl_operator ops[2], size_t nexci, float *tg2e) {
    assert(rh != NULL && ctx != NULL && tg2e != NULL && nexci > 0);

    size_t dims[] = {STDL_OPERATOR_DIM[ops[0]], STDL_OPERATOR_DIM[ops[1]]};
    int show[] = {1, ops[0] != ops[1]};
    size_t biop = _biop(ops[0], ops[1]);
    int show_biop = dims[0] == dims[1];

    // headers
    stdl_log_msg(0, "**    -------- Energy -------");
    for (size_t iop = 0; iop < 2; ++iop) {
        if(show[iop]) {
            stdl_log_msg(0, " ");
            for (size_t cpt = 0; cpt < dims[iop] * 2; ++cpt)
                stdl_log_msg(0, "-");
            stdl_log_msg(0, " <m|%4s|n> ", STDL_OPERATOR_NAME[ops[iop]]);
            for (size_t cpt = 0; cpt < dims[iop] * 7 - 13; ++cpt)
                stdl_log_msg(0, "-");
        }
    }

    if(show_biop) {
        stdl_log_msg(0, " --------");
        if(biop == _biop(STDL_OP_DIPL, STDL_OP_DIPL)
           || biop == _biop(STDL_OP_DIPV, STDL_OP_DIPV)
           || biop == _biop(STDL_OP_DIPL, STDL_OP_ANGM)
           || biop == _biop(STDL_OP_DIPV, STDL_OP_ANGM)
                )
            stdl_log_msg(0, "---------");
    }

    stdl_log_msg(0, "\n");

    stdl_log_msg(0, "       (Eh)     (eV)    (nm) ");
    for (size_t iop = 0; iop < 2; ++iop) {
        if(show[iop]) {
            for (size_t cpt = 0; cpt < dims[iop]; ++cpt)
                stdl_log_msg(0, "   %-3d   ", cpt);
        }
    }

    if(show_biop) {
        stdl_log_msg(0, "    A·B  ");
        if(biop == _biop(STDL_OP_DIPL, STDL_OP_DIPL))
            stdl_log_msg(0, "    fL   ");
        else if(biop == _biop(STDL_OP_DIPV, STDL_OP_DIPV))
            stdl_log_msg(0, "    fV   ");
        else if(biop == _biop(STDL_OP_DIPL, STDL_OP_ANGM))
            stdl_log_msg(0, "    RL   ");
        else if(biop == _biop(STDL_OP_DIPV, STDL_OP_ANGM))
            stdl_log_msg(0, "    RV   ");
    }
    stdl_log_msg(0, "\n");

    for (size_t iexci = 0; iexci < nexci; ++iexci) {
        // print energies
        stdl_log_msg(0, "%4ld %8.5f %7.3f %7.2f", iexci + 1, rh->eexci[iexci], rh->eexci[iexci] * STDL_CONST_AU_TO_EV, STDL_CONST_HC / rh->eexci[iexci]);

        size_t shift = 0;
        for (int iop = 0; iop < 2; ++iop) {
            if(show[iop]) {
                for (size_t cpt = 0; cpt < dims[iop]; ++cpt) {
                    stdl_log_msg(0, " % 8.5f", tg2e[(shift + cpt) * nexci + iexci]);
                }
            }

            shift += dims[iop];
        }

        if(show_biop) {
            float dotp = .0f;
            for (size_t cpt = 0; cpt < dims[1]; ++cpt) {
                dotp += tg2e[cpt * nexci + iexci] * tg2e[(dims[0] + cpt) * nexci + iexci];
            }

            stdl_log_msg(0, " % 8.5f", dotp);

            if(biop == _biop(STDL_OP_DIPL, STDL_OP_DIPL))
                stdl_log_msg(0, " % 8.5f", 2.f/3 * rh->eexci[iexci] * dotp);
            if(biop == _biop(STDL_OP_DIPV, STDL_OP_DIPV))
                stdl_log_msg(0, " % 8.5f", 2.f/3 / rh->eexci[iexci] * dotp);
            if(biop == _biop(STDL_OP_DIPL, STDL_OP_ANGM))
                stdl_log_msg(0, " % 8.5f", dotp);
            if(biop == _biop(STDL_OP_DIPV, STDL_OP_ANGM))
                stdl_log_msg(0, " % 8.5f", -1.f / rh->eexci[iexci] * dotp);

        }

        stdl_log_msg(0, "\n");
    }

    return STDL_ERR_OK;
}

int stdl_log_property_e2e_moments(stdl_responses_handler *rh, stdl_context *ctx, stdl_operator ops[3], size_t nexci, float *te2e) {
    assert(rh != NULL && ctx != NULL && te2e != NULL && nexci > 0);

    size_t dims[] = {STDL_OPERATOR_DIM[ops[0]], STDL_OPERATOR_DIM[ops[1]], STDL_OPERATOR_DIM[ops[2]]};
    int show[] =  {1, ops[1] != ops[0], ops[2] != ops[0] && ops[2] != ops[1]};

    // headers
    stdl_log_msg(0, "**  ");
    for (size_t iop = 0; iop < 3; ++iop) {
        if(show[iop]) {
            stdl_log_msg(0, " ");
            for (size_t cpt = 0; cpt < dims[iop] * 1; ++cpt)
                stdl_log_msg(0, "-");
            stdl_log_msg(0, " <m|bar(%4s)|m> ", STDL_OPERATOR_NAME[ops[iop]], STDL_OPERATOR_NAME[ops[iop]]);
            for (size_t cpt = 0; cpt < dims[iop] * 8 - 18; ++cpt)
                stdl_log_msg(0, "-");
        }
    }

    stdl_log_msg(0, "\n");
    stdl_log_msg(0, "    ");

    for (size_t iop = 0; iop < 3; ++iop) {
        if(show[iop]) {
            for (size_t cpt = 0; cpt < dims[iop]; ++cpt)
                stdl_log_msg(0, "   %-3d   ", cpt);
        }
    }
    stdl_log_msg(0, "\n");

    for(size_t iexci=0; iexci < nexci; iexci++) {
        stdl_log_msg(0, "%4ld", iexci + 1);

        size_t shift = 0;
        for (int iop = 0; iop < 3; ++iop) {
            if(show[iop]) {
                for (size_t cpt = 0; cpt < dims[0]; ++cpt) {
                    stdl_log_msg(0, " % 8.5f", te2e[(shift + cpt) * STDL_MATRIX_SP_SIZE(nexci) + STDL_MATRIX_SP_IDX(iexci, iexci)]);
                }
            }

            shift += dims[iop];
        }

        stdl_log_msg(0, "\n");
    }

    size_t biops[] = {_biop(ops[0], ops[1]), _biop(ops[1], ops[2]), _biop(ops[0], ops[2])};
    int show_biops[] = {dims[0] == dims[1], dims[0] == dims[2] && biops[1] != biops[0], dims[1] == dims[2] && biops[2] != biops[0] && biops[2] != biops[1]};

    // headers
    stdl_log_msg(0, "**         -------- Energy -------");
    for (size_t iop = 0; iop < 3; ++iop) {
        if(show[iop]) {
            stdl_log_msg(0, " ");
            for (size_t cpt = 0; cpt < dims[iop] * 2; ++cpt)
                stdl_log_msg(0, "-");
            stdl_log_msg(0, " <m|%4s|n> ", STDL_OPERATOR_NAME[ops[iop]]);
            for (size_t cpt = 0; cpt < dims[iop] * 7 - 13; ++cpt)
                stdl_log_msg(0, "-");
        }
    }

    for (size_t biop = 0; biop < 3; ++biop) {
        if(show_biops[biop]) {
            stdl_log_msg(0, " --------");
            if(biops[biop] == _biop(STDL_OP_DIPL, STDL_OP_DIPL)
                || biops[biop] == _biop(STDL_OP_DIPV, STDL_OP_DIPV)
                || biops[biop] == _biop(STDL_OP_DIPL, STDL_OP_ANGM)
                || biops[biop] == _biop(STDL_OP_DIPV, STDL_OP_ANGM)
            )
                stdl_log_msg(0, "---------");
        }
    }

    stdl_log_msg(0, "\n");

    stdl_log_msg(0, "   m    n   (Eh)     (eV)    (nm) ");
    for (size_t iop = 0; iop < 3; ++iop) {
        if(show[iop]) {
            for (size_t cpt = 0; cpt < dims[iop]; ++cpt)
                stdl_log_msg(0, "   %-3d   ", cpt);
        }
    }

    for (size_t ibiop = 0; ibiop < 3; ++ibiop) {
        if(show_biops[ibiop]) {
            stdl_log_msg(0, "    %c·%c  ", ibiop == 1 ? 'B' : 'A', ibiop == 0 ? 'B' : 'C');
            if(biops[ibiop] == _biop(STDL_OP_DIPL, STDL_OP_DIPL))
                stdl_log_msg(0, "    fL   ");
            else if(biops[ibiop] == _biop(STDL_OP_DIPV, STDL_OP_DIPV))
                stdl_log_msg(0, "    fV   ");
            else if(biops[ibiop] == _biop(STDL_OP_DIPL, STDL_OP_ANGM))
                stdl_log_msg(0, "    RL   ");
            else if(biops[ibiop] == _biop(STDL_OP_DIPV, STDL_OP_ANGM))
                stdl_log_msg(0, "    RV   ");
        }
    }

    stdl_log_msg(0, "\n");

    for(size_t iexci=0; iexci < nexci - 1; iexci++) {
        for(size_t jexci=iexci + 1; jexci < nexci; jexci++) {
            float ediff = rh->eexci[jexci] - rh->eexci[iexci];
            stdl_log_msg(0, "%4ld %4ld %8.5f %7.3f %7.2f",iexci + 1, jexci + 1, ediff, ediff * STDL_CONST_AU_TO_EV, STDL_CONST_HC / ediff);

            size_t shift = 0;
            for (int iop = 0; iop < 3; ++iop) {
                if(show[iop]) {
                    for (size_t cpt = 0; cpt < dims[iop]; ++cpt) {
                        stdl_log_msg(0, " % 8.5f", te2e[(shift + cpt) * STDL_MATRIX_SP_SIZE(nexci) + STDL_MATRIX_SP_IDX(iexci, jexci)]);
                    }
                }

                shift += dims[iop];
            }

            for (size_t biop = 0; biop < 3; ++biop) {
                if(show_biops[biop]) {
                    size_t xop = biop == 1 ? 1 : 0, yop = biop == 0 ? 1 : 2;
                    float dotp = 0;
                    for (size_t cpt = 0; cpt < dims[xop]; ++cpt) {
                        dotp += te2e[((xop == 1? dims[0]: 0) + cpt) * STDL_MATRIX_SP_SIZE(nexci) + STDL_MATRIX_SP_IDX(iexci, jexci)] *
                                te2e[(dims[0] + (yop == 2? dims[1] : 0) + cpt) * STDL_MATRIX_SP_SIZE(nexci) + STDL_MATRIX_SP_IDX(iexci, jexci)];
                    }

                    stdl_log_msg(0, " % 8.5f", dotp);
                    if(biops[biop] == _biop(STDL_OP_DIPL, STDL_OP_DIPL))
                        stdl_log_msg(0, " % 8.5f", 2.f/3 * ediff * dotp);
                    else if(biops[biop] == _biop(STDL_OP_DIPV, STDL_OP_DIPV))
                        stdl_log_msg(0, " % 8.5f", 2.f/3 / ediff * dotp);
                    else if(biops[biop] == _biop(STDL_OP_DIPL, STDL_OP_ANGM))
                        stdl_log_msg(0, " % 8.5f", dotp);
                    else if(biops[biop] == _biop(STDL_OP_DIPV, STDL_OP_ANGM))
                        stdl_log_msg(0, " % 8.5f", -1.f / ediff * dotp);
                }
            }

            stdl_log_msg(0, "\n");
        }
    }

    return STDL_ERR_OK;

}

