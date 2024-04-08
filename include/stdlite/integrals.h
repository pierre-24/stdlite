#ifndef STDLITE_INTEGRALS_H
#define STDLITE_INTEGRALS_H

#include <assert.h>
#include <stdlite/logging.h>
#include <stdlite/basis.h>

#include <cint_funcs.h>

/**
 * Operators for the linear responses
 * @ingroup ops_integrals
 */
enum stdl_operator_ {
    /// Dipole length
    STDL_OP_DIPL,

    /// Dipole velocity
    STDL_OP_DIPV,

    /// Angular momentum
    STDL_OP_ANGM,

    /// Overlap
    STDL_OP_OVLP,

    STDL_OP_COUNT
};

typedef enum stdl_operator_ stdl_operator;

/**
 * Name of each operator
 * @ingroup ops_integrals
 */
static char* STDL_OPERATOR_NAME[STDL_OP_COUNT] = {
        "dipl",
        "dipv",
        "angm",
        "ovlp",
};

/**
 * Dimensionality of of each operator
 * @ingroup ops_integrals
 */
static size_t STDL_OPERATOR_DIM[STDL_OP_COUNT] = {
        3, // dipl
        3, // dipv
        3, // angm
        1, // ovlp
};

/**
 * Hermicity of each operator
 * @ingroup ops_integrals
 */
static int STDL_OPERATOR_HERMITIAN[STDL_OP_COUNT] = {
        1, // dipl
        0, // dipv
        0, // angm
        1, // ovlp
};

typedef FINT (*_int1e_f)(
        double* /* buff */,
        FINT*,  /* NULL */
        FINT*,  /* shells */
        FINT*, FINT,  /* atm + natm */
        FINT*, FINT,  /* bas + nbas */
        double*,  /* env */
        CINTOpt*, /* optimizer */
        double* /* NULL */
        );

static _int1e_f STDL_OPERATOR_TO_CINT_CART[STDL_OP_COUNT] = {
        int1e_r_cart, // dipl
        int1e_p_cart, // dipv
        int1e_cg_irxp_sph, // angm
        int1e_ovlp_cart, // ovlp
};

static _int1e_f STDL_OPERATOR_TO_CINT_SPH[STDL_OP_COUNT] = {
        int1e_r_sph, // dipl
        int1e_p_sph, // dipv
        int1e_cg_irxp_sph, // angm
        int1e_ovlp_sph, // ovlp
};

/**
 * Compute one-electron ops_integrals values over AO, $\braket{\mu|op|\nu}$.
 *
 * @warning Depending on the operator, the result might be a skew-symmetric matrix
 * @param bs a valid basis set
 * @param op the operator
 * @param fac factor by witch all elements are multiplied (use `1.0` to keep the result unchanged)
 * @param[out] values `double[dim(op), STDL_MATRIX_SP_SIZE(nao)]` the resulting matrix
 * @return error code.
 * @ingroup ops_integrals
 */
int stdl_operator_int1e_dsp(stdl_basis *bs, stdl_operator op, double fac, double *values);

#endif //STDLITE_INTEGRALS_H
