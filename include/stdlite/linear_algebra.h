#ifndef STDLITE_LINEAR_ALGEBRA_H
#define STDLITE_LINEAR_ALGEBRA_H

#include <ctype.h>

#ifdef USE_MKL
#include <mkl.h>
#define STDL_LA_INT MKL_INT
#else
#ifdef USE_OPENBLAS_INC
#include <openblas/cblas.h>
#else
#include <cblas.h>
#endif
#include <lapacke.h>
#ifndef STDL_LA_INT
#define STDL_LA_INT lapack_int
#endif
#endif

#endif //STDLITE_LINEAR_ALGEBRA_H
