#ifndef STDLITE_LINEAR_ALGEBRA_H
#define STDLITE_LINEAR_ALGEBRA_H

#ifdef USE_MKL
#include <mkl.h>
#define STDL_LA_INT MKL_INT
#else
#include <cblas.h>
#include <lapacke.h>
#ifndef STDL_LA_INT
#ifdef USE_OPENBLAS64
#define STDL_LA_INT long
#else
#define STDL_LA_INT int
#endif
#endif
#endif

#endif //STDLITE_LINEAR_ALGEBRA_H
