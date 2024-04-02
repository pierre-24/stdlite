#ifndef STDLITE_LINEAR_ALGEBRA_H
#define STDLITE_LINEAR_ALGEBRA_H

#ifdef USE_MKL
#include <mkl.h>
#define STDL_LA_INT MKL_INT
#else
#include <cblas.h>
#include <lapacke.h>
#define STDL_LA_INT int
#endif

#endif //STDLITE_LINEAR_ALGEBRA_H
