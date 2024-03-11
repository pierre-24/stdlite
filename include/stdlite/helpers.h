#ifndef STDLITE_HELPERS_H
#define STDLITE_HELPERS_H

#include <math.h>

/**
 * Apply a function on a list of objects.
 * From B. Klemens in *21st century C* (O'Reilly).
 * @ingroup helpers
 */
#define STDL_APPLY(type, fn, ...) {                         \
    void* _stper = (int[]) {0};                             \
    type** _lst_to_apply = (type*[]) {__VA_ARGS__, _stper}; \
    for(int i=0; _lst_to_apply[i] != _stper; i++) {         \
        fn(_lst_to_apply[i]);                               \
    }                                                       \
}

/**
 * Free `a` if not `NULL`
 * @ingroup helpers
 */
#define STDL_FREE_IF_USED(a) if((a) != NULL) free(a)

/**
 * Free all objects if they are non-`NULL`.
 * @ingroup helpers
 */
#define STDL_FREE_ALL(...) STDL_APPLY(void, STDL_FREE_IF_USED, __VA_ARGS__)

/// Conversion to eV
/// @ingroup helpers
#ifndef STDL_CONST_AU_TO_EV
#define STDL_CONST_AU_TO_EV 27.21139
#endif

/// Conversion to nanometers
/// @ingroup helpers
#ifndef STDL_CONST_HC
#define STDL_CONST_HC 45.56335
#endif

/// Conversion to Angstrom
/// @ingroup helpers
#ifndef STDL_CONST_AU_TO_ANG
#define STDL_CONST_AU_TO_ANG  0.529177210903
#endif


#define MX(a, b) ((a < b)? b: a)

static inline int stdl_float_equals(float a, float b, float epsilon) {
    return (fabsf(a - b) <= epsilon * MX(fabsf(a), fabsf(b)));
}

#endif //STDLITE_HELPERS_H
