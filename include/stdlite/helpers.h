#ifndef STDLITE_HELPERS_H
#define STDLITE_HELPERS_H

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


#endif //STDLITE_HELPERS_H
