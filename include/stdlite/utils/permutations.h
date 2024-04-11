#ifndef STDLITE_PERMUTATIONS_H
#define STDLITE_PERMUTATIONS_H

#include <stdlib.h>

/**
 * A sequence of permutations.
 * Implemented as a linked list.
 * @ingroup permutations
 */
struct stdl_permutation_ {
    /// `uint8_t[nelements*sz]`, one of the possible permutation of the original set.
    void* perm;

    /// The next element in the set, might be `NULL`.
    struct stdl_permutation_* next;
};

typedef struct stdl_permutation_ stdl_permutations;

/**
 * Build a sequence of `fac(nelm)` permutations from `original_set` (with duplicates if elements of `original_set` are equals).
 *
 * @param original_set `uint8_t[nelm*elemz]` a set of `nelm` elements of `elmsz` bytes each to be permuted.
 * @param nelm number of element in the original_set
 * @param elmsz size of each element of the original_set **in byte** (use `sizeof()`)
 * @param[out] permutations the resulting sequence of permutations.
 * @return error code
 * @ingroup permutations
 */
int stdl_permutations_new(void* original_set, size_t nelm, size_t elmsz, stdl_permutations** permutations);

/**
 * Delete a sequence of permutations.
 *
 * @param permutations a valid sequence of permutations
 * @return error code
 * @ingroup permutations
 */
int stdl_permutations_delete(stdl_permutations* permutations);

/**
 * Remove duplicates in a sequence of permutations, by using `memcmp` in `stdl_permutations_remove_duplicates_with_callback()` (i.e., performing a data comparison).
 *
 * @param permutations a valid sequence of permutations
 * @param nelm number of element in the original set
 * @param elmsz size of each element of the original set **in byte** (use `sizeof()`)
 * @return error code
 * @ingroup permutations
 */
int stdl_permutations_remove_duplicates(stdl_permutations* permutations, size_t nelm, size_t elmsz);

/**
 * Remove duplicates in a sequence of permutations, using `cmp` to compare two permutations
 * Note that this implementation scales as $\mathcal{O}(N^2)$ in the worst case, with $N$ the length of the sequence.
 *
 * @param permutations a valid sequence of permutations
 * @param cmp comparator of the form `int cmp(void* perm1, void* perm2, void* data)`, must return 0 if `perm1` and `perm2` are different, 1 if they are equals.
 * @param data data passed to each call to the callback
 * @return error code
 * @ingroup permutations
 */
int stdl_permutations_remove_duplicates_with_callback(stdl_permutations* permutations, int (*cmp)(void*, void*, void*), void* data);

#endif //STDLITE_PERMUTATIONS_H
