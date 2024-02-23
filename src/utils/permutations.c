#include <string.h>
#include <assert.h>
#include <stdio.h>

#include "stdlite/utils/permutations.h"
#include "stdlite/logging.h"
#include "stdlite/helpers.h"

// create a new element of the permutations, and copy `set` (in its current state) in it.
int _element_new(void* set, size_t szx, stdl_permutations** element) {
    *element = malloc(sizeof(stdl_permutations));
    STDL_ERROR_HANDLE_AND_REPORT(*element == NULL, return STDL_ERR_MALLOC, "malloc");

    (*element)->next = NULL;

    (*element)->perm = malloc(szx);
    STDL_ERROR_HANDLE_AND_REPORT((*element)->perm == NULL, stdl_permutations_delete(*element); return STDL_ERR_MALLOC, "malloc");

    memcpy((*element)->perm, set, szx);

    return STDL_ERR_OK;
}

// swap two elements of `base`.
void _memswap(void* base, size_t i, size_t j, size_t elemsz) {
    assert(base != NULL && elemsz > 0);
    if (i == j)
        return;

    char* byte1 = (char*) base + i * elemsz, *byte2 = (char *) base + j * elemsz;

    while (elemsz--) {

        char byte = *byte1;
        *byte1++ = *byte2;
        *byte2++ = byte;
    }
}

// Populate the set of permutations.
// This is based on the heap's algorithm (https://en.wikipedia.org/wiki/Heap%27s_algorithm)
int _populate(stdl_permutations** last, void* set, size_t nelem, size_t elemsz, size_t k) {

    assert(set != NULL && nelem > 0 && elemsz > 0 && last != NULL);

    int error;

    if(k == 1) {
        stdl_permutations* perm = NULL;
        error = _element_new(set, nelem * elemsz, &perm);
        STDL_ERROR_HANDLE(error, return error);

        (*last)->next = perm;
        *last = perm;

        return STDL_ERR_OK;
    } else {
        error = _populate(last, set, nelem, elemsz, k - 1);
        STDL_ERROR_HANDLE(error, return error);

        for (size_t i = 0; i < k - 1; ++i) {

            if(k % 2 == 0)
                _memswap(set, i, k - 1, elemsz);
            else
                _memswap(set, 0, k - 1, elemsz);

            error = _populate(last, set, nelem, elemsz, k - 1);
            STDL_ERROR_HANDLE(error, return error);
        }
    }

    return STDL_ERR_OK;
}

int stdl_permutations_new(void* original_set, size_t nelm, size_t elmsz, stdl_permutations** permutations) {

    assert(original_set != NULL && nelm > 0 && elmsz > 0 && permutations != NULL);

    stdl_permutations* dummy = malloc(sizeof(stdl_permutations));
    STDL_ERROR_HANDLE_AND_REPORT(dummy == NULL, return STDL_ERR_MALLOC, "malloc");
    dummy->next = NULL;

    stdl_permutations* first = dummy;
    stdl_permutations** last = &dummy;

    int error = _populate(last, original_set, nelm, elmsz, nelm);
    STDL_ERROR_HANDLE(error, free(dummy); return error);

    *permutations = first->next;
    free(first);

    STDL_DEBUG("create permutation %p", *permutations);

    return STDL_ERR_OK;
}

int stdl_permutations_delete(stdl_permutations* permutations) {
    assert(permutations != NULL);

    STDL_DEBUG("delete permutation element %p", permutations);

    if(permutations->next != NULL)
        stdl_permutations_delete(permutations->next);

    STDL_FREE_ALL(permutations->perm, permutations);
    return STDL_ERR_OK;
}

int stdl_permutations_remove_duplicates(stdl_permutations* permutations, size_t nelm, size_t elmsz) {
    size_t szx = nelm * elmsz;

    stdl_permutations* current = permutations;
    while (current != NULL) {
        stdl_permutations *next = current->next, *previous = current;
        while (next != NULL) {
            if(memcmp(current->perm, next->perm, szx) == 0) {
                previous->next = next->next;
                next->next = NULL;
                stdl_permutations_delete(next);
                next = previous->next;
            } else {
                previous = next;
                next = next->next;
            }

        }

        current = current->next;
    }

    return STDL_ERR_OK;
}
