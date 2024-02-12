#include <stdlite/utils/permutations.h>

#include "tests_suite.h"

void test_permutations1() {
    stdl_permutations* set = NULL;
    ASSERT_STDL_OK(stdl_permutations_new((int[]) {0}, 1, sizeof(int), &set));

    TEST_ASSERT_NOT_NULL(set);

    size_t n = 0;
    stdl_permutations* current = set;

    while(current != NULL) {
        current = current->next;
        n++;
    }

    TEST_ASSERT_EQUAL_INT(1, n);

    ASSERT_STDL_OK(stdl_permutations_delete(set));
}

void test_permutations1_remove_duplicates() {
    stdl_permutations* set = NULL;

    ASSERT_STDL_OK(stdl_permutations_new((int[]) {0}, 1, sizeof(int), &set));
    TEST_ASSERT_NOT_NULL(set);

    stdl_permutations_remove_duplicates(set, 1, sizeof(int));

    size_t n = 0;
    stdl_permutations* current = set;

    while(current != NULL) {
        current = current->next;
        n++;
    }

    TEST_ASSERT_EQUAL_INT(1, n);

    ASSERT_STDL_OK(stdl_permutations_delete(set));
}

void test_permutations4() {
    stdl_permutations* set = NULL;
    ASSERT_STDL_OK(stdl_permutations_new((int[]) {0, 1, 2, 3}, 4, sizeof(int), &set));

    TEST_ASSERT_NOT_NULL(set);

    size_t n = 0;
    stdl_permutations* current = set;

    int met[24] = {0};

    while(current != NULL) {
        int* value = (int*) current->perm;
        int cvalue = value[0] * 64 + value[1] * 16 + value[2] * 4 + value[3];
        for (size_t j = 0; j < n; ++j) { // check uniqueness
            TEST_ASSERT_NOT_EQUAL_INT(met[j], cvalue);
        }

        met[n] = cvalue;

        current = current->next;
        n++;
    }

    TEST_ASSERT_EQUAL_INT(24, n);

    ASSERT_STDL_OK(stdl_permutations_delete(set));
}

void test_permutations4_remove_duplicates() {
    stdl_permutations* set = NULL;

    ASSERT_STDL_OK(stdl_permutations_new((int[]) {0, 1, 2, 2}, 4, sizeof(int), &set));
    TEST_ASSERT_NOT_NULL(set);

    size_t n = 0;
    stdl_permutations* current = set;

    // store the 24 values
    int vals[24] = {0};

    while(current != NULL) {
        int* value = (int*) current->perm;
        int cvalue = value[0] * 64 + value[1] * 16 + value[2] * 4 + value[3];
        vals[n] = cvalue;

        current = current->next;
        n++;
    }

    TEST_ASSERT_EQUAL_INT(24, n);

    // check that each value has a duplicate
    for (int i = 0; i < 24; ++i) {
        int val = vals[i];
        int found = -1;
        for (int j = 0; j < 24; ++j) {
            if(i == j)
                continue;

            if(vals[i] == val) {
                found = j;
                break;
            }

            TEST_ASSERT_NOT_EQUAL_INT(found, -1);
        }
    }

    // now remove duplicates:
    stdl_permutations_remove_duplicates(set, 4, sizeof(int));

    n = 0;
    current = set;

    int met[12] = {0};

    while(current != NULL) {
        int* value = (int*) current->perm;
        int cvalue = value[0] * 64 + value[1] * 16 + value[2] * 4 + value[3];
        for (size_t j = 0; j < n; ++j) { // check uniqueness
            TEST_ASSERT_NOT_EQUAL_INT(met[j], cvalue);
        }

        met[n] = cvalue;

        current = current->next;
        n++;
    }

    TEST_ASSERT_EQUAL_INT(12, n);

    ASSERT_STDL_OK(stdl_permutations_delete(set));
}

void test_permutations3_operator() {
    // Quadratic response: permutation over {(coo[0], -f[0]-f[1]), (coo[1], f[0]), (coo[2], f[1])},
    // where coo are coordinates and f are frequencies.

    // if coo = {'x', 'y', 'z'} and f = {1, 1} → 6 possible permutations
    stdl_permutations* set = NULL;
    ASSERT_STDL_OK(stdl_permutations_new((int[]) {'x', -2, 'y', 1, 'z', 1}, 3, 2 * sizeof(int), &set));

    stdl_permutations_remove_duplicates(set, 3, 2 * sizeof(int));

    size_t n = 0;
    stdl_permutations* current = set;

    while(current != NULL) {
        current = current->next;
        n++;
    }

    TEST_ASSERT_EQUAL_INT(6, n);
    ASSERT_STDL_OK(stdl_permutations_delete(set));

    // if coo = {'x', 'y', 'y'} and f = {1, 1} → 3 possible permutations
    ASSERT_STDL_OK(stdl_permutations_new((int[]) {'x', -2, 'y', 1, 'y', 1}, 3, 2 * sizeof(int), &set));

    stdl_permutations_remove_duplicates(set, 3, 2 * sizeof(int));

    n = 0;
    current = set;

    while(current != NULL) {
        current = current->next;
        n++;
    }

    TEST_ASSERT_EQUAL_INT(3, n);

    ASSERT_STDL_OK(stdl_permutations_delete(set));
}
