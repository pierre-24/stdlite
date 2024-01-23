#include <stdio.h>

#include "stdlite/matrix.h"
#include "stdlite/errors.h"


int stdl_matrix_numbers_print(size_t rows, size_t columns, double *matrix, int is_symmetric) {
    size_t i = 0, j, c;
    while (i < columns) {
        // header
        j = 0;
        printf("    ");
        while (j < STDL_MATRIX_MAX_COLS && (i + j) < columns){
            c = i + j;
            printf("    %6ld    ", c);
            j++;
        }
        printf("\n");

        // content
        for(size_t r = (is_symmetric? i:0); r < rows; r++) {
            printf("%4ld", r);
            j = 0;
            while (j < STDL_MATRIX_MAX_COLS && i + j < columns && ((is_symmetric && i + j <= r) || !is_symmetric)){
                c = i + j;
                printf(" % .6e", matrix[r * columns + c]);
                j++;
            }
            printf("\n");
        }


        i += STDL_MATRIX_MAX_COLS;
    }

    return STDL_ERR_OK;
}
