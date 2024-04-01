#include <stdio.h>
#include <assert.h>
#include <stdbool.h>

#define NUMERICAL_IMPLEMENTATION
#include "../numerical.h"

int main()
{
    #define D 3
    float entries[D][D] = {
        {-2.0f, -4.0f, 2.0f},
        {-2.0f, 1.0f, 2.0f},
        {4.0f, 2.0f, 5.0f}
    };
    Matrix M = matrix_from_literal(D, D, entries);
    print_matrix(M, "M");
    
    float* eigenvalues = find_eigenvalues(M);
    printf("Eigenvalues: ");
    for (size_t i = 0; i < D; i++) {
        printf("%.3f ", eigenvalues[i]);
    }
    printf("\n");
    free_matrix(M);
    return 0;
}