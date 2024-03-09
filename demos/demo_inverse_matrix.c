#include <stdio.h>

#define NUMERICAL_IMPLEMENTATION
#include "../numerical.h"

int main()
{
    #define n_dim 3
    float coefs[n_dim][n_dim] = {
        {1.0f, -3.0f, 1.0f},
        {3.0f, -4.0f, 1.0f},
        {1.0f, 1.0f, -1.0f}
    };
    Matrix A = matrix_from_literal(n_dim, n_dim, coefs);
    Matrix A_inv = inverse_matrix(A);
    print_matrix(A, "A");
    print_matrix(A_inv, "A^-1");
}