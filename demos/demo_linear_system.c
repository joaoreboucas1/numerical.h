#include <stdio.h>
#include <math.h>

#define NUMERICAL_IMPLEMENTATION
#include "../numerical.h"


int main()
{
    // Equivalent to the system:
    // x - 3y + z = 2
    // 3x - 4y + z = 0
    // x + y - z = 1
    // Solution is x = -3.5, y = -5, z = -9.5
    #define n_dim 3
    float coefs[n_dim][n_dim] = {
        {1.0f, -3.0f, 1.0f},
        {3.0f, -4.0f, 1.0f},
        {1.0f, 1.0f, -1.0f}
    };
    float bs[n_dim] = {2.0f, 0.0f, 1.0f};
    
    printf("System of equations:\n");
    for (size_t i = 0; i < n_dim; i++) {
        for (size_t j = 0; j < n_dim; j++) {
            printf("%.1fx_%zu", fabsf(coefs[i][j]), j);
            if (j < n_dim - 1) {
                if (coefs[i][j+1] > 0) printf(" + ");
                else printf(" - ");
            }
        }
        printf(" = %.2f\n", bs[i]);
    }

    Matrix A = matrix_from_literal(n_dim, n_dim, coefs);
    Matrix B = column_matrix(n_dim, bs);
    Matrix solution = solve_linear_system(A, B, n_dim);
    printf("Solution: ");
    for (size_t i = 0; i < n_dim; i++) {
        printf("x_%zu = %.2f ", i, matrix_at(solution, i, 0));
    }
    printf("\n");
    return 0;
}