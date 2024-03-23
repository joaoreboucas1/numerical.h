#include <stdio.h>
#include <assert.h>
#include <stdbool.h>

#define NUMERICAL_IMPLEMENTATION
#include "../numerical.h"

ALLOCATES Matrix* QR_decomposition(Matrix M)
{
    Matrix R = zero_matrix(M.rows, M.cols);

    Matrix u = zero_matrix(M.rows, M.cols);
    Matrix e = zero_matrix(M.rows, M.cols);

    float norm_u0 = 0.0f;
    for (size_t i = 0; i < M.rows; i++) {
        norm_u0 += matrix_at(M, i, 0)*matrix_at(M, i, 0);
    }
    norm_u0 = sqrtf(norm_u0);
    for (size_t i = 0; i < M.rows; i++) {
        matrix_at(u, i, 0) = matrix_at(M, i, 0);
        matrix_at(e, i, 0) = matrix_at(M, i, 0)/norm_u0;
    }

    for (size_t j = 1; j < M.cols; j++) {
        float norm_uj = 0.0f;
        for (size_t i = 0; i < M.cols; i++) matrix_at(u, i, j) = matrix_at(M, i, j);
        
        for (size_t k = 0; k < j; k++) {
            float dot_aj_uk = 0.0f;
            float dot_uk_uk = 0.0f;
            for (size_t n = 0; n < M.cols; n++) dot_aj_uk += matrix_at(M, n, j)*matrix_at(u, n, k);
            for (size_t n = 0; n < M.cols; n++) dot_uk_uk += matrix_at(u, n, k)*matrix_at(u, n, k);
            for (size_t i = 0; i < M.cols; i++) matrix_at(u, i, j) -= dot_aj_uk/dot_uk_uk*matrix_at(u, i, k);
        }
        for (size_t i = 0; i < M.cols; i++) norm_uj += matrix_at(u, i, j)*matrix_at(u, i, j);
        norm_uj = sqrtf(norm_uj);
        for (size_t i = 0; i < M.rows; i++) {
            matrix_at(e, i, j) = matrix_at(u, i, j)/norm_uj;
        }
    }

    for (size_t i = 0; i < M.rows; i++) {
        for (size_t j = i; j < M.cols; j++) {
            float dot_ei_aj = 0.0f;
            for (size_t k = 0; k < M.cols; k++) dot_ei_aj += matrix_at(M, k, j)*matrix_at(e, k, i);
            matrix_at(R, i, j) = dot_ei_aj;
        }
    }

    free_matrix(u);
    return (Matrix[2]) { e, R };
}

#define TRIANGULAR_TOL 1e-10
bool is_lower_triangular(Matrix M)
{
    if (M.rows != M.cols) return false;
    for (size_t i = 1; i < M.rows; i++) {
        for (size_t j = 0; j < i; j++) {
            if (fabsf(matrix_at(M, i, j)) > TRIANGULAR_TOL) return false;
        }
    }
    return true;
}

bool is_upper_triangular(Matrix M)
{
    if (M.rows != M.cols) return false;
    for (size_t i = 0; i < M.rows; i++) {
        for (size_t j = i+1; j < M.cols; j++) {
            if (fabsf(matrix_at(M, i, j)) > TRIANGULAR_TOL) return false;
        }
    }
    return true;
}

bool is_triangular(Matrix M)
{
    return is_upper_triangular(M) || is_lower_triangular(M);
}

ALLOCATES float* find_eigenvalues(Matrix M)
{
    bool solved;
    Matrix *QR, *M_current;
    Matrix Q, R, M_next;

    if (M.rows != M.cols) {
        fprintf(stderr, "ERROR: finding eigenvalue of non-square matrix with %zu rows and %zu cols\n", M.rows, M.cols);
        return NULL;
    }
    
    float *eigenvalues = malloc(M.cols*sizeof(float));
    if (eigenvalues == NULL) {
        fprintf(stderr, "ERROR: could not allocate memory for %zu eigenvalues\n", M.rows);
        return NULL;
    }

    size_t iters = 0;
    const size_t max_iters = 10;
    M_current = &M;
    do {
        QR = QR_decomposition(*M_current);
        Q = QR[0];
        R = QR[1];
        M_next = matrix_multiplication(R, Q);
        print_matrix(M_next, "M");
        if (is_triangular(M_next)) {
            solved = true;
            break;
        } else {
            if (iters > 0) free_matrix(*M_current);
            free_matrix(R);
            free_matrix(Q);
            M_current = &M_next;
        }
        iters += 1;
    } while (iters < max_iters);

    if (solved) {
        for (size_t i = 0; i < M.cols; i++) {
            eigenvalues[i] = matrix_at((M_next), i, i);
        }
        free_matrix(Q);
        free_matrix(R);
    }
    free_matrix(M_next);

    return eigenvalues;
    
}

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
    free_matrix(M);
    return 0;
}