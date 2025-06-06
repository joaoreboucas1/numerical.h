#ifndef NUMERICAL_H_
#define NUMERICAL_H_

#include <stddef.h>
#include <stdbool.h>
#include <complex.h>

// Matrix with linear memory layout
typedef struct {
    size_t rows;
    size_t cols;
    float *elements;
} Matrix;
#define matrix_at(M, i, j) (M).elements[(i)*(M).cols + (j)]

// Parameters to be passed into functions for numerical integration
typedef struct {
    int *int_args;
    size_t int_argc;
    float *float_args;
    size_t float_argc;
} Params;

// Linear interpolator
typedef struct {
    float *as;
    float *bs;
    float *xs;
    size_t n_points;
} LinearInterp;

// Cubic Spline following the implementation from https://en.wikipedia.org/wiki/Spline_(mathematics)#Algorithm_for_computing_natural_cubic_splines
typedef struct {
    float *as;
    float *bs;
    float *cs;
    float *ds;
    float *xs;
    size_t n_points;
} CubicSpline;

// Linear algebra
Matrix matrix_from_literal(size_t n_rows, size_t n_cols, float elements[n_rows][n_cols]); // Create a `Matrix` object from matrix literal float[rows][cols]
Matrix column_matrix(size_t n_rows, float *elements); // Create a `Matrix` object with a single column from 1D array
Matrix zero_matrix(size_t n_rows, size_t n_cols); // Create a `Matrix` object with all elements zero
Matrix identity(size_t n_rows); // a new `Matrix` object equivalent to the n_rows x n_rows identity 
void matrix_multiplication(Matrix A, Matrix B, Matrix *AB); // Performs matrix multiplication and returns a new `Matrix` M = AB
bool is_upper_triangular(Matrix A); // Checks if matrix M is upper triangular
bool is_lower_triangular(Matrix A); // Checks if matrix M is lower triangular
bool is_triangular(Matrix A); // Checks if matrix M is triangular
float *find_eigenvalues(Matrix A); // Finds the eigenvalues of the matrix A
void QR_decomposition(Matrix M, Matrix *Q, Matrix *R); // Performs QR decomposition of the matrix M, saving the results into Q and R
void LU_decomposition(Matrix A, Matrix *L, Matrix *U); // Performs LU decomposition of a `Matrix` A, returning two `Matrix` objects L and U
void inverse_matrix(Matrix A, Matrix *A_inv); // Returns a new `Matrix` which is the inverse of the input matrix A
Matrix solve_linear_system(Matrix A, Matrix B, size_t n_dim); // Solves the linear system A*X = B, returning X as a `Matrix` object with a single column
float determinant(Matrix A); // Computes the determinant of a `Matrix` A
void free_matrix(Matrix M); // Frees the dynamic storage of M

// Numerical integration
float trapezoid(float f(float), float a, float b, size_t N); // Integrate f from a to b with N points
float trapezoid_args(float f(float, Params), float a, float b, size_t N, Params params); // Integrate function with parameters

// ODE integration
float* odesolve(float f(float, float), float y_0, float x_0, float x_final, size_t nsteps); // Solve 1D ODE y'(y, x) = f(y, x) with initial conditions (x_0, y_0) and nsteps steps via Euler's algorithm
float **odesolve_nd(void y_prime(float*, float, float*), float* y_initial, float x_initial, float x_final, size_t n_steps, size_t n_dim); // Solve N-dimensional ODE y'(y, x) = f(y, x) (where y is `n_dim`-sized array) with initial conditions (x_0, y_0) and n_steps steps via Euler's algorithm. The user must allocate memory for the derivative array and pass the pointer into the function `y_prime` (so the user must manage the memory for the derivatives)

// Interpolation
LinearInterp get_linear_interpolator(float* xs, float* ys, size_t n_points); // Create a `LinearInterpolator` object from the data, representing the function y(x)
float eval_linear_interpolator(LinearInterp interp, float x); // Evaluate a `LinearInterpolator` object at `x`
void free_linear_interpolator(LinearInterp interp); // Frees the storage of the `LinearInterpolator` (`xs` is not freed)
CubicSpline get_cubic_spline(float* xs, float* ys, size_t n_points); // Create a `CubicSpline` object from the data, representing the function y(x)
float eval_cubic_spline(CubicSpline interp, float x); // Evaluate a `CubicSpline` object at `x`
void free_cubic_spline(CubicSpline interp); // Frees the storage of the `CubicSpline` (`xs` is not freed)

// Root-finding
float find_root(float f(float), float x_min, float x_max, size_t max_iters, float tol); // Finds a solution for the equation f(x) = 0 between x_min and x_max by secant rule. Iterates until abs(f(x)) < tol and by a maximum number of iterations

// Minimization
float minimize(float f(float), float x, float rate, size_t max_iters); // Finds a local minimum of the function f(x) using gradient descent and finite difference. Control descent rate and max number of iterations.

// FFT
void fft(float complex *in, float complex *out, int n); // Returns the fft of the input signal `in` (as float complex)

#ifdef NUMERICAL_IMPLEMENTATION
#include <math.h>
#include <stdlib.h>
#include <string.h>

Matrix matrix_from_literal(size_t n_rows, size_t n_cols, float elements[n_rows][n_cols])
{
    float *m_elements = malloc(n_rows*n_cols*sizeof(float));
    if (m_elements == NULL) {
        fprintf(stderr, "ERROR: could not allocate Matrix");
        return (Matrix) {0};
    }

    Matrix m = {.rows = n_rows, .cols = n_cols, .elements = m_elements};
    for (size_t i = 0; i < n_rows; i++) {
        for (size_t j = 0; j < n_cols; j++) {
            matrix_at(m, i, j) = elements[i][j];
        }
    }
    return m;
}

Matrix column_matrix(size_t n_rows, float *elements)
{
    float *m_elements = malloc(n_rows*sizeof(float));
    if (m_elements == NULL) {
        fprintf(stderr, "ERROR: could not allocate Matrix");
        return (Matrix) {0};
    }
    for (size_t i = 0; i < n_rows; i++) {
        m_elements[i] = elements[i];
    }
    return (Matrix) {.rows = n_rows, .cols = 1, .elements = m_elements};
}

Matrix zero_matrix(size_t n_rows, size_t n_cols)
{
    float *elements = calloc(n_rows*n_cols, sizeof(float));
    if (elements == NULL) {
        fprintf(stderr, "ERROR: could not allocate elements for zero matrix of size %zu x %zu\n", n_rows, n_cols);
        return (Matrix) {0};
    }
    return (Matrix) { .rows = n_rows, .cols = n_cols, .elements = elements };
}

Matrix identity(size_t n_rows)
{
    Matrix M = zero_matrix(n_rows, n_rows);
    for (size_t i = 0; i < n_rows; i++) {
        matrix_at(M, i, i) = 1.0f;
    }
    return M;
}

void matrix_multiplication(Matrix A, Matrix B, Matrix *AB) 
{
    if (A.cols != B.rows) {
        fprintf(stderr, "ERROR: cannot perform matrix multiplication because the number of rows of the first matrix (%zu) is different than the number of columns of the second matrix (%zu)\n", A.cols, B.rows);
        return;
    }

    if (AB->rows != A.rows || AB->cols != B.cols) {
        fprintf(stderr, "ERROR: cannot perform matrix multiplication because the number of rows of the first matrix (%zu) is different than the number of columns of the second matrix (%zu)\n", A.cols, B.rows);
        return;
    }

    for (size_t i = 0; i < AB->rows*AB->cols; i++) AB->elements[i] = 0.0f;

    for (size_t i = 0; i < AB->rows; i++) {
        for (size_t j = 0; j < AB->cols; j++) {
            for (size_t k = 0; k < A.rows; k++) {
                matrix_at(*AB, i, j) += matrix_at(A, i, k)*matrix_at(B, k, j);
            }
        }
    }
}

void QR_decomposition(Matrix M, Matrix *Q, Matrix *R)
{
    Matrix u = zero_matrix(M.rows, M.cols);

    for (size_t i = 0; i < Q->rows*Q->cols; i++) {
        Q->elements[i] = 0.0f;
        R->elements[i] = 0.0f;
    }

    float norm_u0 = 0.0f;
    for (size_t i = 0; i < M.rows; i++) {
        norm_u0 += matrix_at(M, i, 0)*matrix_at(M, i, 0);
    }
    norm_u0 = sqrtf(norm_u0);

    for (size_t i = 0; i < M.rows; i++) {
        matrix_at(u, i, 0) = matrix_at(M, i, 0);
        matrix_at(*Q, i, 0) = matrix_at(M, i, 0)/norm_u0;
    }

    for (size_t k = 1; k < M.cols; k++) {
        for (size_t i = 0; i < M.cols; i++) matrix_at(u, i, k) = matrix_at(M, i, k);
        
        for (size_t j = 0; j < k; j++) {
            float dot_ak_uj = 0.0f;
            float dot_uj_uj = 0.0f;
            for (size_t i = 0; i < M.cols; i++) dot_ak_uj += matrix_at(M, i, k)*matrix_at(u, i, j);
            for (size_t i = 0; i < M.cols; i++) dot_uj_uj += matrix_at(u, i, j)*matrix_at(u, i, j);
            for (size_t i = 0; i < M.cols; i++) matrix_at(u, i, k) -= dot_ak_uj/dot_uj_uj*matrix_at(u, i, j);
        }
        
        float norm_uk = 0.0f;
        for (size_t i = 0; i < M.cols; i++) norm_uk += matrix_at(u, i, k)*matrix_at(u, i, k);
        norm_uk = sqrtf(norm_uk);
        for (size_t i = 0; i < M.rows; i++) {
            matrix_at(*Q, i, k) = matrix_at(u, i, k)/norm_uk;
        }
    }

    for (size_t i = 0; i < M.rows; i++) {
        for (size_t j = i; j < M.cols; j++) {
            float dot_ei_aj = 0.0f;
            for (size_t k = 0; k < M.cols; k++) dot_ei_aj += matrix_at(M, k, j)*matrix_at(*Q, k, i);
            matrix_at(*R, i, j) = dot_ei_aj;
        }
    }

    free_matrix(u);
}

#define TRIANGULAR_TOL 1e-2
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

float* find_eigenvalues(Matrix M)
{
    bool solved;
    Matrix *M_current;
    Matrix Q, R, M_next;

    if (M.rows != M.cols) {
        fprintf(stderr, "ERROR: finding eigenvalue of non-square matrix with %zu rows and %zu cols\n", M.rows, M.cols);
        return NULL;
    }

    Q = zero_matrix(M.rows, M.cols);
    R = zero_matrix(M.rows, M.cols);
    M_next = zero_matrix(M.rows, M.cols);
    
    float *eigenvalues = malloc(M.cols*sizeof(float));
    if (eigenvalues == NULL) {
        fprintf(stderr, "ERROR: could not allocate memory for %zu eigenvalues\n", M.rows);
        return NULL;
    }

    size_t iters = 0;
    const size_t max_iters = 100;
    M_current = &M;
    do {
        QR_decomposition(*M_current, &Q, &R);
        matrix_multiplication(R, Q, &M_next);
        if (is_triangular(M_next)) {
            solved = true;
            break;
        } else {
            M_current = &M_next;
        }
        iters += 1;
    } while (iters < max_iters);

    if (solved) {
        for (size_t i = 0; i < M.cols; i++) {
            eigenvalues[i] = matrix_at(M_next, i, i);
        }
        free_matrix(Q);
        free_matrix(R);
    }
    free_matrix(M_next);

    return eigenvalues;
}

// TODO: partial pivoting
void LU_decomposition(Matrix A, Matrix *L, Matrix *U)
{
    if (A.rows != A.cols) {
        fprintf(stderr, "ERROR: LU decomposition can only be performed for square matrices, found rows = %zu, cols = %zu\n", A.rows, A.cols);
    }

    for (size_t i = 0; i < A.rows; i++) {
        matrix_at(*L, i, i) = 1.0f;
    }
    
    for (size_t j = 0; j < A.rows; j++) {
        for (size_t i = 0; i <= j; i++) {
            matrix_at(*U, i, j) = matrix_at(A, i, j);
            if (i > 0) {
                for (size_t k = 0; k < i; k++) {
                    matrix_at(*U, i, j) -= matrix_at(*L, i, k)*matrix_at(*U, k, j);
                }
            }
        }

        for (size_t i = j; i < A.rows; i++) {
            matrix_at(*L, i, j) = matrix_at(A, i, j);
            for (int k = 0; k < (int) j; k++) {
                matrix_at(*L, i, j) -= matrix_at(*L, i, k)*matrix_at(*U, k, j);
            }
            matrix_at(*L, i, j) /= matrix_at(*U, j, j);
        }
    }
}

void print_matrix(Matrix A, const char* matrix_name)
{
    printf("%s = \n", matrix_name);
    for (size_t i = 0; i < A.rows; i++) {
        for (size_t j = 0; j < A.cols; j++) {
            printf("%f ", matrix_at(A, i, j));
        }
        printf("\n");
    }
}

// Inspired by https://home.cc.umanitoba.ca/~farhadi/Math2120/Inverse%20Using%20LU%20decomposition.pdf
void inverse_matrix(Matrix A, Matrix *A_inv)
{
    Matrix L = zero_matrix(A.rows, A.cols);
    Matrix U = zero_matrix(A.rows, A.cols);
    LU_decomposition(A, &L, &U);

    // Invert U
    Matrix U_inv = identity(U.rows);
    for (size_t j = U.cols; j-- > 0;) {
        for (size_t i = j+1; i-- > 0;) {
            if (i == j) {
                if (matrix_at(U, i, j) == 0.0f) {
                    fprintf(stderr, "ERROR: could not invert matrix A because it's singular\n");
                    free_matrix(L);
                    free_matrix(U);
                    free_matrix(U_inv);
                }
                for (size_t k = 0; k < U.cols; k++) {
                    matrix_at(U_inv, i, k) /= matrix_at(U, i, j);
                }
                for (size_t k = 0; k < U.cols; k++) {
                    matrix_at(U, i, k) /= matrix_at(U, i, j);
                }
            } else {
                for (size_t k = 0; k < A.cols; k++) {
                    matrix_at(U_inv, i, k) -= matrix_at(U, i, j)*matrix_at(U_inv, j, k);
                }
                for (size_t k = 0; k < A.cols; k++) {
                    matrix_at(U, i, k) -= matrix_at(U, i, j)*matrix_at(U, j, k);
                }
            }
        }
    }
    
    // Invert L
    Matrix L_inv = identity(L.rows);
    for (size_t j = 0; j < L.cols; j++) {
        for (size_t i = j; i < L.rows; i++) {
            if (i == j) {
                for (size_t k = 0; k < U.cols; k++) {
                    matrix_at(L_inv, i, k) /= matrix_at(L, i, j);
                }
                for (size_t k = 0; k < U.cols; k++) {
                    matrix_at(L, i, k) /= matrix_at(L, i, j);
                }
            } else {
                for (size_t k = 0; k < A.cols; k++) {
                    matrix_at(L_inv, i, k) -= matrix_at(L, i, j)*matrix_at(L_inv, j, k);
                }
                for (size_t k = 0; k < A.cols; k++) {
                    matrix_at(L, i, k) -= matrix_at(L, i, j)*matrix_at(L, j, k);
                }
            }
        }
    }

    // Multiply (U^-1)(L^-1)
    matrix_multiplication(U_inv, L_inv, A_inv);
    free_matrix(L);
    free_matrix(U);
    free_matrix(U_inv);
}

// Via LU decomposition, from Numerical Recipes in C, 2nd edition, chapter 2.3
Matrix solve_linear_system(Matrix A, Matrix B, size_t n_dim)
{
    if (B.cols != 1) {
        fprintf(stderr, "ERROR: Matrix B must be column matrix, it has %zu rows and %zu cols\n", B.rows, B.cols);
        return (Matrix) {0};
    }

    Matrix X = zero_matrix(B.rows, 1);
    Matrix L = zero_matrix(A.rows, A.cols);
    Matrix U = zero_matrix(A.rows, A.cols);
    LU_decomposition(A, &L, &U);

    for (size_t i = 0; i < U.cols; i++) {
        if (matrix_at(U, i, i) == 0.0f) {
            fprintf(stderr, "ERROR: coefficient matrix of linear system is singular, system cannot be solved\n");
            free_matrix(L);
            free_matrix(U);
            free_matrix(X);
            return (Matrix) {0};
        }
    }

    // Solve L*U*X = B
    float ys[n_dim];
    ys[0] = matrix_at(B, 0, 0)/matrix_at(L, 0, 0);
    for (size_t i = 1; i < n_dim; i++) {
        ys[i] = matrix_at(B, i, 0);
        for (size_t j = 0; j < i; j++) {
            ys[i] -= matrix_at(L, i, j)*ys[j];
        }
        ys[i] /= matrix_at(L, i, i);
    }

    matrix_at(X, n_dim-1, 0) = ys[n_dim-1]/matrix_at(U, n_dim-1, n_dim-1);
    for (int i = n_dim-2; i >= 0; i--) {
        matrix_at(X, i, 0) = ys[i];
        for (size_t j = i+1; j < n_dim; j++) {
            matrix_at(X, i, 0) -= matrix_at(U, i, j)*matrix_at(X, j, 0);
        }
        matrix_at(X, i, 0) /= matrix_at(U, i, i);
    }
    free_matrix(L);
    free_matrix(U);
    return X;
}

float determinant(Matrix A)
{
    if (A.rows != A.cols) {
        fprintf(stderr, "ERROR: determinant is only defined for square matrices, found rows = %zu, cols = %zu\n", A.rows, A.cols);
        return 0.0f;
    }

    Matrix L = zero_matrix(A.rows, A.cols);
    Matrix U = zero_matrix(A.rows, A.cols);
    LU_decomposition(A, &L, &U);
    

    float result = 1;
    for (size_t i = 0; i < U.rows; i++) {
        result *= matrix_at(U, i, i);
    }
    free_matrix(L);
    free_matrix(U);
    return result;

}

float trapezoid(float f(float), float a, float b, size_t N)
{
    float dx = (b - a)/N;
    float result = (f(a) + f(b))/2;
    for (size_t i = 1; i < N; i++) {
        float x_i = a + i*dx;
        result += f(x_i);
    }
    result *= dx;
    return result;
}

float trapezoid_args(float f(float, Params), float a, float b, size_t N, Params params)
{
    float dx = (b - a)/N;
    float result = (f(a, params) + f(b, params))/2;
    for (size_t i = 1; i < N; i++) {
        float x_i = a + i*dx;
        result += f(x_i, params);
    }
    result *= dx;
    return result;
}

float* odesolve(float f(float, float), float y_0, float x_0, float x_final, size_t nsteps)
{
    float *ys = malloc((nsteps + 1) * sizeof(float));
    if (ys == NULL) {
        fprintf(stderr, "ERROR in odesolve: could not allocate array for solution\n");
        return NULL;
    }

    float dx = (x_final - x_0)/nsteps;
    ys[0] = y_0;
    for (size_t i = 1; i <= nsteps; i++) {
        float x_i = x_0 + (i-1)*dx;
        float y_prime = f(ys[i-1], x_i);
        ys[i] = ys[i-1] + y_prime*dx;
    }

    return ys;
}

float **odesolve_nd(void y_prime(float*, float, float*), float* y_initial, float x_initial, float x_final, size_t n_steps, size_t n_dim)
{
    float **result = (float**) malloc((n_steps+1)*sizeof(float*));
    if (result == NULL) {
        fprintf(stderr, "ERROR: could not allocate result of %zu dimensional ODE with %zu\n", n_dim, n_steps);
        return NULL;
    }
    for (size_t i = 0; i < n_steps+1; i++) {
        result[i] = (float*) malloc(n_dim*sizeof(float));
        if (result[i] == NULL) {
            fprintf(stderr, "ERROR: could not allocate result of %zu dimensional ODE with %zu\n", n_dim, n_steps);
            return NULL;
        }
    }

    for (size_t i = 0; i < n_dim; i++) {
        result[0][i] = y_initial[i];
    }

    float x;
    const float dx = (x_final - x_initial)/n_steps;
    float y_primes[n_dim];
    for (size_t i = 0; i < n_steps; i++) {
        x = x_initial + i*dx;
        y_prime(result[i], x, y_primes);
        for (size_t j = 0; j < n_dim; j++) {
            result[i+1][j] = result[i][j] + y_primes[j]*dx;
        }
    }
    return result;
}

LinearInterp get_linear_interpolator(float *xs, float *ys, size_t n_points)
{
    size_t n_segments = n_points - 1;
    float *as = malloc(n_segments*sizeof(float));
    float *bs = malloc(n_segments*sizeof(float));
    if (as == NULL || bs == NULL) {
        fprintf(stderr, "ERROR: could not allocate LinearInterp\n");
        return (LinearInterp) {0};
    }
    for (size_t i = 0; i < n_segments; i++) {
        float dx = xs[i+1] - xs[i];
        if (dx < 0) {
            fprintf(stderr, "ERROR: array `xs` for linear interpolation must be sorted\n");
            free(as);
            free(bs);
            return (LinearInterp) {0};
        }
        float dy = ys[i+1] - ys[i];
        as[i] = dy/dx;
        bs[i] = ys[i] - as[i]*xs[i];
    }
    return (LinearInterp) {.as = as, .bs = bs, .xs = xs, .n_points = n_points};
}

float eval_linear_interpolator(LinearInterp interp, float x)
{
    size_t index = 0;
    if (x < interp.xs[0]) {
        printf("WARNING: interpolation function evaluated at x = %f which is below x_min = %f, extrapolating linearly.\n", x, interp.xs[0]);
    } else if (x > interp.xs[interp.n_points-1]) {
        printf("WARNING: interpolation function evaluated at x = %f which is above x_max = %f, extrapolating linearly.\n", x, interp.xs[interp.n_points-1]);
        index = interp.n_points - 2;
    } else {
        // Find index
        for (size_t i = 1; i < interp.n_points; i++) {
            if (x < interp.xs[i]) {
                index = i-1;
                break;
            }
        }
    }
    return interp.as[index]*x + interp.bs[index];
}

CubicSpline get_cubic_spline(float* xs, float* ys, size_t n_points)
{
    size_t n_segments = n_points - 1;
    float *as = malloc(n_points*sizeof(float));
    float *cs = malloc(n_points*sizeof(float));
    float *bs = malloc(n_segments*sizeof(float));
    float *ds = malloc(n_segments*sizeof(float));

    float hs[n_segments];
    float alphas[n_segments];
    float ls[n_points];
    float mus[n_points];
    float zs[n_points];

    if (as == NULL || bs == NULL || cs == NULL || ds == NULL) {
        fprintf(stderr, "ERROR: could not allocate `CubicSpline` object\n");
        return (CubicSpline) {0};
    }

    for (size_t i = 0; i < n_points; i++) {
        as[i] = ys[i];
    }

    for (size_t i = 0; i < n_segments; i++) {
        hs[i] = xs[i+1] - xs[i];
        if (hs[i] < 0) {
            fprintf(stderr, "ERROR: array `xs` for Cubic Spline must be sorted\n");
            free(as);
            free(bs);
            free(cs);
            free(ds);
            return (CubicSpline) {0};
        }
    }

    for (size_t i = 1; i < n_segments; i++) {
        alphas[i] = 3.0f / hs[i] * (as[i+1] - as[i]) - 3.0f / hs[i-1] * (as[i] - as[i-1]);
    }

    ls[0] = 1.0f;
    mus[0] = 0.0f;
    zs[0] = 0.0f;

    for (size_t i = 1; i < n_segments; i++) {
        ls[i] = 2.0f*(xs[i+1] - xs[i-1]) - hs[i-1]*mus[i-1];
        mus[i] = hs[i]/ls[i];
        zs[i] = (alphas[i] - hs[i-1]*zs[i-1])/ls[i];
    }

    ls[n_segments] = 1.0f;
    cs[n_segments] = 0.0f;
    zs[n_segments] = 0.0f;

    for (int j = n_segments - 1; j >= 0; j--) {
        cs[j] = zs[j] - mus[j]*cs[j+1];
        bs[j] = (as[j+1] - as[j])/hs[j] - hs[j]*(cs[j+1] + 2.0f*cs[j])/3.0f;
        ds[j] = (cs[j+1] - cs[j])/(3.0f*hs[j]);
    }
    return (CubicSpline) { .as = as, .bs = bs, .cs = cs, .ds = ds, .xs = xs, .n_points = n_points };
}

float eval_cubic_spline(CubicSpline interp, float x)
{
    size_t index = 0;
    float delta;
    if (x < interp.xs[0]) {
        printf("WARNING: interpolation function evaluated below x_min = %f, extrapolating with cubic spline (personally i think that's not a good idea).\n", interp.xs[0]);
        delta = (x - interp.xs[0]);
    } else if (x > interp.xs[interp.n_points-1]) {
        printf("WARNING: interpolation function evaluated above x_max = %f, extrapolating with cubic spline (personally i think that's not a good idea).\n", interp.xs[interp.n_points-1]);
        delta = x - interp.xs[interp.n_points-1];
        index = interp.n_points - 2;
    } else {
        // Find index
        for (size_t i = 1; i < interp.n_points; i++) {
            if (x < interp.xs[i]) {
                index = i-1;
                break;
            }
        }
        delta = x - interp.xs[index];
    }
    return interp.as[index] + interp.bs[index]*delta + interp.cs[index]*delta*delta + interp.ds[index]*delta*delta*delta;
}

float find_root(float f(float), float x_min, float x_max, size_t max_iters, float tol)
{
    float y_min = f(x_min);
    float y_max = f(x_max);

    float x_min_in = x_min;
    float x_max_in = x_max;

    if (copysign(1, y_min*y_max) > 0.0f) {
        fprintf(stderr, "ERROR: root finding between x = %f and x = %f expected to have different signs of f(x), found f(%f) = %f and f(%f) = %f, returning zero.\n", x_min, x_max, x_min, y_min, x_max, y_max);
        return 0.0f;
    }

    float root;
    float y_root;
    size_t iters = 0;
    do {
        float a = (y_max - y_min)/(x_max - x_min);
        float b = y_min - a*x_min;
        root = -b/a;
        y_root = f(root);
        if (copysign(1, y_root*y_min) > 0.0f) {
            y_min = y_root;
            x_min = root;
        } else {
            y_max = y_root;
            x_max = root;
        }
    } while (iters++ < max_iters && fabsf(f(root)) > tol);
    
    if (fabsf(f(root)) > tol) printf("WARNING: root finding between %f and %f did not achieve desired tolerance %f\n", x_min_in, x_max_in, tol);
    return root;
}

// Gradient descent with finite difference
float minimize(float f(float), float x, float rate, size_t max_iters)
{
    float x_min = x;
    float derivative;
    static const float tol = 1e-5;
    static const float eps = 1e-4;
    for (size_t i = 0; i < max_iters; i++) {
        derivative = (f(x_min + eps) - f(x_min))/eps;
        if (fabsf(derivative*rate) < fabsf(tol*x_min)) break;
        x_min -= derivative*rate;
    }
    return x_min;
}

// Blatantly copy-pasted from https://rosettacode.org/wiki/Fast_Fourier_transform#C
void _fft(float complex *buf, float complex *aux, int n, int step)
{
	if (step < n) {
		_fft(aux, buf, n, step * 2);
		_fft(aux + step, buf + step, n, step * 2);
 
		for (int i = 0; i < n; i += 2 * step) {
			float complex t = cexp(-I * M_PI * i / n) * aux[i + step];
			buf[i / 2]     = aux[i] + t;
			buf[(i + n)/2] = aux[i] - t;
		}
	}
}

// Blatantly copy-pasted from https://rosettacode.org/wiki/Fast_Fourier_transform#C
void fft(float complex *in, float complex *out, int n)
{
	float complex aux[n];
    
    memcpy(aux, in, n*sizeof(float complex));
    memcpy(out, in, n*sizeof(float complex));
 
	_fft(out, aux, n, 1);
}

void free_matrix(Matrix M)
{
    free(M.elements);
}

void free_linear_interpolator(LinearInterp interp)
{
    free(interp.as);
    free(interp.bs);
}

void free_cubic_spline(CubicSpline interp)
{
    free(interp.as);
    free(interp.bs);
    free(interp.cs);
    free(interp.ds);
}

#endif // INTEGRATE_IMPLEMENTATION
#endif // NUMERICAL_H_