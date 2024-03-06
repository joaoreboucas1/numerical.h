#ifndef INTEGRATE_H_
#define INTEGRATE_H_

#include <stddef.h>
#include <stdlib.h>

typedef struct {
    int *int_args;
    float *float_args;
} Params;

// TODO: LinearInterp doesn't really need to precompute as and bs
// Refactor to { float *xs, float *ys, size_t n_points }
typedef struct {
    float *as;
    float *bs;
    float *xs;
    size_t n_points;
} LinearInterp;


float trapezoid(float f(float), float a, float b, size_t N); // Integrate f from a to b with N points
float trapezoid_args(float f(float, Params), float a, float b, size_t N, Params params); // Integrate function with parameters
float* odesolve(float f(float, float), float y_0, float x_0, float x_final, size_t nsteps); // Solve 1D ODE y'(y, x) = f(y, x) with initial conditions (x_0, y_0) and nsteps steps
LinearInterp get_interpolator(float* xs, float* ys, size_t n_points);
float eval_interpolator(LinearInterp interp, float x);


#ifdef NUMERICAL_IMPLEMENTATION
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


LinearInterp get_interpolator(float *xs, float *ys, size_t n_points)
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
        float dy = ys[i+1] - ys[i];
        as[i] = dy/dx;
        bs[i] = ys[i] - as[i]*xs[i];
    }
    return (LinearInterp) {.as = as, .bs = bs, .xs = xs, .n_points = n_points};
}

float eval_interpolator(LinearInterp interp, float x)
{
    if (x < interp.xs[0]) {
        printf("WARNING: interpolation function evaluated below x_min = %f, extrapolating linearly.\n", interp.xs[0]);
        return interp.as[0]*x + interp.bs[0];
    } else if (x > interp.xs[interp.n_points-1]) {
        printf("WARNING: interpolation function evaluated above x_max = %f, extrapolating linearly.\n", interp.xs[interp.n_points-1]);
        return interp.as[interp.n_points - 2]*x + interp.bs[interp.n_points - 2];
    } else {
        // Find index
        size_t index = 0;
        for (size_t i = 1; i < interp.n_points; i++) {
            if (x < interp.xs[i]) {
                index = i-1;
                printf("a = %f, b = %f\n", interp.as[index], interp.bs[index]);
                break;
            }
        }
        return interp.as[index]*x + interp.bs[index];
    }
}


#endif // INTEGRATE_IMPLEMENTATION

#endif // INTEGRATE_H_