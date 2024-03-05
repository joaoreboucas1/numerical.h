#ifndef INTEGRATE_H_
#define INTEGRATE_H_

#include <stddef.h>
#include <stdlib.h>

typedef struct {
    int *int_args;
    float *float_args;
} Params;

float trapezoid(float f(float), float a, float b, size_t N); // Integrate f from a to b with N points
float trapezoid_args(float f(float, Params), float a, float b, size_t N, Params params); // Integrate function with parameters
float* odesolve(float f(float, float), float y_0, float x_0, float x_final, size_t nsteps); // Solve 1D ODE y'(y, x) = f(y, x) with initial conditions (x_0, y_0) and nsteps steps

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


#endif // INTEGRATE_IMPLEMENTATION

#endif // INTEGRATE_H_