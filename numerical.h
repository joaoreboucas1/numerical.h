#ifndef INTEGRATE_H_
#define INTEGRATE_H_

#include <stddef.h>
#include <stdlib.h>

// Parameters to be passed into functions for numerical integration
typedef struct {
    int *int_args;
    float *float_args;
} Params;

typedef struct {
    float *as;
    float *bs;
    float *xs;
    size_t n_points;
} LinearInterp;

// Following the implementation from https://en.wikipedia.org/wiki/Spline_(mathematics)#Algorithm_for_computing_natural_cubic_splines
typedef struct {
    float *as;
    float *bs;
    float *cs;
    float *ds;
    float *xs;
    size_t n_points;
} CubicSpline;

// Numerical integration
float trapezoid(float f(float), float a, float b, size_t N); // Integrate f from a to b with N points
float trapezoid_args(float f(float, Params), float a, float b, size_t N, Params params); // Integrate function with parameters

// ODE integration
float* odesolve(float f(float, float), float y_0, float x_0, float x_final, size_t nsteps); // Solve 1D ODE y'(y, x) = f(y, x) with initial conditions (x_0, y_0) and nsteps steps

// Interpolation
LinearInterp get_linear_interpolator(float* xs, float* ys, size_t n_points); // Create a `LinearInterpolator` object from the data, representing the function y(x)
float eval_linear_interpolator(LinearInterp interp, float x); // Evaluate a `LinearInterpolator` object at `x`
CubicSpline get_cubic_spline(float* xs, float* ys, size_t n_points); // Create a `CubicSpline` object from the data, representing the function y(x)
float eval_cubic_spline(CubicSpline interp, float x); // Evaluate a `CubicSpline` object at `x`

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
        printf("WARNING: interpolation function evaluated below x_min = %f, extrapolating linearly.\n", interp.xs[0]);
    } else if (x > interp.xs[interp.n_points-1]) {
        printf("WARNING: interpolation function evaluated above x_max = %f, extrapolating linearly.\n", interp.xs[interp.n_points-1]);
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

#endif // INTEGRATE_IMPLEMENTATION
#endif // INTEGRATE_H_