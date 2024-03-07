#include <stdio.h>
#include <math.h>

#define NUMERICAL_IMPLEMENTATION
#include "../numerical.h"

int main()
{
    // Define data on regular grid
    size_t N = 100;
    float x_0 = -5.0f;
    float x_f = 5.0f;
    float dx = (x_f - x_0)/N;
    float xs[N+1];
    float ys[N+1];
    for (size_t i = 0; i <= N; i++) {
        xs[i] = x_0 + i*dx;
        ys[i] = sinf(xs[i]);
    }

    // Get interpolator
    LinearInterp lin_interp = get_linear_interpolator(xs, ys, N+1);
    CubicSpline cs_interp = get_cubic_spline(xs, ys, N+1);

    // Evaluate interpolator
    float x = M_PI/3.0f;
    float result_lin = eval_linear_interpolator(lin_interp, x);
    float result_cs = eval_cubic_spline(cs_interp, x);
    float exact = sinf(x);
    float error_lin = 100*fabsf(result_lin - exact)/exact;
    float error_cs = 100*fabsf(result_cs - exact)/exact;
    printf("Linear interpolator:\n");
    printf("sin(%.2f) = %f (error = %.2f%%)\n", x, result_lin, error_lin);
    printf("Cubic Spline interpolator:\n");
    printf("sin(%.2f) = %f (error = %.2f%%)\n", x, result_cs, error_cs);
}