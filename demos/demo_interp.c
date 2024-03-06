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
    LinearInterp interp = get_interpolator(xs, ys, N+1);

    // Evaluate interpolator
    float x = M_PI/3.0f;
    float result = eval_interpolator(interp, x);
    float exact = sinf(x);
    float error = 100*fabsf(result - exact)/exact;
    printf("sin(%.2f) = %f (error = %.2f%%)\n", x, result, error);
}