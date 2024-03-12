#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NUMERICAL_IMPLEMENTATION
#include "../numerical.h"

float y_prime(float y, float x)
{
    (void) x;
    return y;
}

int main()
{
    float x_0 = 0.0f;
    float x_final = 1.0f;
    float y_0 = 1.0f;
    size_t nsteps = 100;
    float *ys = odesolve(y_prime, y_0, x_0, x_final, nsteps);
    if (ys == NULL) {
        exit(1);
    }
    float dx = (x_final - x_0)/nsteps;
    printf("Numerical solution of y'(x) = y with IC y(0) = 0 and error versus exact solution:\n");
    for (size_t i = 0; i <= nsteps; i++) {
        float x_i = x_0 + i*dx;
        float y_exact = expf(x_i);
        float error = 100*fabsf(ys[i] - y_exact)/y_exact;
        printf("x = %f => y = %f (error = %.2f%%)\n", x_i, ys[i], error);
    }
    return 0;
}