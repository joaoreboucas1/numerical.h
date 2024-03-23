#include <stdio.h>

#define NUMERICAL_IMPLEMENTATION
#include "../numerical.h"

void y_prime(float *y, float x, float* result)
{
    // Equivalent to the system:
    // x' = y
    // y' = -x
    (void) x;
    result[0] = y[1];
    result[1] = -y[0];
}

int main()
{
    float y_initial[2] = { 1.0f, 0.0f };
    const size_t n_steps = 100;
    float **result = odesolve_nd(y_prime, y_initial, 0.0f, 5.0f, n_steps, 2);
    for (size_t i = 0; i < n_steps; i++) {
        float x = i*5.0f/n_steps;
        printf("x_%zu = %f => y_%zu = (%.3f, %.3f)\n", i, x, i, result[i][0], result[i][1]);
    }
    return 0;
}