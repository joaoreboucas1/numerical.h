#include <stdio.h>
#include <math.h>

#define NUMERICAL_IMPLEMENTATION
#include "../numerical.h"

int main()
{
    float x_min = 2.5f;
    float x_max = 3.5f;
    size_t max_iters = 10;
    float tol = 1e-7;
    float root = find_root(sinf, x_min, x_max, max_iters, tol);
    
    printf("Solution of sin(x) = 0 between %f and %f: %.10f\n", x_min, x_max, root);
    return 0;
}