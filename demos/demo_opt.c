#include <stdio.h>
#include <math.h>

#define NUMERICAL_IMPLEMENTATION
#include "../numerical.h"

int main()
{
    float x = 3.0f;
    float x_min = minimize(sinf, x, 1.0f, 10);
    float x_exact = 3*M_PI/2;
    float error = 100*(x_min - x_exact)/x_exact;
    printf("sin(x) has a minimum at value %.3f (error = %f%%)\n", x_min, error);
    return 0;
}