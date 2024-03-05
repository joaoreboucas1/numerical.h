#include <stdio.h>
#include <math.h>

#define NUMERICAL_IMPLEMENTATION
#include "../numerical.h"

// Example functions
float f(float x)
{
    return x*x;
}

float g(float x)
{
    return x*x*x;
}

float h(float x)
{
    return expf(-x*x);
}


// Function with additional parameters
float j(float x, Params p)
{
    float alpha = p.float_args[0];
    int n = p.int_args[0];
    return alpha*pow(x, n);
}


int main()
{
    float a = 0;
    float b = 1;
    size_t N = 100;
    float f_integral = trapezoid(f, a, b, N);
    float g_integral = trapezoid(g, a, b, N);
    float h_integral = trapezoid(h, a, b, N);
    printf("Integral of f(x) = x^2 from %.2f to %.2f: %f\n", a, b, f_integral);
    printf("Integral of g(x) = x^3 from %.2f to %.2f: %f\n", a, b, g_integral);
    printf("Integral of h(x) = exp(-x^2) from %.2f to %.2f: %f\n", a, b, h_integral);
    
    float alpha = 69.0f;
    int n = 2;
    Params p;
    p.float_args = &alpha;
    p.int_args = &n;
    float j_integral = trapezoid_args(j, a, b, N, p);
    printf("Integral of j(x) = %.2f*x^%d from %.2f to %.2f: %f\n", alpha, n, a, b, j_integral);
    return 0;
}