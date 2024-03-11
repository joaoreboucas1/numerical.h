#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <string.h>

#define NUMERICAL_IMPLEMENTATION
#include "../numerical.h"

int main()
{
    const size_t N = pow(2, 10);
    float complex y[N];

    // Create signal as a mixture of sines with different frequencies
    #define n_freqs 6
    float frequencies[n_freqs] = {3.0f, 4.0f, 10.0f, 69.0f, 420.0f, 666.0f};
    for (size_t i = 0; i < N; i++) {
        y[i] = 0;
        for (size_t j = 0; j < n_freqs;j++) {
            y[i] += csinf((float complex) 2*M_PI*frequencies[j]*i/N);
        }
    }

    // Detect the frequencies via FFT
    float complex *y_transf = fft(y, N);
    for (size_t i = 0; i < N; i++) {
        float magnitude = cabsf(y_transf[i]);
        if (magnitude > 1.0f) {
            printf("%zu: y = %f, re = %f, im = %f, magnitude = %f\n", i, crealf(y[i]), crealf(y_transf[i]), cimagf(y_transf[i]), cabsf(y_transf[i]));
        }
    }
    return 0;
}