#ifndef CUDA_FFT_H
#define CUDA_FFT_H

#include "fft.h"

void cuda_fft(complex_t *x, complex_t *x_f, int n, int repetitions, float *cu_fft_time);
void cufft_plan_create(unsigned int* plan, int B, int loops);
#endif
