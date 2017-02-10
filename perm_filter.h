#ifndef PERM_FILTER_H
#define PERM_FILTER_H

#include "fft.h"
#include "cufft.h"
#include "filters.h"

void inner_loop_step_a_plus_b(complex_t *origx, Filter *filter, complex_t *x_sampt, int*permute, int n, int B, double *PF_ALL, double *B_ALL, float *DtoH, float* HtoD, unsigned int plan);


#endif
