#ifndef COMPUTEFOURIER_H
#define COMPUTEFOURIER_H

#include "fft.h"

#include <complex.h>
#include "fftw.h"
#include "filters.h"

#define OPTIMIZE_FFTW 0

extern bool_t WITH_COMB;
extern bool_t ALGORITHM1;
extern bool_t VERBOSE;
extern bool_t TIMING;

//Comments located in the cc file.
  Node *outer_loop(complex_t *origx, int n, Filter *filter, int B2, int num, 
      int B, int W_Comb, int Comb_loops, int loop_threshold, int location_loops, int loops, int *hits_found, int *score, int*hits);

void inner_loop_step_a_plus_b(complex_t *origx, Filter *filter, complex_t *x_sampt, int*permute, int n, int B, int loops, double *PF_ALL, double *B_ALL, float* DtoH, float* HtoD, unsigned int plan);

#endif
