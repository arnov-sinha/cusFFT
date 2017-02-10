#include <cuda.h>
#include "cufft.h"
#include "fft.h"
#include "utils.h"
extern "C"
{
#include "timer.h"
#include "cuda_fft.h"
}

extern "C"
void cuda_fft(complex_t *x, complex_t *x_f, int n, int repetitions, float* time)
{
  
  reset_timer(); 
  double DDD = get_time();
  cufftHandle plan;
  cufftResult err = cufftPlan1d(&plan, n, CUFFT_Z2Z, 1); 
  
  if (err != CUFFT_SUCCESS){
    fprintf(stderr, "CUFFT error: Execution failed, error code is %d\n", err);
    exit(-1); 
  }

  printf("Time to create cuFFT plan: %lf\n", get_time()-DDD);
 
  //cudaEvent_t start, stop;
  //cudaEventCreate(&start);
  //cudaEventCreate(&stop); 

  //cudaEventRecord(start);
  
  complex_t *cufft_x_f = (complex_t *)malloc(n*sizeof(complex_t)); 
  complex_t *d_x, *d_x_f;
  
  cudaMalloc((void**)&d_x, n*sizeof(complex_t));
  cudaMemcpy(d_x, x, n*sizeof(complex_t),cudaMemcpyHostToDevice);
  
  cudaMalloc((void**)&d_x_f, n*sizeof(complex_t));

  for(int i = 0; i < repetitions; i++){
    err = cufftExecZ2Z(plan, (cufftDoubleComplex *)d_x, (cufftDoubleComplex *)d_x_f, CUFFT_FORWARD);
    if (err != CUFFT_SUCCESS){
      fprintf(stderr, "CUFFT error: Execution failed, error code is %d\n", err);
      exit(-1); 
    }
  }
  
  cudaMemcpy(cufft_x_f, d_x_f, n*sizeof(complex_t), cudaMemcpyDeviceToHost);

  cudaFree(d_x);
  cudaFree(d_x_f);
  
  //cudaEventRecord(stop);
  //cudaEventSynchronize(stop);
  //cudaEventElapsedTime(time, start, stop);
  //
  //cudaEventDestroy(start);
  //cudaEventDestroy(stop);
  printf("Time to run cuFFT : %f\n", get_time());
  
  real_t CUFFT_ERROR =0;
  for(int i=0; i< n ; i++){
      CUFFT_ERROR += cabs(cufft_x_f[i]/n- x_f[i]);
  }
  printf("ERROR of CUFFT is %lg\n", CUFFT_ERROR);
  cufftDestroy(plan);
  free(cufft_x_f);
  
}

extern "C"
void cufft_plan_create(unsigned int* plan, int B, int loops)
{
   //cudaFree(0);
   cufftHandle *plan1 = (cufftHandle*)plan;
   cufftResult err;
   err = cufftPlan1d(plan1, B, CUFFT_Z2Z, loops); 
   if (err != CUFFT_SUCCESS){
     fprintf(stderr, "CUFFT error: Plan creation failed");
     exit(-1); 
   }
}
