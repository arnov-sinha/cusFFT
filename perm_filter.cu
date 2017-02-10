#include "fft.h"
#include "cufft.h"
#include "utils.h"
#include "perm_filter.h"
extern "C"
{
#include "timer.h"
}
#include "cuComplex.h"
#include <cuda.h>

__global__ void PermFilterKernel(cuDoubleComplex* d_origx, cuDoubleComplex* d_filter, int* d_permute, cuDoubleComplex* d_x_sampt, int B, int n, int loops, int round)
{
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  if(i < loops*B)
  {
    int i_mod_B = i & (B-1);
    int off = i_mod_B * round;
    int tmp = i/B;
    int ai = d_permute[tmp];
    cuDoubleComplex tmp_value1, tmp_value2;

    for(int j=0; j<round; j++){
      int index = (i_mod_B + B*j)*ai &(n-1);
      
      tmp_value1 = cuCmul(d_origx[index],d_filter[off+j]);
      tmp_value2 = cuCadd(tmp_value1, tmp_value2);
    }
    d_x_sampt[i] = tmp_value2;
  }
}


extern "C"
void inner_loop_step_a_plus_b(complex_t *origx, Filter *filter, complex_t *x_sampt, int*permute, int n, int B, int loops, 
    double *PF_ALL, double *B_ALL, float *DtoH, float *HtoD, unsigned int plan1)
{
  int filter_size = filter->sizet;
  int round = filter_size/B;
  complex_t *d_origx, *d_filter, *d_x_sampt;
  int *d_permute;
  
  cufftHandle plan = (cufftHandle)plan1;
  cufftResult err;
  
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  //Start of timing on HtoD
  cudaEventRecord(start);
  
  //Allocate device memory and copy to device
  cudaMalloc((void**)&d_origx, n*sizeof(complex_t));
  cudaMemcpy(d_origx, origx, n*sizeof(complex_t), cudaMemcpyHostToDevice);

  cudaMalloc((void**)&d_filter, filter_size*sizeof(complex_t));
  cudaMemcpy(d_filter, filter->time, filter_size*sizeof(complex_t), cudaMemcpyHostToDevice);

  cudaMalloc((void**)&d_permute, loops*sizeof(int));
  cudaMemcpy(d_permute, permute, loops*sizeof(int), cudaMemcpyHostToDevice);

  cudaMalloc((void**)&d_x_sampt, loops*B*sizeof(complex_t));
 
  //End of timing on HtoD 
  cudaEventRecord(stop);
  cudaEventSynchronize(stop);
  float time;
  cudaEventElapsedTime(&time,start, stop);
  *HtoD = (float)time/1e3;

  //Start of the kernels
  //Start of timing on kernels
  cudaEventRecord(start);

  dim3 dimBlock(512);
  dim3 dimGrid(loops*B/dimBlock.x);
  PermFilterKernel<<<dimGrid, dimBlock>>>((cuDoubleComplex*)d_origx, (cuDoubleComplex*)d_filter, d_permute, (cuDoubleComplex*)d_x_sampt, B, n,loops,round);
  
  //End of timing on kernerls 
  cudaEventRecord(stop);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&time,start, stop);
  *PF_ALL = (float)time/1e3;
  
  //Step B -- cuFFT of B-dimensional FFT
  double DDD = get_time();
    
  err = cufftExecZ2Z(plan, (cufftDoubleComplex *)d_x_sampt, (cufftDoubleComplex *)d_x_sampt, CUFFT_FORWARD);
  if (err != CUFFT_SUCCESS){
    fprintf(stderr, "CUFFT error: Execution failed, error code is %d\n", err);
    exit(-1); 
  }
  
  *B_ALL = get_time() - DDD;

  
  //Transfer back the d_x_sampt in freq domain
  //Start of timing on DtoH
  cudaEventRecord(start);

  cudaMemcpy(x_sampt, d_x_sampt, loops*B*sizeof(complex_t), cudaMemcpyDeviceToHost); 
  
  //End of timing on DtoH 
  cudaEventRecord(stop);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&time,start, stop);
  *DtoH = (float)time/1e3;
  

  //destroy plan and device memory
  cufftDestroy(plan);
  cudaFree(d_origx);
  cudaFree(d_filter);
  cudaFree(d_x_sampt);
  cudaFree(d_permute);

  cudaEventDestroy(start);
  cudaEventDestroy(stop);
}

