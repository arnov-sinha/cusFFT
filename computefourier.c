#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <sys/time.h>
#include <assert.h>
#include <complex.h>
#include <math.h>
#include <omp.h>
#include "utils.h"
#include "timer.h"
#include "computefourier.h"
#include "filters.h"
#include "fft.h"
#include "wrapper.h"

#include "cuda_fft.h"
bool_t WITH_COMB = false;
bool_t ALGORITHM1   = true;
bool_t VERBOSE      = false;
bool_t TIMING       = true;

#define vprintf(...) if (VERBOSE){ printf(__VA_ARGS__);}



/*use binary search to find the upper bound*/
inline int find_upper_bound(Pair *permuted_approved, int num_Comb, int val)
{ 
  int mid, low = 0, high = num_Comb - 1;
  
  if(val >= permuted_approved[high].first){
    return high;
  }
  
  mid = (low + high) / 2;
  while (high > low) {
    if (permuted_approved[mid].first >=val)
      high  = mid;
    else
      low = mid + 1;
    
    mid = (low + high) / 2;
  }
  return mid;
}

inline int timesmod(const int x, const int a, const int n) {
  return (int)((((long long int)(x))*a)&(n-1));
}


inline int Comb_Filt(complex_t *origx, int n, int num, int W_Comb, 
                      int* Comb_Approved, complex_t *x_sampt, real_t *samples, unsigned int offset, int sigma, fftw_plan p_comb)
{
  
  //Set y_i = x_[i(n/M)+sigma]
  for(int i = 0; i < W_Comb; i++){
    x_sampt[i] = origx[offset + i*sigma];
  }

  //Compute fftw of y_i
  fftw_execute_dft(p_comb, x_sampt, x_sampt); 

  for(int i = 0; i < W_Comb; i ++)
    samples[i] = cabs2(x_sampt[i]);
  
  //Find 2k largest elements of y_i in freq domain
  find_largest_indices(Comb_Approved, num, samples, W_Comb);
  
  return 0;
}



int inner_loop_step_c(int num, int B, int *J, complex_t *x_sampt, real_t *samples)
{
  
  double DDD=get_time();
  
  for(int i = 0; i < B; i++){
    samples[i] = cabs2(x_sampt[i]);
  }

  find_largest_indices(J, num, samples, B);

  if (TIMING) {
    vprintf("Step 1.C (LARGEST BUCKS): %lf\n", get_time() - DDD);
  }

  return 0;
}  

/* 
Find indices that map to J , i.e., lie within n/(2B) of (J * n/B) after permutation.
For each such i, increment score[i] and append to hits if score[i] reaches loop_threshold.
*/ 

int inner_loop_filter_regular(int *J, int n, int num, int B, int a, int ai, int b, int loop_threshold,
                              int *score, int *hits, int *hits_found){
  double DDD = get_time();

  // Given the set of large samples, find the locations in [n] that map there
  // and output them
  int my_hits_found= *hits_found;
  //#pragma omp parallel for
  for(int i = 0; i < num; i++){
    int low = ((int)(ceil((J[i] - 0.5) * n / B)) + n)&(n-1);
    int high = ((int)(ceil((J[i] + 0.5) * n / B)) + n)&(n-1);
    int loc = timesmod(low, a, n);
    for(int j = low; j != high; j = (j + 1)&(n-1)) {
      //#pragma omp atomic
      score[loc]++;
      if(score[loc]==loop_threshold){
        //#pragma omp critical
        hits[my_hits_found++]=loc;
      }
      loc = (loc + a)&(n-1);
    }
  }

  *hits_found = my_hits_found;

  if (TIMING) {
    vprintf("Step 1.D (GROUPING):----------------------------- %lf\n\n", get_time()-DDD);
    vprintf("#####################################################################\n\n");
  }
  
  return 0;
}


/*
 Find indices that (1) map to J under the permutation and (2) lie in Comb_Approved mod W_Comb.

 For each such i, increment hits[i] and append to hits_found if hits[i] reaches loop_threshold.
*/
int inner_loop_filter_Comb(int *J, int n, int num, int B, int a, int ai, int b, int loop_threshold,
                      int *score, int *hits, int *hits_found,
                      int *Comb_Approved,int num_Comb, int W_Comb, Pair *permuted_approved){
  
  
  double DDD = get_time();
  
  
  int my_hits= *hits_found; 
  #pragma omp parallel for
  for(int m=0; m< num_Comb; m++){
    int prev = timesmod(Comb_Approved[m], ai, W_Comb);
    permuted_approved[m].first  = prev;
    permuted_approved[m].second = timesmod(prev, a, n);
  }
  
  qsort(permuted_approved, num_Comb, sizeof(permuted_approved[0]), comp_struct3);

  int nthreads = omp_get_max_threads();

  // compute intersection of permuted_approved and indices close to J * n/B, then invert to get true locations.
  if(nthreads >4)
  {
    #pragma omp parallel for
    for(int i = 0; i < num; i++){
      int low = ((int)(ceil((J[i] - 0.5) * n / B)) + n) & (n-1);
      int high = ((int)(ceil((J[i] + 0.5) * n / B)) + n) & (n-1);
      
      int key = low % W_Comb;
      
      int index = find_upper_bound(permuted_approved, num_Comb, key);
      int location = low - (low % W_Comb);
      int locinv = timesmod(location, a, n);
      for(int j = index; ; j++){
        if (j == num_Comb){
          j -= num_Comb;
          location = (location + W_Comb) & (n-1);
          //locinv = (unsigned long)(location * a) & (n-1);
          locinv = timesmod(location, a, n);
        }
        int approved_loc = location + permuted_approved[j].first;
        if((low < high && (approved_loc >= high || approved_loc < low)) ||
           (low > high && (approved_loc >= high && approved_loc < low)))
          break;
        int loc = (locinv + permuted_approved[j].second) & (n-1);
        #pragma omp atomic
        score[loc]++;
        if(score[loc]==loop_threshold){
          #pragma omp critical 
          hits[my_hits++]=loc;
        }
      }
    }
  }else
  {
    for(int i = 0; i < num; i++){
      int low = ((int)(ceil((J[i] - 0.5) * n / B)) + n) & (n-1);
      int high = ((int)(ceil((J[i] + 0.5) * n / B)) + n) & (n-1);
      
      int key = low % W_Comb;
      
      int index = find_upper_bound(permuted_approved, num_Comb, key);
      int location = low - (low % W_Comb);
      int locinv = timesmod(location, a, n);
      for(int j = index; ; j++){
        if (j == num_Comb){
          j -= num_Comb;
          location = (location + W_Comb) & (n-1);
          locinv = timesmod(location, a, n);
        }
        int approved_loc = location + permuted_approved[j].first;
        if((low < high && (approved_loc >= high || approved_loc < low)) ||
           (low > high && (approved_loc >= high && approved_loc < low)))
          break;
        int loc = (locinv + permuted_approved[j].second) & (n-1);
        score[loc]++;
        if(score[loc]==loop_threshold){
          hits[my_hits++]=loc;
        }
      }
    }
    
  }
  
  *hits_found = my_hits;
  
  if (TIMING) {
    vprintf("Step 1.D (GROUPING):----------------------------- %lf\n\n", get_time()-DDD);
    vprintf("#####################################################################\n\n");
  }

  return 0;
}

/*
  hits contains the indices that we want to estimate.

  x_samp contains a B-dimensional array for each of the `loops`
  iterations of the outer loop.  Every coordinate i of x "hashes to" a
  corresponding coordinate (permute[j] * i) mod B of x_samp[j], which
  gives an estimate of x[i].

  We estimate each coordinate as the median (independently in real and
  imaginary axes) of its images in the rows of x_samp.
 */
void estimate_values(const int *hits, int *hits_found,
                complex_t *x_samp,  const int loops, int n,
                const int *permute,
                const int B, const Filter *filter, 
                int location_loops, Node *ans)
{
  int my_hits_found = *hits_found;
  int n_div_B = n/B;
  int n_div_B_div_2 = n/B/2;
  #pragma omp parallel
  {
    real_t *values[2];
    values[0] = (real_t *) malloc(loops * sizeof(*values[0]));
    values[1] = (real_t *) malloc(loops * sizeof(*values[1]));
    #pragma omp for nowait
    for(int i = 0; i < my_hits_found; i++){
      for(int j = 0; j < loops; j++){
        int permuted_index= timesmod(permute[j], hits[i],  n);
        int hashed_to = permuted_index / n_div_B;
        int dist = permuted_index % n_div_B;
        if (dist > n_div_B_div_2) {
          hashed_to = (hashed_to + 1)%B;
          dist -= n_div_B;
        }
        dist = (n - dist) & (n-1); 
        complex_t filter_value = filter->freq[dist];
        complex_t tmp = x_samp[j*B+hashed_to] / filter_value;
        values[0][j] = creal(tmp);
        values[1][j] = cimag(tmp);
      }
    
      int location = (loops - 1) / 2;
       
      
      #ifndef WRAPPER 
      //Use qsort to find the value in location
      qsort(values[0], loops, sizeof(*values[0]), compare);
      qsort(values[1], loops, sizeof(*values[1]), compare);
      #else
      //Use wrapper to find the nth element.
      nthelement(values[0], location, loops);
      nthelement(values[1], location, loops);
      #endif
      
      real_t realv = values[0][location];
      real_t imagv = values[1][location];
      
      ans[i].key = hits[i];
      ans[i].value = realv + I*imagv;
    }
    
    free(values[0]);
    free(values[1]);
  }
}


/*
  Outer loop of the algorithm.

  If we are performing the Comb heuristic, first we do so.

  Then, `loops` times:
    choose a random permutation
    run inner_loop_locate
    if in the first location_loops loops, also run inner_loop_filter

  at the end, `hits` contains the coordinates that appear at least
  loop_threshold of location_loops times.  We estimate the values at
  these coordinates as the median of the images x_samp[loops].

  Returns a map from coordinates to estimates.
 */
Node *outer_loop(complex_t *origx, int n, Filter *filter, int B2,
  	   int num, int B, int W_Comb, int Comb_loops, int loop_threshold, int location_loops,
	   int loops, int *hits_found, int *score, int *hits){
  
  
  if(n%B){
    fprintf(stderr, "Warning: n is not divisible by B, which algorithm expects.\n");
    exit(-1);
  }
  
  //Variables used for timing
  double PF_ALL =0;
  double G_ALL =0;
  double B_ALL = 0;
  double C_ALL =0;
  float DtoH=0;
  float HtoD=0;


  //Prepare for location+estimation loops
  int permutea[loops];
  int permute[loops];
  int permuteb[loops];
   
  complex_t *x_sampt= (complex_t *) malloc(B*loops*sizeof(complex_t));
  real_t *samples[loops];
  int *J[loops];

  for(int i = 0; i < loops; i++){
    samples[i] = (real_t*)malloc(B* sizeof(*samples[i]));
    J[i] = (int *)malloc(num * sizeof(*J));
  } 
  
  //Prepare for Comb Filt
  complex_t *x_sampt_comb[Comb_loops];
  real_t *samples_comb[Comb_loops];
  unsigned int offset[Comb_loops]; 
  int sigma = n/W_Comb;
  int *Comb_Approved=NULL;
  if(WITH_COMB)
  {
  
    if(n%W_Comb) {
      fprintf(stderr, "Warning: W_Comb is not divisible by N, which algorithm expects.\n");
      assert(n%W_Comb == 0);
    }

    for(int i=0; i<Comb_loops; i++){
      x_sampt_comb[i] = (complex_t*) malloc(W_Comb*sizeof(*x_sampt_comb[i]));
      samples_comb[i] = (real_t*)malloc(W_Comb * sizeof(*samples_comb[i]));
      offset[i] = (unsigned) floor(drand48() * sigma);
    }
    
    Comb_Approved = (int*) malloc(Comb_loops*num*sizeof(*Comb_Approved));
  }
  
  *hits_found = 0;
  
  //Create FFTW plan
  fftw_plan p_comb = fftw_plan_dft_1d(W_Comb, x_sampt_comb[0], x_sampt_comb[0], FFTW_FORWARD, FFTW_ESTIMATE);
  
  //Create cuFFT plan
  unsigned int plan;
  cufft_plan_create(&plan, B, loops);
  
  //Prepare for the stride factor  
  for(int i=0; i<loops; i++){
    int a=0;
    permuteb[i] = 0;
    while(gcd(a,n) != 1){
       a = (int)(random() & (n-1));
    }
    permutea[i] = a;
    permute[i] = mod_inverse(a, n);
  }

  reset_timer();

  //BEGIN Comb
  double DDD = get_time();
  
  int num_Comb = num;
  if(WITH_COMB){
    #pragma omp parallel for
    for(int i=0; i< Comb_loops; i++){
      Comb_Filt(origx, n, num, W_Comb, Comb_Approved+i*num, x_sampt_comb[i], samples_comb[i], offset[i], sigma, p_comb);
    }
    
    if(Comb_loops > 1){
      radix_sort(Comb_Approved, Comb_loops * num);
      int Last =0;
      for(int i = 1; i < Comb_loops * num; i++){
        if(Comb_Approved[i]!=Comb_Approved[Last])
           Comb_Approved[++Last]=Comb_Approved[i];
      }
      num_Comb=Last+1;
      vprintf("Comb:%d----->%d\n\n",num*Comb_loops,num_Comb);
    }
  } 
  

  if(!ALGORITHM1){
    *hits_found = num_Comb * (n/W_Comb);
    for(int j=0; j < n/W_Comb; j++)
      for (int i=0; i<num_Comb; i++)
        hits[j*num_Comb + i] = j*W_Comb + Comb_Approved[i];
  }

  double Comb_time = get_time() - DDD;
  Pair *permuted_approved = (Pair *)malloc(num_Comb * sizeof(*permuted_approved));
  //END Comb
  
  
  //BEGIN INNER LOOPS
  
  //Step A+B permutation+filter+cufft
  inner_loop_step_a_plus_b(origx, filter, x_sampt, permute, n,B,loops, &PF_ALL, &B_ALL, &DtoH, &HtoD, plan);

  
  //Step C -- Find largest indices
  DDD=get_time();
  #pragma omp parallel for
  for(int i=0; i< loops; i++)
  {
      inner_loop_step_c(num, B, J[i], x_sampt+i*B,samples[i]);
  } 
  C_ALL = get_time() - DDD;
  
  
  //Step D -- Grouping 
  DDD=get_time();
  if(!WITH_COMB){
    for (int i=0; i<location_loops; i++){
      inner_loop_filter_regular(J[i], n, num, B,
                                permutea[i], permute[i], permuteb[i], loop_threshold,
                                score, hits, hits_found);
    }
  }
  else{
    for (int i=0; i<location_loops; i++){
      inner_loop_filter_Comb(J[i], n, num, B,
                           permutea[i], permute[i], permuteb[i], loop_threshold,
                           score, hits, hits_found,
                           Comb_Approved, num_Comb, W_Comb, permuted_approved);
    }
  }
  
  G_ALL=get_time() - DDD;
  //END INNER LOOPS
 


  //BEGIN ESTIMATATION LOOPS
  DDD = get_time();
  Node *ans = (Node *) malloc(*hits_found*sizeof(*ans));

  estimate_values(hits, hits_found, x_sampt,  loops, n, permute, B, filter, location_loops, ans);

  //END ESTIMATION
  double E_T = get_time() - DDD;
  DDD = get_time();
  
  
  //Free up allocated arrays
  for(int i = 0; i < loops; i++)
  {
    free(samples[i]);
    free(J[i]);
  }
  free(x_sampt);
  if(WITH_COMB)
  {
    for (int i=0; i<Comb_loops; i++)
    {
      free(x_sampt_comb[i]);
      free(samples_comb[i]);
    }
    free(Comb_Approved);
    fftw_destroy_plan(p_comb);
  }
  
  free(permuted_approved);
  


  if(TIMING){
	  printf("Total sFFT time: %lf\n", DDD); 
	  printf("Time distribution:   comb    perm+filter  group    estimate   stepB    stepC    HtoD     DtoH     other    total\n");
	  printf("                     %6.4f  %6.4f       %6.4f   %6.4f     %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f\n",
			 Comb_time, PF_ALL, G_ALL, E_T, B_ALL, C_ALL, HtoD, DtoH , DDD-PF_ALL-G_ALL-E_T-Comb_time-B_ALL-C_ALL-HtoD-DtoH, DDD);
	  double tott = (DDD)/100;
	  printf("                     %4.1lf%%   %4.1lf%%        %4.1lf%%    %4.1lf%%      %4.1lf%%    %4.1lf%%     %4.1lf%%   %4.1lf%%     %4.1lf%%    %5.1lf%%\n",
			 Comb_time/tott, PF_ALL/tott, G_ALL/tott, E_T/tott, B_ALL/tott, C_ALL/tott, HtoD/tott, DtoH/tott, (DDD-PF_ALL-G_ALL-E_T-Comb_time-B_ALL-C_ALL-DtoH-HtoD)/tott, (DDD)/tott);

	  printf("\n");
  }

  return ans;

}  
