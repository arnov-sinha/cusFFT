/*
 *The C version of the sparse FFT algorithm.
 *File: experiment.c
 *Include the main function and data initialization
 *
 *Copyright (C) and license info
 * 
*/


#include <getopt.h>
#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <complex.h>
#include "fft.h"
#include "filters.h"
#include "computefourier.h"
#include "timer.h"
#include "utils.h"
#include <omp.h>
#include <time.h>
#include "cuda_fft.h"
double evaluate_runtime(int n, double lobefrac, double tolerance,
                        double lobefrac2, double tolerance2,
                        int num, int B, int B2, int location_loops,
                        int est_loops, int loop_threshold, int W_Comb,
                        int Comb_loops);

/*
  Run an experiment.  Parameters are:

  x: the signal
  n: the length of x

  lobefrac_loc:   during location, the main lobe of the filter has half
                  width n*lobefrac_loc
  tolerance_loc:  the linf norm of the residuals of the filter
  b_loc:          the number of adjacent filters to add

  B_loc:          number of samples in subsampling during location
  B_thresh:       number of samples considered "heavy"
  loops_loc:      number of location loops
  loops_thresh:   number of times a coordinate must seem heavy.

  *_est:          ditto as above, but for estimation.

  repetitions:    repeat the experiment this many times for timing
  LARGE_FREQ:     locations of largest coefficients.
  k:              number of HHs, used for evaluation only.
  x_f:            true fft.
 */
  

void run_experiment(complex_t *x, int n,
		    double lobefrac_loc, double tolerance_loc, int b_loc,
		    int B_loc, int B_thresh, int loops_loc, int loops_thresh,
		    double lobefrac_est, double tolerance_est, int b_est,
		    int B_est, int loops_est, int W_Comb, int Comb_loops,
		    int repetitions, bool_t FFTW_OPT,
		    int *LARGE_FREQ, int k, complex_t *x_f, int *hits_found){
  
  printf("sFFT filter parameters for: n=%d, k=%d.\n", n, k);
  printf("******************************************************************************\n");
  
  if (WITH_COMB)
    printf(" Comb Filter: loops: %d mod: %d/%d\n", Comb_loops, B_thresh, W_Comb);
  else
    printf(" Comb Filter: none\n");
  if (ALGORITHM1)
    printf(" Location Filter: (numlobes=%.1lf, tol=%lg, b=%d) B: %d/%d loops: %d/%d\n", 0.5/lobefrac_loc, tolerance_loc, b_loc, B_thresh, B_loc, loops_thresh, loops_loc);
  else
    printf(" Location Filter: none\n");
  
  printf(" Estimation Filter: (numlobes=%.1lf, tol=%lg, b=%d) B: %d loops: %d\n", 0.5/lobefrac_est, tolerance_est, b_est, B_est, loops_est);
  printf("\n");

  if (WITH_COMB)
    assert(B_thresh < W_Comb);
  if (ALGORITHM1) {
    assert(B_thresh < B_loc);
    assert(loops_thresh <= loops_loc);
  }
  
  //Prepare for filter
  int w_loc;
  complex_t *filtert = make_dolphchebyshev_t(lobefrac_loc, tolerance_loc, &w_loc);
  Filter filter = make_multiple_t(filtert, w_loc, n, b_loc);
  
  //Paddled filters
#if 1 
  Filter pad_filter, pad_filter_est;
  complex_t *pad_filter_t=NULL, *pad_filter_est_t=NULL, *pad_filter_trans=NULL;

  if(filter.sizet % B_loc){
    int pad_size = (filter.sizet/B_loc+1)*B_loc;
    
    pad_filter_t = (complex_t *)malloc(pad_size*sizeof(*pad_filter_t));

    memset(pad_filter_t+filter.sizet, 0, (pad_size-filter.sizet)*sizeof(*pad_filter_t));
    
    pad_filter.time = pad_filter_t;
    
    memcpy(pad_filter.time, filter.time, filter.sizet * sizeof(*pad_filter_t));
    
    pad_filter.sizet = pad_size;
    pad_filter.freq = filter.freq;
  }
  else{
    pad_filter.time = filter.time;
    pad_filter.sizet = filter.sizet;
    pad_filter.freq = filter.freq;
  }

  //Transpose the filter
#if 1 
  pad_filter_trans = (complex_t*)malloc(pad_filter.sizet * sizeof(*pad_filter_trans));
  int round  = pad_filter.sizet / B_loc;

  for(int i=0; i<B_loc; i++){
    int off = i*round;
    for(int j=0; j< round; j++){
      pad_filter_trans[off+j] = pad_filter.time[j*B_loc+i];
    }
  }
  pad_filter.time = pad_filter_trans;
#endif
#endif

  int w_est;
  complex_t *filtert_est = make_dolphchebyshev_t(lobefrac_est, tolerance_est, &w_est);
  Filter filter_est = make_multiple_t(filtert_est, w_est, n, b_est);

#if 1 
  //Estimation filter
  if(filter_est.sizet % B_est){
    int pad_size_est = (filter_est.sizet/B_est+1)*B_est;
    
    pad_filter_est_t = (complex_t *)malloc(pad_size_est*sizeof(*pad_filter_est_t));
    memset(pad_filter_est_t+filter_est.sizet, 0, (pad_size_est-filter_est.sizet)*sizeof(*pad_filter_est_t));
    pad_filter_est.time = pad_filter_est_t;
    
    memcpy(pad_filter_est.time, filter_est.time, filter_est.sizet * sizeof(*pad_filter_est_t));
    
    pad_filter_est.sizet = pad_size_est;
    pad_filter_est.freq = filter_est.freq;
  }
  else{
    pad_filter_est.time = filter_est.time;
    pad_filter_est.sizet = filter_est.sizet;
    pad_filter_est.freq = filter_est.freq;
  }
#endif
  printf(" Window size: Location Filter : %d; Estimation Filter : %d;\n", w_loc, w_est);

  real_t filter_noise = 0, filter_noise_est = 0;

  for (int i=0; i<10; i++){
    real_t tmp1 = cabs(pad_filter.freq[n/2+i]) > cabs(pad_filter.freq[n/2-i])? cabs(pad_filter.freq[n/2+i]):cabs(pad_filter.freq[n/2-i]);
    filter_noise = tmp1>filter_noise?tmp1:filter_noise;

    real_t tmp2= cabs(pad_filter_est.freq[n/2+i]) > cabs(pad_filter_est.freq[n/2-i])? cabs(pad_filter_est.freq[n/2+i]):cabs(pad_filter_est.freq[n/2-i]);
    filter_noise_est = tmp2>filter_noise_est? tmp2: filter_noise_est;
  }
  
  printf(" Noise in filter: Location Filter : %lg; Estimation Filter %lg\n", filter_noise, filter_noise_est);
  printf("******************************************************************************\n\n");
 
  //Prepare for array
  // Create the score array
  // calloc is faster if few pages are hit, while malloc/memset is
  // faster if most pages are hit.
  double pages_hit = B_thresh * (n/B_loc) * (WITH_COMB ? B_thresh * 1. / W_Comb : 1.) * loops_loc;
  int *score;
  if (pages_hit > n / 1024) {
    score =(int*)malloc(n*sizeof(*score));
    memset(score, 0, n*sizeof(*score));
  } else {
    score =(int*)calloc(n,sizeof(*score));
  }
  
  int *hits = (int*)malloc(n*sizeof(*hits));

  //Compute the N-dimensional cuFFT
  float cu_fft_time = 0;

  //cuda_fft(x,x_f,n,repetitions, &cu_fft_time);
  
  //Start sFFT computation
  printf("sFFT Results\n");
  printf("******************************************************************************\n");
  
  //struct array to store the output of the outer_loop
  Node *ans=NULL;
  for(int i = 0; i < repetitions; i++){
	  reset_timer();
	  ans = outer_loop(x, n, &pad_filter, B_est, B_thresh, B_loc, W_Comb,
		    Comb_loops, loops_thresh, loops_loc, loops_loc + loops_est, hits_found,score, hits);
  }
 
  //End of the sFFT computation
  
  printf("hits_found = %d\n", *hits_found);
  
  //We only cutoff the k-largest values in the frequency domain 
  //Analyze the results of sFFT

  complex_t *ans_Large = (complex_t *)calloc(n,sizeof(*ans_Large));
  
  
  //sort the ans according to the cabs value of second element
  qsort(ans, *hits_found, sizeof(*ans), comp_struct);
  for(int i=0; i<k; i++){
    ans_Large[ans[i].key] = ans[i].value;
  }

  int large_found = 0;
  int FOUND=0;
  qsort(ans, *hits_found, sizeof(*ans), comp_struct2);
  
  for(int i = 0; i < k; i++){
    Node *item = (Node*)bsearch(&LARGE_FREQ[i], ans, *hits_found, sizeof(*ans), comp_struct2);
    if(item != NULL ){
      FOUND++;
    }
    large_found += (ans_Large[(LARGE_FREQ[i])] != 0);
  }
  
  //Estimate error as the difference of the K largest entries of x_f and ans)
  real_t ERROR =0;
  for(int i=0; i< n ; i++){
      ERROR += cabs(ans_Large[i]- x_f[i]);
  }

   //printf("---------------CONSIDER K LARGEST ONLY ---------------------------------\n");
   printf("ERROR:\n");
   printf("K=%d; MISSED (estimation, result) = (%d, %d); L1 ERROR= %lg  (%lg per large frequency)\n",k, k-FOUND, k-large_found, ERROR, ERROR/k);
   printf("******************************************************************************\n\n");


   printf("CUFFT Results\n");
   printf("******************************************************************************\n");
   
   reset_timer();
   cuda_fft(x,x_f,n,repetitions, &cu_fft_time);
   //printf("Time to run cuFFT : %f\n", cu_fft_time/1e3);
   printf("******************************************************************************\n\n");


  free(ans_Large);
  free(filter.freq);
  free(filter.time);
  free(filter_est.freq);
  free(filter_est.time);
  free(ans);
  free(pad_filter_t);
  
  free(pad_filter_est_t);
  free(pad_filter_trans);
  free(score);  
  free(hits);
}



static void
usage(const char *progname)
{
  const char *p = strrchr(progname, '/');       // drop leading directory path
  if (p)
    p++;
  if (strncmp(p, "lt-", 3) == 0)                // drop lt- libtool prefix
    p += 3;

  fprintf(stderr, "Usage: %s [options]\n\n", p);
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "  -h                   show this message and exit\n");
  fprintf(stderr, "  -N Number of Frequencies [default= 2^22]\n");
  fprintf(stderr, "  -K NUmber of Energetic Frequencies [default=100]\n");
  fprintf(stderr, "  -R Number of Runs [default=1]\n");
  fprintf(stderr, "  -B Constant for number of Buckets [default=2]\n");
  fprintf(stderr, "  -E Constant for number of Buckets for Estimation [default=0.2]\n");
  fprintf(stderr, "  -L Number of estimation loops [default=12]\n");
  fprintf(stderr, "  -l Number of location loops [default=3]\n");
  fprintf(stderr, "  -r Number of required location loops [default=2]\n");
  fprintf(stderr, "  -M Use SFFT2.0 [default=SFFT1.0] and set Constant for number of Comb Buckets [default=12]\n");
  fprintf(stderr, "  -m number of times we run Comb [default=1]\n");
  fprintf(stderr, "  -S SNR as ratio [default=inf]\n");
  fprintf(stderr, "  -A USE Alternative algorithm\n");
  fprintf(stderr, "  -O Use FFTW after optimization\n");
  fprintf(stderr, "  -t Tolerance for the Location Filters [default=1.e-6]\n");
  fprintf(stderr, "  -e Tolerance for the Estimation Filters [default=1.e-8]\n");
  fprintf(stderr, "  -s Calculate the expected runtime and error of sFFT and return.\n");
  fprintf(stderr, "  -v Verbose : prints detailed timing of each step of the code.\n");
}


int main(int argc, char **argv){
  int n = 4*128*8192;
  int k = 100;
  //  int w = 3*(int)(sqrt(k*n));
  int repetitions = 1;
  double Bcst_loc=2;
  double Bcst_est=0.2;
  double Comb_cst=16;
  int loc_loops =3;
  int est_loops =12;
  int threshold_loops =2;
  int Comb_loops = 1;
  int ch;
  int simulate = 0;
  double snr=1000000000;
  double std_noise = 0;
  bool_t FFTW_OPT = false;
  double tolerance_loc = 1.e-6;
  double tolerance_est = 1.e-8;  
  
  while ((ch = getopt(argc, argv, "shvN:K:R:B:E:L:l:r:M:m:S:AOt:e:")) != EOF){
    switch (ch){

    case 'N':
      n =atoi(optarg);
      break;

    case 'K':
      k =atoi(optarg);
      break;

    case 'R':
      repetitions = atoi(optarg);
      break;
    case 'B':
      Bcst_loc = strtod(optarg, 0);
      break;
    case 'E':
      Bcst_est = strtod(optarg, 0);
      break;
  case 'L':
      est_loops = atoi(optarg);
      break;
  case 'l':
      loc_loops =atoi(optarg);
      if (!loc_loops)
        ALGORITHM1 = false;
      break;
  case 'r':
      threshold_loops =atoi(optarg);
      break;
  case 'M':
      Comb_cst =strtod(optarg,0);
      WITH_COMB = true;
      break;
  case 'm':
      Comb_loops =atoi(optarg);
      break;
  case 'S':
      snr =strtod(optarg,0);
      std_noise = sqrt(k/(2*snr));
      break;
  case 'A':
      ALGORITHM1 = false;
      loc_loops = 0;
      break;
  case 'O':
      FFTW_OPT = true;
      break;
  case 'v':
      VERBOSE = true;
      break;
  case 't':
      tolerance_loc = strtod(optarg, 0);
      break;
  case 'e':
      tolerance_est = strtod(optarg, 0);
      break;
  case 's':
    simulate = 1;
    break;
   case 'h':
    default:
      usage(argv[0]);
      exit(1);
    }
  }
  
  if(Bcst_loc!=Bcst_est){
    Bcst_est = Bcst_loc;
    printf("WARNING: E is not equal to B, have assigned B to E.\n");
  }
  
  n = floor_to_pow2(n);
  assert(ALGORITHM1 || WITH_COMB);

  real_t BB_loc = (unsigned) (Bcst_loc*sqrt((double)n*k/(log2(n))));
  real_t BB_est = (unsigned) (Bcst_est*sqrt((double)n*k/(log2(n))));


  double lobefrac_loc = 0.5 / (BB_loc);
  double lobefrac_est = 0.5 / (BB_est);

  int b_loc = (int)(1.2*1.1*((double) n/BB_loc));
  //b_loc = 1;
  //real_t BB2 = (unsigned) (Bcst2*sqrt((double)n*k/(log2(n))));
  int b_est = (int)(1.4*1.1*((double) n/BB_est));
  //int b_est = (int)(1.2*1.1*((double) n/BB_est));

  int B_loc = floor_to_pow2(BB_loc);
  int B_thresh = 2*k;
  int B_est = floor_to_pow2(BB_est);


  int W_Comb = floor_to_pow2(Comb_cst*n/B_loc);

  printf("\n\nRUNNING EXPERIMENT: n=%d, k=%d.\n", n, k);
  
  printf("\n\nSimulation:\n");
  printf("******************************************************************************\n");
  printf("Expected running time: %lg\n",
         evaluate_runtime(n, lobefrac_loc, tolerance_loc, lobefrac_est, tolerance_est,
                          B_thresh, B_loc, B_est, loc_loops, est_loops, threshold_loops, W_Comb, Comb_loops) * 1e-9);
  printf("******************************************************************************\n");  
  if (simulate)
    return 0;
  
  complex_t *x=(complex_t *)malloc(n*sizeof(*x));
  
  //Fix the input, will random Later. 
  //srand(17);
  //srand48( time(NULL) ^ (getpid() * 171717));
  //srand(time(NULL));
  //Randomized the None Zero Bins
  complex_t *x_f = (complex_t *)calloc(n, sizeof(*x_f));
  
  int *LARGE_FREQ = (int *)malloc(k*sizeof(*LARGE_FREQ));

  
  //Initialize the large bins in freq domain
  for (int i=0; i<k; i++)
  {
    LARGE_FREQ[i] = (unsigned)floor(drand48() * n);
    x_f[LARGE_FREQ[i]] = 1.0;
  }

  fftw_dft(x,n,x_f,1);    //Reverse FFT to generate input x

  //ADDED NOISE
  double snr_achieved;
  snr_achieved = AWGN(x,n,std_noise);
  if(std_noise != 0)
    printf("SNR = %g / %.2f dB \n\n", snr_achieved, 10 * log10(snr_achieved));  

  fftw_dft(x_f, n, x, 0);

  for(int i = 0; i < n; i++)
    x_f[i] /= n;

  int hits_found = 0;
  run_experiment(x, n,
                 lobefrac_loc, tolerance_loc, b_loc,
                 B_loc, B_thresh, loc_loops, threshold_loops,
                 lobefrac_est, tolerance_est, b_est,
                 B_est, est_loops, W_Comb, Comb_loops,
                 repetitions, FFTW_OPT, LARGE_FREQ, k, x_f, &hits_found);

  
  
  free(x);
  free(x_f);
  free(LARGE_FREQ);
   
  return 0;
}

double evaluate_runtime(int n, double lobefrac, double tolerance,
                        double lobefrac2, double tolerance2,
                        int num, int B, int B2, int location_loops,
                        int est_loops, int loop_threshold, int W_Comb,
                        int Comb_loops){
  int w = (int)((1 / M_PI) * (1/lobefrac) * acosh(1./tolerance));
  int w2 = (int)((1 / M_PI) * (1/lobefrac2) * acosh(1./tolerance2));

  if (WITH_COMB)
    if(num >= W_Comb)
      return -1;
  if (ALGORITHM1) {
    if(num >= B)
      return -1;
    if (loop_threshold > location_loops)
      return -1;
  }
  if (w > n || w2 > n)
    return -1;

  int loops = location_loops + est_loops;
  double m_frac = WITH_COMB ? num * 1. / W_Comb : 1;
  double projected_hits = ALGORITHM1 ? binomial_cdf(num * (1. / B - 1./n), location_loops, loop_threshold) * n * m_frac + num/2 : n * m_frac;
  //XXX B2 for some, B for some
  int k_est = num/2;
  double projected_noise_on_k = 2 * binomial_cdf(k_est * (1. / B2 - 1./n) / 2, loops, (loops+1)/2) * k_est;
  double projected_error_rate = 2 * binomial_cdf(k_est * (1. / B2 - 1./n) / 4, loops, (loops+1)/2) * n * m_frac + projected_noise_on_k;
  //double projected_error_rate = binomial_cdf((num/2) * (1. / B2 - 1./n), est_loops, (est_loops+1)/2) * (projected_hits - num/2);
  printf("Projected error rate: %lg (%lg per large frequency)\n", projected_error_rate, projected_noise_on_k);


  double pages_to_set = num * (n/B) * m_frac * location_loops * 1024;
  bool_t will_array_memset = pages_to_set > n;

  double const_scorearray = n < (1<<21) ? 0.3 : 1.8;
  double const_permfilt = 38.0;
  double const_Combtime = 90.0;
  double const_estimationtime = WITH_COMB ? 140 : 150;
  double const_grouping = 23;
  double const_group_sort = (WITH_COMB ? 30 : 0);
  double const_bplusctime = 41;

  double time_scorearray = will_array_memset ? const_scorearray * n : 0;
  double time_Comb = const_Combtime * W_Comb * Comb_loops;
  double time_permfilt = const_permfilt*(w *1. * location_loops + w2 * 1. * est_loops);
  double time_grouping = (location_loops * (const_grouping * num * (n/B) * m_frac
                                           + const_group_sort * num * log(num))
                          + const_scorearray * (!will_array_memset) * pages_to_set);
  double time_estimation = const_estimationtime * projected_hits * loops;
  double time_bplusc = const_bplusctime * (location_loops * B + est_loops * B2);
  double time_total = time_scorearray + time_Comb + time_permfilt + time_grouping + time_estimation + time_bplusc;

  return time_total;
}

