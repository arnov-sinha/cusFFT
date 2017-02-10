#ifndef FFT_H
#define FFT_H

#include<math.h>
#include<complex.h>
#include<fftw3.h>

#define WRAPPER

// allow easy change to float or long double
//#define USE_FLOAT
#define USE_DOUBLE

#ifdef USE_FLOAT
typedef float complex complex_t;
typedef float real_t;
#define cexp cexpf
#define exp expf
#endif


#ifdef USE_DOUBLE
//typedef double complex complex_t;
typedef fftw_complex complex_t;
typedef double real_t;
#endif

//#define DEBUG

typedef struct 
{
  unsigned int key;
  complex_t value; 
}Node;

typedef struct 
{
  unsigned int key;
  complex_t value;
  int write; 
}Node2;

typedef struct{
  int first;
  int second;
}Pair;



typedef int bool_t;
#define true 1
#define false 0

#define M_PI		3.14159265358979323846	/* pi */

#endif
