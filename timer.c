#include<time.h>
#include "timer.h"
#include <sys/time.h>

/*
timespec start_time;

void reset_timer(){
  clock_gettime(TIMER_TYPE, &start_time);
}

double get_time(){
  timespec t;
  clock_gettime(TIMER_TYPE, &t);
  return double(t.tv_sec - start_time.tv_sec) + double(t.tv_nsec - start_time.tv_nsec) * 1.e-9;
}

*/

//New timer implementation using gettimeofday() Cheng 03/11/13
//
struct timeval start;

void reset_timer(){
  gettimeofday(&start, NULL);
}

double get_time(){
  struct timeval end;
  gettimeofday(&end, NULL);
  return ((double)(end.tv_sec-start.tv_sec) + (double)(end.tv_usec/1000000.0 - start.tv_usec/1000000.0));
}
