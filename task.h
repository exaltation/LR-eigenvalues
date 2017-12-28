#ifndef task_h
#define task_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

extern int _DEBUG, _ERROR;

int sim(int n, long double* A, long double* tmp, long double precision);
int evc(int n, int max_iterations, long double epsilon, long double* A, long double* E, long double* tmp, long double precision);
size_t sim_memsize(int n);
size_t evc_memsize(int n);
int compare (const void *a, const void *b);

#endif
