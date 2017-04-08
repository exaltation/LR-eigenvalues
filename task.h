#ifndef task_h
#define task_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

extern int _DEBUG, _ERROR;

int sim(int n, double* A, double* tmp, double precision);
int evc(int n, int max_iterations, double epsilon, double* A, double* E, double* tmp, double precision);
size_t sim_memsize(int n);
size_t evc_memsize(int n);
int compare (const void *a, const void *b);

#endif
