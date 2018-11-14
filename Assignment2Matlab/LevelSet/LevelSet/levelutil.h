#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "mex.h"

#define TINY 1.0e-20

/* prototype of functions */

extern double detinv(double **a, int n);
extern void lubksb(double **a, int n, int *indx, double b[]);
extern void ludcmp(double **a, int n, int *indx, double *d);
void add_point(double *I, double *param, int N, int sign);
double cal_uu(double *I, double *param, int N);
void update_param(double *param,int N);
