#include "defines.h"

void direct_convolution_double(PRECISION *x, int nx, PRECISION *h, int nh);
void direct_convolution(REAL *x, int nx, PRECISION *h, int nh);
void direct_convolution2(REAL *x, int nx, PRECISION *h, int nh,REAL * result,int delta);
void convolve(REAL * Signal, size_t SignalLen, REAL * Kernel, size_t KernelLen);
//int convolutionMKL(double  * h, double * x,  double * y, int  start,VSLConvTaskPtr * task);