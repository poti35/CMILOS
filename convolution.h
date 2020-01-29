#include "defines.h"

void direct_convolution_double(PRECISION *x, int nx, PRECISION *h, int nh);
void direct_convolution(REAL *x, int nx, PRECISION *h, int nh);
void direct_convolution2(REAL *x, int nx, PRECISION *h, int nh,REAL * result,int delta);