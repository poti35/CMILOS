#include <complex.h>
#include <fftw3.h> //siempre a continuacion de complex.h
#include <math.h>
#include <stdio.h>


#define FFT_FORWARD -1
#define FFT_BACKWARD +1

double _Complex *fft_d(double *spectra, int nspectra, int direc);
double _Complex *fft_c(double _Complex *spectra, int nspectra, int direc);
void FFTConvolve( double * data,  double * kernel, int n, int h);
void FFTConvolve1D( double * data, int sizeData, int h, fftw_complex * kernel, fftw_complex * in_data, fftw_complex * freq_data, fftw_complex * rev_data, fftw_complex * time_data, fftw_plan p2, fftw_plan rev);