#include "defines.h"

/**
 * 
 */
void spectral_synthesis_convolution(int * nlambda, int * INSTRUMENTAL_CONVOLUTION, int * NMUESTRAS_G);

/**
 * 
 */
void response_functions_convolution(int * nlambda, int * INSTRUMENTAL_CONVOLUTION, int * NMUESTRAS_G);

/**
 * 
 */
void AplicaDelta(Init_Model *model, PRECISION *delta, int *fixed, Init_Model *modelout);
/**
 * 
 */
int check(Init_Model *Model);
/**
 * 
 */
void FijaACeroDerivadasNoNecesarias(PRECISION *d_spectra, int *fixed, int nlambda);

/**
 * 
 */
int mil_svd(PRECISION *h, PRECISION *beta, PRECISION *delta);

/**
 * 
 */
void weights_init(int nlambda, double *sigma, PRECISION *weight, int nweight, PRECISION **wOut, PRECISION **sigOut, double noise);

/*
*
*
* Cálculo de las estimaciones clásicas.
*
*
* lambda_0 :  centro de la línea
* lambda :    vector de muestras
* nlambda :   numero de muesras
* spectro :   vector [I,Q,U,V]
* initModel:  Modelo de atmosfera a ser modificado
*
*
*
* @Author: Juan Pedro Cobos Carrascosa (IAA-CSIC)
*		   jpedro@iaa.es
* @Date:  Nov. 2011
*
*/
void estimacionesClasicas(PRECISION lambda_0, double *lambda, int nlambda, PRECISION *spectro, Init_Model *initModel);


/*
 *
 * nwlineas :   numero de lineas espectrales
 * wlines :		lineas spectrales
 * lambda :		wavelength axis in angstrom
			longitud nlambda
 * spectra : IQUV por filas, longitud ny=nlambda
 */
int lm_mils(Cuantic *cuantic, double *wlines, int nwlines, double *lambda, int nlambda, PRECISION *spectro, int nspectro,
				Init_Model *initModel, PRECISION *spectra, double *chisqrf, int *iterOut,
				double slight, double toplim, int miter, PRECISION *weight, int nweight, int *fix,
				PRECISION *sigma, double filter, double ilambda, double noise, double *pol,
				double getshi, int triplete, int * INSTRUMENTAL_CONVOLUTION, int * NMUESTRAS_G);


/**
 * Generate the guaussian from instrumental profile 
 */
void generateGaussianInstrumentalProfile(PRECISION * G, PRECISION  FWHM, PRECISION DELTA, int INSTRUMENTAL_CONVOLUTION_WITH_PSF, int NMUESTRAS_G);
