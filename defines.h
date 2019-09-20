
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <fitsio.h>
//#include <fftw3.h> //siempre detras de complex.h!
//#include <math.h>
//#include <stdio.h>


#ifndef DEFINES_H_
#define DEFINES_H_

//---------------------------------------------------------
//---------------------------------------------------------
//---------------------------------------------------------
//---------------------------------------------------------
// USER CONFIGURATION

//#define CENTRAL_WL   6173.341000 //6173.341000 (oficial requ.) //6173.3500// 341000  // // 6173.335600 //6173.335400

//#include "SVD_configuration.h"
// #include "SVD_configuration_svdcordic18_norm_Nolimit.h"
// #include "SVD_configuration_svdcordic27.h"
//#include "SVD_configuration_svdcordic36.h"
 #include "SVD_configuration_svdcmp.h"
// #include "SVD_configuration_svdcordic18.h"
// #include "SVD_configuration_svdcordic9.h"

// #include "SVD_configuration_svdcordic18_norm_limit28.h"



//INITIAL MODEL A
/*#define INITIAL_MODEL_B 1200
#define INITIAL_MODEL_GM 170
#define INITIAL_MODEL_AZI 25
#define INITIAL_MODEL_ETHA0 14
#define INITIAL_MODEL_LAMBDADOPP 0.07  //en A
#define INITIAL_MODEL_AA 0.05
#define INITIAL_MODEL_VLOS 0.05 // Km/s
#define INITIAL_MODEL_S0 0.25
#define INITIAL_MODEL_S1 0.75*/

//INITIAL MODEL B
// #define INITIAL_MODEL_B 1200
// #define INITIAL_MODEL_GM 170
// #define INITIAL_MODEL_AZI 20
// #define INITIAL_MODEL_ETHA0 12
// #define INITIAL_MODEL_LAMBDADOPP 0.03  //en A
// #define INITIAL_MODEL_AA 1.2
// #define INITIAL_MODEL_VLOS 0.05 // Km/s
// #define INITIAL_MODEL_S0 0.35
// #define INITIAL_MODEL_S1 0.85


//NumeroS cuanticos
#define CUANTIC_NWL 1
#define CUANTIC_SLOI 2
#define CUANTIC_LLOI 1
#define CUANTIC_JLOI 1
#define CUANTIC_SUPI 2
#define CUANTIC_LUPI 2
#define CUANTIC_JUPI 0


#define NOISE_SIGMA 0.001 

#define CLASSICAL_ESTIMATES_SAMPLE_REF 4 //Muestra referencia para cambio de cuadrante de azimuth. Depende del numero de muestras y posicion Continuo


#define NTERMS 10  //ojo si es mayor q 10 casca el svdCordic (esta version)

//#define	PRECISION double //double or float

// END USER CONFIGURATION
//---------------------------------------------------------
//---------------------------------------------------------
//---------------------------------------------------------
//---------------------------------------------------------

// DONT'T MODIFY ANYTHING BELOW OF THIS LINE

#define PI 	3.14159265358979323846264338327950288419716939937510
 		 	
#define ILAMBDA 10
#define TOPLIM 0.000000000001
#define SLIGHT 0
#define NOISE 1e-10 //0.001

#define RR  0.5641895836

#define VLIGHT 2.99792458e+5 //;light speed (km/s); 

#define CTE4_6_13 4.6686411e-13
#define AH 1.0 //angulo heliocentrico

#define FFT_FORWARD -1 
#define FFT_BACKWARD +1

#define NPARMS 4 //(IQUV)

#define LONG_PUNTERO_CALCULOS_COMPARTIDOS 100  //no se usa...

#define INSTRUMENTAL_CONVOLUTION_INTERPOLACION 0  //realizar interpolacion en la convolucion ?? //No funciona !

//INIT_MODEL=[eta0,magnet,vlos,landadopp,aa,gamma,azi,B1,B2,macro,alfa]
struct INIT_MODEL{
	double eta0; // 0
	double B;//magnetic field    
	double vlos;
	double dopp;
	double aa;
	double gm; //5
	double az;
	double S0;
	double S1;
	double mac; //9
	double alfa;		
};

struct CUANTIC{  
	
	double N_PI;
	double N_SIG;
	double * NUB;//size stored in  n_sig
	double * NUP;//size stored in n_pi
	double * NUR;//size stored in n_sig
	double * WEB;//size stored in n_sig
	double * WEP;//size stored in n_pi
	double * WER;//size stored in n_sig
	double GL;
	double GU;
	double GEFF;
	double FO;	
	
};

typedef struct INIT_MODEL Init_Model;
typedef struct CUANTIC Cuantic;

/******************************************************/

void InitializePointerShareCalculation();
void ResetPointerShareCalculation();
void ReadPointerShareCalculation(int Numero,PRECISION ** a,...);
void AsignPointerShareCalculation(int Numero,PRECISION * a,...);
void DeleteSpectraCalculation();

void AllocateMemoryDerivedSynthesis(int numl);
void FreeMemoryDerivedSynthesis();

/******************************************************/


Cuantic * create_cuantic(double * dat);

int me_der(Cuantic *cuantic,Init_Model *initModel,double * wlines,int nwlines,double *lambda,int nlambda,
			PRECISION *d_spectraOut,double ah,double slight,int triplete,int filter);

int mil_sinrf(Cuantic *cuantic,Init_Model *initModel,double * wlines,int nwlines,double *lambda,int nlambda,PRECISION *spectra,
			double ah,int triplete,int filter);
			

double * fgauss(double MC, double * eje,int neje,double landa,int deriv);

double _Complex * fft_d(double * spectra, int nspectra,int direc);
double _Complex * fft_c(double _Complex * spectra, int nspectra,int direc);


int Guarda(char * nombre,PRECISION *v,int nv);
int GuardaC(char * nombre,double _Complex *v,int nv,int a);

int fvoigt(PRECISION damp,PRECISION *vv,int nvv,PRECISION *h, PRECISION *f);

void direct_convolution(PRECISION * x, int nx,PRECISION * h, int nh,PRECISION delta);
PRECISION * vgauss(PRECISION fwhm,int nmuestras_G,PRECISION delta);

/******************* DEFINITIONS FOR READ FITS FILE *********************/

/* 
	Every fits_image will be store in memory like  4 Vector of PRECISION TYPE using this dimension: rows*cols*nLambdas. 
	If we want access a pixel in the vector: stockesI[ (numRow*numCol) + (rows*cols*nLamba)] (been lambda from 0 to nLambda-1)
*/



#define CTYPE1 "CTYPE1"
#define CTYPE2 "CTYPE2"
#define CTYPE3 "CTYPE3"
#define CTYPE4 "CTYPE4"
#define CUNIT1 "CUNIT1"
#define CUNIT2 "CUNIT2"
#define CUNIT3 "CUNIT3"

#define CUNIT_ANSTROM "Angstrom"
#define CUNIT_ARCSEC "arcsec"
#define CTYPE_WAVE "WAVE-GRI"
#define CTYPE_HPLN_TAN "HPLN-TAN"
#define CTYPE_HPLT_TAN "HPLT-TAN"
#define CTYPE_STOKES "STOKES"


/* This Sttructure will store information relative a pixel in the image to process by fuction lmils in C, MPI and CUDA. 
	In order to do more flexible the structure and pass this structure in the message to MPI and KERNEL of CUDA, we have decided 
	put an attribute for number of pixels to store in spectro. We will play with this parameter to pass more or less pixels through MPI and to the 
	KERNEL of CUDA. 
 */
struct VPIXEL {
	PRECISION * vLambda;
	PRECISION * spectro;
	int nLambda;
};

typedef struct VPIXEL vpixels;



struct FITS_IMAGE{  

	int rows;  // number of rows in the image
	int cols;  // number of cols in the image
	int nLambdas; // number of lambdas in the image 
	int numStokes; // number of stokes paramters, normally is 4

	/* POSITION OF EACH DATA IN THE DIMENSIONS OF FITS */
	int pos_lambda;
	int pos_row;
	int pos_col;
	int pos_stokes_parameters; 

	char ctype_1[FLEN_CARD];
	char ctype_2[FLEN_CARD];
	char ctype_3[FLEN_CARD];
	char ctype_4[FLEN_CARD];
	char cunit_1[FLEN_CARD];
	char cunit_2[FLEN_CARD];
	char cunit_3[FLEN_CARD];
	char cunit_4[FLEN_CARD];

	int numPixels;
	vpixels * pixels;
	PRECISION * vLambdaImagen;
	PRECISION * spectroImagen;
};

typedef struct FITS_IMAGE FitsImage;

#define NUMBER_PARAM_MODELS 11

#define PATH_MAX 4096



#endif /*DEFINES_H_*/
