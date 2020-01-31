#include <math.h>
#include <stdio.h>
#include <math.h>
#include <stdio.h>
#include "convolution.h"
//#include "mkl_vsl.h"

/*

	@autor: Juan Pedro Cobos
	@Date: 31 Marzo 2011
	@loc: IAA- CSIC
	
	Convolucion para el caso Sophi: convolucion central de x con h.
	
	direct_convolution(x,h,delta)
	x spectro
	h gaussiana -perfil instrumental- � ojo, solo con longitud impar!
	delta anchura de muestreo

	--nota: como h es simetrico no se invierte su orden 
	
	//result=calloc(nresult,sizeof(PRECISION*));
	//free();

	_juanp
*/
extern PRECISION *dirConvPar;


void direct_convolution_double(PRECISION *x, int nx, PRECISION *h, int nh)
{

	
	int nx_aux;
	int k, j;

	nx_aux = nx + nh - 1; // tamano de toda la convolucion
	

	int mitad_nh = nh / 2;
   
	// rellenamos el vector auxiliar
	for (k = 0; k < nx_aux; k++)
	{
		dirConvPar[k] = 0;
	}

	for (k = 0; k < nx; k++)
	{
		dirConvPar[k + mitad_nh] = x[k];
	}

	// vamos a tomar solo la convolucion central
	
	
	for (k = 0; k < nx; k++)
	{
		//x[k] = 0;
		double aux = 0;	
		for (j = 0; j < nh; j++)
		{
			aux += h[j] * dirConvPar[j + k];
			
		}
		x[k] = aux;
	}

}

void direct_convolution(REAL *x, int nx, PRECISION *h, int nh)
{

	
	int nx_aux;
	int k, j;

	nx_aux = nx + nh - 1; // tamano de toda la convolucion
	

	int mitad_nh = nh / 2;
   
	// rellenamos el vector auxiliar
	for (k = 0; k < nx_aux; k++)
	{
		dirConvPar[k] = 0;
	}

	for (k = 0; k < nx; k++)
	{
		dirConvPar[k + mitad_nh] = x[k];
	}

	// vamos a tomar solo la convolucion central
	
	
	for (k = 0; k < nx; k++)
	{
		//x[k] = 0;
		double aux = 0;	
		for (j = 0; j < nh; j++)
		{
			aux += h[j] * dirConvPar[j + k];
			
		}
		x[k] = aux;
	}

}


void direct_convolution2(REAL *x, int nx, PRECISION *h, int nh,REAL * result,int delta)
{

	
	int nx_aux;
	int k, j;

	nx_aux = nx + nh - 1; // tamano de toda la convolucion
	

	int mitad_nh = nh / 2;
   
	// rellenamos el vector auxiliar
	for (k = 0; k < nx_aux; k++)
	{
		dirConvPar[k] = 0;
	}

	for (k = 0; k < nx; k++)
	{
		dirConvPar[k + mitad_nh] = x[k];
	}

	// vamos a tomar solo la convolucion central
	
	for (k = 0; k < nx; k++)
	{
		
		double aux = 0;
		for (j = 0; j < nh; j++)
		{
			aux += h[j] * dirConvPar[j + k];
		}
		result[k] = aux*delta;
	}

}

void convolve(REAL * Signal, size_t SignalLen, double * Kernel, size_t KernelLen, REAL * Result , int delta)
{
  size_t n;

  for (n = 0; n < SignalLen + KernelLen - 1; n++)
  {
    size_t kmin, kmax, k;

    dirConvPar[n] = 0;
	double aux = 0;
    kmin = (n >= KernelLen - 1) ? n - (KernelLen - 1) : 0;
    kmax = (n < SignalLen - 1) ? n : SignalLen - 1;

    for (k = kmin; k <= kmax; k++)
    {
      aux += Signal[k] * Kernel[n - k];
    }
	dirConvPar[n] = aux;
  }

  int mitad_nh = SignalLen / 2;
  for(n=0;n<SignalLen;n++){
	  Result[n] = dirConvPar[n+mitad_nh];
  }
}


