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

void direct_convolution2(REAL *x, int nx, PRECISION *h, int nh, REAL *result, int delta)
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
		result[k] = aux * delta;
	}
}

void convolve(REAL *Signal, size_t SignalLen, double *Kernel, size_t KernelLen, REAL *Result, int delta)
{
	size_t n, i;

	/*for (n = 0; n < SignalLen + KernelLen - 1; n++)
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
  */
	REAL *signalAux = calloc(SignalLen + SignalLen - 1, sizeof(REAL));
	double *ketnelAux = calloc(KernelLen + KernelLen - 1, sizeof(double));
	int mitad_nh = SignalLen / 2;
	for (n = mitad_nh; n < SignalLen + mitad_nh; n++)
	{
		signalAux[n] = Signal[n];
		ketnelAux[n] = Kernel[n];
	}
	for (n = 0; n < SignalLen + KernelLen - 1; n++)
	{
		dirConvPar[n] = 0;
		for (i = 0; i <= n; i++)
		{
			dirConvPar[n] = dirConvPar[n] + (signalAux[i] * ketnelAux[n - i]);
		}
	}

	for (n = 0; n < SignalLen; n++)
	{
		Result[n] = dirConvPar[n + mitad_nh];
	}
}

void convCircular(REAL *x, int m, double *h, int n, REAL *result)
{
	int i,j, k,ishift;
	double aux;
	double x2[m];
	double a[m];
	int odd=(m%2);		
	int startShift = m/2;
	if(odd) startShift+=1;	
	ishift = startShift;

	//result[0]=0;
    a[0]=h[0];

    for(j=1;j<n;j++)            /*folding h(n) to h(-n)*/
    	a[j]=h[n-j];

    /*Circular convolution*/
	aux = 0;
	for(i=0;i<n;i++)
        aux+=x[i]*a[i];
	result[ishift++] = aux;

	for(k=1;k<n;k++)
    {
    	aux=0;
        /*circular shift*/

        for(j=1;j<n;j++)
        	x2[j]=a[j-1];
        x2[0]=a[n-1];
        for(i=0;i<n;i++)
        {
        	a[i]=x2[i];
            aux+=x[i]*x2[i];
        }
		if(k <m/2)
			result[ishift++] = aux;
		else
			result[k-(n/2)] = aux;
		
    }

	/*for(i=0,ishift=startShift;i<numl/2;i++,ishift++){
		d_spectra[ishift+9*numl]=results[i];
	}
	for(i=(numl/2),ishift=0;i<numl;i++,ishift++){
		d_spectra[ishift+9*numl]=results[i];
	}	*/

}
