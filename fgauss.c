
#include "defines.h"

/*
 * 
 * deriv : 1 true, 0 false
 */

//;this function builds a gauss function
//;landa(amstrong) ;Central wavelength
//;eje(amstrong) ;Wavelength axis
//;macro ;Macroturbulence in km/s
extern REAL *GMAC;

void fgauss(PRECISION MC, double *eje, int neje, double landa, int deriv)
{
	double centro;
	double ild;
	double term[neje];
	int i;
	double cte;

	centro = eje[(int)neje / 2];		  //center of the axis
	ild = (landa * MC) / 2.99792458e5; //Sigma

	//	printf("ild-> %f  ...\n",ild);


	for (i = 0; i < neje; i++)
	{
		double aux = ((eje[i] - centro) / ild);
		term[i] = ( aux * aux) / 2; //exponent
		//printf("term (%d) %f  ...\n",i,term[i]);
	}

	for (i = 0; i < neje; i++)
	{
		GMAC[i] = exp(-term[i]);
	}

	cte = 0;
	//normalization
	for (i = 0; i < neje; i++)
	{
		cte += GMAC[i];
	}
	for (i = 0; i < neje; i++)
	{
		GMAC[i] /= cte;
	}

	//In case we need the deriv of f gauss /deriv
	if (deriv == 1)
	{
		for (i = 0; i < neje; i++)
		{
			//mtb2=mtb/macro*(((eje-centro)/ILd)^2d0-1d0)
			GMAC[i] = GMAC[i] / MC * ((((eje[i] - centro) / ild) * ((eje[i] - centro) / ild)) - 1.0);			
		}
	}

}



/*
 * 
 * deriv : 1 true, 0 false
 */

//;this function builds a gauss function
//;landa(amstrong) ;Central wavelength
//;eje(amstrong) ;Wavelength axis
//;macro ;Macroturbulence in km/s

REAL * fgauss_WL(REAL FWHM, REAL step_between_lw, REAL lambda0, REAL lambdaCentral, int nLambda, int * sizeG)
{
	//int fgauss(PRECISION MC, PRECISION * eje,int neje,PRECISION landa,int deriv,PRECISION * mtb,int nmtb){

	REAL *mtb;
	REAL *term, *loai;
	int i;
	int nloai, nmtb;
	REAL cte;
	

	//int even = (nLambda%2);
	*sizeG = nLambda;

	/*if(even) {
		*sizeG = nLambda;
	}
	else{
		*sizeG = nLambda + 1;
	}*/

	//printf("\nCENTRO DE CREACIÓN DE LA GAUSSIANA: %d\n",centro);
	//lambda0 = vLambda[centro];		  //center of the axis
	//step_between_lw = vLambda[1]-vLambda[0]; // step to create the gaussian 
	//ild = (landa * MC) / 2.99792458e5; //Sigma

	///Conversion from FWHM to Gaussian sigma (1./(2*sqrt(2*alog2)))
	REAL sigma=FWHM*0.42466090/1000.0; // in Angstroms
	//PRECISION sigma = FWHM * (2 * sqrt(2 * log(2)))/1000;

	//printf("lambda0-> %f  ...sigma %lf\n",lambda0,sigma);

	/*
	for(i=0;i<neje;i++){
		printf("eje (%d) %f  ...\n",i,eje[i]);
	}
*/
	//int half = nLambda/2 +1;
	term = (REAL *)calloc(*sizeG, sizeof(REAL));

	for (i = 0; i < *sizeG; i++)
	{
		REAL lambdaX = lambda0 +i*step_between_lw;
		REAL aux = ((lambdaX - lambdaCentral) / sigma);
		term[i] = ( aux * aux) / 2; //exponent
		//printf("term (%d) %f  ...\n",i,term[i]);
	}

	nloai = 0;
	loai = calloc(*sizeG, sizeof(REAL));
	for (i = 0; i < *sizeG; i++)
	{
		if (term[i] < 1e30)
		{
			nloai++;
			loai[i] = 1;
		}
	}

	if (nloai > 0)
	{
		nmtb = nloai;
		mtb = calloc(nmtb, sizeof(REAL));
		for (i = 0; i < *sizeG; i++)
		{
			if (loai[i])
			{
				mtb[i] = EXP(-term[i]);
				//printf("term (%d) %f  ...\n",i,mtb[i]);
			}
		}
	}
	else
	{

		nmtb = *sizeG;
		mtb = calloc(nmtb, sizeof(REAL));
		for (i = 0; i < *sizeG; i++)
		{
			mtb[i] = EXP(-term[i]);
			//printf("term (%d) %f  ...\n",i,mtb[i]);
		}
	}

	cte = 0;
	//normalization
	for (i = 0; i < nmtb; i++)
	{
		cte += mtb[i];
	}
	for (i = 0; i < *sizeG; i++)
	{
		mtb[i] /= cte;
	}

	free(loai);
	free(term);

	return mtb;
}