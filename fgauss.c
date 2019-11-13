
#include "defines.h"

/*
 * 
 * deriv : 1 true, 0 false
 */

//;this function builds a gauss function
//;landa(amstrong) ;Central wavelength
//;eje(amstrong) ;Wavelength axis
//;macro ;Macroturbulence in km/s

double * fgauss(double MC, double *eje, int neje, double landa, int deriv)
{
	//int fgauss(double MC, double * eje,int neje,double landa,int deriv,double * mtb,int nmtb){

	double centro, *mtb;
	double ild;
	double term[neje];
	double * loai;
	//double *term, *loai;
	int i;
	int nloai, nmtb;
	double cte;

	centro = eje[(int)neje / 2];		  //center of the axis
	ild = (landa * MC) / 2.99792458e5; //Sigma

	//	printf("ild-> %f  ...\n",ild);

	/*
	for(i=0;i<neje;i++){
		printf("eje (%d) %f  ...\n",i,eje[i]);
	}
*/
	//term = (double *) malloc (neje * sizeof(double));


	for (i = 0; i < neje; i++)
	{
		double aux = ((eje[i] - centro) / ild);
		term[i] = ( aux * aux) / 2; //exponent
		//printf("term (%d) %f  ...\n",i,term[i]);
	}

	nloai = 0;
	loai = calloc(neje, sizeof(double));
	for (i = 0; i < neje; i++)
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
		mtb = calloc(nmtb, sizeof(double));
		for (i = 0; i < neje; i++)
		{
			if (loai[i])
			{
				mtb[i] = exp(-term[i]);
			}
		}
	}
	else
	{

		nmtb = neje;
		mtb = malloc ( nmtb * sizeof(double));
		for (i = 0; i < neje; i++)
		{
			mtb[i] = exp(-term[i]);
		}
	}

	cte = 0;
	//normalization
	for (i = 0; i < nmtb; i++)
	{
		cte += mtb[i];
	}
	for (i = 0; i < neje; i++)
	{
		mtb[i] /= cte;
	}

	free(loai);
	//free(term);

	//In case we need the deriv of f gauss /deriv
	if (deriv == 1)
	{
		for (i = 0; i < nmtb; i++)
		{
			//mtb2=mtb/macro*(((eje-centro)/ILd)^2d0-1d0)
			mtb[i] = mtb[i] / MC * ((((eje[i] - centro) / ild) * ((eje[i] - centro) / ild)) - 1.0);			
		}
	}

	return mtb;
}







/*
 * 
 * deriv : 1 true, 0 false
 */

//;this function builds a gauss function
//;landa(amstrong) ;Central wavelength
//;eje(amstrong) ;Wavelength axis
//;macro ;Macroturbulence in km/s

double * fgauss_WL(double FWHM, double step_between_lw, double lambda0, double lambdaCentral, int nLambda, int * sizeG)
{
	//int fgauss(double MC, double * eje,int neje,double landa,int deriv,double * mtb,int nmtb){

	double *mtb;
	double *term, *loai;
	int i;
	int nloai, nmtb;
	double cte;
	

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
	//double sigma=FWHM*0.42466090/1000.0; // in Angstroms
	double sigma = FWHM / (2 * sqrt(2 * log(2)));

	//printf("lambda0-> %f  ...sigma %lf\n",lambda0,sigma);

	/*
	for(i=0;i<neje;i++){
		printf("eje (%d) %f  ...\n",i,eje[i]);
	}
*/
	//int half = nLambda/2 +1;
	term = (double *)calloc(*sizeG, sizeof(double));

	for (i = 0; i < *sizeG; i++)
	{
		double lambdaX = lambda0 +i*step_between_lw;
		double aux = ((lambdaX - lambdaCentral) / sigma);
		term[i] = ( aux * aux) / 2; //exponent
		//printf("term (%d) %f  ...\n",i,term[i]);
	}

	nloai = 0;
	loai = calloc(*sizeG, sizeof(double));
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
		mtb = calloc(nmtb, sizeof(double));
		for (i = 0; i < *sizeG; i++)
		{
			if (loai[i])
			{
				mtb[i] = exp(-term[i]);
				//printf("term (%d) %f  ...\n",i,mtb[i]);
			}
		}
	}
	else
	{

		nmtb = *sizeG;
		mtb = calloc(nmtb, sizeof(double));
		for (i = 0; i < *sizeG; i++)
		{
			mtb[i] = exp(-term[i]);
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