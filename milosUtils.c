#include "defines.h"
#include "milosUtils.h"
#include "nrutil.h"
#include "svdcordic.c"
#include "svdcmp.h"
#include "lib.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>


#define tiempo(ciclos) asm volatile("rdtsc \n\t" \
												: "=A"(ciclos))

// PSF obtenida desde los datos teoricos de CRISP CRISP_6173_28mA.psf
// Se usa el scrip interpolar_psf.m
const PRECISION crisp_psf[141] = {0.0004, 0.0004, 0.0005, 0.0005, 0.0005, 0.0005, 0.0006, 0.0006, 0.0006, 0.0007, 0.0007, 0.0008, 0.0008, 0.0009, 0.0009, 0.0010, 0.0010, 0.0011, 0.0012, 0.0012, 0.0013, 0.0014,
											 0.0015, 0.0016, 0.0017, 0.0019, 0.0020, 0.0021, 0.0023, 0.0025, 0.0027, 0.0029, 0.0031, 0.0034, 0.0037, 0.0040, 0.0044, 0.0048, 0.0052, 0.0057, 0.0062, 0.0069, 0.0076, 0.0083,
											 0.0092, 0.0102, 0.0114, 0.0127, 0.0142, 0.0159, 0.0178, 0.0202, 0.0229, 0.0261, 0.0299, 0.0343, 0.0398, 0.0465, 0.0546, 0.0645, 0.0765, 0.0918, 0.1113, 0.1363, 0.1678, 0.2066,
											 0.2501, 0.2932, 0.3306, 0.3569, 0.3669, 0.3569, 0.3306, 0.2932, 0.2501, 0.2066, 0.1678, 0.1363, 0.1113, 0.0918, 0.0765, 0.0645, 0.0546, 0.0465, 0.0398, 0.0343, 0.0299, 0.0261,
											 0.0229, 0.0202, 0.0178, 0.0159, 0.0142, 0.0127, 0.0114, 0.0102, 0.0092, 0.0083, 0.0076, 0.0069, 0.0062, 0.0057, 0.0052, 0.0048, 0.0044, 0.0040, 0.0037, 0.0034, 0.0031, 0.0029,
											 0.0027, 0.0025, 0.0023, 0.0021, 0.0020, 0.0019, 0.0017, 0.0016, 0.0015, 0.0014, 0.0013, 0.0012, 0.0012, 0.0011, 0.0010, 0.0010, 0.0009, 0.0009, 0.0008, 0.0008, 0.0007, 0.0007,
											 0.0006, 0.0006, 0.0006, 0.0005, 0.0005, 0.0005, 0.0005, 0.0004, 0.0004};

extern long long int c1, c2, cd, semi, c1a, c2a, cda; //variables de 64 bits para leer ciclos de reloj
extern long long int c1total, cdtotal;

extern PRECISION **PUNTEROS_CALCULOS_COMPARTIDOS;
extern int POSW_PUNTERO_CALCULOS_COMPARTIDOS;
extern int POSR_PUNTERO_CALCULOS_COMPARTIDOS;

extern PRECISION *gp4_gp2_rhoq, *gp5_gp2_rhou, *gp6_gp2_rhov;

extern PRECISION *gp1, *gp2, *dt, *dti, *gp3, *gp4, *gp5, *gp6, *etai_2;
extern PRECISION *dgp1, *dgp2, *dgp3, *dgp4, *dgp5, *dgp6, *d_dt;
extern PRECISION *d_ei, *d_eq, *d_eu, *d_ev, *d_rq, *d_ru, *d_rv;
extern PRECISION *dfi, *dshi;
extern PRECISION *fi_p, *fi_b, *fi_r, *shi_p, *shi_b, *shi_r;
extern PRECISION *spectra, *d_spectra;
extern PRECISION *etain, *etaqn, *etaun, *etavn, *rhoqn, *rhoun, *rhovn;
extern PRECISION *etai, *etaq, *etau, *etav, *rhoq, *rhou, *rhov;
extern PRECISION *parcial1, *parcial2, *parcial3;
extern PRECISION *nubB, *nupB, *nurB;
extern PRECISION *G;

extern Cuantic *cuantic; // Variable global, está hecho así, de momento,para parecerse al original

void spectral_synthesis_convolution(int * nlambda, int * INSTRUMENTAL_CONVOLUTION, int * NMUESTRAS_G)
{

	int i;
	//int nlambda = NLAMBDA;

	//convolucionamos los perfiles IQUV (spectra)
	if (*INSTRUMENTAL_CONVOLUTION)
	{

		PRECISION Ic;

		if (!INSTRUMENTAL_CONVOLUTION_INTERPOLACION)
		{
			//convolucion de I
			Ic = spectra[*nlambda - 1];

			for (i = 0; i < *nlambda - 1; i++)
				spectra[i] = Ic - spectra[i];

			direct_convolution(spectra, *nlambda - 1, G, *NMUESTRAS_G, 1); //no convolucionamos el ultimo valor Ic

			for (i = 0; i < *nlambda - 1; i++)
				spectra[i] = Ic - spectra[i];

			//convolucion QUV
			for (i = 1; i < NPARMS; i++)
				direct_convolution(spectra + *nlambda * i, *nlambda - 1, G, *NMUESTRAS_G, 1); //no convolucionamos el ultimo valor
		}
		else
		{
			if ( *nlambda == 6)
			{

				//convolucion de I
				Ic = spectra[ *nlambda - 1];

				for (i = 0; i < *nlambda - 1; i++)
					spectra[i] = Ic - spectra[i];

				PRECISION *spectra_aux;
				spectra_aux = (PRECISION *)calloc(*nlambda * 2 - 2, sizeof(PRECISION));

				int j = 0;
				for (i = 0, j = 0; i < *nlambda * 2 - 2; i = i + 2, j++)
					spectra_aux[i] = spectra[j];

				for (i = 1, j = 0; i < *nlambda * 2 - 2; i = i + 2, j++)
					spectra_aux[i] = (spectra[j] + spectra[j + 1]) / 2;

				direct_convolution(spectra_aux, *nlambda * 2 - 2 - 1, G, *NMUESTRAS_G, 1); //no convolucionamos el ultimo valor Ic

				for (i = 0, j = 0; i < *nlambda * 2 - 2; i = i + 2, j++)
					spectra[j] = spectra_aux[i];

				for (i = 0; i < *nlambda - 1; i++)
					spectra[i] = Ic - spectra[i];

				free(spectra_aux);

				//convolucion QUV
				int k;
				for (k = 1; k < NPARMS; k++)
				{

					PRECISION *spectra_aux;
					spectra_aux = (PRECISION *)calloc(*nlambda * 2 - 2, sizeof(PRECISION));

					int j = 0;
					for (i = 0, j = 0; i < *nlambda * 2 - 2; i = i + 2, j++)
						spectra_aux[i] = spectra[j + *nlambda * k];

					for (i = 1, j = 0; i < *nlambda * 2 - 2; i = i + 2, j++)
						spectra_aux[i] = (spectra[j + *nlambda * k] + spectra[j + 1 + *nlambda * k]) / 2;

					direct_convolution(spectra_aux, *nlambda * 2 - 2 - 1, G, *NMUESTRAS_G, 1); //no convolucionamos el ultimo valor Ic

					for (i = 0, j = 0; i < *nlambda * 2 - 2; i = i + 2, j++)
						spectra[j + *nlambda * k] = spectra_aux[i];

					free(spectra_aux);
				}
			}
		}
	}
}

void response_functions_convolution(int * nlambda, int * INSTRUMENTAL_CONVOLUTION, int * NMUESTRAS_G)
{

	int i, j;
	//int nlambda = NLAMBDA;

	//convolucionamos las funciones respuesta ( d_spectra )
	if (*INSTRUMENTAL_CONVOLUTION)
	{
		if (!INSTRUMENTAL_CONVOLUTION_INTERPOLACION)
		{

			for (j = 0; j < NPARMS; j++)
			{
				for (i = 0; i < NTERMS; i++)
				{
					if (i != 7)																															 //no convolucionamos S0
						direct_convolution(d_spectra + *nlambda * i + *nlambda * NTERMS * j, *nlambda - 1, G, *NMUESTRAS_G, 1); //no convolucionamos el ultimo valor
				}
			}
		}
		else
		{

			int k, m;
			for (k = 0; k < NPARMS; k++)
			{
				for (m = 0; m < NTERMS; m++)
				{

					PRECISION *spectra_aux;
					spectra_aux = (PRECISION *)calloc(*nlambda * 2 - 2, sizeof(PRECISION));

					int j = 0;
					for (i = 0, j = 0; i < *nlambda * 2 - 2; i = i + 2, j++)
						spectra_aux[i] = d_spectra[j + *nlambda * m + *nlambda * NTERMS * k];

					for (i = 1, j = 0; i < *nlambda * 2 - 2; i = i + 2, j++)
						spectra_aux[i] = (d_spectra[j + *nlambda * m + *nlambda * NTERMS * k] + d_spectra[j + *nlambda * m + *nlambda * NTERMS * k]) / 2;

					direct_convolution(spectra_aux, *nlambda * 2 - 2 - 1, G, *NMUESTRAS_G, 1); //no convolucionamos el ultimo valor Ic

					for (i = 0, j = 0; i < *nlambda * 2 - 2; i = i + 2, j++)
						d_spectra[j + *nlambda * m + *nlambda * NTERMS * k] = spectra_aux[i];

					free(spectra_aux);
				}
			}
		}
	}
}



void AplicaDelta(Init_Model *model, PRECISION *delta, int *fixed, Init_Model *modelout)
{

	//INIT_MODEL=[eta0,magnet,vlos,landadopp,aa,gamma,azi,B1,B2,macro,alfa]

	if (fixed[0])
	{
		modelout->eta0 = model->eta0 - delta[0]; // 0
	}
	if (fixed[1])
	{
		if (delta[1] < -800) //300
			delta[1] = -800;
		else if (delta[1] > 800)
			delta[1] = 800;
		modelout->B = model->B - delta[1]; //magnetic field
	}
	if (fixed[2])
	{

		if (delta[2] > 2)
			delta[2] = 2;

		if (delta[2] < -2)
			delta[2] = -2;

		modelout->vlos = model->vlos - delta[2];
	}

	if (fixed[3])
	{

		if (delta[3] > 1e-2)
			delta[3] = 1e-2;
		else if (delta[3] < -1e-2)
			delta[3] = -1e-2;

		modelout->dopp = model->dopp - delta[3];
	}

	if (fixed[4])
		modelout->aa = model->aa - delta[4];

	if (fixed[5])
	{
		if (delta[5] < -15) //15
			delta[5] = -15;
		else if (delta[5] > 15)
			delta[5] = 15;

		modelout->gm = model->gm - delta[5]; //5
	}
	if (fixed[6])
	{
		if (delta[6] < -15)
			delta[6] = -15;
		else if (delta[6] > 15)
			delta[6] = 15;

		modelout->az = model->az - delta[6];
	}
	if (fixed[7])
		modelout->S0 = model->S0 - delta[7];
	if (fixed[8])
		modelout->S1 = model->S1 - delta[8];
	if (fixed[9])
		modelout->mac = model->mac - delta[9]; //9
	if (fixed[10])
		modelout->alfa = model->alfa - delta[10];
}


int check(Init_Model *model)
{

	double offset = 0;
	double inter;

	//Magnetic field
	if (model->B < 0)
	{
		//model->B = 190;
		model->B = -(model->B);
		model->gm = 180.0 - (model->gm);
	}
	if (model->B > 5000)
		model->B = 5000;

	//Inclination
	if (model->gm < 0)
		model->gm = -(model->gm);
	if (model->gm > 180)
	{
		model->gm = 360.0 - model->gm;
		// model->gm = 179; //360.0 - model->gm;
	}

	//azimuth
	if (model->az < 0)
		model->az = 180 + (model->az); //model->az= 180 + (model->az);
	if (model->az > 180)
	{
		model->az = model->az - 180.0;
		// model->az = 179.0;
	}

	//RANGOS
	//Eta0
	if (model->eta0 < 1)
		model->eta0 = 1;

	// if(model->eta0 >8)
	// model->eta0=8;
	if (model->eta0 > 2500) //idl 2500
		model->eta0 = 2500;

	//velocity
	if (model->vlos < (-20)) //20
		model->vlos = (-20);
	if (model->vlos > 20)
		model->vlos = 20;

	//doppler width ;Do NOT CHANGE THIS
	if (model->dopp < 0.0001)
		model->dopp = 0.0001;

	if (model->dopp > 0.6) // idl 0.6
		model->dopp = 0.6;

	if (model->aa < 0.0001) // idl 1e-4
		model->aa = 0.0001;
	if (model->aa > 10) //10
		model->aa = 10;

	//S0
	if (model->S0 < 0.0001)
		model->S0 = 0.0001;
	if (model->S0 > 1.500)
		model->S0 = 1.500;

	//S1
	if (model->S1 < 0.0001)
		model->S1 = 0.0001;
	if (model->S1 > 2.000)
		model->S1 = 2.000;

	//macroturbulence
	if (model->mac < 0)
		model->mac = 0;
	if (model->mac > 4)
		model->mac = 4;

	return 1;
}


void FijaACeroDerivadasNoNecesarias(PRECISION *d_spectra, int *fixed, int nlambda)
{

	int In, j, i;
	for (In = 0; In < NTERMS; In++)
		if (fixed[In] == 0)
			for (j = 0; j < 4; j++)
				for (i = 0; i < nlambda; i++)
					d_spectra[i + nlambda * In + j * nlambda * NTERMS] = 0;
}


/*
	Tamaño de H es 	 NTERMS x NTERMS
	Tamaño de beta es 1xNTERMS

	return en delta tam 1xNTERMS
*/

int mil_svd(PRECISION *h, PRECISION *beta, PRECISION *delta)
{

	double epsilon, top;
	static PRECISION v2[TAMANIO_SVD][TAMANIO_SVD], w2[TAMANIO_SVD], v[NTERMS * NTERMS], w[NTERMS];
	static PRECISION h1[NTERMS * NTERMS], h_svd[TAMANIO_SVD * TAMANIO_SVD];
	static PRECISION aux[NTERMS * NTERMS];
	int i, j;
	//	static double aux2[NTERMS*NTERMS];
	static PRECISION aux2[NTERMS];
	int aux_nf, aux_nc;
	PRECISION factor, maximo, minimo;
	int posi, posj;

	epsilon = 1e-12;
	top = 1.0;

	factor = 0;
	maximo = 0;
	minimo = 1000000000;

	/**/
	for (j = 0; j < NTERMS * NTERMS; j++)
	{
		h1[j] = h[j];
	}

	if (USE_SVDCMP)
	{

		svdcmp(h1, NTERMS, NTERMS, w, v);
	}
	else
	{
		//printf(" NORMALIZACION y CORDIC######################################\n");
		//	NORMALIZACION
		for (j = 0; j < NTERMS * NTERMS; j++)
		{
			if (fabs(h[j]) > maximo)
			{
				maximo = fabs(h[j]);
			}
		}

		factor = maximo;

		//printf("maximo : %.12e \n",maximo);
		//exit(-1);

		if (!NORMALIZATION_SVD)
			factor = 1;

		for (j = 0; j < NTERMS * NTERMS; j++)
		{
			h1[j] = h[j] / (factor);
		}

		for (i = 0; i < TAMANIO_SVD - 1; i++)
		{
			for (j = 0; j < TAMANIO_SVD; j++)
			{
				if (j < NTERMS)
					h_svd[i * TAMANIO_SVD + j] = h1[i * NTERMS + j];
				else
					h_svd[i * TAMANIO_SVD + j] = 0;
			}
		}

		for (j = 0; j < TAMANIO_SVD; j++)
		{
			h_svd[(TAMANIO_SVD - 1) * TAMANIO_SVD + j] = 0;
		}

		svdcordic(h_svd, TAMANIO_SVD, TAMANIO_SVD, w2, v2, NUM_ITER_SVD_CORDIC);

		for (i = 0; i < TAMANIO_SVD - 1; i++)
		{
			for (j = 0; j < TAMANIO_SVD - 1; j++)
			{
				v[i * NTERMS + j] = v2[i][j];
			}
		}

		for (j = 0; j < TAMANIO_SVD - 1; j++)
		{
			w[j] = w2[j] * factor;
		}
	}

	static PRECISION vaux[NTERMS * NTERMS], waux[NTERMS];

	for (j = 0; j < NTERMS * NTERMS; j++)
	{
		vaux[j] = v[j]; //*factor;
	}

	for (j = 0; j < NTERMS; j++)
	{
		waux[j] = w[j]; //*factor;
	}

	multmatrix(beta, 1, NTERMS, vaux, NTERMS, NTERMS, aux2, &aux_nf, &aux_nc);

	for (i = 0; i < NTERMS; i++)
	{
      aux2[i]= aux2[i]*((fabs(waux[i]) > epsilon) ? (1/waux[i]): 0.0);
	}

	multmatrix(vaux, NTERMS, NTERMS, aux2, NTERMS, 1, delta, &aux_nf, &aux_nc);

	return 1;
}



void weights_init(int nlambda, double *sigma, PRECISION *weight, int nweight, PRECISION **wOut, PRECISION **sigOut, double noise)
{
	int i, j;
	PRECISION *w, *sig;

	sig = calloc(4, sizeof(PRECISION));
	if (sigma == NULL)
	{
		for (i = 0; i < 4; i++)
			sig[i] = noise * noise;
	}
	else
	{

		for (i = 0; i < 4; i++)
			sig[i] = (*sigma); // * (*sigma);
	}

	*wOut = w;
	*sigOut = sig;
}


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
void estimacionesClasicas(PRECISION lambda_0, double *lambda, int nlambda, PRECISION *spectro, Init_Model *initModel)
{

	PRECISION x, y, aux, LM_lambda_plus, LM_lambda_minus, Blos, beta_B, Ic, Vlos;
	PRECISION *spectroI, *spectroQ, *spectroU, *spectroV;
	PRECISION L, m, gamma, gamma_rad, tan_gamma, maxV, minV, C, maxWh, minWh;
	int i, j;

	//Es necesario crear un lambda en FLOAT para probar como se hace en la FPGA
	PRECISION *lambda_aux = lambda;
	//lambda_aux = (PRECISION *)calloc(nlambda, sizeof(PRECISION));
	PRECISION lambda0, lambda1, lambda2, lambda3, lambda4;

	/*lambda0 = 6.1732012e+3 + 0;		// RTE_WL_0
	lambda1 = lambda0 + 0.070000000; //RTE_WL_STEP
	lambda2 = lambda1 + 0.070000000;
	lambda3 = lambda2 + 0.070000000;
	lambda4 = lambda3 + 0.070000000;

	lambda_aux[0] = lambda0;
	lambda_aux[1] = lambda1;
	lambda_aux[2] = lambda2;
	lambda_aux[3] = lambda3;
	lambda_aux[4] = lambda4;*/

	spectroI = spectro;
	spectroQ = spectro + nlambda;
	spectroU = spectro + nlambda * 2;
	spectroV = spectro + nlambda * 3;

	Ic = spectro[nlambda - 1]; // Continuo ultimo valor de I

	x = 0;
	y = 0;
	for (i = 0; i < nlambda - 1; i++)
	{
		aux = (Ic - (spectroI[i] + spectroV[i]));
		x = x + aux * (lambda_aux[i] - lambda_0);
		y = y + aux;
	}

	//Para evitar nan
	if (fabs(y) > 1e-15)
		LM_lambda_plus = x / y;
	else
		LM_lambda_plus = 0;

	x = 0;
	y = 0;
	for (i = 0; i < nlambda - 1; i++)
	{
		aux = (Ic - (spectroI[i] - spectroV[i]));
		x = x + aux * (lambda_aux[i] - lambda_0);
		y = y + aux;
	}

	if (fabs(y) > 1e-15)
		LM_lambda_minus = x / y;
	else
		LM_lambda_minus = 0;

	C = (CTE4_6_13 * lambda_0 * lambda_0 * cuantic->GEFF);
	beta_B = 1 / C;

	Blos = beta_B * ((LM_lambda_plus - LM_lambda_minus) / 2);
	Vlos = (VLIGHT / (lambda_0)) * ((LM_lambda_plus + LM_lambda_minus) / 2);


	//------------------------------------------------------------------------------------------------------------
	// //Para probar fórmulación propuesta por D. Orozco (Junio 2017)
	//La formula es la 2.7 que proviene del paper:
	// Diagnostics for spectropolarimetry and magnetography by Jose Carlos del Toro Iniesta and Valent´ýn Mart´ýnez Pillet
	//el 0.08 Es la anchura de la línea en lugar de la resuloción del etalón.

	//Vlos = ( 2*(VLIGHT)*0.08 / (PI*lambda_0)) * atan((spectroI[0]+spectroI[1]-spectroI[3]-spectroI[4])/(spectroI[0]-spectroI[1]-spectroI[3]+spectroI[4]));

	//------------------------------------------------------------------------------------------------------------

	Blos = Blos * 1; //factor de correción x campo debil
	Vlos = Vlos * 1; //factor de correción ...

	//inclinacion
	x = 0;
	y = 0;
	for (i = 0; i < nlambda - 1; i++)
	{
		L = fabs(sqrtf(spectroQ[i] * spectroQ[i] + spectroU[i] * spectroU[i]));
		m = fabs((4 * (lambda_aux[i] - lambda_0) * L)); // / (3*C*Blos) ); //2*3*C*Blos mod abril 2016 (en test!)

		x = x + fabs(spectroV[i]) * m;
		y = y + fabs(spectroV[i]) * fabs(spectroV[i]);

	}

	y = y * fabs((3 * C * Blos));

	tan_gamma = fabs(sqrtf(x / y));

	gamma_rad = atan(tan_gamma); //gamma en radianes

	gamma = gamma_rad * (180 / PI); //gamma en grados

	//correccion
	//utilizamos el signo de Blos para ver corregir el cuadrante
	PRECISION gamma_out = gamma;

	if (Blos < 0)
		gamma = (180) - gamma;

	//azimuth

	PRECISION tan2phi, phi;
	int muestra;

	if (nlambda == 6)
		muestra = CLASSICAL_ESTIMATES_SAMPLE_REF;
	else
		muestra = nlambda * 0.75;

	tan2phi = spectroU[muestra] / spectroQ[muestra];

	phi = (atan(tan2phi) * 180 / PI) / 2; //atan con paso a grados

	if (spectroU[muestra] > 0 && spectroQ[muestra] > 0)
		phi = phi;
	else if (spectroU[muestra] < 0 && spectroQ[muestra] > 0)
		phi = phi + 180;
	else if (spectroU[muestra] < 0 && spectroQ[muestra] < 0)
		phi = phi + 90;
	else if (spectroU[muestra] > 0 && spectroQ[muestra] < 0)
		phi = phi + 90;

	PRECISION B_aux;

	B_aux = fabs(Blos / cos(gamma_rad)) * 2; // 2 factor de corrección

	//Vlos = Vlos * 1.5;
	if (Vlos < (-20))
		Vlos = -20;
	if (Vlos > (20))
		Vlos = (20);


	initModel->B = (B_aux > 4000 ? 4000 : B_aux);
	initModel->vlos = Vlos; //(Vlos*1.5);//1.5;
	initModel->gm = gamma;
	initModel->az = phi;
	initModel->S0 = Blos;

	//Liberar memoria del vector de lambda auxiliar
	//free(lambda_aux);
}




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
				double getshi, int triplete, int * INSTRUMENTAL_CONVOLUTION, int * NMUESTRAS_G)
{

	int *diag;
	int iter;
	int i, j, In, *fixed, nfree;
	static PRECISION delta[NTERMS];
	double max[3], aux;
	int repite, pillado, nw, nsig;
	double *landa_store, flambda;
	static PRECISION beta[NTERMS], alpha[NTERMS * NTERMS];
	double chisqr, ochisqr;
	int nspectra, nd_spectra, clanda, ind;
	Init_Model model;

	//Genera aleatoriamente los componentes del vector
	tiempo(semi); //semilla para  generar los valores de la lista de forma aleatoria con srand
	srand((char)semi);

	iter = 0;

	nfree = CalculaNfree(spectro, nspectro);


	if (nfree == 0)
	{
		return -1; //'NOT ENOUGH POINTS'
	}

	flambda = ilambda;

	if (fix == NULL)
	{
		fixed = calloc(NTERMS, sizeof(double));
		for (i = 0; i < NTERMS; i++)
		{
			fixed[i] = 1;
		}
	}
	else
	{
		fixed = fix;
	}

	clanda = 0;
	iter = 0;
	repite = 1;
	pillado = 0;

	static PRECISION covar[NTERMS * NTERMS];
	static PRECISION betad[NTERMS];

	PRECISION chisqr_mem;
	int repite_chisqr = 0;

	mil_sinrf(cuantic, initModel, wlines, nwlines, lambda, nlambda, spectra, AH, triplete, filter);


	//convolucionamos los perfiles IQUV (spectra)
	spectral_synthesis_convolution(&nlambda,INSTRUMENTAL_CONVOLUTION,NMUESTRAS_G);

	me_der(cuantic, initModel, wlines, nwlines, lambda, nlambda, d_spectra, AH, slight, triplete, filter);

	//convolucionamos las funciones respuesta ( d_spectra )
	response_functions_convolution(&nlambda,INSTRUMENTAL_CONVOLUTION,NMUESTRAS_G);

	//FijaACeroDerivadasNoNecesarias(d_spectra,fixed,nlambda);

	covarm(weight, sigma, nsig, spectro, nlambda, spectra, d_spectra, beta, alpha);

	for (i = 0; i < NTERMS; i++)
		betad[i] = beta[i];

	for (i = 0; i < NTERMS * NTERMS; i++)
		covar[i] = alpha[i];

	/**************************************************************************/

	ochisqr = fchisqr(spectra, nspectro, spectro, weight, sigma, nfree);

	model = *initModel;
	do
	{
		chisqr_mem = (PRECISION)ochisqr;

		/**************************************************************************/
		for (i = 0; i < NTERMS; i++)
		{
			ind = i * (NTERMS + 1);
			covar[ind] = alpha[ind] * (1.0 + flambda);
		}

		mil_svd(covar, betad, delta);

		AplicaDelta(initModel, delta, fixed, &model);

		check(&model);

		/**************************************************************************/

		mil_sinrf(cuantic, &model, wlines, nwlines, lambda, nlambda, spectra, AH, triplete, filter);

		//convolucionamos los perfiles IQUV (spectra)
		spectral_synthesis_convolution(&nlambda,INSTRUMENTAL_CONVOLUTION,NMUESTRAS_G);

		chisqr = fchisqr(spectra, nspectro, spectro, weight, sigma, nfree);

		/**************************************************************************/
		if (chisqr - ochisqr < 0)
		{

			flambda = flambda / 10.0;

			*initModel = model;

			//printf("iteration=%d , chisqr = %f CONVERGE	- lambda= %e \n",iter,chisqr,flambda);

			me_der(cuantic, initModel, wlines, nwlines, lambda, nlambda, d_spectra, AH, slight, triplete, filter);

			//convolucionamos las funciones respuesta ( d_spectra )
			response_functions_convolution(&nlambda,INSTRUMENTAL_CONVOLUTION,NMUESTRAS_G);

			//FijaACeroDerivadasNoNecesarias(d_spectra,fixed,nlambda);

			covarm(weight, sigma, nsig, spectro, nlambda, spectra, d_spectra, beta, alpha);
			for (i = 0; i < NTERMS; i++)
				betad[i] = beta[i];

			for (i = 0; i < NTERMS * NTERMS; i++)
				covar[i] = alpha[i];

			ochisqr = chisqr;
		}
		else
		{
			flambda = flambda * 10; //10;

			//printf("iteration=%d , chisqr = %f NOT CONVERGE	- lambda= %e \n",iter,ochisqr,flambda);
		}

		/**************************************************************************/
		iter++;


	} while (iter <= miter); // && !clanda);

	*iterOut = iter;

	*chisqrf = ochisqr;

	if (fix == NULL)
		free(fixed);

	return 1;
}


void generateGaussianInstrumentalProfile(PRECISION * G, PRECISION  FWHM, PRECISION DELTA, int NMUESTRAS_G){
		G = vgauss(FWHM, NMUESTRAS_G, DELTA);

		/*if (INSTRUMENTAL_CONVOLUTION_WITH_PSF)
		{
			//if you wish to convolution with other instrumental profile you have to declare here and to asign it to "G"
			free(G);
			NMUESTRAS_G = 9;

			G = vgauss(FWHM, NMUESTRAS_G, DELTA); //solo para reservar memoria

			int kk;
			PRECISION sum = 0;
			for (kk = 0; kk < NMUESTRAS_G; kk++)
			{
				int pos = 70 - (int)((DELTA / 0.005) * (((int)(NMUESTRAS_G / 2)) - kk));

				if (pos < 0 || pos > 140) //140 length de crisp_psf
					G[kk] = 0;
				else
					G[kk] = crisp_psf[pos]; //70 es el centro de psf

				sum += G[kk];
			}

			for (kk = 0; kk < NMUESTRAS_G; kk++)
				G[kk] /= sum;
		}
		*/
}



/**
 * Make the interpolation between deltaLambda and PSF where deltaLambda es x and PSF f(x)
 *  Return the array with the interpolation. 
 * */
int interpolationPSF(PRECISION *deltaLambda, PRECISION * PSF, PRECISION * lambdasSamples, PRECISION centralLambda, size_t N_PSF, PRECISION * fInterpolated, size_t NSamples){

	size_t i;
	gsl_interp_accel *acc = gsl_interp_accel_alloc();
  	gsl_spline *spline_cubic = gsl_spline_alloc(gsl_interp_cspline, N_PSF);
	//gsl_spline *spline_akima = gsl_spline_alloc(gsl_interp_akima, NSamples);
	//gsl_spline *spline_steffen = gsl_spline_alloc(gsl_interp_steffen, NSamples);

	gsl_spline_init(spline_cubic, deltaLambda, PSF, N_PSF);
	//gsl_spline_init(spline_akima, deltaLambda, PSF, N_PSF);
	//gsl_spline_init(spline_steffen, deltaLambda, PSF, N_PSF);

	for (i = 0; i < NSamples; ++i){
   	PRECISION xi = lambdasSamples[i]-centralLambda;
      fInterpolated[i] = gsl_spline_eval(spline_cubic, xi, acc);
      //double yi_akima = gsl_spline_eval(spline_akima, xi, acc);
      //double yi_steffen = gsl_spline_eval(spline_steffen, xi, acc);
    }

  gsl_spline_free(spline_cubic);
  //gsl_spline_free(spline_akima);
  //gsl_spline_free(spline_steffen);
  gsl_interp_accel_free(acc);

	return 1;
}