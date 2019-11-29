
//    _______             _______ _________ _        _______  _______
//   (  ____ \           (       )\__   __/( \      (  ___  )(  ____ \
//   | (    \/           | () () |   ) (   | (      | (   ) || (    \/
//   | |         _____   | || || |   | |   | |      | |   | || (_____
//   | |        (_____)  | |(_)| |   | |   | |      | |   | |(_____  )
//   | |                 | |   | |   | |   | |      | |   | |      ) |
//   | (____/\           | )   ( |___) (___| (____/\| (___) |/\____) |
//   (_______/           |/     \|\_______/(_______/(_______)\_______)
//
//
// CMILOS v1.0 (2019)
// RTE INVERSION C code for SOPHI (based on the ILD code MILOS by D. Orozco)
// Juanp && Manu (IAA-CSIC)
//
// How to use:
//
//  >> milos_Fits parameters.txt
//
//

#include <time.h>
#include "defines.h"
//#include "nrutil.h"

#include <string.h>
#include <stdio.h>
#include "fitsio.h"
#include "utilsFits.h"
#include "milosUtils.h"
#include "lib.h"

#ifdef _OPENMP
   #include <omp.h>
#else
   #define omp_get_thread_num() 0
#endif




long long int c1, c2, cd, semi, c1a, c2a, cda; //variables de 64 bits para leer ciclos de reloj
long long int c1total, cdtotal;

Cuantic *cuantic; // Variable global, está hecho así, de momento,para parecerse al original


PRECISION **PUNTEROS_CALCULOS_COMPARTIDOS;
int POSW_PUNTERO_CALCULOS_COMPARTIDOS;
int POSR_PUNTERO_CALCULOS_COMPARTIDOS;

PRECISION *gp1, *gp2, *dt, *dti, *gp3, *gp4, *gp5, *gp6, *etai_2;

//PRECISION gp4_gp2_rhoq[NLAMBDA],gp5_gp2_rhou[NLAMBDA],gp6_gp2_rhov[NLAMBDA];
PRECISION *gp4_gp2_rhoq, *gp5_gp2_rhou, *gp6_gp2_rhov;

PRECISION *dgp1, *dgp2, *dgp3, *dgp4, *dgp5, *dgp6, *d_dt;
PRECISION *d_ei, *d_eq, *d_eu, *d_ev, *d_rq, *d_ru, *d_rv;
PRECISION *dfi, *dshi;
PRECISION CC, CC_2, sin_gm, azi_2, sinis, cosis, cosis_2, cosi, sina, cosa, sinda, cosda, sindi, cosdi, sinis_cosa, sinis_sina;
PRECISION *fi_p, *fi_b, *fi_r, *shi_p, *shi_b, *shi_r;
PRECISION *etain, *etaqn, *etaun, *etavn, *rhoqn, *rhoun, *rhovn;
PRECISION *etai, *etaq, *etau, *etav, *rhoq, *rhou, *rhov;
PRECISION *parcial1, *parcial2, *parcial3;
PRECISION *nubB, *nupB, *nurB;
PRECISION **uuGlobalInicial;
PRECISION **HGlobalInicial;
PRECISION **FGlobalInicial;
PRECISION *perfil_instrumental;
PRECISION *G;
int FGlobal, HGlobal, uuGlobal;

PRECISION *d_spectra, *spectra;

//Convolutions values
int NMUESTRAS_G = 0;
PRECISION FWHM = 0;
PRECISION DELTA = 0;

int INSTRUMENTAL_CONVOLUTION = 0;
int INSTRUMENTAL_CONVOLUTION_WITH_PSF = 0;
int CLASSICAL_ESTIMATES = 0;
int PRINT_SINTESIS = 0;

int NUM_OPENMP_THREADS = 12;


int main(int argc, char **argv)
{

   //omp_set_num_threads(NUM_OPENMP_THREADS);
	PRECISION *wlines;
	int nwlines, nlambda, numPixels, iter, nweight, Max_iter;
	Init_Model *vModels;
	PRECISION chisqrf, * vChisqrf;
	PRECISION slight, toplim;
	PRECISION weight[4] = {1., 1., 1., 1.};
	PRECISION CENTRAL_WL;
	clock_t t_ini, t_fin;
	PRECISION secs, total_secs;

	// CONFIGURACION DE PARAMETROS A INVERTIR
	//INIT_MODEL=[eta0,magnet,vlos,landadopp,aa,gamma,azi,B1,B2,macro,alfa]
	int fix[] = {1., 1., 1., 1., 1., 1., 1., 1., 1., 0., 0.}; //Parametros invertidos
	//----------------------------------------------

	PRECISION sigma[NPARMS];
	PRECISION vsig;
	
	PRECISION ilambda;
	PRECISION noise;

	char  nameInputFileSpectra [PATH_MAX];
	char  nameInputFileLambda [PATH_MAX];
	char  nameOutputFileModels [PATH_MAX];
	char  nameOutputFilePerfiles [PATH_MAX];
	char	nameInputFileLines [PATH_MAX];

   FitsImage * fitsImage;

	PRECISION dat[7] = {CUANTIC_NWL, CUANTIC_SLOI, CUANTIC_LLOI, CUANTIC_JLOI, CUANTIC_SUPI, CUANTIC_LUPI, CUANTIC_JUPI};

	/********************* Read data input from file ******************************/

	/* Read data input from file */

	if(!readParametersFileInput(argv[1], &Max_iter,&CLASSICAL_ESTIMATES,&PRINT_SINTESIS,nameInputFileSpectra,nameInputFileLambda,nameInputFileLines, &CENTRAL_WL,nameOutputFileModels,nameOutputFilePerfiles,&INSTRUMENTAL_CONVOLUTION,&FWHM,&DELTA,&NMUESTRAS_G)){
		printf("\n********************* EXITING THE PROGRAM . ERROR READING PARAMETERS FILE ****************************\n");
		return -1;
	}

	/*********************************************** INITIALIZE VARIABLES  *********************************/
	toplim = 1e-18;

	CC = PI / 180.0;
	CC_2 = CC * 2;


	nweight = 4;

	nwlines = 1;
	wlines = (PRECISION *)calloc(2, sizeof(PRECISION));
	wlines[0] = 1;
	wlines[1] = CENTRAL_WL;

	vsig = NOISE_SIGMA; //original 0.001
	sigma[0] = vsig;
	sigma[1] = vsig;
	sigma[2] = vsig;
	sigma[3] = vsig;
	

	noise = NOISE_SIGMA;
	ilambda = ILAMBDA;
	iter = 0;
	numPixels=0;
	/****************************************************************************************************/	


	/******************* APPLY GAUSSIAN, CREATE CUANTINC AND INITIALIZE DINAMYC MEMORY*******************/


	cuantic = create_cuantic(dat,1);
	InitializePointerShareCalculation();	

	/****************************************************************************************************/



	// READ PIXELS FROM IMAGE 
	PRECISION timeReadImage;
	clock_t t;
	t = clock();
	fitsImage = readFitsSpectroImage(nameInputFileSpectra);
	t = clock() - t;
	timeReadImage = ((PRECISION)t)/CLOCKS_PER_SEC; // in seconds 
	printf("\n +++++++++++++++++++++++");
	printf("\n TIME TO READ FITS IMAGE:  %f seconds to execute \n", timeReadImage); 
	printf("\n +++++++++++++++++++++++");
	if(fitsImage!=NULL && readFitsLambdaFile(nameInputFileLambda,fitsImage)){
		// THE NUMBER OF LAMBDAS IS READ FROM INPUT FILES 
		nlambda = fitsImage->nLambdas;
		#pragma omp parallel private (PUNTEROS_CALCULOS_COMPARTIDOS,spectra,d_spectra,gp4_gp2_rhoq,gp5_gp2_rhou,gp6_gp2_rhov,gp1,gp2,gp3,gp4,gp5,gp6,dt,dti,etai_2,dgp1,dgp2,dgp3,dgp4,dgp5,dgp6,d_dt,d_ei,d_eq,d_eu,d_ev,d_rq,d_ru,d_rv,dfi,dshi,fi_p,fi_b,fi_r,shi_p,shi_b,shi_r,etain,etaqn,etaun,etavn,rhoqn,rhoun,rhovn,etai,etaq,etau,etav,rhoq,rhou,rhov,parcial1,parcial2,parcial3,nubB,nurB,nupB,uuGlobalInicial, uuGlobal,HGlobalInicial, HGlobal,FGlobalInicial,FGlobal)
		{

			AllocateMemoryDerivedSynthesis(nlambda);

			//initializing weights
			PRECISION *w, *sig;
			weights_init(sigma, &sig, noise);

			int indexPixel = 0;

			// ALLOCATE MEMORY FOR STORE THE RESULTS 

			vModels = calloc (fitsImage->numPixels , sizeof(Init_Model));
			vChisqrf = calloc (fitsImage->numPixels , sizeof(PRECISION));

			printf("\n START EXECUTION OF INVERSION ");
			printf("\n**********");
			t = clock();
			Init_Model initModel;
			#pragma omp parallel for private(initModel, chisqrf,FGlobal,HGlobal,uuGlobal,POSR_PUNTERO_CALCULOS_COMPARTIDOS,POSW_PUNTERO_CALCULOS_COMPARTIDOS)  num_threads(2)
			for(indexPixel = 0; indexPixel < fitsImage->numPixels; indexPixel++){
				

				//Initial Model
				
				initModel.eta0 = INITIAL_MODEL_ETHA0;
				initModel.B = INITIAL_MODEL_B; //200 700
				initModel.gm = INITIAL_MODEL_GM;
				initModel.az = INITIAL_MODEL_AZI;
				initModel.vlos = INITIAL_MODEL_VLOS; //km/s 0
				initModel.mac = 0.0;
				initModel.dopp = INITIAL_MODEL_LAMBDADOPP;
				initModel.aa = INITIAL_MODEL_AA;
				initModel.alfa = 1; //0.38; //stray light factor
				initModel.S0 = INITIAL_MODEL_S0;
				initModel.S1 = INITIAL_MODEL_S1;

				if (CLASSICAL_ESTIMATES)
				{
					t_ini = clock();
					estimacionesClasicas(wlines[1], fitsImage->pixels[indexPixel].vLambda, fitsImage->pixels[indexPixel].nLambda, fitsImage->pixels[indexPixel].spectro, &initModel);
					t_fin = clock();

					//Se comprueba si el resultado fue "nan" en las CE
					if (isnan(initModel.B))
						initModel.B = 1;
					if (isnan(initModel.vlos))
						initModel.vlos = 1e-3;
					if (isnan(initModel.gm))
						initModel.gm = 1;
					if (isnan(initModel.az))
						initModel.az = 1;

				}
				//inversion
				if (CLASSICAL_ESTIMATES != 2)
				{

					printf("\n inversion rte");
					printf("\n*******************************");
					//Se introduce en S0 el valor de Blos si solo se calculan estimaciones clásicas
					//Aqui se anula esa asignación porque se va a realizar la inversion RTE completa
					initModel.S0 = INITIAL_MODEL_S0;

					lm_mils(cuantic, wlines, fitsImage->pixels[indexPixel].vLambda, fitsImage->pixels[indexPixel].nLambda, fitsImage->pixels[indexPixel].spectro, fitsImage->pixels[indexPixel].nLambda, &initModel, spectra, &chisqrf, NULL, toplim, Max_iter,
							weight, fix, sig, filter, ilambda,&INSTRUMENTAL_CONVOLUTION,&NMUESTRAS_G);
					printf("\nESCRIBIENDO LOS PRIMEROS VALORES EN LOS VECTORES DE MODELOS Y DE CHI CUADRADAO");
					printf("\n******************************");
				}


				vModels[indexPixel] = initModel;
				vChisqrf[indexPixel] = chisqrf;

			}
		}
		t = clock() - t;
		timeReadImage = ((PRECISION)t)/CLOCKS_PER_SEC; // in seconds 
		printf("\n FINISH EXECUTION OF INVERSION: %f seconds to execute \n", timeReadImage);
		printf("\n**********");
		if(!writeFitsImageModels(nameOutputFileModels,fitsImage->rows,fitsImage->cols,vModels,vChisqrf)){
				printf("\n ERROR WRITING FILE OF MODELS: %s",nameOutputFileModels);
		}

		// PROCESS FILE OF SYNTETIC PROFILES

		if(PRINT_SINTESIS){
			int i;
			for( i=0;i<fitsImage->numPixels;i++)
			{

				Init_Model initModel = vModels[i];
				//PRECISION chisqr = vChisqrf[i];
				int NMODEL = 12; //Numero de parametros del modelo

				mil_sinrf(cuantic, &initModel, wlines, fitsImage->pixels[i].vLambda, nlambda, spectra, AH, 0, filter,NULL);

				me_der(cuantic, &initModel, wlines, fitsImage->pixels[i].vLambda, nlambda, d_spectra, spectra,AH, NULL, 0, filter);
				response_functions_convolution(&nlambda,&INSTRUMENTAL_CONVOLUTION,&NMUESTRAS_G);
				int kk;
				for (kk = 0; kk < (nlambda * NPARMS); kk++)
				{
					fitsImage->pixels[i].spectro[kk] = spectra[kk] ;
				}

			}
			// WRITE SINTHETIC PROFILES TO FITS FILE

			if(!writeFitsImageProfiles(nameOutputFilePerfiles,nameInputFileSpectra,fitsImage)){
				printf("\n ERROR WRITING FILE OF SINTHETIC PROFILES: %s",nameOutputFilePerfiles);
			}
		}
	}
	else{
		printf("\n\n ***************************** FITS FILE WITH THE SPECTRO IMAGE CAN NOT BE READ IT ******************************\n");
	}

	printf(" \n***********************  IMAGE INVERSION DONE, CLEANING MEMORY *********************\n");
	
	freeFitsImage(fitsImage);

	free(cuantic);
	free(wlines);
	FreeMemoryDerivedSynthesis();

	free(G);

	return 0;
}
