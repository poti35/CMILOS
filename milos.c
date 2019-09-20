
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
#include "convolution.c"
#include <string.h>
#include <stdio.h>
#include "fitsio.h"
#include "utilsFits.h"
#include "milosUtils.h"
#include "lib.h"
#include <unistd.h>


long long int c1, c2, cd, semi, c1a, c2a, cda; //variables de 64 bits para leer ciclos de reloj
long long int c1total, cdtotal;

Cuantic *cuantic; // Variable global, está hecho así, de momento,para parecerse al original
char *concatena(char *a, int n, char *b);

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
PRECISION *interpolatedPSF;
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



int main(int argc, char **argv)
{

	double *wlines;
	int nwlines, nlambda, numPixels, iter, nweight, Max_iter;
	Init_Model *vModels;
	double chisqrf, * vChisqrf;
	double slight, toplim;
	PRECISION weight[4] = {1., 1., 1., 1.};

	PRECISION CENTRAL_WL;

	Init_Model INITIAL_MODEL;

	PRECISION * deltaLambda, * PSF;
	int N_SAMPLES_PSF;

	clock_t t_ini, t_fin;
	double secs, total_secs;

	// CONFIGURACION DE PARAMETROS A INVERTIR
	//INIT_MODEL=[eta0,magnet,vlos,landadopp,aa,gamma,azi,B1,B2,macro,alfa]
	int fix[] = {1., 1., 1., 1., 1., 1., 1., 1., 1., 0., 0.}; //Parametros invertidos
	//----------------------------------------------

	double sigma[NPARMS];
	double vsig;
	double filter;
	double ilambda;
	double noise;
	double *pol;
	double getshi;

	char  nameInputFileSpectra [PATH_MAX];
	char  nameInputFileLambda [PATH_MAX];
	char  nameOutputFileModels [PATH_MAX];
	char  nameOutputFilePerfiles [PATH_MAX];
	char	nameInputFileLines [PATH_MAX];
	char 	nameInputFileInitModel [PATH_MAX];
	char	nameInputFilePSF [PATH_MAX];

   FitsImage * fitsImage;

	//double dat[7] = {CUANTIC_NWL, CUANTIC_SLOI, CUANTIC_LLOI, CUANTIC_JLOI, CUANTIC_SUPI, CUANTIC_LUPI, CUANTIC_JUPI};
	double dat[7];

	/********************* Read data input from file ******************************/

	/* Read data input from file */

	if(!readParametersFileInput(argv[1], &Max_iter,&CLASSICAL_ESTIMATES,&PRINT_SINTESIS,nameInputFileSpectra,nameInputFileLambda,nameInputFileLines,nameInputFileInitModel, &CENTRAL_WL, nameOutputFileModels,nameOutputFilePerfiles,&INSTRUMENTAL_CONVOLUTION,nameInputFilePSF,&FWHM,&DELTA,&NMUESTRAS_G)){
		printf("\n********************* EXITING THE PROGRAM . ERROR READING PARAMETERS FILE ****************************\n");
		return -1;
	}

	readFileCuanticLines(nameInputFileLines,dat,CENTRAL_WL);
	readInitialModel(&INITIAL_MODEL,nameInputFileInitModel);

	/****************************************************************************************************/	
	// if parameter name of psf has been apported we read the file, in other case create the gaussian with the parameters
	int usedPSF = 0; // initial false
	if(INSTRUMENTAL_CONVOLUTION){
		if(access(nameInputFilePSF,F_OK) != -1){
			// read the number of lines 
				FILE *fp;
				char ch;
				N_SAMPLES_PSF=0;
				//open file in read more
				fp=fopen(nameInputFilePSF,"r");
				if(fp==NULL)
				{
					printf("File \"%s\" does not exist!!!\n",nameInputFilePSF);
					return 0;
				}

				//read character by character and check for new line	
				while((ch=fgetc(fp))!=EOF)
				{
					if(ch=='\n')
						N_SAMPLES_PSF++;
				}
				
				//close the file
				fclose(fp);
				if(N_SAMPLES_PSF>0){
					deltaLambda = calloc(N_SAMPLES_PSF,sizeof(PRECISION));
					PSF = calloc(N_SAMPLES_PSF,sizeof(PRECISION));
					readPSFFile(deltaLambda,PSF,nameInputFilePSF);
					usedPSF = 1;
				}
				else{
					generateGaussianInstrumentalProfile(G,FWHM,DELTA,NMUESTRAS_G);		
				}
		}else
			generateGaussianInstrumentalProfile(G,FWHM,DELTA,NMUESTRAS_G);
	}
	/*********************************************** INITIALIZE VARIABLES  *********************************/
	toplim = 1e-18;

	CC = PI / 180.0;
	CC_2 = CC * 2;

	filter = 0;
	getshi = 0;
	nweight = 4;

	nwlines = 1;
	wlines = (double *)calloc(2, sizeof(double));
	wlines[0] = 1;
	wlines[1] = CENTRAL_WL;

	vsig = NOISE_SIGMA; //original 0.001
	sigma[0] = vsig;
	sigma[1] = vsig;
	sigma[2] = vsig;
	sigma[3] = vsig;
	pol = NULL;

	noise = NOISE_SIGMA;
	ilambda = ILAMBDA;
	iter = 0;
	numPixels=0;
	/****************************************************************************************************/	


	/******************* APPLY GAUSSIAN, CREATE CUANTINC AND INITIALIZE DINAMYC MEMORY*******************/
	cuantic = create_cuantic(dat);
	InitializePointerShareCalculation();

	/****************************************************************************************************/



	// READ PIXELS FROM IMAGE 
	double timeReadImage;
	clock_t t;
	t = clock();
	fitsImage = readFitsSpectroImage(nameInputFileSpectra);
	t = clock() - t;
	timeReadImage = ((double)t)/CLOCKS_PER_SEC; // in seconds 
	printf("\n +++++++++++++++++++++++");
	printf("\n TIME TO READ FITS IMAGE:  %f seconds to execute \n", timeReadImage); 
	printf("\n +++++++++++++++++++++++");
	if(fitsImage!=NULL && readFitsLambdaFile(nameInputFileLambda,fitsImage)){
		// THE NUMBER OF LAMBDAS IS READ FROM INPUT FILES 
		nlambda = fitsImage->nLambdas;
		AllocateMemoryDerivedSynthesis(nlambda);

		// INTERPOLATE PSF WITH ARRAY OF LAMBDA READ
		if(usedPSF){
			PRECISION * fInterpolated = calloc(nlambda,sizeof(PRECISION));
			interpolationPSF(deltaLambda,  PSF, fitsImage->pixels[0].vLambda ,CENTRAL_WL, N_SAMPLES_PSF,fInterpolated, nlambda);
		}

		//initializing weights
		PRECISION *w, *sig;
		weights_init(nlambda, sigma, weight, nweight, &w, &sig, noise);

		int indexPixel = 0;

		// ALLOCATE MEMORY FOR STORE THE RESULTS 

		vModels = calloc (fitsImage->numPixels , sizeof(Init_Model));
		vChisqrf = calloc (fitsImage->numPixels , sizeof(double));

		printf("\n START EXECUTION OF INVERSION ");
		printf("\n**********");
		t = clock();
		for(indexPixel = 0; indexPixel < fitsImage->numPixels; indexPixel++){
			

			//Initial Model
			Init_Model initModel;
			initModel.eta0 = INITIAL_MODEL.eta0;
			initModel.B = INITIAL_MODEL.B; //200 700
			initModel.gm = INITIAL_MODEL.gm;
			initModel.az = INITIAL_MODEL.az;
			initModel.vlos = INITIAL_MODEL.vlos; //km/s 0
			initModel.mac = INITIAL_MODEL.mac;
			initModel.dopp = INITIAL_MODEL.dopp;
			initModel.aa = INITIAL_MODEL.aa;
			initModel.alfa = INITIAL_MODEL.alfa; //0.38; //stray light factor
			initModel.S0 = INITIAL_MODEL.S0;
			initModel.S1 = INITIAL_MODEL.S1;

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

			int indexAux = 0;
			/*printf("\n Valores del spectro para el pixel %d ", indexPixel);
			for(indexAux = 0; indexAux < (nlambda*NPARMS);indexAux++){
				printf(" %fd ", fitsImage->pixels[indexPixel].spectro[indexAux]);
			}
			printf("\n *****************");*/
			/*printf("\n Valores del lambda para el pixel %d ", indexPixel);
			for(indexAux = 0; indexAux < (nlambda);indexAux++){
				printf(" %fd ", fitsImage->pixels[indexPixel].vLambda[indexAux]);
			}
			printf("\n *****************");*/
			//inversion
			if (CLASSICAL_ESTIMATES != 2)
			{
				//Se introduce en S0 el valor de Blos si solo se calculan estimaciones clásicas
				//Aqui se anula esa asignación porque se va a realizar la inversion RTE completa
				initModel.S0 = INITIAL_MODEL.S0;

				lm_mils(cuantic, wlines, nwlines, fitsImage->pixels[indexPixel].vLambda, fitsImage->pixels[indexPixel].nLambda, fitsImage->pixels[indexPixel].spectro, fitsImage->pixels[indexPixel].nLambda, &initModel, spectra, &chisqrf, &iter, slight, toplim, Max_iter,
						weight, nweight, fix, sig, filter, ilambda, noise, pol, getshi, 0,&INSTRUMENTAL_CONVOLUTION,&NMUESTRAS_G);
			}

			vModels[indexPixel] = initModel;
			vChisqrf[indexPixel] = chisqrf;

		}
		t = clock() - t;
		timeReadImage = ((double)t)/CLOCKS_PER_SEC; // in seconds 
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
				//double chisqr = vChisqrf[i];
				int NMODEL = 12; //Numero de parametros del modelo

				mil_sinrf(cuantic, &initModel, wlines, nwlines, fitsImage->pixels[i].vLambda, nlambda, spectra, AH, 0, filter);

				me_der(cuantic, &initModel, wlines, nwlines, fitsImage->pixels[i].vLambda, nlambda, d_spectra, AH, slight, 0, filter);
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
