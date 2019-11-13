
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
#include "readConfig.h"
#include <unistd.h>
#include <complex.h>
#include <fftw3.h> //siempre a continuacion de complex.h
#include "fftw.h"
#include <omp.h>
//#include "mkl_vsl.h"

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
PRECISION *interpolatedPSF;
int FGlobal, HGlobal, uuGlobal;

PRECISION *d_spectra, *spectra, *spectra_mac;


// GLOBAL variables to use for FFT calculation 
//VSLConvTaskPtr taskConv;
fftw_complex * inSpectraFwPSF, *inSpectraBwPSF, *outSpectraFwPSF, *outSpectraBwPSF;
fftw_complex * inSpectraFwMAC, *inSpectraBwMAC, *outSpectraFwMAC, *outSpectraBwMAC;
fftw_plan planForwardPSF, planBackwardPSF;
fftw_plan planForwardMAC, planBackwardMAC;
fftw_complex * inFilterMAC, * inFilterMAC_DERIV, * outFilterMAC, * outFilterMAC_DERIV;
fftw_plan planFilterMAC, planFilterMAC_DERIV;


fftw_complex * fftw_G_PSF, * fftw_G_MAC_PSF, * fftw_G_MAC_DERIV_PSF;
fftw_complex * inPSF_MAC, * inMulMacPSF, * inPSF_MAC_DERIV, *inMulMacPSFDeriv, *outConvFilters, * outConvFiltersDeriv;
fftw_plan planForwardPSF_MAC, planForwardPSF_MAC_DERIV,planBackwardPSF_MAC, planBackwardPSF_MAC_DERIV;


//Convolutions values
int sizeG = 0;
PRECISION FWHM = 0;
int KIND_CONVOLUTION = 0; // 0 --> Discreet convolution ; 1 --> Convolution with FFT
ConfigControl configCrontrolFile;

int main(int argc, char **argv)
{
	double *wlines;
	int nlambda;
	Init_Model *vModels;
	double chisqrf, * vChisqrf;
	
	int posCENTRAL_WL; // position central wl in file of LINES
	Init_Model INITIAL_MODEL;
	PRECISION * deltaLambda, * PSF;
	int N_SAMPLES_PSF;
	clock_t t_ini;
	

	//----------------------------------------------
	
	
	
	
	PRECISION * slight = NULL;
	int dimStrayLight;

	const char  * nameInputFileSpectra ;
	const char  * nameInputFileLambda ;
	const char  * nameOutputFileModels;
	const char  * nameOutputFilePerfiles;
	const char	* nameInputFileLines;
	const char 	* nameInputFileInitModel ;
	const char	* nameInputFilePSF ;	

   FitsImage * fitsImage;
	double  dat[7];

	/********************* Read data input from file ******************************/

	/* Read data input from file */

	loadInitialValues(&configCrontrolFile);
	readConfigControl(argv[1],&configCrontrolFile,1);

	nameInputFileSpectra = configCrontrolFile.ObservedProfiles;
	nameInputFileLambda = configCrontrolFile.WavelengthFile;
	nameInputFileLines = configCrontrolFile.AtomicParametersFile;
	nameInputFileInitModel = configCrontrolFile.InitialGuessModel;
	nameOutputFileModels = configCrontrolFile.OutputModelFile;
	nameOutputFilePerfiles = configCrontrolFile.OutputSynthesisFile;
	
	nameInputFilePSF = configCrontrolFile.PSFFile;
	FWHM = configCrontrolFile.FWHM;
	if(strcmp(configCrontrolFile.TypeConvolution, CONVOLUTION_FFT)==0)
		KIND_CONVOLUTION = 1;
	else if(strcmp(configCrontrolFile.TypeConvolution, CONVOLUTION_DIRECT)==0)
	{
		 KIND_CONVOLUTION = 0;
	}
	
	/*if(!readParametersFileInput(argv[1], &configCrontrolFile.NumberOfCycles,&CLASSICAL_ESTIMATES,&PRINT_SINTESIS,nameInputFileSpectra,nameInputFileLambda,nameInputFileLines,nameInputFileInitModel, &CENTRAL_WL, nameOutputFileModels,nameOutputFilePerfiles,&configCrontrolFile.ConvolveWithPSF,nameInputFilePSF,&FWHM,&KIND_CONVOLUTION)){
		printf("\n********************* EXITING THE PROGRAM . ERROR READING PARAMETERS FILE ****************************\n");
		return -1;
	}*/

	posCENTRAL_WL = readFileCuanticLines(nameInputFileLines,dat,configCrontrolFile.CentralWaveLenght,1);
	if(!posCENTRAL_WL){
		printf("\n CUANTIC LINE NOT FOUND, REVIEW IT. INPUT CENTRAL WAVE LENGHT: %f",configCrontrolFile.CentralWaveLenght);
		exit(1);
	}
	readInitialModel(&INITIAL_MODEL,nameInputFileInitModel);
	
	/*********************************************** INITIALIZE VARIABLES  *********************************/

	CC = PI / 180.0;
	CC_2 = CC * 2;

	wlines = (double *)calloc(2, sizeof(double));
	wlines[0] = 1;
	wlines[1] = configCrontrolFile.CentralWaveLenght;
			
	/******************* CREATE CUANTINC AND INITIALIZE DINAMYC MEMORY*******************/

	cuantic = create_cuantic(dat);
	InitializePointerShareCalculation();

	/****************************************************************************************************/

	// READ PIXELS FROM IMAGE 
	double timeReadImage, timeExecuteClassicalEstimates;
	clock_t t;
	t = clock();
	fitsImage = readFitsSpectroImage(nameInputFileSpectra);
	t = clock() - t;
	timeReadImage = ((double)t)/CLOCKS_PER_SEC; // in seconds 
	
	printf("\n\n TIME TO READ FITS IMAGE:  %f seconds to execute \n", timeReadImage); 
	

	if(fitsImage!=NULL && readFitsLambdaFile(nameInputFileLambda,fitsImage)){

		// check if read stray light
		if(strcmp(configCrontrolFile.StrayLightFile,"")){ //  IF NOT EMPTY READ stray light file 
			slight = readFitsStrayLightFile(configCrontrolFile.StrayLightFile,&dimStrayLight,fitsImage->nLambdas,fitsImage->rows,fitsImage->cols);
		}
		// INTERPOLATE PSF WITH ARRAY OF LAMBDA READ
		/****************************************************************************************************/	
		// if parameter name of psf has been apported we read the file, in other case create the gaussian with the parameters	
		//initializing weights
		PRECISION *w, *sig;
		// THE NUMBER OF LAMBDAS IS READ FROM INPUT FILES 
		nlambda = fitsImage->nLambdas;
		//***************************************** INIT MEMORY WITH SIZE OF LAMBDA ****************************************************//
		AllocateMemoryDerivedSynthesis(nlambda);
		weights_init(configCrontrolFile.sigma, &w, &sig, configCrontrolFile.noise);

		// ************************** DEFINE PLANS TO EXECUTE MACROTURBULENCE IF NECESSARY **********************************************//
		//fftw_init_threads();
		int edge=(nlambda%2)+1;
		int numln=nlambda-edge+1;
		// MACROTURBULENCE PLANS
		inFilterMAC = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * numln);
		outFilterMAC = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * numln);
		planFilterMAC = fftw_plan_dft_1d(numln, inFilterMAC, outFilterMAC, FFT_FORWARD, FFTW_EXHAUSTIVE);
		inFilterMAC_DERIV = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * numln);
		outFilterMAC_DERIV = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * numln);
		planFilterMAC_DERIV = fftw_plan_dft_1d(numln, inFilterMAC_DERIV, outFilterMAC_DERIV, FFT_FORWARD, FFTW_EXHAUSTIVE);


		inSpectraFwMAC = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * numln);
		outSpectraFwMAC = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * numln);
		planForwardMAC = fftw_plan_dft_1d(numln, inSpectraFwMAC, outSpectraFwMAC, FFT_FORWARD, FFTW_EXHAUSTIVE);
		inSpectraBwMAC = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * numln);
		outSpectraBwMAC = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * numln);		
		planBackwardMAC = fftw_plan_dft_1d(numln, inSpectraBwMAC, outSpectraBwMAC, FFT_BACKWARD, FFTW_EXHAUSTIVE);

		// ********************************************* IF PSF HAS BEEN SELECTEC IN TROL READ PSF FILE OR CREATE GAUSSIAN FILTER ***********//
		if(configCrontrolFile.ConvolveWithPSF){
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
						PRECISION * fInterpolated = calloc(nlambda,sizeof(PRECISION));
						interpolationLinearPSF(deltaLambda,  PSF, fitsImage->pixels[0].vLambda ,configCrontrolFile.CentralWaveLenght, N_SAMPLES_PSF,fInterpolated, nlambda);						
					}
					else{
						//G = vgauss(FWHM, NMUESTRAS_G, DELTA);
						//double * fgauss_WL(double FWHM, double step_between_lw, double lambda0, double lambdaCentral, int nLambda, int * sizeG)
						G = fgauss_WL(FWHM,fitsImage->pixels[0].vLambda[1]-fitsImage->pixels[0].vLambda[0],fitsImage->pixels[0].vLambda[0],fitsImage->pixels[0].vLambda[nlambda/2],fitsImage->pixels[0].nLambda,&sizeG);
					}
			}else{
				//G = vgauss(FWHM, NMUESTRAS_G, DELTA);
				G = fgauss_WL(FWHM,fitsImage->pixels[0].vLambda[1]-fitsImage->pixels[0].vLambda[0],fitsImage->pixels[0].vLambda[0],fitsImage->pixels[0].vLambda[nlambda/2],fitsImage->pixels[0].nLambda,&sizeG);
			}


			
			//PSF FILTER PLANS 
			inSpectraFwPSF = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * numln);
			outSpectraFwPSF = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * numln);
			planForwardPSF = fftw_plan_dft_1d(numln, inSpectraFwPSF, outSpectraFwPSF, FFT_FORWARD, FFTW_EXHAUSTIVE);
			inSpectraBwPSF = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * numln);
			outSpectraBwPSF = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * numln);		
			planBackwardPSF = fftw_plan_dft_1d(numln, inSpectraBwPSF, outSpectraBwPSF, FFT_BACKWARD, FFTW_EXHAUSTIVE);

			fftw_complex * in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * numln);
			int i;
			for (i = 0; i < numln; i++)
			{
				in[i] = G[i] + 0 * _Complex_I;
			}
			fftw_G_PSF = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * numln);
			fftw_plan p = fftw_plan_dft_1d(numln, in, fftw_G_PSF, FFT_FORWARD, FFTW_ESTIMATE);
			fftw_execute(p);
			for (i = 0; i < numln; i++)
			{
				fftw_G_PSF[i] = fftw_G_PSF[i] / numln;
			}
			fftw_destroy_plan(p);
			fftw_free(in);
			
			inPSF_MAC = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * numln);
			fftw_G_MAC_PSF = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * numln);
			planForwardPSF_MAC = fftw_plan_dft_1d(numln, inPSF_MAC, fftw_G_MAC_PSF, FFT_FORWARD, FFTW_EXHAUSTIVE);
			inMulMacPSF = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * numln);
			outConvFilters = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * numln);
			planBackwardPSF_MAC = fftw_plan_dft_1d(numln, inMulMacPSF, outConvFilters, FFT_BACKWARD, FFTW_EXHAUSTIVE);


			inPSF_MAC_DERIV = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * numln);
			fftw_G_MAC_DERIV_PSF = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * numln);
			planForwardPSF_MAC_DERIV = fftw_plan_dft_1d(numln, inPSF_MAC_DERIV, fftw_G_MAC_DERIV_PSF, FFT_FORWARD, FFTW_EXHAUSTIVE);
			inMulMacPSFDeriv = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * numln);
			outConvFiltersDeriv = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * numln);
			planBackwardPSF_MAC_DERIV = fftw_plan_dft_1d(numln, inMulMacPSFDeriv, outConvFiltersDeriv, FFT_BACKWARD, FFTW_EXHAUSTIVE);			

			
		}		

		int indexPixel = 0;

		// ALLOCATE MEMORY FOR STORE THE RESULTS 

		vModels = calloc (fitsImage->numPixels , sizeof(Init_Model));
		vChisqrf = calloc (fitsImage->numPixels , sizeof(double));

		t = clock();
		
		printf("\n***********************  PROGRESS INVERSION *******************************\n\n");

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
			

			if (configCrontrolFile.UseClassicalEstimates)
			{

				t_ini = clock();
				estimacionesClasicas(wlines[1], fitsImage->pixels[indexPixel].vLambda, fitsImage->pixels[indexPixel].nLambda, fitsImage->pixels[indexPixel].spectro, &initModel);
				
				timeExecuteClassicalEstimates += (clock() - t_ini);  
				

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
			if (configCrontrolFile.UseRTEInversion)
			{
				//Se introduce en S0 el valor de Blos si solo se calculan estimaciones clásicas
				//Aqui se anula esa asignación porque se va a realizar la inversion RTE completa
				initModel.S0 = INITIAL_MODEL.S0;
				PRECISION * slightPixel;
				if(slight==NULL) 
					slightPixel = NULL;
				else{
					if(dimStrayLight==nlambda) 
						slightPixel = slight;
					else 
						slightPixel = slight+nlambda*indexPixel;
				}
				lm_mils(cuantic, wlines, fitsImage->pixels[indexPixel].vLambda, fitsImage->pixels[indexPixel].nLambda, fitsImage->pixels[indexPixel].spectro, fitsImage->pixels[indexPixel].nLambda, &initModel, spectra, &chisqrf, slightPixel, configCrontrolFile.toplim, configCrontrolFile.NumberOfCycles,
						configCrontrolFile.WeightForStokes, configCrontrolFile.fix, sig, configCrontrolFile.InitialDiagonalElement,0,&configCrontrolFile.ConvolveWithPSF);						
			}

			vModels[indexPixel] = initModel;
			vChisqrf[indexPixel] = chisqrf;
			
			//printf ("\t\t %.2f seconds -- %.2f %%\r",  ((double)(clock() - t)/CLOCKS_PER_SEC) , ((indexPixel*100.)/fitsImage->numPixels));
		}
		t = clock() - t;

		printf("\n\n TIME EXECUTIN CLASSICAL ESTIMATES: %f seconds to execute \n", ((double)timeExecuteClassicalEstimates)/CLOCKS_PER_SEC);
		timeReadImage = ((double)t)/CLOCKS_PER_SEC; // in seconds 
		printf("\n FINISH EXECUTION OF INVERSION: %f seconds to execute \n", timeReadImage);
		printf("\n**********");
		if(!writeFitsImageModels(nameOutputFileModels,fitsImage->rows,fitsImage->cols,vModels,vChisqrf,configCrontrolFile.fix,configCrontrolFile.saveChisqr)){
				printf("\n ERROR WRITING FILE OF MODELS: %s",nameOutputFileModels);
		}

		// PROCESS FILE OF SYNTETIC PROFILES

		if(configCrontrolFile.SaveSynthesisProfile){
			FreeMemoryDerivedSynthesis();
			InitializePointerShareCalculation();
			AllocateMemoryDerivedSynthesis(nlambda);

			weights_init(configCrontrolFile.sigma,&w, &sig, configCrontrolFile.noise);
			int i;
			for( i=0;i<fitsImage->numPixels;i++)
			{

				Init_Model initModel = vModels[i];
				mil_sinrf(cuantic, &initModel, wlines, fitsImage->pixels[i].vLambda, nlambda, spectra, AH, 0,slight,NULL,configCrontrolFile.ConvolveWithPSF);
				spectral_synthesis_convolution(&nlambda);
				me_der(cuantic, &initModel, wlines, fitsImage->pixels[i].vLambda, nlambda, d_spectra, spectra, AH, slight, 0,1,configCrontrolFile.ConvolveWithPSF);
				response_functions_convolution(&nlambda);
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
		
		//vslConvDeleteTask(&taskConv);
		// FREE MACROTURBULENCE PLANS AND MEMORY
		fftw_free(inFilterMAC);
		fftw_free(outFilterMAC);
		fftw_destroy_plan(planFilterMAC);
		fftw_free(inFilterMAC_DERIV);
		fftw_free(outFilterMAC_DERIV);
		fftw_destroy_plan(planFilterMAC_DERIV);
		fftw_free(inSpectraFwMAC);
		fftw_free(outSpectraFwMAC);
		fftw_destroy_plan(planForwardMAC);
		fftw_free(inSpectraBwMAC);
		fftw_free(outSpectraBwMAC);
		fftw_destroy_plan(planBackwardMAC);

		if(configCrontrolFile.ConvolveWithPSF && KIND_CONVOLUTION==1){
			fftw_free(inSpectraFwPSF);
			fftw_free(outSpectraFwPSF);
			fftw_destroy_plan(planForwardPSF);
			fftw_free(inSpectraBwPSF);
			fftw_free(outSpectraBwPSF);
			fftw_destroy_plan(planBackwardPSF);
		}

		free(vModels);
		free(vChisqrf);
	}
	else{
		printf("\n\n ***************************** FITS FILE WITH THE SPECTRO IMAGE CAN NOT BE READ IT ******************************\n");
	}

	printf(" \n***********************  IMAGE INVERSION DONE, CLEANING MEMORY *********************\n");


	freeFitsImage(fitsImage);

	free(cuantic);
	free(wlines);
	FreeMemoryDerivedSynthesis();
	if(configCrontrolFile.ConvolveWithPSF)
		free(fftw_G_PSF);
	free(G);

	return 0;
}
