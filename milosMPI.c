
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
// RTE INVERSION C - MPI code for SOPHI (based on the ILD code MILOS by D. Orozco)
// Juanp && Manu (IAA-CSIC)


#include "mpi.h"
#include <time.h>
#include "defines.h"
//#include "nrutil.h"
#include <string.h>
#include <stdio.h>
#include <stddef.h>
#include "fitsio.h"
#include "utilsFits.h"
#include "readConfig.h"
#include "lib.h"
#include "milosUtils.h"
#include <unistd.h>
#include <complex.h>
#include <fftw3.h> //always after complex.h




// ***************************** FUNCTIONS TO READ FITS FILE *********************************************************


long long int c1, c2, cd, semi, c1a, c2a, cda; //variables of 64 bits to read whatch cycles 
long long int c1total, c2total, cdtotal, ctritotal;

Cuantic *cuantic; // Global variable with cuantic information 


PRECISION **PUNTEROS_CALCULOS_COMPARTIDOS;
int POSW_PUNTERO_CALCULOS_COMPARTIDOS;
int POSR_PUNTERO_CALCULOS_COMPARTIDOS;

PRECISION *gp1, *gp2, *dt, *dti, *gp3, *gp4, *gp5, *gp6, *etai_2;


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

fftw_complex * inSpectraFwPSF, *inSpectraBwPSF, *outSpectraFwPSF, *outSpectraBwPSF;
fftw_complex * inSpectraFwMAC, *inSpectraBwMAC, *outSpectraFwMAC, *outSpectraBwMAC;
fftw_plan planForwardPSF, planBackwardPSF;
fftw_plan planForwardMAC, planBackwardMAC;
fftw_complex * inFilterMAC, * inFilterMAC_DERIV, * outFilterMAC, * outFilterMAC_DERIV;
fftw_plan planFilterMAC, planFilterMAC_DERIV;
fftw_complex * fftw_G_PSF;

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
	int i, j;  // indexes 
	// INIT MPI  PROGRAM 
	int numProcs, idProc;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
   MPI_Comm_rank(MPI_COMM_WORLD, &idProc);
	MPI_Datatype mpiInitModel;

	const int root=0;	
		
	// FINISH STARTING PROGRAM 

	double *wlines;
	int nlambda, numPixels, indexPixel;


	
	int posCENTRAL_WL; // position central wl in file of LINES
	Init_Model INITIAL_MODEL;
	PRECISION * deltaLambda, * PSF;
	int N_SAMPLES_PSF;	
	
	

	// CONFIGURACION DE PARAMETROS A INVERTIR
	//INIT_MODEL=[eta0,magnet,vlos,landadopp,aa,gamma,azi,B1,B2,macro,alfa]
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


   FitsImage * fitsImage = NULL;
	double  dat[7];

	/********************* Read data input from file ******************************/

	loadInitialValues(&configCrontrolFile);
	readConfigControl(argv[1],&configCrontrolFile,(idProc==root));

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

	/*if(!readParametersFileInput(argv[1], &Max_iter,&CLASSICAL_ESTIMATES,&PRINT_SINTESIS,nameInputFileSpectra,nameInputFileLambda,nameInputFileLines,nameInputFileInitModel, &CENTRAL_WL, nameOutputFileModels,nameOutputFilePerfiles,&INSTRUMENTAL_CONVOLUTION,nameInputFilePSF,&FWHM,&KIND_CONVOLUTION)){
		printf("\n********************* EXITING THE PROGRAM . ERROR READING PARAMETERS FILE ****************************\n");
		return -1;
	}*/
	
	posCENTRAL_WL = readFileCuanticLines(nameInputFileLines,dat,configCrontrolFile.CentralWaveLenght,(idProc==root));
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

	
	numPixels=0;
	
	/******************* APPLY GAUSSIAN, CREATE CUANTINC AND INITIALIZE DINAMYC MEMORY*******************/

	cuantic = create_cuantic(dat);
	InitializePointerShareCalculation();

	/****************************************************************************************************/


	// ROOT PROCESS READ IMAGE FROM FILE 
	if(idProc==root){
		double timeReadImage;
		clock_t t;
		t = clock();
		fitsImage = readFitsSpectroImage(nameInputFileSpectra);
		t = clock() - t;
		timeReadImage = ((double)t)/CLOCKS_PER_SEC; // in seconds 
		printf("\n TIME TO READ FITS IMAGE:  %f seconds to execute \n", timeReadImage); 
		//READ LAMBDA VALUES FROM FITS FILE
		t = clock();
		if(fitsImage!=NULL && readFitsLambdaFile(nameInputFileLambda,fitsImage)){
			nlambda = fitsImage->nLambdas;
			//NLAMBDA = nlambda;
			numPixels = fitsImage->numPixels;
		}
		else{ /* EXIT IF THE IMAGE HAS NOT BEEN READ CORRECTLY */
			printf("**** ERROR READING FITS IMAGE ********************");
			return -1; 
		}
		t = clock() - t;
		timeReadImage = ((double)t)/CLOCKS_PER_SEC; // in seconds 
		printf("\n TIME TO READ FITS IMAGE LAMBDA:  %f seconds to execute \n", timeReadImage); 
				// check if read stray light
		if(strcmp(configCrontrolFile.StrayLightFile,"")){ //  IF NOT EMPTY READ stray light file 
			slight = readFitsStrayLightFile(configCrontrolFile.StrayLightFile,&dimStrayLight,fitsImage->nLambdas,fitsImage->rows,fitsImage->cols);
		}
	}
	// TO-DO send fInterpolated by broadcast 	
	/****************************************************************************************************/	
	// if parameter name of psf has been apported we read the file, in other case create the gaussian with the parameters
	
	if(configCrontrolFile.ConvolveWithPSF && idProc == root){
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
					G = fgauss_WL(FWHM,fitsImage->pixels[0].vLambda[1]-fitsImage->pixels[0].vLambda[0],fitsImage->pixels[0].vLambda[0],fitsImage->pixels[0].vLambda[nlambda/2],fitsImage->pixels[0].nLambda,&sizeG);
				}
		}else
			//G = vgauss(FWHM, NMUESTRAS_G, DELTA);
			G = fgauss_WL(FWHM,fitsImage->pixels[0].vLambda[1]-fitsImage->pixels[0].vLambda[0],fitsImage->pixels[0].vLambda[0],fitsImage->pixels[0].vLambda[nlambda/2],fitsImage->pixels[0].nLambda,&sizeG);
	}


	MPI_Barrier(MPI_COMM_WORLD); // Wait UNTIL THE IMAGE HAS BEEN READED COMPLETELY
	//  BROADCAST THE NUMBER OF LAMBDAS READS FROM THE FILE AND THE NUMBER OF PIXELS
	MPI_Bcast(&nlambda, 1, MPI_INT, root , MPI_COMM_WORLD);
	MPI_Bcast(&numPixels, 1, MPI_INT, root , MPI_COMM_WORLD);
	MPI_Bcast(&sizeG, 1, MPI_INT, root , MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	if(idProc!=root)
		G = calloc(sizeG,sizeof(double));
	MPI_Bcast(G, sizeG, MPI_DOUBLE, root , MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD); // WAIT UNTIL G HAS BEEN READ
	
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

	if(configCrontrolFile.ConvolveWithPSF){
		inSpectraFwPSF = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * numln);
		outSpectraFwPSF = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * numln);
		planForwardPSF = fftw_plan_dft_1d(numln, inSpectraFwPSF, outSpectraFwPSF, FFT_FORWARD, FFTW_EXHAUSTIVE);
		inSpectraBwPSF = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * numln);
		outSpectraBwPSF = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * numln);		
		planBackwardPSF = fftw_plan_dft_1d(numln, inSpectraBwPSF, outSpectraBwPSF, FFT_BACKWARD, FFTW_EXHAUSTIVE);

		// CALCULATE FFT OF GAUSSIAN FILTER 
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

	/*if(idProc == root){
		printf("\n\n VALORES GAUSSIANA: \n");
		int i;
		for(i=0;i<sizeG;i++){
			printf("%lf\n", G[i]);
		}
	}
	printf("\n*********\n");*/

	//vsldConvNewTask1D(&taskConv,VSL_CONV_MODE_AUTO,nlambda,sizeG,nlambda);
	// IF THE NUMBER OF PIXELS IS NOT GREATER THAN 0 WE DON'T CONITUNUE 

	if(numPixels > 0){
		/*printf("\n************************************************");		
		printf("\n************************************************");		
		printf("\n ESTOY EN EL PROCESO %d . Número de pixeles leidos: %d", idProc, numPixels);
		printf("\n************************************************");*/
		double local_start, local_finish, local_elapsed, elapsed;
		double local_start_execution, local_finish_execution, local_elapsed_execution, elapsed_execution;
		double local_start_scatter, local_finish_scatter, local_elapsed_scatter, elapsed_scatter;
		double local_start_gather, local_finish_gather, local_elapsed_gather, elapsed_gather;
		
		Init_Model * resultsInitModel;
		Init_Model * resultsInitModelTotal;
		double * chisqrfTotal, * chisqrf;
		if(idProc == root){
			resultsInitModelTotal = calloc (numPixels , sizeof(Init_Model));
			chisqrfTotal = calloc (numPixels , sizeof(double));
		}

		// allocate memory in all processes 
		AllocateMemoryDerivedSynthesis(nlambda);
		//initializing weights
		PRECISION *w, *sig;
		weights_init(configCrontrolFile.sigma, &w, &sig, configCrontrolFile.noise);
		MPI_Request mpiRequestSpectro, mpiRequestLambda;
		

	   const int nitemsStructInitModel = 11;
		int blocklenghtInitModel [11] = {1,1,1,1,1,1,1,1,1,1,1};
		MPI_Datatype typesInitModel [11] = {MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE};
			
		MPI_Aint offsetsInitModel [11];
		offsetsInitModel[0] = offsetof(Init_Model, eta0);
		offsetsInitModel[1] = offsetof(Init_Model, B);
		offsetsInitModel[2] = offsetof(Init_Model, vlos);
		offsetsInitModel[3] = offsetof(Init_Model, dopp);
		offsetsInitModel[4] = offsetof(Init_Model, aa);
		offsetsInitModel[5] = offsetof(Init_Model, gm);
		offsetsInitModel[6] = offsetof(Init_Model, az);
		offsetsInitModel[7] = offsetof(Init_Model, S0);
		offsetsInitModel[8] = offsetof(Init_Model, S1);
		offsetsInitModel[9] = offsetof(Init_Model, mac);
		offsetsInitModel[10] = offsetof(Init_Model, alfa);

		
		MPI_Type_create_struct(nitemsStructInitModel, blocklenghtInitModel, offsetsInitModel, typesInitModel, &mpiInitModel);
		MPI_Type_commit(&mpiInitModel);

		int numPixelsProceso = numPixels/numProcs;
		int resto = numPixels % numProcs;
		int sum = 0;                // Sum of counts. Used to calculate displacements
		int sumSpectro = 0;
		int sumLambda = 0;
		
		
		
		int * sendcountsPixels = calloc(numProcs, sizeof(int)); // array describing how many elements to send to each process
		int * sendcountsSpectro = calloc(numProcs, sizeof(int));
		int * sendcountsLambda = calloc(numProcs, sizeof(int));
		int * displsPixels = calloc(numProcs, sizeof(int));  // array describing the displacements where each segment begins
		int * displsSpectro = calloc(numProcs, sizeof(int));
		int * displsLambda = calloc(numProcs, sizeof(int));

		for ( i = 0; i < numProcs; i++) {
			sendcountsPixels[i] = numPixelsProceso;
			if (resto > 0) {
					sendcountsPixels[i]++;
					resto--;
			}
			sendcountsSpectro[i] = sendcountsPixels[i]*nlambda*NPARMS;
			sendcountsLambda[i] = sendcountsPixels[i]*nlambda;
			displsPixels[i] = sum;
			displsSpectro[i] = sumSpectro;
			displsLambda[i] = sumLambda;
			sum += sendcountsPixels[i];
			sumSpectro += sendcountsSpectro[i];
			sumLambda += sendcountsLambda[i];
		}

		MPI_Barrier(MPI_COMM_WORLD); // Wait until all processes have their vlambda
		local_start = MPI_Wtime();

		/*if (root == idProc) {
			for ( i = 0; i < numProcs; i++) {
				printf("\n sendcounts[%d] sendcountsSpectro[%d] sendcountsLambda[%d] \n", sendcountsPixels[i], sendcountsSpectro[i] , sendcountsLambda[i] );
				printf("\n displcounts[%d] displSpectro[%d] displLambda[%d] \n", displsPixels[i], displsSpectro[i] , displsLambda[i] );
				printf("\n**********************************************");
			}
		}*/

		// SCATTER VPIXELS 
		PRECISION  * vSpectraSplit = calloc(sendcountsSpectro[idProc],sizeof(PRECISION));
		PRECISION  * vLambdaSplit = calloc(sendcountsLambda[idProc],sizeof(PRECISION));

		
		local_start_scatter = MPI_Wtime();
		
		if( root == idProc){
			MPI_Iscatterv(fitsImage->spectroImagen, sendcountsSpectro, displsSpectro, MPI_DOUBLE, vSpectraSplit, sendcountsSpectro[idProc], MPI_DOUBLE, root, MPI_COMM_WORLD,&mpiRequestSpectro);
			MPI_Iscatterv(fitsImage->vLambdaImagen, sendcountsLambda, displsLambda, MPI_DOUBLE,vLambdaSplit, sendcountsLambda[idProc], MPI_DOUBLE, root, MPI_COMM_WORLD,&mpiRequestLambda);
			MPI_Wait(&mpiRequestSpectro, MPI_STATUS_IGNORE);
			MPI_Wait(&mpiRequestLambda, MPI_STATUS_IGNORE);
		}
		else{
			MPI_Iscatterv(NULL, NULL,NULL, MPI_DOUBLE, vSpectraSplit, sendcountsSpectro[idProc], MPI_DOUBLE, root, MPI_COMM_WORLD,&mpiRequestSpectro);
			MPI_Iscatterv(NULL, NULL,NULL, MPI_DOUBLE, vLambdaSplit, sendcountsLambda[idProc], MPI_DOUBLE, root, MPI_COMM_WORLD,&mpiRequestLambda);
			MPI_Wait(&mpiRequestSpectro, MPI_STATUS_IGNORE);
			MPI_Wait(&mpiRequestLambda, MPI_STATUS_IGNORE);
		}		
		local_finish_scatter = MPI_Wtime();

		resultsInitModel = calloc(sendcountsPixels[idProc], sizeof(Init_Model));
		chisqrf = calloc(sendcountsPixels[idProc], sizeof(double));
		local_start_execution = MPI_Wtime();
		for(indexPixel = 0; indexPixel < sendcountsPixels[idProc]; indexPixel++){
			// get spectro and lambdas for the current pixel
			PRECISION * spectroPixel = calloc(nlambda*NPARMS,sizeof(PRECISION));
			PRECISION * lambdaPixel = calloc(nlambda,sizeof(PRECISION));

			for( j=0;j<(nlambda*NPARMS);j++){
				spectroPixel[j] = vSpectraSplit[(indexPixel*(nlambda*NPARMS))+j];
			}
		
			for( j=0;j<(nlambda);j++){
				lambdaPixel[j] = vLambdaSplit[(indexPixel*(nlambda))+j];
			}

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

			double auxChisqrf;

			if (configCrontrolFile.UseClassicalEstimates)
			{
				estimacionesClasicas(wlines[1], lambdaPixel, nlambda, spectroPixel, &initModel);
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
				/*lm_mils(cuantic, wlines, nwlines, lambdaPixel, nlambda, spectroPixel, nlambda, &initModel, spectra, &auxChisqrf, NULL, toplim, miter,
						weight, fix, sig, ilambda,0,&INSTRUMENTAL_CONVOLUTION,&sizeG,&KIND_CONVOLUTION);*/
				lm_mils(cuantic, wlines, lambdaPixel, nlambda, spectroPixel, nlambda, &initModel, spectra, &auxChisqrf, slightPixel, configCrontrolFile.toplim, configCrontrolFile.NumberOfCycles,
						configCrontrolFile.WeightForStokes, configCrontrolFile.fix, sig, configCrontrolFile.InitialDiagonalElement, 0,&configCrontrolFile.ConvolveWithPSF);												
			}
			free(spectroPixel);
			free(lambdaPixel);
			resultsInitModel[indexPixel] = initModel;
			chisqrf[indexPixel] = auxChisqrf;
		}
		local_finish_execution = MPI_Wtime();

		local_start_gather = MPI_Wtime();
		/*if(idProc==root){
			MPI_Gatherv(resultsInitModel, sendcountsPixels[idProc], mpiInitModel, resultsInitModelTotal, sendcountsPixels, displsPixels, mpiInitModel, root, MPI_COMM_WORLD);
			MPI_Gatherv(chisqrf, sendcountsPixels[idProc], MPI_DOUBLE, chisqrfTotal, sendcountsPixels, displsPixels, MPI_DOUBLE, root, MPI_COMM_WORLD);
		}
		else{
			MPI_Gatherv(resultsInitModel, sendcountsPixels[idProc], mpiInitModel, NULL, NULL, NULL, mpiInitModel, root, MPI_COMM_WORLD);
			MPI_Gatherv(chisqrf, sendcountsPixels[idProc], MPI_DOUBLE, NULL, NULL, NULL, MPI_DOUBLE, root, MPI_COMM_WORLD);
		}*/

		MPI_Igatherv(resultsInitModel, sendcountsPixels[idProc], mpiInitModel, resultsInitModelTotal, sendcountsPixels, displsPixels, mpiInitModel, root, MPI_COMM_WORLD,&mpiRequestSpectro);
		MPI_Igatherv(chisqrf, sendcountsPixels[idProc], MPI_DOUBLE, chisqrfTotal, sendcountsPixels, displsPixels, MPI_DOUBLE, root, MPI_COMM_WORLD,&mpiRequestLambda);		
		MPI_Wait(&mpiRequestSpectro, MPI_STATUS_IGNORE);
		MPI_Wait(&mpiRequestLambda, MPI_STATUS_IGNORE);		
		local_finish_gather = MPI_Wtime();
		
	
		//printf("\n PIXEL CALCULATION FINISHED IN PROCESS: %d \n", idProc);
		//MPI_Barrier(MPI_COMM_WORLD);  // WAIT UNTIL RECEIVED ALL INFORMATION
		local_finish = MPI_Wtime();
		local_elapsed = local_finish - local_start;
		local_elapsed_execution = local_finish_execution - local_start_execution;
		local_elapsed_scatter = local_finish_scatter - local_start_scatter;
		local_elapsed_gather = local_finish_gather - local_start_gather;
		MPI_Reduce(&local_elapsed, &elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		MPI_Reduce(&local_elapsed_execution, &elapsed_execution, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		MPI_Reduce(&local_elapsed_scatter, &elapsed_scatter, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		MPI_Reduce(&local_elapsed_gather, &elapsed_gather, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

		if(idProc==root){		 
			printf("\n Elapsed SCATTER time = %lf seconds\n", elapsed_scatter);
			printf("\n***********************************\n");
			printf("\n Elapsed GATHER time = %lf seconds\n", elapsed_gather);
			printf("\n***********************************\n");			
			printf("\n Elapsed EXECUTION time = %lf seconds\n", elapsed_execution);
			printf("\n***********************************\n");
			printf("\n Elapsed TOTAL time = %lf seconds\n", elapsed);
			printf("\n***********************************\n");
			double timeWriteImage;
			clock_t t;
			t = clock();
			if(!writeFitsImageModels(nameOutputFileModels,fitsImage->rows,fitsImage->cols,resultsInitModelTotal,chisqrfTotal,configCrontrolFile.fix,configCrontrolFile.saveChisqr)){
					printf("\n ERROR WRITING FILE OF MODELS: %s",nameOutputFileModels);
			}
			t = clock() - t;
			timeWriteImage = ((double)t)/CLOCKS_PER_SEC; // in seconds 
			printf("\n TIME TO WRITE FITS IMAGE:  %f seconds to execute \n", timeWriteImage); 
			// PROCESS FILE OF SYNTETIC PROFILES
			if(configCrontrolFile.SaveSynthesisProfile){
				for(indexPixel=0;indexPixel<fitsImage->numPixels;indexPixel++)
				{

					Init_Model initModel = resultsInitModelTotal[indexPixel];
					//double chisqr = vChisqrf[i];

					mil_sinrf(cuantic, &initModel, wlines, fitsImage->pixels[indexPixel].vLambda, nlambda, spectra, AH, 0,slight,NULL,configCrontrolFile.ConvolveWithPSF);
					spectral_synthesis_convolution(&nlambda);
					me_der(cuantic, &initModel, wlines, fitsImage->pixels[indexPixel].vLambda, nlambda, d_spectra,spectra, AH, NULL, 0 ,1,configCrontrolFile.ConvolveWithPSF);
					response_functions_convolution(&nlambda);
					int kk;
					for (kk = 0; kk < (nlambda * NPARMS); kk++)
					{
						fitsImage->pixels[indexPixel].spectro[kk] = spectra[kk] ;
					}
				}
				// WRITE SINTHETIC PROFILES TO FITS FILE
				if(!writeFitsImageProfiles(nameOutputFilePerfiles,nameInputFileSpectra,fitsImage)){
					printf("\n ERROR WRITING FILE OF SINTHETIC PROFILES: %s",nameOutputFilePerfiles);
				}
			}

			//vslConvDeleteTask(&taskConv);
			free(resultsInitModelTotal);		
			free(chisqrfTotal);
		}
		else{
			free(vSpectraSplit);
			free(vLambdaSplit);
			free(resultsInitModel);				
			free(chisqrf);
		}
		free(sendcountsPixels); // array describing how many elements to send to each process
		free(sendcountsSpectro);
		free(sendcountsLambda);
		free(displsPixels);  // array describing the displacements where each segment begins
		free(displsSpectro);
		free(displsLambda);

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


	}
	else{
		printf("\n\n ***************************** FITS FILE CAN NOT BE READ IT ******************************");
	}
	 if(idProc==root){
		freeFitsImage(fitsImage);
	}
	//printf("\n\n TOTAL sec : %.16g segundos\n", total_secs);
	free(cuantic);
	free(wlines);
	FreeMemoryDerivedSynthesis();
	// FREE TYPE OF MPI
	MPI_Type_free(&mpiInitModel);
	MPI_Finalize() ;
	free(G);
	return 0;
}

