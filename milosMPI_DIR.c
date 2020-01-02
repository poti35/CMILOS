
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
#include <sys/types.h>
#include <dirent.h>
#include <sys/stat.h>



// ***************************** FUNCTIONS TO READ FITS FILE *********************************************************


long long int c1, c2, cd, semi, c1a, c2a, cda; //variables of 64 bits to read whatch cycles 
long long int c1total, c2total, cdtotal, ctritotal;

Cuantic *cuantic; // Global variable with cuantic information 


PRECISION **PUNTEROS_CALCULOS_COMPARTIDOS;
int POSW_PUNTERO_CALCULOS_COMPARTIDOS;
int POSR_PUNTERO_CALCULOS_COMPARTIDOS;

PRECISION *dtaux, *etai_gp3, *ext1, *ext2, *ext3, *ext4;
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
PRECISION *G,*GMAC;
PRECISION *interpolatedPSF;
PRECISION AP[NTERMS*NTERMS*NPARMS],BT[NPARMS*NTERMS];
PRECISION * opa;
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

ConfigControl configCrontrolFile;

_Complex PRECISION *z,* zden, * zdiv;

int main(int argc, char **argv)
{
	int i, j;  // indexes 
	// INIT MPI  PROGRAM 
	int numProcs, idProc;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
   MPI_Comm_rank(MPI_COMM_WORLD, &idProc);

	/** DEFINE MPI TYPE TO SEND MODELS **/
	MPI_Datatype mpiInitModel;
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


	MPI_Datatype mpiName;
	const int nItemsStructName = 1;
	int blocklenghName [1] = {PATH_MAX};
	MPI_Datatype typesName [1] = {MPI_CHAR};
	MPI_Aint offsetName [1];
	offsetName[0] = offsetof(nameFile,name) ;
	MPI_Type_create_struct(nItemsStructName,blocklenghName,offsetName,typesName,&mpiName);
	MPI_Type_commit(&mpiName);

	MPI_Request mpiRequestSpectro, mpiRequestLambda,mpiRequestName;

	/************************************/
	const int root=0;	
	
	// FINISH STARTING PROGRAM 

	PRECISION *wlines;
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
	int numberOfFileSpectra;

   nameFile * vInputFileSpectra;
	nameFile * vInputFileSpectraParalell = NULL;
	nameFile * vOutputNameModels;
	nameFile * vOutputNameModelsParalell = NULL;
	nameFile * vInputFileSpectraLocal;
	nameFile * vOutputNameModelsLocal;

	const char  * nameInputFileSpectra ;
	const char  * nameInputFileLambda ;
	const char  * nameOutputFileModels;
	const char  * nameOutputFilePerfiles;
	const char	* nameInputFileLines;
	
	const char	* nameInputFilePSF ;


   FitsImage * fitsImage = NULL;
	PRECISION  dat[7];

	double local_start, local_finish, local_elapsed, elapsed;
	double local_start_execution, local_finish_execution, local_elapsed_execution, elapsed_execution;
	double local_start_scatter, local_finish_scatter, local_elapsed_scatter, elapsed_scatter;
	double local_start_gather, local_finish_gather, local_elapsed_gather, elapsed_gather;
	
	Init_Model * resultsInitModel;
	Init_Model * resultsInitModelTotal;
	PRECISION * chisqrfTotal, *vChisqrf, chisqrf;
	int * vNumIter, * vNumIterTotal; // to store the number of iterations used to converge for each pixel
	PRECISION  * vSpectraSplit, * vLambdaSplit;

	int sendcountsPixels [numProcs] ; // array describing how many elements to send to each process
	int sendcountsSpectro [numProcs];
	int sendcountsLambda [numProcs];
	int sendcountsNameInputFiles [numProcs];  // how many files per process
	int displsPixels [numProcs];  // array describing the displacements where each segment begins
	int displsSpectro [numProcs];
	int displsLambda [numProcs];
	int displsNameInputFiles [numProcs]; // how many 




	/********************* Read data input from file ******************************/

	loadInitialValues(&configCrontrolFile);
	readParametersFileInput(argv[1],&configCrontrolFile,(idProc==root));

	nameInputFileSpectra = configCrontrolFile.ObservedProfiles;
	nameInputFileLambda = configCrontrolFile.WavelengthFile;
	nameInputFileLines = configCrontrolFile.AtomicParametersFile;
	
	nameOutputFileModels = configCrontrolFile.OutputModelFile;  // IN THIS CASE MUST BE A DIRECTORY 
	if(!isDirectory(nameOutputFileModels)){
		printf("\n ERROR: Name Output Models must be a directory: %s",nameOutputFileModels);
		exit(EXIT_FAILURE);
	}
	nameOutputFilePerfiles = configCrontrolFile.OutputSynthesisFile;
	
	nameInputFilePSF = configCrontrolFile.PSFFile;
	FWHM = configCrontrolFile.FWHM;
	

	posCENTRAL_WL = readFileCuanticLines(nameInputFileLines,dat,configCrontrolFile.CentralWaveLenght,(idProc==root));
	if(!posCENTRAL_WL){
		printf("\n CUANTIC LINE NOT FOUND, REVIEW IT. INPUT CENTRAL WAVE LENGHT: %f",configCrontrolFile.CentralWaveLenght);
		exit(1);
	}
	readInitialModel(&INITIAL_MODEL,configCrontrolFile.InitialGuessModel);

	/*********************************************** INITIALIZE VARIABLES  *********************************/

	CC = PI / 180.0;
	CC_2 = CC * 2;
	
	wlines = (PRECISION *)calloc(2, sizeof(PRECISION));
	wlines[0] = 1;
	wlines[1] = configCrontrolFile.CentralWaveLenght;

	numPixels=0;	
	/******************* APPLY GAUSSIAN, CREATE CUANTINC AND INITIALIZE DINAMYC MEMORY*******************/

	cuantic = create_cuantic(dat,(idProc==root));
	
	/****************************************************************************************************/


	/**************************************** READ FITS LAMBDA EACH PROCESS AND STRAY LIGHT ******************************/
	PRECISION * vGlobalLambda = calloc(configCrontrolFile.nliobs, sizeof(PRECISION));
	int readOK = readFitsLambdaToArray(configCrontrolFile.WavelengthFile,0,0,configCrontrolFile.nliobs,vGlobalLambda);

	if(!readOK){
		printf("\n FILE WITH WAVELENGHT HAS NOT BEEN READ PROPERLY, please check it.\n");
		free(vGlobalLambda);
		exit(EXIT_FAILURE);
	}

	if(access(configCrontrolFile.StrayLightFile,F_OK)!=-1){ //  IF NOT EMPTY READ stray light file 
		slight = readFitsStrayLightFile(configCrontrolFile.StrayLightFile,&dimStrayLight,configCrontrolFile.nliobs);
	}
	
	/*****************************************************************************************************/

	/************************************ CREARE GAUSSIAN AND RESERVE MEMORY *****************************/

	
	// ************************** DEFINE PLANS TO EXECUTE MACROTURBULENCE IF NECESSARY **********************************************//
	// MACROTURBULENCE PLANS
	nlambda = configCrontrolFile.nliobs;
	inFilterMAC = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nlambda);
	outFilterMAC = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nlambda);
	planFilterMAC = fftw_plan_dft_1d(nlambda, inFilterMAC, outFilterMAC, FFT_FORWARD, FFTW_EXHAUSTIVE);
	inFilterMAC_DERIV = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nlambda);
	outFilterMAC_DERIV = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nlambda);
	planFilterMAC_DERIV = fftw_plan_dft_1d(nlambda, inFilterMAC_DERIV, outFilterMAC_DERIV, FFT_FORWARD, FFTW_EXHAUSTIVE);


	inSpectraFwMAC = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nlambda);
	outSpectraFwMAC = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nlambda);
	planForwardMAC = fftw_plan_dft_1d(nlambda, inSpectraFwMAC, outSpectraFwMAC, FFT_FORWARD, FFTW_EXHAUSTIVE);
	inSpectraBwMAC = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nlambda);
	outSpectraBwMAC = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nlambda);		
	planBackwardMAC = fftw_plan_dft_1d(nlambda, inSpectraBwMAC, outSpectraBwMAC, FFT_BACKWARD, FFTW_EXHAUSTIVE);

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
					//slog_error(0,"File \"%s\" does not exist!!!\n",nameInputFilePSF);
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
					interpolationLinearPSF(deltaLambda,  PSF, vGlobalLambda ,configCrontrolFile.CentralWaveLenght, N_SAMPLES_PSF,fInterpolated, nlambda);						
				}
				else{
					//G = vgauss(FWHM, NMUESTRAS_G, DELTA);
					//PRECISION * fgauss_WL(PRECISION FWHM, PRECISION step_between_lw, PRECISION lambda0, PRECISION lambdaCentral, int nLambda, int * sizeG)
					G = fgauss_WL(FWHM,vGlobalLambda[1]-vGlobalLambda[0],vGlobalLambda[0],vGlobalLambda[nlambda/2],nlambda,&sizeG);
				}
		}else{
			//G = vgauss(FWHM, NMUESTRAS_G, DELTA);
			G = fgauss_WL(FWHM,vGlobalLambda[1]-vGlobalLambda[0],vGlobalLambda[0],vGlobalLambda[nlambda/2],nlambda,&sizeG);
		}


		
		//PSF FILTER PLANS 
		inSpectraFwPSF = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nlambda);
		outSpectraFwPSF = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nlambda);
		planForwardPSF = fftw_plan_dft_1d(nlambda, inSpectraFwPSF, outSpectraFwPSF, FFT_FORWARD, FFTW_EXHAUSTIVE);
		inSpectraBwPSF = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nlambda);
		outSpectraBwPSF = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nlambda);		
		planBackwardPSF = fftw_plan_dft_1d(nlambda, inSpectraBwPSF, outSpectraBwPSF, FFT_BACKWARD, FFTW_EXHAUSTIVE);

		fftw_complex * in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * nlambda);
		int i;
		for (i = 0; i < nlambda; i++)
		{
			in[i] = G[i] + 0 * _Complex_I;
		}
		fftw_G_PSF = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * nlambda);
		fftw_plan p = fftw_plan_dft_1d(nlambda, in, fftw_G_PSF, FFT_FORWARD, FFTW_ESTIMATE);
		fftw_execute(p);
		for (i = 0; i < nlambda; i++)
		{
			fftw_G_PSF[i] = fftw_G_PSF[i] / nlambda;
		}
		fftw_destroy_plan(p);
		fftw_free(in);
		
		inPSF_MAC = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nlambda);
		fftw_G_MAC_PSF = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nlambda);
		planForwardPSF_MAC = fftw_plan_dft_1d(nlambda, inPSF_MAC, fftw_G_MAC_PSF, FFT_FORWARD, FFTW_EXHAUSTIVE);
		inMulMacPSF = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nlambda);
		outConvFilters = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nlambda);
		planBackwardPSF_MAC = fftw_plan_dft_1d(nlambda, inMulMacPSF, outConvFilters, FFT_BACKWARD, FFTW_EXHAUSTIVE);


		inPSF_MAC_DERIV = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nlambda);
		fftw_G_MAC_DERIV_PSF = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nlambda);
		planForwardPSF_MAC_DERIV = fftw_plan_dft_1d(nlambda, inPSF_MAC_DERIV, fftw_G_MAC_DERIV_PSF, FFT_FORWARD, FFTW_EXHAUSTIVE);
		inMulMacPSFDeriv = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nlambda);
		outConvFiltersDeriv = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nlambda);
		planBackwardPSF_MAC_DERIV = fftw_plan_dft_1d(nlambda, inMulMacPSFDeriv, outConvFiltersDeriv, FFT_BACKWARD, FFTW_EXHAUSTIVE);			

	}

	/*****************************************************************************************************/

	// ROOT PROCESS READ IMAGE FROM FILE 
	if(idProc==root){
		clock_t t;
		t = clock();
      	// CHECK IF INPUT OBSERVED PROFILES COMES IN A DIRECTORY OR IS A FILE

      numberOfFileSpectra = 0;
      if(isDirectory(nameInputFileSpectra)){
         DIR * folder = opendir(nameInputFileSpectra);
         
         struct dirent *entry;
         struct stat filestat;
         if(folder == NULL)
         {
            perror("Unable to read directory");
            return(1);
         }
         while ((entry=readdir(folder)) )
         {
            stat(entry->d_name,&filestat);
            if( !isDirectory(entry->d_name) && strcmp(entry->d_name,".")!=0 && strcmp(entry->d_name,"..")!=0 ){
               numberOfFileSpectra++;
            }
         }
         
         closedir(folder);
         vInputFileSpectra = (nameFile *)malloc(numberOfFileSpectra*sizeof(nameFile));
			vOutputNameModels = (nameFile *) malloc(numberOfFileSpectra*sizeof(nameFile));
         int i;

         
         folder = opendir(nameInputFileSpectra);
         int nameIndex = 0;
         while ((entry=readdir(folder)) )
         {
            stat(entry->d_name,&filestat);
            
            
            if( !isDirectory(entry->d_name) && strcmp(entry->d_name,".")!=0 && strcmp(entry->d_name,"..")!=0 ){
					// IN 
               strcpy( vInputFileSpectra[nameIndex].name, nameInputFileSpectra);
               if( nameInputFileSpectra[strlen(nameInputFileSpectra)-1]!='/')
                  strcat( vInputFileSpectra[nameIndex].name, "/");
               strcat( vInputFileSpectra[nameIndex].name, entry->d_name);
					// OUT 
					strcpy( vOutputNameModels[nameIndex].name, nameOutputFileModels);
					if( nameOutputFileModels[strlen(nameOutputFileModels)-1]!='/')
                  strcat( vOutputNameModels[nameIndex].name, "/");
					strcat(vOutputNameModels[nameIndex].name, configCrontrolFile.outputPrefix);
               strcat( vOutputNameModels[nameIndex].name, entry->d_name);
					/*strip_ext(vOutputNameModels[nameIndex].name);
					strcat (vOutputNameModels[nameIndex].name, "_result_inversion.fits");*/
               nameIndex++;
            }
         }
         closedir(folder);
      }
      else{
         vInputFileSpectra = (nameFile *)malloc(numberOfFileSpectra*sizeof(nameFile));
         strcpy(vInputFileSpectra[0].name,nameInputFileSpectra);
         vOutputNameModels = (nameFile *)malloc(numberOfFileSpectra*sizeof(nameFile));
         strcpy(vOutputNameModels[0].name,nameOutputFileModels);			
      }
	}

	MPI_Barrier(MPI_COMM_WORLD); // Wait UNTIL THE IMAGE HAS BEEN READED COMPLETELY
	//  BROADCAST THE NUMBER OF LAMBDAS READS FROM THE FILE AND THE NUMBER OF PIXELS
	MPI_Bcast(&numberOfFileSpectra, 1, MPI_INT, root , MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

	/**
	 * What we do know is the following: 
	 * --- If the number of files is greater than number of available processors , then scatter the list of files proporcionaly between N processores. 
	 * --- If the number of files is less than number of available processors, then we process each fits file sequentially doing a scatter of N pixels in the each Fits File. 
	 * */

	if(idProc == root){
		
		printf("\n \n NUMBER OF FILES FOUNDED IN THE DIRECTORY  %d \n\n",numberOfFileSpectra);
		for(i=0;i<numberOfFileSpectra;i++){
			printf("\n %s",vInputFileSpectra[i].name);
		}
	}
	
	int numFilesPerProcess = numberOfFileSpectra / numProcs;
	int numFilesPerProcessParallel = numberOfFileSpectra % numProcs;
	int sum = 0;                // Sum of counts. Used to calculate displacements

	for ( i = 0; i < numProcs; i++) {
		sendcountsNameInputFiles[i] = numFilesPerProcess;
		displsNameInputFiles[i] = sum;
		sum += sendcountsNameInputFiles[i];
	}

	/*MPI_Barrier(MPI_COMM_WORLD); // Wait UNTIL THE IMAGE HAS BEEN READED COMPLETELY
	//  BROADCAST THE NUMBER OF LAMBDAS READS FROM THE FILE AND THE NUMBER OF PIXELS
	MPI_Bcast(&numFilesPerProcessParallel, 1, MPI_INT, root , MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);*/

	//printf("\n IDPROC: %d NÚMERO DE FICHEROS POR PROCESO PARALELO %d , NUMERO DE FICHEROS POR PROCESO %d\n",idProc,numFilesPerProcessParallel,numFilesPerProcess);

	
	if(idProc == root && numFilesPerProcess>=1){
		vInputFileSpectraParalell = (nameFile *)malloc(numFilesPerProcessParallel*sizeof(nameFile));
		vOutputNameModelsParalell = (nameFile *)malloc(numFilesPerProcessParallel*sizeof(nameFile));
		
		for(i=0;i<numFilesPerProcessParallel;i++){
			strcpy(vInputFileSpectraParalell[i].name,vInputFileSpectra[(numFilesPerProcess*numProcs)+i].name);
			strcpy(vOutputNameModelsParalell[i].name,vOutputNameModels[(numFilesPerProcess*numProcs)+i].name);
		}
		nameFile * auxInput  = (nameFile *)malloc((numFilesPerProcess*numProcs)*sizeof(nameFile));
		nameFile * auxOutput = (nameFile *)malloc((numFilesPerProcess*numProcs)*sizeof(nameFile));
		for(i=0;i<(numFilesPerProcess*numProcs);i++){
			strcpy(auxInput[i].name,vInputFileSpectra[i].name);
			strcpy(auxOutput[i].name,vOutputNameModels[i].name);
		}		
		free(vInputFileSpectra);
		free(vOutputNameModels);
		vInputFileSpectra = auxInput;
		vOutputNameModels = auxOutput;
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if(numberOfFileSpectra>numProcs){ // Scatter name of files 
		Init_Model *vModels;
		

		/*int numFilesPerProcess = numberOfFileSpectra / numProcs;
		int resto = numberOfFileSpectra % numProcs;
		int sum = 0;                // Sum of counts. Used to calculate displacements
		int sumSpectro = 0;

		for ( i = 0; i < numProcs; i++) {
			sendcountsNameInputFiles[i] = numFilesPerProcess;
			if (resto > 0) {
					sendcountsNameInputFiles[i]++;
					resto--;
			}
			displsNameInputFiles[i] = sum;
			sum += sendcountsNameInputFiles[i];
		}*/

		vInputFileSpectraLocal = (nameFile *) malloc(sendcountsNameInputFiles[idProc]*sizeof(nameFile));
		vOutputNameModelsLocal = (nameFile *) malloc(sendcountsNameInputFiles[idProc]*sizeof(nameFile));
		
		if( root == idProc){
			MPI_Scatterv(vInputFileSpectra, sendcountsNameInputFiles, displsNameInputFiles, mpiName, vInputFileSpectraLocal, sendcountsNameInputFiles[idProc], mpiName, root, MPI_COMM_WORLD);
			MPI_Scatterv(vOutputNameModels, sendcountsNameInputFiles, displsNameInputFiles, mpiName, vOutputNameModelsLocal, sendcountsNameInputFiles[idProc], mpiName, root, MPI_COMM_WORLD);
			
		}
		else{
			MPI_Scatterv(NULL, NULL,NULL, mpiName, vInputFileSpectraLocal, sendcountsNameInputFiles[idProc], mpiName, root, MPI_COMM_WORLD);
			MPI_Scatterv(NULL, NULL,NULL, mpiName, vOutputNameModelsLocal, sendcountsNameInputFiles[idProc], mpiName, root, MPI_COMM_WORLD);
		}
		MPI_Barrier(MPI_COMM_WORLD);

		// milos over each file of the processor
		int indexInputFits;
		for(indexInputFits=0;indexInputFits<sendcountsNameInputFiles[idProc];indexInputFits++){
			/****************************************************************************************************/
			// READ PIXELS FROM IMAGE 
			PRECISION timeReadImage,timeExecuteClassicalEstimates;
			clock_t t;
			t = clock();
			fitsImage = readFitsSpectroImage(vInputFileSpectraLocal[indexInputFits].name);
			t = clock() - t;
			timeReadImage = ((PRECISION)t)/CLOCKS_PER_SEC; // in seconds 
			
			printf("\n\n IDPROC: %d -->  TIME TO READ FITS IMAGE (%s):  %f seconds to execute \n",idProc, vInputFileSpectraLocal[indexInputFits].name, timeReadImage); 
			//slog_info(0,"\n\n TIME TO READ FITS IMAGE:  %f seconds to execute \n", timeReadImage);

			//if(fitsImage!=NULL && readFitsLambdaFile(nameInputFileLambda,fitsImage)){
			if(fitsImage!=NULL){



				// INTERPOLATE PSF WITH ARRAY OF LAMBDA READ
				/****************************************************************************************************/	
				// THE NUMBER OF LAMBDAS IS READ FROM INPUT FILES 
				//nlambda = fitsImage->nLambdas;
				// check if read stray light
			// COPY LAMBDA READ IN THE TOP OF FILE 
				int contLambda = 0;
				for( i=0;i<fitsImage->numPixels;i++){
					for( j=0;j<nlambda;j++){
						fitsImage->pixels[i].vLambda[j]=vGlobalLambda[j];
						fitsImage->vLambdaImagen[contLambda++] = vGlobalLambda[j];
					}
				}

				//***************************************** INIT MEMORY WITH SIZE OF LAMBDA ****************************************************//
				InitializePointerShareCalculation();
				AllocateMemoryDerivedSynthesis(nlambda);
				
				int indexPixel = 0;

				// ALLOCATE MEMORY FOR STORE THE RESULTS 

				vModels = calloc (fitsImage->numPixels , sizeof(Init_Model));
				vChisqrf = calloc (fitsImage->numPixels , sizeof(PRECISION));
				vNumIter = calloc (fitsImage->numPixels, sizeof(int));
				t = clock();
				
				printf("\nIDPROC: %d -->  PROCESSING INVERSION: %s  \n\n",idProc,vInputFileSpectraLocal[indexInputFits].name);
				//slog_info(0,"\n***********************  PROGRESS INVERSION *******************************\n\n");

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

						clock_t t_ini = clock();
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
						lm_mils(cuantic, wlines, fitsImage->pixels[indexPixel].vLambda, fitsImage->pixels[indexPixel].nLambda, fitsImage->pixels[indexPixel].spectro, fitsImage->pixels[indexPixel].nLambda, &initModel, spectra, &vChisqrf[indexPixel], slightPixel, configCrontrolFile.toplim, configCrontrolFile.NumberOfCycles,
								configCrontrolFile.WeightForStokes, configCrontrolFile.fix, configCrontrolFile.sigma, configCrontrolFile.InitialDiagonalElement,&configCrontrolFile.ConvolveWithPSF,&vNumIter[indexPixel]);						
					}

					vModels[indexPixel] = initModel;
					//vChisqrf[indexPixel] = chisqrf;
					
					//printf ("\t\t %.2f seconds -- %.2f %%\r",  ((PRECISION)(clock() - t)/CLOCKS_PER_SEC) , ((indexPixel*100.)/fitsImage->numPixels));
				}
				t = clock() - t;
				if (configCrontrolFile.UseClassicalEstimates)
					printf("\n\n TIME EXECUTING CLASSICAL ESTIMATES: %f seconds to execute \n", ((PRECISION)timeExecuteClassicalEstimates)/CLOCKS_PER_SEC);
				//slog_info(0,"\n\n TIME EXECUTIN CLASSICAL ESTIMATES: %f seconds to execute \n", ((PRECISION)timeExecuteClassicalEstimates)/CLOCKS_PER_SEC);
				timeReadImage = ((PRECISION)t)/CLOCKS_PER_SEC; // in seconds 
				printf("\nIDPROC --> %d  FINISH EXECUTION OF INVERSION %s : %f seconds to execute \n",idProc, vInputFileSpectraLocal[indexInputFits].name,timeReadImage);
				//slog_info(0,"\n FINISH EXECUTION OF INVERSION: %f seconds to execute \n", timeReadImage);
				//printf("\n**********");
				if(!writeFitsImageModels(vOutputNameModelsLocal[indexInputFits].name,fitsImage->rows,fitsImage->cols,vModels,vChisqrf,vNumIter,configCrontrolFile.saveChisqr)){
						printf("\n ERROR WRITING FILE OF MODELS: %s",nameOutputFileModels);
						//slog_error(0,"\n ERROR WRITING FILE OF MODELS: %s",nameOutputFileModels);
				}

				// PROCESS FILE OF SYNTETIC PROFILES

				if(configCrontrolFile.SaveSynthesisProfile){
					FreeMemoryDerivedSynthesis();
					InitializePointerShareCalculation();
					AllocateMemoryDerivedSynthesis(nlambda);

					//weights_init(configCrontrolFile.sigma, &sig, configCrontrolFile.noise);
					int i;
					for( i=0;i<fitsImage->numPixels;i++)
					{

						Init_Model initModel = vModels[i];
						mil_sinrf(cuantic, &initModel, wlines, fitsImage->pixels[i].vLambda, nlambda, spectra, AH,slight,NULL,configCrontrolFile.ConvolveWithPSF);
						spectral_synthesis_convolution(&nlambda);
						me_der(cuantic, &initModel, wlines, fitsImage->pixels[i].vLambda, nlambda, d_spectra, spectra, spectra, AH, slight,1,configCrontrolFile.ConvolveWithPSF);
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
						//slog_error(0,"\n ERROR WRITING FILE OF SINTHETIC PROFILES: %s",nameOutputFilePerfiles);
					}
				}

				free(vModels);
				free(vChisqrf);
				free(vNumIter);
			}
			else{
				printf("\n\n IDPROC: %d --> FITS FILE: %s WITH THE SPECTRO IMAGE CAN NOT BE READ IT ******************************\n",idProc, vInputFileSpectraLocal[indexInputFits].name);
				//slog_error(0,"\n\n ***************************** FITS FILE WITH THE SPECTRO IMAGE CAN NOT BE READ IT ******************************\n");
			}

			printf(" \n IDPROC: %d --> IMAGE INVERSION FOR IMAGE %s DONE, *********************\n", idProc, vInputFileSpectraLocal[indexInputFits].name);
			//slog_info(0," \n***********************  IMAGE INVERSION DONE, CLEANING MEMORY *********************\n");

			freeFitsImage(fitsImage);
			FreeMemoryDerivedSynthesis();

		}
		
		free(vInputFileSpectraLocal);
		free(vOutputNameModelsLocal);
	}
	else{ // all files as parallel 
		if(idProc==root){
			if(vInputFileSpectraParalell!=NULL)
				free(vInputFileSpectraParalell);
			if(vOutputNameModelsParalell!=NULL)
				free(vOutputNameModelsParalell);
			
			numFilesPerProcessParallel = numberOfFileSpectra;
			vInputFileSpectraParalell = (nameFile *)malloc(numFilesPerProcessParallel*sizeof(nameFile));
			vOutputNameModelsParalell = (nameFile *)malloc(numFilesPerProcessParallel*sizeof(nameFile));
			for(i=0;i<numFilesPerProcessParallel;i++){
				strcpy(vInputFileSpectraParalell[i].name,vInputFileSpectra[i].name);
				strcpy(vOutputNameModelsParalell[i].name,vOutputNameModels[i].name);
			}
			free(vInputFileSpectra);
			free(vOutputNameModels);
		}
		MPI_Barrier(MPI_COMM_WORLD); // Wait UNTIL THE IMAGE HAS BEEN READED COMPLETELY
		MPI_Bcast(&numFilesPerProcessParallel, 1, MPI_INT, root , MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);		
	}
	

	

	int indexInputFits;
	for(indexInputFits=0;indexInputFits<numFilesPerProcessParallel;indexInputFits++){
			//  IF NOT EMPTY READ stray light file
		numPixels=0;
		if(idProc==root && (access(vInputFileSpectraParalell[indexInputFits].name,F_OK)!=-1)){
			printf("\n\n FILE TO READ %s \n", vInputFileSpectraParalell[indexInputFits].name);
			
			clock_t t = clock();
			fitsImage = readFitsSpectroImage(vInputFileSpectraParalell[indexInputFits].name);
			t = clock() - t;
			PRECISION timeReadImage = ((PRECISION)t)/CLOCKS_PER_SEC; // in seconds 
			printf("\n TIME TO READ FITS IMAGE:  %f seconds to execute \n", timeReadImage); 
			
			// COPY LAMBDA READ IN THE TOP OF FILE 
			int contLambda = 0;
			for( i=0;i<fitsImage->numPixels;i++){
				for( j=0;j<nlambda;j++){
					fitsImage->pixels[i].vLambda[j]=vGlobalLambda[j];
					fitsImage->vLambdaImagen[contLambda++] = vGlobalLambda[j];
				}
			}
			numPixels = fitsImage->numPixels;
			/*if(fitsImage!=NULL && readFitsLambdaFile(nameInputFileLambda,fitsImage)){
				nlambda = fitsImage->nLambdas;
				//NLAMBDA = nlambda;
				numPixels = fitsImage->numPixels;
			}
			else{ // EXIT IF THE IMAGE HAS NOT BEEN READ CORRECTLY 
				printf("**** ERROR READING FITS IMAGE ********************");
				return -1; 
			}
			t = clock() - t;
			timeReadImage = ((PRECISION)t)/CLOCKS_PER_SEC; // in seconds 
			printf("\n TIME TO READ FITS IMAGE LAMBDA:  %f seconds to execute \n", timeReadImage); */
			
		}


		MPI_Barrier(MPI_COMM_WORLD); // Wait UNTIL THE IMAGE HAS BEEN READED COMPLETELY
		//  BROADCAST THE NUMBER OF PIXELS
		MPI_Bcast(&numPixels, 1, MPI_INT, root , MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD); // WAIT UNTIL G HAS BEEN READ

		// IF THE NUMBER OF PIXELS IS NOT GREATER THAN 0 WE DON'T CONITUNUE 

		if(numPixels > 0){

			if(idProc == root){
				printf("\nPROCESSING INVERSION: %s  \n\n",vInputFileSpectraParalell[indexInputFits].name);
				resultsInitModelTotal = calloc (numPixels , sizeof(Init_Model));
				chisqrfTotal = calloc (numPixels , sizeof(PRECISION));
				vNumIterTotal = calloc (numPixels, sizeof(int));
			}
			// allocate memory in all processes 
			InitializePointerShareCalculation();
			AllocateMemoryDerivedSynthesis(nlambda);
			
			
			int numPixelsProceso = numPixels/numProcs;
			int resto = numPixels % numProcs;
			int sum = 0;                // Sum of counts. Used to calculate displacements
			int sumSpectro = 0;
			int sumLambda = 0;
			

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

			// SCATTER VPIXELS 
			vSpectraSplit = calloc(sendcountsSpectro[idProc],sizeof(PRECISION));
			vLambdaSplit = calloc(sendcountsLambda[idProc],sizeof(PRECISION));

			
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
			vChisqrf = calloc(sendcountsPixels[idProc], sizeof(PRECISION));
			vNumIter = calloc(sendcountsPixels[idProc], sizeof(int));

			local_start_execution = MPI_Wtime();
			for(indexPixel = 0; indexPixel < sendcountsPixels[idProc]; indexPixel++){
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

				PRECISION auxChisqrf;

				if (configCrontrolFile.UseClassicalEstimates)
				{
					estimacionesClasicas(wlines[1], vLambdaSplit+(indexPixel*(nlambda)), nlambda, vSpectraSplit+(indexPixel*(nlambda*NPARMS)), &initModel);
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
					lm_mils(cuantic, wlines, vLambdaSplit+(indexPixel*(nlambda)), nlambda, vSpectraSplit+(indexPixel*(nlambda*NPARMS)), nlambda, &initModel, spectra, &vChisqrf[indexPixel], slightPixel, configCrontrolFile.toplim, configCrontrolFile.NumberOfCycles,
						configCrontrolFile.WeightForStokes, configCrontrolFile.fix, configCrontrolFile.sigma, configCrontrolFile.InitialDiagonalElement,&configCrontrolFile.ConvolveWithPSF,&vNumIter[indexPixel]);																		
				}
				resultsInitModel[indexPixel] = initModel;
				//vChisqrf[indexPixel] = auxChisqrf;
			}
			local_finish_execution = MPI_Wtime();

			local_start_gather = MPI_Wtime();


			MPI_Igatherv(resultsInitModel, sendcountsPixels[idProc], mpiInitModel, resultsInitModelTotal, sendcountsPixels, displsPixels, mpiInitModel, root, MPI_COMM_WORLD,&mpiRequestSpectro);
			MPI_Igatherv(vChisqrf, sendcountsPixels[idProc], MPI_DOUBLE, chisqrfTotal, sendcountsPixels, displsPixels, MPI_DOUBLE, root, MPI_COMM_WORLD,&mpiRequestLambda);		
			MPI_Gatherv(vNumIter, sendcountsPixels[idProc], MPI_INT, vNumIterTotal, sendcountsPixels, displsPixels, MPI_INT, root, MPI_COMM_WORLD);		
			MPI_Wait(&mpiRequestSpectro, MPI_STATUS_IGNORE);
			MPI_Wait(&mpiRequestLambda, MPI_STATUS_IGNORE);		
			local_finish_gather = MPI_Wtime();
			

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
				/*time(&t_now);
				struct tm *  t_local = = localtime(&t_now);
				char s_date[256];
				sprintf(s_date, "%02d/%02d/%d",t_local->tm_mday , (t_local->tm_mon + 1), (t_local->tm_year + 1900));*/
				if(!writeFitsImageModels(vOutputNameModelsParalell[indexInputFits].name,fitsImage->rows,fitsImage->cols,resultsInitModelTotal,chisqrfTotal,vNumIterTotal,configCrontrolFile.saveChisqr)){
						printf("\n ERROR WRITING FILE OF MODELS: %s",vOutputNameModelsParalell[indexInputFits].name);
				}
				t = clock() - t;
				timeWriteImage = ((double)t)/CLOCKS_PER_SEC; // in seconds 
				printf("\n TIME TO WRITE FITS IMAGE:  %f seconds to execute \n", timeWriteImage); 
				// PROCESS FILE OF SYNTETIC PROFILES
				if(configCrontrolFile.SaveSynthesisProfile){
					FreeMemoryDerivedSynthesis();
					InitializePointerShareCalculation();
					AllocateMemoryDerivedSynthesis(nlambda);						
					for(indexPixel=0;indexPixel<fitsImage->numPixels;indexPixel++)
					{

						Init_Model initModel = resultsInitModelTotal[indexPixel];
						//double chisqr = vChisqrf[i];

						
						mil_sinrf(cuantic, &initModel, wlines, fitsImage->pixels[indexPixel].vLambda, nlambda, spectra, AH,slight,NULL,configCrontrolFile.ConvolveWithPSF);
						spectral_synthesis_convolution(&nlambda);							
						me_der(cuantic, &initModel, wlines, fitsImage->pixels[indexPixel].vLambda, nlambda, d_spectra,spectra, spectra, AH, slight,1,configCrontrolFile.ConvolveWithPSF);
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

				free(resultsInitModelTotal);		
				free(chisqrfTotal);
				free(vNumIterTotal);
				printf(" \n IMAGE INVERSION OF %s DONE\n", vInputFileSpectraParalell[indexInputFits].name);
			}
			else{
				free(vSpectraSplit);
				free(vLambdaSplit);
				free(resultsInitModel);				
				free(vChisqrf);
				free(vNumIter);
			}
			
			FreeMemoryDerivedSynthesis();
			
		}
		else{
			if(idProc==root){
				printf("\n\n ***************************** FITS FILE CAN NOT BE READ IT %s ******************************",vInputFileSpectraParalell[indexInputFits].name);
			}
		}
		if(idProc==root){
			freeFitsImage(fitsImage);
		}
	}		
		

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

	if(configCrontrolFile.ConvolveWithPSF){
		fftw_free(inSpectraFwPSF);
		fftw_free(outSpectraFwPSF);
		fftw_destroy_plan(planForwardPSF);
		fftw_free(inSpectraBwPSF);
		fftw_free(outSpectraBwPSF);
		fftw_destroy_plan(planBackwardPSF);

		fftw_free(fftw_G_PSF);
		fftw_free(fftw_G_MAC_PSF);
		fftw_free(fftw_G_MAC_DERIV_PSF);

		fftw_free(inPSF_MAC);
		fftw_free(inMulMacPSF);
		fftw_free(inPSF_MAC_DERIV);
		fftw_free(inMulMacPSFDeriv);
		fftw_free(outConvFilters);
		fftw_free(outConvFiltersDeriv);	

		fftw_destroy_plan(planForwardPSF_MAC);
		fftw_destroy_plan(planForwardPSF_MAC_DERIV);
		fftw_destroy_plan(planBackwardPSF_MAC);
		fftw_destroy_plan(planBackwardPSF_MAC_DERIV);		
	}
	//printf("\n\n TOTAL sec : %.16g segundos\n", total_secs);
	free(cuantic);
	free(wlines);
	free(vGlobalLambda);
	//FreeMemoryDerivedSynthesis();
	// FREE TYPE OF MPI
	MPI_Type_free(&mpiInitModel);
	MPI_Finalize() ;
	free(G);
	return 0;
}

