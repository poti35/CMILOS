
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
//
// How to use:
//
//  >> milos_Fits parameters.txt
//
//

#include "mpi.h"
#include <time.h>
#include "defines.h"
//#include "nrutil.h"
#include "convolution.c"
#include <string.h>
#include <stdio.h>
#include <stddef.h>
#include "fitsio.h"
#include "utilsFits.h"
#include "mpiTools.h"
#include "lib.h"
#include "milosUtils.h"


// ***************************** FUNCTIONS TO READ FITS FILE *********************************************************


long long int c1, c2, cd, semi, c1a, c2a, cda; //variables de 64 bits para leer ciclos de reloj
long long int c1total, c2total, cdtotal, ctritotal;

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
int FGlobal, HGlobal, uuGlobal;

PRECISION *d_spectra, *spectra;

//Number of lambdas in the input profiles
int NLAMBDA = 0;

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

	// INIT MPI  PROGRAM 
	int numProcs, idProc;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
   MPI_Comm_rank(MPI_COMM_WORLD, &idProc);
	MPI_Datatype mpiInitModel;
	MPI_Datatype mpiVPixels;	
	MPI_Datatype mpiCuantic;
	const int root=0;	
		
	// FINISH STARTING PROGRAM 

	double *wlines;
	int nwlines, nlambda, numPixels, iter, miter, nweight, indexPixel;
	Init_Model initModel;
	double * chisqrf;
	double slight, toplim;

	PRECISION weight[4] = {1., 1., 1., 1.};

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
	int Max_iter;
   FitsImage * fitsImage;

	double dat[7] = {CUANTIC_NWL, CUANTIC_SLOI, CUANTIC_LLOI, CUANTIC_JLOI, CUANTIC_SUPI, CUANTIC_LUPI, CUANTIC_JUPI};

	/********************* Read data input from file ******************************/

	if(!readParametersFileInput(argv[1], &Max_iter,&CLASSICAL_ESTIMATES,&PRINT_SINTESIS,nameInputFileSpectra,nameInputFileLambda,nameOutputFileModels,nameOutputFilePerfiles,&INSTRUMENTAL_CONVOLUTION,&FWHM,&DELTA,&NMUESTRAS_G)){
		printf("\n********************* EXITING THE PROGRAM . ERROR READING PARAMETERS FILE ****************************\n");
		return -1;
	}
	/******************************************************************************/


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
	miter = Max_iter;
	
	numPixels=0;

	/****************************************************************************************************/

	
	/******************* APPLY GAUSSIAN, CREATE CUANTINC AND INITIALIZE DINAMYC MEMORY*******************/
	if (INSTRUMENTAL_CONVOLUTION)
	{
		generateGaussianInstrumentalProfile(G,FWHM,DELTA,INSTRUMENTAL_CONVOLUTION_WITH_PSF,NMUESTRAS_G);
	}

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
		if(fitsImage!=NULL && readFitsLambdaFile(nameInputFileLambda,fitsImage)){
			nlambda = fitsImage->nLambdas;
			NLAMBDA = nlambda;
			numPixels = fitsImage->numPixels;
		}
		else{ /* EXIT IF THE IMAGE HAS NOT BEEN READ CORRECTLY */
			printf("**** ERROR READING FITS IMAGE ********************");
			return -1; 
		}
	}

	MPI_Barrier(MPI_COMM_WORLD); // Wait UNTIL THE IMAGE HAS BEEN READED COMPLETELY

	//  BROADCAST THE NUMBER OF LAMBDAS READS FROM THE FILE AND THE NUMBER OF PIXELS
	MPI_Bcast(&nlambda, 1, MPI_INT, root , MPI_COMM_WORLD);
	MPI_Bcast(&numPixels, 1, MPI_INT, root , MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD); // Wait UNTIL THE IMAGE HAS BEEN READED COMPLETELY	

	// IF THE NUMBER OF PIXELS IS NOT GREATER THAN 0 WE DON'T CONITUNUE 

	if(numPixels > 0){
		printf("\n************************************************");		
		printf("\n************************************************");		
		printf("\n ESTOY EN EL PROCESO %d . Número de pixeles leidos: %d", idProc, numPixels);
		printf("\n************************************************");
		Init_Model * resultsInitModel;
		Init_Model * resultsInitModelTotal;
		double * chisqrfTotal;
		double * vLambda;

		// allocate memory in all processes 
		AllocateMemoryDerivedSynthesis(nlambda);
		//initializing weights
		PRECISION *w, *sig;
		weights_init(nlambda, sigma, weight, nweight, &w, &sig, noise);

		//buildMpiType(&mpiInitModel,nlambda, &mpiVPixels,cuantic[0].N_SIG,cuantic[0].N_PI,&mpiCuantic);
		printf("\n************************************************");		
		printf("\n************************************************");		
		printf("\n ESTOY EN EL PROCESO %d . TIPO DE DATOS CONTRUIDOS", idProc);
		printf("\n************************************************");
		printf("\n************************************************");


		int numPixelsProceso = numPixels/numProcs;
		int resto = numPixels % numProcs;
		int sum = 0;                // Sum of counts. Used to calculate displacements
		int maxSendCount = 0;
		int *sendcounts = malloc(sizeof(int)*numProcs); // array describing how many elements to send to each process
		int *displs = malloc(sizeof(int)*numProcs);  // array describing the displacements where each segment begins
		for (int i = 0; i < numProcs; i++) {
			sendcounts[i] = numPixelsProceso;
			if (resto > 0) {
					sendcounts[i]++;
					resto--;
			}
			displs[i] = sum;
			sum += sendcounts[i];
			if(sendcounts[i]>maxSendCount) maxSendCount = sendcounts[i];
		}

		// BROADCAST VECTOR OF LAMBDAS 
		if(root==idProc){
			vLambda = calloc (nlambda, sizeof(double));
			for(int i=0;i<nlambda;i++){
				vLambda[i]=fitsImage->pixels[0].vLambda[i];
			}
		}
		MPI_Bcast(vLambda, nlambda, MPI_DOUBLE, root , MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD); // Wait until all processes have theis vlambda

		for(int i=0;i<nlambda;i++){
			if(i==0){
				printf("\n IMPRIMIENDO EL VECTOR DE LAMBDA DEL PROCESO %d: %f", idProc, vLambda[i]);
			}
			else
			{
				printf(" %f", vLambda[i]);
			}
		}
		printf("\n***********************************************\n");
		const int nitemsStructVPixel = 3;
		int blocklenghtVPixel [3] = {nlambda,nlambda*NPARMS,1};
		MPI_Datatype typesVPixel [3] = {MPI_DOUBLE,MPI_DOUBLE,MPI_INT};
		const MPI_Aint offsetsVPixels [3] = { 0, sizeof(PRECISION)*nlambda, (sizeof(PRECISION)*nlambda) + (sizeof(PRECISION)*nlambda*4) };
		/* offsetsVPixels[0] = offsetof(vpixels, vLambda);
		offsetsVPixels[1] = offsetof(vpixels, spectro);
		offsetsVPixels[2] = offsetof(vpixels, nLambda);*/
		MPI_Type_create_struct(nitemsStructVPixel, blocklenghtVPixel, offsetsVPixels, typesVPixel, &mpiVPixels);
		MPI_Type_commit(&mpiVPixels);

		if (root == idProc) {
			for (int i = 0; i < numProcs; i++) {
				printf("\n sendcounts[%d] = %d\tdispls[%d] = %d\n", i, sendcounts[i], i, displs[i]);
				printf("\n**********************************************");
			}
		}
		// SCATTER VPIXELS 
		vpixels * vPixelsProcess = calloc(maxSendCount, sizeof(vpixels));
		MPI_Barrier(MPI_COMM_WORLD); // Wait until all processes have their pixels
		if(idProc==root){
			MPI_Scatterv(fitsImage->pixels, sendcounts, displs, mpiVPixels, vPixelsProcess, sendcounts[idProc], mpiVPixels, root, MPI_COMM_WORLD);
		}
		else{
			MPI_Scatterv(NULL, NULL, NULL, mpiVPixels, vPixelsProcess, sendcounts[idProc], mpiVPixels, root, MPI_COMM_WORLD);
		}
		printf("\n************************************************");
		printf("\n SCATTER HECHO EN PROCESO %d",idProc);
		printf("\n************************************************");		
		MPI_Barrier(MPI_COMM_WORLD); // Wait until all processes have their pixels


		

		// READ EACH PIXEL FROM THE IMAGE 

		// USE BROADCASTING TO SEND ALL THE DATA NEEDED FOR PROCESS THE PIXELS. 
		resultsInitModel = calloc(sendcounts[idProc],sizeof(Init_Model));
		chisqrf = calloc(sendcounts[idProc], sizeof(double));

		
		for(indexPixel = 0; indexPixel < sendcounts[idProc]; indexPixel++){
			

			//Initial Model

			resultsInitModel[indexPixel].eta0 = INITIAL_MODEL_ETHA0;
			resultsInitModel[indexPixel].B = INITIAL_MODEL_B; //200 700
			resultsInitModel[indexPixel].gm = INITIAL_MODEL_GM;
			resultsInitModel[indexPixel].az = INITIAL_MODEL_AZI;
			resultsInitModel[indexPixel].vlos = INITIAL_MODEL_VLOS; //km/s 0
			resultsInitModel[indexPixel].mac = 0.0;
			resultsInitModel[indexPixel].dopp = INITIAL_MODEL_LAMBDADOPP;
			resultsInitModel[indexPixel].aa = INITIAL_MODEL_AA;
			resultsInitModel[indexPixel].alfa = 1; //0.38; //stray light factor
			resultsInitModel[indexPixel].S0 = INITIAL_MODEL_S0;
			resultsInitModel[indexPixel].S1 = INITIAL_MODEL_S1;

			if (CLASSICAL_ESTIMATES)
			{
				estimacionesClasicas(wlines[1], vLambda, nlambda, vPixelsProcess[indexPixel].spectro, &resultsInitModel[indexPixel]);
				//Se comprueba si el resultado fue "nan" en las CE
				if (isnan(resultsInitModel[indexPixel].B))
					resultsInitModel[indexPixel].B = 1;
				if (isnan(resultsInitModel[indexPixel].vlos))
					resultsInitModel[indexPixel].vlos = 1e-3;
				if (isnan(resultsInitModel[indexPixel].gm))
					resultsInitModel[indexPixel].gm = 1;
				if (isnan(resultsInitModel[indexPixel].az))
					resultsInitModel[indexPixel].az = 1;

			}

			//inversion
			if (CLASSICAL_ESTIMATES != 2)
			{
				//Se introduce en S0 el valor de Blos si solo se calculan estimaciones clásicas
				//Aqui se anula esa asignación porque se va a realizar la inversion RTE completa
				resultsInitModel[indexPixel].S0 = INITIAL_MODEL_S0;

				lm_mils(cuantic, wlines, nwlines, vLambda, nlambda, vPixelsProcess[indexPixel].spectro, nlambda, &initModel, spectra, &chisqrf[indexPixel], &iter, slight, toplim, miter,
						weight, nweight, fix, sig, filter, ilambda, noise, pol, getshi, 0,&INSTRUMENTAL_CONVOLUTION,&NMUESTRAS_G);
			}

		}
		if(idProc ==0){
			resultsInitModelTotal = calloc(numPixels,sizeof(Init_Model));
			chisqrfTotal = calloc(numPixels, sizeof(double));
		}
		MPI_Gatherv(resultsInitModel, sendcounts[idProc], mpiInitModel, resultsInitModelTotal, sendcounts, displs, mpiInitModel, root, MPI_COMM_WORLD);
		MPI_Gatherv(chisqrf, sendcounts[idProc], MPI_DOUBLE, chisqrfTotal, sendcounts, displs, MPI_DOUBLE, root, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);  // FINISH UNTIL EVERY PROCESS FINISHED
 		
		free(resultsInitModel);
		free(chisqrf);
		free(vLambda);
		if(idProc==0){		 	  
			if(!writeFitsImageModels(nameOutputFileModels,fitsImage->rows,fitsImage->cols,resultsInitModelTotal,chisqrfTotal)){
					printf("\n ERROR WRITING FILE OF MODELS: %s",nameOutputFileModels);
			}
			// PROCESS FILE OF SYNTETIC PROFILES
			if(PRINT_SINTESIS){
				for(indexPixel=0;indexPixel<fitsImage->numPixels;indexPixel++)
				{

					Init_Model initModel = resultsInitModelTotal[indexPixel];
					//double chisqr = vChisqrf[i];
					int NMODEL = 12; //Numero de parametros del modelo

					mil_sinrf(cuantic, &initModel, wlines, nwlines, fitsImage->pixels[indexPixel].vLambda, nlambda, spectra, AH, 0, filter);

					me_der(cuantic, &initModel, wlines, nwlines, fitsImage->pixels[indexPixel].vLambda, nlambda, d_spectra, AH, slight, 0, filter);
					response_functions_convolution(&nlambda,&INSTRUMENTAL_CONVOLUTION,&NMUESTRAS_G);
					int kk;
					for (kk = 0; kk < (nlambda * NPARMS); kk++)
					{
						fitsImage->pixels[indexPixel].spectro[kk] = spectra[kk] ;
						//fprintf(fOutput,"%lf %le %le %le %le \n", lambda[kk], spectra[kk], spectra[kk + NLAMBDA], spectra[kk + NLAMBDA * 2], spectra[kk + NLAMBDA * 3]);
					}
				}
				// WRITE SINTHETIC PROFILES TO FITS FILE
				if(!writeFitsImageProfiles(nameOutputFilePerfiles,nameInputFileSpectra,fitsImage)){
					printf("\n ERROR WRITING FILE OF SINTHETIC PROFILES: %s",nameOutputFilePerfiles);
				}
			}
			free(resultsInitModelTotal);
			free(chisqrfTotal);
		}

	}
	else{
		printf("\n\n ***************************** FITS FILE CAN NOT BE READ IT ******************************");
	}
	if(idProc==0){
		freeFitsImage(fitsImage);
	}
	

	//printf("\n\n TOTAL sec : %.16g segundos\n", total_secs);

	free(cuantic);
	free(wlines);


	FreeMemoryDerivedSynthesis();
	// FREE TYPE OF MPI
	MPI_Type_free(&mpiInitModel);
	MPI_Type_free(&mpiVPixels);
	MPI_Type_free(&mpiCuantic);
	MPI_Finalize() ;

	free(G);

	return 0;
}

