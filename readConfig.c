#include "readConfig.h"
#include "defines.h"
#include <time.h>
#include <string.h>
#include <stdio.h>
#include <stddef.h>
#include <locale.h>
#include <stdlib.h>
#include <libconfig.h>



int readConfigControl(char * configFile, ConfigControl * trolConfig, int printLog){


	config_t cfg;
	config_init(&cfg);
	  /* Read the file. If there is an error, report it and exit. */
	if(! config_read_file(&cfg, configFile))
	{
		if(printLog) printf("%s:%d - %s\n", config_error_file(&cfg),config_error_line(&cfg), config_error_text(&cfg));
		config_destroy(&cfg);
		return(EXIT_FAILURE);
	}

	//  number of iterations 
	
	if(config_lookup_int(&cfg, NUMBER_OF_CYCLES, &trolConfig->NumberOfCycles)){
    	if(printLog) printf("NumberOfCycles to apply: %d\n", trolConfig->NumberOfCycles);
	}
  	else
    	if(printLog) printf("No 'NumberOfIterations' setting in configuration file.\n");

	//  file with spectro observed 
	if(config_lookup_string(&cfg, OBSERVED_PROFILES, &trolConfig->ObservedProfiles)){
    	if(printLog) printf("ObservedProfiles to apply: %s\n", trolConfig->ObservedProfiles);
	}
  	else
    	if(printLog) printf("No 'ObservedProfiles' setting in configuration file.\n");

	//  file with stray light 
	if(config_lookup_string(&cfg, STRAY_LIGHT_FILE, &trolConfig->StrayLightFile)){
    	if(printLog) printf("StrayLightFile to apply: %s\n", trolConfig->StrayLightFile);
	}
  	else
    	if(printLog) printf("No 'StrayLightFile' setting in configuration file.\n");

	//  file with psf of instrument 
	if(config_lookup_string(&cfg, PSF_FILE, &trolConfig->PSFFile)){
    	if(printLog) printf("PSFFile to apply: %s\n", trolConfig->PSFFile);
	}
  	else
    	if(printLog) printf("No 'PSFFile' setting in configuration file.\n");		

	//  fits file with wave lenght 
	if(config_lookup_string(&cfg, WAVE_LENGHT_FILE, &trolConfig->WavelengthFile)){
    	if(printLog) printf("WavelengthFile to apply: %s\n", trolConfig->WavelengthFile);
	}
  	else
    	if(printLog) printf("No 'WavelengthFile' setting in configuration file.\n");		

	//  file with atomic lines 
	if(config_lookup_string(&cfg, ATOMIC_PARAMETERS_FILE, &trolConfig->AtomicParametersFile)){
    	if(printLog) printf("AtomicParametersFile to apply: %s\n", trolConfig->AtomicParametersFile);
	}
  	else
    	if(printLog) printf("No 'AtomicParametersFile' setting in configuration file.\n");

	//  READ INITIAL GUESS MODEL 
	if(config_lookup_string(&cfg, INITIAL_GUESS_MODEL, &trolConfig->InitialGuessModel)){
    	if(printLog) printf("InitialGuessModel to apply: %s\n", trolConfig->InitialGuessModel);
	}
  	else
    	if(printLog) printf("No 'InitialGuessModel' setting in configuration file.\n");		

	//  Weight for stokes I  
	if(config_lookup_float(&cfg, WEIGHT_FOR_STOKESI, &trolConfig->WeightForStokes[0])){
    	if(printLog) printf("WeightForStokesI to apply: %f\n", trolConfig->WeightForStokes[0]);
	}
  	else
    	if(printLog) printf("No 'WeightForStokesI' setting in configuration file.\n");

	//  Weight for stokes Q  
	if(config_lookup_float(&cfg, WEIGHT_FOR_STOKESQ, &trolConfig->WeightForStokes[1])){
    	if(printLog) printf("WeightForStokesQ to apply: %f\n", trolConfig->WeightForStokes[1]);
	}
  	else
    	if(printLog) printf("No 'WeightForStokesQ' setting in configuration file.\n");		

	//  Weight for stokes U  
	if(config_lookup_float(&cfg, WEIGHT_FOR_STOKESU, &trolConfig->WeightForStokes[2])){
    	if(printLog) printf("WeightForStokesU to apply: %f\n", trolConfig->WeightForStokes[2]);
	}
  	else
    	if(printLog) printf("No 'WeightForStokesU' setting in configuration file.\n");		

	//  Weight for stokes V 
	if(config_lookup_float(&cfg, WEIGHT_FOR_STOKESV, &trolConfig->WeightForStokes[3])){
    	if(printLog) printf("WeightForStokesV to apply: %f\n", trolConfig->WeightForStokes[3]);
	}
  	else
    	if(printLog) printf("No 'WeightForStokesV' setting in configuration file.\n");

	//  invert macroturbulence
	if(config_lookup_bool(&cfg, INVERT_MACROTURBULENCE, &trolConfig->InvertMacroturbulence)){
    	if(printLog) printf("InvertMacroturbulence to apply: %d\n", trolConfig->InvertMacroturbulence);
	}
  	else
    	if(printLog) printf("No 'InvertMacroturbulence' setting in configuration file.\n");		

	//  invert filling factor 
	if(config_lookup_bool(&cfg, INVERT_FILLING_FACTOR, &trolConfig->InvertFillingFactor)){
    	if(printLog) printf("InvertFillingFactor to apply: %d\n", trolConfig->InvertFillingFactor);
	}
  	else
    	if(printLog) printf("No 'InvertFillingFactor' setting in configuration file.\n");		

	//  invert stray light factor
	if(config_lookup_bool(&cfg, INVERT_STRAY_LIGHT_FACTOR, &trolConfig->InvertStrayLightFactor)){
    	if(printLog) printf("InvertStrayLightFactor to apply: %d\n", trolConfig->InvertStrayLightFactor);
	}
  	else
    	if(printLog) printf("No 'InvertStrayLightFactor' setting in configuration file.\n");

	//  Read mu : Scalar containing the cosine of the heliocentric angle 
	if(config_lookup_float(&cfg, MU, &trolConfig->mu)){
    	if(printLog) printf("mu to apply: %f\n", trolConfig->mu);
	}
  	else
    	if(printLog) printf("No 'mu' setting in configuration file.\n");		

	//  Estimated S/N FOR I 
	if(config_lookup_int(&cfg, ESTIMATEDSNFORI, &trolConfig->EstimatedSNForI)){
    	if(printLog) printf("EstimatedSNForI to apply: %d\n", trolConfig->EstimatedSNForI);
	}
  	else
    	if(printLog) printf("No 'EstimatedSNForI' setting in configuration file.\n");		

	//  Continuum Contrast
	if(config_lookup_int(&cfg, CONTINUUM_CONTRAST, &trolConfig->ContinuumContrast)){
    	if(printLog) printf("ContinuumContrast to apply: %d\n", trolConfig->ContinuumContrast);
	}
  	else
    	if(printLog) printf("No 'ContinuumContrast' setting in configuration file.\n");			

	//  Tolerance for SVD 
	if(config_lookup_float(&cfg, TOLERANCE_FOR_SVD, &trolConfig->ToleranceForSVD)){
    	if(printLog) printf("ToleranceForSVD to apply: %f\n", trolConfig->ToleranceForSVD);
	}
  	else
    	if(printLog) printf("No 'ToleranceForSVD' setting in configuration file.\n");	

	//  initial diagonal element
	if(config_lookup_float(&cfg, INITIAL_DIAGONAL_ELEMENT, &trolConfig->InitialDiagonalElement)){
    	if(printLog) printf("InitialDiagonalElement to apply: %f\n", trolConfig->InitialDiagonalElement);
	}
  	else
    	if(printLog) printf("No 'InitialDiagonalElement' setting in configuration file.\n");	
			

	//  Use interpolar with splines or linear
	if(config_lookup_int(&cfg, USE_INTERPOLAR_SPLINES_OR_LINEAR, &trolConfig->useInterpolarSplinesOrLinear)){
    	if(printLog) printf("useInterpolarSplinesOrLinear to apply: %d\n", trolConfig->useInterpolarSplinesOrLinear);
	}
  	else
    	if(printLog) printf("No 'useInterpolarSplinesOrLinear' setting in configuration file.\n");					

	//  Indicate if convolve with PSF
	if(config_lookup_bool(&cfg, CONVOLVE_WITH_PSF, &trolConfig->ConvolveWithPSF)){
    	if(printLog) printf("ConvolveWithPSF to apply: %d\n", trolConfig->ConvolveWithPSF);
	}
  	else
    	if(printLog) printf("No 'ConvolveWithPSF' setting in configuration file.\n");							

	//  FWHM
	if(config_lookup_float(&cfg, FWHM_FILE, &trolConfig->FWHM)){
    	if(printLog) printf("FWHM to apply: %f\n", trolConfig->FWHM);
	}
  	else
    	if(printLog) printf("No 'FWHM' setting in configuration file.\n");			

	//  TYPE CONVOLUTION TO USE 
	if(config_lookup_string(&cfg, TYPE_CONVOLUTION, &trolConfig->TypeConvolution)){
    	if(printLog) printf("TypeConvolution to apply: %s\n", trolConfig->TypeConvolution);
	}
  	else
    	if(printLog) printf("No 'TypeConvolution' setting in configuration file.\n");			

	//  Gas pressure at surface 1
	if(config_lookup_float(&cfg, GAS_PRESSURE_AT_SURFACE_1, &trolConfig->GasPressureAtSurface1)){
    	if(printLog) printf("GasPressureAtSurface1to apply: %f\n", trolConfig->GasPressureAtSurface1);
	}
  	else
    	if(printLog) printf("No 'GasPressureAtSurface1' setting in configuration file.\n");					

	//  Gas pressure at surface 2
	if(config_lookup_float(&cfg, GAS_PRESSURE_AT_SURFACE_2, &trolConfig->GasPressureAtSurface2)){
    	if(printLog) printf("GasPressureAtSurface2 to apply: %f\n", trolConfig->GasPressureAtSurface2);
	}
  	else
    	if(printLog) printf("No 'GasPressureAtSurface2' setting in configuration file.\n");			

	//  Magnetic pressure term
	if(config_lookup_float(&cfg, MAGNETIC_PRESSURE_TERM, &trolConfig->MagneticPressureTerm)){
    	if(printLog) printf("MagneticPressureTerm to apply: %f\n", trolConfig->MagneticPressureTerm);
	}
  	else
    	if(printLog) printf("No 'MagneticPressureTerm' setting in configuration file.\n");		

	//  Numbe of lines tu use
	if(config_lookup_int(&cfg, NTL, &trolConfig->ntl)){
    	if(printLog) printf("ntl to apply: %d\n", trolConfig->ntl);
	}
  	else
    	if(printLog) printf("No 'ntl' setting in configuration file.\n");		

	//  Number of wavelenght observed
	if(config_lookup_int(&cfg, NLIOBS, &trolConfig->nliobs)){
    	if(printLog) printf("nliobs to apply: %d\n", trolConfig->nliobs);
	}
  	else
    	if(printLog) printf("No 'nliobs' setting in configuration file.\n");	

	//  Central wavelenght 
	if(config_lookup_float(&cfg, CENTRAL_WAVE_LENGHT, &trolConfig->CentralWaveLenght)){
    	if(printLog) printf("CentralWaveLenght to apply: %f\n", trolConfig->CentralWaveLenght);
	}
  	else
    	if(printLog) printf("No 'CentralWaveLenght' setting in configuration file.\n");	

	//  indicate if calculate eta0 for compute output model
	if(config_lookup_int(&cfg, ETA0_LINE_TO_CONTINUUM_ABSORPTION, &trolConfig->fix[0])){
    	if(printLog) printf("ETA0 to apply: %d\n", trolConfig->fix[0]);
	}
  	else
    	if(printLog) printf("No 'ETA0' setting in configuration file.\n");	

	//  indicate if calculate B for compute output model
	if(config_lookup_int(&cfg, B_MAGNETIC_FIELD_STRENGTH, &trolConfig->fix[1])){
    	if(printLog) printf("B to apply: %d\n", trolConfig->fix[1]);
	}
  	else
    	if(printLog) printf("No 'B' setting in configuration file.\n");	

	//  Indicate if calculate VLOS for output model
	if(config_lookup_int(&cfg, VLOS_LINE_OF_SIGHT_VELOCITY, &trolConfig->fix[2])){
    	if(printLog) printf("VLOS to apply: %d\n", trolConfig->fix[2]);
	}
  	else
    	if(printLog) printf("No 'VLOS' setting in configuration file.\n");			

	//  Indicate if calculate DOPP 
	if(config_lookup_int(&cfg, DOPP_DOOPLER_WIDTH, &trolConfig->fix[3])){
    	if(printLog) printf("DOPP to apply: %d\n", trolConfig->fix[3]);
	}
  	else
    	if(printLog) printf("No 'DOPP' setting in configuration file.\n");		

	//  indicate if calculate aa  
	if(config_lookup_int(&cfg, AA_DAMPING_PARAMETER, &trolConfig->fix[4])){
    	if(printLog) printf("AA to apply: %d\n", trolConfig->fix[4]);
	}
  	else
    	if(printLog) printf("No 'AA' setting in configuration file.\n");		

	//  indicate if calculate gm
	if(config_lookup_int(&cfg, GM_MAGNETIC_FIELD_INCLINATION, &trolConfig->fix[5])){
    	if(printLog) printf("GM to apply: %d\n", trolConfig->fix[5]);
	}
  	else
    	if(printLog) printf("No 'GM' setting in configuration file.\n");		

	//  indicate if calculate azimuth 
	if(config_lookup_int(&cfg, AZ_MAGNETIC_FIELD_AZIMUTH, &trolConfig->fix[6])){
    	if(printLog) printf("AZ to apply: %d\n", trolConfig->fix[6]);
	}
  	else
    	if(printLog) printf("No 'AZ' setting in configuration file.\n");		

	//  indicate if calculate s0
	if(config_lookup_int(&cfg, S0_SOURCE_FUNCTION_CONSTANT, &trolConfig->fix[7])){
    	if(printLog) printf("S0 to apply: %d\n", trolConfig->fix[7]);
	}
  	else
    	if(printLog) printf("No 'S0' setting in configuration file.\n");	

	//  indicate if calculate s1 
	if(config_lookup_int(&cfg, S1_SOURCE_FUNCTION_GRADIENT, &trolConfig->fix[8])){
    	if(printLog) printf("S1 to apply: %d\n", trolConfig->fix[8]);
	}
  	else
    	if(printLog) printf("No 'S1' setting in configuration file.\n");						

	//  indicate if calculate macroturbulence
	if(config_lookup_int(&cfg, MAC_MACROTURBULENT_VELOCITY, &trolConfig->fix[9])){
    	if(printLog) printf("MAC to apply: %d\n", trolConfig->fix[9]);
	}
  	else
    	if(printLog) printf("No 'MAC' setting in configuration file.\n");			

	//  indicate if calculate filling factor 
	if(config_lookup_int(&cfg, ALPHA_FILLING_FACTOR, &trolConfig->fix[10])){
    	if(printLog) printf("ALPHA to apply: %d\n", trolConfig->fix[10]);
	}
  	else
    	if(printLog) printf("No 'ALPHA' setting in configuration file.\n");	

	// if save chisqr 
	if(config_lookup_bool(&cfg, SAVE_CHISQR, &trolConfig->saveChisqr)){
    	if(printLog) printf("saveChisqr to apply: %d\n", trolConfig->saveChisqr);
	}
  	else{
    	if(printLog) printf("No 'saveChisqr' setting in configuration file. Default = FALSE\n");	
		trolConfig->saveChisqr=0;
	}

	//  if use classical estimates
	if(config_lookup_bool(&cfg, USE_CLASSICAL_ESTIMATES, &trolConfig->UseClassicalEstimates)){
    	if(printLog) printf("UseClassicalEstimates to apply: %d\n", trolConfig->UseClassicalEstimates);
	}
  	else
    	if(printLog) printf("No 'UseClassicalEstimates' setting in configuration file.\n");	

	//  if use rte inversion 
	if(config_lookup_bool(&cfg, USE_RTE_INVERSION, &trolConfig->UseRTEInversion)){
    	if(printLog) printf("UseRTEInversion to apply: %d\n", trolConfig->UseRTEInversion);
	}
  	else
    	if(printLog) printf("No 'UseRTEInversion' setting in configuration file.\n");	

	//  if save synthesis profile
	if(config_lookup_bool(&cfg, SAVE_SYNTHESIS_PROFILE, &trolConfig->SaveSynthesisProfile)){
    	if(printLog) printf("SaveSynthesisProfile to apply: %d\n", trolConfig->SaveSynthesisProfile);
	}
  	else
    	if(printLog) printf("No 'SaveSynthesisProfile' setting in configuration file.\n");	

	//  Output model file 
	if(config_lookup_string(&cfg, OUTPUT_MODEL_FILE, &trolConfig->OutputModelFile)){
    	if(printLog) printf("OutputModelFile to apply: %s\n", trolConfig->OutputModelFile);
	}
  	else
    	if(printLog) printf("No 'OutputModelFile' setting in configuration file.\n");	

	//  Output synthesis file 
	if(config_lookup_string(&cfg, OUTPUT_SYNTHESIS_FILE, &trolConfig->OutputSynthesisFile)){
    	if(printLog) printf("OutputSynthesisFile to apply: %s\n", trolConfig->OutputSynthesisFile);
	}
  	else
    	if(printLog) printf("No 'OutputSynthesisFile' setting in configuration file.\n");	

	//  noise to apply in all stoke parameter, only will be applied if sigma is empty
	if(config_lookup_float(&cfg, NOISE_FILE, &trolConfig->noise)){
    	if(printLog) printf("noise to apply: %f\n", trolConfig->noise);
	}
  	else
    	if(printLog) printf("No 'noise' setting in configuration file.\n");

	//  sigma to apply in each stokes parameter 
	if(config_lookup_float(&cfg, SIGMA_FILE, trolConfig->sigma)){
    	if(printLog) printf("Sigma to apply\n");
	}
  	else{
   	if(printLog) printf("No 'Sigma' setting in configuration file. Assigning noise to each parameter of the vector. \n");			
		trolConfig->sigma[0] = trolConfig->noise;
		trolConfig->sigma[1] = trolConfig->noise;
		trolConfig->sigma[2] = trolConfig->noise;
		trolConfig->sigma[3] = trolConfig->noise;
	}

		

	//  Optional minimum relative difference between two succesive merit-function values
	if(config_lookup_float(&cfg, TOPLIM_FILE, &trolConfig->toplim)){
    	if(printLog) printf("toplim to apply: %e\n", trolConfig->toplim);
	}
  	else
    	if(printLog) printf("No 'toplim' setting in configuration file.\n");					

	return 1;

}

/**
 * Read Cuantic data from a file with the Cuantic Lines. 
 * 
 * 
 * */
int readFileCuanticLines(const char * inputLineFile, PRECISION * cuanticDat, PRECISION centralLambda, int printLog){
	// try open the file with the 
	FILE * fp;
	char * line = NULL;

   size_t len = 0;
   ssize_t read;
	char atomo [2];
	fp = fopen(inputLineFile, "r");
   if(fp == NULL)	return 0;

	int indexLine, ionicState;
	int found = 0;
	double damping, potentialExcitation, logGf;
	PRECISION lambdaLine;
	int SLOI, SUPI;
	PRECISION LLOI,JLOI,LUPI,JUPI;
	char levelD,levelU;

	while ((read = getline(&line, &len, fp)) != -1  && !found ) {
		sscanf(line,"%i=%s %i %lf %lf %lf -%lf %i%c %lf- %i%c %lf",&indexLine,atomo,&ionicState,&lambdaLine,&damping,&potentialExcitation,&logGf,&SLOI,&levelD,&JLOI,&SUPI,&levelU,&JUPI);
		if(lambdaLine==centralLambda){ // read the rest of the line, else read next line
			switch (levelD)
			{
			case 'S':
				LLOI = 0;
				break;
			case 'P':
				LLOI = 1;
				break;
			case 'D':
				LLOI = 2;
				break;				
			case 'F':
				LLOI = 3;
				break;
			case 'G':
				LLOI = 4;
				break;				
			case 'H':
				LLOI = 5;
				break;
			case 'J':
				LLOI = 6;
				break;
			default:
				break;
			}
			switch (levelU)
			{
			case 'S':
				LUPI = 0;
				break;
			case 'P':
				LUPI = 1;
				break;
			case 'D':
				LUPI = 2;
				break;				
			case 'F':
				LUPI = 3;
				break;
			case 'G':
				LUPI = 4;
				break;				
			case 'H':
				LUPI = 5;
				break;
			case 'J':
				LUPI = 6;
				break;
			default:
				break;
			}
			if(SLOI==5) SLOI=2;
			if(SUPI==5) SUPI=2;

			found = 1; 
		}
   }
	cuanticDat[0] =1 ; // LINE NUMBER 1
	cuanticDat[1] = SLOI;
	cuanticDat[2] = LLOI;
	cuanticDat[3] = JLOI;
	cuanticDat[4] = SUPI;
	cuanticDat[5] = LUPI;
	cuanticDat[6] = JUPI;
	
	if(printLog){
		printf("\n***********************************");
		printf("\n\n CUANTIC NUMBERS READ FROM FILE: \n");
		printf("\n\tSLOI: %fd",cuanticDat[1]);
		printf("\n\tLLOI: %fd",cuanticDat[2]);
		printf("\n\tJLOI: %fd",cuanticDat[3]);
		printf("\n\tSUPI: %fd",cuanticDat[4]);
		printf("\n\tLUPI: %fd",cuanticDat[5]);
		printf("\n\tJUPI: %fd",cuanticDat[6]);
		printf("\n\n***********************************");
	}
	if(!found)
		return found;
	else
		return indexLine;

}


int readInitialModel(Init_Model * INIT_MODEL,const char * fileInitModel){
	
	config_t cfg;
  	
	config_init(&cfg);
	if(! config_read_file(&cfg, fileInitModel))
	{
		printf("%s:%d - %s\n", config_error_file(&cfg),config_error_line(&cfg), config_error_text(&cfg));
		config_destroy(&cfg);
		return(EXIT_FAILURE);
	}

	//  ETA0
	if(!config_lookup_float(&cfg, INITIAL_MODEL_ETHA0, &INIT_MODEL->eta0))
    	printf( "No 'ETA0' setting in configuration file. Used by default \n");

	//  B
	if(!config_lookup_float(&cfg, INITIAL_MODEL_B, &INIT_MODEL->B))
    	printf("No 'B' setting in configuration file. Used by default \n");

	//  vlos
	if(!config_lookup_float(&cfg, INITIAL_MODEL_VLOS, &INIT_MODEL->vlos))
    	printf("No 'VLOS' setting in configuration file. Used by default \n");

	//  LAMBDADOPP
	if(!config_lookup_float(&cfg, INITIAL_MODEL_LAMBDADOPP, &INIT_MODEL->dopp))
    	printf("No 'LAMBDADOPP' setting in configuration file. Used by default \n");

	//  AA
	if(!config_lookup_float(&cfg, INITIAL_MODEL_AA, &INIT_MODEL->aa))
    	printf("No 'AA' setting in configuration file. Used by default \n");		

	//  GM
	if(!config_lookup_float(&cfg, INITIAL_MODEL_GM, &INIT_MODEL->gm))
    	printf("No 'GM' setting in configuration file. Used by default \n");		

	//  AZI
	if(!config_lookup_float(&cfg, INITIAL_MODEL_AZI, &INIT_MODEL->az))
    	printf("No 'AZI' setting in configuration file. Used by default \n");		

	//  S0
	if(!config_lookup_float(&cfg, INITIAL_MODEL_S0, &INIT_MODEL->S0))
    	printf("No 'S0' setting in configuration file. Used by default \n");				

	//  S1
	if(!config_lookup_float(&cfg, INITIAL_MODEL_S1, &INIT_MODEL->S1))
    	printf("No 'S1' setting in configuration file. Used by default \n");				

	//  MAC
	if(!config_lookup_float(&cfg, INITIAL_MODEL_MAC, &INIT_MODEL->mac))
    	printf("No 'MAC' setting in configuration file. Used by default \n");				

	//  ALFA
	if(!config_lookup_float(&cfg, INITIAL_MODEL_ALFA, &INIT_MODEL->alfa))
    	printf("No 'ALFA' setting in configuration file. Used by default \n");						

	config_destroy(&cfg);

	return EXIT_SUCCESS;
}


/**
 * 
 * */
int readPSFFile(PRECISION * deltaLambda, PRECISION * PSF, const char * nameInputPSF){

	// first of all read the number of lines to give size for arrays deltaLambda and PSF
	FILE *fp;

	// alloc memory 

	char * line = NULL;
	size_t len = 0;
   ssize_t read;
	fp=fopen(nameInputPSF,"r");
	if(fp==NULL)
	{
	printf("File \"%s\" does not exist!!!\n",nameInputPSF);
			return 0;
	}	
	int index =0;
	while ((read = getline(&line, &len, fp)) != -1) {
		double delta, psf;
		sscanf(line,"%lf  %lf", &delta, &psf);
		deltaLambda[index] = delta;
		PSF[index] = psf;
		index++;
	}

	fclose(fp);
	return 1;
}


void loadInitialValues(ConfigControl * configControlFile){

	// array of weight 
	configControlFile->WeightForStokes[0]=1.;
	configControlFile->WeightForStokes[1]=1.;
	configControlFile->WeightForStokes[2]=1.;
	configControlFile->WeightForStokes[3]=1.;

	//INIT_MODEL=[eta0,magnet,vlos,landadopp,aa,gamma,azi,B1,B2,macro,alfa]
	int i;
	for(i=0;i<10;i++){
		configControlFile->fix[i]= 1;
	}
	configControlFile->fix[10] = 0;

	configControlFile->noise = NOISE_SIGMA;
	configControlFile->sigma[0] = NOISE_SIGMA;
	configControlFile->sigma[1] = NOISE_SIGMA;
	configControlFile->sigma[2] = NOISE_SIGMA;
	configControlFile->sigma[3] = NOISE_SIGMA;

	configControlFile->ToleranceForSVD = 1e-4;
	configControlFile->InitialDiagonalElement = ILAMBDA;
	configControlFile->toplim = TOPLIM;
	configControlFile->mu = AH;
	configControlFile->TypeConvolution = "FFT";

}