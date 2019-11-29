#include "readConfig.h"
#include "defines.h"
#include <time.h>
#include <string.h>
#include <stdio.h>
#include <stddef.h>
#include <locale.h>
#include <stdlib.h>
//#include <libconfig.h>



/*int readConfigControl(char * configFile, ConfigControl * trolConfig, int printLog){


	config_t cfg;
	config_init(&cfg);
	  
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

}*/

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
			SLOI= (SLOI-1)/2;
			SUPI= (SUPI-1)/2;

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
		
	}
	if(!found)
		return found;
	else
		return indexLine;

}


int readInitialModel(Init_Model * INIT_MODEL, char * fileInitModel){
	
	FILE * fReadInitModel;
	char * line = NULL;
	size_t len = 0;
   ssize_t read;
	char comment[200], name[100];
	int rfscanf;
	fReadInitModel = fopen(fileInitModel, "r");
	if (fReadInitModel == NULL)
	{
		printf("Error opening the file of parameters, it's possible that the file doesn't exist. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileInitModel);
		fclose(fReadInitModel);
		return 0;
	}
	
	while ((read = getline(&line, &len, fReadInitModel)) != -1) {
		double aux_value;
		rfscanf = sscanf(line,"%99[^:]:%lf%99[^!]!",name, &aux_value,comment);
		if(strstr(name,"INITIAL_MODEL_B")!=NULL){ // B
			INIT_MODEL->B = aux_value;
		}
		if(strstr(name,"INITIAL_MODEL_GM")!=NULL){ // GM
			INIT_MODEL->gm = aux_value;
		}
		if(strstr(name,"INITIAL_MODEL_AZI")!=NULL){ // AZI
			INIT_MODEL->az = aux_value;
		}
		if(strstr(name,"INITIAL_MODEL_ETHA0")!=NULL){ // ETHA0
			INIT_MODEL->eta0 = aux_value;
		}
		if(strstr(name,"INITIAL_MODEL_LAMBDADOPP")!=NULL){ // LAMBDADOPP
			INIT_MODEL->dopp = aux_value;
		}
		if(strstr(name,"INITIAL_MODEL_AA")!=NULL){ // AA
			INIT_MODEL->aa = aux_value;
		}
		if(strstr(name,"INITIAL_MODEL_ALFA")!=NULL){ // ALFA
			INIT_MODEL->alfa = aux_value;
		}
		if(strstr(name,"INITIAL_MODEL_MAC")!=NULL){ // MAC
			INIT_MODEL->mac = aux_value;
		}		
		if(strstr(name,"INITIAL_MODEL_VLOS")!=NULL){ // VLOS
			INIT_MODEL->vlos = aux_value;
		}
		if(strstr(name,"INITIAL_MODEL_S0")!=NULL){ // S0
			INIT_MODEL->S0 = aux_value;
		}
		if(strstr(name,"INITIAL_MODEL_S1")!=NULL){ // S1
			INIT_MODEL->S1 = aux_value;
		}				
	}
	fclose(fReadInitModel);
	return 1;
}

/*int readInitialModel(Init_Model * INIT_MODEL,const char * fileInitModel){
	
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
*/

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
	

}

int readParametersFileInput(char * fileParameters,  ConfigControl * trolConfig, int printLog){

	// try open the file with the 
	FILE * fReadParameters;
	char LINE [4096], * returnLine;
	char comment[200], name[100];
	fReadParameters = fopen(fileParameters, "r");
	if (fReadParameters == NULL)
	{
		printf("Error opening the file of parameters, it's possible that the file doesn't exist. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;
	}
	int rfscanf; 
	
	/***************************  NUMBER OF CYCLES  ********************************************/

	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->NumberOfCycles,comment);
	
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param NumberOfCycles. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if( trolConfig->NumberOfCycles < 0){
		printf("milos: Error in NumberOfCycles parameter. review it. Not accepted: %d\n", trolConfig->NumberOfCycles);
		return 0;
	}
	if(printLog) printf("NumberOfCycles to apply: %i\n", trolConfig->NumberOfCycles);

	/***************************  OBSERVED PROFILES  ********************************************/

	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;
	rfscanf = sscanf(LINE,"%99[^:]:%s%99[^!]!",name, trolConfig->ObservedProfiles,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Observed Profiles. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Observed profiles to apply: %s\n", trolConfig->ObservedProfiles);

	/***************************  STRAY LIGHT FILE ********************************************/

	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;
	rfscanf = sscanf(LINE,"%99[^:]:%s%99[^!]!",name, trolConfig->StrayLightFile,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Stray light file. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Stray light file to apply: %s\n", trolConfig->StrayLightFile);


	/***************************  PSF FILE ********************************************/

	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;
	rfscanf = sscanf(LINE,"%99[^:]:%s%99[^!]!",name, trolConfig->PSFFile,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param PSF file. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("PSF file to apply: %s\n", trolConfig->PSFFile);

	/*************************** WAVELENGHT FILE ********************************************/

	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;
	rfscanf = sscanf(LINE,"%99[^:]:%s%99[^!]!",name, trolConfig->WavelengthFile,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param wavelength file. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("wavelength file to apply: %s\n", trolConfig->WavelengthFile);

	/*************************** ATOMIC PARAMETER  FILE ********************************************/

	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;
	rfscanf = sscanf(LINE,"%99[^:]:%s%99[^!]!",name, trolConfig->AtomicParametersFile,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Atomic parameters file. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Atomic parameters file to apply: %s\n", trolConfig->AtomicParametersFile);

	/*************************** INITIAL GUESS MODEL   FILE ********************************************/

	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;
	rfscanf = sscanf(LINE,"%99[^:]:%s%99[^!]!",name, trolConfig->InitialGuessModel,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Initial guess model 1. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Initial guess model 1 to apply: %s\n", trolConfig->InitialGuessModel);

	/*************************** INITIAL GUESS MODEL  2  FILE ********************************************/

	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;
	rfscanf = sscanf(LINE,"%99[^:]:%s%99[^!]!",name, trolConfig->InitialGuessModel_2,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Initial guess model 2. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Initial guess model 2 to apply: %s\n", trolConfig->InitialGuessModel_2);

	/*************************** WEIGHT FOT STOKES I ********************************************/

	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;
	rfscanf = sscanf(LINE,"%99[^:]:%lf%99[^!]!",name, &trolConfig->WeightForStokes[0],comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Weight for Stokes I. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Weight for Stokes I to apply: %lf\n", trolConfig->WeightForStokes[0]);

	/*************************** WEIGHT FOT STOKES Q ********************************************/

	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;
	rfscanf = sscanf(LINE,"%99[^:]:%lf%99[^!]!",name, &trolConfig->WeightForStokes[1],comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Weight for Stokes Q. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Weight for Stokes Q to apply: %lf\n", trolConfig->WeightForStokes[1]);

	/*************************** WEIGHT FOT STOKES U ********************************************/

	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;
	rfscanf = sscanf(LINE,"%99[^:]:%lf%99[^!]!",name, &trolConfig->WeightForStokes[2],comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Weight for Stokes U. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Weight for Stokes U to apply: %lf\n", trolConfig->WeightForStokes[2]);

	/*************************** WEIGHT FOT STOKES V ********************************************/

	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;
	rfscanf = sscanf(LINE,"%99[^:]:%lf%99[^!]!",name, &trolConfig->WeightForStokes[3],comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Weight for Stokes V. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Weight for Stokes V to apply: %lf\n", trolConfig->WeightForStokes[3]);


	/*************************** MU ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%lf%99[^!]!",name, &trolConfig->mu,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param mu=cos (theta). Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("mu=cos (theta) to apply: %f\n", trolConfig->mu);

	/*************************** EstimatedSNForI ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->EstimatedSNForI,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Estimated S/N for I. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog){
		printf("Estimated S/N for I  to apply: %i\n", trolConfig->EstimatedSNForI);
		printf("Estimated noise for I  to apply: %lf\n", 1.0/trolConfig->EstimatedSNForI);
	} 
	trolConfig->noise = 1.0/trolConfig->EstimatedSNForI;
	trolConfig->sigma[0] = trolConfig->noise*trolConfig->noise;
	trolConfig->sigma[1] = trolConfig->sigma[0];
	trolConfig->sigma[2] = trolConfig->sigma[0];
	trolConfig->sigma[3] = trolConfig->sigma[0];
	// PUT VALUES IN ARRAY OF SIGMA 


 	/*************************** CONTINIUM CONTRAST ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->ContinuumContrast,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Continuum contrast. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Continuum contrast to apply: %i\n", trolConfig->ContinuumContrast);

	/*************************** TOLERANCE_FOR_SVD ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%lf%99[^!]!",name, &trolConfig->ToleranceForSVD,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Tolerance for SVD. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Tolerance for SVD to apply: %le\n", trolConfig->ToleranceForSVD);

	/*************************** INITIAL_DIAGONAL_ELEMENT ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%lf%99[^!]!",name, &trolConfig->InitialDiagonalElement,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Initial diagonal element. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Initial diagonal element  to apply: %le\n", trolConfig->InitialDiagonalElement);

	/*************************** USE_INTERPOLAR_SPLINES_OR_LINEAR ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->useInterpolarSplinesOrLinear,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Splines/Linear Interpolation. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Splines/Linear Interpolation to apply: %d\n", trolConfig->useInterpolarSplinesOrLinear);

	/*************************** USE PSF FILTER GAUSSIAN ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->ConvolveWithPSF,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Use Gaussian PSF . Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Use Gaussian PSF  to apply: %d\n", trolConfig->ConvolveWithPSF);

	/*************************** FWHM ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%lf%99[^!]!",name, &trolConfig->FWHM,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param FWHM. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("FWHM to apply: %f\n", trolConfig->FWHM);

	/*************************** NTL ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->ntl,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param NTL. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("NTL to apply: %d\n", trolConfig->ntl);

	/*************************** NLIOBS ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->nliobs,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param NLIOBS. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("NLIOBS to apply: %d\n", trolConfig->nliobs);


	/*************************** CENTRAL_WAVE_LENGHT ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%lf%99[^!]!",name, &trolConfig->CentralWaveLenght,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Central Wave Lenght. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Central Wave Lenght  to apply: %f\n", trolConfig->CentralWaveLenght);

	/*************************** ETA0_LINE_TO_CONTINUUM_ABSORPTION ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->fix[0],comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Line To Continiuum Absorption. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Line To Continiuum Absorption to apply: %d\n", trolConfig->fix[0]);

	/*************************** B_MAGNETIC_FIELD_STRENGTH ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->fix[1],comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Magnetic Field Strength. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Magnetic Field Strength to apply: %d\n", trolConfig->fix[1]);

	/*************************** VLOS_LINE_OF_SIGHT_VELOCITY ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->fix[2],comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Line Of Sight Velocity. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Line Of Sight Velocity to apply: %d\n", trolConfig->fix[2]);

	/*************************** DOPP_DOOPLER_WIDTH ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->fix[3],comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Doopler Width. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Doopler Width to apply: %d\n", trolConfig->fix[3]);

	/*************************** AA_DAMPING_PARAMETER ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->fix[4],comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Damping Parameter. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Damping Parameter to apply: %d\n", trolConfig->fix[4]);

	/*************************** GM_MAGNETIC_FIELD_INCLINATION ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->fix[5],comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Magnetic FieldInclination. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Magnetic FieldInclination to apply: %d\n", trolConfig->fix[5]);

	/*************************** AZ_MAGNETIC_FIELD_AZIMUTH ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->fix[6],comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Magnetic Field Azimuth. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Magnetic Field Azimuth to apply: %d\n", trolConfig->fix[6]);

	/*************************** S0_SOURCE_FUNCTION_CONSTANT ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->fix[7],comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Source Function Constant. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Source Function Constant to apply: %d\n", trolConfig->fix[7]);

	/*************************** S1_SOURCE_FUNCTION_GRADIENT ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->fix[8],comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Source Function Gradient . Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Source Function Gradient  to apply: %d\n", trolConfig->fix[8]);

	/*************************** MAC_MACROTURBULENT_VELOCITY ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->fix[9],comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Macroturbulent Velocity. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Macroturbulent Velocity to apply: %d\n", trolConfig->fix[9]);

	/*************************** ALPHA_FILLING_FACTOR ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->fix[10],comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Filling Factor. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Filling Factor  to apply: %d\n", trolConfig->fix[10]);

	/*************************** SAVE_CHISQR ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->saveChisqr,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Save Chisqr. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Save Chisqr  to apply: %d\n", trolConfig->saveChisqr);

	/*************************** USE_CLASSICAL_ESTIMATES ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->UseClassicalEstimates,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Use Classical Estimates. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Use Classical Estimates to apply: %d\n", trolConfig->UseClassicalEstimates);

	/*************************** USE_RTE_INVERSION ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->UseRTEInversion,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Use RTE Inversion. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Use RTE Inversion to apply: %d\n", trolConfig->UseRTEInversion);


	/*************************** SAVE SYNTHESIS PROFILE ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, &trolConfig->SaveSynthesisProfile,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Save Synthesis Profile. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Save Synthesis Profile to apply: %d\n", trolConfig->SaveSynthesisProfile);

	/*************************** OUTPUT_MODEL_FILE ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%s%99[^!]!",name, trolConfig->OutputModelFile,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Output Model File. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Output Model File  to apply: %s\n", trolConfig->OutputModelFile);

	/*************************** OUTPUT_SYNTHESIS_FILE ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%s%99[^!]!",name, trolConfig->OutputSynthesisFile,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param Output Synthesis File. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("Output Synthesis File to apply: %s\n", trolConfig->OutputSynthesisFile);

	/*************************** SIGMA_FILE ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%lf%99[^!]!",name, &trolConfig->noise,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param NOISE. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("NOISE to apply: %le\n", trolConfig->noise);


	/*************************** TOPLIM_FILE ********************************************/
	
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%lf%99[^!]!",name, &trolConfig->toplim,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param TOPLIM. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if(printLog) printf("TOPLIM to apply: %le\n", trolConfig->toplim);


	return 1;

}