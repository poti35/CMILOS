
#include "defines.h"

int readConfigControl(char * configFile, ConfigControl * trolConfig, int printLog);

int readFileCuanticLines(const char * inputLineFile, PRECISION * cuanticDat, PRECISION centralLambda, int printLog);

int readInitialModel(Init_Model * INIT_MODEL, const char * fileInitModel);

int readPSFFile(PRECISION * deltaLambda, PRECISION * PSF, const char * nameInputPSF);

void loadInitialValues(ConfigControl * configControlFile);

