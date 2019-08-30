#include <string.h>
#include <stdio.h>
#include "fitsio.h"
#include "utilsFits.h"
#include "defines.h"

int main(int argc, char **argv){

   FitsImage * fitsImage;

   printf("\n NOMBRE DEL FICHERO FITS: %s", argv[1]);
   printf("\n lambda file: %s", argv[2]);
   printf("\n fichero a escribir %s", argv[3]);
   printf("\n**********************************");
   fitsImage = readFitsSpectroImage(argv[1]);
   printf("\n Número de pixeles: %d: ", fitsImage->numPixels);   
   printf("\n**********************************");
   
   writeFitsImageProfiles(argv[3],argv[1],fitsImage);

   return 1;
}
