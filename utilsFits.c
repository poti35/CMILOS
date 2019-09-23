#include "utilsFits.h"
#include "fitsio.h"
#include <time.h>
#include <string.h>
#include <stdio.h>
#include <stddef.h>
#include <locale.h>

/**
 * image --> 	Pointer to store the image 
 * fitsFile --> name of file to read (included the path to the file)
 * return NULL in case of error, 1 if everything was fine
 */
vpixels * readImagePixels (char * fitsFile, int * numPixels){

   fitsfile *fptr;   /* FITS file pointer, defined in fitsio.h */
   int bitpix, naxis;
   long naxes[4] = {1,1,1,1}, fpixel[4] = {1,1,1,1};
	char card[FLEN_CARD];   /* Standard string lengths defined in fitsio.h */
   int status = 0;   /* CFITSIO status value MUST be initialized to zero! . Status 0 indicate that everything was fine. */
	int contPixels, currentLambda;
	int nLambda;
	vpixels * image;
	PRECISION * spectro;
	PRECISION * vlambda;
	
	printf("\n ***********************START READING FITS FILE***********************");

   // OPEN THE FITS FILE TO READ THE DEPTH OF EACH DIMENSION
   
	if (!fits_open_file(&fptr, fitsFile, READONLY, &status)){
	
	   printf("\n ***********************FITS FILE OPENED***********************");
   
	   // READ THE HDU PARAMETER FROM THE FITS FILE
   
	   int hdutype;
      fits_get_hdu_type(fptr, &hdutype, &status);

		// We want only fits image 
		if(hdutype==IMAGE_HDU){

			printf("\n ***********************FITS FILE TYPE IMAGE***********************");
			
			// READ IMAGE AND STORAGE IN STRUCTURE IMAGE 
			// TODO --> READ IN WITCH NAXES IS STORED LAMBDA
			nLambda = 100000;
			if (!fits_get_img_param(fptr, 4, &bitpix, &naxis, naxes, &status) ){
								
				// GET THE CURRENT POSITION OF EVERY PARAMETER
				int pos_lambda; 
				int pos_row;
				int pos_col;
				int pos_stokes_parameters; // assume that we have 4 stokes parameters
				int ii;
				// Stokes paramter position , assume that there are 4 stokes parameters. Positions for other we will assume other that can be wrong but it's not important por us. 
				for(ii = 0; ii< naxis; ii++){
					if(naxes[ii]==4) pos_stokes_parameters = ii;
					else{
						if(naxes[ii]<nLambda){
							nLambda = naxes[ii];
							pos_lambda = ii;
						}
					}
					printf("\n naxes %4d : %4ld", ii, naxes[ii]);
				}

				if( (pos_stokes_parameters == 0 && pos_lambda == 1) || (pos_stokes_parameters == 1 && pos_lambda == 0) ){ 
					pos_row = 2;
					pos_col = 3;
				}
				if( (pos_stokes_parameters == 0 && pos_lambda == 2) || (pos_stokes_parameters == 2 && pos_lambda == 0) ){ 
					pos_row = 1;
					pos_col = 3;
				} 
				if( (pos_stokes_parameters == 0 && pos_lambda == 3) || (pos_stokes_parameters == 3 && pos_lambda == 0) ){ 
					pos_row = 1;
					pos_col = 2;					
				}
				if( (pos_stokes_parameters == 1 && pos_lambda == 2) || (pos_stokes_parameters == 2 && pos_lambda == 1) ){ 
					pos_row = 0;
					pos_col = 3;					
				}
				if( (pos_stokes_parameters == 1 && pos_lambda == 3) || (pos_stokes_parameters == 3 && pos_lambda == 1) ){ 
					pos_row = 0;
					pos_col = 2;					
				}
				if( (pos_stokes_parameters == 2 && pos_lambda == 3) || (pos_stokes_parameters == 3 && pos_lambda == 2) ){ 
					pos_row = 0;
					pos_col = 1;					
				}				
				printf(" \n Pos stokes %4d , pos lambda %4d , pos row %4d, pos col %4d", pos_stokes_parameters, pos_lambda, pos_row, pos_col);
				printf("\n***************************");				
				// ALLOCATE MEMORY FOR STORE THE IMAGE. 
				printf(" \n Reservo memoria para %4ld x %4ld", naxes[pos_row], naxes[pos_col]);
				image = malloc ( ((naxes[pos_row]*naxes[pos_col])+1)*sizeof(vpixels));
				*numPixels = (naxes[pos_row]*naxes[pos_col]);
				
				contPixels = 0;
				
				for(fpixel[pos_row] = 1; (fpixel[pos_row] <= naxes[pos_row]) && contPixels<(*numPixels) ; fpixel[pos_row]++){
					for(fpixel[pos_col] = 1; (fpixel[pos_col] <= naxes[pos_col]) && contPixels<(*numPixels) ; fpixel[pos_col]++){
						
						//printf(" \n Reservo memoria para nLambda %4d", nLambda);
						spectro = calloc( 4 * nLambda, sizeof(PRECISION));
						vlambda = calloc(nLambda, sizeof(PRECISION));
						currentLambda = 0;
						for(fpixel[pos_lambda] = 1; fpixel[pos_lambda] <= naxes[pos_lambda]; fpixel[pos_lambda]++){
							vlambda[currentLambda] = 6173+ (currentLambda+1)*0.03;
							// allocat memory for every stoke parameter
							// READ naxes[4] elements, assume naxes[4] is 4 	
							
							for(fpixel[pos_stokes_parameters] = 1; fpixel[pos_stokes_parameters] <= naxes[pos_stokes_parameters]; fpixel[pos_stokes_parameters]++){
								//stokeParameter = malloc(sizeof(PRECISION));
								PRECISION stokeParameter;
								//printf(" \n Leyendo pixel [%4ld][%4ld][%4ld][%4ld]", fpixel[0],fpixel[1],fpixel[2],fpixel[3]);
								if (fits_read_pix(fptr, TFLOAT, fpixel, 1 , NULL, &stokeParameter , NULL, &status) ) 
									printf("\n PIXEL NO SE PUDE LEER");
								
								//printf(" \n Pixel leido: %f", stokeParameter);
								if(fpixel[pos_stokes_parameters] == 1) spectro[currentLambda] = stokeParameter;  // I
								if(fpixel[pos_stokes_parameters] == 2) spectro[currentLambda+ nLambda] = stokeParameter; // Q
								if(fpixel[pos_stokes_parameters] == 3) spectro[currentLambda+ (nLambda *2)] = stokeParameter; // U 
								if(fpixel[pos_stokes_parameters] == 4) spectro[currentLambda+ (nLambda *3)] = stokeParameter; // V
							}
							currentLambda++;
						}
						//printf("\n Escribiendo pixel: %4d", contPixels);
						image[contPixels].nLambda = nLambda;
						image[contPixels].spectro = spectro;
						image[contPixels].vLambda = vlambda;
						contPixels++;
					}
				}
				printf(" \n  IMAGEN LEIDA ");
				printf("\n*************");
			}


		}
		else{
			return NULL;  // we are interested only in FITS image
		}
	}
   else{ // IN CASE AN ERROR OPENING THE FILE RETURN THE ERROR CODE
      if (status) fits_report_error(stderr, status); /* print any error message */
      return NULL;
   }

	fits_close_file(fptr, &status);
	return image; 
}


/**
 * Clean memory from fits image
 */
void freeVpixels(vpixels * image, int numPixels){
	int i;
	for(i=9;i<numPixels;i++){
		free(image[i].spectro);
		free(image[i].vLambda);
	}
	free(image);
}


/**
 * This function read the spectro image from the file "fitsFileSpectra" and store it into a struct of FitsImage
 * fitsFileSpectra --> name of the fits file to read 
 * Return the image read or NULL if something was wrong during the lecture. 
 */
FitsImage *  readFitsSpectroImage (char * fitsFileSpectra){
   fitsfile *fptr;   /* FITS file pointer, defined in fitsio.h */
	FitsImage * image =  malloc(sizeof(FitsImage));
   int status = 0;   /* CFITSIO status value MUST be initialized to zero! */
	PRECISION bscale = 1.0, bzero = 0.0, nulval = 0.; // define null value to 0 because the performance to read from fits file is better doing this. 
   int bitpix, bytepix, naxis, anynul,first, numPixelsFitsFile;
   long naxes [4] = {1,1,1,1}; /* The maximun number of dimension that we will read is 4*/
	char comment[FLEN_CARD];   /* Standard string lengths defined in fitsio.h */
	int ii,kk, hh;
	int i, j, k, h;
   // OPEN THE FITS FILE TO READ THE DEPTH OF EACH DIMENSION
   if (!fits_open_file(&fptr, fitsFileSpectra, READONLY, &status)){
      // READ THE HDU PARAMETER FROM THE FITS FILE
      int hdutype;
      fits_get_hdu_type(fptr, &hdutype, &status);

		// We want only fits image 
		if(hdutype==IMAGE_HDU){
			// We assume that we have only on HDU as primary 
			if(fits_read_key(fptr, TSTRING, CTYPE1, image->ctype_1, comment, &status)) return 0;
			if(fits_read_key(fptr, TSTRING, CTYPE2, image->ctype_2, comment, &status)) return 0;
			if(fits_read_key(fptr, TSTRING, CTYPE3, image->ctype_3, comment, &status)) return 0;
			if(fits_read_key(fptr, TSTRING, CTYPE4, image->ctype_4, comment, &status)) return 0;
			/* if(fits_read_card(fptr, CUNIT1, image->cunit_1, &status)) return 0;
			if(fits_read_card(fptr, CUNIT2, image->cunit_2, &status)) return 0;
			if(fits_read_card(fptr, CUNIT3, image->cunit_3, &status)) return 0;*/

			// GET THE CURRENT POSITION OF EVERY PARAMETER
			int pos_lambda; 
			int pos_row;
			int pos_col;
			int pos_stokes_parameters;
			// LAMBDA POSITION
			if(strcmp(image->ctype_1,CTYPE_WAVE)==0) pos_lambda = 0;
			if(strcmp(image->ctype_2,CTYPE_WAVE)==0) pos_lambda = 1;
			if(strcmp(image->ctype_3,CTYPE_WAVE)==0) pos_lambda = 2;
			if(strcmp(image->ctype_4,CTYPE_WAVE)==0) pos_lambda = 3;

			// HPLN TAN 
			if(strcmp(image->ctype_1,CTYPE_HPLN_TAN)==0) pos_row = 0;
			if(strcmp(image->ctype_2,CTYPE_HPLN_TAN)==0) pos_row = 1;
			if(strcmp(image->ctype_3,CTYPE_HPLN_TAN)==0) pos_row = 2;
			if(strcmp(image->ctype_4,CTYPE_HPLN_TAN)==0) pos_row = 3;

			// HPLT TAN 
			if(strcmp(image->ctype_1,CTYPE_HPLT_TAN)==0) pos_col = 0;
			if(strcmp(image->ctype_2,CTYPE_HPLT_TAN)==0) pos_col = 1;
			if(strcmp(image->ctype_3,CTYPE_HPLT_TAN)==0) pos_col = 2;
			if(strcmp(image->ctype_4,CTYPE_HPLT_TAN)==0) pos_col = 3;			

			// Stokes paramter position , 
			if(strcmp(image->ctype_1,CTYPE_STOKES)==0) pos_stokes_parameters = 0;
			if(strcmp(image->ctype_2,CTYPE_STOKES)==0) pos_stokes_parameters = 1;
			if(strcmp(image->ctype_3,CTYPE_STOKES)==0) pos_stokes_parameters = 2;
			if(strcmp(image->ctype_4,CTYPE_STOKES)==0) pos_stokes_parameters = 3;

			// READ IMAGE AND STORAGE IN STRUCTURE IMAGE 
			if (!fits_get_img_param(fptr, 4, &bitpix, &naxis, naxes, &status) ){

				int datatype = 0;
				switch(bitpix) {
              case BYTE_IMG:
                  datatype = TBYTE;
                  break;
              case SHORT_IMG:
                  datatype = TSHORT;
                  break;
              case LONG_IMG:
                  datatype = TINT;
                  break;
              case FLOAT_IMG:
                  datatype = TFLOAT;
                  break;
              case DOUBLE_IMG:
                  datatype = TDOUBLE;
                  break;
          	}
				image->rows=naxes[pos_row];
				image->cols=naxes[pos_col];
				image->nLambdas=naxes[pos_lambda];
				image->numStokes=naxes[pos_stokes_parameters];
				image->numPixels = naxes[pos_col] * naxes[pos_row]; // we will read the image by columns 
				image->pos_lambda = pos_lambda;
				image->pos_col = pos_col;
				image->pos_row = pos_row;
				image->pos_stokes_parameters = pos_stokes_parameters;
				numPixelsFitsFile = naxes[pos_row]*naxes[pos_col]*naxes[pos_lambda]*naxes[pos_stokes_parameters];
				printf("\n NÚMERO DE PIXELES EN LA IMAGEN %d", numPixelsFitsFile);
				printf("\n**********************");
				// allocate memory to read all pixels in the same array 
				float * imageTemp = calloc(numPixelsFitsFile, sizeof(float));
				if (!imageTemp)  {
					printf("ERROR ALLOCATION MEMORY FOR TEMP IMAGE");
					return NULL;
          	}
				first = 1;
				//fits_read_img(fptr, datatype, first, numPixelsFitsFile, &nulval, imageTemp, &anynul, &status);
				long fpixel [4] = {1,1,1,1}; 
				fits_read_pix(fptr, datatype, fpixel, numPixelsFitsFile, &nulval, imageTemp, &anynul, &status);
				if(status){
					fits_report_error(stderr, status);
               return NULL;	
				}

				// allocate memory for reorder the image
				image->pixels = calloc(image->numPixels, sizeof(vpixels));
				image->vLambdaImagen = calloc(image->numPixels*image->nLambdas, sizeof(PRECISION));
				image->spectroImagen = calloc(image->numPixels*image->nLambdas*image->numStokes, sizeof(PRECISION));
				printf("\n Número de pixeles: %d", image->numPixels);
				printf("\n ***********************************************");
				for( i=0;i<image->numPixels;i++){
					image->pixels[i].spectro = calloc ((image->numStokes*image->nLambdas),sizeof(PRECISION));
					image->pixels[i].vLambda = calloc (image->nLambdas, sizeof(PRECISION));
					image->pixels[i].nLambda = image->nLambdas;
				}
				int currentLambda = 0, currentRow = 0, currentStokeParameter=0, currentCol = 0, currentPixel;
				//PRECISION pixel;
				if(naxis==4){ // image with 4 dimension 
					
					for( i=0; i<naxes[3];i++){
						for( j=0; j<naxes[2];j++){
							for( k=0;k<naxes[1];k++){
								for( h=0;h<naxes[0];h++){
//					for( fpixel[3] = 1; fpixel[3]<=naxes[3];fpixel[3]++){
//						for( fpixel[2] = 1; fpixel[2]<=naxes[2];fpixel[2]++){
//							for( fpixel[1] = 1; fpixel[1]<=naxes[1];fpixel[1]++){
//								for( fpixel[0] = 1; fpixel[0]<=naxes[0]; fpixel[0]++){	
									float pixel = 0.0;
									//fits_read_pix(fptr, datatype, fpixel, 1, &nulval, &pixel, &anynul, &status);
									// I NEED TO KNOW THE CURRENT POSITION OF EACH ITERATOR 
									switch (pos_lambda)
									{
										case 0:
											currentLambda = h;
											//currentLambda = fpixel[0]-1;
											break;
										case 1:
											currentLambda = k;
											//currentLambda = fpixel[1]-1;
											break;
										case 2:
											currentLambda = j;
											//currentLambda = fpixel[2]-1;
											break;
										case 3:
											currentLambda = i;
											//currentLambda = fpixel[3]-1;
											break;																						
									}
									switch (pos_stokes_parameters)
									{
										case 0:
											currentStokeParameter = h;
											//currentStokeParameter = fpixel[0]-1;
											break;
										case 1:
											currentStokeParameter = k;
											//currentStokeParameter = fpixel[1]-1;
											break;
										case 2:
											currentStokeParameter = j;
											//currentStokeParameter = fpixel[2]-1;
											break;
										case 3:
											currentStokeParameter = i;
											//currentStokeParameter = fpixel[3]-1;
											break;																						
									}
									switch (pos_row)
									{
										case 0:
											currentRow = h;
											//currentRow = fpixel[0]-1;
											break;
										case 1:
											currentRow = k;
											//currentRow = fpixel[1]-1;
											break;
										case 2:
											currentRow = j;
											//currentRow = fpixel[2]-1;
											break;
										case 3:
											currentRow = i;
											//currentRow = fpixel[3]-1;
											break;																						
									}
									switch (pos_col)
									{
										case 0:
											currentCol = h;
											//currentCol = fpixel[0]-1;
											break;
										case 1:
											currentCol = k;
											//currentCol = fpixel[1]-1;;
											break;
										case 2:
											currentCol = j;
											//currentCol = fpixel[2]-1;
											break;
										case 3:
											currentCol = i;
											//currentCol = fpixel[3]-1;
											break;																						
									}			
									pixel = imageTemp [(i*naxes[2]*naxes[1]*naxes[0]) + (j*naxes[1]*naxes[0]) + (k*naxes[0]) + h];
									currentPixel = (currentCol*naxes[pos_row]) + currentRow;
									//currentPixel = (currentRow*naxes[pos_col]) + currentCol;
									//printf("\n CURRENTLAMBDA %d CURRENTSTOKEPARAMETER %d CURRENTROW %d CURRENTCOL %d NUMERO DE SPECTRO %d NÚMERO DE ITER %d -- NUMERO DE PIXEL -- %d  VALOR PIXEL: %lf",currentLambda, currentStokeParameter,currentRow, currentCol,(currentLambda + image->nLambdas * currentStokeParameter), numiter, currentPixel, pixel);									
									//printf("\n %d %d %d %d %d %d %d %lf",currentLambda, currentStokeParameter,currentRow, currentCol,(currentLambda + image->nLambdas * currentStokeParameter), numiter, currentPixel, pixel);									
									//printf("\n*");
									image->pixels[currentPixel].spectro[currentLambda + (image->nLambdas * currentStokeParameter)] = pixel;  // I =0, Q = 1, U = 2, V = 3
								}
							}
						}
					}
				}
				int contSpectro = 0;
				for( i=0;i<image->numPixels;i++){
					for( j=0;j<(image->nLambdas*image->numStokes);j++){
						image->spectroImagen[contSpectro++] = image->pixels[i].spectro[j];
					}
				}
				printf("\n IMAGEN LEIDA size spectro %d ", contSpectro);
				printf("**********");
				free(imageTemp);
				fits_close_file(fptr, &status);
				if (status){
					fits_report_error(stderr, status);
					return NULL;
				}				
			}
		}
		else{
			return NULL;  // we are interested only in FITS image
		}
	}
   else{ // IN CASE AN ERROR OPENING THE FILE RETURN THE ERROR CODE
      if (status) fits_report_error(stderr, status); /* print any error message */
      return NULL;
   }
	
	return image; 
}


/**
 * This function read the lambda values for the image from the file "fitsFileLambda" and store it into a struct of FitsImage. The file of spectro must
 * be read it before call this method. 
 * fitsFileLambda --> name of the fits file to read with lambda values 
 * fitsImage --> struct of image 
 * Return 1 If the image has been read corectly if not return 0 
 */

int  readFitsLambdaFile (char * fitsFileLambda, FitsImage * fitsImage){
	int i, j, k;
	fitsfile *fptr;   /* FITS file pointer, defined in fitsio.h */
	int status = 0;   /* CFITSIO status value MUST be initialized to zero! */
	PRECISION bscale = 1.0, bzero = 0.0, nulval = 0.; // define null value to 0 because the performance to read from fits file is better doing this. 
	int bitpix, bytepix, naxis, anynul,first, numPixelsFitsFile;
	long naxes [4] = {1,1,1,1}; /* The maximun number of dimension that we will read is 4*/
	char comment[FLEN_CARD];   /* Standard string lengths defined in fitsio.h */
	printf("\n READING IMAGE WITH LAMBDA ");
	printf("\n**********");
	if (!fits_open_file(&fptr, fitsFileLambda, READONLY, &status)){
		if (!fits_get_img_param(fptr, 4, &bitpix, &naxis, naxes, &status) ){
			int datatype = 0;
			switch(bitpix) {
				case BYTE_IMG:
					datatype = TBYTE;
					break;
				case SHORT_IMG:
					datatype = TSHORT;
					break;
				case LONG_IMG:
					datatype = TINT;
					break;
				case FLOAT_IMG:
					datatype = TFLOAT;
					break;
				case DOUBLE_IMG:
					datatype = TDOUBLE;
					break;
			}

			if(naxis!=1  || naxis!=3){
				if(naxis == 1){ // array of lambads 
					if(naxes[0]!=fitsImage->nLambdas){ // image of lambas has different size of spectra image 
						printf("\n IMAGE OF LAMBAS HAS DIFFERENT SIZE OF SPECTRA IMAGE. NUMBER OF LAMBDAS IN SPECTRA IMAGE %d NUMBER OF LAMBDA IMAGE %ld ", fitsImage->nLambdas,naxes[0]);
						freeFitsImage(fitsImage);
						return 0;
					}
					PRECISION  * vAuxLambdas = calloc(naxes[0], sizeof(PRECISION));
					long fpixel [1] = {1};
					fits_read_pix(fptr, TDOUBLE, fpixel, naxes[0], &nulval, vAuxLambdas, &anynul, &status) ;
					printf("\n*********************");
					if(status){
						fits_report_error(stderr, status);
						return 0;	
					}
					int contLambda = 0;
					for( i=0;i<fitsImage->numPixels;i++){
						for( j=0;j<naxes[0];j++){
							fitsImage->pixels[i].vLambda[j]=vAuxLambdas[j];
							fitsImage->vLambdaImagen[contLambda++] = fitsImage->pixels[i].vLambda[j];
						}
					}								
				}
				else if(naxis == 3){  // matrix of lambdas  
					if( naxes[0]!= fitsImage->rows || naxes[1]!= fitsImage->cols || naxes[2]!=fitsImage->nLambdas){ // image of lambas has different size of spectra image 
						printf("\n IMAGE OF LAMBAS HAS DIFFERENT SIZE OF SPECTRA IMAGE. SIZE SPECTRA %d X %d X %d. SIZE LAMBDA IMAGE %ld X %ld X %ld", fitsImage->rows, fitsImage->cols, fitsImage->nLambdas,naxes[0], naxes[1], naxes[2]);
						freeFitsImage(fitsImage);
						return 0;
					}
					// READ ALL FILE IN ONLY ONE ARRAY 
					// WE ASSUME THAT DATA COMES IN THE FORMAT ROW x COL x LAMBDA
					int numLambdas2Read = naxes[0]*naxes[1]*naxes[2];
					PRECISION  * vAuxLambdas = calloc(numLambdas2Read, sizeof(PRECISION));
					first = 1;
					//fits_read_img(fptr, datatype, first, numLambdas2Read, &nulval, vAuxLambdas, &anynul, &status);
					long fpixel [3] = {1,1,1};
					fits_read_pix(fptr, TDOUBLE, fpixel, naxes[0]*naxes[1]*naxes[2], &nulval, vAuxLambdas, &anynul, &status);
					if(status){
						fits_report_error(stderr, status);
						return 0;	
					}
					int offset_1 = naxes[1] * naxes[0];
					int offset_2 = naxes[0];
					int contLambda = 0;
					for( i=0;i<naxes[2];i++){ // LAMBDA
						for( j=0;j<naxes[1];j++){ // COLS 
							for( k=0;j<naxes[0];k++){ // ROWS
								fitsImage->pixels[ ((j* offset_2) + k) ].vLambda[i] = vAuxLambdas [ (i*offset_1) + (j*offset_2)+ k];
								fitsImage->vLambdaImagen[contLambda++] = fitsImage->pixels[ ((j* offset_2) + k) ].vLambda[i];
							}
						}
					}
				}
			}
			else{
				printf("\n NAXIS FROM LAMBA FILE IS NOT VALID %d ** \n", naxis);
				freeFitsImage(fitsImage);
				return 0;
			}
			// CLOSE FILE FITS LAMBDAS
			fits_close_file(fptr, &status);
			if (status){
				fits_report_error(stderr, status);
				return 0;
			}
		}
		else {
			printf("\n WE CAN NOT OPEN FILE OF LAMBAS ** \n");
			if (status) fits_report_error(stderr, status); /* print any error message */
			freeFitsImage(fitsImage);
			return 0;
		}
	}
	else {
		printf("\n WE CAN NOT READ PARAMETERS FROM THE FILE  ** \n");
		if (status) fits_report_error(stderr, status); /* print any error message */
		freeFitsImage(fitsImage);
		return 0;
	}
	printf("\n LAMBDA IMAGE READ");
	printf("\n**********");

	return 1;

}




void freeFitsImage(FitsImage * image){
	int i;
	for( i=0;i<image->numPixels;i++){
		free(image->pixels[i].spectro);
		free(image->pixels[i].vLambda);
	}
	free(image->pixels);
	free(image->spectroImagen);
	free(image->vLambdaImagen);
	free(image);
}



int writeFitsImageModels(char * fitsFile, int numRows, int numCols, Init_Model * vInitModel, double * vChisqrf){

	fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
   int status, ii, jj;
	int i, j, h; // indexes for loops
   long  fpixel, nelements, exposure;

	int bitpix =  DOUBLE_IMG; 
   long naxis =   3;  /* 2-dimensional image                            */    
   long naxes[3] = { numRows, numCols, NUMBER_PARAM_MODELS };   /* Image of numRows X numCols x 10 parameters of model and chisqrf */

   remove(fitsFile);               /* Delete old file if it already exists */
   status = 0;         /* initialize status before calling fitsio routines */
   if (fits_create_file(&fptr, fitsFile, &status)) /* create new FITS file */
   	printerror( status );           /* call printerror if error occurs */
	
	 /* write the required keywords for the primary array image.     */
    /* Since bitpix = FLOAT_IMG, this will cause cfitsio to create */
    /* a FITS image with BITPIX = -32 (float) .Note that the BSCALE  */
    /* and BZERO keywords will be automatically written by cfitsio  */
    /* in this case.                                                */
	if ( fits_create_img(fptr,  bitpix, naxis, naxes, &status) ){
		printerror( status );
		return 0;
	}

	double * vModel = calloc(numRows * numCols * NUMBER_PARAM_MODELS, sizeof(double));

	int indexModel = 0;
	for( i=0;i<NUMBER_PARAM_MODELS;i++){
		for( j=0;j<numCols;j++){
			for( h=0; h<numRows;h++){
				switch (i)
				{
				case 0:
					vModel[indexModel++] = vInitModel[( j*numCols) + h].B;
					break;
				case 1:
					vModel[indexModel++] = vInitModel[( j*numCols) + h].gm;
					break;
				case 2:
					vModel[indexModel++] = vInitModel[( j*numCols) + h].az;
					break;
				case 3:
					vModel[indexModel++] = vInitModel[( j*numCols) + h].eta0;
					break;
				case 4:
					vModel[indexModel++] = vInitModel[( j*numCols) + h].dopp;
					break;
				case 5:
					vModel[indexModel++] = vInitModel[( j*numCols) + h].aa;
					break;					
				case 6:
					vModel[indexModel++] = vInitModel[( j*numCols) + h].vlos;
					break;					
				case 7:
					vModel[indexModel++] = vInitModel[( j*numCols) + h].alfa;
					break;					
				case 8:
					vModel[indexModel++] = vInitModel[( j*numCols) + h].S0;
					break;					
				case 9:
					vModel[indexModel++] = vInitModel[( j*numCols) + h].S1;
					break;					
				case 10: // READ FROM CHISQR
					vModel[indexModel++] = vChisqrf[( j*numCols) + h];
					break;					
				default:
					break;
				}
			}
		}
	}

   fpixel = 1;                               /* first pixel to write      */
   //nelements = naxes[0] * naxes[1] * naxes[2];          /* number of pixels to write */


   if ( fits_write_img(fptr, TDOUBLE, fpixel, indexModel, vModel, &status) ){
		printerror( status );
		free(vModel);
		return 0;
	}

	// CLEAN MEMORY 
	free(vModel);

	    /* write another optional keyword to the header */
    /* Note that the ADDRESS of the value is passed in the routine */
   exposure = 1500;
	if ( fits_update_key(fptr, TLONG, "EXPOSURE", &exposure,
		"Total Exposure Time", &status) ){
		printerror( status );           
		return 0;
	}
	
	if ( fits_close_file(fptr, &status) ){        
		printerror( status );
		return 0;
	}
	
	return 1;

}


int writeFitsImageProfiles(char * fitsProfileFile, char * fitsFileOrigin, FitsImage * image){

	fitsfile *infptr, *outfptr;   /* FITS file pointers defined in fitsio.h */
	int status = 0, tstatus, ii = 1, iteration = 0, single = 0, hdupos;
	int i, j, k, h; // indexes for loops
	int hdutype, bitpix, bytepix, naxis = 0, nkeys, datatype = 0, anynul;
	long naxes [4] = {1,1,1,1}, fpixel[4] = {1,1,1,1}; /* The maximun number of dimension that we will read is 4*/
	char card[FLEN_CARD];
	char keyname [FLEN_CARD];
	char value [FLEN_CARD];
	remove(fitsProfileFile);               /* Delete old file if it already exists */
	/* Open the input file and create output file */
   fits_open_file(&infptr, fitsFileOrigin, READONLY, &status);
   fits_create_file(&outfptr, fitsProfileFile, &status);
	if (status != 0) {    
		fits_report_error(stderr, status);
		return(status);
	}

	// read a maximun of 4 dimensions 
	fits_get_img_param(infptr, 4, &bitpix, &naxis, naxes, &status);
	// CREATE A NEW IMAGE 
   fits_create_img(outfptr, bitpix, naxis, naxes, &status);
	if (status) {
		fits_report_error(stderr, status);
		return(status);
	}
	/* copy all the user keywords (not the structural keywords) */
	fits_get_hdrspace(infptr, &nkeys, NULL, &status); 

	for (ii = 1; ii <= nkeys; ii++) {
		fits_read_record(infptr, ii, card, &status);
		fits_read_keyn(infptr,ii,keyname, value, NULL, &status);
		fits_update_card(outfptr, keyname,card, &status);
	}

	// CLOSE THE ORIGIN FILE, HE HAVE ALREADY THE INFORMATION OF KEYWORDS. 

	fits_close_file(infptr, &status);
	if (status){
		fits_report_error(stderr, status);
		return 0;
	}


	switch(bitpix) {
	case BYTE_IMG:
		datatype = TBYTE;
		break;
	case SHORT_IMG:
		datatype = TSHORT;
		break;
	case LONG_IMG:
		datatype = TINT;
		break;
	case FLOAT_IMG:
		datatype = TFLOAT;
		break;
	case DOUBLE_IMG:
		datatype = TDOUBLE;
		break;
	}

	// ALLOCATE MEMORY TO WRITE THE IMAGE
	int numElemWrite = naxes[3]*naxes[2]*naxes[1]*naxes[0];

	float * outputImage = calloc(numElemWrite, sizeof(float));
	int currentLambda = 0, currentRow = 0, currentStokeParameter=0, currentCol = 0, currentPixel;
	 
	int pos_lambda = image->pos_lambda;
	int pos_col = image->pos_col;
	int pos_row = image->pos_row;
	int pos_stokes_parameters = image->pos_stokes_parameters;
	
	int numiter = 0;	
	for( i=0; i <naxes[3]; i++){
		for( j=0;j <naxes[2]; j++){
			for( k=0; k<naxes[1]; k++){
				for( h=0; h<naxes[0]; h++){
//	for( fpixel[3] = 1; fpixel[3]<=naxes[3];fpixel[3]++){
//		for( fpixel[2] = 1; fpixel[2]<=naxes[2];fpixel[2]++){
//			for( fpixel[1] = 1; fpixel[1]<=naxes[1];fpixel[1]++){
//				for( fpixel[0] = 1; fpixel[0]<=naxes[0]; fpixel[0]++){	
					// I NEED TO KNOW THE CURRENT POSITION OF EACH ITERATOR 
					switch (pos_lambda)
					{
						case 0:
							//currentLambda = h;
							currentLambda = fpixel[0]-1;
							break;
						case 1:
							//currentLambda = k;
							currentLambda = fpixel[1]-1;
							break;
						case 2:
							//currentLambda = j;
							currentLambda = fpixel[2]-1;
							break;
						case 3:
							//currentLambda = i;
							currentLambda = fpixel[3]-1;
							break;																						
					}
					switch (pos_stokes_parameters)
					{
						case 0:
							//currentStokeParameter = h;
							currentStokeParameter = fpixel[0]-1;
							break;
						case 1:
							//currentStokeParameter = k;
							currentStokeParameter = fpixel[1]-1;
							break;
						case 2:
							//currentStokeParameter = j;
							currentStokeParameter = fpixel[2]-1;
							break;
						case 3:
							//currentStokeParameter = i;
							currentStokeParameter = fpixel[3]-1;
							break;																						
					}
					switch (pos_row)
					{
						case 0:
							//currentRow = h;
							currentRow = fpixel[0]-1;
							break;
						case 1:
							//currentRow = k;
							currentRow = fpixel[1]-1;
							break;
						case 2:
							//currentRow = j;
							currentRow = fpixel[2]-1;
							break;
						case 3:
							//currentRow = i;
							currentRow = fpixel[3]-1;
							break;																						
					}
					switch (pos_col)
					{
						case 0:
							//currentCol = h;
							currentCol = fpixel[0]-1;
							break;
						case 1:
							//currentCol = k;
							currentCol = fpixel[1]-1;;
							break;
						case 2:
							//currentCol = j;
							currentCol = fpixel[2]-1;
							break;
						case 3:
							//currentCol = i;
							currentCol = fpixel[3]-1;
							break;																						
					}
					//double pixel = image->pixels[(currentRow*image->cols) + currentCol].spectro[currentLambda+(image->nLambdas * currentStokeParameter)];
					//fits_write_pix(outfptr, datatype, fpixel, 1, &pixel, &status);			
					outputImage[(i*naxes[2]*naxes[1]*naxes[0]) + (j*naxes[1]*naxes[0]) + (k*naxes[0]) + h] = image->pixels[(currentRow*image->cols) + currentCol].spectro[currentLambda+(image->nLambdas * currentStokeParameter)];
					//outputImage[numiter++] = image->pixels[(currentCol*image->rows) + currentRow].spectro[currentLambda+(image->nLambdas * currentStokeParameter)];
				}
			}
		}
	}

    /* write the array of unsigned integers to the FITS file */
	if ( fits_write_img(outfptr, datatype, 1, numElemWrite, outputImage, &status) ){
		printerror( status );
		free(outputImage);
		return 0;
	}

	free(outputImage);
	fits_close_file(outfptr,  &status);
	if(status){
		printerror( status );
		return 0;
	}
	return 1;
}

int writeFitsImageModelsWithArray(char * fitsFile, int numRows, int numCols, double * eta0, double * B, double * vlos, double * dopp, double * aa, double * gm, double * az, double * S0, double * S1, double * mac, double * alfa, double * vChisqrf){

	fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
   int status, ii, jj;
	int i, j, h; // indexes for loops
   long  fpixel, nelements, exposure;

	int bitpix =  FLOAT_IMG; /* 16-bit unsigned short pixel values       */
   long naxis =   3;  /* 2-dimensional image                            */    
   long naxes[3] = { numRows, numCols, NUMBER_PARAM_MODELS };   /* Image of numRows X numCols x 10 parameters of model and chisqrf */

   remove(fitsFile);               /* Delete old file if it already exists */
   status = 0;         /* initialize status before calling fitsio routines */
   if (fits_create_file(&fptr, fitsFile, &status)) /* create new FITS file */
   	printerror( status );           /* call printerror if error occurs */
	
	 /* write the required keywords for the primary array image.     */
    /* Since bitpix = FLOAT_IMG, this will cause cfitsio to create */
    /* a FITS image with BITPIX = -32 (float) .Note that the BSCALE  */
    /* and BZERO keywords will be automatically written by cfitsio  */
    /* in this case.                                                */
	if ( fits_create_img(fptr,  bitpix, naxis, naxes, &status) ){
		printerror( status );
		return 0;
	}

	double * vModel = malloc(numRows * numCols * NUMBER_PARAM_MODELS);

	int indexModel = 0;
	for( i=0;i<NUMBER_PARAM_MODELS;i++){
		for( j=0;j<numCols;j++){
			for( h=0; h<numRows;h++){
				switch (i)
				{
				case 0:
					vModel[indexModel++] = B[( j*numRows) + h];
					break;
				case 1:
					vModel[indexModel++] = gm[( j*numRows) + h];
					break;
				case 2:
					vModel[indexModel++] = az[( j*numRows) + h];
					break;
				case 3:
					vModel[indexModel++] = eta0[( j*numRows) + h];
					break;
				case 4:
					vModel[indexModel++] = dopp[( j*numRows) + h];
					break;
				case 5:
					vModel[indexModel++] = aa[( j*numRows) + h];
					break;					
				case 6:
					vModel[indexModel++] = vlos[( j*numRows) + h];
					break;					
				case 7:
					vModel[indexModel++] = alfa[( j*numRows) + h];
					break;					
				case 8:
					vModel[indexModel++] = S0[( j*numRows) + h];
					break;					
				case 9:
					vModel[indexModel++] = S1[( j*numRows) + h];
					break;					
				case 10: // READ FROM CHISQR
					vModel[indexModel++] = vChisqrf[( j*numRows) + h];
					break;					
				default:
					break;
				}
			}
		}
	}

   fpixel = 1;                               /* first pixel to write      */
   nelements = naxes[0] * naxes[1] * naxes[2];          /* number of pixels to write */

    /* write the array of unsigned integers to the FITS file */
   if ( fits_write_img(fptr, TUSHORT, fpixel, nelements, vModel, &status) ){
		printerror( status );
		free(vModel);
		return 0;
	}

	// CLEAN MEMORY 
	free(vModel);

	    /* write another optional keyword to the header */
    /* Note that the ADDRESS of the value is passed in the routine */
    exposure = 1500;
	if ( fits_update_key(fptr, TLONG, "EXPOSURE", &exposure,
		"Total Exposure Time", &status) ){
		printerror( status );           
		return 0;
	}

	if ( fits_close_file(fptr, &status) ){                /* close the file */
		printerror( status );
		return 0;
	}
	
	return 1;

}


/*--------------------------------------------------------------------------*/
void printerror( int status)
{
    /*****************************************************/
    /* Print out cfitsio error messages and exit program */
    /*****************************************************/
    if (status)
    {
       fits_report_error(stderr, status); /* print error report */
       //exit( status );    /* terminate the program, returning error status */
    }
    return;
}


int readParametersFileInput(char * fileParameters, int * maxIter, int * clasicalEstimate, int * printSintesis, char * nameInputFileSpectra, char * nameInputFileLambda,char * nameInputFileLines, char * nameInputFileInitModel,  PRECISION *  centralLambda, char *nameOutputFileModels, char * nameOutputFileProfiles, int * useConvolution, char * nameInputFilePSF, PRECISION * FWHM, PRECISION * DELTA, int * NMUESTRAS_G){

	// try open the file with the 
	FILE * fReadParameters;
	char LINE [4096], * returnLine;
	char comment[200], name[100], name2[100];
	fReadParameters = fopen(fileParameters, "r");
	if (fReadParameters == NULL)
	{
		printf("Error opening the file of parameters, it's possible that the file doesn't exist. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;
	}
	int rfscanf; 
	// READ NUM_ITER
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, maxIter,comment);
	
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param MAX_ITER. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if( *maxIter <= 0){
		printf("milos: Error in MAX_ITER parameter. review it. Not accepted: %d\n", *maxIter);
		return 0;
	}


	// READ CLASICAL ESTIMATE
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!", name, clasicalEstimate, comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param CLASICAL_ESTIMATE. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	if (*clasicalEstimate != 0 && *clasicalEstimate != 1 && *clasicalEstimate != 2)
	{
		printf("milos: Error in CLASSICAL_ESTIMATES parameter. [0,1,2] are valid values. Not accepted: %d\n", *clasicalEstimate);
		return 0;
	}


	// READ if print printSintesis
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!",name, printSintesis, comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param PRINT_SINTESIS. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}

	// READ nameInputFileSpectra
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;
	rfscanf = sscanf(LINE,"%99[^:]:%s%99[^!]!",name, nameInputFileSpectra,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param INPUT FILE SPECTRO. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}

	// READ nameInputFileLambda
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;	
	rfscanf = sscanf(LINE,"%99[^:]:%s%99[^!]!", name, nameInputFileLambda,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param INPUT FILE OF LAMBDAS. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}

	// READ nameInputFileLines
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;
	rfscanf = sscanf(LINE,"%99[^:]:%s%99[^!]!",name, nameInputFileLines,comment);

	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param INUT FILE LINES. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}

	// READ nameInputFileInitModel
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;
	rfscanf = sscanf(LINE,"%99[^:]:%s%99[^!]!",name, nameInputFileInitModel,comment);

	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param INUT FILE LINES. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}


	// READ CENTRAL LAMBDA
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;						
	rfscanf = sscanf(LINE,"%99[^:]:%lf%99[^!]!",name, centralLambda,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param DELTA. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}


	// READ nameOutputFileModels
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;		
	rfscanf = sscanf(LINE,"%99[^:]:%s%99[^!]!", name, nameOutputFileModels, comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param OUTPUT FILE OF MODELS. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}

	//READ nameOutputFileProfiles
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;			
	rfscanf = sscanf(LINE,"%99[^:]:%s%99[^!]!", name, nameOutputFileProfiles,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param OUTPUT FILE PROFILES Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	// READ if use Convolution
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;				
	rfscanf = sscanf(LINE,"%99[^:]:%i%99[^!]!", name, useConvolution,comment);
	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param IF USE CONVOLUTION. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}


	// READ nameInputFilePSF
	returnLine = fgets(LINE,4096,fReadParameters);
	if(returnLine == NULL) return 0;
	rfscanf = sscanf(LINE,"%99[^:]:%s%99[^!]!",name, nameInputFilePSF,comment);

	if(rfscanf ==0 || rfscanf == EOF){
		printf("Error reading the file of parameters, param INUT FILE LINES. Please verify it. \n");
		printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
		return 0;		
	}
	

	if(*useConvolution==1){
		// READ FWHM
		returnLine = fgets(LINE,4096,fReadParameters);
		if(returnLine == NULL) return 0;						
		rfscanf = sscanf(LINE,"%99[^:]:%lf%99[^!]!", name, FWHM,comment);
		if(rfscanf ==0 || rfscanf == EOF){
			printf("Error reading the file of parameters, param FWHM. Please verify it. \n");
			printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
			return 0;		
		}
		// READ DELTA
		returnLine = fgets(LINE,4096,fReadParameters);
		if(returnLine == NULL) return 0;						
		rfscanf = sscanf(LINE,"%99[^:]:%lf%99[^!]!",name, DELTA,comment);
		if(rfscanf ==0 || rfscanf == EOF){
			printf("Error reading the file of parameters, param DELTA. Please verify it. \n");
			printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
			return 0;		
		}
		// READ NMUESTRASG
		returnLine = fgets(LINE,4096,fReadParameters);
		if(returnLine == NULL) return 0;						
		rfscanf = sscanf(LINE,"%99[^:]:%i[^!]%s", name, NMUESTRAS_G,comment);
		if(rfscanf ==0 || rfscanf == EOF){
			printf("Error reading the file of parameters, param NMUESTRAS_G. Please verify it. \n");
			printf("\n ******* THIS IS THE NAME OF THE FILE RECEVIED : %s \n", fileParameters);
			return 0;		
		}		
	}
	return 1;

}

/**
 * 
 * 
 * 
 * */
int readFileCuanticLines(char * inputLineFile, PRECISION * cuanticDat, PRECISION centralLambda){
	// try open the file with the 
	FILE * fp;
	char * line = NULL;

   size_t len = 0;
   ssize_t read;
	char atomo [2];
	fp = fopen(inputLineFile, "r");
   if (fp == NULL)
   	return 0;

	int indexLine, ionicState;
	double damping, potentialExcitation, logGf;
	PRECISION lambdaLine;
	int SLOI, SUPI;
	PRECISION LLOI,JLOI,LUPI,JUPI;
	char levelD,levelU;
	while ((read = getline(&line, &len, fp)) != -1) {
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

		}
   }
	cuanticDat[0] =1 ; // LINE NUMBER 1
	cuanticDat[1] = SLOI;
	cuanticDat[2] = LLOI;
	cuanticDat[3] = JLOI;
	cuanticDat[4] = SUPI;
	cuanticDat[5] = LUPI;
	cuanticDat[6] = JUPI;

	printf("\n***********************************");
	printf("\n\n CUANTIC NUMBERS READ FROM FILE: ");
	printf("\n********* SLOI: %fd",cuanticDat[1]);
	printf("\n********* LLOI: %fd",cuanticDat[2]);
	printf("\n********* JLOI: %fd",cuanticDat[3]);
	printf("\n********* SUPI: %fd",cuanticDat[4]);
	printf("\n********* LUPI: %fd",cuanticDat[5]);
	printf("\n********* JUPI: %fd",cuanticDat[6]);
	printf("\n***********************************");


	
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
		if(strstr(name,"B")!=NULL){ // B
			INIT_MODEL->B = aux_value;
		}
		if(strstr(name,"GM")!=NULL){ // GM
			INIT_MODEL->gm = aux_value;
		}
		if(strstr(name,"AZI")!=NULL){ // AZI
			INIT_MODEL->az = aux_value;
		}
		if(strstr(name,"ETHA0")!=NULL){ // ETHA0
			INIT_MODEL->eta0 = aux_value;
		}
		if(strstr(name,"LAMBDADOPP")!=NULL){ // LAMBDADOPP
			INIT_MODEL->dopp = aux_value;
		}
		if(strstr(name,"AA")!=NULL){ // AA
			INIT_MODEL->aa = aux_value;
		}
		if(strstr(name,"ALFA")!=NULL){ // ALFA
			INIT_MODEL->alfa = aux_value;
		}
		if(strstr(name,"MAC")!=NULL){ // MAC
			INIT_MODEL->mac = aux_value;
		}		
		if(strstr(name,"VLOS")!=NULL){ // VLOS
			INIT_MODEL->vlos = aux_value;
		}
		if(strstr(name,"S0")!=NULL){ // S0
			INIT_MODEL->S0 = aux_value;
		}
		if(strstr(name,"S1")!=NULL){ // S1
			INIT_MODEL->S1 = aux_value;
		}				
	}
	fclose(fReadInitModel);


	return 1;
}


/**
 * 
 * */
int readPSFFile(PRECISION * deltaLambda, PRECISION * PSF, char * nameInputPSF){

	// first of all read the number of lines to give size for arrays deltaLambda and PSF
	FILE *fp;

	// alloc memory 

	char * line = NULL;
	int rfscanf;
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
		rfscanf = sscanf(line,"%lf  %lf", &delta, &psf);
		deltaLambda[index] = delta;
		PSF[index] = psf;
		index++;
	}

	fclose(fp);
	return 1;
}