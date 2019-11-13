#include "defines.h"

/**
 * 
 */
int covarm(PRECISION *w, PRECISION *sig, PRECISION *spectro, int nspectro, PRECISION *spectra, PRECISION *d_spectra,
			  PRECISION *beta, PRECISION *alpha);
/**
 * 
 */
double * leeVector(char *nombre,int tam);
/**
 * 
 */
int multmatrixIDL(double *a,int naf,int nac, double *b,int nbf,int nbc,double **resultOut,int *fil,int *col);
/**
 * 
 */
int multmatrix2(double *a,int naf,int nac, PRECISION *b,int nbf,int nbc,double **result,int *fil,int *col);
/**
 * 
 */
double * transpose(double *mat,int fil,int col);
/**
 * 
 */
double *totalParcial(double * A, int f,int c,int dire);
/**
 * 
 */
double *totalParcialMatrix(double * A, int f,int c,int p);
/**
 * 
 */
double total(double * A, int f,int c);
/**
 * 
 */
int multmatrix(PRECISION *a,int naf,int nac, PRECISION *b,int nbf,int nbc,PRECISION *result,int *fil,int *col);
/**
 * 
 * */
int multmatrixCblas(PRECISION *a,int naf,int nac, PRECISION *b,int nbf,int nbc,PRECISION *result,int *fil,int *col);
/**
 * 
 */
double fchisqr(PRECISION * spectra,int nspectro,PRECISION *spectro,PRECISION *w,PRECISION *sig,double nfree);
/**
 * 
 */
PRECISION * transposef(PRECISION *mat,int fil,int col);
/**
 * 
 */
int multmatrixIDLf(PRECISION *a,int naf,int nac,PRECISION *b,int nbf,int nbc,PRECISION **resultOut,int *fil,int *col);
/**
 * 
 */
int multmatrixIDLValue(PRECISION *a,int naf,int nac,PRECISION *b,int nbf,int nbc,PRECISION *result,int *fil,int *col,PRECISION value);
/**
 * 
 */
void totalParcialMatrixf(PRECISION * A, int f,int c,int p,PRECISION *result);
/**
 * 
 */
void totalParcialf(PRECISION * A, int f,int c,PRECISION * result);
/**
 * 
 */
int multmatrix_transpose(PRECISION *a,int naf,int nac, PRECISION *b,int nbf,int nbc,PRECISION *result,int *fil,int *col,PRECISION value);

int multmatrix_transpose_omp(PRECISION *a,int naf,int nac, PRECISION *b,int nbf,int nbc,PRECISION *result,int *fil,int *col,PRECISION value);

int multmatrix_transpose_cblas(PRECISION *a,int naf,int nac, PRECISION *b,int nbf,int nbc,PRECISION *result,int *fil,int *col,PRECISION value);
/**
 * 
 */
int multmatrix_transposeD(double *a,int naf,int nac, double *b,int nbf,int nbc,double *result,int *fil,int *col);
/**
 * 
 */
int multmatrix3(PRECISION *a,int naf,int nac,double *b,int nbf,int nbc,double **result,int *fil,int *col);
/**
 * 
 */
int CalculaNfree(int nspectro);


/**
 * 
 * */
PRECISION mean(double *dat, int numl);

/**
 * PRINT 
 * */
void printProgress (double percentage);

/*
* Check if path is a directory or not. 
*/
int isDirectory(const char *path);