
//##############################################
//SVD CONFIGURATION
#define USE_SVDCMP 0    //1 for using SVDCMP and 0 for using SVD_CORDIC  -->  Note: the SVDCMP doesn't work in float! only double

#define NORMALIZATION_SVD 1 //1 for using normalization matrixes ONLY  in the SVD_CORDIC

#define NUM_ITER_SVD_CORDIC 9 //9,18,27,36  --> 18 parece ok!

#define LIMITE_INFERIOR_PRECISION_SVD pow(2.0,-54)
#define LIMITE_INFERIOR_PRECISION_TRIG pow(2.0,-23)
#define LIMITE_INFERIOR_PRECISION_SINCOS pow(2.0,-23)
//#############################################