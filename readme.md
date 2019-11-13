# C-MILOS


## Description 

Intructions to compile and execute  C-MILOS. 

Look **Deployment** to understand how to execute the program.


### Pre-conditions 

The libraries CFITSIO AND OpenMPI must be installed in the version 3.3.4 and 1.4-4, respectively. 


### Instalation

_After install the pre-conditions libraries, execute: make , to compile the programs milos and milosMPI_



## Deployment

To execute Milos in sequential way, you must to the following call: 
_milos parameters.trol_

The file  parameters.trol  will have the parameters of entry with the following order 

```
NumberOfCycles                  : 15        // NUMBER_OF_CYCLES 
ObservedProfiles                : "/home/mcabrera/SolarImages/fits_CRISP_spot/2014.09.28_09_18_00_t098.fits"         // OBSERVED_PROFILES
StrayLightFile                  : ""                                                                                //STRAY LIGHT FILE 
PSFFile                         : ""                                                                  //PSF_FILE 
WavelengthFile                  : "/home/mcabrera/MILOS/utilsFits/lambda_30.fits"                   //WAVE_LENGHT_FILE 
AtomicParametersFile            : "/home/mcabrera/MILOS/initFiles/LINES_6173"                       //ATOMIC_PARAMETERS_FILE 
InitialGuessModel               : "/home/mcabrera/MILOS/initFiles/initModel.mod"                    //INITIAL_GUESS_MODEL
WeightForStokesI                : 1.              //  (DEFAULT=1; 0=not inverted)
WeightForStokesQ                : 1.              //  (DEFAULT=1; 0=not inverted)
WeightForStokesU                : 1.              //  (DEFAULT=1; 0=not inverted)
WeightForStokesV                : 1.              //  (DEFAULT=1; 0=not inverted)
InvertMacroturbulence           : FALSE          //  (FALSE=no, TRUE=yes)
InvertFillingFactor             : FALSE          //  (FALSE=no, TRUE=yes)
InvertStrayLightFactor          : FALSE          //  (FALSE=no, TRUE=yes)
mu                              : 1.              // cos (theta), (DEFAULT: mu=1. mu<0 => West)
EstimatedSNForI                 : 3000           // (DEFAULT: 1000)
ContinuumContrast               : -1               // (DEFAULT: not used)
ToleranceForSVD                 : 1e-3               // (DEFAULT value: 1e-3)
InitialDiagonalElement          : 0.1               // (DEFAULT value: 1.e-1)
UseInterpolarSplinesOrLinear    : 0               // (0=splines, 1=linear)
ConvolveWithPSF                 : true               // (FALSE=no convolve with psf, TRUE=use psf for convolution)
FWHM                            : 0.075               // Full Width at Half Maximum for create a Gaussian PSF Filter
TypeConvolution                 : "FFT"                // possible values: DIRECT or FFT 
GasPressureAtSurface1           : 0.              // (0 =Pe boundary cond.
GasPressureAtSurface2           : 0.              // (0 =Pe boundary cond. 
MagneticPressureTerm            : 0.              // (0 =no, 1=yes
ntl                             : 1              // (number of spectral lines)
nliobs                          : 30               // (number of observed wavelengths)
CentralWaveLenght               : 6173.3356      // Lambda central to use from file LINES 
ETA0_LineToContiniuumAbsorption : 1     // (0 or blank=do not generate models for different ETA0)
B_MagneticFieldStrength         : 1     // (0 or blank=do not generate models for different B)
VLOS_LineOfSightVelocity        : 1     // (0 or blank=do not generate models for different VLOS)
DOPP_DooplerWidth               : 1     // (0 or blank=do not generate models for different DOPP)
AA_DampingParameter             : 1     // (0 or blank=do not generate models for different AA)
GM_MagneticFieldInclination     : 1     // (0 or blank=do not generate models for different GM)
AZ_MagneticFieldAzimuth         : 1     // (0 or blank=do not generate models for different AZ)
S0_SourceFunctionConstant       : 1     // (0 or blank=do not generate models for different S0)
S1_SourceFunctionGradient       : 1     // (0 or blank=do not generate models for different S1)
MAC_MacroturbulentVelocity      : 1     // (0 or blank=do not generate models for different MAC)
ALPHA_FillingFactor             : 0     // (0 or blank=do not generate models for different ALPHA)
SaveChisqr                      : TRUE      // (FALSE = don't store the resultant chisqr of inversion in model output file,  TRUE= Store chisqr in the last dimension of model output file.)
UseClassicalEstimates           : TRUE     // (FALSE or blank= don't use classical estimates, TRUE= Use classical estimates)
UseRTEInversion                 : TRUE      // (FALSE or blank= don't use rte inversion, TRUE= Use rte inversion)
SaveSynthesisProfile            : FALSE     // (FALSE or blank= don't save synthesis profile from result model of inversion, TRUE= Create synthesis from result model of inversion )
OutputModelFile                 : "/home/mcabrera/MILOS/MODELOS/INSERSION_IMAGE_30_SI_PSF_NO_MAC.fits"        // model result file
OutputSynthesisFile             : "/home/mcabrera/MILOS/MODELOS/sintesis_secuencial_without_convolution.fits"   // synthesis result file
sigma                           : [0.001,0.001,0.001,0.001]             // Array with noise values for the profiles. Dimensions: (4). It can't be blank for not apply it. 
noise                           : 1e-3           // Scalar with the noise values for the Stokes profiles (Ignored when SIGMA is set)
toplim                          : 1e-12          // Optional minimum relative difference between two succesive merit-function values. Default: 1e-12
```





