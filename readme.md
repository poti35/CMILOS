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
Number of cycles          (*):15          ! 
Observed profiles         (*):/home/mcabrera/SolarImages/fits_CRISP_spot/2014.09.28_09_18_00_t098.fits !
Stray light file             :    ! (none=no stray light contam)
PSF file                     :    ! (none=no convolution with PSF)
wavelength file           (*):/home/mcabrera/MILOS/utilsFits/lambda_30.fits ! (none=automatic selection)
Atomic parameters file       :/home/mcabrera/MILOS/initFiles/LINES      ! (none=DEFAULT LINES file)
Initial guess model 1     (*):/home/mcabrera/MILOS/initFiles/initModel.mod  !
Initial guess model 2        :                                           !
Weight for Stokes I          :1               ! (DEFAULT=1; 0=not inverted)
Weight for Stokes Q          :1                 ! (DEFAULT=1; 0=not inverted)
Weight for Stokes U          :1                 ! (DEFAULT=1; 0=not inverted)
Weight for Stokes V          :1                 ! (DEFAULT=1; 0=not inverted)
mu=cos (theta)               :1                 ! (DEFAULT: mu=1. mu<0 => West)
Estimated S/N for I          :1000              ! (DEFAULT: 1000) 
Continuum contrast           :                  ! (DEFAULT: not used)
Tolerance for SVD            :1e-3                ! (DEFAULT value: 1e-3)
Initial diagonal element     :0.1                 ! (DEFAULT value: 1.e-1)
Splines/Linear Interpolation :0                ! (0 or blank=splines, 1=linear)
Use Gaussian PSF             :1               ! (0=not apply gaussian PSF filter, 1= apply gaussian PSF filter)
FWHM                         :0.075          ! Full Width at Half Maximum for create a Gaussian PSF Filter 
ntl                          :1              ! (number of spectral lines)
nliobs                       :30             ! (number of observed wavelengths)
Central Wave Lenght          :6173.3356      ! Lambda central to use from file LINES 
Line To Continiuum Absorption :1     ! (0 or blank=do not generate models for different ETA0)
Magnetic Field Strength       :1     ! (0 or blank=do not generate models for different B)
Line Of Sight Velocity        :1     ! (0 or blank=do not generate models for different VLOS)
Doopler Width                 :1     ! (0 or blank=do not generate models for different DOPP)
Damping Parameter             :1     ! (0 or blank=do not generate models for different AA)
Magnetic FieldInclination     :1     ! (0 or blank=do not generate models for different GM)
Magnetic Field Azimuth        :1     ! (0 or blank=do not generate models for different AZ)
Source Function Constant      :1     ! (0 or blank=do not generate models for different S0)
Source Function Gradient      :1     ! (0 or blank=do not generate models for different S1)
Macroturbulent Velocity       :1     ! (0 or blank=do not generate models for different MAC)
Filling Factor                :0     ! (0 or blank=do not generate models for different ALPHA)
Save Chisqr                   :1     ! (0 = don't store the resultant chisqr of inversion in model output file,  1= Store chisqr in the last dimension of model output file.)
Use Classical Estimates       :0     !(0 or blank= don't use classical estimates, 1= Use classical estimates)
Use RTE Inversion             :1     !( or blank= don't use rte inversion, TRUE= Use rte inversion)
Save Synthesis Profile        :0     !(0 or blank= don't save synthesis profile from result model of inversion, 1= Create synthesis from result model of inversion )
Output Model File             :/home/mcabrera/MILOS/MODELOS/2014.09.28_09_18_00_t098_inversion.fits       ! model result file
Output Synthesis File         :/home/mcabrera/MILOS/MODELOS/sintesis_secuencial_without_convolution.fits ! synthesis result file
noise                         :1e-3          ! Scalar with the noise values for the Stokes profiles (Ignored when SIGMA is set)
toplim                        :1e-12         ! Optional minimum relative difference between two succesive merit-function values. Default: 1e-12
```





