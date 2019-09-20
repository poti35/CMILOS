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
Number of Iterations                    :15                             ! by default is 15 
Use classic estimates and RTE           :1                              ! use classical estimates only
Save synthesis profiles                 :0                              ! don't save synthesis profiles
Path of spectro file                    :2014.05.16_18_08_37_t000.fits  ! file with lambda vector
Path of lambda file                     :lambda_12_bueno.fits           ! file with lambda vector
LINE file                               :LINES_6173                     ! file with atomic lines
File with Init Model                    :initModel.mod                  ! File with init model for inversion
Central lambda                          :6173.3356                      ! Lambda central to use from file LINES 
Path of output model file               :result_ce_secuencial_2.fits    ! model result file
Path of output synthesis file           :sintesis_secuencial.fits       ! synthesis result file
Indicate if use convolution             :1                              ! value 0 don't use it , value 1 then use it 
File with PSF                           :ps.psf                         ! Psf file, if not specify we will use Gauss samples
FWHM                                    :0.076                          ! In Armstrong
DELTA                                   :0.069                          ! In Armstrong
Number of samples to create gaussian    :13                             ! Must be odd
```

The first row explains what is the value of that row, after ":" appears the value and after "!" you can put any comment. 



