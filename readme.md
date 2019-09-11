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
_milos parameters.txt_

The file  parameters.txt  will have the parameters of entry with the following order 

```
15  // Number of max iterations of the inversion algorithm 
1   // Option 0: not use classic estimations and use RTE inversion, Option 1: use classic estimations and RTE inversion, Option 2: use only RTE inversion 
0   // Option 0 to not write syntesis file and 1 to write it 
/home/mcabrera/MILOS/2014.09.28_09_18_00_t012.fits  // source file with spectro 
/home/mcabrera/MILOS/lambda.fits                    // source file with lambda 
/home/mcabrera/MILOS/MODELOS/result.fits            // file to store the results models 
/home/mcabrera/MILOS/MODELOS/sintesis.fits          // file to store synthesis file 
0   // Option 0 not use convolution, Option 1 use convolution. 
NULL  // FWHM
NULL  // DELTA 
NULL  // NMUESTRASG

```



