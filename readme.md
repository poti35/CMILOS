# C-MILOS


## Comenzando 

Instrucción para compilar y ejecutar C-MILOS. 

Mira **Deployment** para conocer como ejecutar el programa.


### Pre-requisitos 

Tener instalado la libreria CFITSIO en su versión 3.3.4  y OpenMPI 1.4-4


### Instalación 

_Una vez instaladas las librerias especificadas como prerequisitos ejecutar: make , para compilar tanto el programa milos como milosMPI_



## Deployment

Para ejecutar Milos en forma secuencial realizar la siguiente llamada: 
_milos parameters.txt_

El fichero parameters.txt contendrá los parametros de entrada al programa con el siguiente orden 

```
15  // Número de iteraciones máxismas del algoritmo de inversion 
1   // Indica si usar estimaciones clásicas: 1 si e inversion, 0 no, 2 solo estimaciones clásicas. 
0   // 0 no se escribe fichero de sintesis, 1 si se escribe. 
/home/mcabrera/MILOS/2014.09.28_09_18_00_t012.fits  // Fichero origen con el espectro 
/home/mcabrera/MILOS/lambda.fits                    // Fichero origen con las longitudes de onda
/home/mcabrera/MILOS/MODELOS/result.fits            // Fichero donde almacenar los modelos de resultados. 
/home/mcabrera/MILOS/MODELOS/sintesis.fits          // Fichero de sintesis en caso de ser seleccionada su escritua. 
0   // Indica si usar convolucion o no, 0 no, 1 si . 
NULL  // FWHM
NULL  // DELTA 
NULL  // NMUESTRASG

```



