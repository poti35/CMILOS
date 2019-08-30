#!/bin/sh
#PBS –N milos_openmpi
#PBS –e ~/outputMilos/errfile
#PBS –o ~/outputMilos/outputfile

#display date
date

NSLOTS=2
mpirun -np $NSLOTS ~/CMilos_OPENMPI/milosMPI ~/CMilos_OPENMPI/parameters.txt 
