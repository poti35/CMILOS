#by manu 2019 (IAA-CSIC)
mil: print svdcmp.o calculosCompartidos.o create_cuantic.o fvoigt.o fgauss.o me_der.o mil_sinrf.o lib.o mpiTools.o utilsFits.o milosUtils.o copyImage milos milosMPI fin

print:
	#clear
	#make clear
	@echo Compilando CMILOS...

create_cuantic.o: create_cuantic.c 
	mpicc -g -c -O1   create_cuantic.c

fvoigt.o: fvoigt.c
	mpicc -g -O1 -c     fvoigt.c

fgauss.o: fgauss.c
	mpicc -g -O1 -c     fgauss.c

me_der.o: me_der.c 
	mpicc -g -O1 -c    me_der.c

mil_sinrf.o: mil_sinrf.c 
	mpicc -g -O1 -c     mil_sinrf.c 

milosUtils.o: milosUtils.c svdcordic.c
	mpicc -g -O1 -c  milosUtils.c

calculosCompartidos.o: calculosCompartidos.c
	mpicc -g -c -O1   calculosCompartidos.c

svdcmp.o: svdcmp.c
	mpicc -g -O1 -c svdcmp.c

lib.o: lib.c
	mpicc -g -O1 -c   lib.c

mpiTools.o: mpiTools.c
	mpicc -g -O1 -c  mpiTools.c

utilsFits.o: utilsFits.c
	mpicc -g -O1 -c  utilsFits.c

milos.o: milos.c 
	gcc -g -O1 -c    milos.c

milosMPI.o: milosMPI.c 
	mpicc -g -O1 -c    milosMPI.c

copyImage: utilsFits.o
	gcc -g -o copyImage copyImage.c utilsFits.o  -L. -lcfitsio -lm -lnsl

milos:  calculosCompartidos.o fgauss.o fvoigt.o  milos.o me_der.o mil_sinrf.o lib.o create_cuantic.o utilsFits.o milosUtils.o svdcmp.o
	gcc -o milos milos.o calculosCompartidos.o create_cuantic.o fgauss.o fvoigt.o me_der.o mil_sinrf.o utilsFits.o milosUtils.o svdcmp.o lib.o -L. -lm -lcfitsio -lnsl

milosMPI:  calculosCompartidos.o fgauss.o fvoigt.o  milosMPI.o me_der.o mil_sinrf.o lib.o create_cuantic.o utilsFits.o mpiTools.o milosUtils.o svdcmp.o
	mpicc -g -o  milosMPI milosMPI.o calculosCompartidos.o create_cuantic.o fgauss.o fvoigt.o me_der.o mil_sinrf.o  lib.o utilsFits.o mpiTools.o milosUtils.o  svdcmp.o -L. -lm -lcfitsio -lnsl

fin:
	@echo --
	@echo All done.

clear:
	rm *.o
	rm milos
	rm milosMPI
