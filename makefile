CC=mpicc

CFLAGS=
CFLAGS+=-Ofast -pg -g 
CFLAGS+=-fno-omit-frame-pointer 
#CFLAGS+=-DGLOBAL_DEBUG
#CFLAGS+=-Wall -Wextra
#CFLAGS+=-Wconversion
#CFLAGS+=-Wno-unused-but-set-variable -Wno-unused-parameter

#LDFLAGS=-L/home/bsc15/bsc15557/root/usr/lib/
LDLIBS= -pg -lm -lcfitsio -lnsl -lgsl -lgslcblas -lfftw3 -ldl -lconfig #-shared-intel
BIN= milos milosMPI 

all: $(BIN)

milos: calculosCompartidos.o fftw.o fgauss.o fvoigt.o  milos.o me_der.o mil_sinrf.o lib.o create_cuantic.o utilsFits.o milosUtils.o convolution.o readConfig.o

milosMPI:  calculosCompartidos.o fgauss.o fvoigt.o  milosMPI.o me_der.o mil_sinrf.o lib.o create_cuantic.o utilsFits.o milosUtils.o convolution.o readConfig.o


clean:
	rm -f *.o $(BIN)
