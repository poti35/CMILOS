CC=mpicc

CFLAGS=
CFLAGS+=-O3 
CFLAGS+=-g 
#CFLAGS+=-D USE_DOUBLE_PRECISION=double
CFLAGS+=-ffast-math
CFLAGS+=-fno-omit-frame-pointer
#CFLAGS+=-Wall -Wextra
#CFLAGS+=-Wconversion
#CFLAGS+=-Wno-unused-but-set-variable -Wno-unused-parameter

LDLIBS= -lm -lcfitsio -lnsl -lgsl -lgslcblas -lfftw3 -ldl -lpthread # -shared-intel
BIN= milos milosMPI 


all: $(BIN)

milos: calculosCompartidos.o fftw.o fgauss.o fvoigt.o  milos.o me_der.o mil_sinrf.o lib.o create_cuantic.o utilsFits.o milosUtils.o convolution.o readConfig.o

milosMPI:  calculosCompartidos.o fgauss.o fvoigt.o  milosMPI.o me_der.o mil_sinrf.o lib.o create_cuantic.o utilsFits.o milosUtils.o convolution.o readConfig.o

clean:
	rm -f *.o $(BIN)
