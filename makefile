CC=mpicc

CFLAGS=
CFLAGS+=-O3 
CFLAGS+=-g 
#CFLAGS+=-floop-block
#CFLAGS+=-ffinite-math-only
#CFLAGS+=-floop-parallelize-all
#CFLAGS+=-floop-strip-mine
#CFLAGS+=-fno-rounding-math -fno-trapping-math -fno-signaling-nans -ffinite-math-only -fno-signed-zeros -funsafe-math-optimizations
CFLAGS+=-ffast-math
CFLAGS+=-fno-omit-frame-pointer
#CFLAGS+=-DGLOBAL_DEBUG
#CFLAGS+=-Wall -Wextra
#CFLAGS+=-Wconversion
#CFLAGS+=-Wno-unused-but-set-variable -Wno-unused-parameter

#LDFLAGS=-L/home/bsc15/bsc15557/root/usr/lib/
LDLIBS= -fopenmp -lm -lcfitsio -lnsl -lgsl -lgslcblas  -lfftw3_omp -lfftw3 -ldl -lpthread  -pg  # -shared-intel
BIN= milos milosMPI milosMPI_DIR 


all: $(BIN)

milos: calculosCompartidos.o fftw.o fgauss.o fvoigt.o  milos.o me_der.o mil_sinrf.o lib.o create_cuantic.o utilsFits.o milosUtils.o convolution.o readConfig.o slog.o

milosMPI:  calculosCompartidos.o fgauss.o fvoigt.o  milosMPI.o me_der.o mil_sinrf.o lib.o create_cuantic.o utilsFits.o milosUtils.o convolution.o readConfig.o slog.o

milosMPI_DIR:  calculosCompartidos.o fgauss.o fvoigt.o  milosMPI_DIR.o me_der.o mil_sinrf.o lib.o create_cuantic.o utilsFits.o milosUtils.o convolution.o readConfig.o slog.o

clean:
	rm -f *.o $(BIN)
