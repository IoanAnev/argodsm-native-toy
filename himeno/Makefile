export ARGO_INSTALL_DIRECTORY=${HOME}/argodsm

CC      = mpic++
CFLAGS  =-O3 -fopenmp -std=c++11
CFLAGS +=-L${ARGO_INSTALL_DIRECTORY}/build/lib
CFLAGS +=-I${ARGO_INSTALL_DIRECTORY}/src
CFLAGS +=-largo -largobackend-mpi -lnuma -lrt
CFLAGS +=-Wl,-rpath=${ARGO_INSTALL_DIRECTORY}/build/lib

all: himenobmtxpa

himenobmtxpa: himenobmtxpa.cpp
	$(CC) $(CFLAGS) himenobmtxpa.cpp -o himenobmtxpa

clean:
	rm -f himenobmtxpa *.o

