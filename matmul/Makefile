export ARGO_INSTALL_DIRECTORY=${HOME}/argodsm

CC      = mpic++
CFLAGS  =-O3 -fopenmp -std=c++17
CFLAGS +=-L${ARGO_INSTALL_DIRECTORY}/build/lib
CFLAGS +=-I${ARGO_INSTALL_DIRECTORY}/src
CFLAGS +=-largo -largobackend-mpi -lnuma -lrt
CFLAGS +=-Wl,-rpath=${ARGO_INSTALL_DIRECTORY}/build/lib

all: matmul

matmul: matmul.cpp
	$(CC) $(CFLAGS) matmul.cpp -o matmul

clean:
	rm -f matmul *.o

