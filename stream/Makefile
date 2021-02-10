export ARGO_INSTALL_DIRECTORY=${HOME}/argodsm

CC      = mpic++
CFLAGS  =-O2 -fopenmp -std=c++11 -DSTREAM_ARRAY_SIZE=10000 -DRANDACC=0 -DTEAMINIT=1
CFLAGS +=-L${ARGO_INSTALL_DIRECTORY}/build/lib
CFLAGS +=-I${ARGO_INSTALL_DIRECTORY}/src
CFLAGS +=-largo -largobackend-mpi -lnuma -lrt
CFLAGS +=-Wl,-rpath=${ARGO_INSTALL_DIRECTORY}/build/lib

all: stream

stream: stream.cpp
	$(CC) $(CFLAGS) stream.cpp -o stream

clean:
	rm -f stream *.o

