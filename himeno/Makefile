# Set your ArgoDSM installation directory
export ARGO_INSTALL_DIRECTORY=${HOME}/argodsm/install

CXX      =mpic++
CXXFLAGS =-O3 -fopenmp -std=c++11

IDFLAGS +=-I${ARGO_INSTALL_DIRECTORY}/include/argo
IDFLAGS +=-I${ARGO_INSTALL_DIRECTORY}/include/argo/synchronization
LDFLAGS +=-L${ARGO_INSTALL_DIRECTORY}/lib
LDFLAGS +=-Wl,-rpath=${ARGO_INSTALL_DIRECTORY}/lib
LDLIBS  +=-largo -largobackend-mpi -lnuma -lrt

all: himenobmtxpa

himenobmtxpa: himenobmtxpa.cpp
	$(CXX) $(CXXFLAGS) -o himenobmtxpa $(IDFLAGS) $(LDFLAGS) himenobmtxpa.cpp $(LDLIBS)

clean:
	rm -f himenobmtxpa *.o

