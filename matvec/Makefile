# Set your ArgoDSM installation directory
export ARGO_INSTALL_DIRECTORY=${HOME}/argodsm/install

CXX      =mpic++
CXXFLAGS =-O3 -fopenmp -std=c++11

IDFLAGS +=-I${ARGO_INSTALL_DIRECTORY}/include/argo
LDFLAGS +=-L${ARGO_INSTALL_DIRECTORY}/lib
LDFLAGS +=-Wl,-rpath=${ARGO_INSTALL_DIRECTORY}/lib
LDLIBS  +=-largo -largobackend-mpi -lnuma -lrt

all: matvec

matvec: matvec.cpp
	$(CXX) $(CXXFLAGS) -o matvec $(IDFLAGS) $(LDFLAGS) matvec.cpp $(LDLIBS)

clean:
	rm -f matvec *.o

