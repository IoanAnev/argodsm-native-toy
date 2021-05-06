# Set your ArgoDSM installation directory
export ARGO_INSTALL_DIRECTORY=${HOME}/argodsm/install

CXX      =mpic++
CXXFLAGS =-O3 -fopenmp -std=c++17

IDFLAGS +=-I${ARGO_INSTALL_DIRECTORY}/include/argo
LDFLAGS +=-L${ARGO_INSTALL_DIRECTORY}/lib
LDFLAGS +=-Wl,-rpath=${ARGO_INSTALL_DIRECTORY}/lib
LDLIBS  +=-largo -largobackend-mpi -lnuma -lrt

all: matmul

matmul: matmul.cpp
	$(CXX) $(CXXFLAGS) -o matmul $(IDFLAGS) $(LDFLAGS) matmul.cpp $(LDLIBS)

clean:
	rm -f matmul *.o

