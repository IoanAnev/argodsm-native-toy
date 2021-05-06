# Set your ArgoDSM installation directory
export ARGO_INSTALL_DIRECTORY=${HOME}/argodsm/install

CXX      =mpic++
CXXFLAGS =-O2 -fopenmp -std=c++11
CPPFLAGS =-DSTREAM_ARRAY_SIZE=10000 -DRANDACC=0 -DTEAMINIT=1

IDFLAGS +=-I${ARGO_INSTALL_DIRECTORY}/include/argo
LDFLAGS +=-L${ARGO_INSTALL_DIRECTORY}/lib
LDFLAGS +=-Wl,-rpath=${ARGO_INSTALL_DIRECTORY}/lib
LDLIBS  +=-largo -largobackend-mpi -lnuma -lrt

all: stream

stream: stream.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o stream $(IDFLAGS) $(LDFLAGS) stream.cpp $(LDLIBS)

clean:
	rm -f stream *.o

