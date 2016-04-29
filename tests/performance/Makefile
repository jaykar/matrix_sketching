CXX = g++
CXXFLAGS = -std=c++14 -Wall -g -v

BOOST_LIB_PATH = -I /usr/local/include/

BOOST_COMPILE_FLAG = -O3 -D NDEBUG -DBOOST_UBLAS_NDEBUG
INTEL_COMPILE_FLAG = -m64 -I$(MKLROOT)/include

LAPACK_LIB_PATH = -L /usr/local/lib/ -llapack -L /usr/local/lib/ -lblas
INTEL_LIB_PATH =  -L$(MKLROOT)/lib -Wl,-rpath,$(MKLROOT)/lib -L$(MKLROOT)/../compiler//lib -Wl,-rpath,$(MKLROOT)/../compiler/lib -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -lmkl_scalapack_lp64 -lmkl_blacs_mpich_lp64
ARMA_LIB_PATH = -larmadillo -framework Accelerate

INTERFACE = -I ../../src/interface

time: time.o
	$(CXX) $(CXXFLAGS) $(LAPACK_LIB_PATH) $(INTEL_LIB_PATH) $(ARMA_LIB_PATH) $(BOOST_COMPILE_FLAG) $(INTEL_COMPILE_FLAG) -o time time.o

time.o: time.cpp
	$(CXX) $(CXXFLAGS) $(BOOST_LIB_PATH) $(INTERFACE) -c time.cpp

.PHONY: clean
clean:
	rm -f *.o a.out core main time

.PHONY: all
all: clean main