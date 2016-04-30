CXX = g++
CXXFLAGS = -std=c++14 -Wall -g

BOOST_LIB_PATH = -I /usr/local/include/
BOOST_COMPILE_OPTION = -O3 -DNDEBUG -DBOOST_UBLAS_NDEBUG

LAPACK_LIB_PATH = -L /usr/local/lib/ -llapack -L /usr/local/lib/ -lblas

MKL_FLAGS =  -L$(MKLROOT)/lib -Wl,-rpath,$(MKLROOT)/lib -L$(MKLROOT)/../compiler//lib -Wl,-rpath,$(MKLROOT)/../compiler/lib -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -lmkl_scalapack_lp64 -lmkl_blacs_mpich_lp64
MKL_INCLUDE = -I$(MKLROOT)/include

INTERFACE = -I ../../sketchy/

ARMADILLO_LIB = -L../../src/libraries/armadillo -larmadillo

time: time.o
	$(CXX) $(CXXFLAGS) $(LAPACK_LIB_PATH) $(BOOST_COMPILE_OPTION) -o time time.o #$(MKL_FLAGS) $(ARMADILLO_LIB)

time.o: time.cpp
	$(CXX) $(CXXFLAGS) $(BOOST_LIB_PATH) $(INTERFACE) -c time.cpp

.PHONY: clean
clean:
	rm -f *.o a.out core main time *.csv

.PHONY: all
all: clean main
