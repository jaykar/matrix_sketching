CXX = clang++
CXXFLAGS = -std=c++14 -Wall -g

BOOST_LIB_PATH = -I /usr/local/include/

BOOST_COMPILE_FLAG = -O3 -D NDEBUG -DBOOST_UBLAS_NDEBUG
INTEL_COMPILE_FLAG = -m64 -I$(MKLROOT)/include

# LAPACK_LIB_PATH = -L /usr/local/lib/ -llapack -L /usr/local/lib/ -lblas
INTEL_LIB_PATH =  -L$(MKLROOT)/lib -Wl,-rpath,$(MKLROOT)/lib -L$(MKLROOT)/../compiler//lib -Wl,-rpath,$(MKLROOT)/../compiler/lib -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -lmkl_scalapack_lp64 -lmkl_blacs_mpich_lp64
ARMA_LIB_PATH = -larmadillo -framework Accelerate

INTERFACE = -I ../../sketchy/

time: time.o
	$(CXX) $(CXXFLAGS) $(BOOST_COMPILE_FLAG) $(INTEL_LIB_PATH) $(ARMA_LIB_PATH) -o time time.o

time.o: time.cpp
	$(CXX) $(CXXFLAGS) $(BOOST_LIB_PATH) $(INTERFACE) -c time.cpp

.PHONY: clean
clean:
	rm -f *.o a.out core main time

.PHONY: all
all: clean main