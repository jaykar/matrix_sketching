CXX = g++
CXXFLAGS = -std=c++14 -Wall -g -I../interface/ -I../libraries/armadillo/include
LDFLAGS =  -L../libraries/armadillo -larmadillo -L$(MKLROOT)/lib -Wl,-rpath,$(MKLROOT)/lib -L$(MKLROOT)/../compiler//lib -Wl,-rpath,$(MKLROOT)/../compiler/lib -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -lmkl_scalapack_lp64 -lmkl_blacs_mpich_lp64 

arm_interface.o: Arm_interface.hpp 
	$(CXX) $(CXXFLAGS) -c Arm_interface.hpp ../interface/SKMatrix.hpp

.PHONY: clean
clean:
	rm -f *.o a.out core main arm_interface
.PHONY: all
all: clean main

