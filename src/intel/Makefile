CXX = g++
CXXFLAGS = -std=c++14 -Wall -g -I../interface/  -L$(MKLROOT)/../compiler/lib  -m64 -I$(MKLROOT)/include

LDFLAGS =  -L$(MKLROOT)/lib -Wl,-rpath,$(MKLROOT)/lib -L$(MKLROOT)/../compiler//lib -Wl,-rpath,$(MKLROOT)/../compiler/lib -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -lmkl_scalapack_lp64 -lmkl_blacs_mpich_lp64 
#LDFLAGS = -L$(MKLROOT)/lib -Wl,-rpath,$(MKLROOT)/lib -lmkl_scalapack_lp64 -lmkl_blacs_mpich_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm

OBJECTS = main.o #sk_intel.o

main: $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o main $(OBJECTS) $(LDFLAGS) 

#~/mkl/mac/mkl/lib/libmkl_core.dylib 

main.o: main.cpp sk_intel.hpp  ../interface/SKMatrix.hpp
	$(CXX) main.cpp -c $(CXXFLAGS) 

#sk_intel.o: sk_intel.cpp sk_intel.h
#	$(CXX) $(CXXFLAGS) -c sk_intel.cpp

.PHONY: clean
clean:
	rm -rf *.o main *.dSYM

.PHONY: all
all: clean main

valgrind: all
	valgrind 
