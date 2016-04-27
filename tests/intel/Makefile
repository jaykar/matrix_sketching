CXX = g++
CXXFLAGS = -std=c++14 -Wall -g -I../interface/  -L$(MKLROOT)/../compiler/lib  -m64 -I$(MKLROOT)/include

LDFLAGS =  -L$(MKLROOT)/lib -Wl,-rpath,$(MKLROOT)/lib -L$(MKLROOT)/../compiler//lib -Wl,-rpath,$(MKLROOT)/../compiler/lib -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm

OBJECTS = main.o #intel.o

main: $(OBJECTS) intel.h ../interface/Matrix.hpp
	$(CXX) $(CXXFLAGS) -o main $(OBJECTS) $(LDFLAGS) 

#~/mkl/mac/mkl/lib/libmkl_core.dylib 

main.o: main.cpp 
	$(CXX) main.cpp -c $(CXXFLAGS) 

#intel.o: intel.cpp intel.h
#	$(CXX) $(CXXFLAGS) -c intel.cpp

.PHONY: clean
clean:
	rm -rf *.o main *.dSYM

.PHONY: all
all: clean main

valgrind: all
	valgrind 
