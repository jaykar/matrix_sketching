CXX = g++
CXXFLAGS = -std=c++14 -Wall -g -I../interface/ -I/Users/adam/mkl/mac/mkl/include  -I/Users/adam/mkl/mac/mkl/lib
LDFLAGS = /Users/adam/mkl/mac/mkl/lib/*.a -rpath /Users/adam/mkl/mac/mkl/lib 

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
