CXX = g++
CXXFLAGS = -std=c++11 -g -I../../sketchy/ -I/home/adam/NVIDIA_CUDA-7.5_Samples/common/inc/ -I/usr/local/cuda-7.5/include -m64 #-I$(MKLROOT)/include
CUDA_PATH ?= /usr/local/cuda-7.5
NVCC = $(CUDA_PATH)/bin/nvcc -ccbin $(CXX)
LDFLAGS = -L/usr/local/cuda-7.5/lib64/ -lcublas

OBJECTS = main.o 

main: $(OBJECTS)
	export LD_LIBRARY_PATH=/usr/local/cuda-7.5/lib64
	$(NVCC) $(CXXFLAGS) -o main $(OBJECTS) $(LDFLAGS)

#~/mkl/mac/mkl/lib/libmkl_core.dylib

main.o: main.cpp ../../sketchy/cuda.hpp  ../../sketchy/SKMatrix.hpp
	$(NVCC) main.cpp -c $(CXXFLAGS)

#sk_intel.o: sk_intel.cpp sk_intel.h
#	$(CXX) $(CXXFLAGS) -c sk_intel.cpp

.PHONY: clean
clean:
	rm -rf *.o main *.dSYM

.PHONY: all
all: clean main

valgrind: all
	valgrind
