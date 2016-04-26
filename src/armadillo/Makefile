CXX = g++
CXXFLAGS = -std=c++14 -Wall -g -I../interface/ -I../libraries/include
LDLIBS = -larmadillo

arm_interface.o: sk_arm.hpp 
	$(CXX) $(CXXFLAGS) -c sk_arm.hpp ../interface/SKMatrix.hpp

.PHONY: clean
clean:
	rm -f *.o a.out core main arm_interface
.PHONY: all
all: clean main

