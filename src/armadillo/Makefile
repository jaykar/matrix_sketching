CXX = g++
CXXFLAGS = -std=c++14 -Wall -g -I../interface/
LDLIBS = -larmadillo

arm_interface.o: Arm_interface.hpp 
	$(CXX) $(CXXFLAGS) -c Arm_interface.hpp ../interface/SKMatrix.hpp

.PHONY: clean
clean:
	rm -f *.o a.out core main arm_interface
.PHONY: all
all: clean main

