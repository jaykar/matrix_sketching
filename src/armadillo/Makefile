CXX = g++
CXXFLAGS = -std=c++14 -Wall -g -I../interface/
LDLIBS = -larmadillo

arm_interface: arm_interface.o
	$(CXX) $(CXXFLAGS) -o arm_interface arm_interface.o $(LDLIBS)

arm_interface.o: Arm_interface.cpp 
	$(CXX) $(CXXFLAGS) -c Arm_interface.cpp ../interface/SKMatrix.hpp

.PHONY: clean
clean:
	rm -f *.o a.out core main arm_interface
.PHONY: all
all: clean main

