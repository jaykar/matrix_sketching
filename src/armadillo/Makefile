CXX = g++
CXXFLAGS = -std=c++14 -Wall -g -I../interface/

arm_interface: arm_interface.o
	$(CXX) $(CXXFLAGS) -o arm_interface arm_interface.o

arm_interface.o: Arm_interface.cpp 
	$(CXX) $(CXXFLAGS) -c Arm_interface.cpp
.PHONY: clean
clean:
	rm -f *.o a.out core main arm_interface
.PHONY: all
all: clean main

