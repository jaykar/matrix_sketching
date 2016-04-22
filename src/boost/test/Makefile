CXX = g++
CXXFLAGS = -std=c++14 -Wall -g -I /usr/local/include/

boost_wrapper_test: boost_wrapper_test.o
	$(CXX) $(CXXFLAGS) -o boost_wrapper_test boost_wrapper_test.o

boost_wrapper_test.o: boost_wrapper_test.cpp
	$(CXX) $(CXXFLAGS) -c boost_wrapper_test.cpp

.PHONY: clean
clean:
	rm -f *.o a.out core main boost_wrapper_test
.PHONY: all
all: clean main