CXX = g++
CXXFLAGS = -std=c++14 -Wall -g

BOOST_LIB_PATH = -I /usr/local/include/

# INTERFACE = -I ../interface
BOOST_INTERFACE = -I ../../src/

#-O3 -D -DNDEBUG -DBOOST_UBLAS_NDEBUG
boost_wrapper_test: boost_wrapper_test.o
	$(CXX) $(CXXFLAGS) -o boost_wrapper_test boost_wrapper_test.o

boost_wrapper_test.o: boost_wrapper_test.cpp
	$(CXX) $(CXXFLAGS) $(BOOST_LIB_PATH) $(BOOST_INTERFACE) -c boost_wrapper_test.cpp

.PHONY: clean
clean:
	rm -f *.o a.out core main boost_wrapper_test
.PHONY: all
all: clean main