CXX = g++
CXXFLAGS = -std=c++14 -Wall -g

BOOST_LIB_PATH = -I /usr/local/include/
BOOST_COMPILE_OPTION = -O3 -DNDEBUG -DBOOST_UBLAS_NDEBUG
LAPACK_LIB_PATH = -L /usr/local/lib/ -llapack -L /usr/local/lib/ -lblas
# INTERFACE = -I ../interface
INTERFACE = -I ../../sketchy/

#-O3 -D -DNDEBUG -DBOOST_UBLAS_NDEBUG
boost_wrapper_test: boost_wrapper_test.o
	$(CXX) $(CXXFLAGS) $(LAPACK_LIB_PATH) $(BOOST_COMPILE_OPTION) -o boost_wrapper_test boost_wrapper_test.o

boost_wrapper_test.o: boost_wrapper_test.cpp
	$(CXX) $(CXXFLAGS) $(BOOST_LIB_PATH) $(INTERFACE) -c boost_wrapper_test.cpp

.PHONY: clean
clean:
	rm -f *.o a.out core main boost_wrapper_test
.PHONY: all
all: clean main