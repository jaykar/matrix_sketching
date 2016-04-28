CXX = g++
CXXFLAGS = -std=c++14 -Wall -g

BOOST_LIB_PATH = -I /usr/local/include/
BOOST_COMPILE_OPTION = -O3 -D NDEBUG -DBOOST_UBLAS_NDEBUG
LAPACK_LIB_PATH = -L /usr/local/lib/ -llapack -L /usr/local/lib/ -lblas
# INTERFACE = -I ../interface
BOOST_INTERFACE = -I ../../src/

#-O3 -D -DNDEBUG -DBOOST_UBLAS_NDEBUG
time: time.o
	$(CXX) $(CXXFLAGS) $(LAPACK_LIB_PATH) $(BOOST_COMPILE_OPTION) -o time time.o

time.o: time.cpp
	$(CXX) $(CXXFLAGS) $(BOOST_LIB_PATH) $(BOOST_INTERFACE) -c time.cpp

.PHONY: clean
clean:
	rm -f *.o a.out core main time
.PHONY: all
all: clean main