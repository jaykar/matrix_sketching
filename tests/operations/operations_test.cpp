#include <iostream>
#include "operations.hpp"
//#include "boost.hpp"
#include "armadillo.hpp"

using namespace std;

//namespace bnu = boost::numeric::ublas;

template<typename T>
T rand_n(int total_size){
    auto input = T();
    input.rand_n(total_size, total_size);
    return input;
}

template<typename T>
T transpose_mult(T input){
    auto b = input;
    b.transpose();
    return input.mult(b);
}

int main(){
    int total_size = 10;
    int sketchyetch_size = 4;
    //sketchy::arm X_train = sketchy::arm("../../X_train.mat"); 
    sktechy::arm X; 
    X.rand_n(30,30); 
    cout << X << endl; 
    return 0;
}
