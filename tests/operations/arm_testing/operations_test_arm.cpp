#include <iostream>
#include "operations.hpp"
//#include "boost.hpp"
#include "armadillo.hpp"
#include <string>
#include <chrono>

using namespace std;

//namespace bnu = boost::numeric::ublas;

int main(){

    string x_label = "comple1x_g"; 
    //string x_label = "data_image"; 
    cout << "loading x train" << endl; 
    sketchy::armadillo x_train = sketchy::armadillo(x_label); 
    int sketch_size = 400; 

    cout << "doing k-svd" << endl; 
    auto t1 = chrono::high_resolution_clock::now();
    auto vec_ans = sketchy::ops::k_svd(x_train, sketch_size); 
    auto t2 = chrono::high_resolution_clock::now();
    auto time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
    std::cout << "seconds: " << time_span.count() << " seconds.";
    std::cout << std::endl;

    cout << "saving" << endl; 
    vec_ans[0].save("u"); 
    vec_ans[1].save("s"); 
    vec_ans[2].save("v"); 
    return 0;
}
