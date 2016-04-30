#include <iostream>
#include "operations.hpp"
//#include "boost.hpp"
#include "armadillo.hpp"
#include <string>
#include <chrono>

using namespace std;

int main(){
    //load file
    string x_label = "paris"; 
    cout << "loading x train" << endl; 
    sketchy::armadillo x_train = sketchy::armadillo(x_label); 
    //sketchy::armadillo ans = sketchy::ops::count_sketch(x_train, 10); 
    //parameters for sketching; 
    
    int sketch_size = 700; 
    sketchy::armadillo U;
    sketchy::armadillo S;
    sketchy::armadillo V;

    cout << "doing k-svd" << endl; 
    auto t1 = chrono::high_resolution_clock::now();
    sketchy::ops::k_svd(x_train, U, S, V, sketch_size); 
    auto t2 = chrono::high_resolution_clock::now();
    auto time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
    cout << "seconds: " << time_span.count() << " seconds." << endl;

    cout << "saving" << endl; 
    U.save("u");
    S.save("s");
    V.save("v");
           
    return 0;
}
