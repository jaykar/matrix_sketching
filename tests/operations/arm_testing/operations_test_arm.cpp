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
    sketchy::armadillo x_train = sketchy::armadillo(x_label); 

    //parameters for sketching; 
    sketchy::armadillo U;
    sketchy::armadillo S;
    sketchy::armadillo V;
    int sketch_size = 0; 
    vector<int> sizes = {1, 5, 10, 15, 20}; 
    int resolution = 20; 
    int tenth = x_train.num_cols() / resolution; 
    for(int i=1; i< resolution ; i++){
        int size = i*tenth; 
        sizes.push_back(size); 
    }
    
    sizes.push_back(x_train.num_cols());
    int n_rep = 5; 

    for (auto i: sizes){
        auto avg_time = 0.0; 
        auto avg_diff = 0.0; 
        for (int j = 0; j<n_rep; j++){
            sketch_size = i; 
            auto t1 = chrono::high_resolution_clock::now();
            sketchy::ops::k_svd(x_train, U, S, V, sketch_size); 
            auto t2 = chrono::high_resolution_clock::now();
            avg_time += chrono::duration_cast<chrono::duration<double>>(t2 - t1).count(); 
            
            sketchy::armadillo temp = U.mult(sketchy::armadillo(arma::diagmat(S.data()))); 
            V.transpose(); 
            sketchy::armadillo reconstruct = temp.mult(V); 
            sketchy::armadillo diff = x_train.subtract(reconstruct); 
            avg_diff += arma::norm(diff.data(), "fro")/arma::norm(x_train.data(), "fro"); 

        }

        avg_time /= n_rep; 
        avg_diff /= n_rep; 
        cout <<i << " " << avg_time << " " << avg_diff << endl;
    }

   
    
    auto avg_time = 0.0;  
    auto avg_diff = 0.0; 
    for(int i=0; i < n_rep; i ++){
        auto t1 = chrono::high_resolution_clock::now();
        x_train.svd(U, S, V, sketch_size); 
        auto t2 = chrono::high_resolution_clock::now();
        auto time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
        avg_time += time_span.count(); 


        sketchy::armadillo temp = U.mult(sketchy::armadillo(arma::diagmat(S.data()))); 
        V.transpose(); 
        sketchy::armadillo reconstruct = temp.mult(V); 
        sketchy::armadillo diff = x_train.subtract(reconstruct); 
        avg_diff += arma::norm(diff.data(), "fro")/arma::norm(x_train.data(), "fro"); 
    }
    
    avg_time /= n_rep; 
    cout <<"svd " << avg_time << " " << avg_diff << endl; 
    return 0;
}
