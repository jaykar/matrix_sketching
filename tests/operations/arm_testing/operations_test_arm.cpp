#include <iostream>
#include "operations.hpp"
//#include "boost.hpp"
#include "armadillo.hpp"
#include <string>
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

    string x_label = "data_image"; 
    cout << "loading x train" << endl; 
    sketchy::armadillo x_train = sketchy::armadillo(x_label); 
    auto vec_ans = sketchy::ops::k_svd(x_train, 5, 20); 

    vec_ans[0].save("u"); 
    vec_ans[1].save("s"); 
    vec_ans[2].save("v"); 

    /*
    cout << "loading y train" << endl; 
    string y_label = "../../../../mnist/regularized_label_small"; 
    sketchy::armadillo y_train = sketchy::armadillo(y_label); 
    
    cout << "solving for b" << endl;  
    sketchy::armadillo b = x_train.solve_x(y_train);  
    
    cout << "solving for xtilde" << endl;  
    sketchy::armadillo x_tilde = sketchy::ops::lin_regress<sketchy::armadillo>(x_train, y_train, 600); 
            
    //sketchy::armadillo y_pred = x_train.mult(x_tilde); 
    //sketchy::armadillo y_pred2 = x_train.mult(b); 

    b.save("b"); 
    x_tilde.save("x_tilde"); 
    */
    //cout << (y_pred.subtract(y_pred2)).accumulate()  << endl; 

    return 0;
}
