#include "../armadillo/Arm_interface.hpp"
#include <iostream>
#include <math.h>

using namespace std; 

template <class T>
T gaussian_projection(T input, int s){
    int n_col = input.dimensions()[1]; 
    auto a = T(); 
    a.rand_n(n_col, s); 
    a.elem_div(sqrt(s*1.0));
    return input.mult(a); 
}

int main(){
    auto input = Armadillo_Matrix(); 
    input.rand_n(10, 10); 
    auto a = gaussian_projection<Armadillo_Matrix>(input, 4); 
    cout << a.data() << endl;
    return 0; 
}
