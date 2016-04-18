#include <iostream>
#include "../interface/SKMatrix.hpp"
#include <armadillo>

using namespace arma; 
class Armadillo_Matrix: SKMatrix<Armadillo_Matrix, arma::mat>{
    mat matrix_data; 
    public:
    int row; 
    int col; 

    Armadillo_Matrix(){}
    Armadillo_Matrix(int r, int c){
        matrix_data = mat(r,c); 
    }

    mat data() const{
        return mat(matrix_data); 
    }

    int get_row(){
        return row; 
    }

};

