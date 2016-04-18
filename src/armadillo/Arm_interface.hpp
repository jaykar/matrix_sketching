#include <iostream>
#include "../interface/SKMatrix.hpp"
#include <armadillo>

using namespace arma; 
class Armadillo_Matrix: SKMatrix<Armadillo_Matrix, arma::mat>{
    private:
        mat matrix_data; 
    public:

        Armadillo_Matrix(){
            //std::cout << "constructor" << std::endl; 
            matrix_data = mat(); 
        }

        ~Armadillo_Matrix(){
        }

        Armadillo_Matrix(int r, int c){
            matrix_data = mat(r,c); 
        }

        mat data() const{
            return mat(matrix_data); 
        }


        Armadillo_Matrix(mat other){
            this->matrix_data = other; 
        } 

        Armadillo_Matrix(const Armadillo_Matrix& other){
            //std::cout << "using copy operator" << std::endl; 
            auto temp = other.matrix_data; 
            this->matrix_data = mat(temp); 
        }

        Armadillo_Matrix& operator=(const Armadillo_Matrix& other){
            std::cout << "using copy operator" << std::endl; 
            this->matrix_data = mat(other.matrix_data); 
            return *this; 
        }

        Armadillo_Matrix(Armadillo_Matrix&& other){
            std::cout << "using move operator" << std::endl; 
            this->matrix_data = other.matrix_data; 

            //might not need to do this
            other.matrix_data = mat(); 
        }

        Armadillo_Matrix& operator=(Armadillo_Matrix&& other){
            std::cout << "using move operator" << std::endl; 
            this->matrix_data = other.matrix_data;
            other.matrix_data = mat(); 
            return *this; 
        }

        Armadillo_Matrix rand_n(int row, int col, int mean, int std); 

        Armadillo_Matrix mult(Armadillo_Matrix& rhs){
            mat a = this->matrix_data * rhs.matrix_data;
            return std::move(Armadillo_Matrix(a)); 
        }
};

Armadillo_Matrix Armadillo_Matrix::rand_n(int row, int col, int mean, int std){
    std::cout << "filling with gaussian noise" << std::endl; 
    mat a(row, col); 
    a.randn(); 
    a = a*std + mean; 
    auto b = Armadillo_Matrix(a); 
    return std::move(b); 
}


