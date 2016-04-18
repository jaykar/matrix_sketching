#include <iostream>
#include "../interface/SKMatrix.hpp"
#include <armadillo>

using namespace arma; 
class Armadillo_Matrix: SKMatrix<Armadillo_Matrix, arma::mat>{
    private:
        mat *matrix_data; 
    public:

        Armadillo_Matrix(){
            matrix_data = new mat(); 
        }
        
        ~Armadillo_Matrix(){
            if (matrix_data){
                delete [] matrix_data; 
            }
        }

        Armadillo_Matrix(int r, int c){
            matrix_data = new mat(r,c); 
        }
        
        mat data() const{
            return mat(*matrix_data); 
        }
        
        Armadillo_Matrix(Armadillo_Matrix& other){
            auto temp = *other.matrix_data; 
            this->matrix_data = new mat(size(temp)); 
            *(this->matrix_data) = temp; 
        }
        
        Armadillo_Matrix(mat other){
            auto temp = other; 
            this->matrix_data = new mat(size(temp)); 
            *(this->matrix_data) = temp; 
        } 

        Armadillo_Matrix& operator=(const Armadillo_Matrix& other){
            std::cout << "using copy operator" << std::endl; 
            if (this != &other){
                if (this->matrix_data){
                    delete this->matrix_data; 
                }
                auto temp = *other.matrix_data; 
                this->matrix_data = new mat(size(temp)); 
                *(this->matrix_data) = temp; 
            }
            return *this; 
        }
        
        Armadillo_Matrix(Armadillo_Matrix&& other){
            this->matrix_data = other.matrix_data; 
            other.matrix_data = NULL; 
        }
        
        Armadillo_Matrix& operator=(Armadillo_Matrix&& other){
            std::cout << "using move operator" << std::endl; 
            if (this != &other){
                delete this->matrix_data; 
                this->matrix_data = other.matrix_data;
                other.matrix_data = NULL; 
            }
            return *this; 
        }

        Armadillo_Matrix rand_n(int row, int col, int mean, int std){
            mat a(row, col); 
            a.randn(); 
            a = a*std + mean; 
            auto b = Armadillo_Matrix(a); 
            return b; 
        }
};


