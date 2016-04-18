#include <iostream>
#include "SKMatrix.hpp"

class Armadillo_Matrix: SKMatrix<Armadillo_Matrix>{
    int row; 
    int col; 
    int data; 
    public:
        Armadillo_Matrix(int r, int c){
            row = r; 
            col = c; 
        }
        
        Armadillo_Matrix& mult(Armadillo_Matrix &a) const{
            return new Armadillo_Matrix(3,3); 
        }

        int get_row(){
            return row; 
        }

};

using namespace std;
int main(){

    auto a = Armadillo_Matrix(3,4); 
    cout << a.get_row() << endl;
    return 0;
}
