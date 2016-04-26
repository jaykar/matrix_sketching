#ifndef __ARM_INTERFACE_H__
#define __ARM_INTERFACE_H__

#include "../interface/SKMatrix.hpp"
#include <armadillo>
#include <vector>
using namespace arma;
class sk_arm: SKMatrix<sk_arm, arma::mat>{
    private:
        mat matrix_data;
    public:
        sk_arm(){
            matrix_data = mat();
        }
        
        std::vector<int> dimensions() const{
            auto a = std::vector<int>();
            a.push_back(matrix_data.n_rows);
            a.push_back(matrix_data.n_cols);
            return a;
        }
        ~sk_arm() = default;

        sk_arm(int r, int c){
            matrix_data = mat(r,c);
        }
        
        mat data() const{
            return mat(matrix_data); 
        }

        mat& data(){
            return matrix_data; 
        }

        sk_arm(mat other){
            this->matrix_data = other;
        }
        

        sk_arm(const sk_arm& other){
            std::cout << "using copy constructor" << std::endl;
            auto temp = other.matrix_data;
            //this->matrix_data = mat(temp);
            this->matrix_data = temp;
        }
        
        sk_arm& operator=(const sk_arm& other){
            std::cout << "using copy assignment" << std::endl;
            this->matrix_data = mat(other.matrix_data);
            return *this;
        }

        sk_arm& operator=(const mat& other){
            std::cout << "using copy assignmnet for mat" << std::endl;
            this->matrix_data = mat(other);
            return *this;
        }

        sk_arm(sk_arm&& other){
            std::cout << "using move constructor" << std::endl;
            this->matrix_data = other.matrix_data;
            //other.matrix_data = mat();
        }

        sk_arm& operator=(sk_arm&& other){
            std::cout << "using move assignment" << std::endl;
            this->matrix_data = other.matrix_data;
            //other.matrix_data = mat();
            return *this;
        }
        

        void clear(){
            this->matrix_data = mat(); 

        }

        int size() const{
            return this->matrix_data.n_elem;
        }

        int num_rows() const{
            return this->matrix_data.n_rows; 
        }

        int num_cols() const{
            return this->matrix_data.n_cols; 
        }

        sk_arm mult(const sk_arm& rhs) const{
            mat a = this->matrix_data * rhs.matrix_data;
            return a; 
        }

        sk_arm rand_n(int row, int col){
            mat a; 
            a.randn(row, col);
            //a = a*std + mean;
            this->matrix_data = a;
            return *this;
        }

        sk_arm elem_div(const double b) const{
           mat a = data();
           if (b != 0.0){
               a = a/b;
           }
           else{
               std::cout << "cannot divide by 0" << std::endl;
           }
           return sk_arm(a); 
        }

        
        sk_arm concat(const sk_arm& column) const{
            mat a = data(); 
            a.insert_cols(a.n_cols-1, column.matrix_data); 
            return sk_arm(a); 
        }

        sk_arm solve_x(const sk_arm& B) const{
                auto X = solve(matrix_data, B.matrix_data);
                return sk_arm(X); 
        }

        sk_arm get_cols(int start, int end) const{
            mat a = matrix_data.cols(start, end);
            return sk_arm(a); 
        }


        sk_arm get_col(int col_n) const{
            mat a = matrix_data.col(col_n);
            return sk_arm(a); 
        }

        void transpose(){
            matrix_data = matrix_data.t();
        }

        sk_arm subtract(const sk_arm& rhs) const{
            mat a = matrix_data - rhs.matrix_data;
            return sk_arm(a); 
        }

        double accumulate() const{
            return accu(matrix_data);
        }

        void qr_decompose(sk_arm& a, sk_arm& b) const{
                mat Q;
                mat R; 
                qr(Q, R, matrix_data); 
                a.matrix_data = Q; 
                b.matrix_data = R; 
        }
        
        friend std::ostream& operator<<(std::ostream&os, const sk_arm& out);

        std::vector<sk_arm> svds(int k){
                mat U;
                vec s;
                mat V;
                arma::svds(U, s, V, sp_mat(matrix_data), k);
                auto ans = std::vector<sk_arm>(3); 
                ans[0] = sk_arm(U); 
                ans[1] = sk_arm(mat(s)); 
                ans[2] = sk_arm(V); 
                return ans; 
        }
};

std::ostream& operator<<(std::ostream& os, const sk_arm&out)
{
    return os << out.matrix_data ;
}
#endif
