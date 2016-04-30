#ifndef __ARMADILLO_H__
#define __ARMADILLO_H__

#include <armadillo>
#include "SKMatrix.hpp"
#include <string>
using namespace arma;

namespace sketchy {
    class armadillo: public SKMatrix<armadillo, arma::mat>{
        public:
            armadillo(){
                matrix_data = mat();
            }

            armadillo(std::string filename){
                mat A;
                A.load(filename);
                matrix_data = A;
            }

            void save(std::string fn){
                matrix_data.save(fn, raw_ascii);
            }

<<<<<<< HEAD
            /*
                          void _eye(const int n){
                          matrix_data = eye(n);
=======
            /*            
                          void _eye(const int n){
                          matrix_data = eye(n); 
>>>>>>> fd6e2ac697c2d4288085bb24cd3d6cd844e6ef6e

                          }
                          */

            std::vector<int> dimensions() const{
                std::vector<int> a;
                a.push_back(matrix_data.n_rows);
                a.push_back(matrix_data.n_cols);
                return a;
            }

            ~armadillo() = default;

            armadillo(int r, int c){
                if(r < 0 || c < 0) {
                    throw std::invalid_argument("Column and row lengths must be non-negative integers");
                } else {
                    matrix_data = mat(r,c);
                }
            }

            mat data() const{
                return mat(matrix_data);
            }

            armadillo(mat other){
                this->matrix_data = other;
            }

            armadillo(const armadillo& other){
                mat temp = other.matrix_data;
                this->matrix_data = temp;
            }

            armadillo& operator=(const armadillo& other){
                this->matrix_data = mat(other.matrix_data);
                return *this;
            }

            armadillo& operator=(const mat& other){
                this->matrix_data = mat(other);
                return *this;
            }

            armadillo(armadillo&& other){
                this->matrix_data = other.matrix_data;
            }

            armadillo& operator=(armadillo&& other){
                this->matrix_data = other.matrix_data;
                return *this;
            }


            void clear(){
                this->matrix_data = mat();

            }

            int size() const{
                return this->matrix_data.n_elem;
            }

            int num_rows(void) const{
                return this->matrix_data.n_rows;
            }

            int num_cols(void) const{
                return this->matrix_data.n_cols;
            }

            armadillo mult(const armadillo& rhs) const{
                if (rhs.num_rows() != num_cols())
                    throw std::invalid_argument("Column of left matrix does not match row of right matrix");
                mat a = this->matrix_data * rhs.matrix_data;
                return a;
            }

            armadillo rand_n(int row, int col){
                if (row < 0 || col < 0)
                    throw std::invalid_argument("Column and row lengths must be non-negative integers");
                mat a;
                a.randn(row, col);
                this->matrix_data = a;
                return *this;
            }

            armadillo elem_div(const float a) const{
                if(a == 0)
                    throw std::overflow_error("Cannot divide by 0" );
                mat b = data();
                b = b/a;
                return armadillo(b);
            }


            armadillo concat(const armadillo& m) const{
                if (m.num_rows() != this->num_rows())
                    throw std::invalid_argument("Number of rows do not match");
                mat a(data());
                a.insert_cols(a.n_cols-1, m.matrix_data);
                return armadillo(a);
            }

            armadillo solve_x(const armadillo& B) const{
                if(this->num_rows() != this->num_cols())
                    throw std::invalid_argument("A must be square in Ax = B");
                if(this->num_cols() != B.num_rows())
                    throw std::invalid_argument("B.rows must equal A.cols in Ax = B");
                if(B.num_cols() != 1)
                    throw std::invalid_argument("B must be nx1 vector in Ax = B");

                mat X = solve(matrix_data, B.matrix_data);
                return armadillo(X);
            }

            armadillo get_cols(int start, int end) const{
                if (start < 0 || end > num_cols())
                    throw std::range_error("Column index out of bound");
                throw;
                if (start > end)
                    throw std::range_error("Start column greater than end column");
                mat a = matrix_data.cols(start, end);
                return armadillo(a);
            }


            armadillo get_col(int col_n) const{
                if(col_n < 0 || col_n >= num_cols())
                    throw std::range_error("Column index out of bound");
                mat a = matrix_data.col(col_n);
                return armadillo(a);
            }

            void transpose(){
                matrix_data = matrix_data.t();
            }

            armadillo subtract(const armadillo& rhs) const{
                if(rhs.num_rows() != this->num_rows())
                    throw std::invalid_argument("Number of rows do not match");
                mat a = matrix_data - rhs.matrix_data;
                return armadillo(a);
            }

            float accumulate() const{
                return accu(matrix_data);
            }

            void qr_decompose(armadillo& Q, armadillo& R) const{
                mat q1;
                mat r1;
                qr_econ(q1, r1, matrix_data);
                Q.matrix_data = q1;
                R.matrix_data = r1;
            }

            friend std::ostream& operator<<(std::ostream&os, const armadillo& out);

            void svd(armadillo& U, armadillo& S, armadillo& V, int k) const{
                mat u;
                vec s;
                mat v;

                //arma::svds(u, s, v, sp_mat(matrix_data), k);
                svd_econ(u, s, v, matrix_data);

                U.matrix_data = u;
<<<<<<< HEAD
                S.matrix_data = s; //diagmat(s.rows(0, k-1));
=======
                S.matrix_data = s; //diagmat(s.rows(0, k-1)); 
>>>>>>> fd6e2ac697c2d4288085bb24cd3d6cd844e6ef6e
                V.matrix_data = v;
                //S.matrix_data = mat(s);
            }
    };
    std::ostream& operator<<(std::ostream& os, const armadillo&out){
        return os << out.matrix_data ;
    }

}

#endif