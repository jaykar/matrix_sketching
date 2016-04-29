#ifndef __ARM_INTERFACE_H__
#define __ARM_INTERFACE_H__

#include "SKMatrix.hpp"
#include <armadillo>

using namespace arma;

namespace sketchy {
    class armadillo: SKMatrix<armadillo, arma::mat>{
        public:
            armadillo(){
                matrix_data = mat();
            }

            std::vector<int> dimensions() const{
                auto a = std::vector<int>();
                a.push_back(matrix_data.n_rows);
                a.push_back(matrix_data.n_cols);
                return a;
            }

            ~armadillo() = default;

            armadillo(int r, int c){
                matrix_data = mat(r,c);
            }

            mat data() const{
                return mat(matrix_data);
            }

<<<<<<< HEAD:src/sk_arm.hpp
            /*
            mat& data(){
                return matrix_data;
            }
            */

            arm(mat other){
=======
            armadillo(mat other){
>>>>>>> 723f20fd38bec2016d27aec829af4fcb0c49d443:sketchy/armadillo.hpp
                this->matrix_data = other;
            }

            armadillo(const armadillo& other){
                auto temp = other.matrix_data;
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

            int num_rows() const{
                return this->matrix_data.n_rows;
            }

            int num_cols() const{
                return this->matrix_data.n_cols;
            }

            armadillo mult(const armadillo& rhs) const{
                if (rhs.num_rows() != num_cols()) {
                    throw std::range_error("Column of left matrix does not match row of right matrix");
                } else {
                    mat a = this->matrix_data * rhs.matrix_data;
                    return a;
                }
            }

            armadillo rand_n(int row, int col){
                if (row < 0 || col < 0) {
                    throw std::range_error("Column and row lengths must be non-negative integers");
                } else {
                    mat a;
                    a.randn(row, col);
                    this->matrix_data = a;
                    return *this;
                }
            }

            armadillo elem_div(const double a) const{
                if(a == 0) {
                    throw std::overflow_error("Cannot divide by 0" );
                } else{
                    mat b = data();
                    b = b/a;
                    return armadillo(b);
                }
            }


<<<<<<< HEAD:src/sk_arm.hpp
            arm concat(const arm& m) const{
                if (m.num_rows() != this->num_rows()) {
                    throw std::range_error("Number of rows do not match");
                } else {
                    mat a(data());
                    a.insert_cols(a.n_cols-1, m.matrix_data);
                    return arm(a);
=======
            armadillo concat(const armadillo& m) const{
                if (m.num_rows() != this->num_rows()) {
                    throw std::range_error("Number of rows do not match");
                } else {
                    m a = data();
                    a.insert_cols(a.n_cols-1, m.matrix_data);
                    return armadillo(a);
>>>>>>> 723f20fd38bec2016d27aec829af4fcb0c49d443:sketchy/armadillo.hpp
                }
            }

            armadillo solve_x(const armadillo& B) const{
                auto X = solve(matrix_data, B.matrix_data);
                return armadillo(X);
            }

            armadillo get_cols(int start, int end) const{
                mat a = matrix_data.cols(start, end);
                return armadillo(a);
            }


            armadillo get_col(int col_n) const{
                mat a = matrix_data.col(col_n);
                return armadillo(a);
            }

            void transpose(){
                matrix_data = matrix_data.t();
            }

            armadillo subtract(const armadillo& rhs) const{
                if(rhs.num_rows() != this->num_rows()){
                    throw std::range_error("Number of rows do not match");
                } else {
                    mat a = matrix_data - rhs.matrix_data;
                    return armadillo(a);
                }
            }

            float accumulate() const{
                return accu(matrix_data);
            }

            void qr_decompose(armadillo& Q, armadillo& R) const{
                mat q;
                mat r;
                qr(q, r, matrix_data);
                Q.matrix_data = q;
                R.matrix_data = r;
            }

            friend std::ostream& operator<<(std::ostream&os, const armadillo& out);

<<<<<<< HEAD:src/sk_arm.hpp
            void svd(arm& U, arm& S, arm& V, const int k) const{
=======
            void svd(armadillo& U, armadillo& S, armadillo& V, const int k){
>>>>>>> 723f20fd38bec2016d27aec829af4fcb0c49d443:sketchy/armadillo.hpp
                mat u;
                vec s;
                mat v;
                arma::svds(u, s, v, sp_mat(matrix_data), k);
                auto ans = std::vector<armadillo>(3);
                U.matrix_data = u;
                S.matrix_data = mat(s);
                V.matrix_data = v;
            }
    };

    std::ostream& operator<<(std::ostream& os, const armadillo&out){
        return os << out.matrix_data ;
    }
}
#endif
