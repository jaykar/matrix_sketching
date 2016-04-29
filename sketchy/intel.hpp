#ifndef __INTEL_H__
#define __INTEL_H__

#include <sstream>
#include <chrono>
#include "SKMatrix.hpp"
#include "mkl.h"

namespace sketchy {
    class intel : SKMatrix<intel, float *>  {
        private:
            int rows;
            int cols;

        public:
            intel() {
                matrix_data = NULL;
                rows = 0;
                cols = 0;
            }

            intel(const int row, const int col) {
                matrix_data = (float *) mkl_malloc(row * col * sizeof(float), sizeof(float));
                memset(matrix_data, 0, row * col * sizeof(float));
                rows = row;
                cols = col;
            }

            void identity(const int n) {
                *this = intel(n, n);
                int i;
                for(i = 0; i < n; i++)
                    matrix_data[i*cols + i] = 1.0;
            }

            void set(int r, int c, float a) {
                this->matrix_data[r*this->rows + c] = a;
            }

            ~intel() {
                mkl_free(matrix_data);
            }

            intel(const intel& im) {
                matrix_data = (float *) mkl_malloc(im.size() * sizeof(float), sizeof(float));
                cblas_scopy(im.size(), im.matrix_data, 1, matrix_data, 1);
                rows = im.rows;
                cols = im.cols;
            }

            intel& operator=(const intel& rhs) {
                if(this->matrix_data)
                    mkl_free(this->matrix_data);
                this->rows = rhs.rows;
                this->cols = rhs.cols;
                this->matrix_data = (float *) mkl_malloc(rhs.size() * sizeof(float), sizeof(float));
                cblas_scopy(rhs.size(), rhs.matrix_data, 1, this->matrix_data, 1);
                return *this;
            }

            void clear() {
                if(this->matrix_data)
                    memset(this->matrix_data, 0, this->rows * this->cols * sizeof(float));
            }

            int size() const {
                return this->rows * this->cols;
            }

            int num_rows(void) const {
                return this->rows;
            }

            int num_cols(void) const {
                return this->cols;
            }

            float *data() const {
                return this->matrix_data;
            }

            intel rand_n(const int row, const int col) {
                 if (row < 0 || col < 0) {
                    throw std::range_error("Column and row lengths must be non-negative integers");
                } else {
                    if(this->matrix_data)
                        mkl_free(this->matrix_data);
                    this->matrix_data = (float *) mkl_malloc(this->size() * sizeof(float), sizeof(float));
                    int i;
                    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
                    std::mt19937 gen(seed);
                    std::normal_distribution<float> n;

                    for(i = 0; i < this->size(); i++) {
                        this->matrix_data[i] = n(gen);
                    }

                    return *this;
                }
            }

            intel mult(const intel& rhs) const {
                if (rhs.num_rows() != num_cols()) {
                    throw std::range_error("Column of left matrix does not match row of right matrix");
                } else {
                    intel product(this->rows, rhs.cols);
                    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, this->rows,
                            rhs.cols, this->cols, 1, this->matrix_data, rhs.rows,
                            rhs.matrix_data, rhs.cols, 0, product.matrix_data, rhs.cols);
                    return product;
                }
            }

            intel elem_div(const float a) const {
                if(a == 0) {
                    throw std::overflow_error("Cannot divide by 0" );
                } else{
                    int i;
                    intel temp(*this);
                    for(i = 0; i < this->size(); i++)
                        temp.matrix_data[i] = this->matrix_data[i] / a;
                    return temp;
                }
            }

            intel concat(const intel& col) const {
                if (col.num_rows() != this->num_rows()) {
                    throw std::range_error("Number of rows do not match");
                } else {
                    intel ret(this->rows, this->cols + col.cols);
                    float *current = ret.matrix_data;
                    int i;
                    for(i = 0; i < this->rows; i++) {
                        cblas_scopy(this->cols, this->matrix_data + i * this->cols, 1, current, 1);
                        current += this->cols;
                        cblas_scopy(col.cols, col.matrix_data + i * col.cols, 1, current, 1);
                        current += col.cols;
                    }
                    return ret;
                }
            }

            intel solve_x(const intel& B) const {
                if(this->rows != this->cols)
                    throw std::range_error("A must be square in Ax = B");
                if(this->cols != B.rows)
                    throw std::range_error("B.rows must equal A.cols in Ax = B");
                if(B.cols != 1)
                    throw std::range_error("B must be nx1 vector in Ax = B");
                intel X(B.rows, 1);
                MKL_INT lda = this->rows;
                MKL_INT n = this->cols;
                MKL_INT nrhs = 1;
                MKL_INT ldb = nrhs;
                MKL_INT ipiv[n];

                intel copy(B);
                LAPACKE_sgesv(LAPACK_ROW_MAJOR, n, nrhs, this->matrix_data, lda, ipiv, copy.matrix_data, ldb);

                return copy;
            }

            intel get_col(const int n) const {
                if(col_n < 0 || col_n >= num_cols()) {
                    throw std::range_error("Column index out of bound");
                } else {
                    return get_cols(n, n);
                }
            }

            intel get_cols(const int start, const int end) const {
                if (start < 0 || end > num_cols()) {
                    throw std::range_error("Column index out of bound");
                    throw;
                } else if (start > end){
                    throw std::range_error("Start column greater than end column");
                } else {
                    intel ret(this->rows, end - start + 1);
                    int i;
                    float *current = ret.matrix_data;
                    for(i = 0; i < this->rows; i++) {
                        cblas_scopy(end - start + 1, this->matrix_data + i * this->cols + start, 1, current, 1);
                        current += end - start + 1;
                    }
                    return ret;
                }
            }

            void transpose() {
                mkl_simatcopy('r', 't', this->rows, this->cols, 1.0, this->matrix_data,
                        this->cols, this->rows);
                int i = this->rows;
                this->rows = this->cols;
                this->cols = i;
            }

            intel subtract(const intel& rhs) const {
                if(rhs.num_rows() != this->num_rows()){
                    throw std::range_error("Number of rows do not match");
                } else {
                    intel temp(*this);
                    cblas_saxpy(temp.size(), -1, rhs.matrix_data, 1, temp.matrix_data, 1);
                    return temp;
                }
            }

            float accumulate() const {
                return cblas_sasum(this->size(), this->matrix_data, 1);
            }

            void qr_decompose(intel& Q, intel& R) const {
            }

            void svd(intel& U, intel& S, intel& V, const int k) const {
                // return std::vector<intel>();
            }

            std::vector<int> dimensions() const {
                std::vector<int> ret = {this->rows, this->cols};
                return ret;
            }

            intel add(const intel& rhs) const {
                if(this->dimensions() != rhs.dimensions()) {
                    throw std::range_error("mismatched dimensions");
                } else {
                    intel temp(*this);
                    cblas_saxpy(temp.size(), 1, rhs.matrix_data, 1, temp.matrix_data, 1);
                    return temp;
                }
            }

            intel operator+(const intel& rhs) const {
                return this->add(rhs);
            }

            intel& operator+=(const intel& rhs) {
                *this = *this + rhs;
                return *this;
            }

            intel operator-(const intel& rhs) const {
                return this->subtract(rhs);
            }

            intel& operator-=(const intel& rhs) {
                *this = *this - rhs;
                return *this;
            }

            intel operator*(const intel& rhs) const {
                return this->mult(rhs);
            }

            intel& operator*=(const intel& rhs) {
                *this = this->mult(rhs);
                return *this;
            }

            intel operator/(const float a) {
                return this->elem_div(a);
            }

            intel& operator/=(const float a) {
                *this = this->elem_div(a);
                return *this;
            }

            friend bool operator==(const intel& lhs, const intel& rhs);
            friend std::ostream& operator<<(std::ostream&os, const intel& im);
    };

    bool operator==(const intel& lhs, const intel& rhs) {
        if(rhs.rows != lhs.rows || rhs.cols != lhs.cols)
            return false;
        int i;
        for(i = 0; i < lhs.size(); i++)
            if(lhs.matrix_data[i] != rhs.matrix_data[i])
                return false;
        return true;
    }

    std::ostream& operator<<(std::ostream&os, const intel& im) {
        int i;
        std::ostringstream out;
        for(i = 0; i < im.size(); i++) {
            if(i > 0 && i % im.cols == 0)
                out << '\n';
            out << std::to_string(im.matrix_data[i]);
            out << '\t';
        }
        os << out.str();
        return os;
    }
}

#endif