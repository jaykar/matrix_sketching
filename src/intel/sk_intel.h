#include "SKMatrix.hpp"
#include "mkl.h"
#include <sstream>
#include <chrono>

class sk_intel : SKmatrix<sk_intel, float *, float>  {
    private:
        int rows;
        int cols;

    public:
        sk_intel() {
            matrix_data = NULL;
            rows = 0;
            cols = 0;
        }

        sk_intel(const int row, const int col) {
            matrix_data = (float *) mkl_malloc(row * col * sizeof(float), sizeof(float));
            memset(matrix_data, 0, row * col * sizeof(float));
            rows = row;
            cols = col;
        }

        void identity(const int n) {
            *this = sk_intel(n, n);
            int i;
            for(i = 0; i < n; i++)
                matrix_data[i*cols + i] = 1.0;
        }

        ~sk_intel() {
            mkl_free(matrix_data);
        }

        sk_intel(const sk_intel& im) {
            matrix_data = (float *) mkl_malloc(im.size() * sizeof(float), sizeof(float));
            cblas_scopy(im.size(), im.matrix_data, 1, matrix_data, 1);
            rows = im.rows;
            cols = im.cols;
        }

        sk_intel& operator=(const sk_intel& rhs) {
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
                mkl_free(this->matrix_data);
            this->rows = 0;
            this->cols = 0;
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

        sk_intel rand_n(const int row = -1, const int col = -1) {

            if(row != -1)
                this->rows = row;
            if(col != -1)
                this->cols = col;
            if(this->matrix_data)
                mkl_free(this->matrix_data);
            this->matrix_data = (float *) mkl_malloc(this->size() * sizeof(float), sizeof(float));
            int i;
            unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
            std::mt19937 gen(seed);
            std::normal_distribution<float> n;
            
            for(i = 0; i < this->size(); i++) {
                this->matrix_data[i] = n(gen);
                //std::cout << i << ": " << this->matrix_data[i] << std::endl;
            }
            
            return *this;
        }

        std::vector<int> dimensions() const {
            std::vector<int> ret = {this->rows, this->cols};
            return ret;
        }

        //gpu parallelizable
        sk_intel add(const sk_intel& rhs) const {
            if(this->dimensions() != rhs.dimensions())
                throw std::invalid_argument("mismatched dimensions");
            sk_intel temp(*this);
            std::cout << "temp: " << temp << std::endl;
            cblas_saxpy(temp.size(), 1, rhs.matrix_data, 1, temp.matrix_data, 1);
            return temp;
        }

        sk_intel subtract(const sk_intel& rhs) const {
            if(this->dimensions() != rhs.dimensions())
                throw std::invalid_argument("mismatched dimensions");
            sk_intel temp(*this);
            cblas_saxpy(temp.size(), -1, rhs.matrix_data, 1, temp.matrix_data, 1);
            return temp;
        }

        sk_intel mult(const sk_intel& rhs) const {
            if(this->cols != rhs.rows)
                throw std::invalid_argument("mismatched dimensions");
            sk_intel product(this->rows, rhs.cols);
            cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, this->rows, 
                    rhs.cols, this->cols, 1, this->matrix_data, rhs.rows, 
                    rhs.matrix_data, rhs.cols, 0, product.matrix_data, rhs.cols);
            std::cout << "p: " << product << std::endl;
            return product;
        }

        sk_intel elem_div(const float a) const {
            if(a == 0.0)
                throw std::invalid_argument("can't divide by zero");
            int i;
            sk_intel temp;
            for(i = 0; i < this->size(); i++)
                temp.matrix_data[i] = this->matrix_data[i] / a;
            return temp;
        }

        sk_intel operator+(const sk_intel& rhs) const {
            return this->add(rhs);
        }

        sk_intel& operator+=(const sk_intel& rhs) {
            *this = *this + rhs;
            return *this;
        }

        sk_intel operator-(const sk_intel& rhs) const {
            return this->subtract(rhs);
        }

        sk_intel& operator-=(const sk_intel& rhs) {
            *this = *this - rhs;
            return *this;
        }

        sk_intel operator*(const sk_intel& rhs) const {
            return this->mult(rhs);
        }

        sk_intel& operator*=(const sk_intel& rhs) {
            *this = this->mult(rhs); 
            return *this;
        }

        sk_intel operator/(const float a) {
            return this->elem_div(a);
        }

        sk_intel& operator/=(const float a) {
            *this = this->elem_div(a);
            return *this;
        }

        /*
        // Count Sketch
        std::vector<int>& flip_signs(const int col...) const;
        std::vector<int>& bucket(const int num_buckets) const;
        sk_intel& count_sketch() const;

        // Regression 
        sk_intel& concat(const sk_intel& col) const ;
        sk_intel& solve_x(const sk_intel& A, const sk_intel& B) const;

        // TODO: K-SVD 
        sk_intel& override_col(const int col, const sk_intel& B) const;
        std::vector<sk_intel> qr_decompose() const;
        sk_intel& svd() const;
         */

        bool operator==(const sk_intel& a) {
            if(a.rows != this->rows || a.cols != this->cols)
                return false;
            int i;
            for(i = 0; i < this->size(); i++)
                if(this->matrix_data[i] != a.matrix_data[i])
                    return false;
            return true;
        }

        friend std::ostream& operator<<(std::ostream&os, const sk_intel& im);

};

std::ostream& operator<<(std::ostream&os, const sk_intel& im) {
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
