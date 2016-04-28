#include "SKMatrix.hpp"
#include "mkl.h"
//#include "mkl_scalapack.h"
//#include "mkl_lapacke.h"
#include <sstream>
#include <chrono>

class sk_intel : SKMatrix<sk_intel, float *>  {
    private:
        int rows;
        int cols;

    public:
        sk_intel() {
            matrix_data = NULL;
            rows = 0;
            cols = 0;
        }

        /*
        sk_intel& operator=(const float * rhs) {
            this->matrix_data = rhs;
        }
        */

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

        void set(int r, int c, float a) {
            this->matrix_data[r*this->rows + c] = a;
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

        sk_intel mult(const sk_intel& rhs) const {
            if(this->cols != rhs.rows)
                throw std::invalid_argument("mismatched dimensions");
            sk_intel product(this->rows, rhs.cols);
            cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, this->rows, 
                    rhs.cols, this->cols, 1, this->matrix_data, rhs.rows, 
                    rhs.matrix_data, rhs.cols, 0, product.matrix_data, rhs.cols);
            //std::cout << "p: " << product << std::endl;
            return product;
        }

        sk_intel elem_div(const double a) const {
            if(a == 0.0)
                throw std::invalid_argument("can't divide by zero");
            int i;
            sk_intel temp;
            for(i = 0; i < this->size(); i++)
                temp.matrix_data[i] = this->matrix_data[i] / a;
            return temp;
        }

        //TODO
        sk_intel concat(const sk_intel& col) const {
            if(this->rows != col.rows)
                throw std::invalid_argument("mismatched dimensions in concat");

            sk_intel ret(this->rows, this->cols + col.cols);
            float *current = ret.matrix_data;
            int i;
            for(i = 0; i < this->rows; i++) {
                cblas_scopy(this->cols, current, 1, this->matrix_data, 1);
                current += this->cols;
                cblas_scopy(col.cols, current, 1, col.matrix_data, 1);
                current += col.cols;
            }
            return ret;
        }

        sk_intel solve_x(const sk_intel& B) const {
            if(this->rows != this->cols)
                throw std::invalid_argument("A must be square in Ax = B");
            if(this->cols != B.rows)
                throw std::invalid_argument("B.rows must equal A.cols in Ax = B");
            if(B.cols != 1)
                throw std::invalid_argument("B must be nx1 vector in Ax = B");
            sk_intel X(B.rows, 1);
            MKL_INT lda = this->rows; // rows in A
            MKL_INT n = this->cols; // cols in A
            MKL_INT nrhs = 1; // cols in B
            MKL_INT ldb = nrhs;
            MKL_INT ipiv[n];

            sk_intel copy(B);
            LAPACKE_sgesv(LAPACK_ROW_MAJOR, n, nrhs, this->matrix_data, lda, ipiv, copy.matrix_data, ldb);

            return copy;
        }

        //TODO
        sk_intel get_col(const int n) const {
            return get_cols(n, n);
        }

        //TODO
        sk_intel get_cols(const int start, const int end) const {
            return sk_intel();
        }

        void transpose() {
            mkl_simatcopy('r', 't', this->rows, this->cols, 1.0, this->matrix_data, 
                    this->cols, this->rows);
            int i = this->rows;
            this->rows = this->cols;
            this->cols = i;
        }

        sk_intel subtract(const sk_intel& rhs) const {
            if(this->dimensions() != rhs.dimensions())
                throw std::invalid_argument("mismatched dimensions");
            sk_intel temp(*this);
            cblas_saxpy(temp.size(), -1, rhs.matrix_data, 1, temp.matrix_data, 1);
            return temp;
        }

        //TODO
        double accumulate() const {
            return 0.0;
        }

        //TODO
        void qr_decompose(sk_intel& Q, sk_intel& R) const {
        }

        //TODO
        std::vector<sk_intel> svds(const int k) const {
            return std::vector<sk_intel>();
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

        // Regression 
        */


        /*
        // TODO: K-SVD 
        sk_intel& override_col(const int col, const sk_intel& B) const;
         */

        friend bool operator==(const sk_intel& lhs, const sk_intel& rhs);
        friend std::ostream& operator<<(std::ostream&os, const sk_intel& im);
};

bool operator==(const sk_intel& lhs, const sk_intel& rhs) {
    if(rhs.rows != lhs.rows || rhs.cols != lhs.cols)
        return false;
    int i;
    for(i = 0; i < lhs.size(); i++)
        if(lhs.matrix_data[i] != rhs.matrix_data[i])
            return false;
    return true;
}

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
