#include "SKMatrix.hpp"
#include "mkl.h"
#include "mkl_blas.h"
#include <iostream>
#include <sstream>
#include <string>
#include <chrono>
#include <random>

class IntelMatrix : SKmatrix<IntelMatrix, float *, float>  {
    private:
        int rows;
        int cols;

    public:
        float *data;
        IntelMatrix() {
            //std::cout << "intel constructor" << std::endl;
            data = NULL;
            rows = 0;
            cols = 0;
        }

        IntelMatrix(const int row, const int col) {
            //std::cout << "intel parametric constructor" << std::endl;
            data = (float *) mkl_malloc(row * col * sizeof(float), sizeof(float));
            memset(data, 0, row * col * sizeof(float));
            rows = row;
            cols = col;
        }

        ~IntelMatrix() {
            //std::cout << "intel destructor: " << (long) data << std::endl;
            mkl_free(data);
        }

        IntelMatrix(const IntelMatrix& im) {
            //std::cout << "copy constructor" << std::endl;
            data = (float *) mkl_malloc(im.size() * sizeof(float), sizeof(float));
            cblas_scopy(im.size(), im.data, 1, data, 1);
            rows = im.rows;
            cols = im.cols;
        }

        /*
        void init() {
            if(!this->data) {
                this->data = (float *) mkl_malloc(sizeof(float), sizeof(float));
                this->rows = 1;
                this->cols = 1;
            }
            memset(this->data, 0, this->size() * sizeof(float));
        }
        */

        IntelMatrix& operator=(const IntelMatrix& rhs) {
            if(this->data)
                mkl_free(this->data);
            this->rows = rhs.rows;
            this->cols = rhs.cols;
            this->data = (float *) mkl_malloc(rhs.size() * sizeof(float), sizeof(float));
            cblas_scopy(rhs.size(), rhs.data, 1, this->data, 1);
            return *this;
        }

        void clear() {
            if(this->data)
                mkl_free(this->data);
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

        IntelMatrix rand_n(const int row = -1, const int col = -1) {

            if(row != -1)
                this->rows = row;
            if(col != -1)
                this->cols = col;
            if(this->data)
                mkl_free(this->data);
            this->data = (float *) mkl_malloc(this->size() * sizeof(float), sizeof(float));
            int i;
            unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
            std::mt19937 gen(seed);
            std::normal_distribution<float> n;
            
            for(i = 0; i < this->size(); i++) {
                this->data[i] = n(gen);
                //std::cout << i << ": " << this->data[i] << std::endl;
            }
            
            return *this;
        }

        std::vector<int> dimensions() const {
            std::vector<int> ret = {this->rows, this->cols};
            return ret;
        }

        /*
        //IntelMatrix(IntelMatrix&& other);
        */

        //gpu parallelizable
        IntelMatrix operator+(const IntelMatrix& rhs) const {
            if(this->dimensions() != rhs.dimensions())
                throw std::invalid_argument("mismatched dimensions");
            IntelMatrix temp(*this);
            std::cout << "temp: " << temp << std::endl;
            cblas_saxpy(temp.size(), 1, rhs.data, 1, temp.data, 1);
            return temp;
        }


        IntelMatrix subtract(const IntelMatrix& rhs) const {
            if(this->dimensions() != rhs.dimensions())
                throw std::invalid_argument("mismatched dimensions");
            IntelMatrix temp(*this);
            //std::cout << "temp: " << temp << std::endl;
            /*
            int i;
            for(i = 0; i < temp.size(); i++)
                std::cout << temp.data[i] - rhs.data[i] << '\t';
            std::cout <<std::endl;
            */
            cblas_saxpy(temp.size(), -1, rhs.data, 1, temp.data, 1);
            /*
            for(i = 0; i < temp.size(); i++)
                std::cout << temp.data[i] << '\t';
            std::cout << std::endl;
            */
            return temp;
        }

        IntelMatrix operator-(const IntelMatrix& rhs) const {
            return this->subtract(rhs);
        }

        IntelMatrix& operator+=(const IntelMatrix& rhs) {
            *this = *this + rhs;
            return *this;
        }

        IntelMatrix& operator-=(const IntelMatrix& rhs) {
            *this = *this - rhs;
            return *this;
        }

        /*
        IntelMatrix& operator*(const IntelMatrix& rhs); 
        IntelMatrix& operator*=(const IntelMatrix& rhs); 

        IntelMatrix& operator/(const IntelMatrix& rhs); 
        IntelMatrix& operator/=(const IntelMatrix& rhs); 

        IntelMatrix& mult(IntelMatrix& m) const;

        // Gaussian projection 
        IntelMatrix& rand_n(const int row, const int col) const;
        IntelMatrix& elem_div(const float a) const;

        // Count Sketch
        std::vector<int>& flip_signs(const int col...) const;
        std::vector<int>& bucket(const int num_buckets) const;
        IntelMatrix& count_sketch() const;

        // Regression 
        IntelMatrix& concat(const IntelMatrix& col) const ;
        IntelMatrix& solve_x(const IntelMatrix& A, const IntelMatrix& B) const;

        // TODO: K-SVD 
        IntelMatrix& override_col(const int col, const IntelMatrix& B) const;
        std::vector<IntelMatrix> qr_decompose() const;
        IntelMatrix& svd() const;
         */

        friend std::ostream& operator<<(std::ostream&os, const IntelMatrix& im);

};

std::ostream& operator<<(std::ostream&os, const IntelMatrix& im) {
    int i;
    std::ostringstream out;
    for(i = 0; i < im.size(); i++) {
        /*
        if(i > 0 && i % im.cols == 0)
            out << '\n';
            */
        out << std::to_string(im.data[i]);
        out << '\t';
    }
    os << out.str();
    return os;
}
