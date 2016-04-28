#include "Matrix.hpp"
#include "mkl.h"
#include "mkl_blas.h"
#include <iostream>
#include <sstream>
#include <string>
#include <chrono>
#include <random>

namespace sketchy {
    class intel : SKMatrix<intel>  {
        private:
            int rows;
            int cols;

        public:
            double *data;
            intel() {
                //std::cout << "intel constructor" << std::endl;
                data = NULL;
                rows = 0;
                cols = 0;
            }

            intel(const int row, const int col) {
                //std::cout << "intel parametric constructor" << std::endl;
                data = (double *) mkl_malloc(row * col * sizeof(double), sizeof(double));
                memset(data, 0, row * col * sizeof(double));
                rows = row;
                cols = col;
            }

            ~intel() {
                //std::cout << "intel destructor: " << (long) data << std::endl;
                mkl_free(data);
            }

            intel(const intel& im) {
                //std::cout << "copy constructor" << std::endl;
                data = (double *) mkl_malloc(im.size() * sizeof(double), sizeof(double));
                cblas_dcopy(im.size(), im.data, 1, data, 1);
                rows = im.rows;
                cols = im.cols;
            }

            void init() {
                int i, j;
                if(!this->data) {
                    this->data = (double *) mkl_malloc(sizeof(double), sizeof(double));
                    this->rows = 1;
                    this->cols = 1;
                }
                memset(this->data, 0, this->size() * sizeof(double));
            }

            intel& operator=(const intel& rhs) {
                if(this->data)
                    mkl_free(this->data);
                this->rows = rhs.rows;
                this->cols = rhs.cols;
                this->data = (double *) mkl_malloc(rhs.size() * sizeof(double), sizeof(double));
                cblas_dcopy(rhs.size(), rhs.data, 1, this->data, 1);
                return *this;
            }

            int size() const {
                return this->rows * this->cols;
            }

            intel rand_n(int row = -1, int col = -1) {

                if(row != -1)
                    this->rows = row;
                if(col != -1)
                    this->cols = col;
                if(this->data)
                    mkl_free(this->data);
                this->data = (double *) mkl_malloc(this->size() * sizeof(double), sizeof(double));
                int i;
                unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
                std::mt19937 gen(seed);
                std::normal_distribution<double> n;

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
            //intel(intel&& other);
            */

            //gpu parallelizable
            intel operator+(const intel& rhs) const {
                if(this->dimensions() != rhs.dimensions())
                    throw std::invalid_argument("mismatched dimensions");
                intel temp(*this);
                std::cout << "temp: " << temp << std::endl;
                cblas_daxpy(temp.size(), 1, rhs.data, 1, temp.data, 1);
                return temp;
            }


            /*

            intel operator-(const intel& rhs) const;

            intel& operator+=(const intel& rhs);
            intel& operator-=(const intel& rhs);

            intel& operator*(const intel& rhs);
            intel& operator*=(const intel& rhs);

            intel& operator/(const intel& rhs);
            intel& operator/=(const intel& rhs);

            intel& mult(intel& m) const;

            // Gaussian projection
            intel& rand_n(const int row, const int col) const;
            intel& elem_div(const float a) const;

            // Count Sketch
            std::vector<int>& flip_signs(const int col...) const;
            std::vector<int>& bucket(const int num_buckets) const;
            intel& count_sketch() const;

            // Regression
            intel& concat(const intel& col) const ;
            intel& solve_x(const intel& A, const intel& B) const;

            // TODO: K-SVD
            intel& override_col(const int col, const intel& B) const;
            std::vector<intel> qr_decompose() const;
            intel& svd() const;
             */

            friend std::ostream& operator<<(std::ostream&os, const intel& im);
    };

    std::ostream& operator<<(std::ostream&os, const intel& im) {
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
}