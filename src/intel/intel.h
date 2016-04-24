#include "Matrix.hpp"
#include "mkl.h"
#include <iostream>

class IntelMatrix : SKMatrix<IntelMatrix>  {
    private:
        double *data;
        int rows;
        int cols;

    public:
        IntelMatrix() {
            std::cout << "intel constructor" << std::endl;
            data = NULL;
            rows = 0;
            cols = 0;
        }

        IntelMatrix(const int row, const int col) {
            std::cout << "intel parametric constructor" << std::endl;
            data = (double *) mkl_malloc(row * col * sizeof(double), sizeof(double));
            rows = row;
            cols = col;
        }

        ~IntelMatrix() {
            std::cout << "intel destructor" << std::endl;
            mkl_free(data);
        }

        /*
           IntelMatrix(const IntelMatrix& im);
           IntelMatrix& operator=(const IntelMatrix& rhs); 
        //IntelMatrix(IntelMatrix&& other);

        //gpu parallaizable
        IntelMatrix operator+(const IntelMatrix& rhs) const;
        IntelMatrix operator-(const IntelMatrix& rhs) const;

        IntelMatrix& operator+=(const IntelMatrix& rhs); 
        IntelMatrix& operator-=(const IntelMatrix& rhs);

        IntelMatrix& operator*(const IntelMatrix& rhs); 
        IntelMatrix& operator*=(const IntelMatrix& rhs); 

        IntelMatrix& operator/(const IntelMatrix& rhs); 
        IntelMatrix& operator/=(const IntelMatrix& rhs); 

        int size() const;
        std::vector<int>& dimensions() const;

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
};
