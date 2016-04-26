#ifndef __SKMATRIX_H__
#define __SKMATRIX_H__

#include <vector>
#include <random>
#include <cstdarg>
#include <time.h>
#include <iostream>

template <typename T, typename C>
class SKMatrix {
    protected:
        C matrix_data;

    public:
        SKMatrix(){}

        ~SKMatrix(){}
        T& operator=(const T& rhs){};
        T& operator=(const C& rhs){};
        // to do
        // T operator+(const T& rhs){};
        // T operator+(const C& rhs){};
        // T& operator+=(const T& rhs){};
        // T& operator+=(const C& rhs){};

        // T operator-(const T& rhs){};
        // T operator-(const C& rhs){};
        // T& operator-=(const T& rhs){};
        // T& operator-=(const C& rhs){};

        // T operator/(const T& rhs){};
        // T operator/(const C& rhs){};
        // T& operator/=(const T& rhs){};
        // T& operator/=(const C& rhs){};

        // T operator*(const T& rhs){};
        // T operator*(const C& rhs){};
        // T& operator*=(const T& rhs){};
        // T& operator*=(const C& rhs){};
        
        virtual void clear(void) = 0;
        virtual int size() const = 0;
        virtual int num_rows(void) const = 0;
        virtual int num_cols(void) const = 0;

        virtual C data(void) const = 0; //if we want this function then fix the other SKMatrix<T> instances
        virtual C& data(void) = 0; //if we want this function then fix the other SKMatrix<T> instances

        // Gaussian projection
        virtual T rand_n(const int row, const int col) = 0;
        virtual T mult(const T& rhs) const = 0;
        virtual T elem_div(const double a) const= 0;

        // Count Sketch
        std::vector<bool> flip_signs(){
            std::vector<bool> indices(this->num_cols());
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_int_distribution<> dis(1, 2);
            bool n_cols = this->matrix_data.n_cols;
            for(int i=0; i<n_cols; i++){
                int num = dis(gen);
                if (num == 2){
                    indices.push_back(true);
                }
                else{
                    indices.push_back(false);
                }
            }
            return indices;
        }

        std::vector<std::vector<int> > bucket(const int num_buckets) const {
            if(num_buckets > this->num_cols()){
                std::cout << "Number of buckets must be less than or equal to the number of columns";
                std::cout << '\n';
                throw;
            } else {
                std::vector<std::vector<int> > buckets(num_buckets);
                std::random_device rd;
                std::mt19937 gen(rd());
                std::uniform_int_distribution<> dis(0, num_buckets-1);

                for(int i = 0; i < this->num_cols(); i++){
                    int bucket = dis(gen);
                    buckets[bucket].push_back(i);
                }

                return buckets;
            }
        }

        // Regression
        virtual T concat(const T& col) const = 0;
        virtual T solve_x(const T& B) const = 0;

        virtual T get_cols(const int start, const int end) const = 0;
        virtual T get_col(const int col_n) const = 0;

        virtual void transpose() = 0;

        virtual T subtract(const T& rhs) const = 0;
        virtual double accumulate() const = 0;

        // TODO: K-SVD
        virtual void qr_decompose(T& a, T& b) const = 0;
};

#endif
