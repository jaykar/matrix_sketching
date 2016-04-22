#ifndef __SKMATRIX_H__
#define __SKMATRIX_H__

#include <vector>
#include <random>
#include <cstdarg>
#include <time.h>
#include <iostream>

template <typename T, class C>
class SKMatrix {
    public:
        SKMatrix(){}
        ~SKMatrix(){}

        virtual int size() const = 0;
        virtual int num_rows(void) const = 0;
        virtual int num_cols(void) const = 0;
        virtual C data(void) const = 0; //if we want this function then fix the other SKMatrix<T> instances

        // Gaussian projection
        virtual T rand_n(const int row, const int col) const = 0;
        virtual T mult(const T& rhs) const = 0;
        virtual T elem_div(const double a) const= 0;

        // Count Sketch
        std::vector<int> flip_signs(){
            std::vector<int> indices(this->num_cols());
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_int_distribution<> dis(1, 2);
            int n_cols = this->matrix_data.n_cols;
            for(int i=0; i<n_cols; i++){
                int num = dis(gen);
                if (num == 2){
                    indices.push_back(-1);
                }
                else{
                    indices.push_back(1);
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

                // std::vector<int> indices(num_buckets);
                // std::random_device rd;
                // std::mt19937 gen(rd());
                // std::uniform_int_distribution<> dis(0, num_buckets-1);
                // int n_cols = this->matrix_data.n_cols;
                // for(int i=0; i<n_cols; i++){
                //     int num = dis(gen);
                //     indices.push_back(num);
                // }
                // return indices;
            }
        }

        // Regression
        virtual T concat(const T& col) const = 0;
        virtual T solve_x(const T& B) const = 0;

        virtual T get_cols(int start, int end) const = 0;
        virtual T get_col(int col_n) const = 0;

        virtual void transpose() = 0;

        virtual T subtract(const T& rhs) const = 0;
        virtual double accumulate() const = 0;

        // TODO: K-SVD
        virtual T override_col(const int col, const T& B) const = 0;
        virtual void qr_decompose(T& a, T& b) const = 0;
};

#endif