#ifndef __SKMATRIX_H__
#define __SKMATRIX_H__

#include <vector>
#include <random>
#include <cstdarg>
#include <time.h>
#include <iostream>

template <typename t, typename c, typename f>
class SKmatrix {
    protected:
        c matrix_data;
    public:
        SKmatrix(){}
        ~SKmatrix(){}
        t& operator=(const t& rhs){};
        t& operator=(const c& rhs){};

        // to do
        t operator+(const t& rhs){};
        // t operator+(const c& rhs){};
        t& operator+=(const t& rhs){};
        // t& operator+=(const c& rhs){};

        t operator-(const t& rhs){};
        // t operator-(const c& rhs){};
        t& operator-=(const t& rhs){};
        // t& operator-=(const c& rhs){};

        // t operator/(const t& rhs){};
        // t operator/(const c& rhs){};
        // t& operator/=(const t& rhs){};
        // t& operator/=(const c& rhs){};

        // t operator*(const t& rhs){};
        // t operator*(const c& rhs){};
        // t& operator*=(const t& rhs){};
        // t& operator*=(const c& rhs){};

        virtual void clear(void) = 0;
        virtual int size() const = 0;
        virtual int num_rows(void) const = 0;
        virtual int num_cols(void) const = 0;

        // virtual c data(void) const = 0; //if we want this function then fix the other SKMatrix<T> instances
        // virtual c& data(void) = 0; //if we want this function then fix the other SKMatrix<T> instances

        // gaussian projection
        virtual t rand_n(const int row, const int col) = 0;
        // virtual t mult(const t& rhs) const = 0;
        //virtual t elem_div(const double a) const= 0;

        /*
        // count sketch
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
                std::cout << "number of buckets must be less than or equal to the number of columns";
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
        */

        // regression
        //virtual t concat(const t& col) const = 0;
        //virtual t solve_x(const t& b) const = 0;

        //virtual t get_cols(int start, int end) const = 0;
        //virtual t get_col(int col_n) const = 0;

        //virtual void transpose() = 0;

        virtual t subtract(const t& rhs) const = 0;
        //virtual double accumulate() const = 0;

        /*
        // TODO: K-SVD
        virtual void qr_decompose(t& a, t& b) const = 0;
        */
};

#endif
