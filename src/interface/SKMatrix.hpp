#import <vector>
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
        virtual C data() = 0; //if we want this function then fix the other SKMatrix<T> instances

        // Gaussian projection
        //virtual T rand_n(int row, int col, int mean, int std) = 0;
        virtual T rand_n(int row, int col) = 0;
        virtual T mult(T& rhs) = 0;
        virtual T elem_div(const double a) = 0;

        // Count Sketch
        virtual std::vector<int> flip_signs()  = 0;
        virtual std::vector<int> s(const int num_buckets) = 0;

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

        //virtual T count_sketch() = 0;

        // Regression
        virtual T concat(const T& col)  = 0;
        virtual T solve_x(const T& B) = 0;

        virtual T get_cols(int start, int end) = 0;
        virtual T get_col(int col_n) = 0;

        virtual void t() = 0;

        virtual T subtract(const T& rhs) = 0;
        virtual double accumulate() = 0;

        // TODO: K-SVD
        virtual T overridce_col(const int col, const SKMatrix& B) const = 0;
        virtual void qr_decompose(T& a, T& b) const = 0;
        virtual std::vector<T> qr_decompose() const = 0;
};

