#ifndef __SKMATRIX_H__
#define __SKMATRIX_H__

#include <vector>
#include <random>
#include <cstdarg>
#include <time.h>
#include <iostream>
#include <exception>

template <typename T, typename C>
class SKMatrix {
    protected:
       /**
        * internal data structure to hold a matrix
        */
        C matrix_data;

    public:
       /**
        * Constructor
        */
        SKMatrix(){}

       /**
        * Destructor
        */
        ~SKMatrix(){}

       /**
        * Assignment operators
        * Allows assignments from both another SKMatrix object
        * @param rhs another SKMatrix instance
        * @see matrix_data
        */
        T& operator=(const T& rhs){};

        /**
        * Assignment operators
        * Allows assignments from other objects of same type as matrix_data
        * @param rhs another object of tyoe matrix_data
        * @see matrix_data
        */
        //T& operator=(const C& rhs){};

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

       /**
        * Clears matrix_data and sets its values to 0
        * @see matrix_data
        */
        virtual void clear(void) = 0;

       /**
        * Returns the number of elements in matrix_data
        * @see matrix_data
        * @return size number of elements
        */
        virtual int size() const = 0;

       /**
        * Returns the number of rows on matrix_data
        * @see matrix_data
        * @return number of rows
        */
        virtual int num_rows(void) const = 0;

       /**
        * Returns the number of columns on matrix_data
        * @see matrix_data
        * @return number of columns
        */
        virtual int num_cols(void) const = 0;

       /**
        * Returns a copy of matrix_data
        * for non-modifying purposes
        * @see matrix_data
        * @return copy of matrix_data
        */
        virtual C data(void) const = 0;

       /**
        * Create matrix filled with random normal values
        * with mean 0 and standard deviation 1
        * @param row number of rows for a new random matrix
        * @param columns number of columns for a new random matrix
        * @return a random row x col matrix
        */
        virtual T rand_n(const int row, const int col) = 0;

       /**
        * Multiplies current matrix with a new matrix
        * and returns the product
        * @param rhs right hand side of multiplication
        * @return product of two matrices
        */
        virtual T mult(const T& rhs) const = 0;

       /**
        * Divides internal matrix by a scalar
        * @param a divisor
        * @return matrix of quotients of scalar division
        */
        virtual T elem_div(const double a) const= 0;

       /**
        * Indicates which column needs to be flipped
        * with 50% chance
        * @return vector indicating indices of vectors to be flipped
        */
        std::vector<bool> flip_signs() const {
            std::vector<bool> indices(this->num_cols());
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_int_distribution<> dis(1, 2);
            int n_cols = num_cols();
            for(int i=0; i< n_cols; i++){
                int num = dis(gen);
                if (num == 2){
                    indices[i] = true;
                }
                else{
                    indices[i] = false;
                }
            }
            return indices;
        }

       /**
        * Buckets/hashes each column based on random uniform probability
        * @param num_buckets number of partitions
        * @return vector representing hashing of columns
        */
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

       /**
        * Concatenates/appends one matrix to the end of another
        * @param mat column to be concatenated
        * @return copy of current matrix with new matrix appended
        */
        virtual T concat(const T& mat) const = 0;

       /**
        * Solves a system of linear equation of form Ax=B
        * @param B scalar matrix
        * @return solution to this linear system
        */
        virtual T solve_x(const T& B) const = 0;

       /**
        * Retrieves a partiular column
        * @param col_n column index
        * @return copy of a particular column at given index
        */
        virtual T get_col(const int col_n) const = 0;

       /**
        * Retrieves columns in range
        * @param start left index
        * @param end right index
        * @return copy of columns in specified range
        */
        virtual T get_cols(const int start, const int end) const = 0;

       /**
        * Transposes matrix_data
        * @see matrix_data
        */
        virtual void transpose() = 0;

       /**
        * Subtracts matrix_data by that of right hand side
        * @param rhs right hand side matrix of subtraction
        * @see matrix_data
        * @return resultant matrix of subtraction
        */
        virtual T subtract(const T& rhs) const = 0;

       /**
        * Sums up all the values in matrix_data
        * @see matrix_data
        * @return sum of elements
        */
        virtual float accumulate() const = 0;

       /**
        * Performs QR decomposition of matrix_data
        * @param Q n by n Q matrix to be overwritten
        * @param R n by m R matrix to be overwritten
        */
        virtual void qr_decompose(T& Q, T& R) const = 0;

        virtual void svd(T& U, T& S, T& V, const int k) const = 0;
};

#endif
