#include <vector>

template <typename T, class C>
class SKMatrix {
    public:
        SKMatrix(){}
        ~SKMatrix(){}
        virtual int size() const = 0;
        virtual std::vector<int> dimensions() const = 0;

        virtual C data() const = 0; //if we want this function then fix the other SKMatrix<T> instances

        // Gaussian projection
        //virtual T rand_n(int row, int col, int mean, int std) = 0; 
        virtual T rand_n(int row, int col) = 0; 
        virtual T mult(T& rhs) = 0; 
        virtual T elem_div(const double a) = 0;

        // Count Sketch 
        virtual std::vector<int> flip_signs()  = 0;
        virtual std::vector<int> bucket(const int num_buckets) = 0;

        //virtual T count_sketch() = 0;

        // Regression 
        virtual T concat(const T& col) const  = 0;
        virtual T solve_x(const T& B) = 0;
        virtual T get_cols(int start, int end) = 0; 
        virtual T get_col(int col_n) = 0; 
        virtual void t() = 0; 
        
        virtual T subtract(const T& rhs) = 0; 
        virtual double accumulate() = 0; 
        /*
        // TODO: K-SVD 
        virtual T& overridce_col(const int col, const SKMatrix& B) const = 0;
        virtual std::vector<T> qr_decompose() const = 0;
        */

};

