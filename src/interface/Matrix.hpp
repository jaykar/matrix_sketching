#include <tuple>
#include <vector>

template <class T>
class SKMatrix {
    public:
        SKMatrix();
        SKMatrix(const int row, const int col);
        ~SKMatrix();
        SKMatrix(const T& sk);
        T& operator=(const T& rhs); 
        //T(SKMatrix&& other);
        
        //gpu parallaizable
        virtual T operator+(const T& rhs) const = 0;
        virtual T operator-(const T& rhs) const = 0;

        virtual T& operator+=(const T& rhs) = 0; 
        virtual T& operator-=(const T& rhs) = 0;

        virtual T& operator*(const T& rhs) = 0; 
        virtual T& operator*=(const T& rhs) = 0; 

        virtual T& operator/(const T& rhs) = 0; 
        virtual T& operator/=(const T& rhs) = 0; 

        int size() const = 0;
        std::vector<int>& dimensions() const = 0;

        virtual T& mult(T& m) const = 0;

        //template<typename F>
        //virtual T& col_op(const F& lambda, const int col) const = 0;

        // Gaussian projection 
        virtual T& rand_n(const int row, const int col) const = 0;
        virtual T& elem_div(const float a) const = 0;

        // Count Sketch
        virtual std::vector<int>& flip_signs(const int col...) const = 0;
        virtual std::vector<int>& bucket(const int num_buckets) const = 0;
        virtual T& count_sketch() const = 0;

        // Regression 
        virtual T& concat(const T& col) const  = 0;
        virtual T& solve_x(const T& A, const T& B) const = 0;

        // TODO: K-SVD 
        virtual T& override_col(const int col, const T& B) const = 0;
        virtual std::vector<T> qr_decompose() const = 0;
        virtual T& svd() const = 0;
};
