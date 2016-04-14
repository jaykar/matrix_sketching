#include <tuple>

template <class T>
class SKMatrix {
    public:
        SKMatrix();
        SKMatrix(const int row, const int col);
        ~SKMatrix();
        SKMatrix(const SKMatrix& sk);

        SKMatrix& operator=(const SKMatrix& rhs);
        SKMatrix(SKMatrix&& other);

        /* gpu parallaizable */
        virtual SKMatrix operator+(const SKMatrix& lhs, const SKMatrix& rhs) const = 0;
        virtual SKMatrix operator-(const SKMatrix& lhs, const SKMatrix& rhs) const = 0;

        virtual SKMatrix& operator+=(const SKMatrix& rhs)
        virtual SKMatrix& operator-=(const SKMatrix& rhs)

        virtual SKMatrix& operator*(const SKMatrix& lhs, const SKMatrix& rhs)
        virtual SKMatrix& operator*=(const SKMatrix& rhs)

        virtual SKMatrix& operator/(const SKMatrix& lhs, const SKMatrix& rhs)
        virtual SKMatrix& operator/=(const SKMatrix& rhs)

        int size() const = 0;
        vector<int>& dimensions() const = 0;

        virtual SKMatrix<T>& mult(SKMatrix<T>& m) const = 0;

        template<typename F>
        virtual SKMatrix<T>& col_op(const F& lambda, const int col) const = 0;

        /* Gaussian projection */
        virtual SKMatrix<T>& rand_n(const Container& dimensions) const = 0;
        virtual SKMatrix<T>& elem_div(const T a) const = 0;

        /* Count Sketch */
        virtual vector<int>& flip_signs(const int col...) const = 0;
        virtual vector<int>& bucket(const int num_buckets) const = 0;

        virtual SKMatrix<T>& count_sketch() const = 0;

        /* Regression */
        virtual SKMatrix<T>& concat(const SKMatrix<T>& col) const  = 0;
        virtual SKMatrix<T>& solve_x(const SKMatrix<T>& A, const SKMatrix& B) const = 0;

        /* TODO: K-SVD */
        virtual SKMatrix<T>& overridce_col(const int col, const SKMatrix& B) const = 0;
        virtual vector<SKMatrix<T> > qr_decompose() const = 0;
        virtual SKMatrix<T>& svd() const = 0;
};