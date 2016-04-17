#import <vector>
template <typename T, typename C>
class SKMatrix {
    public:
        virtual int size(void) const = 0;
        virtual vector<int>& dimensions(void) const = 0;
        virtual SKMatrix<T> data() const = 0;
        virtual SKMatrix<T>& mult(SKMatrix<T>& rhs) const = 0;

        template<typename F>
        virtual SKMatrix<T>& col_op(const F& lambda, const int col) const = 0;

        /* Gaussian projection */
        virtual SKMatrix<double>& rand_n(int row, int col, int mean, int std) const = 0;
        virtual SKMatrix<T>& elem_div(const T a) const = 0;

        /* Count Sketch */
        virtual vector<int>& flip_signs(const int col...) const = 0;
        virtual vector<int>& bucket(const int num_buckets) const = 0;

        virtual SKMatrix<T>& count_sketch(void) const = 0;

        /* Regression */
        virtual SKMatrix<T>& concat(const SKMatrix<T>& col) const  = 0;
        virtual SKMatrix<T>& solve_x(const SKMatrix<T>& A, const SKMatrix& B) const = 0;

        /* TODO: K-SVD */
        virtual SKMatrix<T>& overridce_col(const int col, const SKMatrix& B) const = 0;
        virtual vector<SKMatrix<T> > qr_decompose() const = 0;
};

