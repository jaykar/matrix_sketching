template <class T>  
class Matrix {
    public:
        virtual int size() const = 0;
        virtual vector<int>& dimensions() const = 0;
        virtual Matrix<T>& mult(Matrix<T>& m) const = 0;

        template<typename F>
            virtual Matrix<T>& col_op(const F& lambda, const int col) const = 0;

        /* Gaussian projection */
        virtual Matrix<T>& rand_n(const Container& dimensions) const = 0;
        virtual Matrix<T>& elem_div(const T a) const = 0;

        /* Count Sketch */
        virtual vector<int>& flip_signs(const int col...) const = 0;
        virtual vector<int>& bucket(const int num_buckets) const = 0;

        virtual Matrix<T>& count_sketch() const = 0;

        /* Regression */
        virtual Matrix<T>& concat(const Matrix<T>& col) const  = 0;
        virtual Matrix<T>& solve_x(const Matrix<T>& A, const Matrix& B) const = 0;
        
        /* TODO: K-SVD */

}
