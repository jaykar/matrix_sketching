#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <vector>
#include <iostream>
#include <random>

using namespace boost::numeric::ublas;

template <typename T, typename C>
class Boost_Wrapper public SKMatrix {
    private:
        matrix<T> elems;
    public:
        Boost_Wrapper(int row, int col) {
            elems = matrix(row, col);
        }

        Boost_Wrapper(matrix<T>& m){
            this.elems = m;
        }

        int size() const {
            return matrix.size1 * matrix.size2
        }

        SKMatrix<T> data() const { return matrix(elems) };

        std::vector<int>& dimensions() const {
            return std::vector<int>(elems.size1, elems.size2)
        }

        Boost_Wrapper<T>& mult(Boost_Wrapper<T>& rhs) const;

        template<typename F>
        SKMatrix<T>& col_op(const F& lambda, const int col) const;

        SKMatrix<T>& rand_n(int row, int col, int mean, int std) const;
        SKMatrix<T>& elem_div(const T a) const;

        /* Count Sketch */
        vector<int>& flip_signs(const int col...) const;
        vector<int>& bucket(const int num_buckets) const;

        SKMatrix<T>& count_sketch() const;

        /* Regression */
        SKMatrix<T>& concat(const SKMatrix<T>& col) const;
        SKMatrix<T>& solve_x(const SKMatrix<T>& A, const SKMatrix& B) const;

        /* TODO: K-SVD */
        SKMatrix<T>& overridce_col(const int col, const SKMatrix& B) const;
        vector<SKMatrix<T> > qr_decompose();
};

template <typename T>
Boost_Wrapper<T>& mult(const Boost_Wrapper<T>& rhs) const {
    if(rhs.size1 != matrix.size_2) {
        cout << "Column of left matrix: " << elems.size_2 << " does not match row of right matrix: " << rhs.data.size_2 << "";
        throw;
    } else {
        matrix<T> prod(elems.size_1, rhs.data.size_2);
        noalias(C) = prod(elems, rhs.data);
        return Boost_Wrapper(prod);
    }
}

SKMatrix<double>& rand_n(int row, int col, int mean, int std) const {
    if(row < 0 || col < 0) {
        cout << "Column and row lengths must be non-negative integers";
        throw;
    } else {
        std::default_random_engine generator;
        std::normal_distribution<double> distribution(mean, std);

        matrix<double> rand_matrix(row, col);
        for (unsigned i = 0; i < rand_matrix.size1 (); ++ i){
            for (unsigned j = 0; j < rand_matrix.size2 (); ++ j){
                rand_matrix(i, j) = distribution(generator);
            }
        }
    }
}

SKMatrix<T>& elem_div(const T a) const {
    matrix<T> prod(elems);
    return prod / a;
}