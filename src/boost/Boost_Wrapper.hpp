#ifndef __BOOST_WRAPPER_H__
#define __BOOST_WRAPPER_H__

#include "../interface/SKMatrix.hpp"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation.hpp>

namespace bnu = boost::numeric::ublas;

template <typename F>
class Boost_Wrapper: public SKMatrix<Boost_Wrapper<F>, bnu::matrix<F> >{
    private:
        bnu::matrix<F> matrix_data;

    public:
        Boost_Wrapper<F>(int row, int col){
            this->matrix_data = bnu::matrix<F>(row, col);
        }

        Boost_Wrapper<F>(bnu::matrix<F>& m){
            this->matrix_data = m;
        }

        ~Boost_Wrapper<F>() = default;

        Boost_Wrapper<F>(Boost_Wrapper<F>&& other){
            matrix_data = std::move(other.data());
            return *this;
        }

        Boost_Wrapper<F>& operator=(Boost_Wrapper<F>&& other){
            matrix_data = std::move(other.data());
            return *this;
        }

        int size() const { return matrix_data.size1 * matrix_data.size2; }

        int num_rows(void) const { return matrix_data.size1; }

        int num_cols(void) const { return matrix_data.size2; }

        bnu::matrix<F> data(void) const { return bnu::matrix<F>(matrix_data); };

        Boost_Wrapper<F> mult(const Boost_Wrapper<F>& rhs) const;

        Boost_Wrapper<F> rand_n(const int row, const int col, const int mean, const int std) const;
        Boost_Wrapper<F> elem_div(const F a) const;

        /* Count Sketch */
        Boost_Wrapper<F> flip_signs(const std::vector<int> cols);

        /* Regression */
        Boost_Wrapper<F> concat(const Boost_Wrapper<F>& col) const;
        Boost_Wrapper<F> solve_x(const Boost_Wrapper<F>& B) const;
        Boost_Wrapper<F> get_cols(int start, int end) = 0;
        Boost_Wrapper<F> get_col(int col_n) = 0;

        void t() = 0;

        Boost_Wrapper<F> subtract(const Boost_Wrapper<F>& rhs) = 0;
        double accumulate() = 0;

        /* FODO: K-SVD */

        Boost_Wrapper<F> override_col(const int col, const Boost_Wrapper<F>& new_col) const;
        void qr_decompose(Boost_Wrapper<F>& Q, Boost_Wrapper<F>& R) const;
};

template <typename F>
Boost_Wrapper<F> Boost_Wrapper<F>::mult(const Boost_Wrapper<F>& rhs) const {
    if (rhs.size1 != this->matrix_data.size_2) {
        std::cout << "Column of left matrix: " << this->matrix_data.size_2 << " does not match row of right matrix: " << rhs.data.size_2 << "\n";
        throw;
    } else {
        F prod(this->matrix_data.size_1, rhs.data.size_2);
        bnu::noalias(prod) = bnu::prod(this->matrix_data, rhs.data);
        return std::move(Boost_Wrapper<F>(prod));
    }
}

template <typename F>
Boost_Wrapper<F> Boost_Wrapper<F>::rand_n(const int row, const int col, const int mean, const int std) const {
    if (row < 0 || col < 0) {
        std::cout << "Column and row lengths must be non-negative integers" << std::endl;
        throw;
    } else {
        std::default_random_engine generator;
        std::normal_distribution<F> distribution(mean, std);

        bnu::matrix<F> rand_matrix(row, col);
        for (unsigned i = 0; i < rand_matrix.size1 (); ++ i){
            for (unsigned j = 0; j < rand_matrix.size2 (); ++ j){
                rand_matrix(i, j) = distribution(generator);
            }
        }
        return std::move(Boost_Wrapper<F>(rand_matrix));
    }
}

template <typename F>
Boost_Wrapper<F> Boost_Wrapper<F>::elem_div(const F a) const {
    bnu::matrix<F> prod(this->matrix_data);
    return std::move(Boost_Wrapper<F>(prod / a));
}

template <typename F>
Boost_Wrapper<F> Boost_Wrapper<F>::flip_signs(const std::vector<int> cols) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(1, 2);
    for(int col : cols){
        if (dis(gen) % 2){
            bnu::column(this->matrix_data, col) *= -1;
        }
    }
    return std::move(this);
}

template <typename F>
Boost_Wrapper<F> Boost_Wrapper<F>::concat(const Boost_Wrapper<F>& new_col) const {
    if (new_col.num_rows != this->num_rows) {
        std::cout << "Number of rows don't match" << std::endl;
        throw;
    } else {
        bnu::matrix<F> concat_mat(this->matrix_data);
        return std::move(Boost_Wrapper<F>(concat_mat += new_col));
    }
}

template <typename F>
Boost_Wrapper<F> Boost_Wrapper<F>::solve_x(const Boost_Wrapper<F>& B) const {
    bnu::matrix<F> x(this->matrix_data);
    bnu::permutation_matrix<size_t> pm(x.size1());
    lu_factorize(x, pm);
    lu_substitute(x, pm, B->matrix_data);
    return std::move(Boost_Wrapper<F>(x));
}

// http://www.keithlantz.net/2012/05/qr-decomposition-using-householder-transformations/
//
template <typename F>
void Boost_Wrapper<F>::qr_decompose(Boost_Wrapper<F>& Q, Boost_Wrapper<F>& R) const {
    F mag;
    F alpha;

    bnu::matrix<F> u(this->num_rows, 1);
    bnu::matrix<F> v(this->num_rows, 1);

    bnu::identity_matrix<F> I(this->num_rows);
    R = this->matrix_data;

    for (int i = 0; i < this->num_rows; i++) {
        bnu::matrix<F> P(this->num_rows, this->num_rows);

        u *= 0;
        v *= 0;

        mag = 0.0;

        for (int j = i; j < this->num_cols; j++) {
            u(j,1) = R(j, i);
            mag += u(j,1) * u(j,1);
        }

        mag = std::sqrt(mag);

        alpha = u(i,1) < 0 ? mag : -1 * mag;

        mag = 0.0;

        for (int j = i; j < this->num_cols; j++) {
            v(j,1) = j == i ? u(j,1) + alpha : u(j,1);
            mag += v(j,1) * v(j,1);
        }
        mag = std::sqrt(mag); // norm

        if  (mag < 0.0000000001) continue;

        v /= mag;

        P = I - (v * bnu::trans(v)) * 2.0;

        Q *= P;
        R = P * R;
    }
}

template <typename F>
Boost_Wrapper<F> Boost_Wrapper<F>::override_col(const int col_index, const Boost_Wrapper<F>& new_col) const {
    if (col_index < 0 || col_index >= this->num_cols) {
        std::cout << "Column index out of bound" << std::endl;
        throw;
    } else if  (new_col.num_rows != this->num_rows) {
        std::cout << "Number of rows don't match" << std::endl;
        throw;
    } else {
        bnu::column(this->matrix_data, col_index) = new_col.data();
        return std::move(this);
    }
}

#endif