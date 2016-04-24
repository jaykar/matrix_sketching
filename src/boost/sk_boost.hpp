#ifndef __BOOST_WRAPPER_H__
#define __BOOST_WRAPPER_H__

#include "../interface/SKMatrix.hpp"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/storage.hpp>

namespace bnu = boost::numeric::ublas;

// use assert instead of throw
template <typename F>
class sk_boost: public SKMatrix<sk_boost<F>, bnu::matrix<F>, F>{
    private:
        bnu::matrix<F> matrix_data;

    public:
        sk_boost<F>(const int row, const int col){
            this->matrix_data = bnu::matrix<F>(row, col);
        }

        sk_boost<F>(const bnu::matrix<F>& m){
            this->matrix_data = m;
        }

        ~sk_boost<F>() = default;

        sk_boost<F>& operator=(const sk_boost<F>& rhs){
            matrix_data = rhs.data();
            return *this;
        }

        sk_boost<F>& operator=(const bnu::matrix<F>& rhs){
            matrix_data = bnu::matrix<F>(rhs);
            return *this;
        }


        void clear() { matrix_data.clear(); }
        int size() const { return matrix_data.size1() * matrix_data.size2(); }

        int num_rows(void) const { return matrix_data.size1(); }
        int num_cols(void) const { return matrix_data.size2(); }

        bnu::matrix<F> data(void) const { return bnu::matrix<F>(matrix_data); };
        bnu::matrix<F>& data(void) { return matrix_data; };

        sk_boost<F> mult(const sk_boost<F>& rhs) const;

        sk_boost<F> rand_n(const int row, const int col) const;
        sk_boost<F> elem_div(const double a) const;


        /* Regression */
        sk_boost<F> concat(const sk_boost<F>& col) const;
        sk_boost<F> solve_x(const sk_boost<F>& B) const;

        sk_boost<F> get_cols(int start, int end) const;
        sk_boost<F> get_col(int col_n) const;

        void transpose() {
            matrix_data = bnu::trans(matrix_data);
        };

        sk_boost<F> subtract(const sk_boost<F>& rhs) const;

        double accumulate() const {
            bnu::matrix<F> temp(data());
            return bnu::sum(bnu::prod(bnu::scalar_vector<F>(temp.size1()), temp));
        }

        /* FODO: K-SVD */
        void qr_decompose(sk_boost<F>& Q, sk_boost<F>& R) const;
};

template <typename F>
sk_boost<F> sk_boost<F>::mult(const sk_boost<F>& rhs) const {
    if (rhs.num_rows() != num_cols()) {
        std::cout << "Column of left matrix: " << num_cols() << " does not match row of right matrix: " << rhs.num_rows() << "\n";
        throw;
    } else {
        bnu::matrix<F> prod(num_rows(), rhs.num_cols());
        bnu::noalias(prod) = bnu::prod(this->data(), rhs.data());
        return std::move(sk_boost<F>(prod));
    }
}

template <typename F>
sk_boost<F> sk_boost<F>::rand_n(const int row, const int col) const {
    if (row < 0 || col < 0) {
        std::cout << "Column and row lengths must be non-negative integers" << std::endl;
        throw;
    } else {
        std::default_random_engine generator;
        std::normal_distribution<F> distribution(0, 1);

        bnu::matrix<F> rand_matrix(row, col);
        for (unsigned i = 0; i < rand_matrix.size1 (); ++ i){
            for (unsigned j = 0; j < rand_matrix.size2 (); ++ j){
                rand_matrix(i, j) = distribution(generator);
            }
        }
        return std::move(sk_boost<F>(rand_matrix));
    }
}

template <typename F>
sk_boost<F> sk_boost<F>::elem_div(const double a) const {
    bnu::matrix<F> mat(data());
    bnu::matrix<F> result = mat/a;
    return std::move(sk_boost<F>(result));
}


template <typename F>
sk_boost<F> sk_boost<F>::concat(const sk_boost<F>& new_col) const {
    if (new_col.num_rows() != num_rows()) {
        std::cout << "Number of rows don't match" << std::endl;
        throw;
    } else {
        bnu::matrix<F> concat_mat(data());
        int end_index1 = concat_mat.size2();
        concat_mat.resize(concat_mat.size1(), end_index1 + new_col.data().size2(), true);
        std::cout << concat_mat.size2() << std::endl;

        for(int i = end_index1, j = 0; i < concat_mat.size2(); i++, j++){
            bnu::column(concat_mat, i) = bnu::column(new_col.data(), j);
        }

        return std::move(sk_boost<F>(concat_mat));
    }
}

template <typename F>
sk_boost<F> sk_boost<F>::solve_x(const sk_boost<F>& B) const {
    bnu::matrix<F> A(this->data());
    bnu::matrix<F> y = bnu::trans(B.data());

    bnu::vector<F> b(B.num_rows());
    std::copy(y.begin1(), y.end1(), b.begin());

    bnu::permutation_matrix<> piv(b.size());
    bnu::lu_factorize(A, piv);
    bnu::lu_substitute(A, piv, b);

    bnu::matrix<F> x(this->num_rows(), 1);
    std::copy(b.begin(), b.end(), x.begin1());
    return std::move(sk_boost<F>(x));
}

// http://www.keithlantz.net/2012/05/qr-decomposition-using-householder-transformations/
template <typename F>
void sk_boost<F>::qr_decompose(sk_boost<F>& Q, sk_boost<F>& R) const {
    F mag;
    F alpha;

    bnu::matrix<F> u(this->num_rows(), 1);
    bnu::matrix<F> v(this->num_rows(), 1);

    bnu::identity_matrix<F> I(this->num_rows());

    bnu::matrix<F> q = bnu::identity_matrix<F>(this->num_rows());;
    bnu::matrix<F> r(data());

    for (int i = 0; i < num_rows(); i++) {
        bnu::matrix<F> p(num_rows(), num_rows());

        u.clear();
        v.clear();

        mag = 0.0;

        for (int j = i; j < num_cols(); j++) {
            u(j,0) = r(j, i);
            mag += u(j,0) * u(j,0);
        }

        mag = std::sqrt(mag);

        alpha = -1 * mag;

        mag = 0.0;

        for (int j = i; j < num_cols(); j++) {
            v(j,0) = j == i ? u(j,0) + alpha : u(j,0);
            mag += v(j,0) * v(j,0);
        }

        mag = std::sqrt(mag); // norm

        if  (mag < 0.0000000001) continue;

        v /= mag;

        p = I - (bnu::prod(v, bnu::trans(v))) * 2.0;
        q = bnu::prod(q, p);
        r = bnu::prod(p, r);
    }

    Q = sk_boost<F>(q);
    R = sk_boost<F>(r);
}

template <typename F>
sk_boost<F> sk_boost<F>::subtract(const sk_boost<F>& rhs) const {
    bnu::matrix<F> diff = data() - rhs.data();
    return std::move(sk_boost<F>(diff));
}

template <typename F>
sk_boost<F> sk_boost<F>::get_col(int col_n) const {
    if(col_n < 0 || col_n >= num_cols()) {
        std::cout << "Column index out of bound" << std::endl;
        throw;
    } else {
        bnu::matrix<F> column = bnu::subrange(data(), 0, num_rows(), col_n, col_n+1);
        return std::move(sk_boost<F>(column));
    }
};

template <typename F>
sk_boost<F> sk_boost<F>::get_cols(int start, int end) const {
    if (start < 0 || end > num_cols()) {
        std::cout << "Column indices out of bound" << std::endl;
        throw;
    } else if (start > end){
        std::cout << "Start column greater than end column" << std::endl;
        throw;
    }else {
        bnu::matrix<F> columns = bnu::subrange(data(), 0, num_rows(), start, end);
        return sk_boost<F>(columns);
    }
}

#endif