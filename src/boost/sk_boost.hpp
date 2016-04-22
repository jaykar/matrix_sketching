#ifndef __BOOST_WRAPPER_H__
#define __BOOST_WRAPPER_H__

#include "../interface/SKMatrix.hpp"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation.hpp>

namespace bnu = boost::numeric::ublas;

// use assert instead of throw

template <typename F>
class sk_boost: public SKMatrix<sk_boost<F>, bnu::matrix<F> >{
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

        int size() const { return matrix_data.size1() * matrix_data.size2(); }

        int num_rows(void) const { return matrix_data.size1(); }

        int num_cols(void) const { return matrix_data.size2(); }

        bnu::matrix<F> data(void) const { return bnu::matrix<F>(matrix_data); };

        sk_boost<F> mult(const sk_boost<F>& rhs) const;

        sk_boost<F> rand_n(const int row, const int col) const;
        sk_boost<F> elem_div(const F a) const;

        /* Count Sketch */
        // sk_boost<F> flip_signs(const std::vector<int> cols);

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
            bnu::matrix<F> temp(matrix_data);
            return bnu::sum(bnu::prod(bnu::scalar_vector<F>(temp.size1()), temp));
        }

        /* FODO: K-SVD */

        sk_boost<F> override_col(const int col, const sk_boost<F>& new_col) const;
        void qr_decompose(sk_boost<F>& Q, sk_boost<F>& R) const;
};

template <typename F>
sk_boost<F> sk_boost<F>::mult(const sk_boost<F>& rhs) const {
    if (rhs.num_rows() != this->num_cols()) {
        std::cout << "Column of left matrix: " << this->num_cols() << " does not match row of right matrix: " << rhs.num_rows() << "\n";
        throw;
    } else {
        bnu::matrix<F> prod(this->num_rows(), rhs.num_cols());
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
sk_boost<F> sk_boost<F>::elem_div(const F a) const {
    bnu::matrix<F> mat(this->matrix_data);
    bnu::matrix<F> result = mat/a;
    return std::move(sk_boost<F>(result));
}

// template <typename F>
// sk_boost<F> sk_boost<F>::flip_signs(const std::vector<int> cols) {
//     std::random_device rd;
//     std::mt19937 gen(rd());
//     std::uniform_int_distribution<> dis(1, 2);
//     for(int col : cols){
//         if (dis(gen) % 2){
//             bnu::column(this->matrix_data, col) *= -1;
//         }
//     }
//     return std::move(this);
// }

template <typename F>
sk_boost<F> sk_boost<F>::concat(const sk_boost<F>& new_col) const {
    if (new_col.num_rows() != this->num_rows()) {
        std::cout << "Number of rows don't match" << std::endl;
        throw;
    } else {
        bnu::matrix<F> concat_mat(this->matrix_data);
        concat_mat = new_col.data();
        return std::move(sk_boost<F>(concat_mat));
    }
}

template <typename F>
sk_boost<F> sk_boost<F>::solve_x(const sk_boost<F>& B) const {
    bnu::matrix<F> Aiv(this->num_rows(), this->num_cols());
    bnu::matrix<F> x(this->data());
    bnu::permutation_matrix<size_t> pm(x.size1());
    bnu::lu_factorize(x, pm);
    bnu::lu_substitute(x, pm, Aiv);
    return std::move(sk_boost<F>(x));
}

// http://www.keithlantz.net/2012/05/qr-decomposition-using-householder-transformations/
template <typename F>
void sk_boost<F>::qr_decompose(sk_boost<F>& Q, sk_boost<F>& R) const {
// template<class T>
// void TransposeMultiply(const ublas::vector<T>& vector, ublas::matrix<T>& result, size_t size) {
//     result.resize(size, size);
//     result.clear();

//     for(int row=0; row < vector.size(); ++row) {
//         for(int col=0; col < vector.size(); ++col) {
//             result(row,col) = vector(col) * vector(row);
//         }
//     }
// }

// template<class T>
// void HouseholderCornerSubstraction(ublas::matrix<T>& LeftLarge, const ublas::matrix<T>& RightSmall?) {
//     if(!((LeftLarge.size1() >= RightSmall.size1()) && (LeftLarge.size2() >= RightSmall.size2()))) {
//       cerr << "invalid matrix dimensions" << endl;
//       return;
//     }

//     size_t row_offset = LeftLarge.size2() - RightSmall.size2();
//     size_t col_offset = LeftLarge.size1() - RightSmall.size1();

//     for(unsigned int row = 0; row < RightSmall.size2(); ++row ) {
//         for(unsigned int col = 0; col < RightSmall.size1(); ++col )
//             LeftLarge(col_offset+col,row_offset+row) -= RightSmall?(col,row);
//         }
//     }
// }

// void HouseholderQR? (const ublas::matrix<T>& M, ublas::matrix<T>& Q, ublas::matrix<T>& R) {
//     if(matrix_data.num_rows != matrix_data.num_cols){
//         std::cout << "Matrix is not square" << endl;
//         throw;
//     }

//     int size = this->num_rows();
//     // init Matrices
//     bnu::matrix<F> H, HTemp;
//     HTemp = bnu::identity_matrix<F>(size);
//     bnu::matrix<F> q = bnu::identity_matrix<T>(size);
//     bnu::matrix<F> r(matrix_data);

//     // find Householder reflection matrices
//     for(int col = 0; col < size; col++) {
//       // create X vector
//         ublas::vector<T> RRowView = bnu::column(R,col);
//         std::vector_range<ublas::vector<T> > X2 (RRowView, range(col, size));
//         ublas::vector<T> X = X2;

//         if(X(0) >= 0) {
//             X(0) += norm_2(X);
//         }
//         else {
//             X(0) += -1*norm_2(X);
//         }

//         HTemp.resize(X.size(),X.size(),true);
//         TransposeMultiply(X, HTemp, X.size());
//         // HTemp = the 2UUt part of H
//         HTemp *= ( 2 / inner_prod(X,X) );
//         // H = I - 2UUt
//         H = identity_matrix<T>(size);
//         HouseholderCornerSubstraction(H,HTemp);
//         // add H to Q and R
//         q = prod(q, H);
//         r = prod(H, r);
//     }
// }

    F mag;
    F alpha;

    bnu::matrix<F> u(this->num_rows(), 1);
    bnu::matrix<F> v(this->num_rows(), 1);

    bnu::identity_matrix<F> I(this->num_rows());

    bnu::matrix<F> q = bnu::identity_matrix<F>(this->num_rows());;
    bnu::matrix<F> r(this->matrix_data);

    for (int i = 0; i < this->num_rows(); i++) {
        bnu::matrix<F> P(this->num_rows(), this->num_rows());

        u *= 0;
        v *= 0;

        mag = 0.0;

        for (int j = i; j < this->num_cols(); j++) {
            u(j,1) = r(j, i);
            mag += u(j,1) * u(j,1);
        }

        mag = std::sqrt(mag);

        alpha = u(i,1) < 0 ? mag : -1 * mag;

        mag = 0.0;

        for (int j = i; j < this->num_cols(); j++) {
            v(j,1) = j == i ? u(j,1) + alpha : u(j,1);
            mag += v(j,1) * v(j,1);
        }
        mag = std::sqrt(mag); // norm

        if  (mag < 0.0000000001) continue;

        v /= mag;

        P = I - (bnu::prod(v,bnu::trans(v))) * 2.0;

        q = bnu::prod(q, P);
        r = bnu::prod(P, r);
    }

    Q = sk_boost<F>(q);
    R = sk_boost<F>(r);
}

template <typename F>
sk_boost<F> sk_boost<F>::override_col(const int col_index, const sk_boost<F>& new_col) const {
    if (col_index < 0 || col_index >= this->num_cols()) {
        std::cout << "Column index out of bound" << std::endl;
        throw;
    } else if  (new_col.num_rows() != this->num_rows()) {
        std::cout << "Number of rows don't match" << std::endl;
        throw;
    } else {
        // bnu::column(this->matrix_data, col_index) = bnu::column(new_col.data());
        return *this;
    }
}


template <typename F>
sk_boost<F> sk_boost<F>::subtract(const sk_boost<F>& rhs) const {
    bnu::matrix<F> diff = this->matrix_data - rhs.data();
    return std::move(sk_boost<F>(diff));
}

template <typename F>
sk_boost<F> sk_boost<F>::get_col(int col_n) const {
    if(col_n < 0 || col_n >= this->num_cols()) {
        std::cout << "Column index out of bound" << std::endl;
        throw;
    } else {
        bnu::matrix<F> column(1, 2);
        bnu::column(this->matrix_data, col_n);
        return std::move(sk_boost<F>(column));
    }
};

#endif