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
class sk_boost: public SKMatrix<sk_boost, bnu::matrix<float> >{
    private:
        bnu::matrix<float> matrix_data;

    public:
        sk_boost(){
            this->matrix_data = bnu::matrix<float>();
        }

        sk_boost(const int row, const int col){
            this->matrix_data = bnu::matrix<float>(row, col, 0);
        }

        sk_boost(const bnu::matrix<float>& mat){
            this->matrix_data = mat;
        }

        ~sk_boost() = default;

        sk_boost& operator=(const sk_boost& rhs){
            matrix_data = rhs.data();
            return *this;
        }

        sk_boost& operator=(const bnu::matrix<float>& rhs){
            matrix_data = bnu::matrix<float>(rhs);
            return *this;
        }

        friend std::ostream& operator<<(std::ostream&os, const sk_boost& mat);

        void clear() { matrix_data.clear(); }
        int size() const { return matrix_data.size1() * matrix_data.size2(); }

        int num_rows(void) const { return matrix_data.size1(); }
        int num_cols(void) const { return matrix_data.size2(); }

        bnu::matrix<float> data(void) const { return bnu::matrix<float>(matrix_data); };
        bnu::matrix<float>& data(void) { return matrix_data; };

        sk_boost mult(const sk_boost& rhs) const;
        sk_boost mult_scalar(const sk_boost& rhs) const;

        sk_boost rand_n(const int row, const int col);
        sk_boost elem_div(const double a) const;


        /* Regression */
        sk_boost concat(const sk_boost& col) const;
        sk_boost solve_x(const sk_boost& B) const;

        sk_boost get_cols(const int start, const int end) const;
        sk_boost get_col(const int col_n) const;

        void transpose() {
            matrix_data = bnu::trans(matrix_data);
        };

        sk_boost subtract(const sk_boost& rhs) const;

        double accumulate() const {
            bnu::matrix<float> temp(data());
            return bnu::sum(bnu::prod(bnu::scalar_vector<float>(temp.size1()), temp));
        }

        /* floatODO: K-SVD */
        void qr_decompose(sk_boost& Q, sk_boost& R) const;
};


std::ostream& operator<<(std::ostream&os, const sk_boost& mat){
    std::ostringstream out;
    out << mat.data();
    os << out.str();
    return os;
}

sk_boost sk_boost::mult(const sk_boost& rhs) const {
    if (rhs.num_rows() != num_cols()) {
        std::cout << "Column of left matrix: " << num_cols() << " does not match row of right matrix: " << rhs.num_rows() << "\n";
        throw;
    } else {
        bnu::matrix<float> prod(num_rows(), rhs.num_cols());
        bnu::noalias(prod) = bnu::prod(this->data(), rhs.data());
        return std::move(sk_boost(prod));
    }
}

sk_boost sk_boost::rand_n(const int row, const int col) {
    if (row < 0 || col < 0) {
        std::cout << "Column and row lengths must be non-negative integers" << std::endl;
        throw;
    } else {
        std::default_random_engine generator;
        std::normal_distribution<float> distribution(0, 1);

        bnu::matrix<float> rand_matrix(row, col);
        for (unsigned i = 0; i < rand_matrix.size1 (); ++ i){
            for (unsigned j = 0; j < rand_matrix.size2 (); ++ j){
                rand_matrix(i, j) = distribution(generator);
            }
        }
        matrix_data = rand_matrix;
        return *this;
    }
}

sk_boost sk_boost::elem_div(const double a) const {
    bnu::matrix<float> mat(data());
    bnu::matrix<float> result = mat/a;
    return std::move(sk_boost(result));
}

sk_boost sk_boost::concat(const sk_boost& new_col) const {
    if (new_col.num_rows() != num_rows()) {
        std::cout << "Number of rows don't match" << std::endl;
        throw;
    } else {
        bnu::matrix<float> concat_mat(data());
        int end_index1 = concat_mat.size2();
        concat_mat.resize(concat_mat.size1(), end_index1 + new_col.data().size2(), true);

        for(unsigned int i = end_index1, j = 0; i < concat_mat.size2(); i++, j++){
            bnu::column(concat_mat, i) = bnu::column(new_col.data(), j);
        }

        return std::move(sk_boost(concat_mat));
    }
}

sk_boost sk_boost::solve_x(const sk_boost& B) const {
    bnu::matrix<float> A(this->data());
    bnu::matrix<float> y = bnu::trans(B.data());

    bnu::vector<float> b(B.num_rows());
    std::copy(y.begin1(), y.end1(), b.begin());

    bnu::permutation_matrix<> piv(b.size());
    bnu::lu_factorize(A, piv);
    bnu::lu_substitute(A, piv, b);

    bnu::matrix<float> x(this->num_rows(), 1);
    std::copy(b.begin(), b.end(), x.begin1());
    return std::move(sk_boost(x));
}

// http://www.keithlantz.net/2012/05/qr-decomposition-using-householder-transformations/
void sk_boost::qr_decompose(sk_boost& Q, sk_boost& R) const {
    float mag;
    float alpha;

    bnu::matrix<float> u(this->num_rows(), 1);
    bnu::matrix<float> v(this->num_rows(), 1);

    bnu::identity_matrix<float> I(this->num_rows());

    bnu::matrix<float> q = bnu::identity_matrix<float>(this->num_rows());;
    bnu::matrix<float> r(data());

    for (int i = 0; i < num_rows(); i++) {
        bnu::matrix<float> p(num_rows(), num_rows());

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

    Q = sk_boost(q);
    R = sk_boost(r);
}

sk_boost sk_boost::subtract(const sk_boost& rhs) const {
    bnu::matrix<float> diff = data() - rhs.data();
    return std::move(sk_boost(diff));
}

sk_boost sk_boost::get_col(const int col_n) const {
    if(col_n < 0 || col_n >= num_cols()) {
        std::cout << "Column index out of bound" << std::endl;
        throw;
    } else {
        bnu::matrix<float> column = bnu::subrange(data(), 0, num_rows(), col_n, col_n+1);
        return std::move(sk_boost(column));
    }
};

sk_boost sk_boost::get_cols(const int start, const int end) const {
    if (start < 0 || end > num_cols()) {
        std::cout << "Column indices out of bound" << std::endl;
        throw;
    } else if (start > end){
        std::cout << "Start column greater than end column" << std::endl;
        throw;
    }else {
        bnu::matrix<float> columns = bnu::subrange(data(), 0, num_rows(), start, end);
        return sk_boost(columns);
    }
}

#endif