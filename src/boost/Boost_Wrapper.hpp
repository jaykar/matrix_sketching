#include "../interface/SKMatrix.hpp"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation.hpp>

namespace bnu = boost::numeric::ublas;

template <typename T, typename F>
class Boost_Wrapper : public SKMatrix<Boost_Wrapper<T, F>, T>{
    private:
        bnu::matrix<F> elems;

    public:
        Boost_Wrapper(bnu::matrix<F>& m){
            this->elems = m;
        }

        int size() const {
            return elems.size1 * elems.size2;
        }

        int num_rows(void) const {
            return elems.size1;
        }

        int num_cols(void) const {
            return elems.size2;
        }

        T data(void) const { return bnu::matrix<F>(elems); };

        Boost_Wrapper<T, F>& mult(const Boost_Wrapper<T, F>& rhs) const;

        Boost_Wrapper<T, F>& rand_n(const int row, const int col, const int mean, const int std) const;
        Boost_Wrapper<T, F>& elem_div(const F a) const;

        /* Count Sketch */
        Boost_Wrapper<T, F>& flip_signs(const std::vector<int> cols);

        /* Regression */
        Boost_Wrapper<T, F>& concat(const Boost_Wrapper<T, F>& col) const;
        Boost_Wrapper<T, F>& solve_x(const Boost_Wrapper<T, F>& B) const;

        /* TODO: K-SVD */
        Boost_Wrapper<T, F>& override_col(const int col, const Boost_Wrapper<T, F>& B) const;

        void qr_decompose(Boost_Wrapper<T, F>& Q, Boost_Wrapper<T, F>& R) const;
};

template <typename T, typename F>
Boost_Wrapper<T, F>& Boost_Wrapper<T, F>::mult(const Boost_Wrapper<T, F>& rhs) const {
    if(rhs.size1 != this->elems.size_2) {
        std::cout << "Column of left matrix: " << this->elems.size_2 << " does not match row of right matrix: " << rhs.data.size_2 << "\n";
        throw;
    } else {
        T prod(this->elems.size_1, rhs.data.size_2);
        bnu::noalias(prod) = bnu::prod(this->elems, rhs.data);
        return Boost_Wrapper<T, F>(prod);
    }
}

template <typename T, typename F>
Boost_Wrapper<T, F>& Boost_Wrapper<T, F>::rand_n(const int row, const int col, const int mean, const int std) const {
    if(row < 0 || col < 0) {
        std::cout << "Column and row lengths must be non-negative integers" << std::endl;
        throw;
    } else {
        std::default_random_engine generator;
        std::normal_distribution<float> distribution(mean, std);

        bnu::matrix<float> rand_matrix(row, col);
        for (unsigned i = 0; i < rand_matrix.size1 (); ++ i){
            for (unsigned j = 0; j < rand_matrix.size2 (); ++ j){
                rand_matrix(i, j) = distribution(generator);
            }
        }
        return Boost_Wrapper<bnu::matrix<T>, F>(rand_matrix);
    }
}

template <typename T, typename F>
Boost_Wrapper<T, F>& Boost_Wrapper<T, F>::elem_div(const F a) const {
    bnu::matrix<F> prod(this->elems);
    return Boost_Wrapper<bnu::matrix<F>, F>(prod / a);
}

template <typename T, typename F>
Boost_Wrapper<T, F>& Boost_Wrapper<T, F>::flip_signs(const std::vector<int> cols) {
    srand( time(NULL) ); //Randomize seed initialization
    for(int col : cols){
        if(rand() % 2){
            bnu::column(this->elems, col) *= -1;
        }
    }
    return *this;
}

template <typename T, typename F>
Boost_Wrapper<T, F>& Boost_Wrapper<T, F>::concat(const Boost_Wrapper<T, F>& new_col) const {
    if(new_col.num_rows != this->num_rows) {
        std::cout << "Number of rows don't match" << std::endl;
        throw;
    } else {
        bnu::matrix<T> concat_mat(this->elems);
        return concat_mat += new_col;
    }
}

template <typename T, typename F>
Boost_Wrapper<T, F>& Boost_Wrapper<T, F>::solve_x(const Boost_Wrapper<T, F>& B) const {
    bnu::matrix<T> Q(this->num_rows, this->num_cols);
    bnu::matrix<T> R(this->num_cols, this->num_cols);

    this->qr_decompose(Q, R);
}

// http://www.keithlantz.net/2012/05/qr-decomposition-using-householder-transformations/
template <typename T, typename F>
void Boost_Wrapper<T, F>::qr_decompose(Boost_Wrapper<T, F>& Q, Boost_Wrapper<T, F>& R) const {
    F mag;
    F alpha;

    bnu::matrix<F> u(this->num_rows, 1);
    bnu::matrix<F> v(this->num_rows, 1);

    bnu::identity_matrix<F> I(this->num_rows);
    R = this->elems;

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

        if (mag < 0.0000000001) continue;

        v /= mag;

        P = I - (v * bnu::trans(v)) * 2.0;

        Q *= P;
        R = P * R;
    }
}

template <typename T, typename F>
Boost_Wrapper<T, F>& Boost_Wrapper<T, F>::override_col(const int col_index, const Boost_Wrapper<T, F>& new_col) const {
    if(col_index < 0 || col_index >= this->num_cols) {
        std::cout << "Column index out of bound" << std::endl;
        throw;
    } else if (new_col.num_rows != this->num_rows) {
        std::cout << "Number of rows don't match" << std::endl;
        throw;
    } else {
        bnu::column(this->elems, col_index) = new_col.data();
    }
}