#include "../interface/SKMatrix.hpp"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation.hpp>

namespace bnu = boost::numeric::ublas;

template <typename F>
class Boost_Wrapper: public SKMatrix<Boost_Wrapper<F>, bnu::matrix<F> >{
    private:
        bnu::matrix<F> elems;

    public:
        Boost_Wrapper<F>(bnu::matrix<F>& m){
            this->elems = m;
        }
        ~Boost_Wrapper<F>() = default;

        // Boost_Wrapper<F>(Boost_Wrapper<F>&& other){
        //     std::cout << "using move operator" << std::endl;
        //     this->matrix_data = other.matrix_data;
        //     //might not need to do this
        //     other.matrix_data = mat();
        // }

        // Boost_Wrapper<F>& operator=(Boost_Wrapper<F>&& other){
        //     std::cout << "using move operator" << std::endl;
        //     this->matrix_data = other.matrix_data;
        //     other.matrix_data = mat();
        //     return *this;
        // }

        int size() const {
            return elems.size1 * elems.size2;
        }

        int num_rows(void) const {
            return elems.size1;
        }

        int num_cols(void) const {
            return elems.size2;
        }

        bnu::matrix<F> data(void) const { return bnu::matrix<F>(elems); };

        Boost_Wrapper<F> mult(const Boost_Wrapper<F>& rhs) const;

        Boost_Wrapper<F> rand_n(const int row, const int col, const int mean, const int std) const;
        Boost_Wrapper<F> elem_div(const F a) const;

        /* Count Sketch */
        Boost_Wrapper<F> flip_signs(const std::vector<int> cols);

        /* Regression */
        Boost_Wrapper<F> concat(const Boost_Wrapper<F>& col) const;
        Boost_Wrapper<F> solve_x(const Boost_Wrapper<F>& B) const;

        /* FODO: K-SVD */
        Boost_Wrapper<F> override_col(const int col, const Boost_Wrapper<F>& new_col) const;

        void qr_decompose(Boost_Wrapper<F>& Q, Boost_Wrapper<F>& R) const;
};

template <typename F>
Boost_Wrapper<F> Boost_Wrapper<F>::mult(const Boost_Wrapper<F>& rhs) const {
    if(rhs.size1 != this->elems.size_2) {
        std::cout << "Column of left matrix: " << this->elems.size_2 << " does not match row of right matrix: " << rhs.data.size_2 << "\n";
        throw;
    } else {
        F prod(this->elems.size_1, rhs.data.size_2);
        bnu::noalias(prod) = bnu::prod(this->elems, rhs.data);
        return std::move(Boost_Wrapper<F>(prod));
    }
}

template <typename F>
Boost_Wrapper<F> Boost_Wrapper<F>::rand_n(const int row, const int col, const int mean, const int std) const {
    if(row < 0 || col < 0) {
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
    bnu::matrix<F> prod(this->elems);
    return std::move(Boost_Wrapper<F>(prod / a));
}

template <typename F>
Boost_Wrapper<F> Boost_Wrapper<F>::flip_signs(const std::vector<int> cols) {
    srand( time(NULL) ); //Randomize seed initialization
    for(int col : cols){
        if(rand() % 2){
            bnu::column(this->elems, col) *= -1;
        }
    }
    return std::move(this);
}

template <typename F>
Boost_Wrapper<F> Boost_Wrapper<F>::concat(const Boost_Wrapper<F>& new_col) const {
    if(new_col.num_rows != this->num_rows) {
        std::cout << "Number of rows don't match" << std::endl;
        throw;
    } else {
        bnu::matrix<F> concat_mat(this->elems);
        return std::move(Boost_Wrapper<F>(concat_mat += new_col));
    }
}

template <typename F>
Boost_Wrapper<F> Boost_Wrapper<F>::solve_x(const Boost_Wrapper<F>& B) const {
    bnu::matrix<F> Q(this->num_rows, this->num_cols);
    bnu::matrix<F> R(this->num_cols, this->num_cols);

    this->qr_decompose(Q, R);
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

template <typename F>
Boost_Wrapper<F> Boost_Wrapper<F>::override_col(const int col_index, const Boost_Wrapper<F>& new_col) const {
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