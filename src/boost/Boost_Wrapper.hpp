#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation.hpp>

using namespace boost::numeric::ublas;

template <typename T, typename V>
class Boost_Wrapper public SKMatrix<Boost_Wrapper>{
    private:
        T elems;
    public:
        Boost_Wrapper(const int row, const int col, ) {
            elems = scalar_matrix(row, col);
        }

        Boost_Wrapper(scalar_matrix<T>& m){
            this.elems = m;
        }

        int size() const {
            return elems.size1 * elems.size2
        }

        int num_rows(void) const {
            return elems.size1
        }

        int num_cols(void) const {
            return elems.size2
        }

        T data() const { return elems };

        std::vector<int>& dimensions() const {
            return std::vector<int>(elems.size1, elems.size2)
        }

        Boost_Wrapper<T>& mult(Boost_Wrapper<T>& rhs) const;

        Boost_Wrapper<float>& rand_n(const int row, const int col, const int mean, const int std) const;
        Boost_Wrapper<T>& elem_div(const T a) const;

        /* Count Sketch */
        Boost_Wrapper<T>& flip_signs(const std::vector<int> cols);

        /* Regression */
        Boost_Wrapper<T>& concat(const Boost_Wrapper<T>& col) const;
        Boost_Wrapper<T>& solve_x(const Boost_Wrapper<T>& A, const SKMatrix& B) const;

        /* TODO: K-SVD */
        Boost_Wrapper<T>& overridce_col(const int col, const SKMatrix& B) const;

        std::vector<Boost_Wrapper<T> > qr_decompose(Boost_Wrapper<V>& Q, Boost_Wrapper<T>& R) const;
};

template <typename T>
Boost_Wrapper<T>& mult(const Boost_Wrapper<T>& rhs) const {
    if(rhs.size1 != matrix.size_2) {
        std::cout << "Column of left matrix: " << elems.size_2 << " does not match row of right matrix: " << rhs.data.size_2 << "\n";
        throw;
    } else {
        T prod(elems.size_1, rhs.data.size_2);
        noalias(C) = prod(elems, rhs.data);
        return Boost_Wrapper<T>(prod);
    }
}

Boost_Wrapper<scalar_matrix<float> >& rand_n(const int row, const int col, const int mean, const int std) const {
    if(row < 0 || col < 0) {
        std::cout << "Column and row lengths must be non-negative integers" << '\n';
        throw;
    } else {
        std::default_random_engine generator;
        std::normal_distribution<float> distribution(mean, std);

        scalar_matrix<dfloat> rand_matrix(row, col);
        for (unsigned i = 0; i < rand_matrix.size1 (); ++ i){
            for (unsigned j = 0; j < rand_matrix.size2 (); ++ j){
                rand_matrix(i, j) = distribution(generator);
            }
        }
        return Boost_Wrapper<scalar_matrix<float> >(rand_matrix)
    }
}

template <typename T>
Boost_Wrapper<T>& elem_div(const T a) const {
    scalar_matrix<T> prod(elems);
    return Boost_Wrapper<scalar_matrix<T> >(prod / a);
}

template <typename T>
Boost_Wrapper<T>& flip_signs(const std::vector<int> cols) {
    srand( time(NULL) ); //Randomize seed initialization
    for(int col : cols){
        for(int row : elems.rows) {
            if(rand() % 2) {
                elems(row, col) = -1 * elems(row, col);
            }
        }
    }
    return *this;
}

// http://www.keithlantz.net/2012/05/qr-decomposition-using-householder-transformations/
template <typename T, typename V>
std::vector<Boost_Wrapper<T> > qr_decompose(Boost_Wrapper<V>& Q, Boost_Wrapper<T>& R) const {
    V mag;
    V alpha;

    scalar_matrix<T> u(this.num_rows, 1);
    scalar_matrix<T> v(this.num_rows, 1);

    scalar_matrix<V> P(this.num_rows, this.num_rows);
    identity_matrix I(this.num_rows);

    for (int i = 0; i < n; i++) {
        u *= 0;
        v *= 0;

        mag = 0.0;

        for (int j = i; j < this.num_cols; j++) {
            u(j,1) = R(j, i);
            mag += u(j,1) * u(j,1);
        }

        mag = std::sqrt(mag);

        alpha = u(i,1) < 0 ? mag : -1 * mag;

        mag = 0.0;

        for (int j = i; j < this.num_cols; j++) {
            v(j,1) = j == i ? u(j,1) + alpha : u(j,1);
            mag += v(j,1) * v(j,1);
        }
        mag = std::sqrt(mag);

        if (mag < 0.0000000001) continue;

        v /= mag;

        P = I - (v * boost::numeric::ublas::trans(v)) * 2.0;

        R = P * R;
        Q *= P;
    }
}