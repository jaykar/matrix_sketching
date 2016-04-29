
# Example
~~~{.c++}
#include "skecthy/boost.hpp"
#include "skecthy/armadillo.hpp"
#include "skecthy/intel.hpp"
#include "skecthy/nvblas.hpp"
#include "skecthy/nvblas.hpp"

#include <armadillo>

int main() {
    sketchy::intel intel_matrix();
    sketchy::boost boost_matrix(1, 2);

    armadillo::mat arma_mat;
    sketchy::armadillo arma_matrix(arma_mat);

    sketchy::ops::count_sketch<sketchy::boost>(boost_matrix, 10);
}
~~~