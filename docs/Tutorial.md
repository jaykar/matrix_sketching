
# Libraries Supported
1. [Boost](http://www.boost.org/doc/libs/1_60_0/libs/numeric/ublas/doc/)
2. [Armadillo](http://arma.sourceforge.net/docs.html)
3. [Intel MKL](https://software.intel.com/en-us/intel-mkl)
4. [CUDA](http://docs.nvidia.com/cuda/nvblas/)

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