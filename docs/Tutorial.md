# Tutorial

First download which ever of these implementations you would like to use 

1. [Boost](http://www.boost.org/doc/libs/1_60_0/libs/numeric/ublas/doc/)
2. [Armadillo](http://arma.sourceforge.net/docs.html)
3. [Intel MKL](https://software.intel.com/en-us/intel-mkl)
4. [NVBLAS](http://docs.nvidia.com/cuda/nvblas/)

## Linking
All of our interfaces are written as header files. As such, assuming
above libraries are installed, you just have to link our project directory 
along with whatever is needed to implement the matrix of your choice
```
gcc/clang++ -I ./sketchy $(BOOST_Path) $(ARMADILLO_PATH) $(INTEL_PATH) $(ETC)
```

## Example
Include header files as such
~~~{.c++}
#include "boost.hpp"
#include "armadillo.hpp"
#include "intel.hpp"
#include "nvblas.hpp"
#include "operations.hpp"
~~~

`SKMatrices` can be instantiated in three ways
~~~{.c++}
void sample {
    sketchy::intel intel_matrix();
    sketchy::boost boost_matrix(1, 2);

    armadillo::mat arma_mat;
    sketchy::armadillo arma_matrix(arma_mat);
}
~~~

Outputs of methods can be assigned and chained
~~~{.c++}
void sample {
    sketchy::boost mat1(5, 4);
    sketchy::boost mat2(4, 10);

    sketchy::boost mat3 = mat1.mult(mat2);
    sketchy::boost mat4 = mat1.subtract(mat2);

    sketchy::boost mat5 = mat1.mult(mat2.rand_n(4, 30));
}
~~~

Sketching operations are in `operations.hpp` header under namespace `sketchy::ops`. 

They are templated so explicit type needs to be provided

~~~{.c++}
#include "sketchy/operations.hpp"

void sample {
    sketchy::armadillo mat(20, 50);
    sketchy::ops::count_sketch<sketchy::armadillo>(mat, 10);
}
~~~

