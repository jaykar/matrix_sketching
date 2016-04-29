# Architecture
![Architecture](./images/architecture.png)

`SKMatrix`, which stands for SKetching Matrix, is the core interface that allows different matrix libraries to share methods in `sketchy/operations.hpp` for matrix sketching, such as

#### Organization / Structure
```bash
.
├── docs             # Doxygen and Markdown Documentations
├── sketchy          # Sketchy interface and wrappers
│   ├── SKMatrix     # primary Interface
│   ├── armadillo    # wrapper for Armadillo
│   ├── boost        # wrapper for oost
│   ├── intel        # wrapper for Intel MKL 
│   ├── cuda         # wrapper for Nvidea cuda
├── tests            # testing directory
└── Doxyfile         # configuration for Doxygen
```

1. [Gaussian Projection](https://en.wikipedia.org/wiki/Random_projection)
2. [K-SVD](http://www.cs.technion.ac.il/~elad/publications/journals/2004/32_KSVD_IEEE_TSP.pdf)
3. [Linear Regression](http://researcher.watson.ibm.com/researcher/files/us-dpwoodru/journal.pdf)
4. [Count Sketch](https://www.cs.rutgers.edu/~farach/pubs/FrequentStream.pdf)

This library has four wrapper classes based on `SKMatrix`. It currently supports 3 CPU matrix libraries

1. [Boost](http://www.boost.org/doc/libs/1_60_0/libs/numeric/ublas/doc/)
2. [Armadillo](http://arma.sourceforge.net/docs.html)
3. [Intel MKL](https://software.intel.com/en-us/intel-mkl)

and one GPU backed implementation

4. [NVBLAS](http://docs.nvidia.com/cuda/nvblas/)

One can easily create a wrapper for other matrix libraries by extending `SKMatrix` interface.