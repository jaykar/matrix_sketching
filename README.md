# Sketchy

Sketchy is a high level library that acts as a wrapper for various
matrix implementations in C++ designed to support matrix sketching operations
on both CPU and GPU backed libraries.

Matrix sketching is a new and exciting topic in Machine Learning. The basic concept in sketching is that we can perform matrix operations such as inversions, svd, and eigenvalue decomposition much faster by forgoing some accuracy. However there are theoretical guarantees that our approximation is within some set threshold from the true value.

Many of these normal matrix operations perform in O(n^3) time, and as such do not scale well. By parallelizing them and using the new algorithms on a GPU, we will hopefully drastically speed up matrix decomposition time. 

# Libraries Supported
1. [Boost](http://www.boost.org/doc/libs/1_60_0/libs/numeric/ublas/doc/)
2. [Armadillo](http://arma.sourceforge.net/docs.html)
3. [Intel MKL](https://software.intel.com/en-us/intel-mkl)
4. [CUDA](http://docs.nvidia.com/cuda/nvblas/)

# Linking
All of our interfaces are written as header files. As such, assuming
above libraries are installed, you just have to link our sketchy directory 

Release 0.8: Matrix Sketching, Wrappers for cuBLAS, easy interface to load csv files onto GPU
Release 1.0: Matrix Inversion using sketching, Linear Regression
Release 1.2: SVD using sketching, maybe Kernel PCA (we need to figure what that is)

# Collaborators
Adam Incera (aji2112) </br>
Jaykar Nayeck (jan2150) </br>
Howon Byun (hb2458) </br>