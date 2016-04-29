Adam Incera (aji2112)
Jaykar Nayeck (jan2150)
Howon Byun (hb2458)

Needs
1. boost
2. armadillo
3. intel mkl
4. CUDA

# Sketchy

Writing GPU code is difficult. We want to write a library that allows C++ developers to easily perform matrix sketching on a GPU. Matrix sketching is a new and exciting topic in Machine Learning. The basic concept in sketching is that we can perform matrix operations such as inversions, svd, and eigenvalue decomposition much faster by forgoing some accuracy. However there are theoretical guarantees that our approximation is within some set threshold from the true value. Many of these normal matrix operations perform in O(n3) time, and as such do not scale well. By parallelizing them and using the new algorithms on a GPU, we will hopefully drastically speed up matrix decomposition time. 

# What we will Build
We will first create a small matrix library in C++ that wraps NVIDIA’s CUDA toolkit. This library will wrap specific functions in the cuBLAS library. The cuBLAS library is NVIDIA’s library for matrix operations on the GPU. It is based on the more general LAPACK library which has years of development focusing on matrix operations. With GPU programming, we need to allocate memory on the CPU, send it to the GPU, perform some operation in GPU, then finally return it back to CPU. The explicit memory allocation and deallocation is a constant source of bugs. Our sketching library will use RAII principles to abstract the memory allocations and deallocations. Our sketching library will call these internal classes that wrap cuBLAS. 

Our first tests will be to compare other matrix libraries with our wrappers. So we might perform matrix multiplication using Boost and compare the answer with ours. We will also make sure that we’re not leaking memory on the GPU, and that we have successfully implemented the move operator for our matrices. We will be using github to host our code. We will have a stable branch, a dev branch, and each of the team members will have their own branch. 

To demonstrate this library, our “hello world” program is linear regression on a large matrix using the Ordinary Least Squares solution. The solution is given by (X^T * X)^-1 * (X^T Y). The inverse operation will be much more costly without sketching, so hopefully our library will be able to reflect the theory. The program will look something like this: 

X = gpu_matrix(from_csv_file(‘features.csv’));
Y = gpu_matrix(from_csv_file(‘labels.csv’));
Regression_coef = pinv(X) * (X.T * Y);



Tools and Libraries we will Use
Our primary dependency will be NVIDIA’s CUDA toolkit. From this, we will primarily depend on cuBLAS library and cuRAND. 

We will be doing most of our development on Adam’s desktop, which has an NVIDIA GTX 960 graphics card. This card has 2GB of memory, meaning we will be able to test on matrices up to approximately 700x700 floats. If we find that we need more GPU memory, we can also test on a computer in Jaykar’s lab which has an NVIDIA Tesla K80 with 24GB of memory.

Adam will write the GPU wrappers and Jaykar will write the sketching algorithm. Because our library is so dependent on previous features, we will set milestones to have completed core parts: 

Week 1: Initial framework set up so we can perform basic CUDA operations such as transfer memory, and perform basic matrix operations. We will split writing the wrappers, Adam might write subtract and divide, while Jaykar writes add and multiply. We can both work on the wrapper for the more complex functions such as inv, and svd. We will be continually testing. 

Week 2: Matrix sketching algorithm

Week 3: Matrix Inversion using sketching

Week 4: Performance analysis

Release 0.8: Matrix Sketching, Wrappers for cuBLAS, easy interface to load csv files onto GPU
Release 1.0: Matrix Inversion using sketching, Linear Regression
Release 1.2: SVD using sketching, maybe Kernel PCA (we need to figure what that is)


Matrix: general, interface for other matrix libraries to plug into

## Dependencies
http://arma.sourceforge.net/
https://github.com/google/googletest
