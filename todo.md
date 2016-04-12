#Functions needed

##Sketching
* Gaussian Projection:
  * Access to size information, e.g. size(A)
  * Create matrix filled with random normal values of size m,n
  * Elementwise division
  * Matrix multiplication

* Count Sketch:
  * Flip the signs of entire columns in the matrix
  * Hashing function to discrete set (modulo is NOT fine)
  * Add columns in matrix
  * 

##Regression
* Approximate Matrix Solver, i.e. Ax = b
  * Concatenate matrix with vector
  * Call sketching algorithms
  * Built-in solver for Ax = b, Lapack also has built in solver, most matrix libraries will have one.

##K-SVD
* Approximate K-SVD
  * QR decomposition
  * Overwrite entire columns, e.g. A(:,5) = B
  * Overwrite columns in range, e.g. A(:, i-5 : i+5) = C
  * SVD
 
##Miscellaneous
 * +
 * -
 * /
 * *
 * =
 * Move
 * Constructors

