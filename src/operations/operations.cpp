#include "SKMatrix"
#include "operations.hpp"
#include <cstdarg>

template <typename SKMatrix>
SKMatrix& gaussian_projection(const SKMatrix& input, const int s){
    int n_col = input.num_cols();
    SKMatrix a = SKMatrix();
    a.rand_n(n_col, s);
    a.elem_div(sqrt(s*1.0));
    return input.mult(a);
}

template <typename SKMatrix>
SKMatrix& solve_x_sketch(const SKMatrix& A, const SKMatrix& B, const int s){
    int n_col = A.num_cols();
    auto concat_mat = A.concat(B);
    concat_mat.transpose();
    auto sketch = gaussian_projection<SKMatrix>(concat_mat, s);
    sketch.transpose();
    auto A_sketch = sketch.get_cols(0, n_col -1);
    auto B_sketch = sketch.get_col(n_col);
    auto X_tilde  = A_sketch.solve_x(B_sketch);
    return X_tilde;
}

template <typename SKMatrix>
SKMatrix count_sketch(const SKMatrix& A, const int num_buckets) {
    std::vector<std::vector<int> > buckets = bucket(num_buckets);
    for(int hash : buckets){

    }
}
