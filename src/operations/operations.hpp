#ifndef __OPERATIONS_H__
#define __OPERATIONS_H__

#include <math.h>
#include <vector>

template <typename SK>
SK gaussian_projection(const SK& input, const int s){
    int n_col = input.num_cols();
    SK a = SK();
    a.rand_n(n_col, s);
    a.elem_div(sqrt(s*1.0));
    return input.mult(a);
}

template <typename SK>
SK& lin_regress(const SK& A, const SK& B, const int s){
    int n_col = A.num_cols();

    auto concat_mat = A.concat(B);
    concat_mat.transpose();

    auto sketch = gaussian_projection<SK>(concat_mat, s);
    sketch.transpose();

    auto A_sketch = sketch.get_cols(0, n_col -1);
    auto B_sketch = sketch.get_col(n_col);
    auto X_tilde  = A_sketch.solve_x(B_sketch);

    return X_tilde;
}

template <typename SK>
SK count_sketch(const SK& A, const int num_buckets) {
    std::vector<std::vector<int> > buckets = A.bucket(num_buckets);
    std::vector<bool> flipped = A.flip_signs();

    SK sum(A.num_rows(), 0);

    for(auto hash_set : buckets){
        SK set(A.num_rows(), 1);
        for(int index : hash_set) {
            SK col = A.get_col(index);
            if(flipped[index]) {
                col.data() *= -1;
            }
            set.data() += col.data();
        }
        sum = sum.concat(set);
    }
    return sum;
}


template <typename SK>
std::vector<SK> k_svd(const SK&A, int k, int s){
    auto ans = std::vector<SK>(3); 
    //auto sketch = count_sketch(A, s); 
    auto sketch = gaussian_projection(A, s); 
    std::cout << sketch.size() << std::endl; 
    SK Q; 
    SK R; 
    sketch.qr_decompose(Q,R); 

    std::cout << Q.size() << std::endl; 
    std::cout << R.size() << std::endl; 
    auto Q_temp = Q; 
    Q_temp.transpose(); 
    Q_temp = Q_temp.mult(A); 
    auto svd_vec = Q_temp.svds(k); 
    auto U_tilde = Q.mult(svd_vec[0]); 
    ans[0] = U_tilde; 
    ans[1] = svd_vec[1]; 
    ans[2] = svd_vec[2]; 
    return ans; 
}
#endif
