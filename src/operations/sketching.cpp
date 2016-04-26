#include "../armadillo/sk_arm.hpp"
// #include "../boost/sk_boost.hpp"
#include <iostream>
#include <math.h>

using namespace std;

template <class T>
T gaussian_projection(T input, int s){
    int n_col = input.dimensions()[1];
    auto a = T();
    a.rand_n(n_col, s);
    a.elem_div(sqrt(s*1.0));
    return input.mult(a);
}

template <class T>
T solve_x_sketch(T A, T B, int s){
    int n_col = A.dimensions()[1];
    auto concat_mat = A.concat(B);
    concat_mat.t();
    auto sketch = gaussian_projection<T>(concat_mat, s);
    sketch.t();
    auto A_sketch = sketch.get_cols(0, n_col -1);
    auto B_sketch = sketch.get_col(n_col);
    auto X_tilde  = A_sketch.solve_x(B_sketch);
    return X_tilde;
}

int main(){
    int total_size = 10;
    int sketch_size = 4;
    auto input = sk_boost<float>();
    input.rand_n(total_size, total_size);


    auto b = input;
    b.transpose();
    auto c = input.mult(b);
    cout << c.data() << endl;

    auto d = gaussian_projection(input, sketch_size);
    auto d_t = d;
    d.t();
    auto skt_d = d_t.mult(d);
    cout << skt_d.data() << endl;

    auto output = c.subtract(skt_d);
    cout << output.accumulate() << endl;
    /*
       auto b = Armadillo_Matrix();
       b.rand_n(100, 1);

       auto x = solve_x_sketch(input, b, 20);
       cout << x.data() << endl;
     */
    return 0;
}
