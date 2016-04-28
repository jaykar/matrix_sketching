#include <iostream>
#include "operations.hpp"
#include "sketchy_boost.hpp"
// #include "sketchy_arm.hpp"

using namespace std;

namespace bnu = boost::numeric::ublas;

template<typename T>
T rand_n(int total_size){
    auto input = T();
    input.rand_n(total_size, total_size);
    return input;
}

template<typename T>
T transpose_mult(T input){
    auto b = input;
    b.transpose();
    return input.mult(b);
}

int main(){
    int total_size = 10;
    int sketchyetch_size = 4;

    // sketchy::arm input_arm;

    sketchy::intel intel_matrix();

    sketchy::arm arm_matrix(3,4);

    boost::numeric::ublas::matrix<float> mat;
    sketchy::boost boost_matrix(mat);


    // = rand_n<sketchy_boost>(total_size);
    // input_arm = rand_n<sketchy::arm>(total_size);

    cout << input_boost << endl;
    // cout << input_arm << endl;

    // sketchy::arm t_arm;
    // auto t_boost = transpose_mult(input_boost);
    // t_arm = transpose_mult(input_arm);

    cout << input_boost << endl;
    // cout << t_arm << endl;

    // input_arm = rand_n<sketchy::arm>(total_size*5);
    // auto ans = k_svd(input_arm, 10, 15);
    auto d = sketchy::ops::gaussian_projection<sketchy::boost>(input_boost, sketchyetch_size);
    /*

    cout << d << endl;
    auto d_t = d;
    d.transpose();

    auto sketchyt_d = d_t.mult(d);
    cout << sketchyt_d.data() << endl;

    auto output = c.subtract(sketchyt_d);
    cout << output.accumulate() << endl;

    auto x = count_sketchyetch<sketchy_boost>(input, sketchyetch_size);
    */
    return 0;
}
