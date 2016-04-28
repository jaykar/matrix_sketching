#include <iostream>
#include "operations.hpp"
#include "sk_boost.hpp"
#include "sk_arm.hpp"

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
    int sketch_size = 4;

    sk_arm input_arm;
    sk_boost input_boost;

    input_boost = rand_n<sk_boost>(total_size);
    input_arm = rand_n<sk_arm>(total_size);

    cout << input_boost << endl;
    cout << input_arm << endl;

    sk_arm t_arm;
    auto t_boost = transpose_mult(input_boost);
    t_arm = transpose_mult(input_arm);

    cout << t_boost << endl;
    cout << t_arm << endl;

    input_arm = rand_n<sk_arm>(total_size*5);
    auto ans = k_svd(input_arm, 10, 15);
    /*

    auto d = gaussian_projection<sk_boost>(input, sketch_size);
    cout << d << endl;
    auto d_t = d;
    d.transpose();

    auto skt_d = d_t.mult(d);
    cout << skt_d.data() << endl;

    auto output = c.subtract(skt_d);
    cout << output.accumulate() << endl;

    auto x = count_sketch<sk_boost>(input, sketch_size);
    */
    return 0;
}
