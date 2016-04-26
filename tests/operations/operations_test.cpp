#include <iostream>
#include "operations.hpp"
#include "sk_boost.hpp"

using namespace std;

namespace bnu = boost::numeric::ublas;

int main(){
    int total_size = 10;
    int sketch_size = 4;
    auto input = sk_boost();
    input.rand_n(total_size, total_size);

    auto b = input;
    b.transpose();
    auto c = input.mult(b);

    auto d = gaussian_projection<sk_boost>(input, sketch_size);
    cout << d << endl;
    auto d_t = d;
    d.transpose();

    auto skt_d = d_t.mult(d);
    cout << skt_d.data() << endl;

    auto output = c.subtract(skt_d);
    cout << output.accumulate() << endl;

    auto x = count_sketch<sk_boost>(input, sketch_size);
    return 0;
}