#include "cuda.hpp"

int main() {
    sketchy::cuda a(3,4);
    a.rand_n(3,4);
    std::cout << a << std::endl;
    return 1;
}
