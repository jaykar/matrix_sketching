#include "cuda.hpp"

int main() {
    sketchy::cuda a(3,4);
    sketchy::cuda b = a;
    std::cout << a << std::endl;
    std::cout << b << std::endl;
    return 1;
}
