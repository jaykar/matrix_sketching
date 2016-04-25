#include "intel.h"

int main() {
    IntelMatrix a;
    a.init();
    std::cout << "a: ";
    std::cout << a << std::endl;
    IntelMatrix b(4, 3);
    std::cout << "b: ";
    std::cout << b << std::endl;
    IntelMatrix c(b);
    std::cout << "c: ";
    std::cout << c << std::endl;
    IntelMatrix d = c;
    std::cout << "d: ";
    std::cout << d << std::endl;
    a.rand_n(3, 4);
    b.rand_n(3, 4);
    std::cout << "a: ";
    std::cout << a << std::endl;
    std::cout << "b: ";
    std::cout << b << std::endl;
    c = a + b;
    std::cout << "c: ";

    std::cout << c << std::endl;

    //std::cout << c.data[2] << " " << a.data[2] + b.data[2] << std::endl;
    return 0;
}
