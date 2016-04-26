#include "sk_intel.h"

int main() {
    sk_intel a;
    std::cout << "a: ";
    std::cout << a << std::endl;
    sk_intel b(4, 3);
    std::cout << "b: ";
    std::cout << b << std::endl;
    sk_intel c(b);
    std::cout << "c: ";
    std::cout << c << std::endl;
    sk_intel d = c;
    std::cout << "d: ";
    std::cout << d << std::endl;
    a.rand_n(4, 4);
    b.rand_n(4, 4);
    std::cout << "a: ";
    std::cout << a << std::endl;
    std::cout << "b: ";
    std::cout << b << std::endl;
    a -= b;
    std::cout << "a: ";
    a.identity(4);


    std::cout << a << std::endl;

    sk_intel product;
    product = a * b;
    
    std::cout << (product == b) << std::endl;

    //std::cout << c.data[2] << " " << a.data[2] + b.data[2] << std::endl;
    return 0;
}
