#include "sk_intel.hpp"

int main() {
    sk_intel a(2,2);
    a.set(0,0,1.0);
    a.set(0,1,1.0);
    a.set(1,0,0.09);
    a.set(1,1,0.11);
    std::cout << "a: ";
    std::cout << a << std::endl;
    sk_intel b(2, 1);
    b.set(0, 0, 12000.0);
    b.set(0, 1, 1180.0);
    std::cout << "b: ";
    std::cout << b << std::endl;
    sk_intel c(b);
    std::cout << "c: ";
    std::cout << c << std::endl;
    sk_intel d = c;
    std::cout << "d: ";
    std::cout << d << std::endl;
    //a.rand_n(4, 4);
    std::cout << "a * X:" << std::endl;
    sk_intel X = a.solve_x(b);
    std::cout << a * b << std::endl;
//b.rand_n(4, 4);
    std::cout << "a: ";
    std::cout << a << std::endl;
    std::cout << "b: ";
    std::cout << X << std::endl;
    /*
    a -= b;
    std::cout << "a: ";
    a.identity(4);

    std::cout << a << std::endl;

    sk_intel product;
    product = a * b;
    
    std::cout << (product == b) << std::endl;

    //std::cout << c.data[2] << " " << a.data[2] + b.data[2] << std::endl;
    */
    return 0;
}
