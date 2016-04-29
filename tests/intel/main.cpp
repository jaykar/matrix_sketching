#include "../../src/sk_intel.hpp"

int main() {
    sketchy::intel a(2,2);
    a.set(0,0,1.0);
    a.set(0,1,1.0);
    a.set(1,0,0.09);
    a.set(1,1,0.11);
    std::cout << "a: ";
    std::cout << a << std::endl;
    sketchy::intel b(2, 1);
    b.set(0, 0, 12000.0);
    b.set(0, 1, 1180.0);
    std::cout << "b: ";
    std::cout << b << std::endl;
    sketchy::intel c(b);
    std::cout << "c: ";
    std::cout << c << std::endl;
    sketchy::intel d = c;
    std::cout << "d: ";
    std::cout << d << std::endl;
    //a.rand_n(4, 4);
    std::cout << "a * X:" << std::endl;
    sketchy::intel X = a.solve_x(b);
    std::cout << a * b << std::endl;
//b.rand_n(4, 4);
    std::cout << "a: ";
    std::cout << a << std::endl;
    std::cout << "X: ";
    std::cout << a << std::endl;
    std::cout << a.accumulate() << std::endl;
    std::cout << "Xt: ";
    a.transpose();
    std::cout << a << std::endl;
    std::cout << a.accumulate() << std::endl;
    std::cout << "[a b]: ";
    sketchy::intel f = a.concat(b);
    std::cout << f << std:: endl;
    std::cout << f.get_cols(1,2) << std::endl;

    /*
    a -= b;
    std::cout << "a: ";
    a.identity(4);

    std::cout << a << std::endl;

    sketchy::intel product;
    product = a * b;
    
    std::cout << (product == b) << std::endl;

    //std::cout << c.data[2] << " " << a.data[2] + b.data[2] << std::endl;
    */
    return 0;
}
