#include <iostream>
#include <time.h>
#include <fstream>
#include <string>

#include "operations.hpp"
#include "sk_boost.hpp"
#include "sk_arm.hpp"
#include "sk_intel.hpp"

using namespace std;

template <typename T>
void test(ofstream& myfile) {
    string header = "Size, Instantiation, Rand_Fill, Mult, Elem_Div, Concat, Subtract, Accumulate, Gaussian Decomposition, Count Sketch\n";
    cout << header;
    myfile << header;

    chrono::time_point<chrono::system_clock> start, end;

    for(int matrix_size = 100; matrix_size <= 1000; matrix_size+=100){
        start = chrono::system_clock::now();
        T matrix(matrix_size, matrix_size);
        std::cout << matrix_size << " " << (chrono::system_clock::now() - start).count() << " ";
        myfile << matrix_size << "," << (chrono::system_clock::now() - start).count() << ',';

        start = chrono::system_clock::now();
        matrix.rand_n(matrix_size, matrix_size);
        std::cout << (chrono::system_clock::now() - start).count() << " ";
        myfile << (chrono::system_clock::now() - start).count() << ',';

        T matrix2(matrix_size, matrix_size);
        start = chrono::system_clock::now();
        matrix.mult(matrix2);
        std::cout << (chrono::system_clock::now() - start).count() << " ";
        myfile << (chrono::system_clock::now() - start).count() << ',';

        start = chrono::system_clock::now();
        matrix.elem_div(matrix_size/2);
        std::cout << (chrono::system_clock::now() - start).count() << " ";
        myfile << (chrono::system_clock::now() - start).count() << ',';

        start = chrono::system_clock::now();
        matrix.concat(matrix2);
        std::cout << (chrono::system_clock::now() - start).count() << " ";
        myfile << (chrono::system_clock::now() - start).count() << ',';

        start = chrono::system_clock::now();
        matrix.subtract(matrix2);
        std::cout << (chrono::system_clock::now() - start).count() << " ";
        myfile << (chrono::system_clock::now() - start).count() << ',';

        start = chrono::system_clock::now();
        matrix.accumulate();
        std::cout << (chrono::system_clock::now() - start).count() << " ";
        myfile << (chrono::system_clock::now() - start).count() << ',';

        // T Q(1, 1);
        // T R(1, 1);
        // start = chrono::system_clock::now();
        // matrix.qr_decompose(Q, R);
        // std::cout << (chrono::system_clock::now() - start).count() << ' ';
        // myfile << (chrono::system_clock::now() - start).count() << ',';

        // start = chrono::system_clock::now();
        // sketchy::ops::count_sketch<T>(matrix, matrix_size/2);
        // std::cout << (chrono::system_clock::now() - start).count() << '\n';
        // myfile << (chrono::system_clock::now() - start).count() << '\n';


        // start = chrono::system_clock::now();
        // sketchy::ops::gaussian_projection<T>(matrix, matrix_size/2);
        // std::cout << (chrono::system_clock::now() - start).count() << '\n';
        // myfile << (chrono::system_clock::now() - start).count() << '\n';
    }
}

int main(int argc, char *argv[]){
    ofstream myfile;

    if(argc != 2 ) {
        myfile.open ("construction_test1.csv");
    } else {
        myfile.open (argv[1]);
    }

    myfile << "Boost" << '\n';
    test<sketchy::boost>(myfile) ;
    myfile << '\n';

    myfile << "Armadillo" << '\n';
    test<sketchy::arm>(myfile) ;
    myfile << '\n';

    myfile << "Intel" << '\n';
    test<sketchy::intel>(myfile) ;
    myfile << '\n';

    myfile.close();

    return 0;
}
