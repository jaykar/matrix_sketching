#include <iostream>
#include <time.h>
#include <fstream>
#include <string>

#include "operations.hpp"
#include "boost.hpp"
// #include "armadillo.hpp"
// #include "intel.hpp"

using namespace std;

template <typename T>
void test_perf(ofstream& myfile){
    string header = "Size, Instantiation, Rand_Fill, Mult, Elem_Div, Concat, Subtract, Accumulate\n";
    cout << header;
    myfile << header;

    chrono::time_point<chrono::system_clock> start, end;

    for(int matrix_size = 100; matrix_size <= 5000; matrix_size+=100){
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
        std::cout << (chrono::system_clock::now() - start).count() << "\n";
        myfile << (chrono::system_clock::now() - start).count() << '\n';
    }
}

int main(int argc, char *argv[]){
    ofstream myfile;
    if(argc != 2 ){
        myfile.open ("construction_test1.csv");
    } else {
        myfile.open (argv[1]);
    }

    myfile << "Boost" << '\n';
    test_perf<sketchy::boost>(myfile);
    myfile << '\n';

    // myfile << "Armadillo" << '\n';
    // test_perf<sketchy::armadillo>(myfile);
    // myfile << '\n';

    // myfile << "Intel" << '\n';
    // test_perf<sketchy::intel>(myfile);
    // myfile << '\n';

    myfile.close();

    return 0;
}
