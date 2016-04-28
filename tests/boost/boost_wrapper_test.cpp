#include <iostream>
#include <boost/numeric/ublas/matrix.hpp>
#include "sk_boost.hpp"

namespace bnu = boost::numeric::ublas;
using namespace std;

int main(){
  bnu::matrix<float> m1(5,5, 1.0);
  sketchy::boost mat1(m1);

  bnu::matrix<float> m2(5,5, 2.0);
  sketchy::boost mat2(m2);

  cout << "mat 1" << endl;
  cout << mat1 << endl;

  cout << "mat 2" << endl;
  cout << mat2 << endl;

  cout << "mat 1 subtract mat 2" << endl;
  cout << mat1.subtract(mat2) << endl;

  cout << "mat 1 concat mat 2" << endl;
  cout << mat1.concat(mat2) << endl;
  cout << mat1 <<endl;
  cout << "mat 1 mult mat 2" << endl;
  cout << mat1.mult(mat2) << endl;

  auto buckets = mat1.bucket(3);
  cout << buckets[0].size() << endl;
  cout << buckets[1].size() << endl;
  cout << buckets[2].size() << endl;

  cout << "rand_n(2, 2)" << endl;
  cout << mat1.rand_n(2, 2) << endl;

  cout << "mat1 elem_div 1.5" << endl;
  cout << mat1.elem_div(1.5) << endl;

  cout<< mat1 << endl;
  cout << "mat1 get_col 1" << endl;
  cout << mat1.get_col(1) << endl;

  cout << "mat1 get_col 2, 4" << endl;
  cout << mat1.get_cols(0,1) << endl;

  cout << "mat1 accumulate" << endl;
  cout << mat1.accumulate() << endl;
  bnu::matrix<double> B(3, 3);
  B(0,0) = 12;
  B(0,1) = -51;
  B(0,2) = 4;

  B(1,0) = 6;
  B(1,1) = 167;
  B(1,2) = -68;

  B(2,0) = -4;
  B(2,1) = 24;
  B(2,2) = -41;
  cout << B << endl;

  bnu::matrix<double> q(3, 3);
  bnu::matrix<double> r(3, 3);

  sketchy::boost matB(B);
  sketchy::boost Q(q);
  sketchy::boost R(r);
  cout << endl;
  matB.qr_decompose(Q, R);

  cout << Q << endl;
  cout << R << endl;

  bnu::matrix<float> D(2,1);
  D(0,0) = 5;
  D(1,0) = 7;
  sketchy::boost matD(D);

  bnu::matrix<float> E(2,2);
  E(0,0) = 3;
  E(0,1) = 4;

  E(1,0) = 2;
  E(1,1) = -1;
  sketchy::boost matE(E);

  cout << "Ax = B" << endl;

  cout << matE.solve_x(matD) << endl;

  cout << matB << endl;
  matB.svd(matE, matD, mat2, 1);

  return 0;
}